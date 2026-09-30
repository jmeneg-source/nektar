///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCSpecies.cpp
//
// See PEMFCSpecies.h. Part of the nektar++_fuel_cells project (EPUSP/RCGI).
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <PEMFCSolver/EquationSystems/PEMFCSpecies.h>

namespace Nektar
{
std::string PEMFCSpecies::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "PEMFCSpecies", PEMFCSpecies::create);

PEMFCSpecies::PEMFCSpecies(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph)
{
}

void PEMFCSpecies::v_InitObject(bool DeclareFields)
{
    AdvectionSystem::v_InitObject(DeclareFields);

    // Transport parameters (SI unless the case is nondimensional).
    m_session->LoadParameter("epsilon", m_epsilon, 1.0e-5);
    m_session->LoadParameter("rho", m_rho, 1.0);
    m_session->LoadParameter("Htot", m_Htot, 1.38e-3);
    m_session->LoadParameter("MO2", m_MO2, 0.032);
    m_session->LoadParameter("MH2O", m_MH2O, 0.018);
    int ne;
    m_session->LoadParameter("ne", ne, 4);
    m_ne = ne;

    // Kinetics parameters (defaults: Razmara et al. 2024, Table 2.1).
    m_session->LoadParameter("i0", m_kinetics.m_i0, 0.17);
    m_session->LoadParameter("alphaC", m_kinetics.m_alphaC, 0.5);
    m_session->LoadParameter("Tcell", m_kinetics.m_T, 343.15);
    m_session->LoadParameter("Voc", m_kinetics.m_Voc, 1.22);
    m_session->LoadParameter("Vcell", m_kinetics.m_Vcell, 0.4);
    m_session->LoadParameter("Rlump", m_kinetics.m_Rlump, 7.0e-5);
    m_session->LoadParameter("CO2ref", m_kinetics.m_CO2ref, 30.0);

    // Variable indices.
    m_idxO2  = -1;
    m_idxH2O = -1;
    for (int i = 0; i < m_session->GetVariables().size(); ++i)
    {
        if (m_session->GetVariable(i) == "wO2")
            m_idxO2 = i;
        if (m_session->GetVariable(i) == "wH2O")
            m_idxH2O = i;
    }
    ASSERTL0(m_idxO2 >= 0, "Variable wO2 not found in session.");
    ASSERTL0(m_idxH2O >= 0, "Variable wH2O not found in session.");

    // Frozen advection velocity from session function.
    m_velocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    std::vector<std::string> vel;
    vel.push_back("Vx");
    vel.push_back("Vy");
    vel.push_back("Vz");
    vel.resize(m_spacedim);
    GetFunction("AdvectionVelocity")->Evaluate(vel, m_velocity);

    m_iloc = Array<OneD, NekDouble>(GetNpoints(), 0.0);

    // Optional design-dependent diffusivity factor f_D(x) (thesis
    // Eq. 2.26: 0.19703 in ribs, 0.92167 in channels). If the session
    // defines a function "DiffusivityFactor", the implicit diffusion
    // solves div(f_D grad w) with variable coefficients; otherwise
    // f_D = 1 and the constant-coefficient path is used.
    m_hasVarDiff = m_session->DefinesFunction("DiffusivityFactor");
    if (m_hasVarDiff)
    {
        const int nq = GetNpoints();
        m_diffFactor = Array<OneD, NekDouble>(nq, 1.0);
        GetFunction("DiffusivityFactor")
            ->Evaluate(m_session->GetVariable(0), m_diffFactor);
        ASSERTL0(Vmath::Vmin(nq, m_diffFactor, 1) > 0.0,
                 "DiffusivityFactor must be positive.");
    }

    // Advection operator (CG, non-conservative), diffusion implicit.
    ASSERTL0(m_projectionType == MultiRegions::eGalerkin,
             "PEMFCSpecies stage 1 supports continuous Galerkin only.");
    ASSERTL0(!m_explicitDiffusion,
             "PEMFCSpecies requires implicit diffusion (IMEX).");

    std::string advName;
    m_session->LoadSolverInfo("AdvectionType", advName, "NonConservative");
    m_advObject =
        SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
    m_advObject->SetFluxVector(&PEMFCSpecies::GetFluxVector, this);

    m_ode.DefineImplicitSolve(&PEMFCSpecies::DoImplicitSolve, this);
    m_ode.DefineProjection(&PEMFCSpecies::DoOdeProjection, this);
    m_ode.DefineOdeRhs(&PEMFCSpecies::DoOdeRhs, this);
}

PEMFCSpecies::~PEMFCSpecies()
{
}

void PEMFCSpecies::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    ASSERTL1(flux[0].size() == m_velocity.size(),
             "Dimension of flux array and velocity array do not match");

    const int nq = m_fields[0]->GetNpoints();

    for (int i = 0; i < flux.size(); ++i)
    {
        for (int j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1, flux[i][j], 1);
        }
    }
}

void PEMFCSpecies::ComputeReaction(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    const int nq        = GetNpoints();
    const double F      = 96485.3;
    const double volO2  = -1.0 * m_MO2 / (m_ne * F) / (m_rho * m_Htot);
    const double volH2O = +2.0 * m_MH2O / (m_ne * F) / (m_rho * m_Htot);

    for (int q = 0; q < nq; ++q)
    {
        // O2 molar concentration from mass fraction (constant rho).
        const double CO2 = m_rho * inarray[m_idxO2][q] / m_MO2;
        const double il  = m_kinetics.Solve(CO2);
        m_iloc[q]        = il;
        outarray[m_idxO2][q] += volO2 * il;
        outarray[m_idxH2O][q] += volH2O * il;
    }
}

void PEMFCSpecies::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    const int nVariables   = inarray.size();
    const int nSolutionPts = GetNpoints();

    m_advObject->Advect(nVariables, m_fields, m_velocity, inarray, outarray,
                        time);

    for (int i = 0; i < nVariables; ++i)
    {
        Vmath::Neg(nSolutionPts, outarray[i], 1);
    }

    ComputeReaction(inarray, outarray);
}

void PEMFCSpecies::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    SetBoundaryConditions(time);

    Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
    for (int i = 0; i < inarray.size(); ++i)
    {
        m_fields[i]->FwdTrans(inarray[i], coeffs);
        m_fields[i]->BwdTrans(coeffs, outarray[i]);
    }
}

void PEMFCSpecies::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    const int nvariables = inarray.size();
    const int nq         = m_fields[0]->GetNpoints();

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;

    StdRegions::VarCoeffMap varCoeffMap = StdRegions::NullVarCoeffMap;
    if (m_hasVarDiff)
    {
        // div( f_D(x) grad w ) via diagonal diffusion varcoeffs.
        varCoeffMap[StdRegions::eVarCoeffD00] = m_diffFactor;
        varCoeffMap[StdRegions::eVarCoeffD11] = m_diffFactor;
    }

    // Helmholtz solve for each species:
    //   (div(f_D grad) - lambda) out = -1/(lambda*eps) * in
    for (int i = 0; i < nvariables; ++i)
    {
        Array<OneD, NekDouble> rhs(nq);
        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1, rhs,
                    1);

        m_fields[i]->HelmSolve(rhs, m_fields[i]->UpdateCoeffs(), factors,
                               varCoeffMap);
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
    }
}

void PEMFCSpecies::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    AdvectionSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Kinetics", "Butler-Volmer ORR");
    SolverUtils::AddSummaryItem(s, "Vcell", m_kinetics.m_Vcell);
    SolverUtils::AddSummaryItem(s, "i0", m_kinetics.m_i0);
    SolverUtils::AddSummaryItem(s, "Rlump", m_kinetics.m_Rlump);
    SolverUtils::AddSummaryItem(s, "epsilon", m_epsilon);
}

void PEMFCSpecies::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
    m_fields[0]->FwdTransLocalElmt(m_iloc, coeffs);
    variables.push_back("iloc");
    fieldcoeffs.push_back(coeffs);
}
} // namespace Nektar
