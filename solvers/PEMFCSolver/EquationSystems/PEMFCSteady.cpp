///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCSteady.cpp
//
// See PEMFCSteady.h. Part of the nektar++_fuel_cells project (EPUSP/RCGI).
//
///////////////////////////////////////////////////////////////////////////////

#include <PEMFCSolver/EquationSystems/PEMFCSteady.h>

// Portability: official v5.0.0 has ContField2D; the unified ContField
// appeared in later trees. PEMFC_CONTFIELD2D is set by the CMake probe.
#ifdef PEMFC_CONTFIELD2D
#include <MultiRegions/ContField2D.h>
typedef Nektar::MultiRegions::ContField2D PemfcContField;
typedef Nektar::MultiRegions::ContField2DSharedPtr PemfcContFieldSharedPtr;
#else
#include <MultiRegions/ContField.h>
typedef Nektar::MultiRegions::ContField PemfcContField;
typedef Nektar::MultiRegions::ContFieldSharedPtr PemfcContFieldSharedPtr;
#endif

namespace Nektar
{
std::string PEMFCSteady::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "PEMFCSteady", PEMFCSteady::create);

PEMFCSteady::PEMFCSteady(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph)
{
}

PEMFCSteady::~PEMFCSteady()
{
}

void PEMFCSteady::v_InitObject(bool DeclareFields)
{
    EquationSystem::v_InitObject(DeclareFields);

    ASSERTL0(m_fields.size() == 1 && m_session->GetVariable(0) == "wO2",
             "PEMFCSteady expects a single variable named wO2.");

    // Transport and closure parameters (SI).
    m_session->LoadParameter("epsilon", m_epsilon, 2.2e-5);
    m_session->LoadParameter("rho", m_rho, 1.0);
    m_session->LoadParameter("Htot", m_Htot, 1.38e-3);
    m_session->LoadParameter("MO2", m_MO2, 0.032);
    int ne;
    m_session->LoadParameter("ne", ne, 4);
    m_ne = ne;
    m_session->LoadParameter("winlet", m_winlet, 0.228);

    m_session->LoadParameter("i0", m_kinetics.m_i0, 0.17);
    m_session->LoadParameter("alphaC", m_kinetics.m_alphaC, 0.5);
    m_session->LoadParameter("Tcell", m_kinetics.m_T, 343.15);
    m_session->LoadParameter("Voc", m_kinetics.m_Voc, 1.22);
    m_session->LoadParameter("Vcell", m_kinetics.m_Vcell, 0.4);
    m_session->LoadParameter("Rlump", m_kinetics.m_Rlump, 7.0e-5);
    m_session->LoadParameter("CO2ref", m_kinetics.m_CO2ref, 30.0);

    // Stefan/dilution correction: net gas-phase mass source (water
    // production + electro-osmotic drag + membrane flux) dilutes the
    // oxygen mass fraction beyond consumption. Modeled as an extra
    // linear sink k_dil [1/s] = (net mass source)/(rho*Htot), taken
    // proportional to i_loc via the same secant scaling. Ground truth
    // from the thesis COMSOL model (R_wH2O = production + F_eo + J_H2O).
    m_session->LoadParameter("DilutionMassRatio", m_dilRatio, 0.0);

    int iters;
    m_session->LoadParameter("PicardIters", iters, 8);
    m_picardIters = iters;
    m_session->LoadParameter("PicardTol", m_picardTol, 1.0e-6);
    m_session->LoadParameter("PicardRelax", m_picardRelax, 1.0);

    const int nq = GetNpoints();

    // Physical advection velocity (unscaled).
    m_velocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    std::vector<std::string> vel;
    vel.push_back("Vx");
    vel.push_back("Vy");
    vel.push_back("Vz");
    vel.resize(m_spacedim);
    GetFunction("AdvectionVelocity")->Evaluate(vel, m_velocity);

    // Variable diffusion coefficient eps * f_D(x).
    m_diffCoeff = Array<OneD, NekDouble>(nq, m_epsilon);
    m_gamma     = Array<OneD, NekDouble>(nq, 1.0);
    if (m_session->DefinesSolverInfo("DesignRaster"))
    {
        // Design-field mask gamma(r,theta) from a polar ASCII raster
        // (case 05: Sc8 topology). f_D follows the reference SIMP
        // interpolation f_D = fpr + (fpc - fpr) * gamma^simp.
        NekDouble rin, rout, thdeg, fpr, fpc, simp;
        m_session->LoadParameter("RasterRin", rin, 3.8197e-4);
        m_session->LoadParameter("RasterRout", rout, 1.2621e-2);
        m_session->LoadParameter("RasterThetaDeg", thdeg, 60.0);
        m_session->LoadParameter("FDpr", fpr, 0.19703);
        m_session->LoadParameter("FDpc", fpc, 0.92167);
        m_session->LoadParameter("FDsimp", simp, 4.0);
        PEMFC::PolarRaster raster;
        ASSERTL0(raster.Load(m_session->GetSolverInfo("DesignRaster"), rin,
                             rout, 0.0, thdeg * M_PI / 180.0),
                 "PEMFCSteady: could not load DesignRaster file.");
        Array<OneD, NekDouble> xq(nq), yq(nq), zq(nq);
        m_fields[0]->GetCoords(xq, yq, zq);
        for (int q = 0; q < nq; ++q)
        {
            const double g = raster.Eval(xq[q], yq[q]);
            m_gamma[q]     = g;
            m_diffCoeff[q] =
                m_epsilon * (fpr + (fpc - fpr) * std::pow(g, simp));
        }
    }
    else if (m_session->DefinesFunction("DiffusivityFactor"))
    {
        GetFunction("DiffusivityFactor")
            ->Evaluate(m_session->GetVariable(0), m_diffCoeff);
        Vmath::Smul(nq, m_epsilon, m_diffCoeff, 1, m_diffCoeff, 1);
        ASSERTL0(Vmath::Vmin(nq, m_diffCoeff, 1) > 0.0,
                 "DiffusivityFactor must be positive.");
    }

    m_iloc = Array<OneD, NekDouble>(nq, 0.0);
}

void PEMFCSteady::v_DoInitialise([[maybe_unused]] bool dumpInitialConditions)
{
    // Initial guess: inlet composition everywhere (or session-provided
    // InitialConditions if defined).
    if (m_session->DefinesFunction("InitialConditions"))
    {
        GetFunction("InitialConditions")
            ->Evaluate(m_session->GetVariables(), m_fields);
    }
    else
    {
        Vmath::Fill(GetNpoints(), m_winlet, m_fields[0]->UpdatePhys(), 1);
    }
}

void PEMFCSteady::v_DoSolve()
{
    const int nq         = GetNpoints();
    const double F       = 96485.3;
    const double volcoef = m_MO2 / (m_ne * F * m_rho * m_Htot);

    PemfcContFieldSharedPtr cf =
        std::dynamic_pointer_cast<PemfcContField>(m_fields[0]);
    ASSERTL0(cf, "PEMFCSteady requires a continuous Galerkin field.");

    Array<OneD, NekDouble> w = m_fields[0]->UpdatePhys();
    Array<OneD, NekDouble> kfld(nq, 0.0);
    Array<OneD, NekDouble> mass(nq, 1.0);
    Array<OneD, NekDouble> rhs(nq, 0.0);
    Array<OneD, NekDouble> wnew(nq, 0.0);

    // Reference reaction scale for the eFactorLambda / eVarCoeffMass
    // split: lambda * m(x) = -k(x) in the operator convention
    // (operator adds -lambda*m*u; the weak form then contributes
    // +lambda*m*u to the left-hand side, so a SINK k > 0 requires
    // lambda = -k_ref, m(x) = k(x)/k_ref; verified against the
    // ADRSolver SteadyAdvectionDiffusionReaction empirical convention).
    for (int it = 0; it < m_picardIters; ++it)
    {
        // Secant linearization coefficient from the current iterate.
        double kref = 0.0;
        for (int q = 0; q < nq; ++q)
        {
            const double wq  = std::max(w[q], 1.0e-6);
            const double CO2 = m_rho * wq / m_MO2;
            const double il  = m_kinetics.Solve(CO2);
            m_iloc[q]        = il;
            // consumption + dilution: S = -(volcoef*il/w)*w
            //                             -(m_dilRatio*volcoef*il)*w
            kfld[q] = volcoef * il * (1.0 / wq + m_dilRatio);
            kref += kfld[q];
        }
        kref /= nq;

        for (int q = 0; q < nq; ++q)
        {
            mass[q] = kfld[q] / kref;
        }

        StdRegions::VarCoeffMap vc;
        vc[StdRegions::eVarCoeffVelX] = m_velocity[0];
        vc[StdRegions::eVarCoeffVelY] = m_velocity[1];
        vc[StdRegions::eVarCoeffD00]  = m_diffCoeff;
        vc[StdRegions::eVarCoeffD11]  = m_diffCoeff;
        vc[StdRegions::eVarCoeffMass] = mass;

        Vmath::Zero(nq, rhs, 1);

        StdRegions::ConstFactorMap adrFac;
        adrFac[StdRegions::eFactorLambda] = -kref;
        cf->LinearAdvectionDiffusionReactionSolve(
            rhs, m_fields[0]->UpdateCoeffs(), adrFac, vc);

        m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), wnew);

        // Relaxed update and convergence check.
        double dmax = 0.0;
        for (int q = 0; q < nq; ++q)
        {
            const double upd =
                m_picardRelax * wnew[q] + (1.0 - m_picardRelax) * w[q];
            dmax = std::max(dmax, std::fabs(upd - w[q]));
            w[q] = upd;
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            std::cout << "Picard iter " << it << ": kref = " << kref
                      << ", max dw = " << dmax << std::endl;
        }
        if (dmax < m_picardTol)
        {
            break;
        }
    }

    // Final consistent current density.
    for (int q = 0; q < nq; ++q)
    {
        const double CO2 = m_rho * std::max(w[q], 0.0) / m_MO2;
        m_iloc[q]        = m_kinetics.Solve(CO2);
    }
    m_fields[0]->FwdTrans(w, m_fields[0]->UpdateCoeffs());
    m_fields[0]->SetPhysState(true);
}

void PEMFCSteady::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    SolverUtils::AddSummaryItem(s, "EquationType", "PEMFCSteady (Picard)");
    SolverUtils::AddSummaryItem(s, "Kinetics", "Butler-Volmer ORR");
    SolverUtils::AddSummaryItem(s, "Vcell", m_kinetics.m_Vcell);
    SolverUtils::AddSummaryItem(s, "epsilon", m_epsilon);
    SolverUtils::AddSummaryItem(s, "PicardIters", m_picardIters);
}

void PEMFCSteady::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
    m_fields[0]->FwdTransLocalElmt(m_iloc, coeffs);
    variables.push_back("iloc");
    fieldcoeffs.push_back(coeffs);

    if (m_session->DefinesSolverInfo("DesignRaster"))
    {
        Array<OneD, NekDouble> gcoeffs(m_fields[0]->GetNcoeffs());
        m_fields[0]->FwdTransLocalElmt(m_gamma, gcoeffs);
        variables.push_back("gamma");
        fieldcoeffs.push_back(gcoeffs);
    }
}
} // namespace Nektar
