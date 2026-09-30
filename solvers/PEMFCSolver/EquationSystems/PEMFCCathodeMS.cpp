///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCCathodeMS.cpp
//
// See PEMFCCathodeMS.h. Part of the nektar++_fuel_cells project
// (EPUSP / RCGI).
//
///////////////////////////////////////////////////////////////////////////////

#include <PEMFCSolver/EquationSystems/PEMFCCathodeMS.h>

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
std::string PEMFCCathodeMS::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "PEMFCCathodeMS", PEMFCCathodeMS::create);

PEMFCCathodeMS::PEMFCCathodeMS(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph)
{
}

PEMFCCathodeMS::~PEMFCCathodeMS()
{
}

void PEMFCCathodeMS::v_InitObject(bool DeclareFields)
{
    EquationSystem::v_InitObject(DeclareFields);

    ASSERTL0(m_fields.size() == 2 && m_session->GetVariable(0) == "wO2" &&
                 m_session->GetVariable(1) == "wH2O",
             "PEMFCCathodeMS expects variables wO2, wH2O.");

    // Transport parameters: identical calibration to PEMFCSteady.
    m_session->LoadParameter("epsilon", m_epsilon, 2.67e-5);
    m_session->LoadParameter("rho", m_rho, 1.0);
    m_session->LoadParameter("Htot", m_Htot, 1.38e-3);
    m_session->LoadParameter("MO2", m_MO2, 0.032);
    m_session->LoadParameter("MH2O", m_MH2O, 0.018);
    m_session->LoadParameter("MN2", m_MN2, 0.028);
    int ne;
    m_session->LoadParameter("ne", ne, 4);
    m_ne = ne;
    m_session->LoadParameter("WO2In", m_wO2in, 0.228);
    m_session->LoadParameter("WH2OIn", m_wH2Oin, 0.0236);

    // Kinetics: same names and defaults as PEMFCSteady.
    m_session->LoadParameter("i0", m_kinetics.m_i0, 0.17);
    m_session->LoadParameter("alphaC", m_kinetics.m_alphaC, 0.5);
    m_session->LoadParameter("Tcell", m_kinetics.m_T, 343.15);
    m_session->LoadParameter("Voc", m_kinetics.m_Voc, 1.22);
    m_session->LoadParameter("Vcell", m_kinetics.m_Vcell, 0.4);
    m_session->LoadParameter("Rlump", m_kinetics.m_Rlump, 7.0e-5);
    m_session->LoadParameter("CO2ref", m_kinetics.m_CO2ref, 30.0);

    // Multi-species closure parameters.
    m_session->LoadParameter("XO2In", m_xO2in, 0.21);
    m_session->LoadParameter("XH2OIn", m_xH2Oin, 0.038);
    m_session->LoadParameter("DO2H2O", m_DO2H2O, 1.78e-5);
    m_session->LoadParameter("DO2N2", m_DO2N2, 1.40e-5);
    m_session->LoadParameter("DN2H2O", m_DN2H2O, 1.75e-5);
    m_session->LoadParameter("WaterSourceRatio", m_rw, 3.375);
    m_session->LoadParameter("StefanCorrection", m_stefan, 1.0);

    int iters;
    m_session->LoadParameter("PicardIters", iters, 12);
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

    // f_D(x) mask, applied to both species.
    m_fD = Array<OneD, NekDouble>(nq, 1.0);
    if (m_session->DefinesFunction("DiffusivityFactor"))
    {
        GetFunction("DiffusivityFactor")
            ->Evaluate(m_session->GetVariable(0), m_fD);
        ASSERTL0(Vmath::Vmin(nq, m_fD, 1) > 0.0,
                 "DiffusivityFactor must be positive.");
    }

    m_iloc   = Array<OneD, NekDouble>(nq, 0.0);
    m_DO2fac = Array<OneD, NekDouble>(nq, 1.0);
}

void PEMFCCathodeMS::v_DoInitialise([[maybe_unused]] bool dumpInitialConditions)
{
    if (m_session->DefinesFunction("InitialConditions"))
    {
        GetFunction("InitialConditions")
            ->Evaluate(m_session->GetVariables(), m_fields);
        return;
    }
    Vmath::Fill(GetNpoints(), m_wO2in, m_fields[0]->UpdatePhys(), 1);
    Vmath::Fill(GetNpoints(), m_wH2Oin, m_fields[1]->UpdatePhys(), 1);
    for (int f = 0; f < 2; ++f)
    {
        m_fields[f]->FwdTrans(m_fields[f]->GetPhys(),
                              m_fields[f]->UpdateCoeffs());
    }
}

void PEMFCCathodeMS::SolveField(int f,
                                const Array<OneD, Array<OneD, NekDouble>> &vel,
                                const Array<OneD, NekDouble> &diff,
                                const Array<OneD, NekDouble> &kfld,
                                const Array<OneD, NekDouble> &rhs)
{
    const int nq = GetNpoints();
    PemfcContFieldSharedPtr cf =
        std::dynamic_pointer_cast<PemfcContField>(m_fields[f]);
    ASSERTL0(cf, "PEMFCCathodeMS requires continuous Galerkin fields.");

    double kref = 0.0;
    for (int q = 0; q < nq; ++q)
    {
        kref += kfld[q];
    }
    kref /= nq;
    kref = std::max(kref, 1.0e-12);

    Array<OneD, NekDouble> mass(nq);
    for (int q = 0; q < nq; ++q)
    {
        mass[q] = kfld[q] / kref;
    }

    StdRegions::VarCoeffMap vc;
    vc[StdRegions::eVarCoeffVelX] = vel[0];
    vc[StdRegions::eVarCoeffVelY] = vel[1];
    vc[StdRegions::eVarCoeffD00]  = diff;
    vc[StdRegions::eVarCoeffD11]  = diff;
    vc[StdRegions::eVarCoeffMass] = mass;

    // Positive physical source S enters as inarray = -S.
    Array<OneD, NekDouble> rhsNeg(nq);
    Vmath::Smul(nq, -1.0, rhs, 1, rhsNeg, 1);
    StdRegions::ConstFactorMap adrFac;
    adrFac[StdRegions::eFactorLambda] = -kref;
    cf->LinearAdvectionDiffusionReactionSolve(
        rhsNeg, m_fields[f]->UpdateCoeffs(), adrFac, vc);
    m_fields[f]->BwdTrans(m_fields[f]->GetCoeffs(), m_fields[f]->UpdatePhys());
}

void PEMFCCathodeMS::v_DoSolve()
{
    const int nq         = GetNpoints();
    const double F       = 96485.3;
    const double volcoef = m_MO2 / (m_ne * F * m_rho * m_Htot);

    // Inlet-composition D_O2^mix for the epsilon normalization.
    const double xN2in = 1.0 - m_xO2in - m_xH2Oin;
    const double DO2in =
        (1.0 - m_xO2in) / (m_xH2Oin / m_DO2H2O + xN2in / m_DO2N2);

    Array<OneD, NekDouble> kfld(nq), rhs(nq);
    Array<OneD, NekDouble> DO2(nq), DH2O(nq), DN2(nq);
    Array<OneD, NekDouble> diffO2(nq), diffH2O(nq);
    Array<OneD, Array<OneD, NekDouble>> uEff(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> gO2(m_spacedim), gW(m_spacedim);
    for (int d = 0; d < m_spacedim; ++d)
    {
        uEff[d] = Array<OneD, NekDouble>(nq, 0.0);
        gO2[d]  = Array<OneD, NekDouble>(nq, 0.0);
        gW[d]   = Array<OneD, NekDouble>(nq, 0.0);
    }
    Array<OneD, NekDouble> wo_old(nq), ww_old(nq);

    for (int it = 0; it < m_picardIters; ++it)
    {
        Array<OneD, NekDouble> wo = m_fields[0]->UpdatePhys();
        Array<OneD, NekDouble> ww = m_fields[1]->UpdatePhys();
        Vmath::Vcopy(nq, wo, 1, wo_old, 1);
        Vmath::Vcopy(nq, ww, 1, ww_old, 1);

        // Local composition, mixture diffusivities, kinetics.
        for (int q = 0; q < nq; ++q)
        {
            const double o    = std::max(wo[q], 0.0);
            const double w    = std::max(ww[q], 0.0);
            const double n    = std::max(1.0 - o - w, 1.0e-8);
            const double xden = o / m_MO2 + w / m_MH2O + n / m_MN2;
            const double xO2  = (o / m_MO2) / xden;
            const double xH2O = (w / m_MH2O) / xden;
            const double xN2  = (n / m_MN2) / xden;

            DO2[q] =
                (1.0 - xO2) / std::max(xH2O / m_DO2H2O + xN2 / m_DO2N2, 1.0e-8);
            DH2O[q] = (1.0 - xH2O) /
                      std::max(xO2 / m_DO2H2O + xN2 / m_DN2H2O, 1.0e-8);
            DN2[q] =
                (1.0 - xN2) / std::max(xO2 / m_DO2N2 + xH2O / m_DN2H2O, 1.0e-8);

            const double scale = m_epsilon / DO2in;
            m_DO2fac[q]        = DO2[q] / DO2in;
            diffO2[q]          = scale * DO2[q] * m_fD[q];
            diffH2O[q]         = scale * DH2O[q] * m_fD[q];

            const double CO2 = m_rho * std::max(wo[q], 1.0e-6) / m_MO2;
            m_iloc[q]        = m_kinetics.Solve(CO2);
        }

        // Stefan/correction velocity from the current gradients:
        // u_St = fD*scale*[(DO2-DN2) grad wO2 + (DH2O-DN2) grad wH2O]
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vcopy(nq, m_velocity[d], 1, uEff[d], 1);
        }
        if (m_stefan > 0.5)
        {
            m_fields[0]->PhysDeriv(wo, gO2[0], gO2[1]);
            m_fields[1]->PhysDeriv(ww, gW[0], gW[1]);
            const double scale = m_epsilon / DO2in;
            for (int d = 0; d < m_spacedim; ++d)
            {
                for (int q = 0; q < nq; ++q)
                {
                    uEff[d][q] += m_stefan * m_fD[q] * scale *
                                  ((DO2[q] - DN2[q]) * gO2[d][q] +
                                   (DH2O[q] - DN2[q]) * gW[d][q]);
                }
            }
        }

        // O2: consumption (secant) + dilution (linear).
        for (int q = 0; q < nq; ++q)
        {
            const double wq = std::max(wo_old[q], 1.0e-6);
            kfld[q]         = volcoef * m_iloc[q] * (1.0 / wq + (m_rw - 1.0));
            rhs[q]          = 0.0;
        }
        SolveField(0, uEff, diffO2, kfld, rhs);

        // H2O: dilution (linear sink) + production (explicit source).
        for (int q = 0; q < nq; ++q)
        {
            kfld[q] = (m_rw - 1.0) * volcoef * m_iloc[q];
            rhs[q]  = m_rw * volcoef * m_iloc[q];
        }
        SolveField(1, uEff, diffH2O, kfld, rhs);

        // Relaxation + convergence on both fields.
        double dmax = 0.0;
        for (int q = 0; q < nq; ++q)
        {
            double upd =
                m_picardRelax * wo[q] + (1.0 - m_picardRelax) * wo_old[q];
            dmax  = std::max(dmax, std::fabs(upd - wo_old[q]));
            wo[q] = upd;
            upd   = m_picardRelax * ww[q] + (1.0 - m_picardRelax) * ww_old[q];
            dmax  = std::max(dmax, std::fabs(upd - ww_old[q]));
            ww[q] = upd;
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            std::cout << "Picard iter " << it << ": max dw = " << dmax
                      << ", min wO2 = " << Vmath::Vmin(nq, wo, 1)
                      << ", max iloc = " << Vmath::Vmax(nq, m_iloc, 1)
                      << std::endl;
        }

        // Varcoeff maps change every iteration: clear the cached
        // factorizations (Model B memory lesson).
        m_fields[0]->ClearGlobalLinSysManager();

        if (dmax < m_picardTol)
        {
            break;
        }
    }

    // Final consistent current density.
    for (int q = 0; q < nq; ++q)
    {
        const double CO2 =
            m_rho * std::max(m_fields[0]->GetPhys()[q], 0.0) / m_MO2;
        m_iloc[q] = m_kinetics.Solve(CO2);
    }
    for (int f = 0; f < 2; ++f)
    {
        m_fields[f]->FwdTrans(m_fields[f]->GetPhys(),
                              m_fields[f]->UpdateCoeffs());
        m_fields[f]->SetPhysState(true);
    }
}

void PEMFCCathodeMS::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    SolverUtils::AddSummaryItem(s, "EquationType", "PEMFCCathodeMS (gate a)");
    SolverUtils::AddSummaryItem(s, "Kinetics", "Butler-Volmer ORR");
    SolverUtils::AddSummaryItem(s, "Vcell", m_kinetics.m_Vcell);
    SolverUtils::AddSummaryItem(s, "WaterSourceRatio", m_rw);
    SolverUtils::AddSummaryItem(s, "StefanCorrection", m_stefan);
    SolverUtils::AddSummaryItem(s, "PicardIters", m_picardIters);
}

void PEMFCCathodeMS::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    const char *names[2]              = {"iloc", "DO2fac"};
    Array<OneD, NekDouble> *arrays[2] = {&m_iloc, &m_DO2fac};
    for (int i = 0; i < 2; ++i)
    {
        Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
        m_fields[0]->FwdTransLocalElmt(*arrays[i], coeffs);
        variables.push_back(names[i]);
        fieldcoeffs.push_back(coeffs);
    }
}
} // namespace Nektar
