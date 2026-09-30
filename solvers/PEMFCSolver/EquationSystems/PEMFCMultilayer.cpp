///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCMultilayer.cpp
//
// See PEMFCMultilayer.h. Part of the nektar++_fuel_cells project
// (EPUSP / RCGI / BUGA-LAB).
//
///////////////////////////////////////////////////////////////////////////////

#include <PEMFCSolver/EquationSystems/PEMFCMultilayer.h>
#include <algorithm>

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
std::string PEMFCMultilayer::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "PEMFCMultilayer", PEMFCMultilayer::create);

namespace
{
// Layer of each field: 0=cc, 1=cg, 2=ac, 3=ag.
const int kLayerOf[PEMFCMultilayer::eNumFields] = {0, 0, 1, 1, 2, 2, 3, 3};
const char *kLayerFn[4] = {"Velocity_cc", "Velocity_cg", "Velocity_ac",
                           "Velocity_ag"};
} // namespace

PEMFCMultilayer::PEMFCMultilayer(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph)
{
}

PEMFCMultilayer::~PEMFCMultilayer()
{
}

void PEMFCMultilayer::v_InitObject(bool DeclareFields)
{
    EquationSystem::v_InitObject(DeclareFields);

    ASSERTL0((int)m_fields.size() == (int)eNumFields,
             "PEMFCMultilayer expects 8 variables: wO2_cc, wH2O_cc, "
             "wO2_cg, wH2O_cg, wH2_ac, wH2O_ac, wH2_ag, wH2O_ag.");
    const char *names[eNumFields] = {"wO2_cc", "wH2O_cc", "wO2_cg", "wH2O_cg",
                                     "wH2_ac", "wH2O_ac", "wH2_ag", "wH2O_ag"};
    for (int f = 0; f < eNumFields; ++f)
    {
        ASSERTL0(m_session->GetVariable(f) == names[f],
                 "PEMFCMultilayer: variable order/name mismatch.");
    }

    // Geometry and gas properties (ground-truth defaults).
    m_session->LoadParameter("Hch", m_Hch, 1.0e-3);
    m_session->LoadParameter("Hgdl", m_Hgdl, 380.0e-6);
    m_session->LoadParameter("PressIn", m_pAbs, 2.01e5);
    m_session->LoadParameter("Tcell", m_T, 353.15);
    // Mixture densities (ideal gas at inlet composition; constant
    // in B1). Cathode humid air at 2 atm, anode humid H2.
    m_session->LoadParameter("RhoCath", m_rhoC, 1.95);
    m_session->LoadParameter("RhoAn", m_rhoA, 0.40);
    m_session->LoadParameter("DO2H2O", m_DO2H2O, 1.78e-5);
    m_session->LoadParameter("DO2N2", m_DO2N2, 1.40e-5);
    m_session->LoadParameter("DN2H2O", m_DN2H2O, 1.75e-5);
    m_session->LoadParameter("DH2H2O", m_DH2H2O, 1.78e-5);
    m_session->LoadParameter("EpsGdl", m_epsGdl, 0.4);

    // Kinetics parameter overrides.
    m_session->LoadParameter("Ecell", m_kin.m_Ecell, 0.8);
    m_session->LoadParameter("Tcell", m_kin.m_T, 353.15);
    m_session->LoadParameter("I0c", m_kin.m_i0c, 1.0e4);
    m_session->LoadParameter("I0a", m_kin.m_i0a, 2.0e8);
    m_session->LoadParameter("AlphaBV", m_kin.m_alpha, 0.5);
    m_session->LoadParameter("EeqRefC", m_kin.m_EeqRefC, 1.19);
    m_session->LoadParameter("EeqRefA", m_kin.m_EeqRefA, 0.0);
    m_kin.m_Hch  = m_Hch;
    m_kin.m_Hgdl = m_Hgdl;

    int iters;
    m_session->LoadParameter("PicardIters", iters, 20);
    m_picardIters = iters;
    m_session->LoadParameter("PicardTol", m_picardTol, 1.0e-7);
    m_session->LoadParameter("PicardRelax", m_picardRelax, 0.7);

    const int nq = GetNpoints();

    // Frozen per-layer velocities (default: quiescent).
    m_vel = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(4);
    for (int l = 0; l < 4; ++l)
    {
        m_vel[l] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_vel[l][d] = Array<OneD, NekDouble>(nq, 0.0);
        }
        if (m_session->DefinesFunction(kLayerFn[l]))
        {
            std::vector<std::string> v;
            v.push_back("Vx");
            v.push_back("Vy");
            v.resize(m_spacedim);
            GetFunction(kLayerFn[l])->Evaluate(v, m_vel[l]);
        }
    }

    // Design masks (default fully open).
    m_designC = Array<OneD, NekDouble>(nq, 1.0);
    m_designA = Array<OneD, NekDouble>(nq, 1.0);
    if (m_session->DefinesFunction("DesignCathode"))
    {
        GetFunction("DesignCathode")
            ->Evaluate(m_session->GetVariable(0), m_designC);
    }
    if (m_session->DefinesFunction("DesignAnode"))
    {
        GetFunction("DesignAnode")
            ->Evaluate(m_session->GetVariable(4), m_designA);
    }

    // v12.3 clamp: file-based design masks over/undershoot [0, 1]
    // after projection; clamp pointwise.
    for (int iq = 0; iq < nq; ++iq)
    {
        m_designC[iq] = std::min(1.0, std::max(0.0, m_designC[iq]));
        m_designA[iq] = std::min(1.0, std::max(0.0, m_designA[iq]));
    }

    m_iloc  = Array<OneD, NekDouble>(nq, 0.0);
    m_etaC  = Array<OneD, NekDouble>(nq, -0.3);
    m_etaAn = Array<OneD, NekDouble>(nq, 0.0);
}

void PEMFCMultilayer::v_DoInitialise(
    [[maybe_unused]] bool dumpInitialConditions)
{
    // Inlet compositions everywhere as the initial guess (or
    // session-provided InitialConditions).
    if (m_session->DefinesFunction("InitialConditions"))
    {
        GetFunction("InitialConditions")
            ->Evaluate(m_session->GetVariables(), m_fields);
        return;
    }
    NekDouble wO2in, wH2Ocin, wH2in, wH2Oain;
    m_session->LoadParameter("WO2In", wO2in, 0.23);
    m_session->LoadParameter("WH2OcIn", wH2Ocin, 0.0236);
    m_session->LoadParameter("WH2In", wH2in, 0.67);
    m_session->LoadParameter("WH2OaIn", wH2Oain, 0.33);
    const NekDouble init[eNumFields] = {wO2in, wH2Ocin, wO2in, wH2Ocin,
                                        wH2in, wH2Oain, wH2in, wH2Oain};
    for (int f = 0; f < eNumFields; ++f)
    {
        Vmath::Fill(GetNpoints(), init[f], m_fields[f]->UpdatePhys(), 1);
        m_fields[f]->FwdTrans(m_fields[f]->GetPhys(),
                              m_fields[f]->UpdateCoeffs());
    }
}

void PEMFCMultilayer::UpdateClosures()
{
    const int nq                               = GetNpoints();
    const Array<OneD, const NekDouble> &wO2g   = m_fields[eO2cg]->GetPhys();
    const Array<OneD, const NekDouble> &wH2Og  = m_fields[eH2Ocg]->GetPhys();
    const Array<OneD, const NekDouble> &wH2g   = m_fields[eH2ag]->GetPhys();
    const Array<OneD, const NekDouble> &wH2Oag = m_fields[eH2Oag]->GetPhys();

    for (int q = 0; q < nq; ++q)
    {
        // Cathode GDL mole fractions (N2 by closure).
        const double wo    = std::max(wO2g[q], 0.0);
        const double ww    = std::max(wH2Og[q], 0.0);
        const double wn    = std::max(1.0 - wo - ww, 1.0e-8);
        double xden        = wo / kMO2 + ww / kMH2O + wn / kMN2;
        const double xO2   = (wo / kMO2) / xden;
        const double xH2Oc = (ww / kMH2O) / xden;
        // Anode GDL mole fractions.
        const double wh  = std::max(wH2g[q], 0.0);
        const double wwa = std::max(wH2Oag[q], 0.0);
        const double wna = std::max(1.0 - wh - wwa, 1.0e-8);
        xden             = wh / kMH2 + wwa / kMH2O + wna / kMN2;
        const double xH2 = (wh / kMH2) / xden;

        double etaC = m_etaC[q], etaAn = m_etaAn[q];
        const double il = m_kin.Solve(xO2 * m_pAbs, xH2Oc * m_pAbs,
                                      xH2 * m_pAbs, etaC, etaAn);
        // Kinetics live where the electrode exists; B1 applies them
        // over the whole domain (the GDL layer spans it).
        m_iloc[q]  = il;
        m_etaC[q]  = etaC;
        m_etaAn[q] = etaAn;
    }
}

void PEMFCMultilayer::SolveField(int f, const Array<OneD, NekDouble> &kfld,
                                 const Array<OneD, NekDouble> &rhs)
{
    const int nq = GetNpoints();
    PemfcContFieldSharedPtr cf =
        std::dynamic_pointer_cast<PemfcContField>(m_fields[f]);
    ASSERTL0(cf, "PEMFCMultilayer requires continuous Galerkin fields.");

    // Reference reaction scale (PEMFCSteady convention: sink k > 0
    // enters as lambda = -kref with mass = k/kref).
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

    // Layer diffusivity: mixture-reduced, Bruggeman in the GDL.
    // B1 uses the inlet-composition constant per layer.
    const int layer   = kLayerOf[f];
    const bool gdl    = (layer == 1 || layer == 3);
    const double brug = gdl ? std::pow(m_epsGdl, 1.5) : 1.0;
    double D;
    if (layer <= 1) // cathode
    {
        // Dm at inlet composition (xO2 0.21, xH2O 0.038, xN2 rest).
        D = (f == eO2cc || f == eO2cg)
                ? (1.0 - 0.21) / (0.038 / m_DO2H2O + 0.752 / m_DO2N2)
                : (1.0 - 0.038) / (0.21 / m_DO2H2O + 0.752 / m_DN2H2O);
    }
    else // anode
    {
        D = m_DH2H2O;
    }
    D *= brug;

    Array<OneD, NekDouble> diff(nq, D);

    StdRegions::VarCoeffMap vc;
    vc[StdRegions::eVarCoeffVelX] = m_vel[layer][0];
    vc[StdRegions::eVarCoeffVelY] = m_vel[layer][1];
    vc[StdRegions::eVarCoeffD00]  = diff;
    vc[StdRegions::eVarCoeffD11]  = diff;
    vc[StdRegions::eVarCoeffMass] = mass;

    // Forcing sign: the solve assembles -IProduct(inarray), so the
    // equation satisfied is  U.grad w - div(D grad w) + k w = -inarray;
    // a POSITIVE physical source S therefore enters as inarray = -S.
    Array<OneD, NekDouble> rhsNeg(nq);
    Vmath::Smul(nq, -1.0, rhs, 1, rhsNeg, 1);
    StdRegions::ConstFactorMap adrFac;
    adrFac[StdRegions::eFactorLambda] = -kref;
    cf->LinearAdvectionDiffusionReactionSolve(
        rhsNeg, m_fields[f]->UpdateCoeffs(), adrFac, vc);
    m_fields[f]->BwdTrans(m_fields[f]->GetCoeffs(), m_fields[f]->UpdatePhys());
}

void PEMFCMultilayer::v_DoSolve()
{
    const int nq   = GetNpoints();
    const double F = 96485.3;

    // Inter-layer diffusive-exchange coefficients [1/s], from the
    // ground-truth Ndiff closure with layer-thickness averaging:
    //   flux = -Dm_eff * 2 * (rho_ch w_ch - rho_gdl w_gdl)/(Hch+Hgdl)
    //   volumetric in the channel: /Hch; in the GDL: /Hgdl.
    // With B1 constant layer densities the exchange is linear.
    const double Hsum = m_Hch + m_Hgdl;

    // Mixture-reduced effective diffusivities (inlet composition),
    // thickness-averaged ch/gdl as in the ground truth.
    const double brug    = std::pow(m_epsGdl, 1.5);
    const double DmO2ch  = (1.0 - 0.21) / (0.038 / m_DO2H2O + 0.752 / m_DO2N2);
    const double DmH2Och = (1.0 - 0.038) / (0.21 / m_DO2H2O + 0.752 / m_DN2H2O);
    const double DmO2eff = (m_Hch * DmO2ch + m_Hgdl * DmO2ch * brug) / Hsum;
    const double DmH2Oeff = (m_Hch * DmH2Och + m_Hgdl * DmH2Och * brug) / Hsum;
    const double DmH2eff = (m_Hch * m_DH2H2O + m_Hgdl * m_DH2H2O * brug) / Hsum;

    // C-coefficients: mass-flux prefactor 2*Dm_eff/Hsum [m/s].
    const double CO2  = 2.0 * DmO2eff / Hsum;
    const double CH2O = 2.0 * DmH2Oeff / Hsum;
    const double CH2  = 2.0 * DmH2eff / Hsum;

    Array<OneD, NekDouble> kfld(nq), rhs(nq);
    std::vector<Array<OneD, NekDouble>> wold(eNumFields);
    for (int f = 0; f < eNumFields; ++f)
    {
        wold[f] = Array<OneD, NekDouble>(nq);
    }

    for (int it = 0; it < m_picardIters; ++it)
    {
        for (int f = 0; f < eNumFields; ++f)
        {
            Vmath::Vcopy(nq, m_fields[f]->GetPhys(), 1, wold[f], 1);
        }
        UpdateClosures();

        for (int f = 0; f < eNumFields; ++f)
        {
            const int layer                    = kLayerOf[f];
            const bool cath                    = (layer <= 1);
            const bool gdl                     = (layer == 1 || layer == 3);
            const double H                     = gdl ? m_Hgdl : m_Hch;
            const double rhoL                  = cath ? m_rhoC : m_rhoA;
            const Array<OneD, NekDouble> &mask = cath ? m_designC : m_designA;

            // Exchange partner field and species C-coefficient.
            int fo; // partner index (other layer)
            double C;
            switch (f)
            {
                case eO2cc:
                    fo = eO2cg;
                    C  = CO2;
                    break;
                case eH2Occ:
                    fo = eH2Ocg;
                    C  = CH2O;
                    break;
                case eO2cg:
                    fo = eO2cc;
                    C  = CO2;
                    break;
                case eH2Ocg:
                    fo = eH2Occ;
                    C  = CH2O;
                    break;
                case eH2ac:
                    fo = eH2ag;
                    C  = CH2;
                    break;
                case eH2Oac:
                    fo = eH2Oag;
                    C  = CH2O;
                    break;
                case eH2ag:
                    fo = eH2ac;
                    C  = CH2;
                    break;
                default:
                    fo = eH2Oac;
                    C  = CH2O;
                    break; // eH2Oag
            }
            const Array<OneD, const NekDouble> &wpart = m_fields[fo]->GetPhys();

            for (int q = 0; q < nq; ++q)
            {
                // Diagonal sink and lagged partner source, per unit
                // layer volume, divided by the layer density.
                const double cex = mask[q] * C / H; // [1/s]
                kfld[q]          = cex;
                rhs[q]           = cex * wpart[q]; // rho cancels

                // Reaction fluxes in the GDL layers (per-area iloc
                // spread over the GDL thickness).
                if (gdl)
                {
                    const double ivol = m_iloc[q] / (m_Hgdl * rhoL * F);
                    if (f == eO2cg)
                    {
                        // consumption: secant-implicit
                        const double w = std::max(wold[f][q], 1.0e-6);
                        kfld[q] += (kMO2 / 4.0) * ivol / w;
                    }
                    else if (f == eH2ag)
                    {
                        const double w = std::max(wold[f][q], 1.0e-6);
                        kfld[q] += (kMH2 / 2.0) * ivol / w;
                    }
                    else if (f == eH2Ocg)
                    {
                        // production: explicit source
                        rhs[q] += 2.0 * (kMH2O / 4.0) * ivol;
                    }
                }
            }
            SolveField(f, kfld, rhs);
        }

        // Relaxation + convergence.
        double dmax = 0.0;
        for (int f = 0; f < eNumFields; ++f)
        {
            Array<OneD, NekDouble> w = m_fields[f]->UpdatePhys();
            for (int q = 0; q < nq; ++q)
            {
                const double upd =
                    m_picardRelax * w[q] + (1.0 - m_picardRelax) * wold[f][q];
                dmax = std::max(dmax, std::fabs(upd - wold[f][q]));
                w[q] = upd;
            }
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            std::cout << "Picard iter " << it << ": max dw = " << dmax
                      << ", max iloc = " << Vmath::Vmax(nq, m_iloc, 1)
                      << std::endl;
        }

        // The varcoeff maps change every iteration, so every solve
        // builds a NEW cached global matrix (8 dense factorizations
        // per iteration). Clear the manager each iteration or the
        // run exhausts memory (observed OOM at P6 after 16 iters).
        m_fields[0]->ClearGlobalLinSysManager();

        if (dmax < m_picardTol)
        {
            break;
        }
    }

    UpdateClosures();
    for (int f = 0; f < eNumFields; ++f)
    {
        m_fields[f]->FwdTrans(m_fields[f]->GetPhys(),
                              m_fields[f]->UpdateCoeffs());
        m_fields[f]->SetPhysState(true);
    }
}

void PEMFCMultilayer::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    SolverUtils::AddSummaryItem(s, "EquationType",
                                "PEMFCMultilayer (B1 skeleton)");
    SolverUtils::AddSummaryItem(s, "Ecell", m_kin.m_Ecell);
    SolverUtils::AddSummaryItem(s, "PicardIters", m_picardIters);
}

void PEMFCMultilayer::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    const char *names[3]              = {"iloc", "etaC", "etaAn"};
    Array<OneD, NekDouble> *arrays[3] = {&m_iloc, &m_etaC, &m_etaAn};
    for (int i = 0; i < 3; ++i)
    {
        Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
        m_fields[0]->FwdTransLocalElmt(*arrays[i], coeffs);
        variables.push_back(names[i]);
        fieldcoeffs.push_back(coeffs);
    }
}
} // namespace Nektar
