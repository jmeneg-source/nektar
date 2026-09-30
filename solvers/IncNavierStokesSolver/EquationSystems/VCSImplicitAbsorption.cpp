///////////////////////////////////////////////////////////////////////////////
//
// File: VCSImplicitAbsorption.cpp
//
// See VCSImplicitAbsorption.h. Part of the nektar++_fuel_cells project.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VCSImplicitAbsorption.h>
#include <PEMFCSolver/Kinetics/PolarRaster.h>
#include <SolverUtils/Core/Misc.h>
#include <algorithm>

namespace Nektar
{
std::string VCSImplicitAbsorption::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "VCSImplicitAbsorption", VCSImplicitAbsorption::create);

// 5.10 port: VerifySolverInfo() rejects any SOLVERTYPE value that is
// not a registered enum value, so the addon must register its name.
// It behaves as a VelocityCorrectionScheme for enum-based dispatch.
std::string VCSImplicitAbsorption::solverTypeLookupId =
    LibUtilities::SessionReader::RegisterEnumValue(
        "SolverType", "VCSImplicitAbsorption", eVelocityCorrectionScheme);

VCSImplicitAbsorption::VCSImplicitAbsorption(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph),
      VelocityCorrectionScheme(pSession, pGraph)
{
}

VCSImplicitAbsorption::~VCSImplicitAbsorption()
{
}

void VCSImplicitAbsorption::v_InitObject(bool DeclareField)
{
    VelocityCorrectionScheme::v_InitObject(DeclareField);

    const int nq = m_fields[0]->GetNpoints();
    m_sigma      = Array<OneD, NekDouble>(nq, 0.0);
    m_massVar    = Array<OneD, NekDouble>(nq, 1.0);

    if (m_session->DefinesSolverInfo("DesignRaster"))
    {
        // Design-field mask gamma(r,theta) (case 05, Sc8): Brinkman
        // drag by the reference RAMP interpolation
        //   sigma(gamma) = AlphaMax * (1 - gamma) / (1 + q * gamma),
        // AlphaMax = mu / kappa_pr (5.447e6 1/s), q = 1 (thesis .mph).
        NekDouble rin, rout, thdeg, amax, qq;
        m_session->LoadParameter("RasterRin", rin, 3.8197e-4);
        m_session->LoadParameter("RasterRout", rout, 1.2621e-2);
        m_session->LoadParameter("RasterThetaDeg", thdeg, 60.0);
        m_session->LoadParameter("AlphaMax", amax, 5.447e6);
        m_session->LoadParameter("AlphaQ", qq, 1.0);
        PEMFC::PolarRaster raster;
        ASSERTL0(raster.Load(m_session->GetSolverInfo("DesignRaster"), rin,
                             rout, 0.0, thdeg * M_PI / 180.0),
                 "VCSImplicitAbsorption: could not load DesignRaster.");
        Array<OneD, NekDouble> xq(nq), yq(nq), zq(nq);
        m_fields[0]->GetCoords(xq, yq, zq);
        for (int q = 0; q < nq; ++q)
        {
            // Solver v10: threshold the bilinear raster value to a
            // sharp 0/1 mask before the RAMP (kills the gray band at
            // channel walls; worth 0.2 Pa on the Sc8 pressure drop).
            const double g = raster.Eval(xq[q], yq[q]) > 0.5 ? 1.0 : 0.0;
            m_sigma[q]     = amax * (1.0 - g) / (1.0 + qq * g);
        }
    }
    else
    {
        ASSERTL0(m_session->DefinesFunction("ImplicitAbsorption") ||
                     m_session->DefinesFunction("AbsorptionDesign"),
                 "VCSImplicitAbsorption requires a session function "
                 "'ImplicitAbsorption' (sigma(x) >= 0, in 1/s) or a "
                 "SOLVERINFO DesignRaster entry.");

        // Evaluate sigma(x) once (steady coefficient field), using the
        // first velocity variable's entry of the function.
        if (m_session->DefinesFunction("ImplicitAbsorption"))
            GetFunction("ImplicitAbsorption")
                ->Evaluate(m_session->GetVariable(0), m_sigma);
    }

    // Stage 4 (v12.3): design-file absorption, RAMP closure applied
    // pointwise at the quadrature points (no projection of sigma).
    if (m_session->DefinesFunction("AbsorptionDesign"))
    {
        NekDouble s0, kmax, qr;
        m_session->LoadParameter("AbsSigma0", s0, 0.0);
        m_session->LoadParameter("AbsKmax", kmax, 0.0);
        m_session->LoadParameter("AbsQ", qr, 0.1);
        Array<OneD, NekDouble> dsg(nq, 1.0);
        GetFunction("AbsorptionDesign")
            ->Evaluate(m_session->GetVariable(0), dsg);
        for (int iq = 0; iq < nq; ++iq)
        {
            const NekDouble d = std::min(1.0, std::max(0.0, dsg[iq]));
            m_sigma[iq]       = s0 + kmax * (1.0 - d * (1.0 + qr) / (d + qr));
        }
    }

    // Sanity: absorption must be non-negative here (the sign is NOT the
    // ForcingAbsorption convention; sigma >= 0 is a sink).
    NekDouble smin = Vmath::Vmin(nq, m_sigma, 1);
    ASSERTL0(smin >= 0.0,
             "ImplicitAbsorption must be >= 0 (do not negate as in "
             "ForcingAbsorption).");
}

void VCSImplicitAbsorption::v_SolveViscous(
    const Array<OneD, const Array<OneD, NekDouble>> &Forcing,
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble aii_Dt)
{
    const int nq = m_fields[0]->GetNpoints();

    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varCoeffMap     = StdRegions::NullVarCoeffMap;
    StdRegions::VarFactorsMap varFactorsMap = StdRegions::NullVarFactorsMap;

    AppendSVVFactors(factors, varFactorsMap);
    ComputeGJPNormalVelocity(inarray, varCoeffMap);

    // Rebuild the mass-variation field only when aii_Dt changes (start-up
    // steps of the IMEX scheme); otherwise reuse the SAME array so the
    // assembled Helmholtz matrix is cached across time steps.
    if (aii_Dt != m_lastAiiDt)
    {
        // m(x) = 1 + sigma(x) * aii_Dt  (per-field nu enters via lambda)
        Vmath::Smul(nq, aii_Dt, m_sigma, 1, m_massVar, 1);
        Vmath::Sadd(nq, 1.0, m_massVar, 1, m_massVar, 1);
        m_lastAiiDt = aii_Dt;
    }
    varCoeffMap[StdRegions::eVarCoeffMass] = m_massVar;

    for (int i = 0; i < m_nConvectiveFields; ++i)
    {
        factors[StdRegions::eFactorLambda] = 1.0 / aii_Dt / m_diffCoeff[i];
        m_fields[i]->HelmSolve(Forcing[i], m_fields[i]->UpdateCoeffs(), factors,
                               varCoeffMap, varFactorsMap);
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
    }
}

void VCSImplicitAbsorption::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    VelocityCorrectionScheme::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Absorption",
                                "Implicit (variable-coefficient)");
    SolverUtils::AddSummaryItem(s, "Sigma max",
                                Vmath::Vmax(m_sigma.size(), m_sigma, 1));
}
} // namespace Nektar
