///////////////////////////////////////////////////////////////////////////////
//
// File: VCSImplicitAbsorption.h
//
// Part of the nektar++_fuel_cells project (EPUSP / RCGI).
//
// Velocity-correction scheme with IMPLICIT spatially-varying absorption
// (Brinkman/Darcy penalization) folded into the viscous Helmholtz solve:
//
//   viscous step:  ( nabla^2 - lambda * m(x) ) u* = F,
//   m(x) = 1 + sigma(x) * aii_Dt,   lambda = 1/(aii_Dt * nu),
//
// equivalent to du/dt = ... - sigma(x) u treated implicitly, which removes
// the explicit-forcing stability limit sigma*dt <~ 1 and admits physical
// rib penalizations sigma = mu/kappa ~ 5e6 1/s (kappa_pr = 3.8e-12 m2).
//
// The coefficient field sigma(x) [1/s] is read once from the session
// function "ImplicitAbsorption" (variables u,v). Select with
//   SOLVERINFO SolverType = "VCSImplicitAbsorption".
// Motivated quantitatively by case 04: the explicit cap sigma = 1e5
// leaked ~4 mm/s through the ribs and washed out the O2 depletion in the
// rib shadows (see the project development log, 15/08/2026).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VCSIMPLICITABSORPTION_H
#define NEKTAR_SOLVERS_VCSIMPLICITABSORPTION_H

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>

namespace Nektar
{
class VCSImplicitAbsorption : public VelocityCorrectionScheme
{
public:
    friend class MemoryManager<VCSImplicitAbsorption>;

    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<VCSImplicitAbsorption>::AllocateSharedPtr(pSession,
                                                                    pGraph);
        p->InitObject();
        return p;
    }

    static std::string className;
    static std::string solverTypeLookupId;

    VCSImplicitAbsorption(const LibUtilities::SessionReaderSharedPtr &pSession,
                          const SpatialDomains::MeshGraphSharedPtr &pGraph);

    virtual ~VCSImplicitAbsorption();

protected:
    void v_InitObject(bool DeclareField = true) override;
    void v_SolveViscous(
        const Array<OneD, const Array<OneD, NekDouble>> &Forcing,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble aii_Dt) override;
    void v_GenerateSummary(SolverUtils::SummaryList &s) override;

private:
    Array<OneD, NekDouble> m_sigma;   // absorption coefficient [1/s]
    Array<OneD, NekDouble> m_massVar; // 1 + sigma*aii_Dt (cached)
    NekDouble m_lastAiiDt = -1.0;
};
} // namespace Nektar

#endif
