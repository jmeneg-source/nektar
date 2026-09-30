///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCSteady.h
//
// PEMFCSolver: STEADY oxygen transport with nonlinear Butler-Volmer
// kinetics, solved by Picard iteration on variable-coefficient linear
// steady advection-diffusion-reaction solves. Part of the
// nektar++_fuel_cells project (EPUSP / RCGI).
//
// Each Picard iterate solves, for w = wO2 (physical form, SI units):
//
//   V.grad(w) - div( eps f_D(x) grad(w) ) + k_m(x) w = 0
//   k_m(x) = (M_O2 / (ne F rho Htot)) * i_BV(w_m(x)) / w_m(x)
//
// with i_BV the full nonlinear closure (Kinetics/ButlerVolmer.h). The
// secant linearization k_m makes each iterate exact for the previous
// field; convergence in a handful of iterations at the operating points
// of interest (weak nonlinearity through the ohmic feedback).
//
// Requires the two platform patches shipped with this project:
//   - StdExpansion::LinearAdvectionDiffusionReactionMatrixOp_MatFree:
//     velocity-only direction counting + eVarCoeffMass reaction;
//   - ContField::LinearAdvectionDiffusionReactionSolveVarCoeff.
//
// Session: EQTYPE PEMFCSteady, variable wO2; functions AdvectionVelocity
// (Vx,Vy from fld, UNSCALED physical velocity), optional
// DiffusivityFactor (f_D); parameters epsilon, rho, Htot, MO2, ne and
// the kinetics set; PicardIters, PicardTol, PicardRelax.
// iloc written as an extra output field.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCSTEADY_H
#define NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCSTEADY_H

#include <PEMFCSolver/Kinetics/ButlerVolmer.h>
#include <PEMFCSolver/Kinetics/PolarRaster.h>
#include <SolverUtils/EquationSystem.h>

namespace Nektar
{
class PEMFCSteady : public SolverUtils::EquationSystem
{
public:
    friend class MemoryManager<PEMFCSteady>;

    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<PEMFCSteady>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    static std::string className;

    virtual ~PEMFCSteady();

protected:
    PEMFCSteady(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void v_InitObject(bool DeclareFields = true) override;
    void v_DoInitialise(bool dumpInitialConditions = true) override;
    void v_DoSolve() override;
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s);
    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        std::vector<std::string> &variables);

private:
    Array<OneD, Array<OneD, NekDouble>> m_velocity;
    Array<OneD, NekDouble> m_diffCoeff; // eps*f_D(x)
    Array<OneD, NekDouble> m_gamma;     // design mask
    Array<OneD, NekDouble> m_iloc;
    PEMFC::ButlerVolmerORR m_kinetics;

    NekDouble m_epsilon;
    NekDouble m_rho;
    NekDouble m_Htot;
    NekDouble m_MO2;
    int m_ne;
    NekDouble m_dilRatio; // net-mass-source / O2-consumption ratio
    NekDouble m_winlet;   // inlet mass fraction (initial guess)
    int m_picardIters;
    NekDouble m_picardTol;
    NekDouble m_picardRelax;
};
} // namespace Nektar

#endif
