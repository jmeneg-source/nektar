///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCSpecies.h
//
// PEMFCSolver: coupled species transport with nonlinear Butler-Volmer
// kinetics on a prescribed (frozen) velocity field. Stage 1 of the
// monolithic PEM fuel cell solver of the nektar++_fuel_cells project.
//
// Solves, for w_k in {w_O2, w_H2O} (w_N2 by closure), CG projection:
//
//   dw_k/dt + V.grad(w_k) = eps * lap(w_k) + S_k(w_O2) / (rho * Htot)
//
//   S_O2  = -M_O2  * i_loc / (ne F)         [ORR: O2 + 4H+ + 4e- -> 2H2O]
//   S_H2O = +2 M_H2O * i_loc / (ne F)
//
// with i_loc from the pointwise nonlinear Butler-Volmer solve (see
// Kinetics/ButlerVolmer.h) and C_O2 = rho w_O2 / M_O2. Advection is
// explicit, diffusion implicit (IMEX), following the pattern of
// ADRSolver's UnsteadyAdvectionDiffusion. i_loc is written to output
// files as an extra field.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCSPECIES_H
#define NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCSPECIES_H

#include <PEMFCSolver/Kinetics/ButlerVolmer.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar
{
class PEMFCSpecies : public SolverUtils::AdvectionSystem
{
public:
    friend class MemoryManager<PEMFCSpecies>;

    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<PEMFCSpecies>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    static std::string className;

    virtual ~PEMFCSpecies();

protected:
    PEMFCSpecies(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void v_InitObject(bool DeclareFields = true) override;
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s);
    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        std::vector<std::string> &variables);

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void GetFluxVector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                       Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

private:
    Array<OneD, Array<OneD, NekDouble>> m_velocity;
    Array<OneD, NekDouble> m_iloc;       // local current density
    Array<OneD, NekDouble> m_diffFactor; // f_D(x)
    bool m_hasVarDiff = false;
    PEMFC::ButlerVolmerORR m_kinetics;

    NekDouble m_epsilon; // effective diffusivity
    NekDouble m_rho;     // mixture density (constant, stage 1)
    NekDouble m_Htot;    // total (channel + GDL) height
    NekDouble m_MO2;     // molar masses
    NekDouble m_MH2O;
    int m_idxO2; // variable indices
    int m_idxH2O;
    int m_ne; // participating electrons

    void ComputeReaction(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);
};
} // namespace Nektar

#endif
