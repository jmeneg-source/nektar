///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCCathodeMS.h
//
// PEMFCSolver: cathode MULTI-SPECIES steady transport on the collapsed
// (pseudo-2D) Behrou/Model-A layer. Gate (a) of Model B stage B1: the
// case-04 configuration E rerun with real two-species transport
// (wO2, wH2O; N2 by closure) in place of the ad hoc dilution parameter
// DilutionMassRatio of PEMFCSteady. Part of the nektar++_fuel_cells
// project (EPUSP / RCGI).
//
// For each species i in {O2, H2O} (physical form, SI units):
//
//   (V + u_St).grad(w_i) - div( eps_i(x,w) f_D(x) grad(w_i) )
//       + k_i(x) w_i = S_i(x)
//
// with, per unit collapsed-layer volume and mixture density
// (ivol = i_BV / (ne F rho Htot), i_BV from Kinetics/ButlerVolmer.h,
// identical closure and calibration to PEMFCSteady):
//
//   O2 : consumption  M_O2 * ivol       (secant-implicit, as PEMFCSteady)
//        dilution     (rw - 1) * M_O2 * ivol * w_O2   (linear sink)
//   H2O: source       rw * M_O2 * ivol                (explicit)
//        dilution     (rw - 1) * M_O2 * ivol * w_H2O  (linear sink)
//
// rw = WaterSourceRatio = (gas-phase H2O source mass)/(O2 sink mass);
// ground truth 3.375 (production 1.125 + electro-osmotic drag and
// membrane flux 2.25, thesis COMSOL R_wH2O). The dilution sink is the
// -w_i S_tot / rho term of the mixture species equation with
// S_tot = (rw - 1) * (O2 sink); with the composition effects switched
// off this degenerates EXACTLY to PEMFCSteady with
// DilutionMassRatio = rw - 1.
//
// Composition effects (the real multi-species content):
//  - eps_i(x,w) = epsilon * D_i^mix(x_local) / D_O2^mix(x_inlet):
//    mixture-averaged diffusivities from the local mole fractions,
//    normalized so the inlet composition reproduces the calibrated
//    Model A epsilon. D_i^mix = (1 - x_i) / sum_j x_j / D_ij.
//  - Stefan/correction velocity (StefanCorrection = 1):
//    u_St = f_D * [ (D_O2 - D_N2) grad wO2 + (D_H2O - D_N2) grad wH2O ]
//    (the mixture-average mass-flux closure does not sum to zero; the
//    correction restores sum J_i = 0 and carries the counter-diffusion
//    drag of outgoing water against incoming oxygen under the ribs).
//
// Session: EQTYPE PEMFCCathodeMS, variables wO2, wH2O (this order);
// functions AdvectionVelocity (Vx,Vy fld, physical), DiffusivityFactor
// (f_D, evaluated for wO2, applied to both species); parameters
// epsilon, rho, Htot, MO2, MH2O, MN2, ne, kinetics set as PEMFCSteady,
// XO2In, XH2OIn (inlet mole fractions for the normalization), DO2H2O,
// DO2N2, DN2H2O (only their ratios matter), WaterSourceRatio,
// StefanCorrection, PicardIters, PicardTol, PicardRelax, WO2In, WH2OIn
// (initial guess). Extra outputs: iloc, DO2fac (local eps_O2/epsilon).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCCATHODEMS_H
#define NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCCATHODEMS_H

#include <PEMFCSolver/Kinetics/ButlerVolmer.h>
#include <SolverUtils/EquationSystem.h>

namespace Nektar
{
class PEMFCCathodeMS : public SolverUtils::EquationSystem
{
public:
    friend class MemoryManager<PEMFCCathodeMS>;

    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<PEMFCCathodeMS>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    static std::string className;

    virtual ~PEMFCCathodeMS();

protected:
    PEMFCCathodeMS(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void v_InitObject(bool DeclareFields = true) override;
    void v_DoInitialise(bool dumpInitialConditions = true) override;
    void v_DoSolve() override;
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s);
    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        std::vector<std::string> &variables);

private:
    void SolveField(int f, const Array<OneD, Array<OneD, NekDouble>> &vel,
                    const Array<OneD, NekDouble> &diff,
                    const Array<OneD, NekDouble> &kfld,
                    const Array<OneD, NekDouble> &rhs);

    Array<OneD, Array<OneD, NekDouble>> m_velocity; // base advection
    Array<OneD, NekDouble> m_fD;                    // f_D(x) mask
    Array<OneD, NekDouble> m_iloc;
    Array<OneD, NekDouble> m_DO2fac; // eps_O2 local/eps
    PEMFC::ButlerVolmerORR m_kinetics;

    NekDouble m_epsilon, m_rho, m_Htot;
    NekDouble m_MO2, m_MH2O, m_MN2;
    int m_ne;
    NekDouble m_xO2in, m_xH2Oin;
    NekDouble m_DO2H2O, m_DO2N2, m_DN2H2O;
    NekDouble m_rw;     // WaterSourceRatio
    NekDouble m_stefan; // 0/1 switch
    NekDouble m_wO2in, m_wH2Oin;
    int m_picardIters;
    NekDouble m_picardTol, m_picardRelax;
};
} // namespace Nektar

#endif
