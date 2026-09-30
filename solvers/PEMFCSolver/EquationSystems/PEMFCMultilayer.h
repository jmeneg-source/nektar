///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCMultilayer.h
//
// PEMFCSolver: model-B stage B1 skeleton. Pseudo-3D multilayer PEMFC:
// four gas layers (cathode channel cc, cathode GDL cg, anode channel
// ac, anode GDL ag) collapsed onto ONE 2D mesh. Isothermal,
// single-phase, dry membrane. Solves EIGHT species mass fractions
//
//   wO2_cc, wH2O_cc, wO2_cg, wH2O_cg, wH2_ac, wH2O_ac, wH2_ag, wH2O_ag
//
// by Picard iteration on variable-coefficient steady ADR solves (same
// machinery as PEMFCSteady), with frozen per-layer velocities and the
// ground-truth inter-layer exchange closures:
//
//   - diffusive z-flux (var1/var5): +-design*Dm_eff*2*
//       (rho_ch w_ch - rho_gdl w_gdl)/(Hch+Hgdl), split into an
//       implicit diagonal sink and a lagged off-diagonal source;
//   - optional convective z-flux (var4/var7) when layer pressure
//     fields are supplied (off by default in B1);
//   - reaction fluxes (var14): Rk = -+ M_s nu_s/(ne F) * iloc in the
//     GDL layers;
//   - electrochemistry closure: Kinetics/MultilayerKinetics.h
//     (ACTIVE var18/var20 model), pointwise per quadrature point.
//
// Session: EQTYPE PEMFCMultilayer, the eight variables above (order
// fixed); functions Velocity_cc/_cg/_ac/_ag (Vx,Vy; default zero),
// DesignCathode/DesignAnode (default 1); parameters from the
// ground-truth set (all defaulted); PicardIters/Tol/Relax.
// Extra output fields: iloc, etaC, etaAn.
//
// Part of the nektar++_fuel_cells project (EPUSP / RCGI / BUGA-LAB).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCMULTILAYER_H
#define NEKTAR_SOLVERS_PEMFCSOLVER_EQUATIONSYSTEMS_PEMFCMULTILAYER_H

#include <PEMFCSolver/Kinetics/MultilayerKinetics.h>
#include <SolverUtils/EquationSystem.h>

namespace Nektar
{
class PEMFCMultilayer : public SolverUtils::EquationSystem
{
public:
    friend class MemoryManager<PEMFCMultilayer>;

    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<PEMFCMultilayer>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    static std::string className;

    virtual ~PEMFCMultilayer();

    // Field indices (fixed order).
    enum FieldIdx
    {
        eO2cc = 0,
        eH2Occ,
        eO2cg,
        eH2Ocg,
        eH2ac,
        eH2Oac,
        eH2ag,
        eH2Oag,
        eNumFields
    };

protected:
    PEMFCMultilayer(const LibUtilities::SessionReaderSharedPtr &pSession,
                    const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void v_InitObject(bool DeclareFields = true) override;
    void v_DoInitialise(bool dumpInitialConditions = true) override;
    void v_DoSolve() override;
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s);
    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        std::vector<std::string> &variables);

private:
    // One steady variable-coefficient ADR solve for field f with
    // pointwise sink kfld [1/s] and source rhs (physical units of
    // w per second), advected by m_vel[layer of f].
    void SolveField(int f, const Array<OneD, NekDouble> &kfld,
                    const Array<OneD, NekDouble> &rhs);

    // Pointwise closure sweep: mixture molar masses, partial
    // pressures, kinetics -> m_iloc/m_etaC/m_etaAn.
    void UpdateClosures();

    // Per-layer (4) frozen velocities, [layer][dim][nq].
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_vel;
    Array<OneD, NekDouble> m_designC; // cathode mask
    Array<OneD, NekDouble> m_designA; // anode mask
    Array<OneD, NekDouble> m_iloc;
    Array<OneD, NekDouble> m_etaC;
    Array<OneD, NekDouble> m_etaAn;

    PEMFC::MultilayerKinetics m_kin;

    // Geometry / properties (ground-truth defaults).
    NekDouble m_Hch, m_Hgdl;
    NekDouble m_pAbs; // absolute pressure [Pa]
    NekDouble m_T;
    NekDouble m_rhoC, m_rhoA;                        // layer densities [kg/m^3]
    NekDouble m_DO2H2O, m_DO2N2, m_DN2H2O, m_DH2H2O; // ref binary D
    NekDouble m_epsGdl;                              // GDL porosity (Bruggeman)
    // Molar masses.
    static constexpr double kMO2  = 0.032;
    static constexpr double kMH2  = 0.002;
    static constexpr double kMN2  = 0.028;
    static constexpr double kMH2O = 0.018;

    int m_picardIters;
    NekDouble m_picardTol;
    NekDouble m_picardRelax;
};
} // namespace Nektar

#endif
