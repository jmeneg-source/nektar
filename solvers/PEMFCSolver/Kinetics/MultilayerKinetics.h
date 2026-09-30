///////////////////////////////////////////////////////////////////////////////
//
// File: MultilayerKinetics.h
//
// Dependency-free pointwise electrochemistry closure for the model-B
// pseudo-3D multilayer PEMFC solver (PEMFCMultilayer). Implements the
// ACTIVE current model of the ground-truth .mph (groups var18/var20):
// surface-style Butler-Volmer on both electrodes with Nernst pressure
// products, coupled through the cell voltage balance
//
//   eta_anode = (Eeq_c - Eeq_a) - E_cell - |eta_c|
//               - deltaV_ohm_c - deltaV_ohm_a,
//   deltaV_ohm_c = deltaV_ohm_a = (Hch/sigma_ch + Hgdl/sigma_gdl)*|iloc|,
//
// and the local current constraint iloc_c(eta_c) + iloc_a(eta_anode) = 0,
// iloc = |iloc_c|. Solved pointwise: outer fixed-point on the ohmic drop,
// inner safeguarded bisection/Newton on eta_c.
//
// The DISABLED agglomerate model (var13) and the volumetric exchange
// currents (i0 in A/m^3, alpha = 0.7) are intentionally NOT used here.
//
// Part of the nektar++_fuel_cells project (EPUSP / RCGI / BUGA-LAB).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PEMFCSOLVER_KINETICS_MULTILAYERKINETICS_H
#define NEKTAR_SOLVERS_PEMFCSOLVER_KINETICS_MULTILAYERKINETICS_H

#include <algorithm>
#include <cmath>

namespace Nektar
{
namespace PEMFC
{

class MultilayerKinetics
{
public:
    // Ground-truth defaults (Paper_2_p3D_v7.mph). All SI.
    double m_T        = 353.15;   // cell temperature [K]
    double m_Ecell    = 0.8;      // cell voltage [V]
    double m_i0c      = 1.0e4;    // i0_ref_c [A/m^2]
    double m_i0a      = 2.0e8;    // i0_ref_a [A/m^2]
    double m_alpha    = 0.5;      // alphaa = alphac (both electrodes)
    double m_EeqRefC  = 1.19;     // Eeq_ref_gdl_c [V]  (TODO: thermo fit)
    double m_EeqRefA  = 0.0;      // Eeq_ref_gdl_a [V]  (TODO: thermo fit)
    double m_Hch      = 1.0e-3;   // channel height [m]
    double m_Hgdl     = 380.0e-6; // GDL height [m]
    double m_sigmaCh  = 330.0;    // graphite through-plane [S/m]
    double m_sigmaGdl = 500.0;    // GDL through-plane [S/m]
    double m_patm     = 101325.0; // [Pa]
    double m_tol      = 1.0e-6;

    static constexpr double kF  = 96485.3;
    static constexpr double kRg = 8.3144;

    // Nernst-corrected equilibrium potentials from GDL partial
    // pressures [Pa]. Cathode: Pred = (p_h2o/atm)^2, Pox = p_o2/atm.
    // Anode: Pred = p_h2/atm, Pox = 1. (var23 conventions, both use
    // the 1/4F log prefactor as in the .mph.)
    double EeqC(double pO2, double pH2Oc) const
    {
        const double pr  = PredC(pH2Oc);
        const double pox = PoxC(pO2);
        return m_EeqRefC - kRg * m_T / (4.0 * kF) * std::log(pr / pox);
    }
    double EeqA(double pH2) const
    {
        const double pr = PredA(pH2);
        return m_EeqRefA - kRg * m_T / (4.0 * kF) * std::log(pr);
    }

    double PredC(double pH2Oc) const
    {
        const double x = pH2Oc / m_patm;
        return (x > 1.0e-6) ? x * x : x * 1.0e-6;
    }
    double PoxC(double pO2) const
    {
        return std::max(pO2 / m_patm, m_tol);
    }
    double PredA(double pH2) const
    {
        return std::max(pH2 / m_patm, m_tol);
    }

    static double SafeExp(double a)
    {
        return std::exp(std::min(std::max(a, -500.0), 500.0));
    }
    double IlocC(double eta, double pO2, double pH2Oc) const
    {
        const double b = m_alpha * kF / (kRg * m_T);
        return m_i0c * (PredC(pH2Oc) * SafeExp(b * eta) -
                        PoxC(pO2) * SafeExp(-b * eta));
    }
    double IlocA(double etaAn, double pH2) const
    {
        const double b = m_alpha * kF / (kRg * m_T);
        return m_i0a *
               (PredA(pH2) * SafeExp(b * etaAn) - 1.0 * SafeExp(-b * etaAn));
    }

    double OhmicDrop(double iloc) const
    {
        return (m_Hch / m_sigmaCh + m_Hgdl / m_sigmaGdl) * std::fabs(iloc);
    }

    // Pointwise solve. Inputs: GDL partial pressures [Pa] on both sides,
    // warm-start eta_c. Outputs through references; returns iloc >= 0.
    double Solve(double pO2, double pH2Oc, double pH2, double &etaC,
                 double &etaAn) const
    {
        const double dEeq = EeqC(pO2, pH2Oc) - EeqA(pH2);
        const double R2   = 2.0 * (m_Hch / m_sigmaCh + m_Hgdl / m_sigmaGdl);

        // Root in x = |eta_c| >= 0 of the strictly decreasing
        //   H(x) = IlocA(etaAn(x)) - |IlocC(-x)|,
        //   etaAn(x) = dEeq - E_cell - x - R2*|IlocC(-x)|,
        // i.e. the ohmic drop is evaluated self-consistently at the
        // cathode current of the trial x (no outer fixed point: the
        // exponential sensitivity makes a lagged ohmic loop diverge).
        // H(0) <= 0 means the point cannot sustain positive current
        // (starved or reversed): iloc = 0.
        auto icath = [&](double x) { return std::fabs(IlocC(-x, pO2, pH2Oc)); };
        auto H     = [&](double x) {
            const double io = icath(x);
            const double ea = dEeq - m_Ecell - x - R2 * io;
            return IlocA(ea, pH2) - io;
        };
        double lo = 0.0, hi = 2.0;
        if (H(lo) <= 0.0)
        {
            etaC  = 0.0;
            etaAn = dEeq - m_Ecell;
            return 0.0;
        }
        while (H(hi) > 0.0 && hi < 32.0)
        {
            hi *= 2.0; // pathological, but keep the bracket sound
        }
        for (int it = 0; it < 100; ++it)
        {
            const double mid = 0.5 * (lo + hi);
            if (H(mid) > 0.0)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
            if (hi - lo < 1.0e-13)
                break;
        }
        const double x    = 0.5 * (lo + hi);
        const double iloc = icath(x);
        etaC              = -x;
        etaAn             = dEeq - m_Ecell - x - R2 * iloc;
        return iloc;
    }
};

} // namespace PEMFC
} // namespace Nektar

#endif
