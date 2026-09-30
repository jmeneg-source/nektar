///////////////////////////////////////////////////////////////////////////////
//
// File: ButlerVolmer.h
//
// PEMFCSolver: pluggable electrochemical kinetics, first module.
// Part of the nektar++_fuel_cells project (EPUSP / RCGI).
//
// Butler-Volmer kinetics for the oxygen reduction reaction with a lumped
// ohmic overpotential closure, after Behrou et al., IJHMT 135 (2019) and
// Razmara et al., SMO 67 (2024):
//
//   i_BV(i)  = i0 * (C_O2 / C_O2ref) * [ exp( k1*eta ) - exp( -k2*eta ) ]
//   eta(i)   = Voc - Vcell - Rlump * i
//   k1 = alpha_c F / (Rg T),   k2 = (1 - alpha_c) F / (Rg T)
//
// The fixed point i = i_BV(i) is solved pointwise by a safeguarded Newton
// iteration on f(i) = i - i_BV(i). f is monotonically increasing in i
// (Rlump > 0), so the root is unique; convergence from i = 0 is robust.
//
// This class is deliberately free of Nektar++ dependencies so that it can
// be unit-tested standalone and, later, templated on the scalar type for
// algorithmic differentiation (CoDiPack) in the topology optimization
// branch of the project.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PEMFCSOLVER_KINETICS_BUTLERVOLMER_H
#define NEKTAR_SOLVERS_PEMFCSOLVER_KINETICS_BUTLERVOLMER_H

#include <algorithm>
#include <cmath>

namespace Nektar
{
namespace PEMFC
{

class ButlerVolmerORR
{
public:
    // Parameters with defaults from Razmara et al. (2024), Table 2.1.
    double m_i0     = 0.17;    // exchange current density [A/m^2]
    double m_alphaC = 0.5;     // cathodic transfer coefficient [-]
    double m_T      = 343.15;  // cell temperature [K]
    double m_Voc    = 1.22;    // open circuit voltage [V]
    double m_Vcell  = 0.4;     // operating cell voltage [V]
    double m_Rlump  = 7.0e-5;  // lumped ohmic resistance [Ohm m^2]
    double m_CO2ref = 30.0;    // reference O2 concentration [mol/m^3]
    double m_F      = 96485.3; // Faraday constant [C/mol]
    double m_Rg     = 8.3144;  // universal gas constant [J/mol/K]

    // Local current density [A/m^2] for a given O2 molar concentration
    // [mol/m^3]. Newton on f(i) = i - i_BV(i), started from i = 0.
    double Solve(double CO2) const
    {
        const double C = std::max(CO2, 0.0) / m_CO2ref;
        if (C <= 0.0)
        {
            return 0.0;
        }
        const double k1 = m_alphaC * m_F / (m_Rg * m_T);
        const double k2 = (1.0 - m_alphaC) * m_F / (m_Rg * m_T);
        const double e0 = m_Voc - m_Vcell;

        double i = 0.0;
        for (int it = 0; it < 50; ++it)
        {
            const double eta = e0 - m_Rlump * i;
            const double ep  = std::exp(k1 * eta);
            const double em  = std::exp(-k2 * eta);
            const double f   = i - m_i0 * C * (ep - em);
            const double fp  = 1.0 + m_i0 * C * m_Rlump * (k1 * ep + k2 * em);
            const double di  = f / fp;
            i -= di;
            if (std::fabs(di) < 1.0e-10 * (1.0 + std::fabs(i)))
            {
                break;
            }
        }
        return std::max(i, 0.0);
    }
};

} // namespace PEMFC
} // namespace Nektar

#endif
