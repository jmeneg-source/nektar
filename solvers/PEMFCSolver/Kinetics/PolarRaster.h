///////////////////////////////////////////////////////////////////////////////
//
// File: PolarRaster.h
//
// Dependency-free polar raster field gamma(r,theta) for annular-sector
// design masks (case 05, circular radial cell). Part of the
// nektar++_fuel_cells project (EPUSP / RCGI / BUGA-LAB).
//
// File format: ASCII matrix, '#' comment lines ignored; rows = radial
// samples (cell centers, Rin..Rout), columns = angular samples (cell
// centers, theta0..theta1). Values in [0,1]: 1 = channel, 0 = rib.
// Evaluation: bilinear interpolation in (r,theta), clamped at the
// raster boundary. Points outside [Rin,Rout]x[theta0,theta1] are
// clamped (the sector mesh should match the raster extent).
//
// Session convention (both PEMFCSteady and VCSImplicitAbsorption):
//   SOLVERINFO  DesignRaster       = <file>
//   PARAMETERS  RasterRin, RasterRout, RasterThetaDeg (default 60)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PEMFCSOLVER_KINETICS_POLARRASTER_H
#define NEKTAR_SOLVERS_PEMFCSOLVER_KINETICS_POLARRASTER_H

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace Nektar
{
namespace PEMFC
{

class PolarRaster
{
public:
    int m_nr     = 0;
    int m_nt     = 0;
    double m_rin = 0.0, m_rout = 1.0;
    double m_th0 = 0.0, m_th1 = 0.0;
    std::vector<double> m_v; // row-major: [i_r * m_nt + j_t]

    bool Load(const std::string &fname, double rin, double rout, double th0,
              double th1)
    {
        std::ifstream f(fname.c_str());
        if (!f.good())
        {
            return false;
        }
        m_rin  = rin;
        m_rout = rout;
        m_th0  = th0;
        m_th1  = th1;
        m_v.clear();
        m_nr = 0;
        m_nt = 0;
        std::string line;
        while (std::getline(f, line))
        {
            if (line.empty() || line[0] == '#')
            {
                continue;
            }
            std::istringstream ss(line);
            double x;
            int cols = 0;
            while (ss >> x)
            {
                m_v.push_back(x);
                ++cols;
            }
            if (cols == 0)
            {
                continue;
            }
            if (m_nt == 0)
            {
                m_nt = cols;
            }
            else if (cols != m_nt)
            {
                return false; // ragged rows
            }
            ++m_nr;
        }
        return (m_nr > 1 && m_nt > 1);
    }

    // gamma at physical point (x,y); theta measured from +x axis.
    double Eval(double x, double y) const
    {
        const double r = std::sqrt(x * x + y * y);
        double th      = std::atan2(y, x);
        // cell-center coordinates
        double fr = (r - m_rin) / (m_rout - m_rin) * m_nr - 0.5;
        double ft = (th - m_th0) / (m_th1 - m_th0) * m_nt - 0.5;
        if (fr < 0.0)
            fr = 0.0;
        if (ft < 0.0)
            ft = 0.0;
        if (fr > m_nr - 1.0)
            fr = m_nr - 1.0;
        if (ft > m_nt - 1.0)
            ft = m_nt - 1.0;
        const int i0    = static_cast<int>(fr);
        const int j0    = static_cast<int>(ft);
        const int i1    = (i0 < m_nr - 1) ? i0 + 1 : i0;
        const int j1    = (j0 < m_nt - 1) ? j0 + 1 : j0;
        const double ai = fr - i0;
        const double aj = ft - j0;
        return (1 - ai) *
                   ((1 - aj) * m_v[i0 * m_nt + j0] + aj * m_v[i0 * m_nt + j1]) +
               ai * ((1 - aj) * m_v[i1 * m_nt + j0] + aj * m_v[i1 * m_nt + j1]);
    }
};

} // namespace PEMFC
} // namespace Nektar

#endif
