////////////////////////////////////////////////////////////////////////////////
//
//  File: RefRegionParallelogram.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Parallelogram surface for the refinement region.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/RefRegionParallelogram.h>

namespace Nektar::SpatialDomains
{

RefRegionParallelogram::RefRegionParallelogram(
    const unsigned int coordim, NekDouble radius, std::vector<NekDouble> coord1,
    std::vector<NekDouble> coord2, std::vector<unsigned int> numModes,
    std::vector<unsigned int> numPoints)
    : RefRegion(coordim, radius, coord1, coord2, numModes, numPoints)
{
}

/**
 * @brief Check if vertex is inside the Parallelogram.
 *
 * @param coords    coordinates of the vertex
 * @return          true or false depending on if the vertex is inside
 *                  or not of the surface defined by the user.
 */
bool RefRegionParallelogram::v_Contains(const Array<OneD, NekDouble> &coords)
{
    // This is simplification for a two-dimenion domain of the algorithm in the
    // RefRegionCylinder::v_Contains method.

    const size_t dim = coords.size();
    Array<OneD, NekDouble> e(dim, 0.0);     // direction: rb - ra
    Array<OneD, NekDouble> m(dim - 1, 0.0); // momentum: ra x rb
    NekDouble d = 0.0;                      // distance

    // Direction
    e[0] = m_coord2[0] - m_coord1[0];
    e[1] = m_coord2[1] - m_coord1[1];

    // Cross product of vectors 'coord1'(ra) and 'coord2' (rb)
    m[0] = m_coord1[0] * m_coord2[1] - m_coord1[1] * m_coord2[0];

    // Distance of P (coords) to line AB (coord1coord2) is equal or less than R
    // d = || e x (rp - ra) || / || e ||

    // rA - rP
    Array<OneD, NekDouble> rpa(dim, 0.0);
    rpa[0] = coords[0] - m_coord1[0];
    rpa[1] = coords[1] - m_coord1[1];

    // || e ||
    NekDouble e_mod = sqrt(e[0] * e[0] + e[1] * e[1]);

    // || e x (rp - ra) ||
    Array<OneD, NekDouble> exrpa(dim - 1, 0.0);
    exrpa[0] = e[0] * rpa[1] - e[1] * rpa[0];

    NekDouble exrpa_mod = sqrt(exrpa[0] * exrpa[0]);

    d = exrpa_mod / e_mod;
    if (d >= m_radius)
    {
        return false;
    }

    // Is P between the regions below?
    // xa, xb plus the radius
    // ya, yb plus the radius
    Array<OneD, bool> insideFlag(dim, false);
    for (int i = 0; i < dim; ++i)
    {
        if (m_coord1[i] < m_coord2[i])
        {
            if ((m_coord1[i] - m_radius < coords[i]) &&
                (m_coord2[i] + m_radius > coords[i]))
            {
                insideFlag[i] = true;
            }
        }
        else
        {
            if ((m_coord1[i] + m_radius > coords[i]) &&
                (m_coord2[i] - m_radius < coords[i]))
            {
                insideFlag[i] = true;
            }
        }
    }

    if ((insideFlag[0] == true) && (insideFlag[1] == true))
    {
        return true;
    }
    else
    {
        return false;
    }
}

} // namespace Nektar::SpatialDomains
