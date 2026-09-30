////////////////////////////////////////////////////////////////////////////////
//
//  File: RefRegionCylinder.cpp
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
//  Description: Cylindrical surface for the refinement region.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/RefRegionCylinder.h>

namespace Nektar::SpatialDomains
{

RefRegionCylinder::RefRegionCylinder(const unsigned int coordim,
                                     NekDouble radius,
                                     std::vector<NekDouble> coord1,
                                     std::vector<NekDouble> coord2,
                                     std::vector<unsigned int> numModes,
                                     std::vector<unsigned int> numPoints)
    : RefRegion(coordim, radius, coord1, coord2, numModes, numPoints)
{
}

/**
 * @brief Check if vertex is inside the cylinder.
 *
 * @param coords    coordinates of the vertex
 * @return          true or false depending on if the vertex is inside
 *                  or not of the surface defined by the user.
 */
bool RefRegionCylinder::v_Contains(const Array<OneD, NekDouble> &coords)
{
    const size_t dim = coords.size();
    Array<OneD, NekDouble> e(dim, 0.0); // direction: rb - ra
    Array<OneD, NekDouble> m(dim, 0.0); // momentum: ra x rb
    NekDouble d = 0.0;                  // distance

    // Direction
    e[0] = m_coord2[0] - m_coord1[0];
    e[1] = m_coord2[1] - m_coord1[1];
    e[2] = m_coord2[2] - m_coord1[2];

    // Cross product of vectors 'coord1' and 'coord2'
    m[0] = m_coord1[1] * m_coord2[2] - m_coord1[2] * m_coord2[1];
    m[1] = m_coord1[2] * m_coord2[0] - m_coord1[0] * m_coord2[2];
    m[2] = m_coord1[0] * m_coord2[1] - m_coord1[1] * m_coord2[0];

    // 1. Distance of P (coords) to line AB (coord1coord2) is equal or less than
    // R
    //  d = || e x (rp - ra) || / || e ||

    // rA - rP
    Array<OneD, NekDouble> rpa(dim, 0.0);
    rpa[0] = coords[0] - m_coord1[0];
    rpa[1] = coords[1] - m_coord1[1];
    rpa[2] = coords[2] - m_coord1[2];

    // || e ||
    NekDouble e_mod = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);

    // || e x (rp - ra) ||
    Array<OneD, NekDouble> exrpa(dim, 0.0);
    exrpa[0] = e[1] * rpa[2] - e[2] * rpa[1];
    exrpa[1] = e[2] * rpa[0] - e[0] * rpa[2];
    exrpa[2] = e[0] * rpa[1] - e[1] * rpa[0];

    NekDouble exrpa_mod =
        sqrt(exrpa[0] * exrpa[0] + exrpa[1] * exrpa[1] + exrpa[2] * exrpa[2]);

    d = exrpa_mod / e_mod;
    if (d >= m_radius)
    {
        return false;
    }

    // 2. Closest point Q on line AB to P
    // rq = rp + (e x (m + e x rp)) / ||e||^{2}

    // (m + e x rp)
    Array<OneD, NekDouble> mpexrp(dim, 0.0);
    mpexrp[0] = m[0] + e[1] * coords[2] - e[2] * coords[1];
    mpexrp[1] = m[1] + e[2] * coords[0] - e[0] * coords[2];
    mpexrp[2] = m[2] + e[0] * coords[1] - e[1] * coords[0];

    // e x (m + e x rp) = numerator
    Array<OneD, NekDouble> numerator(dim, 0.0);
    numerator[0] = e[1] * mpexrp[2] - e[2] * mpexrp[1];
    numerator[1] = e[2] * mpexrp[0] - e[0] * mpexrp[2];
    numerator[2] = e[0] * mpexrp[1] - e[1] * mpexrp[0];

    // rq
    Array<OneD, NekDouble> rq(dim, 0.0);
    rq[0] = coords[0] + (numerator[0] / pow(e_mod, 2));
    rq[1] = coords[1] + (numerator[1] / pow(e_mod, 2));
    rq[2] = coords[2] + (numerator[2] / pow(e_mod, 2));

    // 3. The baricentric coordinates of Q(wa,wb) such that rq = wa*ra+wb*rb
    //  wa = ||rq x rb||/||m||, wb = ||rq x ra||/||m||
    NekDouble wa = 0.0;
    NekDouble wb = 0.0;

    // ||m||
    NekDouble m_mod = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
    // m_mod > 0: condition
    ASSERTL0(m_mod, "The cylinder axis must not go through the origin");

    // ||rq x rb||
    Array<OneD, NekDouble> rqxrb(dim, 0.0);
    rqxrb[0] = rq[1] * m_coord2[2] - rq[2] * m_coord2[1];
    rqxrb[1] = rq[2] * m_coord2[0] - rq[0] * m_coord2[2];
    rqxrb[2] = rq[0] * m_coord2[1] - rq[1] * m_coord2[0];
    NekDouble rqxrb_mod =
        sqrt(rqxrb[0] * rqxrb[0] + rqxrb[1] * rqxrb[1] + rqxrb[2] * rqxrb[2]);

    // ||rq x ra||
    Array<OneD, NekDouble> rqxra(dim, 0.0);
    rqxra[0] = rq[1] * m_coord1[2] - rq[2] * m_coord1[1];
    rqxra[1] = rq[2] * m_coord1[0] - rq[0] * m_coord1[2];
    rqxra[2] = rq[0] * m_coord1[1] - rq[1] * m_coord1[0];
    NekDouble rqxra_mod =
        sqrt(rqxra[0] * rqxra[0] + rqxra[1] * rqxra[1] + rqxra[2] * rqxra[2]);

    // wa
    wa = rqxrb_mod / m_mod;
    // wb
    wb = rqxra_mod / m_mod;

    // 4. Check that point Q between A and B by making sure the baricentric
    // coordinates are between 0 and 1.
    if ((wa >= 0) && (wa <= 1) && (wb >= 0) && (wb <= 1))
    {
        return true;
    }
    else
    {
        return false;
    }
}

} // namespace Nektar::SpatialDomains
