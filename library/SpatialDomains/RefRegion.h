////////////////////////////////////////////////////////////////////////////////
//
//  File: RefRegion.h
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
//  Description: Abstract base class for the refinement region.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_REFREGION_H
#define NEKTAR_SPATIALDOMAINS_REFREGION_H

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <vector>

namespace Nektar::SpatialDomains
{

/**
 * @class RefRegion
 * @brief Abstract base class for the refinement surface region.
 */
class RefRegion
{
public:
    /// Constructor
    SPATIAL_DOMAINS_EXPORT RefRegion(const unsigned int coordim,
                                     NekDouble m_radius,
                                     std::vector<NekDouble> coord1,
                                     std::vector<NekDouble> coord2,
                                     std::vector<unsigned int> numModes,
                                     std::vector<unsigned int> numPoints);
    /// Destructor
    SPATIAL_DOMAINS_EXPORT virtual ~RefRegion() = default;

    /// Pure virtual fuction
    SPATIAL_DOMAINS_EXPORT virtual bool v_Contains(
        const Array<OneD, NekDouble> &coords) = 0;

    /// Get the number of modes to update expansion
    SPATIAL_DOMAINS_EXPORT std::vector<unsigned int> GetNumModes()
    {
        return m_numModes;
    }

    /// Get the number of quadrature points to update expansion
    SPATIAL_DOMAINS_EXPORT std::vector<unsigned int> GetNumPoints()
    {
        return m_numPoints;
    }

protected:
    /// Coordinate 1
    std::vector<NekDouble> m_coord1;
    /// Coordinate 2
    std::vector<NekDouble> m_coord2;
    /// Number of modes
    std::vector<unsigned int> m_numModes;
    /// Number of quadrature points
    std::vector<unsigned int> m_numPoints;
    /// Radius of the surface region
    NekDouble m_radius;
    /// Dimension of the coordinate (space dimension)
    unsigned int m_coordim;
};

} // namespace Nektar::SpatialDomains

#endif // NEKTAR_SPATIALDOMAINS_REFREGION_H
