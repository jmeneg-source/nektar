////////////////////////////////////////////////////////////////////////////////
//
//  File: RefRegionCylinder.h
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

#ifndef NEKTAR_SPATIALDOMAINS_REFREGIONCYLINDER_H
#define NEKTAR_SPATIALDOMAINS_REFREGIONCYLINDER_H

#include <SpatialDomains/RefRegion.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar::SpatialDomains
{

/**
 * @class RefRegionCylinder
 * @brief Derived class for the refinement surface region.
 */
class RefRegionCylinder : public RefRegion
{
public:
    /// Constructor
    SPATIAL_DOMAINS_EXPORT RefRegionCylinder(
        const unsigned int coordim, NekDouble radius,
        std::vector<NekDouble> coord1, std::vector<NekDouble> coord2,
        std::vector<unsigned int> numModes,
        std::vector<unsigned int> numPoints);
    /// Destructor
    SPATIAL_DOMAINS_EXPORT ~RefRegionCylinder() override = default;

protected:
    /// Check if vertex is inside the surface region
    SPATIAL_DOMAINS_EXPORT bool v_Contains(
        const Array<OneD, NekDouble> &coords) override;
};

} // namespace Nektar::SpatialDomains

#endif // NEKTAR_SPATIALDOMAINS_REFREGIONCYLINDER_H
