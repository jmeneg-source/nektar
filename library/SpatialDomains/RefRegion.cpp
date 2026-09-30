////////////////////////////////////////////////////////////////////////////////
//
//  File: RefRegion.cpp
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

#include "RefRegion.h"
#include <SpatialDomains/RefRegion.h>
#include <SpatialDomains/RefRegionCylinder.h>
#include <SpatialDomains/RefRegionLine.h>
#include <SpatialDomains/RefRegionParallelogram.h>
#include <SpatialDomains/RefRegionSphere.h>

namespace Nektar::SpatialDomains
{

RefRegion::RefRegion(const unsigned int coordim, NekDouble radius,
                     std::vector<NekDouble> coord1,
                     std::vector<NekDouble> coord2,
                     std::vector<unsigned int> numModes,
                     std::vector<unsigned int> numPoints)
    : m_coord1(coord1), m_coord2(coord2), m_numModes(numModes),
      m_numPoints(numPoints), m_radius(radius), m_coordim(coordim)
{
}

} // namespace Nektar::SpatialDomains
