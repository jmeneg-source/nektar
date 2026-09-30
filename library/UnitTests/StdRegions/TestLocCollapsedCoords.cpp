///////////////////////////////////////////////////////////////////////////////
//
// File: TestLocCollapsedCoords.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Test Cartesian <-> collapsed coordinates
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SpatialDomains/MeshGraph.h>
#include <StdRegions/StdTetExp.h>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar::StdRegionsTests
{

BOOST_AUTO_TEST_CASE(TestLocCollapsedCoordsTet)
{
    unsigned int nq = 6;
    auto pkey0 =
        LibUtilities::PointsKey(nq, LibUtilities::eGaussLobattoLegendre);
    auto pkey1 =
        LibUtilities::PointsKey(nq - 1, LibUtilities::eGaussRadauMAlpha1Beta0);
    auto pkey2 =
        LibUtilities::PointsKey(nq - 1, LibUtilities::eGaussRadauMAlpha2Beta0);

    unsigned int nm = 4;
    auto bkey0 = LibUtilities::BasisKey(LibUtilities::eModified_A, nm, pkey0);
    auto bkey1 = LibUtilities::BasisKey(LibUtilities::eModified_B, nm, pkey1);
    auto bkey2 = LibUtilities::BasisKey(LibUtilities::eModified_C, nm, pkey2);

    // Create tet expansion
    auto stdTet = StdRegions::StdTetExp(bkey0, bkey1, bkey2);

    // Get some evenly distributed points
    auto pkey =
        LibUtilities::PointsKey(nq, LibUtilities::eNodalTetEvenlySpaced);
    auto evenlySpaced = LibUtilities::PointsManager()[pkey];

    int np = evenlySpaced->GetTotNumPoints();
    Array<OneD, NekDouble> xc(np), yc(np), zc(np);
    evenlySpaced->GetPoints(xc, yc, zc);

    double epsilon = 1.0e-8;

    // Loop over points and ensure that the xi -> eta -> xi mapping returns the
    // original coordinates in the tetrahedron.
    for (auto i = 0; i < np; ++i)
    {
        Array<OneD, NekDouble> xi(3), eta(3), tmp(3);
        xi[0] = xc[i];
        xi[1] = yc[i];
        xi[2] = zc[i];

        stdTet.LocCoordToLocCollapsed(xi, eta);
        stdTet.LocCollapsedToLocCoord(eta, tmp);

        BOOST_CHECK_CLOSE(tmp[0], xi[0], epsilon);
        BOOST_CHECK_CLOSE(tmp[1], xi[1], epsilon);
        BOOST_CHECK_CLOSE(tmp[2], xi[2], epsilon);
    }
}
} // namespace Nektar::StdRegionsTests
