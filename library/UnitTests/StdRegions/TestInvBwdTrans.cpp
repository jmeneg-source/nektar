///////////////////////////////////////////////////////////////////////////////
//
// File: TestInvBwdTrans.cpp
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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// Description: Unit tests for inverse backward transform matrices.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdQuadExp.h>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>

namespace Nektar::StdRegionsUnitTests
{

void CheckInvBwdTrans(const StdRegions::StdExpansionSharedPtr &exp)
{
    const StdRegions::StdMatrixKey bwdKey(StdRegions::eBwdTrans,
                                          exp->DetShapeType(), *exp);
    const StdRegions::StdMatrixKey invKey(StdRegions::eInvBwdTrans,
                                          exp->DetShapeType(), *exp);
    const StdRegions::StdMatrixKey invIProdKey(StdRegions::eInvIProductWRTBase,
                                               exp->DetShapeType(), *exp);

    DNekMatSharedPtr bwd      = exp->GetStdMatrix(bwdKey);
    DNekMatSharedPtr invBwd   = exp->GetStdMatrix(invKey);
    DNekMatSharedPtr invIProd = exp->GetStdMatrix(invIProdKey);

    BOOST_REQUIRE_EQUAL(bwd->GetRows(), exp->GetTotPoints());
    BOOST_REQUIRE_EQUAL(bwd->GetColumns(), exp->GetNcoeffs());
    BOOST_REQUIRE_EQUAL(invBwd->GetRows(), exp->GetNcoeffs());
    BOOST_REQUIRE_EQUAL(invBwd->GetColumns(), exp->GetTotPoints());
    BOOST_REQUIRE_EQUAL(invIProd->GetRows(), exp->GetTotPoints());
    BOOST_REQUIRE_EQUAL(invIProd->GetColumns(), exp->GetNcoeffs());

    const NekDouble epsilon = 1.0e-8;
    DNekMat identity        = (*invBwd) * (*bwd);

    for (unsigned int i = 0; i < identity.GetRows(); ++i)
    {
        for (unsigned int j = 0; j < identity.GetColumns(); ++j)
        {
            const NekDouble expected = i == j ? 1.0 : 0.0;
            BOOST_CHECK_SMALL(std::abs(*identity(i, j) - expected), epsilon);
        }
    }

    for (unsigned int i = 0; i < invBwd->GetRows(); ++i)
    {
        for (unsigned int j = 0; j < invBwd->GetColumns(); ++j)
        {
            BOOST_CHECK_SMALL(std::abs(*(*invIProd)(j, i) - *(*invBwd)(i, j)),
                              epsilon);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestQuadInvBwdTransPseudoInverse)
{
    using namespace LibUtilities;

    const PointsKey pointsKey(6, eGaussLobattoLegendre);
    const BasisKey basisKey(eModified_A, 4, pointsKey);

    StdRegions::StdExpansionSharedPtr exp =
        MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(basisKey,
                                                                 basisKey);
    CheckInvBwdTrans(exp);
}

BOOST_AUTO_TEST_CASE(TestHexInvBwdTransPseudoInverse)
{
    using namespace LibUtilities;

    const PointsKey pointsKey(5, eGaussLobattoLegendre);
    const BasisKey basisKey(eModified_A, 3, pointsKey);

    StdRegions::StdExpansionSharedPtr exp =
        MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(
            basisKey, basisKey, basisKey);
    CheckInvBwdTrans(exp);
}
} // namespace Nektar::StdRegionsUnitTests
