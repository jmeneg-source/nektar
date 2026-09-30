///////////////////////////////////////////////////////////////////////////////
//
// File: TestInterpCoeff.cpp
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
// Description: Unit tests for coefficient interpolation.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>

namespace Nektar::UnitTests
{

Array<OneD, NekDouble> TensorProductInterp2D(
    const LibUtilities::BasisKey &fbasis0,
    const LibUtilities::BasisKey &fbasis1,
    const Array<OneD, const NekDouble> &from,
    const LibUtilities::BasisKey &tbasis0,
    const LibUtilities::BasisKey &tbasis1)
{
    const size_t fnm0 = fbasis0.GetNumModes();
    const size_t fnm1 = fbasis1.GetNumModes();
    const size_t tnm0 = tbasis0.GetNumModes();
    const size_t tnm1 = tbasis1.GetNumModes();

    DNekMatSharedPtr ft0 = LibUtilities::BasisManager()[fbasis0]->GetI(tbasis0);
    DNekMatSharedPtr ft1 = LibUtilities::BasisManager()[fbasis1]->GetI(tbasis1);

    Array<OneD, NekDouble> expected(tnm0 * tnm1, 0.0);
    for (size_t i1 = 0; i1 < tnm1; ++i1)
    {
        for (size_t i0 = 0; i0 < tnm0; ++i0)
        {
            for (size_t j1 = 0; j1 < fnm1; ++j1)
            {
                for (size_t j0 = 0; j0 < fnm0; ++j0)
                {
                    expected[i0 + tnm0 * i1] +=
                        (*ft0)(i0, j0) * from[j0 + fnm0 * j1] * (*ft1)(i1, j1);
                }
            }
        }
    }

    return expected;
}

Array<OneD, NekDouble> TensorProductInterp3D(
    const LibUtilities::BasisKey &fbasis0,
    const LibUtilities::BasisKey &fbasis1,
    const LibUtilities::BasisKey &fbasis2,
    const Array<OneD, const NekDouble> &from,
    const LibUtilities::BasisKey &tbasis0,
    const LibUtilities::BasisKey &tbasis1,
    const LibUtilities::BasisKey &tbasis2)
{
    const size_t fnm0 = fbasis0.GetNumModes();
    const size_t fnm1 = fbasis1.GetNumModes();
    const size_t fnm2 = fbasis2.GetNumModes();
    const size_t tnm0 = tbasis0.GetNumModes();
    const size_t tnm1 = tbasis1.GetNumModes();
    const size_t tnm2 = tbasis2.GetNumModes();

    DNekMatSharedPtr ft0 = LibUtilities::BasisManager()[fbasis0]->GetI(tbasis0);
    DNekMatSharedPtr ft1 = LibUtilities::BasisManager()[fbasis1]->GetI(tbasis1);
    DNekMatSharedPtr ft2 = LibUtilities::BasisManager()[fbasis2]->GetI(tbasis2);

    Array<OneD, NekDouble> expected(tnm0 * tnm1 * tnm2, 0.0);
    for (size_t i2 = 0; i2 < tnm2; ++i2)
    {
        for (size_t i1 = 0; i1 < tnm1; ++i1)
        {
            for (size_t i0 = 0; i0 < tnm0; ++i0)
            {
                for (size_t j2 = 0; j2 < fnm2; ++j2)
                {
                    for (size_t j1 = 0; j1 < fnm1; ++j1)
                    {
                        for (size_t j0 = 0; j0 < fnm0; ++j0)
                        {
                            expected[i0 + tnm0 * (i1 + tnm1 * i2)] +=
                                (*ft0)(i0, j0) * (*ft1)(i1, j1) *
                                (*ft2)(i2, j2) *
                                from[j0 + fnm0 * (j1 + fnm1 * j2)];
                        }
                    }
                }
            }
        }
    }

    return expected;
}

BOOST_AUTO_TEST_CASE(TestInterpCoeff2DDifferentOrdersAndBases)
{
    using namespace LibUtilities;

    const BasisKey fbasis0(eModified_A, 3, PointsKey(5, eGaussLobattoLegendre));
    const BasisKey fbasis1(eOrtho_A, 4, PointsKey(6, eGaussLobattoLegendre));
    const BasisKey tbasis0(eOrtho_A, 5, PointsKey(7, eGaussLobattoLegendre));
    const BasisKey tbasis1(eModified_A, 6, PointsKey(8, eGaussLobattoLegendre));

    Array<OneD, NekDouble> coeffs(fbasis0.GetNumModes() *
                                  fbasis1.GetNumModes());
    for (size_t i = 0; i < coeffs.size(); ++i)
    {
        coeffs[i] = 0.125 + 0.25 * static_cast<NekDouble>(i);
    }

    Array<OneD, NekDouble> actual(tbasis0.GetNumModes() *
                                  tbasis1.GetNumModes());
    Array<OneD, NekDouble> expected =
        TensorProductInterp2D(fbasis0, fbasis1, coeffs, tbasis0, tbasis1);

    LibUtilities::InterpCoeff2D(fbasis0, fbasis1, coeffs, tbasis0, tbasis1,
                                actual);

    BOOST_REQUIRE_EQUAL(expected.size(), actual.size());

    const NekDouble epsilon = 1.0e-10;
    for (size_t i = 0; i < expected.size(); ++i)
    {
        BOOST_CHECK_SMALL(std::abs(expected[i] - actual[i]), epsilon);
    }
}

BOOST_AUTO_TEST_CASE(TestInterpCoeff3DDifferentOrdersAndBases)
{
    using namespace LibUtilities;

    const BasisKey fbasis0(eModified_A, 2, PointsKey(4, eGaussLobattoLegendre));
    const BasisKey fbasis1(eOrtho_A, 3, PointsKey(5, eGaussLobattoLegendre));
    const BasisKey fbasis2(eModified_A, 4, PointsKey(6, eGaussLobattoLegendre));
    const BasisKey tbasis0(eOrtho_A, 4, PointsKey(6, eGaussLobattoLegendre));
    const BasisKey tbasis1(eModified_A, 5, PointsKey(7, eGaussLobattoLegendre));
    const BasisKey tbasis2(eOrtho_A, 6, PointsKey(8, eGaussLobattoLegendre));

    Array<OneD, NekDouble> coeffs(
        fbasis0.GetNumModes() * fbasis1.GetNumModes() * fbasis2.GetNumModes());
    for (size_t i = 0; i < coeffs.size(); ++i)
    {
        coeffs[i] = 0.2 + 0.1 * static_cast<NekDouble>((3 * i) % 11);
    }

    Array<OneD, NekDouble> actual(
        tbasis0.GetNumModes() * tbasis1.GetNumModes() * tbasis2.GetNumModes());
    Array<OneD, NekDouble> expected = TensorProductInterp3D(
        fbasis0, fbasis1, fbasis2, coeffs, tbasis0, tbasis1, tbasis2);

    LibUtilities::InterpCoeff3D(fbasis0, fbasis1, fbasis2, coeffs, tbasis0,
                                tbasis1, tbasis2, actual);

    BOOST_REQUIRE_EQUAL(expected.size(), actual.size());

    const NekDouble epsilon = 1.0e-10;
    for (size_t i = 0; i < expected.size(); ++i)
    {
        BOOST_CHECK_SMALL(std::abs(expected[i] - actual[i]), epsilon);
    }
}
} // namespace Nektar::UnitTests
