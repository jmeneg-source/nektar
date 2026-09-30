///////////////////////////////////////////////////////////////////////////////
//
// File: TestEquiSpaced.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/NodalUtil.h>
#include <StdRegions/StdNodalTetExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdTriExp.h>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar::QuadEquiSpaced
{

BOOST_AUTO_TEST_CASE(TestQuadExpInterpPhysToEquiSpaced)
{

    LibUtilities::PointsType PointsTypeDir1 =
        Nektar::LibUtilities::eGaussLobattoLegendre;
    LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
    unsigned int numPoints                = 6;
    const Nektar::LibUtilities::PointsKey PointsKeyDir1(numPoints,
                                                        PointsTypeDir1);
    const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1, 4,
                                                      PointsKeyDir1);

    StdRegions::StdQuadExpSharedPtr Exp =
        MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

    unsigned int numEQ         = 4;
    Array<OneD, NekDouble> c0  = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c1  = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c2  = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> out = Array<OneD, NekDouble>(numEQ * numEQ);

    double epsilon = 1.0e-8;

    Exp->GetCoords(c0, c1, c2);

    Exp->PhysInterpToSimplexEquiSpaced(c0, out, numEQ);
    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ; ++i)
        {
            BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * i / (numEQ - 1.0),
                              epsilon);
        }
    }

    Exp->PhysInterpToSimplexEquiSpaced(c1, out, numEQ);
    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ; ++i)
        {
            BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * j / (numEQ - 1.0),
                              epsilon);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestTriExpInterpPhysToEquiSpaced)
{
    using namespace LibUtilities;

    PointsType PointsTypeDir1 = eGaussLobattoLegendre;
    PointsType PointsTypeDir2 = eGaussRadauMAlpha1Beta0;
    BasisType basisTypeDir1   = eModified_A;
    BasisType basisTypeDir2   = eModified_B;

    unsigned int numPoints = 6;
    unsigned int numEQ     = 4;

    const PointsKey PointsKeyDir1(numPoints, PointsTypeDir1);
    const BasisKey basisKeyDir1(basisTypeDir1, numEQ, PointsKeyDir1);
    const PointsKey PointsKeyDir2(numPoints - 1, PointsTypeDir2);
    const BasisKey basisKeyDir2(basisTypeDir2, numEQ, PointsKeyDir2);

    StdRegions::StdTriExpSharedPtr Exp =
        MemoryManager<Nektar::StdRegions::StdTriExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2);

    Array<OneD, NekDouble> c0  = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c1  = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c2  = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> out = Array<OneD, NekDouble>(numEQ * numEQ);

    double epsilon = 1.0e-8;

    Exp->GetCoords(c0, c1, c2);

    Exp->PhysInterpToSimplexEquiSpaced(c0, out, numEQ);
    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ - j; ++i)
        {
            BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * i / (numEQ - 1.0),
                              epsilon);
        }
    }

    Exp->PhysInterpToSimplexEquiSpaced(c1, out, numEQ);
    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ - j; ++i)
        {
            BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * j / (numEQ - 1.0),
                              epsilon);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestTetExpInterpPhysToEquiSpaced)
{
    using namespace LibUtilities;

    PointsType PointsTypeDir1 = eGaussLobattoLegendre;
    PointsType PointsTypeDir2 = eGaussRadauMAlpha1Beta0;
    PointsType PointsTypeDir3 = eGaussRadauMAlpha2Beta0;
    BasisType basisTypeDir1   = eModified_A;
    BasisType basisTypeDir2   = eModified_B;
    BasisType basisTypeDir3   = eModified_C;

    unsigned int numPoints = 10;
    unsigned int numEQ     = 8;

    // Set up standard element.
    const PointsKey PointsKeyDir1(numPoints, PointsTypeDir1);
    const PointsKey PointsKeyDir2(numPoints - 1, PointsTypeDir2);
    const PointsKey PointsKeyDir3(numPoints - 1, PointsTypeDir3);
    const BasisKey basisKeyDir2(basisTypeDir2, numEQ, PointsKeyDir2);
    const BasisKey basisKeyDir1(basisTypeDir1, numEQ, PointsKeyDir1);
    const BasisKey basisKeyDir3(basisTypeDir3, numEQ, PointsKeyDir3);

    StdRegions::StdTetExpSharedPtr Exp =
        MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2, basisKeyDir3);

    // define an equispaced Tet points
    PointsType PointsTypeEq = eNodalTetEvenlySpaced;

    StdRegions::StdTetExpSharedPtr ExpEq =
        MemoryManager<StdRegions::StdNodalTetExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2, basisKeyDir3, PointsTypeEq);

    // Get coordinates at quadrature points
    Array<OneD, NekDouble> c0 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c1 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c2 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Exp->GetCoords(c0, c1, c2);

    // Get coordinates at equispaced points
    Array<OneD, NekDouble> ceq0 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> ceq1 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> ceq2 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    ExpEq->GetCoords(ceq0, ceq1, ceq2);

    Array<OneD, NekDouble> out   = Array<OneD, NekDouble>(Exp->GetNcoeffs());
    Array<OneD, NekDouble> outeq = Array<OneD, NekDouble>(Exp->GetNcoeffs());
    double epsilon               = 1.0e-8;

    Array<OneD, int> sorted;
    LibUtilities::NodalUtilTetrahedron::CartesianOrdering(numEQ, sorted);

    // compare x-coordinates
    Exp->PhysInterpToSimplexEquiSpaced(c0, out, numEQ);
    int cnt = 0;
    for (int k = 0; k < numEQ; ++k)
    {
        for (int j = 0; j < numEQ - k; ++j)
        {
            for (int i = 0; i < numEQ - k - j; ++i)
            {
                BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * i / (numEQ - 1.0),
                                  epsilon);
            }
        }
    }

    // compare y-coordinates
    Exp->PhysInterpToSimplexEquiSpaced(c1, out, numEQ);
    cnt = 0;
    for (int k = 0; k < numEQ; ++k)
    {
        for (int j = 0; j < numEQ - k; ++j)
        {
            for (int i = 0; i < numEQ - k - j; ++i)
            {
                BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * j / (numEQ - 1.0),
                                  epsilon);
            }
        }
    }

    // compare z-coordinates
    Exp->PhysInterpToSimplexEquiSpaced(c2, out, numEQ);
    cnt = 0;
    for (int k = 0; k < numEQ; ++k)
    {
        for (int j = 0; j < numEQ - k; ++j)
        {
            for (int i = 0; i < numEQ - k - j; ++i)
            {
                BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * k / (numEQ - 1.0),
                                  epsilon);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(TestPrismExpInterpPhysToEquiSpaced)
{
    using namespace LibUtilities;

    PointsType PointsTypeDir1 = eGaussLobattoLegendre;
    PointsType PointsTypeDir2 = eGaussLobattoLegendre;
    PointsType PointsTypeDir3 = eGaussRadauMAlpha1Beta0;
    BasisType basisTypeDir1   = eModified_A;
    BasisType basisTypeDir2   = eModified_A;
    BasisType basisTypeDir3   = eModified_B;

    unsigned int numPoints = 10;
    unsigned int numEQ     = 8;

    // Set up standard element.
    const PointsKey PointsKeyDir1(numPoints, PointsTypeDir1);
    const PointsKey PointsKeyDir2(numPoints, PointsTypeDir2);
    const PointsKey PointsKeyDir3(numPoints - 1, PointsTypeDir3);
    const BasisKey basisKeyDir1(basisTypeDir1, numEQ, PointsKeyDir1);
    const BasisKey basisKeyDir2(basisTypeDir2, numEQ, PointsKeyDir2);
    const BasisKey basisKeyDir3(basisTypeDir3, numEQ, PointsKeyDir3);

    StdRegions::StdPrismExpSharedPtr Exp =
        MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2, basisKeyDir3);

    // define an equispaced Tet points
    PointsType PointsTypeEq = eNodalPrismEvenlySpaced;

    StdRegions::StdTetExpSharedPtr ExpEq =
        MemoryManager<StdRegions::StdNodalTetExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2, basisKeyDir3, PointsTypeEq);

    // Get coordinates at quadrature points
    Array<OneD, NekDouble> c0 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c1 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c2 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Exp->GetCoords(c0, c1, c2);

    // Get coordinates at equispaced points
    Array<OneD, NekDouble> ceq0 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> ceq1 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> ceq2 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    ExpEq->GetCoords(ceq0, ceq1, ceq2);

    Array<OneD, NekDouble> out   = Array<OneD, NekDouble>(Exp->GetNcoeffs());
    Array<OneD, NekDouble> outeq = Array<OneD, NekDouble>(Exp->GetNcoeffs());
    double epsilon               = 1.0e-8;

    Array<OneD, int> sorted;
    LibUtilities::NodalUtilPrism::CartesianOrdering(numEQ, sorted);

    // compare x-coordinates
    Exp->PhysInterpToSimplexEquiSpaced(c0, out, numEQ);
    int cnt = 0;
    for (int k = 0; k < numEQ; ++k)
    {
        for (int j = 0; j < numEQ; ++j)
        {
            for (int i = 0; i < numEQ - k; ++i)
            {
                BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * i / (numEQ - 1.0),
                                  epsilon);
            }
        }
    }

    // compare y-coordinates
    Exp->PhysInterpToSimplexEquiSpaced(c1, out, numEQ);
    cnt = 0;
    for (int k = 0; k < numEQ; ++k)
    {
        for (int j = 0; j < numEQ; ++j)
        {
            for (int i = 0; i < numEQ - k; ++i)
            {
                BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * j / (numEQ - 1.0),
                                  epsilon);
            }
        }
    }

    // compare z-coordinates
    Exp->PhysInterpToSimplexEquiSpaced(c2, out, numEQ);
    cnt = 0;
    for (int k = 0; k < numEQ; ++k)
    {
        for (int j = 0; j < numEQ; ++j)
        {
            for (int i = 0; i < numEQ - k; ++i)
            {
                BOOST_CHECK_CLOSE(out[cnt++], -1.0 + 2.0 * k / (numEQ - 1.0),
                                  epsilon);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(TestQuadExpEquiSpacedToPhys)
{
    LibUtilities::PointsType PointsTypeDir1 =
        Nektar::LibUtilities::eGaussLobattoLegendre;
    LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
    unsigned int numPoints                = 10;
    unsigned int numModes                 = 8;
    const Nektar::LibUtilities::PointsKey PointsKeyDir1(numPoints,
                                                        PointsTypeDir1);
    const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1, numModes,
                                                      PointsKeyDir1);

    StdRegions::StdQuadExpSharedPtr Exp =
        MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

    unsigned int numEQ          = 4;
    Array<OneD, NekDouble> c0   = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c1   = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c2   = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> phys = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> in   = Array<OneD, NekDouble>(numEQ * numEQ);

    double epsilon = 1.0e-8;

    Exp->GetCoords(c0, c1, c2);

    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ; ++i)
        {
            in[cnt++] = -1.0 + 2.0 * i / (numEQ - 1.0);
        }
    }

    Exp->EquiSpacedToPhys(numEQ, in, phys);
    for (int j = 0, cnt = 0; j < numPoints; ++j)
    {
        for (int i = 0; i < numPoints; ++i)
        {
            BOOST_CHECK_CLOSE(phys[cnt], c0[cnt], epsilon);
            cnt++;
        }
    }

    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ; ++i)
        {
            in[cnt++] = -1.0 + 2.0 * j / (numEQ - 1.0);
        }
    }
    Exp->EquiSpacedToPhys(numEQ, in, phys);
    for (int j = 0, cnt = 0; j < numPoints; ++j)
    {
        for (int i = 0; i < numPoints; ++i)
        {
            BOOST_CHECK_CLOSE(phys[cnt], c1[cnt], epsilon);
            cnt++;
        }
    }
}

BOOST_AUTO_TEST_CASE(TestTriExpEquiSpacedToPhys)
{

    LibUtilities::PointsType PointsTypeDir1 =
        Nektar::LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType PointsTypeDir2 =
        Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
    Nektar::LibUtilities::BasisType basisTypeDir1 =
        Nektar::LibUtilities::eModified_A;
    Nektar::LibUtilities::BasisType basisTypeDir2 =
        Nektar::LibUtilities::eModified_B;

    unsigned int numPoints = 6;

    const Nektar::LibUtilities::PointsKey PointsKeyDir1(numPoints,
                                                        PointsTypeDir1);
    const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1, 4,
                                                      PointsKeyDir1);
    const Nektar::LibUtilities::PointsKey PointsKeyDir2(numPoints - 1,
                                                        PointsTypeDir2);
    const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir2, 4,
                                                      PointsKeyDir2);

    StdRegions::StdTriExpSharedPtr Exp =
        MemoryManager<Nektar::StdRegions::StdTriExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2);

    unsigned int numEQ          = 4;
    Array<OneD, NekDouble> c0   = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c1   = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c2   = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> phys = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> in   = Array<OneD, NekDouble>(numEQ * numEQ);

    double epsilon = 1.0e-8;

    Exp->GetCoords(c0, c1, c2);

    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ - j; ++i)
        {
            in[cnt++] = -1.0 + 2.0 * i / (numEQ - 1.0);
        }
    }
    Exp->EquiSpacedToPhys(numEQ, in, phys);
    for (int j = 0, cnt = 0; j < numPoints - 1; ++j)
    {
        for (int i = 0; i < numPoints; ++i)
        {
            BOOST_CHECK_CLOSE(phys[cnt], c0[cnt], epsilon);
            cnt++;
        }
    }

    for (int j = 0, cnt = 0; j < numEQ; ++j)
    {
        for (int i = 0; i < numEQ - j; ++i)
        {
            in[cnt++] = -1.0 + 2.0 * j / (numEQ - 1.0);
        }
    }
    Exp->EquiSpacedToPhys(numEQ, in, phys);
    for (int j = 0, cnt = 0; j < numPoints - 1; ++j)
    {
        for (int i = 0; i < numPoints; ++i)
        {
            BOOST_CHECK_CLOSE(phys[cnt], c1[cnt], epsilon);
            cnt++;
        }
    }
}

} // namespace Nektar::QuadEquiSpaced
