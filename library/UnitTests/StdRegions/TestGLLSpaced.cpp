///////////////////////////////////////////////////////////////////////////////
//
// File: TestGLLSpaced.cpp
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
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdNodalTetExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar::QuadEquiSpaced
{

BOOST_AUTO_TEST_CASE(TestTriExpInterpGLLSpaced)
{
    using namespace LibUtilities;

    PointsType PointsTypeDir1 = eGaussLobattoLegendre;
    PointsType PointsTypeDir2 = eGaussRadauMAlpha1Beta0;
    BasisType basisTypeDir1   = eModified_A;
    BasisType basisTypeDir2   = eModified_B;

    unsigned int numPoints = 8;
    unsigned int numEQ     = 6;

    const PointsKey PointsKeyDir1(numPoints, PointsTypeDir1);
    const BasisKey basisKeyDir1(basisTypeDir1, numEQ, PointsKeyDir1);
    const PointsKey PointsKeyDir2(numPoints - 1, PointsTypeDir2);
    const BasisKey basisKeyDir2(basisTypeDir2, numEQ, PointsKeyDir2);

    StdRegions::StdTriExpSharedPtr Exp =
        MemoryManager<Nektar::StdRegions::StdTriExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2);

    // define an equispaced Tri points
    PointsType PointsTypeEq = eNodalTriElec;

    StdRegions::StdTriExpSharedPtr ExpEq =
        MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2, PointsTypeEq);

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

    double epsilon = 1.0e-8;

    Array<OneD, int> sorted;
    LibUtilities::NodalUtilTriangle::CartesianOrdering(numEQ, sorted);

    // compare x-coordinates
    Exp->PhysInterpToGLL(c0, out, numEQ);
    ExpEq->FwdTrans(ceq0, outeq);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outeq[sorted[i]], epsilon);
    }

    // compare y-coordinates
    Exp->PhysInterpToGLL(c1, out, numEQ);
    ExpEq->FwdTrans(ceq1, outeq);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outeq[sorted[i]], epsilon);
    }
}

BOOST_AUTO_TEST_CASE(TestTetpInterpPhysToGLLSpaced)
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
    PointsType PointsTypeEq = eNodalTetElec;

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
    Exp->PhysInterpToGLL(c0, out, numEQ);
    ExpEq->FwdTrans(ceq0, outeq);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outeq[sorted[i]], epsilon);
    }

    // compare y-coordinates
    Exp->PhysInterpToGLL(c1, out, numEQ);
    ExpEq->FwdTrans(ceq1, outeq);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outeq[sorted[i]], epsilon);
    }

    // compare z-coordinates
    Exp->PhysInterpToGLL(c2, out, numEQ);
    ExpEq->FwdTrans(ceq2, outeq);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outeq[sorted[i]], epsilon);
    }
}

BOOST_AUTO_TEST_CASE(TestPrismpInterpPhysToGLLSpaced)
{
    using namespace LibUtilities;

    PointsType PointsTypeDir1 = eGaussLobattoLegendre;
    PointsType PointsTypeDir2 = eGaussLobattoLegendre;
    PointsType PointsTypeDir3 = eGaussRadauMAlpha1Beta0;
    BasisType basisTypeDir1   = eModified_A;
    BasisType basisTypeDir2   = eModified_A;
    BasisType basisTypeDir3   = eModified_B;

    unsigned int numPoints = 5;
    unsigned int numPt     = 4;

    // Set up standard element.
    const PointsKey PointsKeyDir1(numPoints, PointsTypeDir1);
    const PointsKey PointsKeyDir2(numPoints, PointsTypeDir2);
    const PointsKey PointsKeyDir3(numPoints - 1, PointsTypeDir3);
    const BasisKey basisKeyDir1(basisTypeDir1, numPt, PointsKeyDir1);
    const BasisKey basisKeyDir2(basisTypeDir2, numPt, PointsKeyDir2);
    const BasisKey basisKeyDir3(basisTypeDir3, numPt, PointsKeyDir3);

    StdRegions::StdPrismExpSharedPtr Exp =
        MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2, basisKeyDir3);

    // define Prism points
    PointsType PointsType = eNodalPrismElec;

    StdRegions::StdPrismExpSharedPtr ExpGll =
        MemoryManager<StdRegions::StdNodalPrismExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir2, basisKeyDir3, PointsType);

    // Get coordinates at quadrature points
    Array<OneD, NekDouble> c0 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c1 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> c2 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Exp->GetCoords(c0, c1, c2);

    // Get coordinates at gll points
    Array<OneD, NekDouble> cgll0 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> cgll1 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    Array<OneD, NekDouble> cgll2 = Array<OneD, NekDouble>(Exp->GetTotPoints());
    ExpGll->GetCoords(cgll0, cgll1, cgll2);

    Array<OneD, NekDouble> out    = Array<OneD, NekDouble>(Exp->GetNcoeffs());
    Array<OneD, NekDouble> outgll = Array<OneD, NekDouble>(Exp->GetNcoeffs());
    double epsilon                = 1.0e-8;

    Array<OneD, int> sorted;
    LibUtilities::NodalUtilPrism::CartesianOrdering(numPt, sorted);

    // compare x-coordinates
    Exp->PhysInterpToGLL(c0, out, numPt);
    ExpGll->FwdTrans(cgll0, outgll);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outgll[sorted[i]], epsilon);
    }

    // compare y-coordinates
    Exp->PhysInterpToGLL(c1, out, numPt);
    ExpGll->FwdTrans(cgll1, outgll);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outgll[sorted[i]], epsilon);
    }

    // compare z-coordinates
    Exp->PhysInterpToGLL(c2, out, numPt);
    ExpGll->FwdTrans(cgll2, outgll);
    for (int i = 0; i < Exp->GetNcoeffs(); ++i)
    {
        BOOST_CHECK_CLOSE(out[i], outgll[sorted[i]], epsilon);
    }
}

} // namespace Nektar::QuadEquiSpaced
