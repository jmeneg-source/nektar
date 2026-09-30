///////////////////////////////////////////////////////////////////////////////
//
// File: BLPoints.cpp
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
// Description: 1D boundary layer points
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/BLPoints.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar::LibUtilities
{
bool BLPoints::initPointsManager[] = {
    PointsManager().RegisterCreator(PointsKey(0, eBoundaryLayerPoints),
                                    BLPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eBoundaryLayerPointsRev),
                                    BLPoints::Create)};

void BLPoints::v_CalculatePoints()
{
    // Allocate the storage for points.
    PointsBaseType::v_CalculatePoints();
    size_t npts = m_pointsKey.GetNumPoints();

    // Derived power coefficient.
    NekDouble r = m_pointsKey.GetFactor();
    ASSERTL0(r != NekConstants::kNekUnsetDouble,
             "Must set factor in BLPoints key");

    if (fabs(r - 1.0) < 1e-6)
    {
        NekDouble tmp = 2.0 / (npts - 1.0);
        for (size_t i = 0; i < npts; ++i)
        {
            m_points[0][i] = -1.0 + i * tmp;
        }
    }
    else
    {
        NekDouble a    = 2.0 * (1.0 - r) / (1.0 - pow(r, (double)(npts - 1)));
        m_points[0][0] = -1.0;

        for (size_t i = 1; i < npts; ++i)
        {
            m_points[0][i] = m_points[0][i - 1] + a * pow(r, (double)(i - 1));
        }

        m_points[0][npts - 1] = 1.0;
    }

    if (m_pointsKey.GetPointsType() == eBoundaryLayerPointsRev)
    {
        std::vector<NekDouble> tmp(npts);
        for (size_t i = 0; i < npts; ++i)
        {
            tmp[i] = -m_points[0][npts - 1 - i];
        }

        for (size_t i = 0; i < npts; ++i)
        {
            m_points[0][i] = tmp[i];
        }
    }
}

void BLPoints::v_CalculateWeights()
{
}

void BLPoints::v_CalculateDerivMatrix()
{
    //// Allocate the derivative matrix.
    Points<NekDouble>::v_CalculateDerivMatrix();
}

std::shared_ptr<Points<NekDouble>> BLPoints::Create(const PointsKey &key)
{
    std::shared_ptr<Points<NekDouble>> returnval(
        MemoryManager<BLPoints>::AllocateSharedPtr(key));

    returnval->Initialize();

    return returnval;
}

std::shared_ptr<NekMatrix<NekDouble>> BLPoints::CreateMatrix(
    [[maybe_unused]] const PointsKey &pkey)
{
    ASSERTL0(false, "CreateMatrix not available for Boundary Layer Points");

    return nullptr;
}

const std::shared_ptr<NekMatrix<NekDouble>> BLPoints::v_GetI(
    [[maybe_unused]] const PointsKey &pkey)
{
    ASSERTL0(false, "Interp not available for Boundary Layer Points");

    return nullptr;
}

const std::shared_ptr<NekMatrix<NekDouble>> BLPoints::v_GetI(
    [[maybe_unused]] const Array<OneD, const NekDouble> &x)
{
    ASSERTL0(false, "Interp not available for Boundary Layer Points");

    return nullptr;
}

const std::shared_ptr<NekMatrix<NekDouble>> BLPoints::v_GetI(
    [[maybe_unused]] size_t numpoints,
    [[maybe_unused]] const Array<OneD, const NekDouble> &x)
{
    ASSERTL0(false, "Interp not available for Boundary Layer Points");

    return nullptr;
}

void BLPoints::CalculateInterpMatrix(
    [[maybe_unused]] size_t npts,
    [[maybe_unused]] const Array<OneD, const NekDouble> &xpoints,
    [[maybe_unused]] Array<OneD, NekDouble> &interp)
{
    ASSERTL0(false,
             "CalculateInterpMatrix not available for Boundary Layer Points");
}
} // namespace Nektar::LibUtilities
