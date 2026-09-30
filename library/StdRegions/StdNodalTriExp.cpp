///////////////////////////////////////////////////////////////////////////////
//
// File: StdNodalTriExp.cpp
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
// Description: Nodal triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h> // for PointsManager, etc
#include <StdRegions/StdNodalTriExp.h>

using namespace std;

namespace Nektar::StdRegions
{
StdNodalTriExp::StdNodalTriExp(const LibUtilities::BasisKey &Ba,
                               const LibUtilities::BasisKey &Bb,
                               LibUtilities::PointsType Ntype)
    : StdExpansion(LibUtilities::StdTriData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes()),
                   2, Ba, Bb),
      StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes()),
                     Ba, Bb),
      StdTriExp(Ba, Bb), m_nodalPointsKey(Ba.GetNumModes(), Ntype)
{
    ASSERTL0(m_base[0]->GetNumModes() == m_base[1]->GetNumModes(),
             "Nodal basis initiated with different orders in the a "
             "and b directions");

    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());

    StdFacKey w1key(eWeights1, Bb);
    // get weights[1] from manager where points are rescaled
    m_weights.push_back(GetStdFac(w1key));
}

bool StdNodalTriExp::v_IsNodalNonTensorialExp()
{
    return true;
}

//-------------------------------
// Nodal basis specific routines
//-------------------------------

void StdNodalTriExp::v_NodalToModal(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    StdMatrixKey Nkey(eInvNBasisTrans, DetShapeType(), *this,
                      NullConstFactorMap, NullVarCoeffMap, NullVarFactorsMap,
                      m_nodalPointsKey.GetPointsType());
    DNekMatSharedPtr inv_vdm = GetStdMatrix(Nkey);

    NekVector<NekDouble> nodal(m_ncoeffs, inarray, eWrapper);
    NekVector<NekDouble> modal(m_ncoeffs, outarray, eWrapper);
    modal = (*inv_vdm) * nodal;
}

// Operate with transpose of NodalToModal transformation
void StdNodalTriExp::NodalToModalTranspose(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    StdMatrixKey Nkey(eInvNBasisTrans, DetShapeType(), *this,
                      NullConstFactorMap, NullVarCoeffMap, NullVarFactorsMap,
                      m_nodalPointsKey.GetPointsType());
    DNekMatSharedPtr inv_vdm = GetStdMatrix(Nkey);

    NekVector<NekDouble> nodal(m_ncoeffs, inarray, eCopy);
    NekVector<NekDouble> modal(m_ncoeffs, outarray, eWrapper);
    modal = Transpose(*inv_vdm) * nodal;
}

void StdNodalTriExp::ModalToNodal(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &outarray)
{
    StdMatrixKey Nkey(eNBasisTrans, DetShapeType(), *this, NullConstFactorMap,
                      NullVarCoeffMap, NullVarFactorsMap,
                      m_nodalPointsKey.GetPointsType());
    DNekMatSharedPtr vdm = GetStdMatrix(Nkey);

    // Multiply out matrix
    NekVector<NekDouble> modal(m_ncoeffs, inarray, eWrapper);
    NekVector<NekDouble> nodal(m_ncoeffs, outarray, eWrapper);
    nodal = (*vdm) * modal;
}

void StdNodalTriExp::GetNodalPoints(Array<OneD, const NekDouble> &x,
                                    Array<OneD, const NekDouble> &y)
{
    LibUtilities::PointsManager()[m_nodalPointsKey]->GetPoints(x, y);
}

DNekMatSharedPtr StdNodalTriExp::GenNBasisTransMatrix()
{
    int i, j;
    Array<OneD, const NekDouble> r, s;
    Array<OneD, NekDouble> c(2);
    DNekMatSharedPtr Mat;

    Mat = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
    GetNodalPoints(r, s);

    // Store the values of m_phys in a temporary array
    int nqtot = GetTotPoints();
    Array<OneD, NekDouble> phys(nqtot);

    for (i = 0; i < m_ncoeffs; ++i)
    {
        // fill physical space with mode i
        StdTriExp::v_FillMode(i, phys);

        // interpolate mode i to the Nodal points 'j' and
        // store in outarray
        for (j = 0; j < m_ncoeffs; ++j)
        {
            c[0]         = r[j];
            c[1]         = s[j];
            (*Mat)(j, i) = StdExpansion2D::v_PhysEvaluate(c, phys);
        }
    }
    return Mat;
}

//---------------------------------------
// Transforms
//---------------------------------------

void StdNodalTriExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp(m_ncoeffs);
    v_NodalToModal(inarray, tmp);
    StdTriExp::v_BwdTrans(tmp, outarray);
}

//---------------------------------------
// Inner product functions
//---------------------------------------

void StdNodalTriExp::v_IProductWRTBase(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    StdTriExp::v_IProductWRTBase(inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

void StdNodalTriExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    StdTriExp::v_IProductWRTDerivBase(dir, inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

//---------------------------------------
// Evaluation functions
//---------------------------------------

void StdNodalTriExp::v_FillMode(const int mode,
                                Array<OneD, NekDouble> &outarray)
{
    ASSERTL2(mode >= m_ncoeffs,
             "calling argument mode is larger than total expansion order");

    Vmath::Zero(m_ncoeffs, outarray, 1);
    outarray[mode] = 1.0;
    v_BwdTrans(outarray, outarray);
}

//---------------------------
// Helper functions
//---------------------------

LibUtilities::ShapeType StdNodalTriExp::v_DetShapeType() const
{
    return LibUtilities::eNodalTri;
}

int StdNodalTriExp::v_NumBndryCoeffs() const
{
    return 3 + (GetBasisNumModes(0) - 2) + 2 * (GetBasisNumModes(1) - 2);
}

const LibUtilities::BasisKey StdNodalTriExp::v_GetTraceBasisKey(
    const int i, [[maybe_unused]] const int j,
    [[maybe_unused]] bool UseGLL) const
{
    ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

    // Get basiskey (0 or 1) according to edge id i
    int dir = (i != 0);

    switch (m_base[dir]->GetBasisType())
    {
        case LibUtilities::eOrtho_A:
        case LibUtilities::eModified_A:
        {
            switch (m_base[dir]->GetPointsType())
            {
                case LibUtilities::eGaussLobattoLegendre:
                {
                    LibUtilities::PointsKey pkey(
                        m_base[dir]
                            ->GetBasisKey()
                            .GetPointsKey()
                            .GetNumPoints(),
                        LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,
                                                  m_base[dir]->GetNumModes(),
                                                  pkey);
                }
                break;
                default:
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Unexpected points distribution " +
                                 LibUtilities::kPointsTypeStr
                                     [m_base[dir]->GetPointsType()] +
                                 " in StdNodalTriExp::v_GetTraceBasisKey");
                }
            }
        }
        break;
        case LibUtilities::eModified_B:
        case LibUtilities::eOrtho_B:
        {
            switch (m_base[dir]->GetPointsType())
            {
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                {
                    LibUtilities::PointsKey pkey(
                        m_base[dir]
                                ->GetBasisKey()
                                .GetPointsKey()
                                .GetNumPoints() +
                            1,
                        LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,
                                                  m_base[dir]->GetNumModes(),
                                                  pkey);
                }
                break;
                default:
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Unexpected points distribution " +
                                 LibUtilities::kPointsTypeStr
                                     [m_base[dir]->GetPointsType()] +
                                 " in StdNodalTriExp::v_GetTraceBasisKey");
                }
            }
        }
        break;
        default:
        {
            NEKERROR(ErrorUtil::efatal,
                     "Information not available to set edge key");
        }
    }
    return LibUtilities::NullBasisKey;
}

//--------------------------
// Mappings
//--------------------------

int StdNodalTriExp::v_GetVertexMap(const int localVertexId,
                                   [[maybe_unused]] bool useCoeffPacking)
{
    ASSERTL0(localVertexId >= 0 && localVertexId <= 2,
             "Local Vertex ID must be between 0 and 2");
    return localVertexId;
}

void StdNodalTriExp::v_GetTraceToElementMap(const int eid,
                                            Array<OneD, unsigned int> &maparray,
                                            Array<OneD, int> &signarray,
                                            Orientation edgeOrient, int P,
                                            [[maybe_unused]] int Q)
{
    ASSERTL0(eid >= 0 && eid <= 2, "Local Edge ID must be between 0 and 2");

    const int nEdgeCoeffs = GetTraceNcoeffs(eid);

    ASSERTL0(P == -1 || P == nEdgeCoeffs,
             "Nodal triangle not set up to deal with variable "
             "polynomial order.");

    if (maparray.size() != nEdgeCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nEdgeCoeffs);
    }

    if (signarray.size() != nEdgeCoeffs)
    {
        signarray = Array<OneD, int>(nEdgeCoeffs, 1);
    }
    else
    {
        fill(signarray.data(), signarray.data() + nEdgeCoeffs, 1);
    }

    Orientation orient = edgeOrient;
    if (eid == 2)
    {
        orient = orient == eForwards ? eBackwards : eForwards;
    }

    maparray[0]               = eid;
    maparray[nEdgeCoeffs - 1] = eid == 2 ? 0 : eid + 1;
    for (int i = 2; i < nEdgeCoeffs; i++)
    {
        maparray[i - 1] = eid * (nEdgeCoeffs - 2) + 1 + i;
    }

    if (orient == eBackwards)
    {
        reverse(maparray.data(), maparray.data() + nEdgeCoeffs);
    }
}

void StdNodalTriExp::v_GetTraceInteriorToElementMap(
    const int eid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation edgeOrient)
{
    ASSERTL0(eid >= 0 && eid <= 2, "Local Edge ID must be between 0 and 2");

    const int nEdgeIntCoeffs = GetTraceNcoeffs(eid) - 2;

    if (maparray.size() != nEdgeIntCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
    }

    if (signarray.size() != nEdgeIntCoeffs)
    {
        signarray = Array<OneD, int>(nEdgeIntCoeffs, 1);
    }
    else
    {
        fill(signarray.data(), signarray.data() + nEdgeIntCoeffs, 1);
    }

    Orientation orient = edgeOrient;
    if (eid == 2)
    {
        orient = orient == eForwards ? eBackwards : eForwards;
    }

    for (int i = 0; i < nEdgeIntCoeffs; i++)
    {
        maparray[i] = eid * nEdgeIntCoeffs + 3 + i;
    }

    if (orient == eBackwards)
    {
        reverse(maparray.data(), maparray.data() + nEdgeIntCoeffs);
    }
}

void StdNodalTriExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    unsigned int i;
    if (outarray.size() != GetNcoeffs() - NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(GetNcoeffs() - NumBndryCoeffs());
    }

    for (i = NumBndryCoeffs(); i < GetNcoeffs(); i++)
    {
        outarray[i - NumBndryCoeffs()] = i;
    }
}

void StdNodalTriExp::v_GetBoundaryMap(Array<OneD, unsigned int> &outarray)
{
    unsigned int i;
    if (outarray.size() != NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
    }

    for (i = 0; i < NumBndryCoeffs(); i++)
    {
        outarray[i] = i;
    }
}

//---------------------------------------
// Wrapper functions
//---------------------------------------

DNekMatSharedPtr StdNodalTriExp::v_GenMatrix(const StdMatrixKey &mkey)
{
    DNekMatSharedPtr Mat;

    switch (mkey.GetMatrixType())
    {
        case eNBasisTrans:
            Mat = GenNBasisTransMatrix();
            break;
        default:
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
            break;
    }

    return Mat;
}

DNekMatSharedPtr StdNodalTriExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return StdNodalTriExp::v_GenMatrix(mkey);
}

//---------------------------------------
// Operator evaluation functions
//---------------------------------------

void StdNodalTriExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdNodalTriExp::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

void StdNodalTriExp::v_LaplacianMatrixOp(
    const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)

{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void StdNodalTriExp::v_WeakDerivMatrixOp(
    const int i, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void StdNodalTriExp::v_HelmholtzMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

//---------------------------------------
// Private helper functions
//---------------------------------------

} // namespace Nektar::StdRegions
