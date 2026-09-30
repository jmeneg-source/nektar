///////////////////////////////////////////////////////////////////////////////
//
// File: StdNodalTetExp.cpp
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
// Description: Nodal tetrahedral routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h> // for PointsManager, etc
#include <StdRegions/StdNodalTetExp.h>

namespace Nektar::StdRegions
{

StdNodalTetExp::StdNodalTetExp(const LibUtilities::BasisKey &Ba,
                               const LibUtilities::BasisKey &Bb,
                               const LibUtilities::BasisKey &Bc,
                               LibUtilities::PointsType Ntype)
    : StdExpansion(LibUtilities::StdTetData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                   3, Ba, Bb, Bc),
      StdExpansion3D(LibUtilities::StdTetData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                     Ba, Bb, Bc),
      StdTetExp(Ba, Bb, Bc), m_nodalPointsKey(Ba.GetNumModes(), Ntype)
{
    ASSERTL0(Ba.GetNumModes() <= Bb.GetNumModes(),
             "order in 'a' direction is higher than order "
             "in 'b' direction");
    ASSERTL0(Ba.GetNumModes() <= Bc.GetNumModes(),
             "order in 'a' direction is higher than order "
             "in 'c' direction");
    ASSERTL0(Bb.GetNumModes() <= Bc.GetNumModes(),
             "order in 'b' direction is higher than order "
             "in 'c' direction");
}

bool StdNodalTetExp::v_IsNodalNonTensorialExp()
{
    return true;
}

//-------------------------------
// Nodal basis specific routines
//-------------------------------

void StdNodalTetExp::NodalToModal(const Array<OneD, const NekDouble> &inarray,
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
void StdNodalTetExp::NodalToModalTranspose(
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

void StdNodalTetExp::ModalToNodal(const Array<OneD, const NekDouble> &inarray,
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

void StdNodalTetExp::GetNodalPoints(Array<OneD, const NekDouble> &x,
                                    Array<OneD, const NekDouble> &y,
                                    Array<OneD, const NekDouble> &z)
{
    LibUtilities::PointsManager()[m_nodalPointsKey]->GetPoints(x, y, z);
}

DNekMatSharedPtr StdNodalTetExp::GenNBasisTransMatrix()
{
    int i, j;
    Array<OneD, const NekDouble> r, s, t;
    Array<OneD, NekDouble> c(3);
    DNekMatSharedPtr Mat;

    Mat = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
    GetNodalPoints(r, s, t);

    // Store the values of m_phys in a temporary array
    int nqtot = GetTotPoints();
    Array<OneD, NekDouble> tmp_phys(nqtot);

    for (i = 0; i < m_ncoeffs; ++i)
    {
        // fill physical space with mode i
        StdTetExp::v_FillMode(i, tmp_phys);

        // interpolate mode i to the Nodal points 'j' and
        // store in outarray
        for (j = 0; j < m_ncoeffs; ++j)
        {
            c[0]         = r[j];
            c[1]         = s[j];
            c[2]         = t[j];
            (*Mat)(j, i) = StdExpansion3D::v_PhysEvaluate(c, tmp_phys);
        }
    }

    return Mat;
}

//---------------------------------------
// Transforms
//---------------------------------------

void StdNodalTetExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp(m_ncoeffs);
    NodalToModal(inarray, tmp);
    StdTetExp::v_BwdTrans(tmp, outarray);
}

//---------------------------------------
// Inner product functions
//---------------------------------------

void StdNodalTetExp::v_IProductWRTBase(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    StdTetExp::v_IProductWRTBase(inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

void StdNodalTetExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    StdTetExp::v_IProductWRTDerivBase(dir, inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

//---------------------------------------
// Evaluation functions
//---------------------------------------

void StdNodalTetExp::v_FillMode(const int mode,
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

LibUtilities::ShapeType StdNodalTetExp::v_DetShapeType() const
{
    return LibUtilities::eNodalTet;
}

//---------------------------------------
// Mapping functions
//---------------------------------------

int StdNodalTetExp::v_GetVertexMap(const int localVertexId,
                                   [[maybe_unused]] bool useCoeffPacking)
{
    ASSERTL0(localVertexId >= 0 && localVertexId <= 3,
             "Local Vertex ID must be between 0 and 3");
    return localVertexId;
}

void StdNodalTetExp::v_GetBoundaryMap(Array<OneD, unsigned int> &outarray)
{
    unsigned int i;
    const unsigned int nBndryCoeff = NumBndryCoeffs();

    if (outarray.size() != nBndryCoeff)
    {
        outarray = Array<OneD, unsigned int>(nBndryCoeff);
    }

    for (i = 0; i < nBndryCoeff; i++)
    {
        outarray[i] = i;
    }
}

void StdNodalTetExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    unsigned int i;
    const unsigned int nBndryCoeff = NumBndryCoeffs();

    if (outarray.size() != m_ncoeffs - nBndryCoeff)
    {
        outarray = Array<OneD, unsigned int>(m_ncoeffs - nBndryCoeff);
    }

    for (i = nBndryCoeff; i < m_ncoeffs; i++)
    {
        outarray[i - nBndryCoeff] = i;
    }
}

//---------------------------------------
// Wrapper functions
//---------------------------------------

DNekMatSharedPtr StdNodalTetExp::v_GenMatrix(const StdMatrixKey &mkey)
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

DNekMatSharedPtr StdNodalTetExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return StdNodalTetExp::v_GenMatrix(mkey);
}
} // namespace Nektar::StdRegions
