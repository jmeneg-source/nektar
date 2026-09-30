///////////////////////////////////////////////////////////////////////////////
//
// File: NodalTetExp.cpp
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
// Description: NodalTetExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/NodalTetExp.h>

using namespace std;

namespace Nektar::LocalRegions
{
NodalTetExp::NodalTetExp(const LibUtilities::BasisKey &Ba,
                         const LibUtilities::BasisKey &Bb,
                         const LibUtilities::BasisKey &Bc,
                         const LibUtilities::PointsType Ntype,
                         SpatialDomains::Geometry3D *geom)
    : StdExpansion(LibUtilities::StdNodalTetData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                   3, Ba, Bb, Bc),
      StdExpansion3D(LibUtilities::StdNodalTetData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                     Ba, Bb, Bc),
      StdTetExp(Ba, Bb, Bc), StdNodalTetExp(Ba, Bb, Bc, Ntype), Expansion(geom),
      Expansion3D(geom), TetExp(Ba, Bb, Bc, geom),
      m_matrixManager(
          std::bind(&Expansion3D::CreateMatrix, this, std::placeholders::_1)),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1))
{
}

NodalTetExp::NodalTetExp(const NodalTetExp &T)
    : StdExpansion(T), StdExpansion3D(T), StdTetExp(T), StdNodalTetExp(T),
      Expansion(T), Expansion3D(T), TetExp(T),
      m_matrixManager(T.m_matrixManager),
      m_staticCondMatrixManager(T.m_staticCondMatrixManager)
{
}

void NodalTetExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp(m_ncoeffs);
    NodalToModal(inarray, tmp);
    StdTetExp::v_BwdTrans(tmp, outarray);
}

void NodalTetExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    TetExp::v_IProductWRTBase(inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

void NodalTetExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    TetExp::v_IProductWRTDerivBase(dir, inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

StdRegions::StdExpansionSharedPtr NodalTetExp::v_GetStdExp(void) const
{

    return MemoryManager<StdRegions::StdNodalTetExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey(), m_base[1]->GetBasisKey(),
        m_base[2]->GetBasisKey(), m_nodalPointsKey.GetPointsType());
}

StdRegions::StdExpansionSharedPtr NodalTetExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());
    LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(), 2,
                                 m_base[1]->GetPointsKey());
    LibUtilities::BasisKey bkey2(m_base[2]->GetBasisType(), 2,
                                 m_base[2]->GetPointsKey());

    return MemoryManager<StdRegions::StdNodalTetExp>::AllocateSharedPtr(
        bkey0, bkey1, bkey2, m_nodalPointsKey.GetPointsType());
}

void NodalTetExp::v_ExtractDataToCoeffs(
    const NekDouble *data, const std::vector<unsigned int> &nummodes,
    const int mode_offset, NekDouble *coeffs,
    [[maybe_unused]] std::vector<LibUtilities::BasisType> &fromType)
{
    Array<OneD, NekDouble> modes(m_ncoeffs);
    Expansion::ExtractDataToCoeffs(data, nummodes, mode_offset, &modes[0],
                                   fromType);

    Array<OneD, NekDouble> nodes(m_ncoeffs, coeffs, eArrayWrapper);
    ModalToNodal(modes, nodes);
}

DNekMatSharedPtr NodalTetExp::v_CreateStdMatrix(
    const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey0   = m_base[0]->GetBasisKey();
    LibUtilities::BasisKey bkey1   = m_base[1]->GetBasisKey();
    LibUtilities::BasisKey bkey2   = m_base[2]->GetBasisKey();
    LibUtilities::PointsType ntype = m_nodalPointsKey.GetPointsType();
    StdRegions::StdNodalTetExpSharedPtr tmp =
        MemoryManager<StdNodalTetExp>::AllocateSharedPtr(bkey0, bkey1, bkey2,
                                                         ntype);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr NodalTetExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}

DNekScalBlkMatSharedPtr NodalTetExp::v_GetLocStaticCondMatrix(
    const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}

void NodalTetExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

void NodalTetExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void NodalTetExp::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

void NodalTetExp::v_LaplacianMatrixOp(
    const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void NodalTetExp::v_WeakDerivMatrixOp(
    const int i, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void NodalTetExp::v_HelmholtzMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

} // namespace Nektar::LocalRegions
