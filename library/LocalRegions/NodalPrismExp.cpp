///////////////////////////////////////////////////////////////////////////////
//
// File: NodalPrismExp.cpp
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
// Description: NodalPrismExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/NodalPrismExp.h>

namespace Nektar::LocalRegions
{
NodalPrismExp::NodalPrismExp(const LibUtilities::BasisKey &Ba,
                             const LibUtilities::BasisKey &Bb,
                             const LibUtilities::BasisKey &Bc,
                             const LibUtilities::PointsType Ntype,
                             SpatialDomains::Geometry3D *geom)
    : StdExpansion(LibUtilities::StdNodalPrismData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                   3, Ba, Bb, Bc),
      StdExpansion3D(LibUtilities::StdNodalPrismData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                     Ba, Bb, Bc),
      StdPrismExp(Ba, Bb, Bc), StdNodalPrismExp(Ba, Bb, Bc, Ntype),
      Expansion(geom), Expansion3D(geom), PrismExp(Ba, Bb, Bc, geom),
      m_matrixManager(
          std::bind(&Expansion3D::CreateMatrix, this, std::placeholders::_1)),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1))
{
}

NodalPrismExp::NodalPrismExp(const NodalPrismExp &T)
    : StdExpansion(T), StdExpansion3D(T), StdPrismExp(T), StdNodalPrismExp(T),
      Expansion(T), Expansion3D(T), PrismExp(T),
      m_matrixManager(T.m_matrixManager),
      m_staticCondMatrixManager(T.m_staticCondMatrixManager)
{
}

void NodalPrismExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp(m_ncoeffs);
    NodalToModal(inarray, tmp);
    StdPrismExp::v_BwdTrans(tmp, outarray);
}

void NodalPrismExp::v_IProductWRTBase(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    PrismExp::v_IProductWRTBase(inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

void NodalPrismExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    PrismExp::v_IProductWRTDerivBase(dir, inarray, outarray);
    NodalToModalTranspose(outarray, outarray);
}

StdRegions::StdExpansionSharedPtr NodalPrismExp::v_GetStdExp(void) const
{

    return MemoryManager<StdRegions::StdNodalPrismExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey(), m_base[1]->GetBasisKey(),
        m_base[2]->GetBasisKey(), m_nodalPointsKey.GetPointsType());
}

StdRegions::StdExpansionSharedPtr NodalPrismExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());
    LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(), 2,
                                 m_base[1]->GetPointsKey());
    LibUtilities::BasisKey bkey2(m_base[2]->GetBasisType(), 2,
                                 m_base[2]->GetPointsKey());

    return MemoryManager<StdRegions::StdNodalPrismExp>::AllocateSharedPtr(
        bkey0, bkey1, bkey2, m_nodalPointsKey.GetPointsType());
}

void NodalPrismExp::v_ExtractDataToCoeffs(
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

DNekMatSharedPtr NodalPrismExp::v_CreateStdMatrix(
    const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey0   = m_base[0]->GetBasisKey();
    LibUtilities::BasisKey bkey1   = m_base[1]->GetBasisKey();
    LibUtilities::BasisKey bkey2   = m_base[2]->GetBasisKey();
    LibUtilities::PointsType ntype = m_nodalPointsKey.GetPointsType();
    StdRegions::StdNodalPrismExpSharedPtr tmp =
        MemoryManager<StdNodalPrismExp>::AllocateSharedPtr(bkey0, bkey1, bkey2,
                                                           ntype);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr NodalPrismExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}

DNekScalBlkMatSharedPtr NodalPrismExp::v_GetLocStaticCondMatrix(
    const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}

void NodalPrismExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

void NodalPrismExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void NodalPrismExp::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

void NodalPrismExp::v_LaplacianMatrixOp(
    const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void NodalPrismExp::v_WeakDerivMatrixOp(
    const int i, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void NodalPrismExp::v_HelmholtzMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

} // namespace Nektar::LocalRegions
