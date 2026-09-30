///////////////////////////////////////////////////////////////////////////////
//
// File: NodalPrismExp.h
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
// Description: Header for NodalPrismExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALPRISMEXP_H
#define NODALPRISMEXP_H

#include <SpatialDomains/PrismGeom.h>
#include <StdRegions/StdNodalPrismExp.h>

#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/PrismExp.h>

namespace Nektar::LocalRegions
{

class NodalPrismExp final : virtual public StdRegions::StdNodalPrismExp,
                            virtual public PrismExp
{
public:
    /** \brief Constructor using BasisKey class for quadrature
    points and order definition */
    LOCAL_REGIONS_EXPORT NodalPrismExp(const LibUtilities::BasisKey &Ba,
                                       const LibUtilities::BasisKey &Bb,
                                       const LibUtilities::BasisKey &Bc,
                                       const LibUtilities::PointsType Ntype,
                                       SpatialDomains::Geometry3D *geom);

    /// Copy Constructor
    LOCAL_REGIONS_EXPORT NodalPrismExp(const NodalPrismExp &T);

    /// Destructor
    LOCAL_REGIONS_EXPORT ~NodalPrismExp() override = default;

protected:
    //---------------------------------------
    // Transforms
    //---------------------------------------
    LOCAL_REGIONS_EXPORT void v_BwdTrans(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------------------
    // Inner product functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT void v_IProductWRTBase(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_IProductWRTDerivBase(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------------------
    // Evaluation functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT void v_GetCoord(
        const Array<OneD, const NekDouble> &Lcoords,
        Array<OneD, NekDouble> &coords) override
    {
        PrismExp::v_GetCoord(Lcoords, coords);
    }

    LOCAL_REGIONS_EXPORT void v_GetCoords(
        Array<OneD, NekDouble> &coords_1, Array<OneD, NekDouble> &coords_2,
        Array<OneD, NekDouble> &coords_3) override
    {
        PrismExp::v_GetCoords(coords_1, coords_2, coords_3);
    }

    LOCAL_REGIONS_EXPORT StdRegions::StdExpansionSharedPtr v_GetStdExp(
        void) const override;
    LOCAL_REGIONS_EXPORT StdRegions::StdExpansionSharedPtr v_GetLinStdExp(
        void) const override;
    LOCAL_REGIONS_EXPORT void v_ExtractDataToCoeffs(
        const NekDouble *data, const std::vector<unsigned int> &nummodes,
        const int mode_offset, NekDouble *coeffs,
        [[maybe_unused]] std::vector<LibUtilities::BasisType> &fromType)
        override;

    //---------------------------------------
    // Matrix creation functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT DNekMatSharedPtr
    v_GenMatrix(const StdRegions::StdMatrixKey &mkey) override
    {
        return PrismExp::v_GenMatrix(mkey);
    }
    LOCAL_REGIONS_EXPORT DNekMatSharedPtr
    v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr
    v_GetLocMatrix(const MatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT DNekScalBlkMatSharedPtr
    v_GetLocStaticCondMatrix(const MatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_DropLocMatrix(const MatrixKey &mkey) override;

    //----------------------------------------
    // Matrix Operators
    //----------------------------------------
    LOCAL_REGIONS_EXPORT void v_MassMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_LaplacianMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_LaplacianMatrixOp(
        const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_WeakDerivMatrixOp(
        const int i, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_HelmholtzMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;

private:
    LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess>
        m_matrixManager;
    LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess>
        m_staticCondMatrixManager;
};

typedef std::shared_ptr<NodalPrismExp> NodalPrismExpSharedPtr;
typedef std::vector<NodalPrismExpSharedPtr> NodalPrismExpVector;

} // namespace Nektar::LocalRegions

#endif // NODAL_PRISMEXP_H
