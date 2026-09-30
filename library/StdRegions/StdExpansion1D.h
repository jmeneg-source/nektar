///////////////////////////////////////////////////////////////////////////////
//
// File: StdExpansion1D.h
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 1d expansion. Typically this inolves physical
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXP1D_H
#define STDEXP1D_H

#include <StdRegions/StdExpansion.h>

namespace Nektar::StdRegions
{
class StdExpansion1D : virtual public StdExpansion
{
public:
    STD_REGIONS_EXPORT StdExpansion1D(int numcoeffs,
                                      const LibUtilities::BasisKey &Ba);
    STD_REGIONS_EXPORT StdExpansion1D()                        = default;
    STD_REGIONS_EXPORT StdExpansion1D(const StdExpansion1D &T) = default;
    STD_REGIONS_EXPORT ~StdExpansion1D() override              = default;

    // find derivative of u (inarray) at all coords points
    STD_REGIONS_EXPORT inline NekDouble BaryTensorDeriv(
        const Array<OneD, NekDouble> &coord,
        const Array<OneD, const NekDouble> &inarray,
        std::array<NekDouble, 3> &firstOrderDerivs)
    {
        return StdExpansion::BaryEvaluate<0, true>(coord[0], &inarray[0],
                                                   firstOrderDerivs[0]);
    }

    // find derivative/2nd Derivative of u (inarray) at all coords points
    STD_REGIONS_EXPORT inline NekDouble BaryTensorDeriv(
        const Array<OneD, NekDouble> &coord,
        const Array<OneD, const NekDouble> &inarray,
        std::array<NekDouble, 3> &firstOrderDerivs,
        std::array<NekDouble, 6> &secondOrderDerivs)
    {
        return StdExpansion::BaryEvaluate<0, true, true>(
            coord[0], &inarray[0], firstOrderDerivs[0], secondOrderDerivs[0]);
    }

    STD_REGIONS_EXPORT void IProductWRTBaseKernel(
        const Array<OneD, const NekDouble> &base0,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const Array<OneD, const NekDouble> &jac, const bool Deformed);

protected:
    /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
     *  physical quadrature points given by \a inarray and return in
     *  \a outarray.
     *
     *  \param inarray array of a function evaluated at the quadrature
     *  points
     *  \param outarray the resulting array of the derivative \f$
     *  du/d_{\xi_1}|_{\xi_{1i}} \f$ will be stored in the array
     *  \a outarray as output of the function
     */
    STD_REGIONS_EXPORT void PhysTensorDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);
    STD_REGIONS_EXPORT void v_PhysDeriv(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    using StdExpansion::v_PhysDeriv;

    STD_REGIONS_EXPORT NekDouble
    v_StdPhysEvaluate(const Array<OneD, const NekDouble> &coords,
                      const Array<OneD, const NekDouble> &physvals) override;

    STD_REGIONS_EXPORT void v_IProductWRTBase(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    STD_REGIONS_EXPORT virtual void v_IProductWRTBaseKernel(
        const Array<OneD, const NekDouble> &base0,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const Array<OneD, const NekDouble> &jac, const bool Deformed) = 0;

    STD_REGIONS_EXPORT void v_PhysInterp(
        std::shared_ptr<StdExpansion> fromExp,
        const Array<OneD, const NekDouble> &fromData,
        Array<OneD, NekDouble> &toData, bool Transpose) override;

    STD_REGIONS_EXPORT void v_MultiplyByStdQuadratureMetric(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    int v_GetShapeDimension() const final
    {
        return 1;
    }
    bool v_IsCollocatedBasis() const final
    {
        return ((m_base[0]->Collocation()));
    }
};

typedef std::shared_ptr<StdExpansion1D> StdExpansion1DSharedPtr;

} // namespace Nektar::StdRegions

#endif // STDEXP1D_H
