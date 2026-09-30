///////////////////////////////////////////////////////////////////////////////
//
// File: StdExpansion1D.cpp
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
// which are common to 1d expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion1D.h>

#include <LibUtilities/BasicUtils/NekInline.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <StdRegions/Operators/SwitchLevel1.h>
#include <StdRegions/Operators/SwitchLevel2.h>

namespace Nektar::StdRegions
{
// Declaration of scalar routine
using vec_t = tinysimd::scalarT<double>;
#include <StdRegions/Operators/PhysDerivSumFacStdKernels.hpp>

StdExpansion1D::StdExpansion1D(
    [[maybe_unused]] int numcoeffs,
    [[maybe_unused]] const LibUtilities::BasisKey &Ba)
{
}

//----------------------------
// Differentiation Methods
//-----------------------------
void StdExpansion1D::PhysTensorDeriv(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad    = GetTotPoints();
    NekDouble *D = m_base[0]->GetD()->GetRawPtr();
    Array<OneD, const NekDouble> intmp;

    // copy inarray data if inarray and outarray are the same.
    if (inarray.data() == outarray.data())
    {
        Array<OneD, NekDouble> wsp(nquad);
        CopyArray(inarray, wsp);
        intmp = wsp;
    }
    else
    {
        intmp = inarray;
    }

    // Switch statment using boost_pp and macros. This unfolls into a
    // nested switch statement which runs from SMIN to SMAX for quadratrure
    // order. If you want to see it unwrapped compile in verbose mode and add
    // --preprocess to the c++ command. Default case
#undef PHYSDERIV_Q
#define PHYSDERIV_Q(r, i)                                                      \
    case NQ1(i):                                                               \
        PhysDerivTensor1DKernel(NQ1(i), (const vec_t *)intmp.data(),           \
                                (const vec_t *)D, (vec_t *)outarray.data());   \
        break;

    // templated cases on  standard quadrature
    // usage where quad order goes from SMIN to SMAX
    switch (nquad)
    {
        BOOST_PP_FOR((SMIN, SMAX), STDLEV1TEST, STDLEV1UPDATE, PHYSDERIV_Q);
        default:
            PhysDerivTensor1DKernel(nquad, (const vec_t *)intmp.data(),
                                    (const vec_t *)D, (vec_t *)outarray.data());
            break;
    }
}

void StdExpansion1D::v_PhysDeriv([[maybe_unused]] const int dir,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray)
{
    ASSERTL1(dir == 0, "input dir is out of range");
    v_PhysDeriv(inarray, outarray, NullNekDouble1DArray, NullNekDouble1DArray);
}

NekDouble StdExpansion1D::v_StdPhysEvaluate(
    const Array<OneD, const NekDouble> &Lcoord,
    const Array<OneD, const NekDouble> &physvals)
{
    ASSERTL2(Lcoord[0] >= -1 - NekConstants::kNekZeroTol, "Lcoord[0] < -1");
    ASSERTL2(Lcoord[0] <= 1 + NekConstants::kNekZeroTol, "Lcoord[0] >  1");

    return StdExpansion::BaryEvaluate<0>(Lcoord[0], &physvals[0]);
}

/** \brief Inner product of \a inarray over region with respect to the
 *  expansion basis (this)->m_base[0] and return in \a outarray
 *
 *  Wrapper call to \a IProductWRTBaseKernel()
 *
 *  @param inarray - Array of function values evaluated at the physical
 *  collocation points
 *  @param outarray - The values of the inner product with respect to
 *  each basis over region will be stored in the array \a outarray as
 *  output of the function
 */
void StdExpansion1D::v_IProductWRTBase(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    if (m_base[0]->Collocation())
    {
        v_MultiplyByStdQuadratureMetric(inarray, outarray);
    }
    else
    {
        const Array<OneD, const NekDouble> one(1, 1.0);
        v_IProductWRTBaseKernel(m_base[0]->GetBdata(), inarray, outarray, one,
                                false);
    }
}

void StdExpansion1D::IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, const NekDouble> &jac,
    const bool Deformed)
{
    v_IProductWRTBaseKernel(base0, inarray, outarray, jac, Deformed);
}

void StdExpansion1D::v_PhysInterp(std::shared_ptr<StdExpansion> fromExp,
                                  const Array<OneD, const NekDouble> &fromData,
                                  Array<OneD, NekDouble> &toData,
                                  [[maybe_unused]] bool Transpose)
{
    LibUtilities::Interp1D(fromExp->GetBasis(0)->GetPointsKey(), fromData,
                           m_base[0]->GetPointsKey(), toData);
}

void StdExpansion1D::v_MultiplyByStdQuadratureMetric(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();

    for (int j = 0; j < nquad0; ++j)
    {
        outarray[j] = inarray[j] * m_weights[0][j];
    }
}
// up to here
} // namespace Nektar::StdRegions
