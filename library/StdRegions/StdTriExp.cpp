///////////////////////////////////////////////////////////////////////////////
//
// File: StdTriExp.cpp
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
// Description: Triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdSegExp.h> // for StdSegExp, etc
#include <StdRegions/StdTriExp.h>

using namespace std;
#include <LibUtilities/BasicUtils/NekInline.hpp>
#include <StdRegions/Operators/SwitchLevel1.h>
#include <StdRegions/Operators/SwitchLevel2.h>

namespace Nektar::StdRegions
{
// Declaration of scalar routine
using vec_t = tinysimd::scalarT<double>;
#include <StdRegions/Operators/BwdTransSumFacStdKernels.hpp>
#include <StdRegions/Operators/IProductWRTBaseSumFacStdKernels.hpp>

StdTriExp::StdTriExp(const LibUtilities::BasisKey &Ba,
                     const LibUtilities::BasisKey &Bb)
    : StdExpansion(LibUtilities::StdTriData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes()),
                   2, Ba, Bb),
      StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes()),
                     Ba, Bb)
{
    ASSERTL0(Ba.GetNumModes() <= Bb.GetNumModes(),
             "order in 'a' direction is higher than order "
             "in 'b' direction");

    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());

    StdFacKey w1key(eWeights1, Bb);
    // get weights[1] from manager where points are rescaled
    m_weights.push_back(GetStdFac(w1key));
}

//-----------------------------
// Differentiation Methods
//-----------------------------
/**
 * \brief Calculate the derivative of the physical points.
 *
 * \f$ \frac{\partial u}{\partial  x_1} =  \left .
 * \frac{2.0}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
 * \right |_{\eta_2}\f$
 *
 * \f$ \frac{\partial u}{\partial  x_2} =  \left .
 * \frac{1+\eta_1}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
 * \right |_{\eta_2}  + \left . \frac{\partial u}{\partial d\eta_2}
 * \right |_{\eta_1}  \f$
 */
void StdTriExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &out_d0,
                               Array<OneD, NekDouble> &out_d1,
                               [[maybe_unused]] Array<OneD, NekDouble> &out_d2)
{
    int i;
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    Array<OneD, NekDouble> wsp(std::max(nquad0, nquad1));

    const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
    const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();

    // set up geometric factor: 2/(1-z1)
    Vmath::Sadd(nquad1, -1.0, z1, 1, wsp, 1);
    Vmath::Sdiv(nquad1, -2.0, wsp, 1, wsp, 1);

    if (out_d0.size() > 0)
    {
        PhysTensorDeriv(inarray, out_d0, out_d1);

        for (i = 0; i < nquad1; ++i)
        {
            Blas::Dscal(nquad0, wsp[i], &out_d0[0] + i * nquad0, 1);
        }

        // if no d1 required do not need to calculate both deriv
        if (out_d1.size() > 0)
        {
            // set up geometric factor: (1+z0)/2
            Vmath::Sadd(nquad0, 1.0, z0, 1, wsp, 1);
            Vmath::Smul(nquad0, 0.5, wsp, 1, wsp, 1);

            for (i = 0; i < nquad1; ++i)
            {
                Vmath::Vvtvp(nquad0, &wsp[0], 1, &out_d0[0] + i * nquad0, 1,
                             &out_d1[0] + i * nquad0, 1,
                             &out_d1[0] + i * nquad0, 1);
            }
        }
    }
    else if (out_d1.size() > 0)
    {
        Array<OneD, NekDouble> diff0(nquad0 * nquad1);
        PhysTensorDeriv(inarray, diff0, out_d1);

        for (i = 0; i < nquad1; ++i)
        {
            Blas::Dscal(nquad0, wsp[i], &diff0[0] + i * nquad0, 1);
        }

        Vmath::Sadd(nquad0, 1.0, z0, 1, wsp, 1);
        Vmath::Smul(nquad0, 0.5, wsp, 1, wsp, 1);

        for (i = 0; i < nquad1; ++i)
        {
            Vmath::Vvtvp(nquad0, &wsp[0], 1, &diff0[0] + i * nquad0, 1,
                         &out_d1[0] + i * nquad0, 1, &out_d1[0] + i * nquad0,
                         1);
        }
    }
}

//---------------------------------------
// Transforms
//---------------------------------------

/**
 * \brief Backward tranform for triangular elements
 *
 * @note 'q' (base[1]) runs fastest in this element.
 */
void StdTriExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &outarray)
{
    ASSERTL2((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B) ||
                 (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
             "Basis[1] is not of general tensor type");

    const Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
    const Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();

    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();

    std::vector<vec_t, tinysimd::allocator<vec_t>> wsp0(nmodes0);
    bool isModified = (m_base[0]->GetBasisType() == LibUtilities::eModified_A);

// Swith statment using boost_pp and macros. This unfolls intwo a
// nested swtich statement where the outer swtich statement runs
// from SMIN to SMAX for modal order and the inner switch
// statemets run from the outer value of the case to 2*SMAX for
// the quadrature order. If you want to see it unwrapped compile
// in verbose mode and add --preprocess to the c++ command.
// Default case
#undef BWDTRANS_DEF
#define BWDTRANS_DEF                                                           \
    BwdTransTriKernel(nmodes0, nmodes1, nquad0, nquad1, isModified,            \
                      (const vec_t *)base0.data(),                             \
                      (const vec_t *)base1.data(), wsp0.data(),                \
                      (const vec_t *)inarray.data(), (vec_t *)outarray.data())

// Inner loop case over quarature points
#undef BWDTRANS_Q
#define BWDTRANS_Q(r, i)                                                       \
    case NQ(i):                                                                \
        BwdTransTriKernel(NM(i), NM(i), NQ(i), NQ_M1(i), isModified,           \
                          (const vec_t *)base0.data(),                         \
                          (const vec_t *)base1.data(), wsp0.data(),            \
                          (const vec_t *)inarray.data(),                       \
                          (vec_t *)outarray.data());                           \
        break;

// outer loop case over modes
#undef BWDTRANS_M
#define BWDTRANS_M(r, i)                                                       \
    case NM(i):                                                                \
    {                                                                          \
        switch (nquad0)                                                        \
        {                                                                      \
            BOOST_PP_FOR_##r((NM(i), NM_P1(i), BOOST_PP_MUL(2, NM(i))),        \
                             STDLEV2TEST1, STDLEV2UPDATE1, BWDTRANS_Q) default \
                : BWDTRANS_DEF;                                                \
            break;                                                             \
        }                                                                      \
    }                                                                          \
    break;

    // templated cases on equi-ordered modes and standard quad
    // usage where quad order goes from mode order to 2(*mode
    // order)
    if ((nmodes0 == nmodes1) && (nquad0 == nquad1 + 1))
    {
        switch (nmodes0)
        {
            BOOST_PP_FOR((SMIN, 0, SMAX), STDLEV2TEST, STDLEV2UPDATE,
                         BWDTRANS_M)
            default:
                BWDTRANS_DEF;
                break;
        }
    }
    else
    {
        BWDTRANS_DEF;
    }
}

void StdTriExp::v_FwdTransBndConstrained(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int i, j;
    int npoints[2] = {m_base[0]->GetNumPoints(), m_base[1]->GetNumPoints()};
    int nmodes[2]  = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes()};

    fill(outarray.data(), outarray.data() + m_ncoeffs, 0.0);

    Array<OneD, NekDouble> physEdge[3];
    Array<OneD, NekDouble> coeffEdge[3];
    for (i = 0; i < 3; i++)
    {
        physEdge[i]  = Array<OneD, NekDouble>(npoints[i != 0]);
        coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i != 0]);
    }

    for (i = 0; i < npoints[0]; i++)
    {
        physEdge[0][i] = inarray[i];
    }

    for (i = 0; i < npoints[1]; i++)
    {
        physEdge[1][i] = inarray[npoints[0] - 1 + i * npoints[0]];
        physEdge[2][i] =
            inarray[(npoints[1] - 1) * npoints[0] - i * npoints[0]];
    }

    StdSegExpSharedPtr segexp[2] = {
        MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
            m_base[0]->GetBasisKey()),
        MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
            m_base[1]->GetBasisKey())};

    Array<OneD, unsigned int> mapArray;
    Array<OneD, int> signArray;
    NekDouble sign;

    for (i = 0; i < 3; i++)
    {
        segexp[i != 0]->FwdTransBndConstrained(physEdge[i], coeffEdge[i]);

        GetTraceToElementMap(i, mapArray, signArray);
        for (j = 0; j < nmodes[i != 0]; j++)
        {
            sign                  = (NekDouble)signArray[j];
            outarray[mapArray[j]] = sign * coeffEdge[i][j];
        }
    }

    Array<OneD, NekDouble> tmp0(m_ncoeffs);
    Array<OneD, NekDouble> tmp1(m_ncoeffs);

    StdMatrixKey masskey(eMass, DetShapeType(), *this);
    MassMatrixOp(outarray, tmp0, masskey);
    v_IProductWRTBase(inarray, tmp1);

    Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

    // get Mass matrix inverse (only of interior DOF)
    // use block (1,1) of the static condensed system
    // note: this block alreay contains the inverse matrix
    DNekMatSharedPtr matsys =
        (m_stdStaticCondMatrixManager[masskey])->GetBlock(1, 1);

    int nBoundaryDofs = v_NumBndryCoeffs();
    int nInteriorDofs = m_ncoeffs - nBoundaryDofs;

    Array<OneD, NekDouble> rhs(nInteriorDofs);
    Array<OneD, NekDouble> result(nInteriorDofs);

    v_GetInteriorMap(mapArray);

    for (i = 0; i < nInteriorDofs; i++)
    {
        rhs[i] = tmp1[mapArray[i]];
    }

    Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, 1.0, &(matsys->GetPtr())[0],
                nInteriorDofs, rhs.data(), 1, 0.0, result.data(), 1);

    for (i = 0; i < nInteriorDofs; i++)
    {
        outarray[mapArray[i]] = result[i];
    }
}

//---------------------------------------
// Inner product functions
//---------------------------------------
/** \brief Inner product of \a inarray over region with respect to the
 *  expansion basis (this)->m_base[0] and return in \a outarray
 *
 *  @param base0 - An array containing the values of the basis in the
 *  0-direction at the quarature poitns
 *  @param base1 - An array containing the values of the basis in the
 *  1-direction at the quarature poitns
 *  @param inarray - Array of values evaluated at the physical
 *  quadrature points
 *  @param outarray the values of the inner product with respect to
 *  each basis over region will be stored in the array \a outarray as
 *  output of the function
 *  @param jac - An array of size 1 if not deformed or the number of
 *  quadrature points if deformed holding the values of the jacobian
 *  @param Deformed - a bool identifying if the inner product is to be
 *  treated as a deformed or regular integration which just relates to
 *  how the \param jac array is treated
 */
void StdTriExp::v_IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &base1,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, NekDouble> &jac,
    const bool Deformed, [[maybe_unused]] const bool CollDir0,
    [[maybe_unused]] const bool CollDir1)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int order0 = m_base[0]->GetNumModes();
    int order1 = m_base[1]->GetNumModes();

    const bool isModified =
        (m_base[0]->GetBasisType() == LibUtilities::eModified_A);

    std::vector<vec_t, tinysimd::allocator<vec_t>> wsp0(nquad1);

    // Swith statment using boost_pp and macros. This unfolls intwo a
    // nested swtich statement where the outer swtich statement runs
    // from SMIN to SMAX for modal order and the inner switch
    // statemets run from the outer value of the case to 2*SMAX for
    // the quadrature order. If you want to see it unwrapped compile
    // in verbose mode and add --preprocess to the c++ command.
    if (Deformed)
    {
        // Default case
#undef IPRODUCTWRTBASE_DEF
#define IPRODUCTWRTBASE_DEF                                                    \
    IProductTriKernel<false, false, true>(                                     \
        order0, order1, nquad0, nquad1, isModified,                            \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)m_weights[0].data(),       \
        (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductTriKernel<false, false, true>(                                 \
            NM(i), NM(i), NQ(i), NQ_M1(i), isModified,                         \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)m_weights[0].data(),   \
            (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)outarray.data());                   \
        break;

        // outer loop case over modes
#undef IPRODUCTWRTBASE_M
#define IPRODUCTWRTBASE_M(r, i)                                                \
    case NM(i):                                                                \
    {                                                                          \
        switch (nquad0)                                                        \
        {                                                                      \
            BOOST_PP_FOR_##r((NM(i), NM_P1(i), BOOST_PP_MUL(2, NM(i))),        \
                             STDLEV2TEST1, STDLEV2UPDATE1,                     \
                             IPRODUCTWRTBASE_Q) default : IPRODUCTWRTBASE_DEF; \
            break;                                                             \
        }                                                                      \
    }                                                                          \
    break;

        // templated cases on equi-ordered modes and standard quad usage
        // where quad order goes from mode order to 2(*mode order)
        if ((order0 == order1) && (nquad0 == nquad1 + 1))
        {
            switch (order0)
            {
                BOOST_PP_FOR((SMIN, 0, SMAX), STDLEV2TEST, STDLEV2UPDATE,
                             IPRODUCTWRTBASE_M)
                default:
                    IPRODUCTWRTBASE_DEF;
                    break;
            }
        }
        else
        {
            IPRODUCTWRTBASE_DEF;
        }
    }
    else // non-deformed case
    {
        // Default case
#undef IPRODUCTWRTBASE_DEF
#define IPRODUCTWRTBASE_DEF                                                    \
    IProductTriKernel<false, false, false>(                                    \
        order0, order1, nquad0, nquad1, isModified,                            \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)m_weights[0].data(),       \
        (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductTriKernel<false, false, false>(                                \
            NM(i), NM(i), NQ(i), NQ_M1(i), isModified,                         \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)m_weights[0].data(),   \
            (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)outarray.data());                   \
        break;

        // outer loop case over modes
#undef IPRODUCTWRTBASE_M
#define IPRODUCTWRTBASE_M(r, i)                                                \
    case NM(i):                                                                \
    {                                                                          \
        switch (nquad0)                                                        \
        {                                                                      \
            BOOST_PP_FOR_##r((NM(i), NM_P1(i), BOOST_PP_MUL(2, NM(i))),        \
                             STDLEV2TEST1, STDLEV2UPDATE1,                     \
                             IPRODUCTWRTBASE_Q) default : IPRODUCTWRTBASE_DEF; \
            break;                                                             \
        }                                                                      \
    }                                                                          \
    break;

        // templated cases on equi-ordered modes and standard quad usage
        // where quad order goes from mode order to 2(*mode order)
        if ((order0 == order1) && (nquad0 == nquad1 + 1))
        {
            switch (order0)
            {
                BOOST_PP_FOR((SMIN, 0, SMAX), STDLEV2TEST, STDLEV2UPDATE,
                             IPRODUCTWRTBASE_M)
                default:
                    IPRODUCTWRTBASE_DEF;
                    break;
            }
        }
        else
        {
            IPRODUCTWRTBASE_DEF;
        }
    }
}

void StdTriExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nqtot  = nquad0 * nquad1;
    Array<OneD, NekDouble> tmpQuad(nqtot);

    // multiply by 2/(1-z1)
    StdFacKey fackey(eTwoOverOneMinusZ1, m_base[1]->GetBasisKey());
    Array<OneD, const NekDouble> gfac = GetStdFac(fackey);
    for (int i = 0; i < nquad1; ++i)
    {
        Vmath::Smul(nquad0, gfac[i], &inarray[0] + i * nquad0, 1,
                    &tmpQuad[0] + i * nquad0, 1);
    }

    const Array<OneD, const NekDouble> one(1, 1.0);
    switch (dir)
    {
        case 0:
            v_IProductWRTBaseKernel(m_base[0]->GetDbdata(),
                                    m_base[1]->GetBdata(), tmpQuad, outarray,
                                    one, false);
            break;
        case 1:
        {
            Array<OneD, NekDouble> tmpCoeff(m_ncoeffs);

            // multiply by 0.5*(1-z0)
            StdFacKey fackey1(eHalfMultOnePlusZ0, m_base[0]->GetBasisKey());
            gfac = GetStdFac(fackey1);
            for (int i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0, &gfac[0], 1, &tmpQuad[0] + i * nquad0, 1,
                            &tmpQuad[0] + i * nquad0, 1);
            }

            v_IProductWRTBaseKernel(m_base[0]->GetDbdata(),
                                    m_base[1]->GetBdata(), tmpQuad, tmpCoeff,
                                    one, false);

            v_IProductWRTBaseKernel(m_base[0]->GetBdata(),
                                    m_base[1]->GetDbdata(), inarray, outarray,
                                    one, false);

            Vmath::Vadd(m_ncoeffs, &tmpCoeff[0], 1, &outarray[0], 1,
                        &outarray[0], 1);
            break;
        }
        default:
        {
            ASSERTL1(false, "input dir is out of range");
            break;
        }
    }
}

//---------------------------------------
// Evaluation functions
//---------------------------------------

void StdTriExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble> &xi,
                                         Array<OneD, NekDouble> &eta)
{
    NekDouble d1 = 1. - xi[1];
    if (fabs(d1) < NekConstants::kNekZeroTol)
    {
        if (d1 >= 0.)
        {
            d1 = NekConstants::kNekZeroTol;
        }
        else
        {
            d1 = -NekConstants::kNekZeroTol;
        }
    }
    eta[0] = 2. * (1. + xi[0]) / d1 - 1.0;
    eta[1] = xi[1];
}

void StdTriExp::v_LocCollapsedToLocCoord(
    const Array<OneD, const NekDouble> &eta, Array<OneD, NekDouble> &xi)
{
    xi[0] = (1.0 + eta[0]) * (1.0 - eta[1]) * 0.5 - 1.0;
    xi[1] = eta[1];
}

void StdTriExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
{
    int i, m;
    int nquad0                         = m_base[0]->GetNumPoints();
    int nquad1                         = m_base[1]->GetNumPoints();
    int order0                         = m_base[0]->GetNumModes();
    int order1                         = m_base[1]->GetNumModes();
    int mode0                          = 0;
    Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
    Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();

    ASSERTL2(mode <= m_ncoeffs, "calling argument mode is larger than "
                                "total expansion order");

    m = order1;
    for (i = 0; i < order0; ++i, m += order1 - i)
    {
        if (m > mode)
        {
            mode0 = i;
            break;
        }
    }

    // deal with top vertex mode in modified basis
    if (mode == 1 && m_base[0]->GetBasisType() == LibUtilities::eModified_A)
    {
        Vmath::Fill(nquad0 * nquad1, 1.0, outarray, 1);
    }
    else
    {
        for (i = 0; i < nquad1; ++i)
        {
            Vmath::Vcopy(nquad0, (NekDouble *)(base0.data() + mode0 * nquad0),
                         1, &outarray[0] + i * nquad0, 1);
        }
    }

    for (i = 0; i < nquad0; ++i)
    {
        Vmath::Vmul(nquad1, (NekDouble *)(base1.data() + mode * nquad1), 1,
                    &outarray[0] + i, nquad0, &outarray[0] + i, nquad0);
    }
}

NekDouble StdTriExp::v_PhysEvaluateBasis(
    const Array<OneD, const NekDouble> &coords, int mode)
{
    Array<OneD, NekDouble> coll(2);
    LocCoordToLocCollapsed(coords, coll);

    // From mode we need to determine mode0 and mode1 in the (p,q)
    // direction. mode1 can be directly inferred from mode.
    const int nm1   = m_base[1]->GetNumModes();
    const double c  = 1 + 2 * nm1;
    const int mode0 = floor(0.5 * (c - sqrt(c * c - 8 * mode)));

    if (mode == 1 && m_base[0]->GetBasisType() == LibUtilities::eModified_A)
    {
        // Account for collapsed vertex.
        return StdExpansion::BaryEvaluateBasis<1>(coll[1], 1);
    }
    else
    {
        return StdExpansion::BaryEvaluateBasis<0>(coll[0], mode0) *
               StdExpansion::BaryEvaluateBasis<1>(coll[1], mode);
    }
}

NekDouble StdTriExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    // Collapse coordinates
    Array<OneD, NekDouble> coll(2, 0.0);
    LocCoordToLocCollapsed(coord, coll);

    // If near singularity do the old interpolation matrix method
    if ((1 - coll[1]) < 1e-5)
    {
        int totPoints = GetTotPoints();
        Array<OneD, NekDouble> EphysDeriv0(totPoints), EphysDeriv1(totPoints);
        v_PhysDeriv(inarray, EphysDeriv0, EphysDeriv1, NullNekDouble1DArray);

        Array<OneD, DNekMatSharedPtr> I(2);
        I[0] = GetBase()[0]->GetI(coll);
        I[1] = GetBase()[1]->GetI(coll + 1);

        firstOrderDerivs[0] = PhysEvaluate(I, EphysDeriv0);
        firstOrderDerivs[1] = PhysEvaluate(I, EphysDeriv1);
        return PhysEvaluate(I, inarray);
    }

    // set up geometric factor: 2.0/(1.0-z1)
    NekDouble fac0 = 2 / (1 - coll[1]);

    NekDouble val = BaryTensorDeriv(coll, inarray, firstOrderDerivs);

    // Copy d0 into temp for d1
    NekDouble temp;
    temp = firstOrderDerivs[0];

    // Multiply by geometric factor
    firstOrderDerivs[0] = firstOrderDerivs[0] * fac0;

    // set up geometric factor: (1+z0)/(1-z1)
    NekDouble fac1 = fac0 * (coll[0] + 1) / 2;

    // Multiply out_d0 by geometric factor and add to out_d1
    firstOrderDerivs[1] += fac1 * temp;

    return val;
}

int StdTriExp::v_GetNverts() const
{
    return 3;
}

int StdTriExp::v_GetNtraces() const
{
    return 3;
}

LibUtilities::ShapeType StdTriExp::v_DetShapeType() const
{
    return LibUtilities::eTriangle;
}

int StdTriExp::v_NumBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
             "BasisType is not a boundary interior form");

    return 3 + (GetBasisNumModes(0) - 2) + 2 * (GetBasisNumModes(1) - 2);
}

int StdTriExp::v_NumDGBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
             "BasisType is not a boundary interior form");

    return GetBasisNumModes(0) + 2 * GetBasisNumModes(1);
}

int StdTriExp::v_GetTraceNcoeffs(const int i) const
{
    ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

    if (i == 0)
    {
        return GetBasisNumModes(0);
    }
    else
    {
        return GetBasisNumModes(1);
    }
}

int StdTriExp::v_GetTraceIntNcoeffs(const int i) const
{
    ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

    if (i == 0)
    {
        return GetBasisNumModes(0) - 2;
    }
    else
    {
        return GetBasisNumModes(1) - 2;
    }
}

int StdTriExp::v_GetTraceNumPoints(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 2), "edge id is out of range");

    if (i == 0)
    {
        return GetNumPoints(0);
    }
    else
    {
        return GetNumPoints(1);
    }
}

int StdTriExp::v_CalcNumberOfCoefficients(
    const std::vector<unsigned int> &nummodes, int &modes_offset)
{
    int nmodes = LibUtilities::StdTriData::getNumberOfCoefficients(
        nummodes[modes_offset], nummodes[modes_offset + 1]);
    modes_offset += 2;

    return nmodes;
}

void StdTriExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                            Array<OneD, NekDouble> &coords_1,
                            [[maybe_unused]] Array<OneD, NekDouble> &coords_2)
{
    Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
    Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
    int nq0                         = GetNumPoints(0);
    int nq1                         = GetNumPoints(1);
    int i, j;

    for (i = 0; i < nq1; ++i)
    {
        for (j = 0; j < nq0; ++j)
        {
            coords_0[i * nq0 + j] = (1 + z0[j]) * (1 - z1[i]) / 2.0 - 1.0;
        }
        Vmath::Fill(nq0, z1[i], &coords_1[0] + i * nq0, 1);
    }
}

bool StdTriExp::v_IsBoundaryInteriorExpansion() const
{
    return m_base[0]->GetBasisType() == LibUtilities::eModified_A &&
           m_base[1]->GetBasisType() == LibUtilities::eModified_B;
}

const LibUtilities::BasisKey StdTriExp::v_GetTraceBasisKey(
    const int i, [[maybe_unused]] const int j,
    [[maybe_unused]] bool UseGLL) const
{
    ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

    // Get basiskey (0 or 1) according to edge id i
    int dir = (i != 0);

    switch (m_base[dir]->GetBasisType())
    {
        case LibUtilities::eGLL_Lagrange:
        case LibUtilities::eModified_A:
        {
            switch (m_base[dir]->GetPointsType())
            {
                case LibUtilities::eGaussLegendreWithMP:
                case LibUtilities::eGaussLobattoLegendre:
                {
                    return m_base[dir]->GetBasisKey();
                }
                break;
                default:
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Unexpected points distribution " +
                                 LibUtilities::kPointsTypeStr
                                     [m_base[dir]->GetPointsType()] +
                                 " in StdTriExp::v_GetTraceBasisKey");
                }
            }
        }
        break;
        case LibUtilities::eModified_B:
        {
            switch (m_base[dir]->GetPointsType())
            {
                case LibUtilities::eGaussLegendreWithMP:
                case LibUtilities::eGaussLobattoLegendre:
                {
                    return LibUtilities::BasisKey(LibUtilities::eModified_A,
                                                  m_base[dir]->GetNumModes(),
                                                  m_base[dir]->GetPointsKey());
                }
                break;
                case LibUtilities::eGaussLegendreWithM:
                {
                    LibUtilities::PointsKey pkey(
                        m_base[dir]
                                ->GetBasisKey()
                                .GetPointsKey()
                                .GetNumPoints() +
                            1,
                        LibUtilities::eGaussLegendreWithMP);
                    return LibUtilities::BasisKey(LibUtilities::eModified_A,
                                                  m_base[dir]->GetNumModes(),
                                                  pkey);
                }
                break;
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                {
                    LibUtilities::PointsKey pkey(
                        m_base[dir]
                                ->GetBasisKey()
                                .GetPointsKey()
                                .GetNumPoints() +
                            1,
                        LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eModified_A,
                                                  m_base[dir]->GetNumModes(),
                                                  pkey);
                }
                break;
                // Currently this does not increase the points by
                // 1 since when using this quadrature we are
                // presuming it is already been increased by one
                // when comopared to
                // GaussRadauMAlpha1Beta0. Currently used in the
                // GJP option
                //
                // Note have put down it back to numpoints +1 to
                // test for use on tri faces and GJP.
                case LibUtilities::eGaussRadauMLegendre:
                {
                    LibUtilities::PointsKey pkey(
                        m_base[dir]
                                ->GetBasisKey()
                                .GetPointsKey()
                                .GetNumPoints() +
                            1,
                        LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eModified_A,
                                                  m_base[dir]->GetNumModes(),
                                                  pkey);
                }
                break;
                case LibUtilities::ePolyEvenlySpaced:
                {
                    LibUtilities::PointsKey pkey(
                        m_base[dir]
                                ->GetBasisKey()
                                .GetPointsKey()
                                .GetNumPoints() +
                            1,
                        LibUtilities::ePolyEvenlySpaced);
                    return LibUtilities::BasisKey(LibUtilities::eModified_A,
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
                                 " in StdTriExp::v_GetTraceBasisKey");
                }
            }
        }
        break;
        case LibUtilities::eOrtho_A:
        {
            switch (m_base[dir]->GetPointsType())
            {
                case LibUtilities::eGaussLegendreWithMP:
                case LibUtilities::eGaussLobattoLegendre:
                {
                    return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,
                                                  m_base[dir]->GetNumModes(),
                                                  m_base[dir]->GetPointsKey());
                }
                break;
                default:
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Unexpected points distribution " +
                                 LibUtilities::kPointsTypeStr
                                     [m_base[dir]->GetPointsType()] +
                                 " in StdTriExp::v_GetTraceBasisKey");
                }
            }
        }
        break;
        case LibUtilities::eOrtho_B:
        {
            switch (m_base[dir]->GetPointsType())
            {
                case LibUtilities::eGaussLegendreWithMP:
                case LibUtilities::eGaussLobattoLegendre:
                {
                    return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,
                                                  m_base[dir]->GetNumModes(),
                                                  m_base[dir]->GetPointsKey());
                }
                break;
                case LibUtilities::eGaussLegendreWithM:
                {
                    LibUtilities::PointsKey pkey(
                        m_base[dir]
                                ->GetBasisKey()
                                .GetPointsKey()
                                .GetNumPoints() +
                            1,
                        LibUtilities::eGaussLegendreWithMP);
                    return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,
                                                  m_base[dir]->GetNumModes(),
                                                  pkey);
                }
                break;
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
                                 " in StdTriExp::v_GetTraceBasisKey");
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

int StdTriExp::v_GetVertexMap(const int localVertexId, bool useCoeffPacking)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eModified_B,
             "Mapping not defined for this type of basis");

    int localDOF = 0;
    if (useCoeffPacking == true)
    {
        switch (localVertexId)
        {
            case 0:
            {
                localDOF = 0;
                break;
            }
            case 1:
            {
                localDOF = 1;
                break;
            }
            case 2:
            {
                localDOF = m_base[1]->GetNumModes();
                break;
            }
            default:
            {
                ASSERTL0(false, "eid must be between 0 and 2");
                break;
            }
        }
    }
    else // follow book format for vertex indexing.
    {
        switch (localVertexId)
        {
            case 0:
            {
                localDOF = 0;
                break;
            }
            case 1:
            {
                localDOF = m_base[1]->GetNumModes();
                break;
            }
            case 2:
            {
                localDOF = 1;
                break;
            }
            default:
            {
                ASSERTL0(false, "eid must be between 0 and 2");
                break;
            }
        }
    }

    return localDOF;
}

void StdTriExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                 GetBasisType(1) == LibUtilities::eModified_B,
             "Expansion not of a proper type");

    int i, j;
    int cnt = 0;
    int nummodes0, nummodes1;
    int startvalue;
    if (outarray.size() != GetNcoeffs() - NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(GetNcoeffs() - NumBndryCoeffs());
    }

    nummodes0 = m_base[0]->GetNumModes();
    nummodes1 = m_base[1]->GetNumModes();

    startvalue = 2 * nummodes1;

    for (i = 0; i < nummodes0 - 2; i++)
    {
        for (j = 0; j < nummodes1 - 3 - i; j++)
        {
            outarray[cnt++] = startvalue + j;
        }
        startvalue += nummodes1 - 2 - i;
    }
}

void StdTriExp::v_GetBoundaryMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                 GetBasisType(1) == LibUtilities::eModified_B,
             "Expansion not of expected type");
    int i;
    int cnt;
    int nummodes0, nummodes1;
    int value;

    if (outarray.size() != NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
    }

    nummodes0 = m_base[0]->GetNumModes();
    nummodes1 = m_base[1]->GetNumModes();

    value = 2 * nummodes1 - 1;
    for (i = 0; i < value; i++)
    {
        outarray[i] = i;
    }
    cnt = value;

    for (i = 0; i < nummodes0 - 2; i++)
    {
        outarray[cnt++] = value;
        value += nummodes1 - 2 - i;
    }
}

void StdTriExp::v_GetTraceCoeffMap(const unsigned int eid,
                                   Array<OneD, unsigned int> &maparray)
{

    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eModified_B,
             "Mapping not defined for this type of basis");

    ASSERTL1(eid < 3, "eid must be between 0 and 2");

    int i;
    int order0   = m_base[0]->GetNumModes();
    int order1   = m_base[1]->GetNumModes();
    int numModes = (eid == 0) ? order0 : order1;

    if (maparray.size() != numModes)
    {
        maparray = Array<OneD, unsigned int>(numModes);
    }

    switch (eid)
    {
        case 0:
        {
            int cnt = 0;
            for (i = 0; i < numModes; cnt += order1 - i, ++i)
            {
                maparray[i] = cnt;
            }
            break;
        }
        case 1:
        {
            maparray[0] = order1;
            maparray[1] = 1;
            for (i = 2; i < numModes; i++)
            {
                maparray[i] = order1 - 1 + i;
            }
            break;
        }
        case 2:
        {
            for (i = 0; i < numModes; i++)
            {
                maparray[i] = i;
            }
            break;
        }
        default:
            ASSERTL0(false, "eid must be between 0 and 2");
            break;
    }
}

void StdTriExp::v_GetTraceInteriorToElementMap(
    const int eid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation edgeOrient)
{
    ASSERTL0(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eModified_B,
             "Mapping not defined for this type of basis");
    int i;
    const int nummodes1      = m_base[1]->GetNumModes();
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

    switch (eid)
    {
        case 0:
        {
            int cnt = 2 * nummodes1 - 1;
            for (i = 0; i < nEdgeIntCoeffs; cnt += nummodes1 - 2 - i, ++i)
            {
                maparray[i] = cnt;
            }
            break;
        }
        case 1:
        {
            for (i = 0; i < nEdgeIntCoeffs; i++)
            {
                maparray[i] = nummodes1 + 1 + i;
            }
            break;
        }
        case 2:
        {
            for (i = 0; i < nEdgeIntCoeffs; i++)
            {
                maparray[i] = 2 + i;
            }
            break;
        }
        default:
        {
            ASSERTL0(false, "eid must be between 0 and 2");
            break;
        }
    }

    if (edgeOrient == eBackwards)
    {
        for (i = 1; i < nEdgeIntCoeffs; i += 2)
        {
            signarray[i] = -1;
        }
    }
}

//---------------------------------------
// Wrapper functions
//---------------------------------------

DNekMatSharedPtr StdTriExp::v_GenMatrix(const StdMatrixKey &mkey)
{

    MatrixType mtype = mkey.GetMatrixType();

    DNekMatSharedPtr Mat;

    switch (mtype)
    {
        case ePhysInterpToEquiSpaced:
        {
            int nq0, nq1, nq;

            nq0 = m_base[0]->GetNumPoints();
            nq1 = m_base[1]->GetNumPoints();

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = max(nq0, nq1);
            }

            int neq = LibUtilities::StdTriData::getNumberOfCoefficients(nq, nq);
            Array<OneD, Array<OneD, NekDouble>> coords(neq);
            Array<OneD, NekDouble> coll(2);
            Array<OneD, DNekMatSharedPtr> I(2);
            Array<OneD, NekDouble> tmp(nq0);

            Mat     = MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1);
            int cnt = 0;

            for (int i = 0; i < nq; ++i)
            {
                for (int j = 0; j < nq - i; ++j, ++cnt)
                {
                    coords[cnt]    = Array<OneD, NekDouble>(2);
                    coords[cnt][0] = -1.0 + 2 * j / (NekDouble)(nq - 1);
                    coords[cnt][1] = -1.0 + 2 * i / (NekDouble)(nq - 1);
                }
            }

            for (int i = 0; i < neq - 1; ++i)
            {
                LocCoordToLocCollapsed(coords[i], coll);

                I[0] = m_base[0]->GetI(coll);
                I[1] = m_base[1]->GetI(coll + 1);

                // interpolate first coordinate direction
                for (int j = 0; j < nq1; ++j)
                {
                    NekDouble fac = (I[1]->GetPtr())[j];
                    Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                    Vmath::Vcopy(nq0, &tmp[0], 1,
                                 Mat->GetRawPtr() + j * nq0 * neq + i, neq);
                }
            }

            // evaluate top vertex as average of all interpolation
            // along collapsed coordinates at singular vertex
            coll[1] = 1.0;
            I[1]    = m_base[1]->GetI(coll + 1);

            // interpolate first coordinate direction
            for (int j = 0; j < nq1; ++j)
            {
                NekDouble fac = (I[1]->GetPtr())[j] / ((NekDouble)nq0);

                Vmath::Fill(nq0, fac,
                            Mat->GetRawPtr() + j * nq0 * neq + neq - 1, neq);
            }

            break;
        }
        case ePhysInterpToGLL:
        {
            int nq0, nq1, nq;

            nq0 = m_base[0]->GetNumPoints();
            nq1 = m_base[1]->GetNumPoints();

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = max(nq0, nq1);
            }

            int neq = LibUtilities::StdTriData::getNumberOfCoefficients(nq, nq);
            Array<OneD, NekDouble> coords(2);
            Array<OneD, NekDouble> coll(2);
            Array<OneD, DNekMatSharedPtr> I(2);
            Array<OneD, NekDouble> tmp(nq0);

            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1);

            const LibUtilities::PointsKey key(nq, LibUtilities::eNodalTriElec);

            Array<OneD, const NekDouble> x, y;
            LibUtilities::PointsManager()[key]->GetPoints(x, y);

            // set up nodal points in a similar mannger to equispaced
            // points startign at bottom and working upwards left to
            // right

            int row = 0;
            // First Vertex
            coords[0] = x[0];
            coords[1] = y[0];
            LocCoordToLocCollapsed(coords, coll);

            I[0] = m_base[0]->GetI(coll);
            I[1] = m_base[1]->GetI(coll + 1);

            // interpolate first coordinate direction
            for (int j = 0; j < nq1; ++j)
            {
                NekDouble fac = (I[1]->GetPtr())[j];
                Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                Vmath::Vcopy(nq0, &tmp[0], 1,
                             Mat->GetRawPtr() + j * nq0 * neq + row, neq);
            }
            row++;

            // First edge
            for (int i = 0; i < nq - 2; ++i)
            {
                coords[0] = x[3 + i];
                coords[1] = y[3 + i];
                LocCoordToLocCollapsed(coords, coll);

                I[0] = m_base[0]->GetI(coll);
                I[1] = m_base[1]->GetI(coll + 1);

                // interpolate first coordinate direction
                for (int j = 0; j < nq1; ++j)
                {
                    NekDouble fac = (I[1]->GetPtr())[j];
                    Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                    Vmath::Vcopy(nq0, &tmp[0], 1,
                                 Mat->GetRawPtr() + j * nq0 * neq + row, neq);
                }
                row++;
            }

            // Second Vertex
            coords[0] = x[1];
            coords[1] = y[1];
            LocCoordToLocCollapsed(coords, coll);

            I[0] = m_base[0]->GetI(coll);
            I[1] = m_base[1]->GetI(coll + 1);

            // interpolate first coordinate direction
            for (int j = 0; j < nq1; ++j)
            {
                NekDouble fac = (I[1]->GetPtr())[j];
                Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                Vmath::Vcopy(nq0, &tmp[0], 1,
                             Mat->GetRawPtr() + j * nq0 * neq + row, neq);
            }
            row++;

            int cnt = 0;
            // other  edges and interior
            for (int i = 0; i < nq - 2; ++i)
            {

                // edge 3 (reversing counter-clockwise ordering)
                coords[0] = x[3 * (nq - 1) - i - 1];
                coords[1] = y[3 * (nq - 1) - i - 1];
                // coords[0] = x[3 + 2*(nq-2)+i]; coords[1] = y[3 + 2*(nq-2)+i];
                LocCoordToLocCollapsed(coords, coll);

                I[0] = m_base[0]->GetI(coll);
                I[1] = m_base[1]->GetI(coll + 1);

                // interpolate first coordinate direction
                for (int j = 0; j < nq1; ++j)
                {
                    NekDouble fac = (I[1]->GetPtr())[j];
                    Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                    Vmath::Vcopy(nq0, &tmp[0], 1,
                                 Mat->GetRawPtr() + j * nq0 * neq + row, neq);
                }
                row++;

                for (int j = 0; j < nq - 3 - i; ++j)
                {
                    coords[0] = x[3 * (nq - 1) + cnt];
                    coords[1] = y[3 * (nq - 1) + cnt];
                    LocCoordToLocCollapsed(coords, coll);

                    I[0] = m_base[0]->GetI(coll);
                    I[1] = m_base[1]->GetI(coll + 1);

                    // interpolate first coordinate direction
                    for (int j = 0; j < nq1; ++j)
                    {
                        NekDouble fac = (I[1]->GetPtr())[j];
                        Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                        Vmath::Vcopy(nq0, &tmp[0], 1,
                                     Mat->GetRawPtr() + j * nq0 * neq + row,
                                     neq);
                    }
                    row++;
                    cnt++;
                }

                // edge 2
                coords[0] = x[3 + (nq - 2) + i];
                coords[1] = y[3 + (nq - 2) + i];
                LocCoordToLocCollapsed(coords, coll);

                I[0] = m_base[0]->GetI(coll);
                I[1] = m_base[1]->GetI(coll + 1);

                // interpolate first coordinate direction
                for (int j = 0; j < nq1; ++j)
                {
                    NekDouble fac = (I[1]->GetPtr())[j];
                    Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                    Vmath::Vcopy(nq0, &tmp[0], 1,
                                 Mat->GetRawPtr() + j * nq0 * neq + row, neq);
                }
                row++;
            }

            // evaluate top vertex as average of all interpolation
            // along collapsed coordinates at singular vertex
            coll[1] = 1.0;
            I[1]    = m_base[1]->GetI(coll + 1);

            // interpolate first coordinate direction
            for (int j = 0; j < nq1; ++j)
            {
                NekDouble fac = (I[1]->GetPtr())[j] / ((NekDouble)nq0);

                Vmath::Fill(nq0, fac,
                            Mat->GetRawPtr() + j * nq0 * neq + neq - 1, neq);
            }

            break;
        }
        case eEquiSpacedToPhys:
        {
            int nm0 = m_base[0]->GetNumPoints();
            int nm1 = m_base[1]->GetNumPoints();
            int neq;

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                neq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                neq = max(nm0, nm1);
            }

            // set up an exansion with the same number of modes as neq;
            LibUtilities::BasisKey ba(m_base[0]->GetBasisType(), neq,
                                      m_base[0]->GetPointsKey());
            LibUtilities::BasisKey bb(m_base[1]->GetBasisType(), neq,
                                      m_base[1]->GetPointsKey());

            StdTriExp Exp2D(ba, bb);
            int ncoeffs = Exp2D.GetNcoeffs();

            // Get hold of equispaced to coeff matrix
            ConstFactorMap cmap;
            cmap[eFactorConst] = neq;
            StdMatrixKey Ikey(eEquiSpacedToCoeffs, DetShapeType(), *this, cmap);
            DNekMatSharedPtr intmat = Exp2D.GetStdMatrix(Ikey);

            int nqtot = GetTotPoints();

            // generate a matrix
            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(nqtot, ncoeffs);
            NekDouble *ptr = Mat->GetRawPtr();

            Array<OneD, NekDouble> qmode(nqtot);

            // Get first mode at quadrature points
            Exp2D.FillMode(0, qmode);

            // first part of matrix-matrix multiply intiailising out matrix
            for (int j = 0; j < ncoeffs; ++j)
            {
                NekDouble val = (*intmat)(0, j);
                Vmath::Smul(nqtot, val, qmode.data(), 1, ptr + j * nqtot, 1);
            }

            for (int i = 1; i < ncoeffs; ++i)
            {
                // Get mode at quadrature points
                Exp2D.FillMode(i, qmode);

                for (int j = 0; j < ncoeffs; ++j)
                {
                    NekDouble val = (*intmat)(i, j);
                    Vmath::Svtvp(nqtot, val, qmode.data(), 1, ptr + j * nqtot,
                                 1, ptr + j * nqtot, 1);
                }
            }
        }
        break;
        default:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
            break;
        }
    }

    return Mat;
}

DNekMatSharedPtr StdTriExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

//---------------------------------------
// Operator evaluation functions
//---------------------------------------

void StdTriExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray,
                               const StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdTriExp::v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdTriExp::v_LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdTriExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void StdTriExp::v_WeakDerivMatrixOp(const int i,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void StdTriExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdTriExp::v_HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdTriExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                     const StdMatrixKey &mkey)
{
    int qa       = m_base[0]->GetNumPoints();
    int qb       = m_base[1]->GetNumPoints();
    int nmodes_a = m_base[0]->GetNumModes();
    int nmodes_b = m_base[1]->GetNumModes();

    // Declare orthogonal basis.
    LibUtilities::PointsKey pa(qa, m_base[0]->GetPointsType());
    LibUtilities::PointsKey pb(qb, m_base[1]->GetPointsType());

    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, nmodes_a, pa);
    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B, nmodes_b, pb);
    StdTriExp OrthoExp(Ba, Bb);

    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());

    // project onto physical space.
    OrthoExp.FwdTrans(array, orthocoeffs);

    if (mkey.ConstFactorExists(
            eFactorSVVPowerKerDiffCoeff)) // Rodrigo's power kern
    {
        NekDouble cutoff = mkey.GetConstFactor(eFactorSVVCutoffRatio);
        NekDouble SvvDiffCoeff =
            mkey.GetConstFactor(eFactorSVVPowerKerDiffCoeff) *
            mkey.GetConstFactor(eFactorSVVDiffCoeff);

        int cnt = 0;
        for (int j = 0; j < nmodes_a; ++j)
        {
            for (int k = 0; k < nmodes_b - j; ++k, ++cnt)
            {
                NekDouble fac = std::max(
                    pow((1.0 * j) / (nmodes_a - 1), cutoff * nmodes_a),
                    pow((1.0 * k) / (nmodes_b - 1), cutoff * nmodes_b));

                orthocoeffs[cnt] *= (SvvDiffCoeff * fac);
            }
        }
    }
    else if (mkey.ConstFactorExists(
                 eFactorSVVDGKerDiffCoeff)) // Rodrigo/mansoor's DG kernel
    {
        NekDouble SvvDiffCoeff = mkey.GetConstFactor(eFactorSVVDGKerDiffCoeff) *
                                 mkey.GetConstFactor(eFactorSVVDiffCoeff);
        int max_ab = max(nmodes_a - kSVVDGFiltermodesmin,
                         nmodes_b - kSVVDGFiltermodesmin);
        max_ab     = max(max_ab, 0);
        max_ab     = min(max_ab, kSVVDGFiltermodesmax - kSVVDGFiltermodesmin);

        int cnt = 0;
        for (int j = 0; j < nmodes_a; ++j)
        {
            for (int k = 0; k < nmodes_b - j; ++k, ++cnt)
            {
                int maxjk = max(j, k);
                maxjk     = min(maxjk, kSVVDGFiltermodesmax - 1);

                orthocoeffs[cnt] *= SvvDiffCoeff * kSVVDGFilter[max_ab][maxjk];
            }
        }
    }
    else
    {
        NekDouble SvvDiffCoeff = mkey.GetConstFactor(eFactorSVVDiffCoeff);

        int cutoff = (int)(mkey.GetConstFactor(eFactorSVVCutoffRatio) *
                           min(nmodes_a, nmodes_b));

        NekDouble epsilon = 1.0;
        int nmodes        = min(nmodes_a, nmodes_b);

        int cnt = 0;

        // apply SVV filter (JEL)
        for (int j = 0; j < nmodes_a; ++j)
        {
            for (int k = 0; k < nmodes_b - j; ++k)
            {
                if (j + k >= cutoff)
                {
                    orthocoeffs[cnt] *=
                        (SvvDiffCoeff *
                         exp(-(j + k - nmodes) * (j + k - nmodes) /
                             ((NekDouble)((j + k - cutoff + epsilon) *
                                          (j + k - cutoff + epsilon)))));
                }
                else
                {
                    orthocoeffs[cnt] *= 0.0;
                }
                cnt++;
            }
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdTriExp::v_ReduceOrderCoeffs(int numMin,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    int n_coeffs = inarray.size();
    int nquad0   = m_base[0]->GetNumPoints();
    int nquad1   = m_base[1]->GetNumPoints();
    Array<OneD, NekDouble> coeff(n_coeffs);
    Array<OneD, NekDouble> coeff_tmp(n_coeffs, 0.0);
    Array<OneD, NekDouble> tmp;
    Array<OneD, NekDouble> tmp2;
    int nqtot = nquad0 * nquad1;
    Array<OneD, NekDouble> phys_tmp(nqtot, 0.0);

    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();
    int numMin2 = nmodes0;
    int i;

    const LibUtilities::PointsKey Pkey0(nmodes0,
                                        LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey Pkey1(nmodes1,
                                        LibUtilities::eGaussLobattoLegendre);

    LibUtilities::BasisKey b0(m_base[0]->GetBasisType(), nmodes0, Pkey0);
    LibUtilities::BasisKey b1(m_base[1]->GetBasisType(), nmodes1, Pkey1);

    LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A, nmodes0, Pkey0);
    LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_B, nmodes1, Pkey1);

    StdRegions::StdTriExpSharedPtr m_OrthoTriExp;
    StdRegions::StdTriExpSharedPtr m_TriExp;

    m_TriExp = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(b0, b1);
    m_OrthoTriExp = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(
        bortho0, bortho1);

    m_TriExp->BwdTrans(inarray, phys_tmp);
    m_OrthoTriExp->FwdTrans(phys_tmp, coeff);

    for (i = 0; i < n_coeffs; i++)
    {
        if (i == numMin)
        {
            coeff[i] = 0.0;
            numMin += numMin2 - 1;
            numMin2 -= 1.0;
        }
    }

    m_OrthoTriExp->BwdTrans(coeff, phys_tmp);
    m_TriExp->FwdTrans(phys_tmp, outarray);
}

//---------------------------------------
// Private helper functions
//---------------------------------------

void StdTriExp::v_GetSimplexEquiSpacedConnectivity(
    Array<OneD, int> &conn, [[maybe_unused]] bool standard)
{
    int np1 = m_base[0]->GetNumPoints();
    int np2 = m_base[1]->GetNumPoints();
    int np  = max(np1, np2);

    conn = Array<OneD, int>(3 * (np - 1) * (np - 1));

    int row   = 0;
    int rowp1 = 0;
    int cnt   = 0;
    for (int i = 0; i < np - 1; ++i)
    {
        rowp1 += np - i;
        for (int j = 0; j < np - i - 2; ++j)
        {
            conn[cnt++] = row + j;
            conn[cnt++] = row + j + 1;
            conn[cnt++] = rowp1 + j;

            conn[cnt++] = rowp1 + j + 1;
            conn[cnt++] = rowp1 + j;
            conn[cnt++] = row + j + 1;
        }

        conn[cnt++] = row + np - i - 2;
        conn[cnt++] = row + np - i - 1;
        conn[cnt++] = rowp1 + np - i - 2;

        row += np - i;
    }
}
} // namespace Nektar::StdRegions
