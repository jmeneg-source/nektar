///////////////////////////////////////////////////////////////////////////////
//
// File: StdPyrExp.cpp
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
// Description: pyramidic routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <StdRegions/StdPyrExp.h>
#include <iomanip>

using namespace std;
#include <LibUtilities/BasicUtils/NekInline.hpp>
#include <StdRegions/Operators/SwitchLevel1.h>
#include <StdRegions/Operators/SwitchLevel2.h>

namespace Nektar::StdRegions
{
// Declaretion of scalar routine
using vec_t = tinysimd::scalarT<double>;
#include <StdRegions/Operators/BwdTransSumFacStdKernels.hpp>
#include <StdRegions/Operators/IProductWRTBaseSumFacStdKernels.hpp>

StdPyrExp::StdPyrExp(const LibUtilities::BasisKey &Ba,
                     const LibUtilities::BasisKey &Bb,
                     const LibUtilities::BasisKey &Bc)
    : StdExpansion(LibUtilities::StdPyrData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                   3, Ba, Bb, Bc),
      StdExpansion3D(LibUtilities::StdPyrData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                     Ba, Bb, Bc)
{

    ASSERTL0(Ba.GetNumModes() <= Bc.GetNumModes(),
             "order in 'a' direction is higher "
             "than order in 'c' direction");
    ASSERTL0(Bb.GetNumModes() <= Bc.GetNumModes(),
             "order in 'b' direction is higher "
             "than order in 'c' direction");
    ASSERTL1(Bc.GetBasisType() == LibUtilities::eModifiedPyr_C ||
                 Bc.GetBasisType() == LibUtilities::eOrthoPyr_C,
             "Expected basis type in 'c' direction to be ModifiedPyr_C or "
             "OrthoPyr_C");

    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());
    m_weights.push_back(m_base[1]->GetW());

    StdFacKey w2key(eWeights2, Bc);
    // get weights[2] from manager where points are rescaled
    m_weights.push_back(GetStdFac(w2key));
}

//---------------------------------------
// Differentiation/integration Methods
//---------------------------------------
/**
 * \brief Calculate the derivative of the physical points
 *
 * The derivative is evaluated at the nodal physical points.
 * Derivatives with respect to the local Cartesian coordinates.
 *
 * \f$\begin{Bmatrix} \frac {\partial} {\partial \xi_1} \\ \frac
 * {\partial} {\partial \xi_2} \\ \frac {\partial} {\partial \xi_3}
 * \end{Bmatrix} = \begin{Bmatrix} \frac 2 {(1-\eta_3)} \frac \partial
 * {\partial \bar \eta_1} \\ \frac {\partial} {\partial \xi_2} \ \
 * \frac {(1 + \bar \eta_1)} {(1 - \eta_3)} \frac \partial {\partial
 * \bar \eta_1} + \frac {\partial} {\partial \eta_3} \end{Bmatrix}\f$
 */
void StdPyrExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &u_physical,
                               Array<OneD, NekDouble> &out_dxi1,
                               Array<OneD, NekDouble> &out_dxi2,
                               Array<OneD, NekDouble> &out_dxi3)
{
    // PhysDerivative implementation based on Spen's book page 152.
    int Qx = m_base[0]->GetNumPoints();
    int Qy = m_base[1]->GetNumPoints();
    int Qz = m_base[2]->GetNumPoints();

    Array<OneD, NekDouble> dEta_bar1(Qx * Qy * Qz, 0.0);
    Array<OneD, NekDouble> dXi2(Qx * Qy * Qz, 0.0);
    Array<OneD, NekDouble> dEta3(Qx * Qy * Qz, 0.0);
    PhysTensorDeriv(u_physical, dEta_bar1, dXi2, dEta3);

    Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
    eta_x = m_base[0]->GetZ();
    eta_y = m_base[1]->GetZ();
    eta_z = m_base[2]->GetZ();

    int i, j, k, n;

    if (out_dxi1.size() > 0)
    {
        for (k = 0, n = 0; k < Qz; ++k)
        {
            NekDouble fac = 2.0 / (1.0 - eta_z[k]);
            for (j = 0; j < Qy; ++j)
            {
                for (i = 0; i < Qx; ++i, ++n)
                {
                    out_dxi1[n] = fac * dEta_bar1[n];
                }
            }
        }
    }

    if (out_dxi2.size() > 0)
    {
        for (k = 0, n = 0; k < Qz; ++k)
        {
            NekDouble fac = 2.0 / (1.0 - eta_z[k]);
            for (j = 0; j < Qy; ++j)
            {
                for (i = 0; i < Qx; ++i, ++n)
                {
                    out_dxi2[n] = fac * dXi2[n];
                }
            }
        }
    }

    if (out_dxi3.size() > 0)
    {
        for (k = 0, n = 0; k < Qz; ++k)
        {
            NekDouble fac = 1.0 / (1.0 - eta_z[k]);
            for (j = 0; j < Qy; ++j)
            {
                NekDouble fac1 = (1.0 + eta_y[j]);
                for (i = 0; i < Qx; ++i, ++n)
                {
                    out_dxi3[n] = (1.0 + eta_x[i]) * fac * dEta_bar1[n] +
                                  fac1 * fac * dXi2[n] + dEta3[n];
                }
            }
        }
    }
}

//---------------------------------------
// Transforms
//---------------------------------------

/**
 * \brief Backward transformation is evaluated at the quadrature
 * points.
 *
 * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
 * u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
 *
 * Backward transformation is three dimensional tensorial expansion
 *
 * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
 *  (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
 *  \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k})
 *  \rbrace} \rbrace}. \f$
 *
 * And sumfactorizing step of the form is as:\ \ \f$ f_{pqr}
 * (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c
 * (\xi_{3k}),\\ g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y}
 * \psi_{p}^a (\xi_{2j}) f_{pqr} (\xi_{3k}),\\ u(\xi_{1i}, \xi_{2j},
 * \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p}
 * (\xi_{2j}, \xi_{3k}).  \f$
 **/
void StdPyrExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &outarray)
{
    const Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
    const Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();
    const Array<OneD, const NekDouble> base2 = m_base[2]->GetBdata();

    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();
    int nmodes2 = m_base[2]->GetNumModes();

    bool isModified = (m_base[0]->GetBasisType() == LibUtilities::eModified_A);

    std::vector<vec_t, tinysimd::allocator<vec_t>> wsp0(nmodes0 * nmodes1),
        wsp1(nmodes0);

    // Switch statment using boost_pp and macros. This unfolls intwo a
    // nested swtich statement where the outer swtich statement runs
    // from SMIN to SMAX for modal order and the inner switch
    // statemets run from the outer value of the case to 2*SMAX for
    // the quadrature order. If you want to see it unwrapped compile
    // in verbose mode and add --preprocess to the c++ command.
    // Default case
#undef BWDTRANS_DEF
#define BWDTRANS_DEF                                                           \
    BwdTransPyrKernel(nmodes0, nmodes1, nmodes2, nquad0, nquad1, nquad2,       \
                      isModified, (const vec_t *)base0.data(),                 \
                      (const vec_t *)base1.data(),                             \
                      (const vec_t *)base2.data(), wsp0.data(), wsp1.data(),   \
                      (const vec_t *)inarray.data(), (vec_t *)outarray.data())

    // Inner loop case over quarature points
#undef BWDTRANS_Q
#define BWDTRANS_Q(r, i)                                                       \
    case NQ(i):                                                                \
        BwdTransPyrKernel(                                                     \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ_M1(i), isModified,           \
            (const vec_t *)base0.data(), (const vec_t *)base1.data(),          \
            (const vec_t *)base2.data(), wsp0.data(), wsp1.data(),             \
            (const vec_t *)inarray.data(), (vec_t *)outarray.data());          \
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
    if ((nmodes0 == nmodes1) && (nmodes1 == nmodes2) && (nquad0 == nquad1) &&
        (nquad1 == nquad2 + 1))
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
 *  @param base2 - An array containing the values of the basis in the
 *  2-direction at the quarature poitns
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
void StdPyrExp::v_IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &base1,
    const Array<OneD, const NekDouble> &base2,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, NekDouble> &jac,
    const bool Deformed, [[maybe_unused]] bool CollDir0,
    [[maybe_unused]] bool CollDir1, [[maybe_unused]] bool CollDir2)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    int order0 = m_base[0]->GetNumModes();
    int order1 = m_base[1]->GetNumModes();
    int order2 = m_base[2]->GetNumModes();

    const bool isModified =
        (m_base[0]->GetBasisType() == LibUtilities::eModified_A);

    std::vector<vec_t, tinysimd::allocator<vec_t>> wsp0(nquad1 * nquad2),
        wsp1(nquad2);

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
    IProductPyrKernel<false, false, true>(                                     \
        order0, order1, order2, nquad0, nquad1, nquad2, isModified,            \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)base2.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(),                                    \
        (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductPyrKernel<false, false, true>(                                 \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ_M1(i), isModified,           \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)base2.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(),                                \
            (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)wsp1.data(),                        \
            (vec_t *)outarray.data());                                         \
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
        if ((order0 == order1) && (order1 == order2) && (nquad0 == nquad1) &&
            (nquad1 == nquad2 + 1))
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
    IProductPyrKernel<false, false, false>(                                    \
        order0, order1, order2, nquad0, nquad1, nquad2, isModified,            \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)base2.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(),                                    \
        (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductPyrKernel<false, false, false>(                                \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ_M1(i), isModified,           \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)base2.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(),                                \
            (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)wsp1.data(),                        \
            (vec_t *)outarray.data());                                         \
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
        if ((order0 == order1) && (order1 == order2) && (nquad0 == nquad1) &&
            (nquad1 == nquad2 + 1))
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

void StdPyrExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    StdFacKey fackey2(eTwoOverOneMinusZ2, m_base[2]->GetBasisKey());
    Array<OneD, const NekDouble> gfac2 = GetStdFac(fackey2);

    Array<OneD, NekDouble> tmp0(nquad0 * nquad1 * nquad2);

    // Derivative in first/second direction is always scaled as follows
    const int nq01 = nquad0 * nquad1;
    for (int i = 0; i < nquad2; ++i)
    {
        Vmath::Smul(nq01, gfac2[i], &inarray[0] + i * nq01, 1,
                    &tmp0[0] + i * nq01, 1);
    }

    const Array<OneD, const NekDouble> one(1, 1.0);

    switch (dir)
    {
        case 0:
        {
            v_IProductWRTBaseKernel(
                m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                m_base[2]->GetBdata(), tmp0, outarray, one, false);
        }
        break;
        case 1:
        {
            v_IProductWRTBaseKernel(
                m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                m_base[2]->GetBdata(), tmp0, outarray, one, false);
        }
        break;
        case 2:
        {
            StdFacKey fackey0(eHalfMultOnePlusZ0, m_base[0]->GetBasisKey());
            Array<OneD, const NekDouble> gfac0 = GetStdFac(fackey0);
            StdFacKey fackey1(eHalfMultOnePlusZ1, m_base[1]->GetBasisKey());
            Array<OneD, const NekDouble> gfac1 = GetStdFac(fackey1);

            Array<OneD, NekDouble> tmp3(m_ncoeffs);
            Array<OneD, NekDouble> tmp4(m_ncoeffs);

            // Scale eta_1 derivative by gfac0
            for (int i = 0; i < nquad1 * nquad2; ++i)
            {
                Vmath::Vmul(nquad0, tmp0.data() + i * nquad0, 1, gfac0.data(),
                            1, tmp0.data() + i * nquad0, 1);
            }
            v_IProductWRTBaseKernel(
                m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                m_base[2]->GetBdata(), tmp0, tmp3, one, false);

            // Scale eta_2 derivative by gfac1*gfac2
            for (int i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nq01, gfac2[i], &inarray[0] + i * nq01, 1,
                            &tmp0[0] + i * nq01, 1);
            }

            for (int i = 0; i < nquad1 * nquad2; ++i)
            {
                Vmath::Smul(nquad0, gfac1[i % nquad1], &tmp0[0] + i * nquad0, 1,
                            &tmp0[0] + i * nquad0, 1);
            }

            v_IProductWRTBaseKernel(
                m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                m_base[2]->GetBdata(), tmp0, tmp4, one, false);

            v_IProductWRTBaseKernel(
                m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                m_base[2]->GetDbdata(), inarray, outarray, one, false);

            Vmath::Vadd(m_ncoeffs, &tmp3[0], 1, &outarray[0], 1, &outarray[0],
                        1);
            Vmath::Vadd(m_ncoeffs, &tmp4[0], 1, &outarray[0], 1, &outarray[0],
                        1);
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

void StdPyrExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble> &xi,
                                         Array<OneD, NekDouble> &eta)
{
    NekDouble d2 = 1.0 - xi[2];
    if (fabs(d2) < NekConstants::kNekZeroTol)
    {
        if (d2 >= 0.)
        {
            d2 = NekConstants::kNekZeroTol;
        }
        else
        {
            d2 = -NekConstants::kNekZeroTol;
        }
    }
    eta[2] = xi[2]; // eta_z = xi_z
    eta[1] = 2.0 * (1.0 + xi[1]) / d2 - 1.0;
    eta[0] = 2.0 * (1.0 + xi[0]) / d2 - 1.0;
}

void StdPyrExp::v_LocCollapsedToLocCoord(
    const Array<OneD, const NekDouble> &eta, Array<OneD, NekDouble> &xi)
{
    xi[0] = (1.0 + eta[0]) * (1.0 - eta[2]) * 0.5 - 1.0;
    xi[1] = (1.0 + eta[1]) * (1.0 - eta[2]) * 0.5 - 1.0;
    xi[2] = eta[2];
}

void StdPyrExp::v_GetCoords(Array<OneD, NekDouble> &xi_x,
                            Array<OneD, NekDouble> &xi_y,
                            Array<OneD, NekDouble> &xi_z)
{
    Array<OneD, const NekDouble> etaBar_x = m_base[0]->GetZ();
    Array<OneD, const NekDouble> eta_y    = m_base[1]->GetZ();
    Array<OneD, const NekDouble> eta_z    = m_base[2]->GetZ();
    int Qx                                = GetNumPoints(0);
    int Qy                                = GetNumPoints(1);
    int Qz                                = GetNumPoints(2);

    // Convert collapsed coordinates into cartesian coordinates: eta --> xi
    for (int k = 0; k < Qz; ++k)
    {
        for (int j = 0; j < Qy; ++j)
        {
            for (int i = 0; i < Qx; ++i)
            {
                int s = i + Qx * (j + Qy * k);

                xi_z[s] = eta_z[k];
                xi_y[s] = (1.0 + eta_y[j]) * (1.0 - eta_z[k]) / 2.0 - 1.0;
                xi_x[s] = (1.0 + etaBar_x[i]) * (1.0 - eta_z[k]) / 2.0 - 1.0;
            }
        }
    }
}

NekDouble StdPyrExp::v_PhysEvaluateBasis(
    const Array<OneD, const NekDouble> &coords, int mode)
{
    Array<OneD, NekDouble> coll(3);
    LocCoordToLocCollapsed(coords, coll);

    const int nm0 = m_base[0]->GetNumModes();
    const int nm1 = m_base[1]->GetNumModes();
    const int nm2 = m_base[2]->GetNumModes();

    int mode0 = 0, mode1 = 0, mode2 = 0, cnt = 0;

    bool found = false;
    for (mode0 = 0; mode0 < nm0; ++mode0)
    {
        for (mode1 = 0; mode1 < nm1; ++mode1)
        {
            int maxpq = max(mode0, mode1);
            for (mode2 = 0; mode2 < nm2 - maxpq; ++mode2, ++cnt)
            {
                if (cnt == mode)
                {
                    found = true;
                    break;
                }
            }

            if (found)
            {
                break;
            }
        }

        if (found)
        {
            break;
        }

        for (int j = nm1; j < nm2; ++j)
        {
            int ijmax = max(mode0, j);
            mode2 += nm2 - ijmax;
        }
    }

    if (mode == 1 && m_base[0]->GetBasisType() == LibUtilities::eModified_A)
    {
        return StdExpansion::BaryEvaluateBasis<2>(coll[2], 1);
    }
    else
    {
        return StdExpansion::BaryEvaluateBasis<0>(coll[0], mode0) *
               StdExpansion::BaryEvaluateBasis<1>(coll[1], mode1) *
               StdExpansion::BaryEvaluateBasis<2>(coll[2], mode2);
    }
}

NekDouble StdPyrExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    // Collapse coordinates
    Array<OneD, NekDouble> coll(3, 0.0);
    LocCoordToLocCollapsed(coord, coll);

    // If near singularity do the old interpolation matrix method
    if ((1 - coll[2]) < 1e-5)
    {
        int totPoints = GetTotPoints();
        Array<OneD, NekDouble> EphysDeriv0(totPoints), EphysDeriv1(totPoints),
            EphysDeriv2(totPoints);
        v_PhysDeriv(inarray, EphysDeriv0, EphysDeriv1, EphysDeriv2);

        Array<OneD, DNekMatSharedPtr> I(3);
        I[0] = GetBase()[0]->GetI(coll);
        I[1] = GetBase()[1]->GetI(coll + 1);
        I[2] = GetBase()[2]->GetI(coll + 2);

        firstOrderDerivs[0] = PhysEvaluate(I, EphysDeriv0);
        firstOrderDerivs[1] = PhysEvaluate(I, EphysDeriv1);
        firstOrderDerivs[2] = PhysEvaluate(I, EphysDeriv2);
        return PhysEvaluate(I, inarray);
    }

    std::array<NekDouble, 3> interDeriv;
    NekDouble val = StdExpansion3D::BaryTensorDeriv(coll, inarray, interDeriv);

    NekDouble fac = 2.0 / (1.0 - coll[2]);

    firstOrderDerivs[0] = fac * interDeriv[0];
    firstOrderDerivs[1] = fac * interDeriv[1];
    firstOrderDerivs[2] = ((1.0 + coll[0]) / (1.0 - coll[2])) * interDeriv[0] +
                          ((1.0 + coll[1]) / (1.0 - coll[2])) * interDeriv[1] +
                          interDeriv[2];

    return val;
}

void StdPyrExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp(m_ncoeffs, 0.0);
    tmp[mode] = 1.0;
    v_BwdTrans(tmp, outarray);
}

void StdPyrExp::v_GetTraceNumModes(const int fid, int &numModes0,
                                   int &numModes1, Orientation faceOrient)
{
    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};
    switch (fid)
    {
        // quad
        case 0:
        {
            numModes0 = nummodes[0];
            numModes1 = nummodes[1];
        }
        break;
        case 1:
        case 3:
        {
            numModes0 = nummodes[0];
            numModes1 = nummodes[2];
        }
        break;
        case 2:
        case 4:
        {
            numModes0 = nummodes[1];
            numModes1 = nummodes[2];
        }
        break;
    }

    if (faceOrient >= 9)
    {
        std::swap(numModes0, numModes1);
    }
}

//---------------------------------------
// Helper functions
//---------------------------------------

int StdPyrExp::v_GetNverts() const
{
    return 5;
}

int StdPyrExp::v_GetNedges() const
{
    return 8;
}

int StdPyrExp::v_GetNtraces() const
{
    return 5;
}

LibUtilities::ShapeType StdPyrExp::v_DetShapeType() const
{
    return LibUtilities::ePyramid;
}

int StdPyrExp::v_NumBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModifiedPyr_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes();
    int Q = m_base[1]->GetNumModes();
    int R = m_base[2]->GetNumModes();

    return LibUtilities::StdPyrData::getNumberOfBndCoefficients(P, Q, R);
}

int StdPyrExp::v_NumDGBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModifiedPyr_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes() - 1;
    int Q = m_base[1]->GetNumModes() - 1;
    int R = m_base[2]->GetNumModes() - 1;

    return (P + 1) * (Q + 1)                    // 1 rect. face on base
           + 2 * (R + 1) + P * (1 + 2 * R - P)  // 2 tri. (P,R) faces
           + 2 * (R + 1) + Q * (1 + 2 * R - Q); // 2 tri. (Q,R) faces
}

int StdPyrExp::v_GetTraceNcoeffs(const int i) const
{
    ASSERTL2(i >= 0 && i <= 4, "face id is out of range");

    if (i == 0)
    {
        return GetBasisNumModes(0) * GetBasisNumModes(1);
    }
    else if (i == 1 || i == 3)
    {
        int P = GetBasisNumModes(0) - 1, Q = GetBasisNumModes(2) - 1;
        return Q + 1 + (P * (1 + 2 * Q - P)) / 2;
    }
    else
    {
        int P = GetBasisNumModes(1) - 1, Q = GetBasisNumModes(2) - 1;
        return Q + 1 + (P * (1 + 2 * Q - P)) / 2;
    }
}

int StdPyrExp::v_GetTraceIntNcoeffs(const int i) const
{
    ASSERTL2(i >= 0 && i <= 4, "face id is out of range");

    int P = m_base[0]->GetNumModes() - 1;
    int Q = m_base[1]->GetNumModes() - 1;
    int R = m_base[2]->GetNumModes() - 1;

    if (i == 0)
    {
        return (P - 1) * (Q - 1);
    }
    else if (i == 1 || i == 3)
    {
        return (P - 1) * (2 * (R - 1) - (P - 1) - 1) / 2;
    }
    else
    {
        return (Q - 1) * (2 * (R - 1) - (Q - 1) - 1) / 2;
    }
}

int StdPyrExp::v_GetTraceNumPoints(const int i) const
{
    ASSERTL2(i >= 0 && i <= 4, "face id is out of range");

    if (i == 0)
    {
        return m_base[0]->GetNumPoints() * m_base[1]->GetNumPoints();
    }
    else if (i == 1 || i == 3)
    {
        return m_base[0]->GetNumPoints() * m_base[2]->GetNumPoints();
    }
    else
    {
        return m_base[1]->GetNumPoints() * m_base[2]->GetNumPoints();
    }
}

int StdPyrExp::v_GetEdgeNcoeffs(const int i) const
{
    ASSERTL2(i >= 0 && i <= 7, "edge id is out of range");

    if (i == 0 || i == 2)
    {
        return GetBasisNumModes(0);
    }
    else if (i == 1 || i == 3)
    {
        return GetBasisNumModes(1);
    }
    else
    {
        return GetBasisNumModes(2);
    }
}

const LibUtilities::BasisKey StdPyrExp::v_GetTraceBasisKey(const int i,
                                                           const int k,
                                                           bool UseGLL) const
{
    ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
    ASSERTL2(k >= 0 && k <= 1, "basis key id is out of range");

    switch (i)
    {
        case 0:
        {
            return EvaluateQuadFaceBasisKey(k, m_base[k]);
        }
        case 1:
        case 3:
        {
            return EvaluateTriFaceBasisKey(k, m_base[2 * k], UseGLL);
        }
        case 2:
        case 4:
        {
            return EvaluateTriFaceBasisKey(k, m_base[k + 1], UseGLL);
        }
    }

    // Should never get here.
    return LibUtilities::NullBasisKey;
}

int StdPyrExp::v_CalcNumberOfCoefficients(
    const std::vector<unsigned int> &nummodes, int &modes_offset)
{
    int nmodes = LibUtilities::StdPyrData::getNumberOfCoefficients(
        nummodes[modes_offset], nummodes[modes_offset + 1],
        nummodes[modes_offset + 2]);

    modes_offset += 3;
    return nmodes;
}

int StdPyrExp::v_GetVertexMap(int vId, bool useCoeffPacking)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eModifiedPyr_C,
             "Mapping not defined for this type of basis");

    int l = 0;

    if (useCoeffPacking == true) // follow packing of coefficients i.e q,r,p
    {
        switch (vId)
        {
            case 0:
                l = GetMode(0, 0, 0);
                break;
            case 1:
                l = GetMode(0, 0, 1);
                break;
            case 2:
                l = GetMode(0, 1, 0);
                break;
            case 3:
                l = GetMode(1, 0, 0);
                break;
            case 4:
                l = GetMode(1, 1, 0);
                break;
            default:
                ASSERTL0(false, "local vertex id must be between 0 and 4");
        }
    }
    else
    {
        switch (vId)
        {
            case 0:
                l = GetMode(0, 0, 0);
                break;
            case 1:
                l = GetMode(1, 0, 0);
                break;
            case 2:
                l = GetMode(1, 1, 0);
                break;
            case 3:
                l = GetMode(0, 1, 0);
                break;
            case 4:
                l = GetMode(0, 0, 1);
                break;
            default:
                ASSERTL0(false, "local vertex id must be between 0 and 4");
        }
    }

    return l;
}

void StdPyrExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModifiedPyr_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes() - 1, p;
    int Q = m_base[1]->GetNumModes() - 1, q;
    int R = m_base[2]->GetNumModes() - 1, r;

    int nIntCoeffs = m_ncoeffs - NumBndryCoeffs();

    if (outarray.size() != nIntCoeffs)
    {
        outarray = Array<OneD, unsigned int>(nIntCoeffs);
    }

    int idx = 0;

    // Loop over all interior modes.
    for (p = 2; p <= P; ++p)
    {
        for (q = 2; q <= Q; ++q)
        {
            int maxpq = max(p, q);
            for (r = 1; r <= R - maxpq; ++r)
            {
                outarray[idx++] = GetMode(p, q, r);
            }
        }
    }
}

void StdPyrExp::v_GetBoundaryMap(Array<OneD, unsigned int> &maparray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModifiedPyr_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P   = m_base[0]->GetNumModes() - 1, p;
    int Q   = m_base[1]->GetNumModes() - 1, q;
    int R   = m_base[2]->GetNumModes() - 1, r;
    int idx = 0;

    int nBnd = NumBndryCoeffs();

    if (maparray.size() != nBnd)
    {
        maparray = Array<OneD, unsigned int>(nBnd);
    }

    // Loop over all boundary modes (in ascending order).
    for (p = 0; p <= P; ++p)
    {
        // First two q-r planes are entirely boundary modes.
        if (p <= 1)
        {
            for (q = 0; q <= Q; ++q)
            {
                int maxpq = max(p, q);
                for (r = 0; r <= R - maxpq; ++r)
                {
                    maparray[idx++] = GetMode(p, q, r);
                }
            }
        }
        else
        {
            // Remaining q-r planes contain boundary modes on the two
            // front and back sides and edges 0 2.
            for (q = 0; q <= Q; ++q)
            {
                if (q <= 1)
                {
                    for (r = 0; r <= R - p; ++r)
                    {
                        maparray[idx++] = GetMode(p, q, r);
                    }
                }
                else
                {
                    maparray[idx++] = GetMode(p, q, 0);
                }
            }
        }
    }
}

void StdPyrExp::v_GetTraceCoeffMap(const unsigned int fid,
                                   Array<OneD, unsigned int> &maparray)
{
    ASSERTL1(GetBasisType(0) == GetBasisType(1),
             "Method only implemented if BasisType is identical"
             "in x and y directions");
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                 GetBasisType(2) == LibUtilities::eModifiedPyr_C,
             "Method only implemented for Modified_A BasisType"
             "(x and y direction) and ModifiedPyr_C BasisType (z "
             "direction)");

    int p, q, r, P = 0, Q = 0, idx = 0;

    int order0 = m_base[0]->GetNumModes();
    int order1 = m_base[1]->GetNumModes();
    int order2 = m_base[2]->GetNumModes();

    switch (fid)
    {
        case 0:
            P = order0;
            Q = order1;
            break;
        case 1:
        case 3:
            P = order0;
            Q = order2;
            break;
        case 2:
        case 4:
            P = order1;
            Q = order2;
            break;
        default:
            ASSERTL0(false, "fid must be between 0 and 4");
    }

    if (maparray.size() != P * Q)
    {
        maparray = Array<OneD, unsigned int>(P * Q);
    }

    // Set up ordering inside each 2D face. Also for triangular faces,
    // populate signarray.
    switch (fid)
    {
        case 0: // Bottom quad

            for (q = 0; q < Q; ++q)
            {
                for (p = 0; p < P; ++p)
                {
                    maparray[q * P + p] = GetMode(p, q, 0);
                }
            }
            break;

        case 1: // Front triangle
            for (p = 0; p < P; ++p)
            {
                for (r = 0; r < Q - p; ++r)
                {
                    maparray[idx++] = GetMode(p, 0, r);
                }
            }
            break;

        case 2: // Right triangle
            maparray[idx++] = GetMode(1, 0, 0);
            maparray[idx++] = GetMode(0, 0, 1);
            for (r = 1; r < Q - 1; ++r)
            {
                maparray[idx++] = GetMode(1, 0, r);
            }

            for (q = 1; q < P; ++q)
            {
                for (r = 0; r < Q - q; ++r)
                {
                    maparray[idx++] = GetMode(1, q, r);
                }
            }
            break;

        case 3: // Rear triangle
            maparray[idx++] = GetMode(0, 1, 0);
            maparray[idx++] = GetMode(0, 0, 1);
            for (r = 1; r < Q - 1; ++r)
            {
                maparray[idx++] = GetMode(0, 1, r);
            }

            for (p = 1; p < P; ++p)
            {
                for (r = 0; r < Q - p; ++r)
                {
                    maparray[idx++] = GetMode(p, 1, r);
                }
            }
            break;

        case 4: // Left triangle
            for (q = 0; q < P; ++q)
            {
                for (r = 0; r < Q - q; ++r)
                {
                    maparray[idx++] = GetMode(0, q, r);
                }
            }
            break;

        default:
            ASSERTL0(false, "Face to element map unavailable.");
    }
}

void StdPyrExp::v_GetElmtTraceToTraceMap(const unsigned int fid,
                                         Array<OneD, unsigned int> &maparray,
                                         Array<OneD, int> &signarray,
                                         Orientation faceOrient, int P, int Q)
{
    ASSERTL1(GetBasisType(0) == GetBasisType(1),
             "Method only implemented if BasisType is identical"
             "in x and y directions");
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                 GetBasisType(2) == LibUtilities::eModifiedPyr_C,
             "Method only implemented for Modified_A BasisType"
             "(x and y direction) and ModifiedPyr_C BasisType (z "
             "direction)");

    int i, j, k, p, r, nFaceCoeffs;
    int nummodesA = 0, nummodesB = 0;

    int order0 = m_base[0]->GetNumModes();
    int order1 = m_base[1]->GetNumModes();
    int order2 = m_base[2]->GetNumModes();

    switch (fid)
    {
        case 0:
            nummodesA = order0;
            nummodesB = order1;
            break;
        case 1:
        case 3:
            nummodesA = order0;
            nummodesB = order2;
            break;
        case 2:
        case 4:
            nummodesA = order1;
            nummodesB = order2;
            break;
        default:
            ASSERTL0(false, "fid must be between 0 and 4");
    }

    if (P == -1)
    {
        P           = nummodesA;
        Q           = nummodesB;
        nFaceCoeffs = GetTraceNcoeffs(fid);
    }
    else if (fid > 0)
    {
        nFaceCoeffs = P * (2 * Q - P + 1) / 2;
    }
    else
    {
        nFaceCoeffs = P * Q;
    }

    // Allocate the map array and sign array; set sign array to ones (+)
    if (maparray.size() != nFaceCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nFaceCoeffs);
    }

    if (signarray.size() != nFaceCoeffs)
    {
        signarray = Array<OneD, int>(nFaceCoeffs, 1);
    }
    else
    {
        fill(signarray.data(), signarray.data() + nFaceCoeffs, 1);
    }

    // triangular faces
    if (fid > 0)
    {
        // zero signmap and set maparray to zero if elemental
        // modes are not as large as face modesl
        int idx   = 0;
        int cnt   = 0;
        int minPA = min(nummodesA, P);
        int minQB = min(nummodesB, Q);

        for (j = 0; j < minPA; ++j)
        {
            // set maparray
            for (k = 0; k < minQB - j; ++k, ++cnt)
            {
                maparray[idx++] = cnt;
            }

            cnt += nummodesB - minQB;

            for (k = nummodesB - j; k < Q - j; ++k)
            {
                signarray[idx]  = 0.0;
                maparray[idx++] = maparray[0];
            }
        }
        for (j = nummodesA; j < P; ++j)
        {
            for (k = 0; k < Q - j; ++k)
            {
                signarray[idx]  = 0.0;
                maparray[idx++] = maparray[0];
            }
        }

        // Triangles only have one possible orientation (base
        // direction reversed); swap edge modes.
        if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
        {
            swap(maparray[0], maparray[Q]);
            for (i = 1; i < Q - 1; ++i)
            {
                swap(maparray[i + 1], maparray[Q + i]);
            }

            idx = 0;
            for (p = 0; p < P; ++p)
            {
                for (r = 0; r < Q - p; ++r, idx++)
                {
                    if (p > 1)
                    {
                        signarray[idx] = p % 2 ? -1 : 1;
                    }
                }
            }
        }
    }
    else
    {

        // Set up an array indexing for quads, since the ordering may need
        // to be transposed.
        Array<OneD, int> arrayindx(nFaceCoeffs, -1);

        for (i = 0; i < Q; i++)
        {
            for (j = 0; j < P; j++)
            {
                if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
                {
                    arrayindx[i * P + j] = i * P + j;
                }
                else
                {
                    arrayindx[i * P + j] = j * Q + i;
                }
            }
        }

        // zero signmap and set maparray to zero if elemental
        // modes are not as large as face modesl
        for (j = 0; j < P; ++j)
        {
            // set up default maparray
            for (k = 0; k < Q; k++)
            {
                maparray[arrayindx[j + k * P]] = j + k * nummodesA;
            }

            for (k = nummodesB; k < Q; ++k)
            {
                signarray[arrayindx[j + k * P]] = 0.0;
                maparray[arrayindx[j + k * P]]  = maparray[0];
            }
        }

        for (j = nummodesA; j < P; ++j)
        {
            for (k = 0; k < Q; ++k)
            {
                signarray[arrayindx[j + k * P]] = 0.0;
                maparray[arrayindx[j + k * P]]  = maparray[0];
            }
        }

        // The code below is exactly the same as that taken from
        // StdHexExp and reverses the 'b' and 'a' directions as
        // appropriate (1st and 2nd if statements respectively) in
        // quadrilateral faces.
        if (faceOrient == eDir1FwdDir1_Dir2BwdDir2 ||
            faceOrient == eDir1BwdDir1_Dir2BwdDir2 ||
            faceOrient == eDir1BwdDir2_Dir2FwdDir1 ||
            faceOrient == eDir1BwdDir2_Dir2BwdDir1)
        {
            if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
            {
                for (i = 3; i < Q; i += 2)
                {
                    for (j = 0; j < P; j++)
                    {
                        signarray[arrayindx[i * P + j]] *= -1;
                    }
                }

                for (i = 0; i < P; i++)
                {
                    swap(maparray[i], maparray[i + P]);
                    swap(signarray[i], signarray[i + P]);
                }
            }
            else
            {
                for (i = 0; i < Q; i++)
                {
                    for (j = 3; j < P; j += 2)
                    {
                        signarray[arrayindx[i * P + j]] *= -1;
                    }
                }

                for (i = 0; i < Q; i++)
                {
                    swap(maparray[i], maparray[i + Q]);
                    swap(signarray[i], signarray[i + Q]);
                }
            }
        }

        if (faceOrient == eDir1BwdDir1_Dir2FwdDir2 ||
            faceOrient == eDir1BwdDir1_Dir2BwdDir2 ||
            faceOrient == eDir1FwdDir2_Dir2BwdDir1 ||
            faceOrient == eDir1BwdDir2_Dir2BwdDir1)
        {
            if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
            {
                for (i = 0; i < Q; i++)
                {
                    for (j = 3; j < P; j += 2)
                    {
                        signarray[arrayindx[i * P + j]] *= -1;
                    }
                }

                for (i = 0; i < Q; i++)
                {
                    swap(maparray[i * P], maparray[i * P + 1]);
                    swap(signarray[i * P], signarray[i * P + 1]);
                }
            }
            else
            {
                for (i = 3; i < Q; i += 2)
                {
                    for (j = 0; j < P; j++)
                    {
                        signarray[arrayindx[i * P + j]] *= -1;
                    }
                }

                for (i = 0; i < P; i++)
                {
                    swap(maparray[i * Q], maparray[i * Q + 1]);
                    swap(signarray[i * Q], signarray[i * Q + 1]);
                }
            }
        }
    }
}

void StdPyrExp::v_GetEdgeInteriorToElementMap(
    const int eid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation edgeOrient)
{
    int i;
    bool signChange;
    const int P              = m_base[0]->GetNumModes() - 1;
    const int Q              = m_base[1]->GetNumModes() - 1;
    const int R              = m_base[2]->GetNumModes() - 1;
    const int nEdgeIntCoeffs = v_GetEdgeNcoeffs(eid) - 2;

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

    // If edge is oriented backwards, change sign of modes which have
    // degree 2n+1, n >= 1.
    signChange = edgeOrient == eBackwards;

    switch (eid)
    {
        case 0:
            for (i = 2; i <= P; ++i)
            {
                maparray[i - 2] = GetMode(i, 0, 0);
            }
            break;

        case 1:
            for (i = 2; i <= Q; ++i)
            {
                maparray[i - 2] = GetMode(1, i, 0);
            }
            break;
        case 2:
            for (i = 2; i <= P; ++i)
            {
                maparray[i - 2] = GetMode(i, 1, 0);
            }
            break;

        case 3:
            for (i = 2; i <= Q; ++i)
            {
                maparray[i - 2] = GetMode(0, i, 0);
            }
            break;
        case 4:
            for (i = 2; i <= R; ++i)
            {
                maparray[i - 2] = GetMode(0, 0, i);
            }
            break;

        case 5:
            for (i = 1; i <= R - 1; ++i)
            {
                maparray[i - 1] = GetMode(1, 0, i);
            }
            break;
        case 6:
            for (i = 1; i <= R - 1; ++i)
            {
                maparray[i - 1] = GetMode(1, 1, i);
            }
            break;

        case 7:
            for (i = 1; i <= R - 1; ++i)
            {
                maparray[i - 1] = GetMode(0, 1, i);
            }
            break;
        default:
            ASSERTL0(false, "Edge not defined.");
            break;
    }

    if (signChange)
    {
        for (i = 1; i < nEdgeIntCoeffs; i += 2)
        {
            signarray[i] = -1;
        }
    }
}

void StdPyrExp::v_GetTraceInteriorToElementMap(
    const int fid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation faceOrient)
{
    const int P              = m_base[0]->GetNumModes() - 1;
    const int Q              = m_base[1]->GetNumModes() - 1;
    const int R              = m_base[2]->GetNumModes() - 1;
    const int nFaceIntCoeffs = v_GetTraceIntNcoeffs(fid);
    int p, q, r, idx = 0;
    int nummodesA = 0;
    int nummodesB = 0;
    int i, j;

    if (maparray.size() != nFaceIntCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nFaceIntCoeffs);
    }

    if (signarray.size() != nFaceIntCoeffs)
    {
        signarray = Array<OneD, int>(nFaceIntCoeffs, 1);
    }
    else
    {
        fill(signarray.data(), signarray.data() + nFaceIntCoeffs, 1);
    }

    // Set up an array indexing for quad faces, since the ordering may
    // need to be transposed depending on orientation.
    Array<OneD, int> arrayindx(nFaceIntCoeffs);
    if (fid == 0)
    {
        nummodesA = P - 1;
        nummodesB = Q - 1;

        for (i = 0; i < nummodesB; i++)
        {
            for (j = 0; j < nummodesA; j++)
            {
                if (faceOrient < 9)
                {
                    arrayindx[i * nummodesA + j] = i * nummodesA + j;
                }
                else
                {
                    arrayindx[i * nummodesA + j] = j * nummodesB + i;
                }
            }
        }
    }

    switch (fid)
    {
        case 0: // Bottom quad
            for (q = 2; q <= Q; ++q)
            {
                for (p = 2; p <= P; ++p)
                {
                    maparray[arrayindx[(q - 2) * nummodesA + (p - 2)]] =
                        GetMode(p, q, 0);
                }
            }
            break;
        case 1: // Front triangle
            for (p = 2; p <= P; ++p)
            {
                for (r = 1; r <= R - p; ++r)
                {
                    if ((int)faceOrient == 7)
                    {
                        signarray[idx] = p % 2 ? -1 : 1;
                    }
                    maparray[idx++] = GetMode(p, 0, r);
                }
            }
            break;
        case 2: // Right triangle
            for (q = 2; q <= Q; ++q)
            {
                for (r = 1; r <= R - q; ++r)
                {
                    if ((int)faceOrient == 7)
                    {
                        signarray[idx] = q % 2 ? -1 : 1;
                    }
                    maparray[idx++] = GetMode(1, q, r);
                }
            }
            break;

        case 3: // Rear triangle
            for (p = 2; p <= P; ++p)
            {
                for (r = 1; r <= R - p; ++r)
                {
                    if ((int)faceOrient == 7)
                    {
                        signarray[idx] = p % 2 ? -1 : 1;
                    }
                    maparray[idx++] = GetMode(p, 1, r);
                }
            }
            break;

        case 4: // Left triangle
            for (q = 2; q <= Q; ++q)
            {
                for (r = 1; r <= R - q; ++r)
                {
                    if ((int)faceOrient == 7)
                    {
                        signarray[idx] = q % 2 ? -1 : 1;
                    }
                    maparray[idx++] = GetMode(0, q, r);
                }
            }
            break;
        default:
            ASSERTL0(false, "Face interior map not available.");
    }

    // Triangular faces are processed in the above switch loop; for
    // remaining quad faces, set up orientation if necessary.
    if (fid > 0)
    {
        return;
    }

    if (faceOrient == 6 || faceOrient == 8 || faceOrient == 11 ||
        faceOrient == 12)
    {
        if (faceOrient < 9)
        {
            for (i = 1; i < nummodesB; i += 2)
            {
                for (j = 0; j < nummodesA; j++)
                {
                    signarray[arrayindx[i * nummodesA + j]] *= -1;
                }
            }
        }
        else
        {
            for (i = 0; i < nummodesB; i++)
            {
                for (j = 1; j < nummodesA; j += 2)
                {
                    signarray[arrayindx[i * nummodesA + j]] *= -1;
                }
            }
        }
    }

    if (faceOrient == 7 || faceOrient == 8 || faceOrient == 10 ||
        faceOrient == 12)
    {
        if (faceOrient < 9)
        {
            for (i = 0; i < nummodesB; i++)
            {
                for (j = 1; j < nummodesA; j += 2)
                {
                    signarray[arrayindx[i * nummodesA + j]] *= -1;
                }
            }
        }
        else
        {
            for (i = 1; i < nummodesB; i += 2)
            {
                for (j = 0; j < nummodesA; j++)
                {
                    signarray[arrayindx[i * nummodesA + j]] *= -1;
                }
            }
        }
    }
}

//---------------------------------------
// Wrapper functions
//---------------------------------------

DNekMatSharedPtr StdPyrExp::v_GenMatrix(const StdMatrixKey &mkey)
{
    return CreateGeneralMatrix(mkey);
}

DNekMatSharedPtr StdPyrExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

/**
 * @brief Compute the mode number in the expansion for a
 * particular tensorial combination.
 *
 * Modes are numbered with the r index travelling fastest,
 * followed by q and then p, and each q-r plane is of size
 *
 * (R+1-p)*(Q+1) - l(l+1)/2 where l = max(0,Q-p)
 *
 * For example, when P=2, Q=3 and R=4 the indexing inside each
 * q-r plane (with r increasing upwards and q to the right)
 * is:
 *
 * p = 0:      p = 1:       p = 2:
 * ----------------------------------
 * 4
 * 3 8         17 21
 * 2 7 11      16 20 24     29 32 35
 * 1 6 10 13   15 19 23 26  28 31 34 37
 * 0 5 9  12   14 18 22 25  27 30 33 36
 *
 * Note that in this element, we must have that \f$ P,Q \leq
 * R\f$.
 */
int StdPyrExp::GetMode(const int I, const int J, const int K)
{
    const int Q = m_base[1]->GetNumModes() - 1;
    const int R = m_base[2]->GetNumModes() - 1;

    int i, l;
    int cnt = 0;

    // Traverse to q-r plane number I
    for (i = 0; i < I; ++i)
    {
        // Size of triangle part
        l = max(0, Q - i);

        // Size of rectangle part
        cnt += (R + 1 - i) * (Q + 1) - l * (l + 1) / 2;
    }

    // Traverse to q column J (Pretend this is a face of width J)
    l = max(0, J - 1 - I);
    cnt += (R + 1 - I) * J - l * (l + 1) / 2;

    // Traverse up stacks to K
    cnt += K;

    return cnt;
}

void StdPyrExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                     const StdMatrixKey &mkey)
{
    // Generate an orthonogal expansion
    int qa       = m_base[0]->GetNumPoints();
    int qb       = m_base[1]->GetNumPoints();
    int qc       = m_base[2]->GetNumPoints();
    int nmodes_a = m_base[0]->GetNumModes();
    int nmodes_b = m_base[1]->GetNumModes();
    int nmodes_c = m_base[2]->GetNumModes();
    // Declare orthogonal basis.
    LibUtilities::PointsKey pa(qa, m_base[0]->GetPointsType());
    LibUtilities::PointsKey pb(qb, m_base[1]->GetPointsType());
    LibUtilities::PointsKey pc(qc, m_base[2]->GetPointsType());

    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, nmodes_a, pa);
    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, nmodes_b, pb);
    LibUtilities::BasisKey Bc(LibUtilities::eOrthoPyr_C, nmodes_c, pc);
    StdPyrExp OrthoExp(Ba, Bb, Bc);

    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
    int i, j, k, cnt = 0;

    // project onto modal  space.
    OrthoExp.FwdTrans(array, orthocoeffs);

    if (mkey.ConstFactorExists(eFactorSVVPowerKerDiffCoeff))
    {
        // Rodrigo's power kernel
        NekDouble cutoff = mkey.GetConstFactor(eFactorSVVCutoffRatio);
        NekDouble SvvDiffCoeff =
            mkey.GetConstFactor(eFactorSVVPowerKerDiffCoeff) *
            mkey.GetConstFactor(eFactorSVVDiffCoeff);

        for (i = 0; i < nmodes_a; ++i)
        {
            for (j = 0; j < nmodes_b; ++j)
            {
                int maxij      = max(i, j);
                NekDouble fac1 = std::max(
                    pow((1.0 * i) / (nmodes_a - 1), cutoff * nmodes_a),
                    pow((1.0 * j) / (nmodes_b - 1), cutoff * nmodes_b));

                for (k = 0; k < nmodes_c - maxij; ++k)
                {
                    NekDouble fac =
                        std::max(fac1, pow((1.0 * k) / (nmodes_c - 1),
                                           cutoff * nmodes_c));

                    orthocoeffs[cnt] *= SvvDiffCoeff * fac;
                    cnt++;
                }
            }
        }
    }
    else if (mkey.ConstFactorExists(
                 eFactorSVVDGKerDiffCoeff)) // Rodrigo/Mansoor's DG Kernel
    {
        NekDouble SvvDiffCoeff = mkey.GetConstFactor(eFactorSVVDGKerDiffCoeff) *
                                 mkey.GetConstFactor(eFactorSVVDiffCoeff);

        int max_abc = max(nmodes_a - kSVVDGFiltermodesmin,
                          nmodes_b - kSVVDGFiltermodesmin);
        max_abc     = max(max_abc, nmodes_c - kSVVDGFiltermodesmin);
        // clamp max_abc
        max_abc = max(max_abc, 0);
        max_abc = min(max_abc, kSVVDGFiltermodesmax - kSVVDGFiltermodesmin);

        for (i = 0; i < nmodes_a; ++i)
        {
            for (j = 0; j < nmodes_b; ++j)
            {
                int maxij = max(i, j);

                for (k = 0; k < nmodes_c - maxij; ++k)
                {
                    int maxijk = max(maxij, k);
                    maxijk     = min(maxijk, kSVVDGFiltermodesmax - 1);

                    orthocoeffs[cnt] *=
                        SvvDiffCoeff * kSVVDGFilter[max_abc][maxijk];
                    cnt++;
                }
            }
        }
    }
    else
    {
        // SVV filter paramaters (how much added diffusion relative
        // to physical one and fraction of modes from which you
        // start applying this added diffusion)
        //
        NekDouble SvvDiffCoeff =
            mkey.GetConstFactor(StdRegions::eFactorSVVDiffCoeff);
        NekDouble SVVCutOff =
            mkey.GetConstFactor(StdRegions::eFactorSVVCutoffRatio);

        // Defining the cut of mode
        int cutoff_a = (int)(SVVCutOff * nmodes_a);
        int cutoff_b = (int)(SVVCutOff * nmodes_b);
        int cutoff_c = (int)(SVVCutOff * nmodes_c);
        // To avoid the fac[j] from blowing up
        NekDouble epsilon = 1;

        int nmodes       = min(min(nmodes_a, nmodes_b), nmodes_c);
        NekDouble cutoff = min(min(cutoff_a, cutoff_b), cutoff_c);

        for (i = 0; i < nmodes_a; ++i) // P
        {
            for (j = 0; j < nmodes_b; ++j) // Q
            {
                int maxij = max(i, j);
                for (k = 0; k < nmodes_c - maxij; ++k) // R
                {
                    if (j + k >= cutoff || i + k >= cutoff)
                    {
                        orthocoeffs[cnt] *=
                            (SvvDiffCoeff *
                             exp(-(i + k - nmodes) * (i + k - nmodes) /
                                 ((NekDouble)((i + k - cutoff + epsilon) *
                                              (i + k - cutoff + epsilon)))) *
                             exp(-(j - nmodes) * (j - nmodes) /
                                 ((NekDouble)((j - cutoff + epsilon) *
                                              (j - cutoff + epsilon)))));
                    }
                    else
                    {
                        orthocoeffs[cnt] *= 0.0;
                    }
                    cnt++;
                }
            }
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdPyrExp::v_ReduceOrderCoeffs(int numMin,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nquad2  = m_base[2]->GetNumPoints();
    int nqtot   = nquad0 * nquad1 * nquad2;
    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();
    int nmodes2 = m_base[2]->GetNumModes();
    int numMax  = nmodes0;

    Array<OneD, NekDouble> coeff(m_ncoeffs);
    Array<OneD, NekDouble> coeff_tmp1(m_ncoeffs, 0.0);
    Array<OneD, NekDouble> phys_tmp(nqtot, 0.0);
    Array<OneD, NekDouble> tmp, tmp2, tmp3, tmp4;

    const LibUtilities::PointsKey Pkey0 = m_base[0]->GetPointsKey();
    const LibUtilities::PointsKey Pkey1 = m_base[1]->GetPointsKey();
    const LibUtilities::PointsKey Pkey2 = m_base[2]->GetPointsKey();

    LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A, nmodes0, Pkey0);
    LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_A, nmodes1, Pkey1);
    LibUtilities::BasisKey bortho2(LibUtilities::eOrthoPyr_C, nmodes2, Pkey2);

    int cnt = 0;
    int u   = 0;
    int i   = 0;
    StdRegions::StdPyrExpSharedPtr OrthoPyrExp;

    OrthoPyrExp = MemoryManager<StdRegions::StdPyrExp>::AllocateSharedPtr(
        bortho0, bortho1, bortho2);

    BwdTrans(inarray, phys_tmp);
    OrthoPyrExp->FwdTrans(phys_tmp, coeff);

    // filtering
    for (u = 0; u < numMin; ++u)
    {
        for (i = 0; i < numMin; ++i)
        {

            int maxui = max(u, i);
            Vmath::Vcopy(numMin - maxui, tmp = coeff + cnt, 1,
                         tmp2 = coeff_tmp1 + cnt, 1);
            cnt += nmodes2 - maxui;
        }

        for (i = numMin; i < nmodes1; ++i)
        {
            int maxui = max(u, i);
            cnt += numMax - maxui;
        }
    }

    OrthoPyrExp->BwdTrans(coeff_tmp1, phys_tmp);
    StdPyrExp::FwdTrans(phys_tmp, outarray);
}

} // namespace Nektar::StdRegions
