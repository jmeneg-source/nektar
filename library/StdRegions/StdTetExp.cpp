///////////////////////////////////////////////////////////////////////////////
//
// File: StdTetExp.cpp
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
// Description: Header field for tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <StdRegions/StdTetExp.h>

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

StdTetExp::StdTetExp(const LibUtilities::BasisKey &Ba,
                     const LibUtilities::BasisKey &Bb,
                     const LibUtilities::BasisKey &Bc)
    : StdExpansion(LibUtilities::StdTetData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                   3, Ba, Bb, Bc),
      StdExpansion3D(LibUtilities::StdTetData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                     Ba, Bb, Bc)
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

    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());

    StdFacKey w1key(eWeights1, Bb);
    // get weights[1] from manager where points are rescaled
    m_weights.push_back(GetStdFac(w1key));

    StdFacKey w2key(eWeights2, Bc);
    // get weights[2] from manager where points are rescaled
    m_weights.push_back(GetStdFac(w2key));
}

//----------------------------
// Differentiation Methods
//----------------------------
/**
 * \brief Calculate the derivative of the physical points
 *
 * The derivative is evaluated at the nodal physical points.
 * Derivatives with respect to the local Cartesian coordinates
 *
 * \f$\begin{Bmatrix} \frac {\partial} {\partial \xi_1} \\ \frac
 * {\partial} {\partial \xi_2} \\ \frac {\partial} {\partial \xi_3}
 * \end{Bmatrix} = \begin{Bmatrix} \frac 4 {(1-\eta_2)(1-\eta_3)}
 * \frac \partial {\partial \eta_1} \ \ \frac {2(1+\eta_1)}
 * {(1-\eta_2)(1-\eta_3)} \frac \partial {\partial \eta_1} + \frac 2
 * {1-\eta_3} \frac \partial {\partial \eta_3} \\ \frac {2(1 +
 * \eta_1)} {2(1 - \eta_2)(1-\eta_3)} \frac \partial {\partial \eta_1}
 * + \frac {1 + \eta_2} {1 - \eta_3} \frac \partial {\partial \eta_2}
 * + \frac \partial {\partial \eta_3} \end{Bmatrix}\f$
 **/
void StdTetExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &out_dxi0,
                               Array<OneD, NekDouble> &out_dxi1,
                               Array<OneD, NekDouble> &out_dxi2)
{
    int Q0   = m_base[0]->GetNumPoints();
    int Q1   = m_base[1]->GetNumPoints();
    int Q2   = m_base[2]->GetNumPoints();
    int Qtot = Q0 * Q1 * Q2;

    // Compute the physical derivative
    Array<OneD, NekDouble> out_dEta0(3 * Qtot, 0.0);
    Array<OneD, NekDouble> out_dEta1 = out_dEta0 + Qtot;
    Array<OneD, NekDouble> out_dEta2 = out_dEta1 + Qtot;

    bool Do_2 = (out_dxi2.size() > 0) ? true : false;
    bool Do_1 = (out_dxi1.size() > 0) ? true : false;

    if (Do_2) // Need all local derivatives
    {
        PhysTensorDeriv(inarray, out_dEta0, out_dEta1, out_dEta2);
    }
    else if (Do_1) // Need 0 and 1 derivatives
    {
        PhysTensorDeriv(inarray, out_dEta0, out_dEta1, NullNekDouble1DArray);
    }
    else // Only need Eta0 derivaitve
    {
        PhysTensorDeriv(inarray, out_dEta0, NullNekDouble1DArray,
                        NullNekDouble1DArray);
    }

    Array<OneD, const NekDouble> eta_0, eta_1, eta_2;
    eta_0 = m_base[0]->GetZ();
    eta_1 = m_base[1]->GetZ();
    eta_2 = m_base[2]->GetZ();

    // calculate 2.0/((1-eta_1)(1-eta_2)) Out_dEta0

    NekDouble *dEta0 = &out_dEta0[0];
    NekDouble fac;
    for (int k = 0; k < Q2; ++k)
    {
        for (int j = 0; j < Q1; ++j, dEta0 += Q0)
        {
            Vmath::Smul(Q0, 2.0 / (1.0 - eta_1[j]), dEta0, 1, dEta0, 1);
        }
        fac = 1.0 / (1.0 - eta_2[k]);
        Vmath::Smul(Q0 * Q1, fac, &out_dEta0[0] + k * Q0 * Q1, 1,
                    &out_dEta0[0] + k * Q0 * Q1, 1);
    }

    if (out_dxi0.size() > 0)
    {
        // out_dxi0 = 4.0/((1-eta_1)(1-eta_2)) Out_dEta0
        Vmath::Smul(Qtot, 2.0, out_dEta0, 1, out_dxi0, 1);
    }

    if (Do_1 || Do_2)
    {
        Array<OneD, NekDouble> Fac0(Q0);
        Vmath::Sadd(Q0, 1.0, eta_0, 1, Fac0, 1);

        // calculate 2.0*(1+eta_0)/((1-eta_1)(1-eta_2)) Out_dEta0
        for (int k = 0; k < Q1 * Q2; ++k)
        {
            Vmath::Vmul(Q0, &Fac0[0], 1, &out_dEta0[0] + k * Q0, 1,
                        &out_dEta0[0] + k * Q0, 1);
        }
        // calculate 2/(1.0-eta_2) out_dEta1
        for (int k = 0; k < Q2; ++k)
        {
            Vmath::Smul(Q0 * Q1, 2.0 / (1.0 - eta_2[k]),
                        &out_dEta1[0] + k * Q0 * Q1, 1,
                        &out_dEta1[0] + k * Q0 * Q1, 1);
        }

        if (Do_1)
        {
            // calculate out_dxi1 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) Out_dEta0
            // + 2/(1.0-eta_2) out_dEta1
            Vmath::Vadd(Qtot, out_dEta0, 1, out_dEta1, 1, out_dxi1, 1);
        }

        if (Do_2)
        {
            // calculate (1 + eta_1)/(1 -eta_2)*out_dEta1
            NekDouble *dEta1 = &out_dEta1[0];
            for (int k = 0; k < Q2; ++k)
            {
                for (int j = 0; j < Q1; ++j, dEta1 += Q0)
                {
                    Vmath::Smul(Q0, (1.0 + eta_1[j]) / 2.0, dEta1, 1, dEta1, 1);
                }
            }

            // calculate out_dxi2 =
            // 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) Out_dEta0 +
            // (1 + eta_1)/(1 -eta_2)*out_dEta1 + out_dEta2
            Vmath::Vadd(Qtot, out_dEta0, 1, out_dEta1, 1, out_dxi2, 1);
            Vmath::Vadd(Qtot, out_dEta2, 1, out_dxi2, 1, out_dxi2, 1);
        }
    }
}

//---------------------------------------
// Transforms
//---------------------------------------

/**
 * @note 'r' (base[2]) runs fastest in this element
 *
 * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
 *  u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
 *
 * Backward transformation is three dimensional tensorial expansion
 * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
 * (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{pq}^b (\xi_{2j})
 * \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k})
 * \rbrace} \rbrace}. \f$ And sumfactorizing step of the form is as:\\
 *
 * \f$ f_{pq} (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c
 * (\xi_{3k}),\\
 *
 * g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y} \psi_{pq}^b
 * (\xi_{2j}) f_{pq} (\xi_{3k}),\\
 *
 * u(\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a
 * (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).  \f$
 */
void StdTetExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &outarray)
{
    ASSERTL1((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B) ||
                 (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
             "Basis[1] is not a general tensor type");

    ASSERTL1((m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                 (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
             "Basis[2] is not a general tensor type");

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
    BwdTransTetKernel(nmodes0, nmodes1, nmodes2, nquad0, nquad1, nquad2,       \
                      isModified, (const vec_t *)base0.data(),                 \
                      (const vec_t *)base1.data(),                             \
                      (const vec_t *)base2.data(), wsp0.data(), wsp1.data(),   \
                      (const vec_t *)inarray.data(), (vec_t *)outarray.data())

    // Inner loop case over quarature points
#undef BWDTRANS_Q
#define BWDTRANS_Q(r, i)                                                       \
    case NQ(i):                                                                \
        BwdTransTetKernel(                                                     \
            NM(i), NM(i), NM(i), NQ(i), NQ_M1(i), NQ_M1(i), isModified,        \
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
    if ((nmodes0 == nmodes1) && (nmodes1 == nmodes2) &&
        (nquad0 == nquad1 + 1) && (nquad1 == nquad2))
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
void StdTetExp::v_IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &base1,
    const Array<OneD, const NekDouble> &base2,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, NekDouble> &jac,
    const bool Deformed, [[maybe_unused]] bool CollDir0,
    [[maybe_unused]] bool CollDir1, [[maybe_unused]] bool CollDir2)
{
    ASSERTL1((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B) ||
                 (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
             "Basis[1] is not a general tensor type");

    ASSERTL1((m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                 (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
             "Basis[2] is not a general tensor type");

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
    IProductTetKernel<false, false, true>(                                     \
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
        IProductTetKernel<false, false, true>(                                 \
            NM(i), NM(i), NM(i), NQ(i), NQ_M1(i), NQ_M1(i), isModified,        \
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
        if ((order0 == order1) && (order1 == order2) &&
            (nquad0 == nquad1 + 1) && (nquad1 == nquad2))
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
    IProductTetKernel<false, false, false>(                                    \
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
        IProductTetKernel<false, false, false>(                                \
            NM(i), NM(i), NM(i), NQ(i), NQ_M1(i), NQ_M1(i), isModified,        \
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
        if ((order0 == order1) && (order1 == order2) &&
            (nquad0 == nquad1 + 1) && (nquad1 == nquad2))
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

/**
 * @param   inarray     Function evaluated at physical collocation
 *                      points.
 * @param   outarray    Inner product with respect to each basis
 *                      function over the element.
 */
void StdTetExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int i;
    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nquad2  = m_base[2]->GetNumPoints();
    int nqtot   = nquad0 * nquad1 * nquad2;
    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();

    Array<OneD, NekDouble> tmp0(max(nqtot, m_ncoeffs));
    Array<OneD, NekDouble> wsp(nquad1 * nquad2 * nmodes0 +
                               nquad2 * nmodes0 * (2 * nmodes1 - nmodes0 + 1) /
                                   2);

    StdFacKey fackey0(eHalfMultOnePlusZ0, m_base[0]->GetBasisKey());
    Array<OneD, const NekDouble> gfac0 = GetStdFac(fackey0);
    StdFacKey fackey1(eTwoOverOneMinusZ1, m_base[1]->GetBasisKey());
    Array<OneD, const NekDouble> gfac1 = GetStdFac(fackey1);
    StdFacKey fackey2(eTwoOverOneMinusZ2, m_base[2]->GetBasisKey());
    Array<OneD, const NekDouble> gfac2 = GetStdFac(fackey2);

    // Derivative in first direction is always scaled as follows
    for (i = 0; i < nquad1 * nquad2; ++i)
    {
        Vmath::Smul(nquad0, gfac1[i % nquad1], &inarray[0] + i * nquad0, 1,
                    &tmp0[0] + i * nquad0, 1);
    }
    for (i = 0; i < nquad2; ++i)
    {
        Vmath::Smul(nquad0 * nquad1, gfac2[i], &tmp0[0] + i * nquad0 * nquad1,
                    1, &tmp0[0] + i * nquad0 * nquad1, 1);
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
            Array<OneD, NekDouble> tmp3(m_ncoeffs);

            for (i = 0; i < nquad1 * nquad2; ++i)
            {
                Vmath::Vmul(nquad0, &gfac0[0], 1, &tmp0[0] + i * nquad0, 1,
                            &tmp0[0] + i * nquad0, 1);
            }

            v_IProductWRTBaseKernel(
                m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                m_base[2]->GetBdata(), tmp0, tmp3, one, false);

            for (i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nquad0 * nquad1, gfac2[i],
                            &inarray[0] + i * nquad0 * nquad1, 1,
                            &tmp0[0] + i * nquad0 * nquad1, 1);
            }

            v_IProductWRTBaseKernel(
                m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                m_base[2]->GetBdata(), tmp0, outarray, one, false);
            Vmath::Vadd(m_ncoeffs, &tmp3[0], 1, &outarray[0], 1, &outarray[0],
                        1);
        }
        break;
        case 2:
        {
            Array<OneD, NekDouble> tmp3(m_ncoeffs);
            Array<OneD, NekDouble> tmp4(m_ncoeffs);
            StdFacKey fackey1a(eHalfMultOnePlusZ1, m_base[1]->GetBasisKey());
            gfac1 = GetStdFac(fackey1a);

            for (i = 0; i < nquad1 * nquad2; ++i)
            {
                Vmath::Vmul(nquad0, &gfac0[0], 1, &tmp0[0] + i * nquad0, 1,
                            &tmp0[0] + i * nquad0, 1);
            }
            v_IProductWRTBaseKernel(
                m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                m_base[2]->GetBdata(), tmp0, tmp3, one, false);

            for (i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nquad0 * nquad1, gfac2[i],
                            &inarray[0] + i * nquad0 * nquad1, 1,
                            &tmp0[0] + i * nquad0 * nquad1, 1);
            }
            for (i = 0; i < nquad1 * nquad2; ++i)
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
        }
        break;
        default:
        {
            ASSERTL1(false, "input dir is out of range");
        }
        break;
    }
}

//---------------------------------------
// Evaluation functions
//---------------------------------------

void StdTetExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble> &xi,
                                         Array<OneD, NekDouble> &eta)
{
    NekDouble d2  = 1.0 - xi[2];
    NekDouble d12 = -xi[1] - xi[2];
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
    if (fabs(d12) < NekConstants::kNekZeroTol)
    {
        if (d12 >= 0.)
        {
            d12 = NekConstants::kNekZeroTol;
        }
        else
        {
            d12 = -NekConstants::kNekZeroTol;
        }
    }
    eta[0] = 2.0 * (1.0 + xi[0]) / d12 - 1.0;
    eta[1] = 2.0 * (1.0 + xi[1]) / d2 - 1.0;
    eta[2] = xi[2];
}

void StdTetExp::v_LocCollapsedToLocCoord(
    const Array<OneD, const NekDouble> &eta, Array<OneD, NekDouble> &xi)
{
    xi[2] = eta[2];
    xi[1] = (1.0 + eta[1]) * (1.0 - xi[2]) * 0.5 - 1.0;
    xi[0] = (1.0 + eta[0]) * (-xi[1] - xi[2]) * 0.5 - 1.0;
}

void StdTetExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp(m_ncoeffs, 0.0);
    tmp[mode] = 1.0;
    StdTetExp::v_BwdTrans(tmp, outarray);
}

NekDouble StdTetExp::v_PhysEvaluateBasis(
    const Array<OneD, const NekDouble> &coords, int mode)
{
    Array<OneD, NekDouble> coll(3);
    LocCoordToLocCollapsed(coords, coll);

    const int nm1 = m_base[1]->GetNumModes();
    const int nm2 = m_base[2]->GetNumModes();

    const int b     = 2 * nm2 + 1;
    const int mode0 = floor(0.5 * (b - sqrt(b * b - 8.0 * mode / nm1)));
    const int tmp =
        mode - nm1 * (mode0 * (nm2 - 1) + 1 - (mode0 - 2) * (mode0 - 1) / 2);
    const int mode1 = tmp / (nm2 - mode0);
    const int mode2 = tmp % (nm2 - mode0);

    if (m_base[0]->GetBasisType() == LibUtilities::eModified_A)
    {
        // Handle the collapsed vertices and edges in the modified
        // basis.
        if (mode == 1)
        {
            // Collapsed top vertex
            return StdExpansion::BaryEvaluateBasis<2>(coll[2], 1);
        }
        else if (mode0 == 0 && mode2 == 1)
        {
            return StdExpansion::BaryEvaluateBasis<1>(coll[1], 0) *
                   StdExpansion::BaryEvaluateBasis<2>(coll[2], 1);
        }
        else if (mode0 == 1 && mode1 == 1 && mode2 == 0)
        {
            return StdExpansion::BaryEvaluateBasis<0>(coll[0], 0) *
                   StdExpansion::BaryEvaluateBasis<1>(coll[1], 1);
        }
    }

    return StdExpansion::BaryEvaluateBasis<0>(coll[0], mode0) *
           StdExpansion::BaryEvaluateBasis<1>(coll[1], mode1) *
           StdExpansion::BaryEvaluateBasis<2>(coll[2], mode2);
}

NekDouble StdTetExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    // Collapse coordinates
    Array<OneD, NekDouble> coll(3, 0.0);
    LocCoordToLocCollapsed(coord, coll);

    // If near singularity do the old interpolation matrix method
    if ((1 - coll[1]) < 1e-5 || (1 - coll[2]) < 1e-5)
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
    NekDouble val = BaryTensorDeriv(coll, inarray, interDeriv);

    // calculate 2.0/((1-eta_1)(1-eta_2)) * Out_dEta0
    NekDouble temp = 2.0 / ((1 - coll[1]) * (1 - coll[2]));
    interDeriv[0] *= temp;

    // out_dxi0 = 4.0/((1-eta_1)(1-eta_2)) * Out_dEta0
    firstOrderDerivs[0] = 2 * interDeriv[0];

    // fac0 = 1 + eta_0
    NekDouble fac0;
    fac0 = 1 + coll[0];

    // calculate 2.0*(1+eta_0)/((1-eta_1)(1-eta_2)) * Out_dEta0
    interDeriv[0] *= fac0;

    // calculate 2/(1.0-eta_2) * out_dEta1
    fac0 = 2 / (1 - coll[2]);
    interDeriv[1] *= fac0;

    // calculate out_dxi1 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2))
    //  * Out_dEta0 + 2/(1.0-eta_2) out_dEta1
    firstOrderDerivs[1] = interDeriv[0] + interDeriv[1];

    // calculate (1 + eta_1)/(1 -eta_2)*out_dEta1
    fac0 = (1 + coll[1]) / 2;
    interDeriv[1] *= fac0;

    // calculate out_dxi2 =
    // 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) Out_dEta0 +
    // (1 + eta_1)/(1 -eta_2)*out_dEta1 + out_dEta2
    firstOrderDerivs[2] = interDeriv[0] + interDeriv[1] + interDeriv[2];

    return val;
}

void StdTetExp::v_GetTraceNumModes(const int fid, int &numModes0,
                                   int &numModes1,
                                   [[maybe_unused]] Orientation faceOrient)
{
    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};
    switch (fid)
    {
        case 0:
        {
            numModes0 = nummodes[0];
            numModes1 = nummodes[1];
        }
        break;
        case 1:
        {
            numModes0 = nummodes[0];
            numModes1 = nummodes[2];
        }
        break;
        case 2:
        case 3:
        {
            numModes0 = nummodes[1];
            numModes1 = nummodes[2];
        }
        break;
    }
}

//---------------------------
// Helper functions
//---------------------------

int StdTetExp::v_GetNverts() const
{
    return 4;
}

int StdTetExp::v_GetNedges() const
{
    return 6;
}

int StdTetExp::v_GetNtraces() const
{
    return 4;
}

LibUtilities::ShapeType StdTetExp::v_DetShapeType() const
{
    return LibUtilities::eTetrahedron;
}

int StdTetExp::v_NumBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes();
    int Q = m_base[1]->GetNumModes();
    int R = m_base[2]->GetNumModes();

    return LibUtilities::StdTetData::getNumberOfBndCoefficients(P, Q, R);
}

int StdTetExp::v_NumDGBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes() - 1;
    int Q = m_base[1]->GetNumModes() - 1;
    int R = m_base[2]->GetNumModes() - 1;

    return (Q + 1) + P * (1 + 2 * Q - P) / 2    // base face
           + (R + 1) + P * (1 + 2 * R - P) / 2  // front face
           + 2 * (R + 1) + Q * (1 + 2 * R - Q); // back two faces
}

int StdTetExp::v_GetTraceNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 3), "face id is out of range");
    int nFaceCoeffs = 0;
    int nummodesA, nummodesB, P, Q;
    if (i == 0)
    {
        nummodesA = GetBasisNumModes(0);
        nummodesB = GetBasisNumModes(1);
    }
    else if ((i == 1) || (i == 2))
    {
        nummodesA = GetBasisNumModes(0);
        nummodesB = GetBasisNumModes(2);
    }
    else
    {
        nummodesA = GetBasisNumModes(1);
        nummodesB = GetBasisNumModes(2);
    }
    P           = nummodesA - 1;
    Q           = nummodesB - 1;
    nFaceCoeffs = Q + 1 + (P * (1 + 2 * Q - P)) / 2;
    return nFaceCoeffs;
}

int StdTetExp::v_GetTraceIntNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 3), "face id is out of range");
    int Pi = m_base[0]->GetNumModes() - 2;
    int Qi = m_base[1]->GetNumModes() - 2;
    int Ri = m_base[2]->GetNumModes() - 2;

    if ((i == 0))
    {
        return Pi * (2 * Qi - Pi - 1) / 2;
    }
    else if ((i == 1))
    {
        return Pi * (2 * Ri - Pi - 1) / 2;
    }
    else
    {
        return Qi * (2 * Ri - Qi - 1) / 2;
    }
}

int StdTetExp::v_GetTraceNumPoints(const int i) const
{
    ASSERTL2(i >= 0 && i <= 3, "face id is out of range");

    if (i == 0)
    {
        return m_base[0]->GetNumPoints() * m_base[1]->GetNumPoints();
    }
    else if (i == 1)
    {
        return m_base[0]->GetNumPoints() * m_base[2]->GetNumPoints();
    }
    else
    {
        return m_base[1]->GetNumPoints() * m_base[2]->GetNumPoints();
    }
}

int StdTetExp::v_GetEdgeNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 5), "edge id is out of range");
    int P = m_base[0]->GetNumModes();
    int Q = m_base[1]->GetNumModes();
    int R = m_base[2]->GetNumModes();

    if (i == 0)
    {
        return P;
    }
    else if (i == 1 || i == 2)
    {
        return Q;
    }
    else
    {
        return R;
    }
}

LibUtilities::PointsKey StdTetExp::v_GetTracePointsKey(const int i,
                                                       const int j) const
{
    ASSERTL2(i >= 0 && i <= 3, "face id is out of range");
    ASSERTL2(j == 0 || j == 1, "face direction is out of range");

    if (i == 0)
    {
        return m_base[j]->GetPointsKey();
    }
    else if (i == 1)
    {
        return m_base[2 * j]->GetPointsKey();
    }
    else
    {
        return m_base[j + 1]->GetPointsKey();
    }
}

int StdTetExp::v_CalcNumberOfCoefficients(
    const std::vector<unsigned int> &nummodes, int &modes_offset)
{
    int nmodes = LibUtilities::StdTetData::getNumberOfCoefficients(
        nummodes[modes_offset], nummodes[modes_offset + 1],
        nummodes[modes_offset + 2]);
    modes_offset += 3;

    return nmodes;
}

const LibUtilities::BasisKey StdTetExp::v_GetTraceBasisKey(const int i,
                                                           const int k,
                                                           bool UseGLL) const
{
    ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
    ASSERTL2(k == 0 || k == 1, "face direction out of range");

    int dir = k;
    switch (i)
    {
        case 0:
            dir = k; // retrun facedir=0-> 0 facedir=1->1
            break;
        case 1:
            dir = 2 * k; // retrun facedir=0-> 0 facedir=1->2
            break;
        case 2:
        case 3:
            dir = k + 1; // retrun facedir=0-> 1 facedir=1->2
            break;
    }

    return EvaluateTriFaceBasisKey(k, m_base[dir], UseGLL);
}

void StdTetExp::v_GetCoords(Array<OneD, NekDouble> &xi_x,
                            Array<OneD, NekDouble> &xi_y,
                            Array<OneD, NekDouble> &xi_z)
{
    Array<OneD, const NekDouble> eta_x = m_base[0]->GetZ();
    Array<OneD, const NekDouble> eta_y = m_base[1]->GetZ();
    Array<OneD, const NekDouble> eta_z = m_base[2]->GetZ();
    int Qx                             = GetNumPoints(0);
    int Qy                             = GetNumPoints(1);
    int Qz                             = GetNumPoints(2);

    // Convert collapsed coordinates into cartesian coordinates: eta
    // --> xi
    for (int k = 0; k < Qz; ++k)
    {
        for (int j = 0; j < Qy; ++j)
        {
            for (int i = 0; i < Qx; ++i)
            {
                int s = i + Qx * (j + Qy * k);
                xi_x[s] =
                    (eta_x[i] + 1.0) * (1.0 - eta_y[j]) * (1.0 - eta_z[k]) / 4 -
                    1.0;
                xi_y[s] = (eta_y[j] + 1.0) * (1.0 - eta_z[k]) / 2 - 1.0;
                xi_z[s] = eta_z[k];
            }
        }
    }
}

bool StdTetExp::v_IsBoundaryInteriorExpansion() const
{
    return (m_base[0]->GetBasisType() == LibUtilities::eModified_A) &&
           (m_base[1]->GetBasisType() == LibUtilities::eModified_B) &&
           (m_base[2]->GetBasisType() == LibUtilities::eModified_C);
}

//--------------------------
// Mappings
//--------------------------
int StdTetExp::v_GetVertexMap(const int localVertexId, bool useCoeffPacking)
{
    ASSERTL0((GetBasisType(0) == LibUtilities::eModified_A) ||
                 (GetBasisType(1) == LibUtilities::eModified_B) ||
                 (GetBasisType(2) == LibUtilities::eModified_C),
             "Mapping not defined for this type of basis");

    int localDOF = 0;
    if (useCoeffPacking == true) // follow packing of coefficients i.e q,r,p
    {
        switch (localVertexId)
        {
            case 0:
            {
                localDOF = GetMode(0, 0, 0);
                break;
            }
            case 1:
            {
                localDOF = GetMode(0, 0, 1);
                break;
            }
            case 2:
            {
                localDOF = GetMode(0, 1, 0);
                break;
            }
            case 3:
            {
                localDOF = GetMode(1, 0, 0);
                break;
            }
            default:
            {
                ASSERTL0(false, "Vertex ID must be between 0 and 3");
                break;
            }
        }
    }
    else
    {
        switch (localVertexId)
        {
            case 0:
            {
                localDOF = GetMode(0, 0, 0);
                break;
            }
            case 1:
            {
                localDOF = GetMode(1, 0, 0);
                break;
            }
            case 2:
            {
                localDOF = GetMode(0, 1, 0);
                break;
            }
            case 3:
            {
                localDOF = GetMode(0, 0, 1);
                break;
            }
            default:
            {
                ASSERTL0(false, "Vertex ID must be between 0 and 3");
                break;
            }
        }
    }

    return localDOF;
}

/**
 * Maps interior modes of an edge to the elemental modes.
 */

/**
 * List of all interior modes in the expansion.
 */
void StdTetExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes();
    int Q = m_base[1]->GetNumModes();
    int R = m_base[2]->GetNumModes();

    int nIntCoeffs = m_ncoeffs - NumBndryCoeffs();

    if (outarray.size() != nIntCoeffs)
    {
        outarray = Array<OneD, unsigned int>(nIntCoeffs);
    }

    int idx = 0;
    for (int i = 2; i < P; ++i)
    {
        for (int j = 1; j < Q - i; ++j)
        {
            for (int k = 1; k < R - i - j; ++k)
            {
                outarray[idx++] = GetMode(i, j, k);
            }
        }
    }
}

/**
 * List of all boundary modes in the the expansion.
 */
void StdTetExp::v_GetBoundaryMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes();
    int Q = m_base[1]->GetNumModes();
    int R = m_base[2]->GetNumModes();

    int i, j, k;
    int idx = 0;

    int nBnd = NumBndryCoeffs();

    if (outarray.size() != nBnd)
    {
        outarray = Array<OneD, unsigned int>(nBnd);
    }

    for (i = 0; i < P; ++i)
    {
        // First two Q-R planes are entirely boundary modes
        if (i < 2)
        {
            for (j = 0; j < Q - i; j++)
            {
                for (k = 0; k < R - i - j; ++k)
                {
                    outarray[idx++] = GetMode(i, j, k);
                }
            }
        }
        // Remaining Q-R planes contain boundary modes on bottom and
        // left edge.
        else
        {
            for (k = 0; k < R - i; ++k)
            {
                outarray[idx++] = GetMode(i, 0, k);
            }
            for (j = 1; j < Q - i; ++j)
            {
                outarray[idx++] = GetMode(i, j, 0);
            }
        }
    }
}

void StdTetExp::v_GetTraceCoeffMap(const unsigned int fid,
                                   Array<OneD, unsigned int> &maparray)
{
    int i, j, k;
    int P = 0, Q = 0, idx = 0;
    int nFaceCoeffs = 0;

    switch (fid)
    {
        case 0:
            P = m_base[0]->GetNumModes();
            Q = m_base[1]->GetNumModes();
            break;
        case 1:
            P = m_base[0]->GetNumModes();
            Q = m_base[2]->GetNumModes();
            break;
        case 2:
        case 3:
            P = m_base[1]->GetNumModes();
            Q = m_base[2]->GetNumModes();
            break;
        default:
            ASSERTL0(false, "fid must be between 0 and 3");
    }

    nFaceCoeffs = P * (2 * Q - P + 1) / 2;

    if (maparray.size() != nFaceCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nFaceCoeffs);
    }

    switch (fid)
    {
        case 0:
            idx = 0;
            for (i = 0; i < P; ++i)
            {
                for (j = 0; j < Q - i; ++j)
                {
                    maparray[idx++] = GetMode(i, j, 0);
                }
            }
            break;
        case 1:
            idx = 0;
            for (i = 0; i < P; ++i)
            {
                for (k = 0; k < Q - i; ++k)
                {
                    maparray[idx++] = GetMode(i, 0, k);
                }
            }
            break;
        case 2:
            idx = 0;
            for (j = 0; j < P - 1; ++j)
            {
                for (k = 0; k < Q - 1 - j; ++k)
                {
                    maparray[idx++] = GetMode(1, j, k);
                    // Incorporate modes from zeroth plane where needed.
                    // Add in top vertex
                    if (j == 0 && k == 0)
                    {
                        maparray[idx++] = GetMode(0, 0, 1);
                    }
                    // Add in bottom  singular vertex  plus singular edge
                    if (j == 0 && k == Q - 2)
                    {
                        for (int r = 0; r < Q - 1; ++r)
                        {
                            maparray[idx++] = GetMode(0, 1, r);
                        }
                    }
                }
            }
            break;
        case 3:
            idx = 0;
            for (j = 0; j < P; ++j)
            {
                for (k = 0; k < Q - j; ++k)
                {
                    maparray[idx++] = GetMode(0, j, k);
                }
            }
            break;
        default:
            ASSERTL0(false, "Element map not available.");
    }
}

void StdTetExp::v_GetElmtTraceToTraceMap(const unsigned int fid,
                                         Array<OneD, unsigned int> &maparray,
                                         Array<OneD, int> &signarray,
                                         Orientation faceOrient, int P, int Q)
{
    int nummodesA = 0, nummodesB = 0, i, j, k, idx;

    ASSERTL1(v_IsBoundaryInteriorExpansion(),
             "Method only implemented for Modified_A BasisType (x "
             "direction), Modified_B BasisType (y direction), and "
             "Modified_C BasisType(z direction)");

    int nFaceCoeffs = 0;

    switch (fid)
    {
        case 0:
            nummodesA = m_base[0]->GetNumModes();
            nummodesB = m_base[1]->GetNumModes();
            break;
        case 1:
            nummodesA = m_base[0]->GetNumModes();
            nummodesB = m_base[2]->GetNumModes();
            break;
        case 2:
        case 3:
            nummodesA = m_base[1]->GetNumModes();
            nummodesB = m_base[2]->GetNumModes();
            break;
        default:
            ASSERTL0(false, "fid must be between 0 and 3");
    }

    if (P == -1)
    {
        P = nummodesA;
        Q = nummodesB;
    }

    nFaceCoeffs = P * (2 * Q - P + 1) / 2;

    // Allocate the map array and sign array; set sign array to ones (+)
    if (maparray.size() != nFaceCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nFaceCoeffs, 1);
    }

    if (signarray.size() != nFaceCoeffs)
    {
        signarray = Array<OneD, int>(nFaceCoeffs, 1);
    }
    else
    {
        fill(signarray.data(), signarray.data() + nFaceCoeffs, 1);
    }

    // zero signmap and set maparray to zero if elemental
    // modes are not as large as face modesl
    idx       = 0;
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

    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
    {
        idx = 0;
        for (i = 0; i < P; ++i)
        {
            for (j = 0; j < Q - i; ++j, idx++)
            {
                if (i > 1)
                {
                    signarray[idx] = (i % 2 ? -1 : 1);
                }
            }
        }

        swap(maparray[0], maparray[Q]);

        for (i = 1; i < Q - 1; ++i)
        {
            swap(maparray[i + 1], maparray[Q + i]);
        }
    }
}

/**
 * Maps interior modes of an edge to the elemental modes.
 */
void StdTetExp::v_GetEdgeInteriorToElementMap(
    const int eid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation edgeOrient)
{
    int i;
    const int P = m_base[0]->GetNumModes();
    const int Q = m_base[1]->GetNumModes();
    const int R = m_base[2]->GetNumModes();

    const int nEdgeIntCoeffs = v_GetEdgeNcoeffs(eid) - 2;

    if (maparray.size() != nEdgeIntCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
    }
    else
    {
        fill(maparray.data(), maparray.data() + nEdgeIntCoeffs, 0);
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
            for (i = 0; i < P - 2; ++i)
            {
                maparray[i] = GetMode(i + 2, 0, 0);
            }
            if (edgeOrient == eBackwards)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
            break;
        case 1:
            for (i = 0; i < Q - 2; ++i)
            {
                maparray[i] = GetMode(1, i + 1, 0);
            }
            if (edgeOrient == eBackwards)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
            break;
        case 2:
            for (i = 0; i < Q - 2; ++i)
            {
                maparray[i] = GetMode(0, i + 2, 0);
            }
            if (edgeOrient == eBackwards)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
            break;
        case 3:
            for (i = 0; i < R - 2; ++i)
            {
                maparray[i] = GetMode(0, 0, i + 2);
            }
            if (edgeOrient == eBackwards)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
            break;
        case 4:
            for (i = 0; i < R - 2; ++i)
            {
                maparray[i] = GetMode(1, 0, i + 1);
            }
            if (edgeOrient == eBackwards)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
            break;
        case 5:
            for (i = 0; i < R - 2; ++i)
            {
                maparray[i] = GetMode(0, 1, i + 1);
            }
            if (edgeOrient == eBackwards)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
            break;
        default:
            ASSERTL0(false, "Edge not defined.");
            break;
    }
}

void StdTetExp::v_GetTraceInteriorToElementMap(
    const int fid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation faceOrient)
{
    int i, j, idx, k;
    const int P = m_base[0]->GetNumModes();
    const int Q = m_base[1]->GetNumModes();
    const int R = m_base[2]->GetNumModes();

    const int nFaceIntCoeffs = v_GetTraceIntNcoeffs(fid);

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

    switch (fid)
    {
        case 0:
            idx = 0;
            for (i = 2; i < P; ++i)
            {
                for (j = 1; j < Q - i; ++j)
                {
                    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
                    {
                        signarray[idx] = (i % 2 ? -1 : 1);
                    }
                    maparray[idx++] = GetMode(i, j, 0);
                }
            }
            break;
        case 1:
            idx = 0;
            for (i = 2; i < P; ++i)
            {
                for (k = 1; k < R - i; ++k)
                {
                    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
                    {
                        signarray[idx] = (i % 2 ? -1 : 1);
                    }
                    maparray[idx++] = GetMode(i, 0, k);
                }
            }
            break;
        case 2:
            idx = 0;
            for (j = 1; j < Q - 1; ++j)
            {
                for (k = 1; k < R - 1 - j; ++k)
                {
                    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
                    {
                        signarray[idx] = ((j + 1) % 2 ? -1 : 1);
                    }
                    maparray[idx++] = GetMode(1, j, k);
                }
            }
            break;
        case 3:
            idx = 0;
            for (j = 2; j < Q; ++j)
            {
                for (k = 1; k < R - j; ++k)
                {
                    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
                    {
                        signarray[idx] = (j % 2 ? -1 : 1);
                    }
                    maparray[idx++] = GetMode(0, j, k);
                }
            }
            break;
        default:
            ASSERTL0(false, "Face interior map not available.");
            break;
    }
}
//---------------------------------------
// Wrapper functions
//---------------------------------------
DNekMatSharedPtr StdTetExp::v_GenMatrix(const StdMatrixKey &mkey)
{

    MatrixType mtype = mkey.GetMatrixType();

    DNekMatSharedPtr Mat;

    switch (mtype)
    {
        case ePhysInterpToEquiSpaced:
        {
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            int nq2 = m_base[2]->GetNumPoints();
            int nq;

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = max(nq0, max(nq1, nq2));
            }

            int neq =
                LibUtilities::StdTetData::getNumberOfCoefficients(nq, nq, nq);
            Array<OneD, Array<OneD, NekDouble>> coords(neq);
            Array<OneD, NekDouble> coll(3);
            Array<OneD, DNekMatSharedPtr> I(3);
            Array<OneD, NekDouble> tmp(nq0);

            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1 * nq2);
            int cnt = 0;

            for (int i = 0; i < nq; ++i)
            {
                for (int j = 0; j < nq - i; ++j)
                {
                    for (int k = 0; k < nq - i - j; ++k, ++cnt)
                    {
                        coords[cnt]    = Array<OneD, NekDouble>(3);
                        coords[cnt][0] = -1.0 + 2 * k / (NekDouble)(nq - 1);
                        coords[cnt][1] = -1.0 + 2 * j / (NekDouble)(nq - 1);
                        coords[cnt][2] = -1.0 + 2 * i / (NekDouble)(nq - 1);
                    }
                }
            }

            for (int i = 0; i < neq; ++i)
            {
                LocCoordToLocCollapsed(coords[i], coll);

                I[0] = m_base[0]->GetI(coll);
                I[1] = m_base[1]->GetI(coll + 1);
                I[2] = m_base[2]->GetI(coll + 2);

                // interpolate first coordinate direction
                NekDouble fac;
                for (int k = 0; k < nq2; ++k)
                {
                    for (int j = 0; j < nq1; ++j)
                    {

                        fac = (I[1]->GetPtr())[j] * (I[2]->GetPtr())[k];
                        Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                        Vmath::Vcopy(nq0, &tmp[0], 1,
                                     Mat->GetRawPtr() + k * nq0 * nq1 * neq +
                                         j * nq0 * neq + i,
                                     neq);
                    }
                }
            }
        }
        break;
        case ePhysInterpToGLL:
        {
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            int nq2 = m_base[2]->GetNumPoints();
            int nq;

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = max(nq0, max(nq1, nq2));
            }

            int neq =
                LibUtilities::StdTetData::getNumberOfCoefficients(nq, nq, nq);
            Array<OneD, NekDouble> coords(3);
            Array<OneD, NekDouble> coll(3);
            Array<OneD, DNekMatSharedPtr> I(3);
            Array<OneD, NekDouble> tmp(nq0);

            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1 * nq2);

            const LibUtilities::PointsKey key(nq, LibUtilities::eNodalTetElec);

            Array<OneD, const NekDouble> x, y, z;
            LibUtilities::PointsManager()[key]->GetPoints(x, y, z);

            Array<OneD, int> sorted;
            LibUtilities::NodalUtilTetrahedron::CartesianOrdering(nq, sorted);

            for (int i = 0; i < neq; ++i)
            {
                coords[0] = x[sorted[i]];
                coords[1] = y[sorted[i]];
                coords[2] = z[sorted[i]];

                LocCoordToLocCollapsed(coords, coll);

                I[0] = m_base[0]->GetI(coll);
                I[1] = m_base[1]->GetI(coll + 1);
                I[2] = m_base[2]->GetI(coll + 2);

                // interpolate first coordinate direction
                NekDouble fac;
                for (int k = 0; k < nq2; ++k)
                {
                    for (int j = 0; j < nq1; ++j)
                    {

                        fac = (I[1]->GetPtr())[j] * (I[2]->GetPtr())[k];
                        Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                        Vmath::Vcopy(nq0, &tmp[0], 1,
                                     Mat->GetRawPtr() + k * nq0 * nq1 * neq +
                                         j * nq0 * neq + i,
                                     neq);
                    }
                }
            }
            // need to set up test?
        }
        break;
        default:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
        }
        break;
    }

    return Mat;
}

DNekMatSharedPtr StdTetExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

//---------------------------------------
// Private helper functions
//---------------------------------------

/**
 * @brief Compute the mode number in the expansion for a particular
 * tensorial combination.
 *
 * Modes are numbered with the r index travelling fastest, followed by
 * q and then p, and each q-r plane is of size
 * (Q+1)*(Q+2)/2+max(0,R-Q-p)*Q. For example, when P=2, Q=3 and R=4 (nm0=3, nm1
 * = 4, nm2 = 5) the indexing inside each q-r plane (with r increasing upwards
 * and q to the right) is:
 *
 * 4
 * 3 8         17
 * 2 7 11      16 20         25
 * 1 6 10 13   15 19 22      24 27
 * 0 5  9 12   14 18 21      23 26
 *
 * Geometrically they can be interpreted as
 * p = 0:      p = 2:       p = 1:
 * ----------------------------------
 * 1
 * 4 8                      17
 * 3 11 7       25          16 20
 * 2 10 13 6    24 27       15 19 22
 * 0 9  12 5    23 26       14 18 21
 *
 * so we have the following breakdown
 *
 * Vertices V[0,1,2,3] = [0, 14, 5, 1]
 * Edges E[0,1,2,3,4,5,6] =[[23],[18, 21],
 * [9, 12], [2,3,4], [15, 16, 17], [6,7,8]]
 * Faces F[0.1,2,3] = [[26], [24,25], [19, 22, 20], [10, 13, 11]
 * Interior [27]
 * Note that in this element, we must have that \f$ P \leq Q \leq
 * R\f$.
 */
int StdTetExp::GetMode(const int I, const int J, const int K)
{
    const int Q = m_base[1]->GetNumModes();
    const int R = m_base[2]->GetNumModes();

    int i, j, q_hat, k_hat;
    int cnt = 0;

    // Traverse to q-r plane number I
    for (i = 0; i < I; ++i)
    {
        // Size of triangle part
        q_hat = Q - i;
        // Size of rectangle part
        k_hat = R - Q;
        cnt += q_hat * (q_hat + 1) / 2 + k_hat * (Q - i);
    }

    // Traverse to q column J
    q_hat = R - I;
    for (j = 0; j < J; ++j)
    {
        cnt += q_hat;
        q_hat--;
    }

    // Traverse up stacks to K
    cnt += K;

    return cnt;
}

void StdTetExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                     const StdMatrixKey &mkey)
{
    // To do : 1) add a test to ensure 0 \leq SvvCutoff \leq 1.
    //        2) check if the transfer function needs an analytical
    //           Fourier transform.
    //        3) if it doesn't : find a transfer function that renders
    //           the if( cutoff_a ...) useless to reduce computational
    //           cost.
    //        4) add SVVDiffCoef to both models!!

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
    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B, nmodes_b, pb);
    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_C, nmodes_c, pc);

    StdTetExp OrthoExp(Ba, Bb, Bc);

    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
    int i, j, k, cnt = 0;

    // project onto physical space.
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
            for (j = 0; j < nmodes_b - j; ++j)
            {
                NekDouble fac1 = std::max(
                    pow((1.0 * i) / (nmodes_a - 1), cutoff * nmodes_a),
                    pow((1.0 * j) / (nmodes_b - 1), cutoff * nmodes_b));

                for (k = 0; k < nmodes_c - i - j; ++k)
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
            for (j = 0; j < nmodes_b - j; ++j)
            {
                int maxij = max(i, j);

                for (k = 0; k < nmodes_c - i - j; ++k)
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

        // SVV filter paramaters (how much added diffusion
        // relative to physical one and fraction of modes from
        // which you start applying this added diffusion)

        NekDouble SvvDiffCoeff =
            mkey.GetConstFactor(StdRegions::eFactorSVVDiffCoeff);
        NekDouble SVVCutOff =
            mkey.GetConstFactor(StdRegions::eFactorSVVCutoffRatio);

        // Defining the cut of mode
        int cutoff_a      = (int)(SVVCutOff * nmodes_a);
        int cutoff_b      = (int)(SVVCutOff * nmodes_b);
        int cutoff_c      = (int)(SVVCutOff * nmodes_c);
        int nmodes        = min(min(nmodes_a, nmodes_b), nmodes_c);
        NekDouble cutoff  = min(min(cutoff_a, cutoff_b), cutoff_c);
        NekDouble epsilon = 1;

        //------"New" Version August 22nd '13--------------------
        for (i = 0; i < nmodes_a; ++i)
        {
            for (j = 0; j < nmodes_b - i; ++j)
            {
                for (k = 0; k < nmodes_c - i - j; ++k)
                {
                    if (i + j + k >= cutoff)
                    {
                        orthocoeffs[cnt] *= ((SvvDiffCoeff)*exp(
                            -(i + j + k - nmodes) * (i + j + k - nmodes) /
                            ((NekDouble)((i + j + k - cutoff + epsilon) *
                                         (i + j + k - cutoff + epsilon)))));
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

void StdTetExp::v_ReduceOrderCoeffs(int numMin,
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
    Array<OneD, NekDouble> coeff_tmp2(m_ncoeffs, 0.0);
    Array<OneD, NekDouble> phys_tmp(nqtot, 0.0);
    Array<OneD, NekDouble> tmp, tmp2, tmp3, tmp4;

    Vmath::Vcopy(m_ncoeffs, inarray, 1, coeff_tmp2, 1);

    const LibUtilities::PointsKey Pkey0 = m_base[0]->GetPointsKey();
    const LibUtilities::PointsKey Pkey1 = m_base[1]->GetPointsKey();
    const LibUtilities::PointsKey Pkey2 = m_base[2]->GetPointsKey();

    LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A, nmodes0, Pkey0);
    LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_B, nmodes1, Pkey1);
    LibUtilities::BasisKey bortho2(LibUtilities::eOrtho_C, nmodes2, Pkey2);

    Vmath::Zero(m_ncoeffs, coeff_tmp2, 1);

    StdRegions::StdTetExpSharedPtr OrthoTetExp;
    OrthoTetExp = MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(
        bortho0, bortho1, bortho2);

    BwdTrans(inarray, phys_tmp);
    OrthoTetExp->FwdTrans(phys_tmp, coeff);

    Vmath::Zero(m_ncoeffs, outarray, 1);

    // filtering
    int cnt = 0;
    for (int u = 0; u < numMin; ++u)
    {
        for (int i = 0; i < numMin - u; ++i)
        {
            Vmath::Vcopy(numMin - u - i, tmp = coeff + cnt, 1,
                         tmp2 = coeff_tmp1 + cnt, 1);
            cnt += numMax - u - i;
        }
        for (int i = numMin; i < numMax - u; ++i)
        {
            cnt += numMax - u - i;
        }
    }

    OrthoTetExp->BwdTrans(coeff_tmp1, phys_tmp);
    FwdTrans(phys_tmp, outarray);
}

void StdTetExp::v_GetSimplexEquiSpacedConnectivity(
    Array<OneD, int> &conn, [[maybe_unused]] bool standard)
{
    int np0 = m_base[0]->GetNumPoints();
    int np1 = m_base[1]->GetNumPoints();
    int np2 = m_base[2]->GetNumPoints();
    int np  = max(np0, max(np1, np2));

    conn = Array<OneD, int>(4 * (np - 1) * (np - 1) * (np - 1));

    int row     = 0;
    int rowp1   = 0;
    int plane   = 0;
    int row1    = 0;
    int row1p1  = 0;
    int planep1 = 0;
    int cnt     = 0;
    for (int i = 0; i < np - 1; ++i)
    {
        planep1 += (np - i) * (np - i + 1) / 2;
        row    = 0; // current plane row offset
        rowp1  = 0; // current plane row plus one offset
        row1   = 0; // next plane row offset
        row1p1 = 0; // nex plane row plus one offset
        for (int j = 0; j < np - i - 1; ++j)
        {
            rowp1 += np - i - j;
            row1p1 += np - i - j - 1;
            for (int k = 0; k < np - i - j - 2; ++k)
            {
                conn[cnt++] = plane + row + k + 1;
                conn[cnt++] = plane + row + k;
                conn[cnt++] = plane + rowp1 + k;
                conn[cnt++] = planep1 + row1 + k;

                conn[cnt++] = plane + row + k + 1;
                conn[cnt++] = plane + rowp1 + k + 1;
                conn[cnt++] = planep1 + row1 + k + 1;
                conn[cnt++] = planep1 + row1 + k;

                conn[cnt++] = plane + rowp1 + k + 1;
                conn[cnt++] = plane + row + k + 1;
                conn[cnt++] = plane + rowp1 + k;
                conn[cnt++] = planep1 + row1 + k;

                conn[cnt++] = planep1 + row1 + k;
                conn[cnt++] = planep1 + row1p1 + k;
                conn[cnt++] = plane + rowp1 + k;
                conn[cnt++] = plane + rowp1 + k + 1;

                conn[cnt++] = planep1 + row1 + k;
                conn[cnt++] = planep1 + row1p1 + k;
                conn[cnt++] = planep1 + row1 + k + 1;
                conn[cnt++] = plane + rowp1 + k + 1;

                if (k < np - i - j - 3)
                {
                    conn[cnt++] = plane + rowp1 + k + 1;
                    conn[cnt++] = planep1 + row1p1 + k + 1;
                    conn[cnt++] = planep1 + row1 + k + 1;
                    conn[cnt++] = planep1 + row1p1 + k;
                }
            }

            conn[cnt++] = plane + row + np - i - j - 1;
            conn[cnt++] = plane + row + np - i - j - 2;
            conn[cnt++] = plane + rowp1 + np - i - j - 2;
            conn[cnt++] = planep1 + row1 + np - i - j - 2;

            row += np - i - j;
            row1 += np - i - j - 1;
        }
        plane += (np - i) * (np - i + 1) / 2;
    }
}

} // namespace Nektar::StdRegions
