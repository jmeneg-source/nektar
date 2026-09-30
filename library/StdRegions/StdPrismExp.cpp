///////////////////////////////////////////////////////////////////////////////
//
// File: StdPrismExp.cpp
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
// Description: Prismatic routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <StdRegions/StdPrismExp.h>

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

StdPrismExp::StdPrismExp(const LibUtilities::BasisKey &Ba,
                         const LibUtilities::BasisKey &Bb,
                         const LibUtilities::BasisKey &Bc)
    : StdExpansion(LibUtilities::StdPrismData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                   3, Ba, Bb, Bc),
      StdExpansion3D(LibUtilities::StdPrismData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                     Ba, Bb, Bc)
{
    ASSERTL0(Ba.GetNumModes() <= Bc.GetNumModes(),
             "order in 'a' direction is higher than order in 'c' direction");

    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());
    m_weights.push_back(m_base[1]->GetW());

    StdFacKey w2key(eWeights2, Bc);
    // get weights[2] from manager where points are rescaled
    m_weights.push_back(GetStdFac(w2key));
}

//---------------------------------------
// Differentiation Methods
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

void StdPrismExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &u_physical,
                                 Array<OneD, NekDouble> &out_dxi1,
                                 Array<OneD, NekDouble> &out_dxi2,
                                 Array<OneD, NekDouble> &out_dxi3)
{
    int Qx   = m_base[0]->GetNumPoints();
    int Qy   = m_base[1]->GetNumPoints();
    int Qz   = m_base[2]->GetNumPoints();
    int Qtot = Qx * Qy * Qz;

    Array<OneD, NekDouble> dEta_bar1(Qtot, 0.0);

    Array<OneD, const NekDouble> eta_x, eta_z;
    eta_x = m_base[0]->GetZ();
    eta_z = m_base[2]->GetZ();

    int i, k;

    bool Do_1 = (out_dxi1.size() > 0) ? true : false;
    bool Do_3 = (out_dxi3.size() > 0) ? true : false;

    // out_dXi2 is just a tensor derivative so is just passed through
    if (Do_3)
    {
        PhysTensorDeriv(u_physical, dEta_bar1, out_dxi2, out_dxi3);
    }
    else if (Do_1)
    {
        PhysTensorDeriv(u_physical, dEta_bar1, out_dxi2, NullNekDouble1DArray);
    }
    else // case if just require 2nd direction
    {
        PhysTensorDeriv(u_physical, NullNekDouble1DArray, out_dxi2,
                        NullNekDouble1DArray);
    }

    if (Do_1)
    {
        for (k = 0; k < Qz; ++k)
        {
            Vmath::Smul(Qx * Qy, 2.0 / (1.0 - eta_z[k]),
                        &dEta_bar1[0] + k * Qx * Qy, 1,
                        &out_dxi1[0] + k * Qx * Qy, 1);
        }
    }

    if (Do_3)
    {
        // divide dEta_Bar1 by (1-eta_z)
        for (k = 0; k < Qz; ++k)
        {
            Vmath::Smul(Qx * Qy, 1.0 / (1.0 - eta_z[k]),
                        &dEta_bar1[0] + k * Qx * Qy, 1,
                        &dEta_bar1[0] + k * Qx * Qy, 1);
        }

        // Multiply dEta_Bar1 by (1+eta_x) and add ot out_dxi3
        for (i = 0; i < Qx; ++i)
        {
            Vmath::Svtvp(Qz * Qy, 1.0 + eta_x[i], &dEta_bar1[0] + i, Qx,
                         &out_dxi3[0] + i, Qx, &out_dxi3[0] + i, Qx);
        }
    }
}

//---------------------------------------
// Transforms
//---------------------------------------

/**
 * @note 'r' (base[2]) runs fastest in this element.
 *
 * Perform backwards transformation at the quadrature points:
 *
 * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
 *  u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
 *
 * In the prism this expansion becomes:
 *
 * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
 *  (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
 *  \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pr}^b (\xi_{3k})
 *  \rbrace} \rbrace}. \f$
 *
 * And sumfactorizing step of the form is as:\\
 *
 * \f$ f_{pr} (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pr}^b
 * (\xi_{3k}),\\
 *
 * g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y} \psi_{p}^a (\xi_{2j})
 * f_{pr} (\xi_{3k}),\ \
 *
 * u(\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a
 *  (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).  \f$
 */
void StdPrismExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
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

    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nquad2  = m_base[2]->GetNumPoints();
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
    BwdTransPrismKernel(                                                       \
        nmodes0, nmodes1, nmodes2, nquad0, nquad1, nquad2, isModified,         \
        (const vec_t *)base0.data(), (const vec_t *)base1.data(),              \
        (const vec_t *)base2.data(), wsp0.data(), wsp1.data(),                 \
        (const vec_t *)inarray.data(), (vec_t *)outarray.data())

    // Inner loop case over quarature points
#undef BWDTRANS_Q
#define BWDTRANS_Q(r, i)                                                       \
    case NQ(i):                                                                \
        BwdTransPrismKernel(                                                   \
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
void StdPrismExp::v_IProductWRTBaseKernel(
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
        wsp1(nquad2), wsp2(order1);

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
    IProductPrismKernel<false, false, true>(                                   \
        order0, order1, order2, nquad0, nquad1, nquad2, isModified,            \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)base2.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(),                                    \
        (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)wsp2.data(),      \
        (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductPrismKernel<false, false, true>(                               \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ_M1(i), isModified,           \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)base2.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(),                                \
            (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)wsp2.data(),  \
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
    IProductPrismKernel<false, false, false>(                                  \
        order0, order1, order2, nquad0, nquad1, nquad2, isModified,            \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)base2.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(),                                    \
        (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)wsp2.data(),      \
        (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductPrismKernel<false, false, false>(                              \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ_M1(i), isModified,           \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)base2.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(),                                \
            (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)wsp2.data(),  \
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

/**
 * \brief Inner product of \a inarray over region with respect to the
 * object's default expansion basis; output in \a outarray.
 */
void StdPrismExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    ASSERTL0(dir >= 0 && dir <= 2, "input dir is out of range");

    int i;
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    Array<OneD, NekDouble> tmp0(nquad0 * nquad1 * nquad2);

    StdFacKey fackey2(eTwoOverOneMinusZ2, m_base[2]->GetBasisKey());
    Array<OneD, const NekDouble> gfac2 = GetStdFac(fackey2);

    const Array<OneD, const NekDouble> one(1, 1.0);

    // Scale first derivative term by gfac2.
    if (dir != 1)
    {
        for (i = 0; i < nquad2; ++i)
        {
            Vmath::Smul(nquad0 * nquad1, gfac2[i],
                        &inarray[0] + i * nquad0 * nquad1, 1,
                        &tmp0[0] + i * nquad0 * nquad1, 1);
        }
    }

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
                m_base[2]->GetBdata(), inarray, outarray, one, false);
        }
        break;
        case 2:
        {
            Array<OneD, NekDouble> tmp1(m_ncoeffs);

            StdFacKey fackey0(eHalfMultOnePlusZ0, m_base[0]->GetBasisKey());
            Array<OneD, const NekDouble> gfac0 = GetStdFac(fackey0);

            // Scale eta_1 derivative with gfac0.
            for (i = 0; i < nquad1 * nquad2; ++i)
            {
                Vmath::Vmul(nquad0, &gfac0[0], 1, &tmp0[0] + i * nquad0, 1,
                            &tmp0[0] + i * nquad0, 1);
            }

            v_IProductWRTBaseKernel(
                m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                m_base[2]->GetBdata(), tmp0, tmp1, one, false);

            v_IProductWRTBaseKernel(
                m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                m_base[2]->GetDbdata(), inarray, outarray, one, false);

            Vmath::Vadd(m_ncoeffs, &tmp1[0], 1, &outarray[0], 1, &outarray[0],
                        1);
            break;
        }
    }
}

//---------------------------------------
// Evaluation functions
//---------------------------------------

void StdPrismExp::v_LocCoordToLocCollapsed(
    const Array<OneD, const NekDouble> &xi, Array<OneD, NekDouble> &eta)
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
    eta[1] = xi[1]; // eta_y = xi_y
    eta[0] = 2.0 * (1.0 + xi[0]) / d2 - 1.0;
}

void StdPrismExp::v_LocCollapsedToLocCoord(
    const Array<OneD, const NekDouble> &eta, Array<OneD, NekDouble> &xi)
{
    xi[0] = (1.0 + eta[0]) * (1.0 - eta[2]) * 0.5 - 1.0;
    xi[1] = eta[1];
    xi[2] = eta[2];
}

void StdPrismExp::v_GetCoords(Array<OneD, NekDouble> &xi_x,
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
                int s   = i + Qx * (j + Qy * k);
                xi_x[s] = (1.0 - eta_z[k]) * (1.0 + etaBar_x[i]) / 2.0 - 1.0;
                xi_y[s] = eta_y[j];
                xi_z[s] = eta_z[k];
            }
        }
    }
}

NekDouble StdPrismExp::v_PhysEvaluateBasis(
    const Array<OneD, const NekDouble> &coords, int mode)
{
    Array<OneD, NekDouble> coll(3);
    LocCoordToLocCollapsed(coords, coll);

    const int nm1 = m_base[1]->GetNumModes();
    const int nm2 = m_base[2]->GetNumModes();
    const int b   = 2 * nm2 + 1;

    const int mode0 = floor(0.5 * (b - sqrt(b * b - 8.0 * mode / nm1)));
    const int tmp =
        mode - nm1 * (mode0 * (nm2 - 1) + 1 - (mode0 - 2) * (mode0 - 1) / 2);
    const int mode1 = tmp / (nm2 - mode0);
    const int mode2 = tmp % (nm2 - mode0);

    if (mode0 == 0 && mode2 == 1 &&
        m_base[0]->GetBasisType() == LibUtilities::eModified_A)
    {
        // handle collapsed top edge to remove mode0 terms
        return StdExpansion::BaryEvaluateBasis<1>(coll[1], mode1) *
               StdExpansion::BaryEvaluateBasis<2>(coll[2], mode2);
    }
    else
    {
        return StdExpansion::BaryEvaluateBasis<0>(coll[0], mode0) *
               StdExpansion::BaryEvaluateBasis<1>(coll[1], mode1) *
               StdExpansion::BaryEvaluateBasis<2>(coll[2], mode2);
    }
}

NekDouble StdPrismExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    // Collapse coordinates
    Array<OneD, NekDouble> coll(3, 0.0);
    LocCoordToLocCollapsed(coord, coll);

    // If near singularity do the old interpolation matrix method
    // @TODO: Dave thinks there might be a way in the Barycentric to
    //        mathematically remove this singularity?
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

    NekDouble val = BaryTensorDeriv(coll, inarray, firstOrderDerivs);

    NekDouble dEta_bar1 = firstOrderDerivs[0];

    NekDouble fac       = 2.0 / (1.0 - coll[2]);
    firstOrderDerivs[0] = fac * dEta_bar1;

    // divide dEta_Bar1 by (1-eta_z)
    fac       = 1.0 / (1.0 - coll[2]);
    dEta_bar1 = fac * dEta_bar1;

    // Multiply dEta_Bar1 by (1+eta_x) and add ot out_dxi3
    fac = 1.0 + coll[0];
    firstOrderDerivs[2] += fac * dEta_bar1;

    return val;
}

void StdPrismExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp(m_ncoeffs, 0.0);
    tmp[mode] = 1.0;
    StdPrismExp::v_BwdTrans(tmp, outarray);
}

void StdPrismExp::v_GetTraceNumModes(const int fid, int &numModes0,
                                     int &numModes1, Orientation faceOrient)
{
    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};
    switch (fid)
    {
        // base quad
        case 0:
        {
            numModes0 = nummodes[0];
            numModes1 = nummodes[1];
        }
        break;
        // front and back quad
        case 2:
        case 4:
        {
            numModes0 = nummodes[1];
            numModes1 = nummodes[2];
        }
        break;
        // triangles
        case 1:
        case 3:
        {
            numModes0 = nummodes[0];
            numModes1 = nummodes[2];
        }
        break;
    }

    if (faceOrient >= eDir1FwdDir2_Dir2FwdDir1)
    {
        std::swap(numModes0, numModes1);
    }
}

int StdPrismExp::v_GetEdgeNcoeffs(const int i) const
{
    ASSERTL2(i >= 0 && i <= 8, "edge id is out of range");

    if (i == 0 || i == 2)
    {
        return GetBasisNumModes(0);
    }
    else if (i == 1 || i == 3 || i == 8)
    {
        return GetBasisNumModes(1);
    }
    else
    {
        return GetBasisNumModes(2);
    }
}

//---------------------------------------
// Helper functions
//---------------------------------------

int StdPrismExp::v_GetNverts() const
{
    return 6;
}

int StdPrismExp::v_GetNedges() const
{
    return 9;
}

int StdPrismExp::v_GetNtraces() const
{
    return 5;
}

/**
 * \brief Return Shape of region, using ShapeType enum list;
 * i.e. prism.
 */
LibUtilities::ShapeType StdPrismExp::v_DetShapeType() const
{
    return LibUtilities::ePrism;
}

int StdPrismExp::v_NumBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes();
    int Q = m_base[1]->GetNumModes();
    int R = m_base[2]->GetNumModes();

    return LibUtilities::StdPrismData::getNumberOfBndCoefficients(P, Q, R);
}

int StdPrismExp::v_NumDGBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int P = m_base[0]->GetNumModes() - 1;
    int Q = m_base[1]->GetNumModes() - 1;
    int R = m_base[2]->GetNumModes() - 1;

    return (P + 1) * (Q + 1)                    // 1 rect. face on base
           + 2 * (Q + 1) * (R + 1)              // other 2 rect. faces
           + 2 * (R + 1) + P * (1 + 2 * R - P); // 2 tri. faces
}

int StdPrismExp::v_GetTraceNcoeffs(const int i) const
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
        return GetBasisNumModes(1) * GetBasisNumModes(2);
    }
}

int StdPrismExp::v_GetTraceIntNcoeffs(const int i) const
{
    ASSERTL2(i >= 0 && i <= 4, "face id is out of range");

    int Pi = GetBasisNumModes(0) - 2;
    int Qi = GetBasisNumModes(1) - 2;
    int Ri = GetBasisNumModes(2) - 2;

    if (i == 0)
    {
        return Pi * Qi;
    }
    else if (i == 1 || i == 3)
    {
        return Pi * (2 * Ri - Pi - 1) / 2;
    }
    else
    {
        return Qi * Ri;
    }
}

int StdPrismExp::v_GetTraceNumPoints(const int i) const
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

LibUtilities::PointsKey StdPrismExp::v_GetTracePointsKey(const int i,
                                                         const int j) const
{
    ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
    ASSERTL2(j == 0 || j == 1, "face direction is out of range");

    if (i == 0)
    {
        return m_base[j]->GetPointsKey();
    }
    else if (i == 1 || i == 3)
    {
        return m_base[2 * j]->GetPointsKey();
    }
    else
    {
        return m_base[j + 1]->GetPointsKey();
    }
}

const LibUtilities::BasisKey StdPrismExp::v_GetTraceBasisKey(const int i,
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
        case 2:
        case 4:
        {
            return EvaluateQuadFaceBasisKey(k, m_base[k + 1]);
        }
        case 1:
        case 3:
        {
            return EvaluateTriFaceBasisKey(k, m_base[2 * k], UseGLL);
        }
        break;
    }

    // Should never get here.
    return LibUtilities::NullBasisKey;
}

int StdPrismExp::v_CalcNumberOfCoefficients(
    const std::vector<unsigned int> &nummodes, int &modes_offset)
{
    int nmodes = LibUtilities::StdPrismData::getNumberOfCoefficients(
        nummodes[modes_offset], nummodes[modes_offset + 1],
        nummodes[modes_offset + 2]);

    modes_offset += 3;
    return nmodes;
}

bool StdPrismExp::v_IsBoundaryInteriorExpansion() const
{
    return (m_base[0]->GetBasisType() == LibUtilities::eModified_A) &&
           (m_base[1]->GetBasisType() == LibUtilities::eModified_A) &&
           (m_base[2]->GetBasisType() == LibUtilities::eModified_B);
}

//---------------------------------------
// Mappings
//---------------------------------------

int StdPrismExp::v_GetVertexMap(const int vId, bool useCoeffPacking)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eModified_B,
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
                l = GetMode(0, 1, 1);
                break;
            case 4:
                l = GetMode(1, 0, 0);
                break;
            case 5:
                l = GetMode(1, 1, 0);
                break;
            default:
                ASSERTL0(false, "local vertex id must be between 0 and 5");
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
            case 5:
                l = GetMode(0, 1, 1);
                break;
            default:
                ASSERTL0(false, "local vertex id must be between 0 and 5");
        }
    }

    return l;
}

void StdPrismExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
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
            for (r = 1; r <= R - p; ++r)
            {
                outarray[idx++] = GetMode(p, q, r);
            }
        }
    }
}

void StdPrismExp::v_GetBoundaryMap(Array<OneD, unsigned int> &maparray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
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
                for (r = 0; r <= R - p; ++r)
                {
                    maparray[idx++] = GetMode(p, q, r);
                }
            }
        }
        else
        {
            // Remaining q-r planes contain boundary modes on the two
            // left-hand sides and bottom edge.
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

void StdPrismExp::v_GetTraceCoeffMap(const unsigned int fid,
                                     Array<OneD, unsigned int> &maparray)
{
    ASSERTL1(GetBasisType(0) == GetBasisType(1),
             "Method only implemented if BasisType is identical"
             "in x and y directions");
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                 GetBasisType(2) == LibUtilities::eModified_B,
             "Method only implemented for Modified_A BasisType"
             "(x and y direction) and Modified_B BasisType (z "
             "direction)");
    int p, q, r, idx = 0;
    int P = 0, Q = 0;

    switch (fid)
    {
        case 0:
            P = m_base[0]->GetNumModes();
            Q = m_base[1]->GetNumModes();
            break;
        case 1:
        case 3:
            P = m_base[0]->GetNumModes();
            Q = m_base[2]->GetNumModes();
            break;
        case 2:
        case 4:
            P = m_base[1]->GetNumModes();
            Q = m_base[2]->GetNumModes();
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
        case 1: // Left triangle
            for (p = 0; p < P; ++p)
            {
                for (r = 0; r < Q - p; ++r)
                {
                    maparray[idx++] = GetMode(p, 0, r);
                }
            }
            break;
        case 2: // Slanted quad
            for (q = 0; q < P; ++q)
            {
                maparray[q] = GetMode(1, q, 0);
            }
            for (q = 0; q < P; ++q)
            {
                maparray[P + q] = GetMode(0, q, 1);
            }
            for (r = 1; r < Q - 1; ++r)
            {
                for (q = 0; q < P; ++q)
                {
                    maparray[(r + 1) * P + q] = GetMode(1, q, r);
                }
            }
            break;
        case 3: // Right triangle
            for (p = 0; p < P; ++p)
            {
                for (r = 0; r < Q - p; ++r)
                {
                    maparray[idx++] = GetMode(p, 1, r);
                }
            }
            break;
        case 4: // Rear quad
            for (r = 0; r < Q; ++r)
            {
                for (q = 0; q < P; ++q)
                {
                    maparray[r * P + q] = GetMode(0, q, r);
                }
            }
            break;
        default:
            ASSERTL0(false, "Face to element map unavailable.");
    }
}

void StdPrismExp::v_GetElmtTraceToTraceMap(const unsigned int fid,
                                           Array<OneD, unsigned int> &maparray,
                                           Array<OneD, int> &signarray,
                                           Orientation faceOrient, int P, int Q)
{
    ASSERTL1(GetBasisType(0) == GetBasisType(1),
             "Method only implemented if BasisType is identical"
             "in x and y directions");
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                 GetBasisType(2) == LibUtilities::eModified_B,
             "Method only implemented for Modified_A BasisType"
             "(x and y direction) and Modified_B BasisType (z "
             "direction)");

    int i, j, k, p, r, nFaceCoeffs, idx = 0;
    int nummodesA = 0, nummodesB = 0;

    switch (fid)
    {
        case 0:
            nummodesA = m_base[0]->GetNumModes();
            nummodesB = m_base[1]->GetNumModes();
            break;
        case 1:
        case 3:
            nummodesA = m_base[0]->GetNumModes();
            nummodesB = m_base[2]->GetNumModes();
            break;
        case 2:
        case 4:
            nummodesA = m_base[1]->GetNumModes();
            nummodesB = m_base[2]->GetNumModes();
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
    else if (fid == 1 || fid == 3)
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

    int minPA = min(nummodesA, P);
    int minQB = min(nummodesB, Q);
    // triangular faces
    if (fid == 1 || fid == 3)
    {
        // zero signmap and set maparray to zero if elemental
        // modes are not as large as face modesl
        idx     = 0;
        int cnt = 0;

        for (j = 0; j < minPA; ++j)
        {
            // set maparray
            for (k = 0; k < minQB - j; ++k, ++cnt)
            {
                maparray[idx++] = cnt;
            }

            cnt += nummodesB - minQB;

            // idx += nummodesB-j;
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

            swap(maparray[0], maparray[Q]);
            for (i = 1; i < Q - 1; ++i)
            {
                swap(maparray[i + 1], maparray[Q + i]);
            }
        }
    }
    else
    {
        // Set up an array indexing for quads, since the
        // ordering may need to be transposed.
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

void StdPrismExp::v_GetEdgeInteriorToElementMap(
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
            // Base quad; reverse direction.
            // signChange = !signChange;
            for (i = 2; i <= P; ++i)
            {
                maparray[i - 2] = GetMode(i, 1, 0);
            }
            break;

        case 3:
            // Base quad; reverse direction.
            // signChange = !signChange;
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
            for (i = 2; i <= R; ++i)
            {
                maparray[i - 2] = GetMode(0, 1, i);
            }
            break;

        case 8:
            for (i = 2; i <= Q; ++i)
            {
                maparray[i - 2] = GetMode(0, i, 1);
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

void StdPrismExp::v_GetTraceInteriorToElementMap(
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
    int i         = 0;
    int j         = 0;

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
    if (fid != 1 && fid != 3)
    {
        if (fid == 0) // Base quad
        {
            nummodesA = P - 1;
            nummodesB = Q - 1;
        }
        else // front and back quad
        {
            nummodesA = Q - 1;
            nummodesB = R - 1;
        }

        for (i = 0; i < nummodesB; i++)
        {
            for (j = 0; j < nummodesA; j++)
            {
                if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
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

        case 1: // Left triangle
            for (p = 2; p <= P; ++p)
            {
                for (r = 1; r <= R - p; ++r)
                {
                    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
                    {
                        signarray[idx] = p % 2 ? -1 : 1;
                    }
                    maparray[idx++] = GetMode(p, 0, r);
                }
            }
            break;

        case 2: // Slanted quad
            for (r = 1; r <= R - 1; ++r)
            {
                for (q = 2; q <= Q; ++q)
                {
                    maparray[arrayindx[(r - 1) * nummodesA + (q - 2)]] =
                        GetMode(1, q, r);
                }
            }
            break;

        case 3: // Right triangle
            for (p = 2; p <= P; ++p)
            {
                for (r = 1; r <= R - p; ++r)
                {
                    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2)
                    {
                        signarray[idx] = p % 2 ? -1 : 1;
                    }
                    maparray[idx++] = GetMode(p, 1, r);
                }
            }
            break;

        case 4: // Back quad
            for (r = 2; r <= R; ++r)
            {
                for (q = 2; q <= Q; ++q)
                {
                    maparray[arrayindx[(r - 2) * nummodesA + (q - 2)]] =
                        GetMode(0, q, r);
                }
            }
            break;

        default:
            ASSERTL0(false, "Face interior map not available.");
    }

    // Triangular faces are processed in the above switch loop; for
    // remaining quad faces, set up orientation if necessary.
    if (fid == 1 || fid == 3)
    {
        return;
    }

    if (faceOrient == eDir1FwdDir1_Dir2BwdDir2 ||
        faceOrient == eDir1BwdDir1_Dir2BwdDir2 ||
        faceOrient == eDir1BwdDir2_Dir2FwdDir1 ||
        faceOrient == eDir1BwdDir2_Dir2BwdDir1)
    {
        if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
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

    if (faceOrient == eDir1BwdDir1_Dir2FwdDir2 ||
        faceOrient == eDir1BwdDir1_Dir2BwdDir2 ||
        faceOrient == eDir1FwdDir2_Dir2BwdDir1 ||
        faceOrient == eDir1BwdDir2_Dir2BwdDir1)
    {
        if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
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

DNekMatSharedPtr StdPrismExp::v_GenMatrix(const StdMatrixKey &mkey)
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
                LibUtilities::StdPrismData::getNumberOfCoefficients(nq, nq, nq);
            Array<OneD, Array<OneD, NekDouble>> coords(neq);
            Array<OneD, NekDouble> coll(3);
            Array<OneD, DNekMatSharedPtr> I(3);
            Array<OneD, NekDouble> tmp(nq0);

            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1 * nq2);
            int cnt = 0;
            for (int i = 0; i < nq; ++i)
            {
                for (int j = 0; j < nq; ++j)
                {
                    for (int k = 0; k < nq - i; ++k, ++cnt)
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
                LibUtilities::StdPrismData::getNumberOfCoefficients(nq, nq, nq);
            Array<OneD, NekDouble> coords(3);
            Array<OneD, NekDouble> coll(3);
            Array<OneD, DNekMatSharedPtr> I(3);
            Array<OneD, NekDouble> tmp(nq0);

            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1 * nq2);

            const LibUtilities::PointsKey key(nq,
                                              LibUtilities::eNodalPrismElec);

            Array<OneD, const NekDouble> x, y, z;
            LibUtilities::PointsManager()[key]->GetPoints(x, y, z);

            Array<OneD, int> sorted;
            LibUtilities::NodalUtilPrism::CartesianOrdering(nq, sorted);

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

DNekMatSharedPtr StdPrismExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

/**
 * @brief Compute the local mode number in the expansion for a
 * particular tensorial combination.
 *
 * Modes are numbered with the r index travelling fastest, followed by
 * q and then p, and each q-r plane is of size (R+1-p). For example,
 * with P=1, Q=2, R=3, the indexing inside each q-r plane (with r
 * increasing upwards and q to the right) is:
 *
 * p = 0:       p = 1:
 * -----------------------
 * 3   7  11
 * 2   6  10    14  17  20
 * 1   5   9    13  16  19
 * 0   4   8    12  15  18
 *
 * Note that in this element, we must have that \f$ P <= R \f$.
 */
int StdPrismExp::GetMode(int p, int q, int r)
{
    int Q = m_base[1]->GetNumModes() - 1;
    int R = m_base[2]->GetNumModes() - 1;

    return r +               // Skip along stacks  (r-direction)
           q * (R + 1 - p) + // Skip along columns (q-direction)
           (Q + 1) * (p * R + 1 -
                      (p - 2) * (p - 1) / 2); // Skip along rows (p-direction)
}

void StdPrismExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
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
    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_B, nmodes_c, pc);
    StdPrismExp OrthoExp(Ba, Bb, Bc);

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
                NekDouble fac1 = std::max(
                    pow((1.0 * i) / (nmodes_a - 1), cutoff * nmodes_a),
                    pow((1.0 * j) / (nmodes_b - 1), cutoff * nmodes_b));

                for (k = 0; k < nmodes_c - i; ++k)
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

                for (k = 0; k < nmodes_c - i; ++k)
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

        //------"New" Version August 22nd '13--------------------
        for (i = 0; i < nmodes_a; ++i) // P
        {
            for (j = 0; j < nmodes_b; ++j) // Q
            {
                for (k = 0; k < nmodes_c - i; ++k) // R
                {
                    if (j >= cutoff || i + k >= cutoff)
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

void StdPrismExp::v_ReduceOrderCoeffs(
    int numMin, const Array<OneD, const NekDouble> &inarray,
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
    LibUtilities::BasisKey bortho2(LibUtilities::eOrtho_B, nmodes2, Pkey2);

    int cnt = 0;
    int u   = 0;
    int i   = 0;
    StdRegions::StdPrismExpSharedPtr OrthoPrismExp;

    OrthoPrismExp = MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(
        bortho0, bortho1, bortho2);

    BwdTrans(inarray, phys_tmp);
    OrthoPrismExp->FwdTrans(phys_tmp, coeff);

    // filtering
    for (u = 0; u < numMin; ++u)
    {
        for (i = 0; i < numMin; ++i)
        {
            Vmath::Vcopy(numMin - u, tmp = coeff + cnt, 1,
                         tmp2 = coeff_tmp1 + cnt, 1);
            cnt += numMax - u;
        }

        for (i = numMin; i < numMax; ++i)
        {
            cnt += numMax - u;
        }
    }

    OrthoPrismExp->BwdTrans(coeff_tmp1, phys_tmp);
    StdPrismExp::FwdTrans(phys_tmp, outarray);
}
} // namespace Nektar::StdRegions
