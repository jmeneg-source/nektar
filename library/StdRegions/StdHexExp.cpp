///////////////////////////////////////////////////////////////////////////////
//
// File: StdHexExp.cpp
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
// Description: Heaxhedral methods
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <StdRegions/StdHexExp.h>

#ifdef max
#undef max
#endif

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

StdHexExp::StdHexExp(const LibUtilities::BasisKey &Ba,
                     const LibUtilities::BasisKey &Bb,
                     const LibUtilities::BasisKey &Bc)
    : StdExpansion(Ba.GetNumModes() * Bb.GetNumModes() * Bc.GetNumModes(), 3,
                   Ba, Bb, Bc),
      StdExpansion3D(Ba.GetNumModes() * Bb.GetNumModes() * Bc.GetNumModes(), Ba,
                     Bb, Bc)
{
    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());

    // cache integration weights for future use
    m_weights.push_back(m_base[1]->GetW());

    // cache integration weights for future use
    m_weights.push_back(m_base[2]->GetW());
}

bool StdHexExp::v_IsBoundaryInteriorExpansion() const
{
    return (m_base[0]->GetBasisType() == LibUtilities::eModified_A &&
            m_base[1]->GetBasisType() == LibUtilities::eModified_A &&
            m_base[2]->GetBasisType() == LibUtilities::eModified_A) ||
           (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange &&
            m_base[1]->GetBasisType() == LibUtilities::eGLL_Lagrange &&
            m_base[1]->GetBasisType() == LibUtilities::eGLL_Lagrange);
}

///////////////////////////////
/// Differentiation Methods
///////////////////////////////
/**
 * For Hexahedral region can use the PhysTensorDeriv function defined
 * under StdExpansion. Following tenserproduct:
 */
void StdHexExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &out_d0,
                               Array<OneD, NekDouble> &out_d1,
                               Array<OneD, NekDouble> &out_d2)
{
    PhysTensorDeriv(inarray, out_d0, out_d1, out_d2);
}

/**
 * Backward transformation is three dimensional tensorial expansion
 * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k})
 *  = \sum_{p=0}^{Q_x} \psi_p^a (\xi_{1i})
 *  \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
 *    \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{r}^a (\xi_{3k})
 *    \rbrace}
 *  \rbrace}. \f$
 * And sumfactorizing step of the form is as:\\
 * \f$ f_{r} (\xi_{3k})
 * = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{r}^a (\xi_{3k}),\\
 * g_{p} (\xi_{2j}, \xi_{3k})
 * = \sum_{r=0}^{Q_y} \psi_{p}^a (\xi_{2j}) f_{r} (\xi_{3k}),\\
 * u(\xi_{1i}, \xi_{2j}, \xi_{3k})
 * = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).
 * \f$
 *
 * @param   inarray     ?
 * @param   outarray    ?
 */
void StdHexExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    if (m_base[0]->Collocation() && m_base[1]->Collocation() &&
        m_base[2]->Collocation())
    {
        std::memcpy(outarray.data(), inarray.data(),
                    nquad0 * nquad1 * nquad2 * sizeof(NekDouble));
    }
    else
    {
        const Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
        const Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();
        const Array<OneD, const NekDouble> base2 = m_base[2]->GetBdata();

        int nmodes0 = m_base[0]->GetNumModes();
        int nmodes1 = m_base[1]->GetNumModes();
        int nmodes2 = m_base[2]->GetNumModes();

        std::vector<vec_t, tinysimd::allocator<vec_t>> wsp0(nmodes1 * nmodes2 *
                                                            nquad0),
            wsp1(nquad1 * nquad0 * nmodes2);

        // Switch statment using boost_pp and macros. This unfolls intwo a
        // nested swtich statement where the outer swtich statement runs
        // from SMIN to SMAX for modal order and the inner switch
        // statemets run from the outer value of the case to 2*SMAX for
        // the quadrature order. If you want to see it unwrapped compile
        // in verbose mode and add --preprocess to the c++ command.
        // Default case
#undef BWDTRANS_DEF
#define BWDTRANS_DEF                                                           \
    BwdTransHexKernel(nmodes0, nmodes1, nmodes2, nquad0, nquad1, nquad2,       \
                      (const vec_t *)base0.data(),                             \
                      (const vec_t *)base1.data(),                             \
                      (const vec_t *)base2.data(), wsp0.data(), wsp1.data(),   \
                      (const vec_t *)inarray.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef BWDTRANS_Q
#define BWDTRANS_Q(r, i)                                                       \
    case NQ(i):                                                                \
        BwdTransHexKernel(                                                     \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ(i),                          \
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
            (nquad0 == nquad1) && (nquad1 == nquad2))
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
}

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
 *  @param CollDir0 - bool to identify if 0-direction basis is a
 *  collocated expansion
 *  @param CollDir1 - bool to identify if 1-direction basis is a
 *  collocated expansion
 *  @param CollDir2 - bool to identify if 2-direction basis is a
 *  collocated expansion
 */
void StdHexExp::v_IProductWRTBaseKernel(
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
    IProductHexKernel<false, false, true>(                                     \
        order0, order1, order2, nquad0, nquad1, nquad2,                        \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)base2.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(),                                    \
        (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)outarray.data(),  \
        1.0, CollDir0, CollDir1, CollDir2)

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductHexKernel<false, false, true>(                                 \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ(i),                          \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)base2.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(),                                \
            (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)wsp1.data(),                        \
            (vec_t *)outarray.data(), 1.0, CollDir0, CollDir1, CollDir2);      \
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
            (nquad1 == nquad2))
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
    IProductHexKernel<false, false, false>(                                    \
        order0, order1, order2, nquad0, nquad1, nquad2,                        \
        (const vec_t *)inarray.data(), (const vec_t *)base0.data(),            \
        (const vec_t *)base1.data(), (const vec_t *)base2.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(),                                    \
        (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)wsp1.data(), (vec_t *)outarray.data(),  \
        1.0, CollDir0, CollDir1, CollDir2)

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductHexKernel<false, false, false>(                                \
            NM(i), NM(i), NM(i), NQ(i), NQ(i), NQ(i),                          \
            (const vec_t *)inarray.data(), (const vec_t *)base0.data(),        \
            (const vec_t *)base1.data(), (const vec_t *)base2.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(),                                \
            (const vec_t *)m_weights[2].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)wsp1.data(),                        \
            (vec_t *)outarray.data(), 1.0, CollDir0, CollDir1, CollDir2);      \
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
            (nquad1 == nquad2))
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

void StdHexExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    ASSERTL0((dir == 0) || (dir == 1) || (dir == 2),
             "input dir is out of range");
    Array<OneD, NekDouble> one(1, 1.0);

    // perform sum-factorisation
    switch (dir)
    {
        case 0:
            v_IProductWRTBaseKernel(
                m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                m_base[2]->GetBdata(), inarray, outarray, one, false, false,
                m_base[1]->Collocation(), m_base[2]->Collocation());
            break;
        case 1:
            v_IProductWRTBaseKernel(
                m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                m_base[2]->GetBdata(), inarray, outarray, one, false,
                m_base[0]->Collocation(), false, m_base[2]->Collocation());
            break;
        case 2:
            v_IProductWRTBaseKernel(
                m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                m_base[2]->GetDbdata(), inarray, outarray, one, false,
                m_base[0]->Collocation(), m_base[1]->Collocation(), false);
            break;
    }
}

void StdHexExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble> &xi,
                                         Array<OneD, NekDouble> &eta)
{
    eta[0] = xi[0];
    eta[1] = xi[1];
    eta[2] = xi[2];
}

void StdHexExp::v_LocCollapsedToLocCoord(
    const Array<OneD, const NekDouble> &eta, Array<OneD, NekDouble> &xi)
{
    xi[0] = eta[0];
    xi[1] = eta[1];
    xi[2] = eta[2];
}

/**
 * @note for hexahedral expansions _base[0] (i.e. p) modes run fastest.
 */
void StdHexExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
    Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();
    Array<OneD, const NekDouble> base2 = m_base[2]->GetBdata();

    int btmp0 = m_base[0]->GetNumModes();
    int btmp1 = m_base[1]->GetNumModes();
    int mode2 = mode / (btmp0 * btmp1);
    int mode1 = (mode - mode2 * btmp0 * btmp1) / btmp0;
    int mode0 = (mode - mode2 * btmp0 * btmp1) % btmp0;

    ASSERTL2(mode == mode2 * btmp0 * btmp1 + mode1 * btmp0 + mode0,
             "Mode lookup failed.");
    ASSERTL2(mode < m_ncoeffs,
             "Calling argument mode is larger than total expansion "
             "order");

    for (int i = 0; i < nquad1 * nquad2; ++i)
    {
        Vmath::Vcopy(nquad0, (NekDouble *)(base0.data() + mode0 * nquad0), 1,
                     &outarray[0] + i * nquad0, 1);
    }

    for (int j = 0; j < nquad2; ++j)
    {
        for (int i = 0; i < nquad0; ++i)
        {
            Vmath::Vmul(nquad1, (NekDouble *)(base1.data() + mode1 * nquad1), 1,
                        &outarray[0] + i + j * nquad0 * nquad1, nquad0,
                        &outarray[0] + i + j * nquad0 * nquad1, nquad0);
        }
    }

    for (int i = 0; i < nquad2; i++)
    {
        Blas::Dscal(nquad0 * nquad1, base2[mode2 * nquad2 + i],
                    &outarray[0] + i * nquad0 * nquad1, 1);
    }
}

NekDouble StdHexExp::v_PhysEvaluateBasis(
    const Array<OneD, const NekDouble> &coords, int mode)
{
    ASSERTL2(coords[0] > -1 - NekConstants::kNekZeroTol, "coord[0] < -1");
    ASSERTL2(coords[0] < 1 + NekConstants::kNekZeroTol, "coord[0] >  1");
    ASSERTL2(coords[1] > -1 - NekConstants::kNekZeroTol, "coord[1] < -1");
    ASSERTL2(coords[1] < 1 + NekConstants::kNekZeroTol, "coord[1] >  1");
    ASSERTL2(coords[2] > -1 - NekConstants::kNekZeroTol, "coord[2] < -1");
    ASSERTL2(coords[2] < 1 + NekConstants::kNekZeroTol, "coord[2] >  1");

    const int nm0   = m_base[0]->GetNumModes();
    const int nm1   = m_base[1]->GetNumModes();
    const int mode2 = mode / (nm0 * nm1);
    const int mode1 = (mode - mode2 * nm0 * nm1) / nm0;
    const int mode0 = (mode - mode2 * nm0 * nm1) % nm0;

    return StdExpansion::BaryEvaluateBasis<0>(coords[0], mode0) *
           StdExpansion::BaryEvaluateBasis<1>(coords[1], mode1) *
           StdExpansion::BaryEvaluateBasis<2>(coords[2], mode2);
}

int StdHexExp::v_GetNverts() const
{
    return 8;
}

int StdHexExp::v_GetNedges() const
{
    return 12;
}

int StdHexExp::v_GetNtraces() const
{
    return 6;
}

LibUtilities::ShapeType StdHexExp::v_DetShapeType() const
{
    return LibUtilities::eHexahedron;
}

int StdHexExp::v_NumBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();
    int nmodes2 = m_base[2]->GetNumModes();

    return (2 * (nmodes0 * nmodes1 + nmodes0 * nmodes2 + nmodes1 * nmodes2) -
            4 * (nmodes0 + nmodes1 + nmodes2) + 8);
}

int StdHexExp::v_NumDGBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();
    int nmodes2 = m_base[2]->GetNumModes();

    return 2 * (nmodes0 * nmodes1 + nmodes0 * nmodes2 + nmodes1 * nmodes2);
}

int StdHexExp::v_GetTraceNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 5), "face id is out of range");
    if ((i == 0) || (i == 5))
    {
        return GetBasisNumModes(0) * GetBasisNumModes(1);
    }
    else if ((i == 1) || (i == 3))
    {
        return GetBasisNumModes(0) * GetBasisNumModes(2);
    }
    else
    {
        return GetBasisNumModes(1) * GetBasisNumModes(2);
    }
}

int StdHexExp::v_GetTraceIntNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 5), "face id is out of range");
    if ((i == 0) || (i == 5))
    {
        return (GetBasisNumModes(0) - 2) * (GetBasisNumModes(1) - 2);
    }
    else if ((i == 1) || (i == 3))
    {
        return (GetBasisNumModes(0) - 2) * (GetBasisNumModes(2) - 2);
    }
    else
    {
        return (GetBasisNumModes(1) - 2) * (GetBasisNumModes(2) - 2);
    }
}

int StdHexExp::v_GetTraceNumPoints(const int i) const
{
    ASSERTL2(i >= 0 && i <= 5, "face id is out of range");

    if (i == 0 || i == 5)
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

LibUtilities::PointsKey StdHexExp::v_GetTracePointsKey(const int i,
                                                       const int j) const
{
    ASSERTL2(i >= 0 && i <= 5, "face id is out of range");
    ASSERTL2(j == 0 || j == 1, "face direction is out of range");

    if (i == 0 || i == 5)
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

int StdHexExp::v_CalcNumberOfCoefficients(
    const std::vector<unsigned int> &nummodes, int &modes_offset)
{
    int nmodes = nummodes[modes_offset] * nummodes[modes_offset + 1] *
                 nummodes[modes_offset + 2];
    modes_offset += 3;

    return nmodes;
}

const LibUtilities::BasisKey StdHexExp::v_GetTraceBasisKey(
    const int i, const int k, [[maybe_unused]] bool UseGLL) const
{
    ASSERTL2(i >= 0 && i <= 5, "face id is out of range");
    ASSERTL2(k >= 0 && k <= 1, "basis key id is out of range");

    int dir = k;
    switch (i)
    {
        case 0:
        case 5:
            dir = k;
            break;
        case 1:
        case 3:
            dir = 2 * k;
            break;
        case 2:
        case 4:
            dir = k + 1;
            break;
    }

    return EvaluateQuadFaceBasisKey(k, m_base[dir]);
}

void StdHexExp::v_GetCoords(Array<OneD, NekDouble> &xi_x,
                            Array<OneD, NekDouble> &xi_y,
                            Array<OneD, NekDouble> &xi_z)
{
    Array<OneD, const NekDouble> eta_x = m_base[0]->GetZ();
    Array<OneD, const NekDouble> eta_y = m_base[1]->GetZ();
    Array<OneD, const NekDouble> eta_z = m_base[2]->GetZ();
    int Qx                             = GetNumPoints(0);
    int Qy                             = GetNumPoints(1);
    int Qz                             = GetNumPoints(2);

    // Convert collapsed coordinates into cartesian coordinates:
    // eta --> xi
    for (int k = 0; k < Qz; ++k)
    {
        for (int j = 0; j < Qy; ++j)
        {
            for (int i = 0; i < Qx; ++i)
            {
                int s   = i + Qx * (j + Qy * k);
                xi_x[s] = eta_x[i];
                xi_y[s] = eta_y[j];
                xi_z[s] = eta_z[k];
            }
        }
    }
}

void StdHexExp::v_GetTraceNumModes(const int fid, int &numModes0,
                                   int &numModes1, Orientation faceOrient)
{
    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};
    switch (fid)
    {
        case 0:
        case 5:
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
        default:
        {
            ASSERTL0(false, "fid out of range");
        }
        break;
    }

    if (faceOrient >= eDir1FwdDir2_Dir2FwdDir1)
    {
        std::swap(numModes0, numModes1);
    }
}

/**
 * Expansions in each of the three dimensions must be of type
 * LibUtilities#eModified_A or LibUtilities#eGLL_Lagrange.
 *
 * @param   localVertexId   ID of vertex (0..7)
 * @returns Position of vertex in local numbering scheme.
 */
int StdHexExp::v_GetVertexMap(const int localVertexId, bool useCoeffPacking)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    ASSERTL1((localVertexId >= 0) && (localVertexId < 8),
             "local vertex id must be between 0 and 7");

    int p = 0;
    int q = 0;
    int r = 0;

    // Retrieve the number of modes in each dimension.
    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};

    if (useCoeffPacking == true) // follow packing of coefficients i.e q,r,p
    {
        if (localVertexId > 3)
        {
            if (GetBasisType(2) == LibUtilities::eGLL_Lagrange)
            {
                r = nummodes[2] - 1;
            }
            else
            {
                r = 1;
            }
        }

        switch (localVertexId % 4)
        {
            case 0:
                break;
            case 1:
            {
                if (GetBasisType(0) == LibUtilities::eGLL_Lagrange)
                {
                    p = nummodes[0] - 1;
                }
                else
                {
                    p = 1;
                }
            }
            break;
            case 2:
            {
                if (GetBasisType(1) == LibUtilities::eGLL_Lagrange)
                {
                    q = nummodes[1] - 1;
                }
                else
                {
                    q = 1;
                }
            }
            break;
            case 3:
            {
                if (GetBasisType(1) == LibUtilities::eGLL_Lagrange)
                {
                    p = nummodes[0] - 1;
                    q = nummodes[1] - 1;
                }
                else
                {
                    p = 1;
                    q = 1;
                }
            }
            break;
        }
    }
    else
    {
        // Right face (vertices 1,2,5,6)
        if ((localVertexId % 4) % 3 > 0)
        {
            if (GetBasisType(0) == LibUtilities::eGLL_Lagrange)
            {
                p = nummodes[0] - 1;
            }
            else
            {
                p = 1;
            }
        }
        // Back face (vertices 2,3,6,7)
        if (localVertexId % 4 > 1)
        {
            if (GetBasisType(1) == LibUtilities::eGLL_Lagrange)
            {
                q = nummodes[1] - 1;
            }
            else
            {
                q = 1;
            }
        }

        // Top face (vertices 4,5,6,7)
        if (localVertexId > 3)
        {
            if (GetBasisType(2) == LibUtilities::eGLL_Lagrange)
            {
                r = nummodes[2] - 1;
            }
            else
            {
                r = 1;
            }
        }
    }
    // Compute the local number.
    return r * nummodes[0] * nummodes[1] + q * nummodes[0] + p;
}

/**
 * @param   outarray    Storage area for computed map.
 */
void StdHexExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int i;
    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};

    int nIntCoeffs = m_ncoeffs - NumBndryCoeffs();

    if (outarray.size() != nIntCoeffs)
    {
        outarray = Array<OneD, unsigned int>(nIntCoeffs);
    }

    const LibUtilities::BasisType Btype[3] = {GetBasisType(0), GetBasisType(1),
                                              GetBasisType(2)};

    int p, q, r;
    int cnt = 0;

    int IntIdx[3][2];

    for (i = 0; i < 3; i++)
    {
        if (Btype[i] == LibUtilities::eModified_A)
        {
            IntIdx[i][0] = 2;
            IntIdx[i][1] = nummodes[i];
        }
        else
        {
            IntIdx[i][0] = 1;
            IntIdx[i][1] = nummodes[i] - 1;
        }
    }

    for (r = IntIdx[2][0]; r < IntIdx[2][1]; r++)
    {
        for (q = IntIdx[1][0]; q < IntIdx[1][1]; q++)
        {
            for (p = IntIdx[0][0]; p < IntIdx[0][1]; p++)
            {
                outarray[cnt++] =
                    r * nummodes[0] * nummodes[1] + q * nummodes[0] + p;
            }
        }
    }
}

/**
 * @param   outarray    Storage for computed map.
 */
void StdHexExp::v_GetBoundaryMap(Array<OneD, unsigned int> &outarray)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    int i;
    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};

    int nBndCoeffs = NumBndryCoeffs();

    if (outarray.size() != nBndCoeffs)
    {
        outarray = Array<OneD, unsigned int>(nBndCoeffs);
    }

    const LibUtilities::BasisType Btype[3] = {GetBasisType(0), GetBasisType(1),
                                              GetBasisType(2)};

    int p, q, r;
    int cnt = 0;

    int BndIdx[3][2];
    int IntIdx[3][2];

    for (i = 0; i < 3; i++)
    {
        BndIdx[i][0] = 0;

        if (Btype[i] == LibUtilities::eModified_A)
        {
            BndIdx[i][1] = 1;
            IntIdx[i][0] = 2;
            IntIdx[i][1] = nummodes[i];
        }
        else
        {
            BndIdx[i][1] = nummodes[i] - 1;
            IntIdx[i][0] = 1;
            IntIdx[i][1] = nummodes[i] - 1;
        }
    }

    for (i = 0; i < 2; i++)
    {
        r = BndIdx[2][i];
        for (q = 0; q < nummodes[1]; q++)
        {
            for (p = 0; p < nummodes[0]; p++)
            {
                outarray[cnt++] =
                    r * nummodes[0] * nummodes[1] + q * nummodes[0] + p;
            }
        }
    }

    for (r = IntIdx[2][0]; r < IntIdx[2][1]; r++)
    {
        for (i = 0; i < 2; i++)
        {
            q = BndIdx[1][i];
            for (p = 0; p < nummodes[0]; p++)
            {
                outarray[cnt++] =
                    r * nummodes[0] * nummodes[1] + q * nummodes[0] + p;
            }
        }

        for (q = IntIdx[1][0]; q < IntIdx[1][1]; q++)
        {
            for (i = 0; i < 2; i++)
            {
                p = BndIdx[0][i];
                outarray[cnt++] =
                    r * nummodes[0] * nummodes[1] + q * nummodes[0] + p;
            }
        }
    }

    sort(outarray.data(), outarray.data() + nBndCoeffs);
}

NekDouble StdHexExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    return BaryTensorDeriv(coord, inarray, firstOrderDerivs);
}

/**
 * Only for basis type Modified_A or GLL_LAGRANGE in all directions.
 */
void StdHexExp::v_GetTraceCoeffMap(const unsigned int fid,
                                   Array<OneD, unsigned int> &maparray)
{
    int i, j;
    int nummodesA = 0, nummodesB = 0;

    ASSERTL1(GetBasisType(0) == GetBasisType(1) &&
                 GetBasisType(0) == GetBasisType(2),
             "Method only implemented if BasisType is indentical in "
             "all directions");
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "Method only implemented for Modified_A or "
             "GLL_Lagrange BasisType");

    const int nummodes0 = m_base[0]->GetNumModes();
    const int nummodes1 = m_base[1]->GetNumModes();
    const int nummodes2 = m_base[2]->GetNumModes();

    switch (fid)
    {
        case 0:
        case 5:
            nummodesA = nummodes0;
            nummodesB = nummodes1;
            break;
        case 1:
        case 3:
            nummodesA = nummodes0;
            nummodesB = nummodes2;
            break;
        case 2:
        case 4:
            nummodesA = nummodes1;
            nummodesB = nummodes2;
            break;
        default:
            ASSERTL0(false, "fid must be between 0 and 5");
    }

    int nFaceCoeffs = nummodesA * nummodesB;

    if (maparray.size() != nFaceCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nFaceCoeffs);
    }

    bool modified = (GetBasisType(0) == LibUtilities::eModified_A);

    int offset = 0;
    int jump1  = 1;
    int jump2  = 1;

    switch (fid)
    {
        case 5:
        {
            if (modified)
            {
                offset = nummodes0 * nummodes1;
            }
            else
            {
                offset = (nummodes2 - 1) * nummodes0 * nummodes1;
                jump1  = nummodes0;
            }
        }
        /* Falls through. */
        case 0:
        {
            jump1 = nummodes0;
            break;
        }
        case 3:
        {
            if (modified)
            {
                offset = nummodes0;
            }
            else
            {
                offset = nummodes0 * (nummodes1 - 1);
                jump1  = nummodes0 * nummodes1;
            }
        }
        /* Falls through. */
        case 1:
        {
            jump1 = nummodes0 * nummodes1;
            break;
        }
        case 2:
        {
            if (modified)
            {
                offset = 1;
            }
            else
            {
                offset = nummodes0 - 1;
                jump1  = nummodes0 * nummodes1;
                jump2  = nummodes0;
            }
        }
        /* Falls through. */
        case 4:
        {
            jump1 = nummodes0 * nummodes1;
            jump2 = nummodes0;
            break;
        }
        default:
            ASSERTL0(false, "fid must be between 0 and 5");
    }

    for (i = 0; i < nummodesB; i++)
    {
        for (j = 0; j < nummodesA; j++)
        {
            maparray[i * nummodesA + j] = i * jump1 + j * jump2 + offset;
        }
    }
}

void StdHexExp::v_GetElmtTraceToTraceMap(const unsigned int fid,
                                         Array<OneD, unsigned int> &maparray,
                                         Array<OneD, int> &signarray,
                                         Orientation faceOrient, int P, int Q)
{
    int i, j;
    int nummodesA = 0, nummodesB = 0;

    ASSERTL1(GetBasisType(0) == GetBasisType(1) &&
                 GetBasisType(0) == GetBasisType(2),
             "Method only implemented if BasisType is indentical in "
             "all directions");
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "Method only implemented for Modified_A or "
             "GLL_Lagrange BasisType");

    const int nummodes0 = m_base[0]->GetNumModes();
    const int nummodes1 = m_base[1]->GetNumModes();
    const int nummodes2 = m_base[2]->GetNumModes();

    switch (fid)
    {
        case 0:
        case 5:
            nummodesA = nummodes0;
            nummodesB = nummodes1;
            break;
        case 1:
        case 3:
            nummodesA = nummodes0;
            nummodesB = nummodes2;
            break;
        case 2:
        case 4:
            nummodesA = nummodes1;
            nummodesB = nummodes2;
            break;
        default:
            ASSERTL0(false, "fid must be between 0 and 5");
    }

    if (P == -1)
    {
        P = nummodesA;
        Q = nummodesB;
    }

    bool modified = (GetBasisType(0) == LibUtilities::eModified_A);

    // check that
    if (modified == false)
    {
        ASSERTL1((P == nummodesA) || (Q == nummodesB),
                 "Different trace space face dimention "
                 "and element face dimention not possible for "
                 "GLL-Lagrange bases");
    }

    int nFaceCoeffs = P * Q;

    if (maparray.size() != nFaceCoeffs)
    {
        maparray = Array<OneD, unsigned int>(nFaceCoeffs);
    }

    // fill default mapping as increasing index
    for (i = 0; i < nFaceCoeffs; ++i)
    {
        maparray[i] = i;
    }

    if (signarray.size() != nFaceCoeffs)
    {
        signarray = Array<OneD, int>(nFaceCoeffs, 1);
    }
    else
    {
        fill(signarray.data(), signarray.data() + nFaceCoeffs, 1);
    }

    // setup indexing to manage transpose directions
    Array<OneD, int> arrayindx(nFaceCoeffs);
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
    // modes are not as large as face models
    for (i = 0; i < nummodesB; i++)
    {
        for (j = nummodesA; j < P; j++)
        {
            signarray[arrayindx[i * P + j]] = 0.0;
            maparray[arrayindx[i * P + j]]  = maparray[0];
        }
    }

    for (i = nummodesB; i < Q; i++)
    {
        for (j = 0; j < P; j++)
        {
            signarray[arrayindx[i * P + j]] = 0.0;
            maparray[arrayindx[i * P + j]]  = maparray[0];
        }
    }

    // zero signmap and set maparray to zero entry if
    // elemental modes are not as large as face modesl
    for (i = 0; i < Q; i++)
    {
        // fill values into map array of trace size
        // for element face index
        for (j = 0; j < P; j++)
        {
            maparray[arrayindx[i * P + j]] = i * nummodesA + j;
        }

        // zero values if P > numModesA
        for (j = nummodesA; j < P; j++)
        {
            signarray[arrayindx[i * P + j]] = 0.0;
            maparray[arrayindx[i * P + j]]  = maparray[0];
        }
    }

    // zero values if Q > numModesB
    for (i = nummodesB; i < Q; i++)
    {
        for (j = 0; j < P; j++)
        {
            signarray[arrayindx[i * P + j]] = 0.0;
            maparray[arrayindx[i * P + j]]  = maparray[0];
        }
    }

    // Now reorientate indices accordign to orientation
    if ((faceOrient == eDir1FwdDir1_Dir2BwdDir2) ||
        (faceOrient == eDir1BwdDir1_Dir2BwdDir2) ||
        (faceOrient == eDir1BwdDir2_Dir2FwdDir1) ||
        (faceOrient == eDir1BwdDir2_Dir2BwdDir1))
    {
        if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
        {
            if (modified)
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
                for (i = 0; i < P; i++)
                {
                    for (j = 0; j < Q / 2; j++)
                    {
                        swap(maparray[i + j * P],
                             maparray[i + P * Q - P - j * P]);
                        swap(signarray[i + j * P],
                             signarray[i + P * Q - P - j * P]);
                    }
                }
            }
        }
        else
        {
            if (modified)
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
            else
            {
                for (i = 0; i < P; i++)
                {
                    for (j = 0; j < Q / 2; j++)
                    {
                        swap(maparray[i * Q + j], maparray[i * Q + Q - 1 - j]);
                        swap(signarray[i * Q + j],
                             signarray[i * Q + Q - 1 - j]);
                    }
                }
            }
        }
    }

    if ((faceOrient == eDir1BwdDir1_Dir2FwdDir2) ||
        (faceOrient == eDir1BwdDir1_Dir2BwdDir2) ||
        (faceOrient == eDir1FwdDir2_Dir2BwdDir1) ||
        (faceOrient == eDir1BwdDir2_Dir2BwdDir1))
    {
        if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
        {
            if (modified)
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
                for (i = 0; i < Q; i++)
                {
                    for (j = 0; j < P / 2; j++)
                    {
                        swap(maparray[i * P + j], maparray[i * P + P - 1 - j]);
                        swap(signarray[i * P + j],
                             signarray[i * P + P - 1 - j]);
                    }
                }
            }
        }
        else
        {
            if (modified)
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
            else
            {
                for (i = 0; i < Q; i++)
                {
                    for (j = 0; j < P / 2; j++)
                    {
                        swap(maparray[i + j * Q],
                             maparray[i + P * Q - Q - j * Q]);
                        swap(signarray[i + j * Q],
                             signarray[i + P * Q - Q - j * Q]);
                    }
                }
            }
        }
    }
}

/**
 * @param   eid         The edge to compute the numbering for.
 * @param   edgeOrient  Orientation of the edge.
 * @param   maparray    Storage for computed mapping array.
 * @param   signarray   ?
 */
void StdHexExp::v_GetEdgeInteriorToElementMap(
    const int eid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation edgeOrient)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    ASSERTL1((eid >= 0) && (eid < 12),
             "local edge id must be between 0 and 11");

    int nEdgeIntCoeffs = GetEdgeNcoeffs(eid) - 2;

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

    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};

    const LibUtilities::BasisType bType[3] = {GetBasisType(0), GetBasisType(1),
                                              GetBasisType(2)};

    bool reverseOrdering = false;
    bool signChange      = false;

    int IdxRange[3][2] = {{0, 0}, {0, 0}, {0, 0}};

    switch (eid)
    {
        case 0:
        case 1:
        case 2:
        case 3:
        {
            IdxRange[2][0] = 0;
            IdxRange[2][1] = 1;
        }
        break;
        case 8:
        case 9:
        case 10:
        case 11:
        {
            if (bType[2] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[2][0] = nummodes[2] - 1;
                IdxRange[2][1] = nummodes[2];
            }
            else
            {
                IdxRange[2][0] = 1;
                IdxRange[2][1] = 2;
            }
        }
        break;
        case 4:
        case 5:
        case 6:
        case 7:
        {
            if (bType[2] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[2][0] = 1;
                IdxRange[2][1] = nummodes[2] - 1;

                if (edgeOrient == eBackwards)
                {
                    reverseOrdering = true;
                }
            }
            else
            {
                IdxRange[2][0] = 2;
                IdxRange[2][1] = nummodes[2];

                if (edgeOrient == eBackwards)
                {
                    signChange = true;
                }
            }
        }
        break;
    }

    switch (eid)
    {
        case 0:
        case 4:
        case 5:
        case 8:
        {
            IdxRange[1][0] = 0;
            IdxRange[1][1] = 1;
        }
        break;
        case 2:
        case 6:
        case 7:
        case 10:
        {
            if (bType[1] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[1][0] = nummodes[1] - 1;
                IdxRange[1][1] = nummodes[1];
            }
            else
            {
                IdxRange[1][0] = 1;
                IdxRange[1][1] = 2;
            }
        }
        break;
        case 1:
        case 9:
        {
            if (bType[1] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[1][0] = 1;
                IdxRange[1][1] = nummodes[1] - 1;

                if (edgeOrient == eBackwards)
                {
                    reverseOrdering = true;
                }
            }
            else
            {
                IdxRange[1][0] = 2;
                IdxRange[1][1] = nummodes[1];

                if (edgeOrient == eBackwards)
                {
                    signChange = true;
                }
            }
        }
        break;
        case 3:
        case 11:
        {
            if (bType[1] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[1][0] = 1;
                IdxRange[1][1] = nummodes[1] - 1;

                if (edgeOrient == eForwards)
                {
                    reverseOrdering = true;
                }
            }
            else
            {
                IdxRange[1][0] = 2;
                IdxRange[1][1] = nummodes[1];

                if (edgeOrient == eForwards)
                {
                    signChange = true;
                }
            }
        }
        break;
    }

    switch (eid)
    {
        case 3:
        case 4:
        case 7:
        case 11:
        {
            IdxRange[0][0] = 0;
            IdxRange[0][1] = 1;
        }
        break;
        case 1:
        case 5:
        case 6:
        case 9:
        {
            if (bType[0] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[0][0] = nummodes[0] - 1;
                IdxRange[0][1] = nummodes[0];
            }
            else
            {
                IdxRange[0][0] = 1;
                IdxRange[0][1] = 2;
            }
        }
        break;
        case 0:
        case 8:
        {
            if (bType[0] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[0][0] = 1;
                IdxRange[0][1] = nummodes[0] - 1;

                if (edgeOrient == eBackwards)
                {
                    reverseOrdering = true;
                }
            }
            else
            {
                IdxRange[0][0] = 2;
                IdxRange[0][1] = nummodes[0];

                if (edgeOrient == eBackwards)
                {
                    signChange = true;
                }
            }
        }
        break;
        case 2:
        case 10:
        {
            if (bType[0] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[0][0] = 1;
                IdxRange[0][1] = nummodes[0] - 1;

                if (edgeOrient == eForwards)
                {
                    reverseOrdering = true;
                }
            }
            else
            {
                IdxRange[0][0] = 2;
                IdxRange[0][1] = nummodes[0];

                if (edgeOrient == eForwards)
                {
                    signChange = true;
                }
            }
        }
        break;
    }

    int cnt = 0;

    for (int r = IdxRange[2][0]; r < IdxRange[2][1]; r++)
    {
        for (int q = IdxRange[1][0]; q < IdxRange[1][1]; q++)
        {
            for (int p = IdxRange[0][0]; p < IdxRange[0][1]; p++)
            {
                maparray[cnt++] =
                    r * nummodes[0] * nummodes[1] + q * nummodes[0] + p;
            }
        }
    }

    if (reverseOrdering)
    {
        reverse(maparray.data(), maparray.data() + nEdgeIntCoeffs);
    }

    if (signChange)
    {
        for (int p = 1; p < nEdgeIntCoeffs; p += 2)
        {
            signarray[p] = -1;
        }
    }
}

/**
 * Generate mapping describing which elemental modes lie on the
 * interior of a given face. Accounts for face orientation.
 */
void StdHexExp::v_GetTraceInteriorToElementMap(
    const int fid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation faceOrient)
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                 GetBasisType(2) == LibUtilities::eGLL_Lagrange,
             "BasisType is not a boundary interior form");

    ASSERTL1((fid >= 0) && (fid < 6), "local face id must be between 0 and 5");

    int nFaceIntCoeffs = v_GetTraceIntNcoeffs(fid);

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

    int nummodes[3] = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes(),
                       m_base[2]->GetNumModes()};

    const LibUtilities::BasisType bType[3] = {GetBasisType(0), GetBasisType(1),
                                              GetBasisType(2)};

    int nummodesA = 0;
    int nummodesB = 0;

    // Determine the number of modes in face directions A & B based
    // on the face index given.
    switch (fid)
    {
        case 0:
        case 5:
        {
            nummodesA = nummodes[0];
            nummodesB = nummodes[1];
        }
        break;
        case 1:
        case 3:
        {
            nummodesA = nummodes[0];
            nummodesB = nummodes[2];
        }
        break;
        case 2:
        case 4:
        {
            nummodesA = nummodes[1];
            nummodesB = nummodes[2];
        }
    }

    Array<OneD, int> arrayindx(nFaceIntCoeffs);

    // Create a mapping array to account for transposition of the
    // coordinates due to face orientation.
    for (int i = 0; i < (nummodesB - 2); i++)
    {
        for (int j = 0; j < (nummodesA - 2); j++)
        {
            if (faceOrient < eDir1FwdDir2_Dir2FwdDir1)
            {
                arrayindx[i * (nummodesA - 2) + j] = i * (nummodesA - 2) + j;
            }
            else
            {
                arrayindx[i * (nummodesA - 2) + j] = j * (nummodesB - 2) + i;
            }
        }
    }

    int IdxRange[3][2];
    int Incr[3];

    Array<OneD, int> sign0(nummodes[0], 1);
    Array<OneD, int> sign1(nummodes[1], 1);
    Array<OneD, int> sign2(nummodes[2], 1);

    // Set the upper and lower bounds, and increment for the faces
    // involving the first coordinate direction.
    switch (fid)
    {
        case 0: // bottom face
        {
            IdxRange[2][0] = 0;
            IdxRange[2][1] = 1;
            Incr[2]        = 1;
        }
        break;
        case 5: // top face
        {
            if (bType[2] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[2][0] = nummodes[2] - 1;
                IdxRange[2][1] = nummodes[2];
                Incr[2]        = 1;
            }
            else
            {
                IdxRange[2][0] = 1;
                IdxRange[2][1] = 2;
                Incr[2]        = 1;
            }
        }
        break;
        default: // all other faces
        {
            if (bType[2] == LibUtilities::eGLL_Lagrange)
            {
                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 2)
                {
                    IdxRange[2][0] = nummodes[2] - 2;
                    IdxRange[2][1] = 0;
                    Incr[2]        = -1;
                }
                else
                {
                    IdxRange[2][0] = 1;
                    IdxRange[2][1] = nummodes[2] - 1;
                    Incr[2]        = 1;
                }
            }
            else
            {
                IdxRange[2][0] = 2;
                IdxRange[2][1] = nummodes[2];
                Incr[2]        = 1;

                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 2)
                {
                    for (int i = 3; i < nummodes[2]; i += 2)
                    {
                        sign2[i] = -1;
                    }
                }
            }
        }
    }

    // Set the upper and lower bounds, and increment for the faces
    // involving the second coordinate direction.
    switch (fid)
    {
        case 1:
        {
            IdxRange[1][0] = 0;
            IdxRange[1][1] = 1;
            Incr[1]        = 1;
        }
        break;
        case 3:
        {
            if (bType[1] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[1][0] = nummodes[1] - 1;
                IdxRange[1][1] = nummodes[1];
                Incr[1]        = 1;
            }
            else
            {
                IdxRange[1][0] = 1;
                IdxRange[1][1] = 2;
                Incr[1]        = 1;
            }
        }
        break;
        case 0:
        case 5:
        {
            if (bType[1] == LibUtilities::eGLL_Lagrange)
            {
                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 2)
                {
                    IdxRange[1][0] = nummodes[1] - 2;
                    IdxRange[1][1] = 0;
                    Incr[1]        = -1;
                }
                else
                {
                    IdxRange[1][0] = 1;
                    IdxRange[1][1] = nummodes[1] - 1;
                    Incr[1]        = 1;
                }
            }
            else
            {
                IdxRange[1][0] = 2;
                IdxRange[1][1] = nummodes[1];
                Incr[1]        = 1;

                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 2)
                {
                    for (int i = 3; i < nummodes[1]; i += 2)
                    {
                        sign1[i] = -1;
                    }
                }
            }
        }
        break;
        default: // case2: case4:
        {
            if (bType[1] == LibUtilities::eGLL_Lagrange)
            {
                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1)
                {
                    IdxRange[1][0] = nummodes[1] - 2;
                    IdxRange[1][1] = 0;
                    Incr[1]        = -1;
                }
                else
                {
                    IdxRange[1][0] = 1;
                    IdxRange[1][1] = nummodes[1] - 1;
                    Incr[1]        = 1;
                }
            }
            else
            {
                IdxRange[1][0] = 2;
                IdxRange[1][1] = nummodes[1];
                Incr[1]        = 1;

                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1)
                {
                    for (int i = 3; i < nummodes[1]; i += 2)
                    {
                        sign1[i] = -1;
                    }
                }
            }
        }
    }

    switch (fid)
    {
        case 4:
        {
            IdxRange[0][0] = 0;
            IdxRange[0][1] = 1;
            Incr[0]        = 1;
        }
        break;
        case 2:
        {
            if (bType[0] == LibUtilities::eGLL_Lagrange)
            {
                IdxRange[0][0] = nummodes[0] - 1;
                IdxRange[0][1] = nummodes[0];
                Incr[0]        = 1;
            }
            else
            {
                IdxRange[0][0] = 1;
                IdxRange[0][1] = 2;
                Incr[0]        = 1;
            }
        }
        break;
        default:
        {
            if (bType[0] == LibUtilities::eGLL_Lagrange)
            {
                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1)
                {
                    IdxRange[0][0] = nummodes[0] - 2;
                    IdxRange[0][1] = 0;
                    Incr[0]        = -1;
                }
                else
                {
                    IdxRange[0][0] = 1;
                    IdxRange[0][1] = nummodes[0] - 1;
                    Incr[0]        = 1;
                }
            }
            else
            {
                IdxRange[0][0] = 2;
                IdxRange[0][1] = nummodes[0];
                Incr[0]        = 1;

                if (((int)(faceOrient - eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1)
                {
                    for (int i = 3; i < nummodes[0]; i += 2)
                    {
                        sign0[i] = -1;
                    }
                }
            }
        }
    }

    int cnt = 0;

    for (int r = IdxRange[2][0]; r != IdxRange[2][1]; r += Incr[2])
    {
        for (int q = IdxRange[1][0]; q != IdxRange[1][1]; q += Incr[1])
        {
            for (int p = IdxRange[0][0]; p != IdxRange[0][1]; p += Incr[0])
            {
                maparray[arrayindx[cnt]] =
                    r * nummodes[0] * nummodes[1] + q * nummodes[0] + p;
                signarray[arrayindx[cnt++]] = sign0[p] * sign1[q] * sign2[r];
            }
        }
    }
}

int StdHexExp::v_GetEdgeNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 11), "edge id is out of range");

    if ((i == 0) || (i == 2) || (i == 8) || (i == 10))
    {
        return GetBasisNumModes(0);
    }
    else if ((i == 1) || (i == 3) || (i == 9) || (i == 11))
    {
        return GetBasisNumModes(1);
    }
    else
    {
        return GetBasisNumModes(2);
    }
}

DNekMatSharedPtr StdHexExp::v_GenMatrix(const StdMatrixKey &mkey)
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
                LibUtilities::StdHexData::getNumberOfCoefficients(nq, nq, nq);
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
                    for (int k = 0; k < nq; ++k, ++cnt)
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
        default:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
        }
        break;
    }

    return Mat;
}

DNekMatSharedPtr StdHexExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

void StdHexExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray,
                               const StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdHexExp::v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdHexExp::v_LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdHexExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void StdHexExp::v_WeakDerivMatrixOp(const int i,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void StdHexExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdHexExp::v_HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdHexExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
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
    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_A, nmodes_c, pc);
    StdHexExp OrthoExp(Ba, Bb, Bc);

    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
    int cnt = 0;

    // project onto modal  space.
    OrthoExp.FwdTrans(array, orthocoeffs);

    if (mkey.ConstFactorExists(eFactorSVVPowerKerDiffCoeff))
    {
        // Rodrigo's power kernel
        NekDouble cutoff = mkey.GetConstFactor(eFactorSVVCutoffRatio);
        NekDouble SvvDiffCoeff =
            mkey.GetConstFactor(eFactorSVVPowerKerDiffCoeff) *
            mkey.GetConstFactor(eFactorSVVDiffCoeff);

        for (int i = 0; i < nmodes_a; ++i)
        {
            for (int j = 0; j < nmodes_b; ++j)
            {
                NekDouble fac1 = std::max(
                    pow((1.0 * i) / (nmodes_a - 1), cutoff * nmodes_a),
                    pow((1.0 * j) / (nmodes_b - 1), cutoff * nmodes_b));

                for (int k = 0; k < nmodes_c; ++k)
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

        for (int i = 0; i < nmodes_a; ++i)
        {
            for (int j = 0; j < nmodes_b; ++j)
            {
                int maxij = max(i, j);

                for (int k = 0; k < nmodes_c; ++k)
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

        int cutoff = (int)(mkey.GetConstFactor(eFactorSVVCutoffRatio) *
                           min(nmodes_a, nmodes_b));
        NekDouble SvvDiffCoeff = mkey.GetConstFactor(eFactorSVVDiffCoeff);
        //  Filter just trilinear space
        int nmodes = max(nmodes_a, nmodes_b);
        nmodes     = max(nmodes, nmodes_c);

        Array<OneD, NekDouble> fac(nmodes, 1.0);
        for (int j = cutoff; j < nmodes; ++j)
        {
            fac[j] = fabs((j - nmodes) / ((NekDouble)(j - cutoff + 1.0)));
            fac[j] *= fac[j]; // added this line to conform with equation
        }

        for (int i = 0; i < nmodes_a; ++i)
        {
            for (int j = 0; j < nmodes_b; ++j)
            {
                for (int k = 0; k < nmodes_c; ++k)
                {
                    if ((i >= cutoff) || (j >= cutoff) || (k >= cutoff))
                    {
                        orthocoeffs[i * nmodes_a * nmodes_b + j * nmodes_c +
                                    k] *=
                            (SvvDiffCoeff * exp(-(fac[i] + fac[j] + fac[k])));
                    }
                    else
                    {
                        orthocoeffs[i * nmodes_a * nmodes_b + j * nmodes_c +
                                    k] *= 0.0;
                    }
                }
            }
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdHexExp::v_ExponentialFilter(Array<OneD, NekDouble> &array,
                                    const NekDouble alpha,
                                    const NekDouble exponent,
                                    const NekDouble cutoff)
{
    // Generate an orthogonal expansion
    int qa      = m_base[0]->GetNumPoints();
    int qb      = m_base[1]->GetNumPoints();
    int qc      = m_base[2]->GetNumPoints();
    int nmodesA = m_base[0]->GetNumModes();
    int nmodesB = m_base[1]->GetNumModes();
    int nmodesC = m_base[2]->GetNumModes();
    int P       = nmodesA - 1;
    int Q       = nmodesB - 1;
    int R       = nmodesC - 1;

    // Declare orthogonal basis.
    LibUtilities::PointsKey pa(qa, m_base[0]->GetPointsType());
    LibUtilities::PointsKey pb(qb, m_base[1]->GetPointsType());
    LibUtilities::PointsKey pc(qc, m_base[2]->GetPointsType());

    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, nmodesA, pa);
    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, nmodesB, pb);
    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_A, nmodesC, pc);
    StdHexExp OrthoExp(Ba, Bb, Bc);

    // Cutoff
    int Pcut = cutoff * P;
    int Qcut = cutoff * Q;
    int Rcut = cutoff * R;

    // Project onto orthogonal space.
    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
    OrthoExp.FwdTrans(array, orthocoeffs);

    //
    NekDouble fac, fac1, fac2, fac3;
    int index = 0;
    for (int i = 0; i < nmodesA; ++i)
    {
        for (int j = 0; j < nmodesB; ++j)
        {
            for (int k = 0; k < nmodesC; ++k, ++index)
            {
                // to filter out only the "high-modes"
                if (i > Pcut || j > Qcut || k > Rcut)
                {
                    fac1 = (NekDouble)(i - Pcut) / ((NekDouble)(P - Pcut));
                    fac2 = (NekDouble)(j - Qcut) / ((NekDouble)(Q - Qcut));
                    fac3 = (NekDouble)(k - Rcut) / ((NekDouble)(R - Rcut));
                    fac  = max(max(fac1, fac2), fac3);
                    fac  = pow(fac, exponent);
                    orthocoeffs[index] *= exp(-alpha * fac);
                }
            }
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdHexExp::v_GetSimplexEquiSpacedConnectivity(
    Array<OneD, int> &conn, [[maybe_unused]] bool standard)
{
    int np0 = m_base[0]->GetNumPoints();
    int np1 = m_base[1]->GetNumPoints();
    int np2 = m_base[2]->GetNumPoints();
    int np  = max(np0, max(np1, np2));

    conn = Array<OneD, int>(6 * (np - 1) * (np - 1) * (np - 1));

    int row   = 0;
    int rowp1 = 0;
    int cnt   = 0;
    int plane = 0;
    for (int i = 0; i < np - 1; ++i)
    {
        for (int j = 0; j < np - 1; ++j)
        {
            rowp1 += np;
            for (int k = 0; k < np - 1; ++k)
            {
                conn[cnt++] = plane + row + k;
                conn[cnt++] = plane + row + k + 1;
                conn[cnt++] = plane + rowp1 + k;

                conn[cnt++] = plane + rowp1 + k + 1;
                conn[cnt++] = plane + rowp1 + k;
                conn[cnt++] = plane + row + k + 1;
            }
            row += np;
        }
        plane += np * np;
    }
}
} // namespace Nektar::StdRegions
