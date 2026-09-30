///////////////////////////////////////////////////////////////////////////////
//
// File: StdQuadExp.cpp
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
// Description: Quadrilateral routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <StdRegions/StdQuadExp.h>

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

/** \brief Constructor using BasisKey class for quadrature
 *  points and order definition
 */
StdQuadExp::StdQuadExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb)
    : StdExpansion(Ba.GetNumModes() * Bb.GetNumModes(), 2, Ba, Bb),
      StdExpansion2D(Ba.GetNumModes() * Bb.GetNumModes(), Ba, Bb)
{
    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());

    // cache integration weights for future use
    m_weights.push_back(m_base[1]->GetW());
}

/////////////////////////////
// Differentiation Methods //
/////////////////////////////
/** \brief Calculate the derivative of the physical points
 *
 *  For quadrilateral region can use the Tensor_Deriv function
 *  defined under StdExpansion.
 */

void StdQuadExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &out_d0,
                                Array<OneD, NekDouble> &out_d1,
                                [[maybe_unused]] Array<OneD, NekDouble> &out_d2)
{
    PhysTensorDeriv(inarray, out_d0, out_d1);
}

////////////////
// Transforms //
////////////////

void StdQuadExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();

    if (m_base[0]->Collocation() && m_base[1]->Collocation())
    {
        std::memcpy(outarray.data(), inarray.data(),
                    nquad0 * nquad1 * sizeof(NekDouble));
    }
    else
    {
        const Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
        const Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();

        int nmodes0 = m_base[0]->GetNumModes();
        int nmodes1 = m_base[1]->GetNumModes();

        std::vector<vec_t, tinysimd::allocator<vec_t>> wsp0(nmodes1 * nquad0);

        // Switch statment using boost_pp and macros. This unfolls intwo a
        // nested swtich statement where the outer swtich statement runs
        // from SMIN to SMAX for modal order and the inner switch
        // statemets run from the outer value of the case to 2*SMAX for
        // the quadrature order. If you want to see it unwrapped compile
        // in verbose mode and add --preprocess to the c++ command.
        // Default case
#undef BWDTRANS_DEF
#define BWDTRANS_DEF                                                           \
    BwdTransQuadKernel(                                                        \
        nmodes0, nmodes1, nquad0, nquad1, (const vec_t *)base0.data(),         \
        (const vec_t *)base1.data(), wsp0.data(),                              \
        (const vec_t *)inarray.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef BWDTRANS_Q
#define BWDTRANS_Q(r, i)                                                       \
    case NQ(i):                                                                \
        BwdTransQuadKernel(                                                    \
            NM(i), NM(i), NQ(i), NQ(i), (const vec_t *)base0.data(),           \
            (const vec_t *)base1.data(), wsp0.data(),                          \
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
        if ((nmodes0 == nmodes1) && (nquad0 == nquad1))
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

void StdQuadExp::v_FwdTransBndConstrained(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    if ((m_base[0]->Collocation()) && (m_base[1]->Collocation()))
    {
        Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
    }
    else
    {
        int i, j;
        int npoints[2] = {m_base[0]->GetNumPoints(), m_base[1]->GetNumPoints()};
        int nmodes[2]  = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes()};

        fill(outarray.data(), outarray.data() + m_ncoeffs, 0.0);

        Array<OneD, NekDouble> physEdge[4];
        Array<OneD, NekDouble> coeffEdge[4];
        for (i = 0; i < 4; i++)
        {
            physEdge[i]  = Array<OneD, NekDouble>(npoints[i % 2]);
            coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i % 2]);
        }

        for (i = 0; i < npoints[0]; i++)
        {
            physEdge[0][i] = inarray[i];
            physEdge[2][i] = inarray[npoints[0] * npoints[1] - 1 - i];
        }

        for (i = 0; i < npoints[1]; i++)
        {
            physEdge[1][i] = inarray[npoints[0] - 1 + i * npoints[0]];
            physEdge[3][i] =
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

        for (i = 0; i < 4; i++)
        {
            segexp[i % 2]->FwdTransBndConstrained(physEdge[i], coeffEdge[i]);

            GetTraceToElementMap(i, mapArray, signArray);
            for (j = 0; j < nmodes[i % 2]; j++)
            {
                sign                  = (NekDouble)signArray[j];
                outarray[mapArray[j]] = sign * coeffEdge[i][j];
            }
        }

        Array<OneD, NekDouble> tmp0(m_ncoeffs);
        Array<OneD, NekDouble> tmp1(m_ncoeffs);

        StdMatrixKey masskey(eMass, DetShapeType(), *this);
        MassMatrixOp(outarray, tmp0, masskey);
        IProductWRTBase(inarray, tmp1);

        Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

        // get Mass matrix inverse (only of interior DOF)
        // use block (1,1) of the static condensed system
        // note: this block alreay contains the inverse matrix
        DNekMatSharedPtr matsys =
            (m_stdStaticCondMatrixManager[masskey])->GetBlock(1, 1);

        int nBoundaryDofs = NumBndryCoeffs();
        int nInteriorDofs = m_ncoeffs - nBoundaryDofs;

        Array<OneD, NekDouble> rhs(nInteriorDofs);
        Array<OneD, NekDouble> result(nInteriorDofs);

        GetInteriorMap(mapArray);

        for (i = 0; i < nInteriorDofs; i++)
        {
            rhs[i] = tmp1[mapArray[i]];
        }

        Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, 1.0,
                    &(matsys->GetPtr())[0], nInteriorDofs, rhs.data(), 1, 0.0,
                    result.data(), 1);

        for (i = 0; i < nInteriorDofs; i++)
        {
            outarray[mapArray[i]] = result[i];
        }
    }
}

/////////////////////////////
// Inner Product Functions //
/////////////////////////////
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
 *  @param CollDir0 - bool to identify if 0-direction basis is a
 *  collocated expansion
 *  @param CollDir1 - bool to identify if 1-direction basis is a
 *  collocated expansion
 */
void StdQuadExp::v_IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &base1,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, NekDouble> &jac,
    const bool Deformed, const bool CollDir0, const bool CollDir1)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int order0 = m_base[0]->GetNumModes();
    int order1 = m_base[1]->GetNumModes();

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
    IProductQuadKernel<false, false, true>(                                    \
        order0, order1, nquad0, nquad1, (const vec_t *)inarray.data(),         \
        (const vec_t *)base0.data(), (const vec_t *)base1.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)outarray.data(), 1.0, CollDir0,         \
        CollDir1)

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductQuadKernel<false, false, true>(                                \
            NM(i), NM(i), NQ(i), NQ(i), (const vec_t *)inarray.data(),         \
            (const vec_t *)base0.data(), (const vec_t *)base1.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)outarray.data(), 1.0, CollDir0,     \
            CollDir1);                                                         \
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
        if ((order0 == order1) && (nquad0 == nquad1))
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
    IProductQuadKernel<false, false, false>(                                   \
        order0, order1, nquad0, nquad1, (const vec_t *)inarray.data(),         \
        (const vec_t *)base0.data(), (const vec_t *)base1.data(),              \
        (const vec_t *)m_weights[0].data(),                                    \
        (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),         \
        (vec_t *)wsp0.data(), (vec_t *)outarray.data(), 1.0, CollDir0,         \
        CollDir1)

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductQuadKernel<false, false, false>(                               \
            NM(i), NM(i), NQ(i), NQ(i), (const vec_t *)inarray.data(),         \
            (const vec_t *)base0.data(), (const vec_t *)base1.data(),          \
            (const vec_t *)m_weights[0].data(),                                \
            (const vec_t *)m_weights[1].data(), (const vec_t *)jac.data(),     \
            (vec_t *)wsp0.data(), (vec_t *)outarray.data(), 1.0, CollDir0,     \
            CollDir1);                                                         \
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
        if ((order0 == order1) && (nquad0 == nquad1))
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

void StdQuadExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    ASSERTL0((dir == 0) || (dir == 1), "input dir is out of range");

    Array<OneD, NekDouble> one(1, 1.0);
    if (dir) // dir == 1
    {
        v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                                inarray, outarray, one, false,
                                m_base[0]->Collocation(), false);
    }
    else // dir == 0
    {
        v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                                inarray, outarray, one, false, false,
                                m_base[1]->Collocation());
    }
}

//////////////////////////
// Evaluation functions //
//////////////////////////

void StdQuadExp::v_LocCoordToLocCollapsed(
    const Array<OneD, const NekDouble> &xi, Array<OneD, NekDouble> &eta)
{
    eta[0] = xi[0];
    eta[1] = xi[1];
}

void StdQuadExp::v_LocCollapsedToLocCoord(
    const Array<OneD, const NekDouble> &eta, Array<OneD, NekDouble> &xi)
{
    xi[0] = eta[0];
    xi[1] = eta[1];
}

/** \brief Fill outarray with mode \a mode of expansion
 *
 *  Note for quadrilateral expansions _base[0] (i.e. p)  modes run
 *  fastest
 */

void StdQuadExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
{
    int i;
    int nquad0                         = m_base[0]->GetNumPoints();
    int nquad1                         = m_base[1]->GetNumPoints();
    Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
    Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();
    int btmp0                          = m_base[0]->GetNumModes();
    int mode0                          = mode % btmp0;
    int mode1                          = mode / btmp0;

    ASSERTL2(mode1 == (int)floor((1.0 * mode) / btmp0),
             "Integer Truncation not Equiv to Floor");

    ASSERTL2(m_ncoeffs > mode,
             "calling argument mode is larger than total expansion order");

    for (i = 0; i < nquad1; ++i)
    {
        Vmath::Vcopy(nquad0, (NekDouble *)(base0.data() + mode0 * nquad0), 1,
                     &outarray[0] + i * nquad0, 1);
    }

    for (i = 0; i < nquad0; ++i)
    {
        Vmath::Vmul(nquad1, (NekDouble *)(base1.data() + mode1 * nquad1), 1,
                    &outarray[0] + i, nquad0, &outarray[0] + i, nquad0);
    }
}

//////////////////////
// Helper functions //
//////////////////////

int StdQuadExp::v_GetNverts() const
{
    return 4;
}

int StdQuadExp::v_GetNtraces() const
{
    return 4;
}

int StdQuadExp::v_GetTraceNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 3), "edge id is out of range");

    if ((i == 0) || (i == 2))
    {
        return GetBasisNumModes(0);
    }
    else
    {
        return GetBasisNumModes(1);
    }
}

int StdQuadExp::v_GetTraceIntNcoeffs(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 4), "edge id is out of range");
    if ((i == 0) || (i == 2))
    {
        return GetBasisNumModes(0) - 2;
    }
    else
    {
        return GetBasisNumModes(1) - 2;
    }
}

int StdQuadExp::v_GetTraceNumPoints(const int i) const
{
    ASSERTL2((i >= 0) && (i <= 3), "edge id is out of range");

    if ((i == 0) || (i == 2))
    {
        return GetNumPoints(0);
    }
    else
    {
        return GetNumPoints(1);
    }
}

const LibUtilities::BasisKey StdQuadExp::v_GetTraceBasisKey(
    const int i, [[maybe_unused]] const int j,
    [[maybe_unused]] bool UseGLL) const
{
    ASSERTL2((i >= 0) && (i <= 3), "edge id is out of range");

    if ((i == 0) || (i == 2))
    {
        switch (GetBasis(0)->GetBasisType())
        {
            case LibUtilities::eOrtho_A:
            {
                return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,
                                              GetBasis(0)->GetNumModes(),
                                              GetBasis(0)->GetPointsKey());
            }
            break;
            default:
            {
                return GetBasis(0)->GetBasisKey();
            }
        }
    }
    else
    {
        switch (GetBasis(1)->GetBasisType())
        {
            case LibUtilities::eOrtho_A:
            {
                return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,
                                              GetBasis(1)->GetNumModes(),
                                              GetBasis(1)->GetPointsKey());
            }
            break;
            default:
            {
                return GetBasis(1)->GetBasisKey();
            }
        }
    }
}

LibUtilities::ShapeType StdQuadExp::v_DetShapeType() const
{
    return LibUtilities::eQuadrilateral;
}

int StdQuadExp::v_NumBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange ||
                 GetBasisType(0) == LibUtilities::eGauss_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange ||
                 GetBasisType(0) == LibUtilities::eGauss_Lagrange,
             "BasisType is not a boundary interior form");

    return 4 + 2 * (GetBasisNumModes(0) - 2) + 2 * (GetBasisNumModes(1) - 2);
}

int StdQuadExp::v_NumDGBndryCoeffs() const
{
    ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                 GetBasisType(0) == LibUtilities::eGLL_Lagrange ||
                 GetBasisType(0) == LibUtilities::eGauss_Lagrange,
             "BasisType is not a boundary interior form");
    ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                 GetBasisType(1) == LibUtilities::eGLL_Lagrange ||
                 GetBasisType(0) == LibUtilities::eGauss_Lagrange,
             "BasisType is not a boundary interior form");

    return 2 * GetBasisNumModes(0) + 2 * GetBasisNumModes(1);
}

int StdQuadExp::v_CalcNumberOfCoefficients(
    const std::vector<unsigned int> &nummodes, int &modes_offset)
{
    int nmodes = nummodes[modes_offset] * nummodes[modes_offset + 1];
    modes_offset += 2;

    return nmodes;
}

bool StdQuadExp::v_IsBoundaryInteriorExpansion() const
{
    bool returnval = false;

    if ((m_base[0]->GetBasisType() == LibUtilities::eModified_A) ||
        (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange))
    {
        if ((m_base[1]->GetBasisType() == LibUtilities::eModified_A) ||
            (m_base[1]->GetBasisType() == LibUtilities::eGLL_Lagrange))
        {
            returnval = true;
        }
    }

    return returnval;
}

void StdQuadExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                             Array<OneD, NekDouble> &coords_1,
                             [[maybe_unused]] Array<OneD, NekDouble> &coords_2)
{
    Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
    Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
    int nq0                         = GetNumPoints(0);
    int nq1                         = GetNumPoints(1);
    int i;

    for (i = 0; i < nq1; ++i)
    {
        Vmath::Vcopy(nq0, z0.data(), 1, &coords_0[0] + i * nq0, 1);
        Vmath::Fill(nq0, z1[i], &coords_1[0] + i * nq0, 1);
    }
}

/**
 * @brief This function evaluates the basis function mode @p mode at a
 * point @p coords of the domain.
 *
 * This function uses barycentric interpolation with the tensor
 * product separation of the basis function to improve performance.
 *
 * @param coord   The coordinate inside the standard region.
 * @param mode    The mode number to be evaluated.
 *
 * @return The value of the basis function @p mode at @p coords.
 */
NekDouble StdQuadExp::v_PhysEvaluateBasis(
    const Array<OneD, const NekDouble> &coords, int mode)
{
    ASSERTL2(coords[0] > -1 - NekConstants::kNekZeroTol, "coord[0] < -1");
    ASSERTL2(coords[0] < 1 + NekConstants::kNekZeroTol, "coord[0] >  1");
    ASSERTL2(coords[1] > -1 - NekConstants::kNekZeroTol, "coord[1] < -1");
    ASSERTL2(coords[1] < 1 + NekConstants::kNekZeroTol, "coord[1] >  1");

    const int nm0 = m_base[0]->GetNumModes();
    const int nm1 = m_base[1]->GetNumModes();

    return StdExpansion::BaryEvaluateBasis<0>(coords[0], mode % nm1) *
           StdExpansion::BaryEvaluateBasis<1>(coords[1], mode / nm0);
}

NekDouble StdQuadExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    return BaryTensorDeriv(coord, inarray, firstOrderDerivs);
}

//////////////
// Mappings //
//////////////

void StdQuadExp::v_GetBoundaryMap(Array<OneD, unsigned int> &outarray)
{
    int i;
    int cnt = 0;
    int nummodes0, nummodes1;
    int value1 = 0, value2 = 0;
    if (outarray.size() != NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
    }

    nummodes0 = m_base[0]->GetNumModes();
    nummodes1 = m_base[1]->GetNumModes();

    const LibUtilities::BasisType Btype0 = GetBasisType(0);
    const LibUtilities::BasisType Btype1 = GetBasisType(1);

    switch (Btype1)
    {
        case LibUtilities::eGLL_Lagrange:
        case LibUtilities::eGauss_Lagrange:
            value1 = nummodes0;
            break;
        case LibUtilities::eModified_A:
            value1 = 2 * nummodes0;
            break;
        default:
            ASSERTL0(0, "Mapping array is not defined for this expansion");
            break;
    }

    for (i = 0; i < value1; i++)
    {
        outarray[i] = i;
    }
    cnt = value1;

    switch (Btype0)
    {
        case LibUtilities::eGLL_Lagrange:
        case LibUtilities::eGauss_Lagrange:
            value2 = value1 + nummodes0 - 1;
            break;
        case LibUtilities::eModified_A:
            value2 = value1 + 1;
            break;
        default:
            ASSERTL0(0, "Mapping array is not defined for this expansion");
            break;
    }

    for (i = 0; i < nummodes1 - 2; i++)
    {
        outarray[cnt++] = value1 + i * nummodes0;
        outarray[cnt++] = value2 + i * nummodes0;
    }

    if (Btype1 == LibUtilities::eGLL_Lagrange ||
        Btype1 == LibUtilities::eGauss_Lagrange)
    {
        for (i = nummodes0 * (nummodes1 - 1); i < GetNcoeffs(); i++)
        {
            outarray[cnt++] = i;
        }
    }
}

void StdQuadExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    int i, j;
    int cnt = 0;
    int nummodes0, nummodes1;
    int startvalue = 0;
    if (outarray.size() != GetNcoeffs() - NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(GetNcoeffs() - NumBndryCoeffs());
    }

    nummodes0 = m_base[0]->GetNumModes();
    nummodes1 = m_base[1]->GetNumModes();

    const LibUtilities::BasisType Btype0 = GetBasisType(0);
    const LibUtilities::BasisType Btype1 = GetBasisType(1);

    switch (Btype1)
    {
        case LibUtilities::eGLL_Lagrange:
            startvalue = nummodes0;
            break;
        case LibUtilities::eModified_A:
            startvalue = 2 * nummodes0;
            break;
        default:
            ASSERTL0(0, "Mapping array is not defined for this expansion");
            break;
    }

    switch (Btype0)
    {
        case LibUtilities::eGLL_Lagrange:
            startvalue++;
            break;
        case LibUtilities::eModified_A:
            startvalue += 2;
            break;
        default:
            ASSERTL0(0, "Mapping array is not defined for this expansion");
            break;
    }

    for (i = 0; i < nummodes1 - 2; i++)
    {
        for (j = 0; j < nummodes0 - 2; j++)
        {
            outarray[cnt++] = startvalue + j;
        }
        startvalue += nummodes0;
    }
}

int StdQuadExp::v_GetVertexMap(int localVertexId, bool useCoeffPacking)
{
    int localDOF = 0;

    if (useCoeffPacking == true)
    {
        switch (localVertexId)
        {
            case 0:
            {
                localDOF = 0;
            }
            break;
            case 1:
            {
                if (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    localDOF = m_base[0]->GetNumModes() - 1;
                }
                else
                {
                    localDOF = 1;
                }
            }
            break;
            case 2:
            {
                if (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    localDOF = m_base[0]->GetNumModes() *
                               (m_base[1]->GetNumModes() - 1);
                }
                else
                {
                    localDOF = m_base[0]->GetNumModes();
                }
            }
            break;
            case 3:
            {
                if (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    localDOF =
                        m_base[0]->GetNumModes() * m_base[1]->GetNumModes() - 1;
                }
                else
                {
                    localDOF = m_base[0]->GetNumModes() + 1;
                }
            }
            break;
            default:
                ASSERTL0(false, "eid must be between 0 and 3");
                break;
        }
    }
    else
    {
        switch (localVertexId)
        {
            case 0:
            {
                localDOF = 0;
            }
            break;
            case 1:
            {
                if (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    localDOF = m_base[0]->GetNumModes() - 1;
                }
                else
                {
                    localDOF = 1;
                }
            }
            break;
            case 2:
            {
                if (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    localDOF =
                        m_base[0]->GetNumModes() * m_base[1]->GetNumModes() - 1;
                }
                else
                {
                    localDOF = m_base[0]->GetNumModes() + 1;
                }
            }
            break;
            case 3:
            {
                if (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
                {
                    localDOF = m_base[0]->GetNumModes() *
                               (m_base[1]->GetNumModes() - 1);
                }
                else
                {
                    localDOF = m_base[0]->GetNumModes();
                }
            }
            break;
            default:
                ASSERTL0(false, "eid must be between 0 and 3");
                break;
        }
    }
    return localDOF;
}

/**  \brief Get the map of the coefficient location to teh
 *    local trace coefficients
 */

void StdQuadExp::v_GetTraceCoeffMap(const unsigned int traceid,
                                    Array<OneD, unsigned int> &maparray)
{
    ASSERTL1(traceid < 4, "traceid must be between 0 and 3");

    unsigned int i;
    unsigned int order0   = m_base[0]->GetNumModes();
    unsigned int order1   = m_base[1]->GetNumModes();
    unsigned int numModes = (traceid % 2) ? order1 : order0;

    if (maparray.size() != numModes)
    {
        maparray = Array<OneD, unsigned int>(numModes);
    }

    const LibUtilities::BasisType bType = GetBasisType(traceid % 2);

    if (bType == LibUtilities::eModified_A)
    {
        switch (traceid)
        {
            case 0:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = i;
                }
            }
            break;
            case 1:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = i * order0 + 1;
                }
            }
            break;
            case 2:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = order0 + i;
                }
            }
            break;
            case 3:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = i * order0;
                }
            }
            break;
            default:
                break;
        }
    }
    else if (bType == LibUtilities::eGLL_Lagrange ||
             bType == LibUtilities::eGauss_Lagrange)
    {
        switch (traceid)
        {
            case 0:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = i;
                }
            }
            break;
            case 1:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = (i + 1) * order0 - 1;
                }
            }
            break;
            case 2:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = order0 * (order1 - 1) + i;
                }
            }
            break;
            case 3:
            {
                for (i = 0; i < numModes; i++)
                {
                    maparray[i] = order0 * i;
                }
            }
            break;
            default:
                break;
        }
    }
    else
    {
        ASSERTL0(false, "Mapping not defined for this type of basis");
    }
}

void StdQuadExp::v_GetTraceInteriorToElementMap(
    const int eid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, const Orientation edgeOrient)
{
    int i;
    const int nummodes0                 = m_base[0]->GetNumModes();
    const int nummodes1                 = m_base[1]->GetNumModes();
    const int nEdgeIntCoeffs            = GetTraceNcoeffs(eid) - 2;
    const LibUtilities::BasisType bType = GetBasisType(eid % 2);

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

    if (bType == LibUtilities::eModified_A)
    {
        switch (eid)
        {
            case 0:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = i + 2;
                }
            }
            break;
            case 1:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = (i + 2) * nummodes0 + 1;
                }
            }
            break;
            case 2:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = nummodes0 + i + 2;
                }
            }
            break;
            case 3:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = (i + 2) * nummodes0;
                }
            }
            break;
            default:
                ASSERTL0(false, "eid must be between 0 and 3");
                break;
        }

        if (edgeOrient == eBackwards)
        {
            for (i = 1; i < nEdgeIntCoeffs; i += 2)
            {
                signarray[i] = -1;
            }
        }
    }
    else if (bType == LibUtilities::eGLL_Lagrange)
    {
        switch (eid)
        {
            case 0:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = i + 1;
                }
            }
            break;
            case 1:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = (i + 2) * nummodes0 - 1;
                }
            }
            break;
            case 2:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = nummodes0 * (nummodes1 - 1) + i + 1;
                }
            }
            break;
            case 3:
            {
                for (i = 0; i < nEdgeIntCoeffs; i++)
                {
                    maparray[i] = nummodes0 * (i + 1);
                }
            }
            break;
            default:
                ASSERTL0(false, "eid must be between 0 and 3");
                break;
        }
        if (edgeOrient == eBackwards)
        {
            reverse(maparray.data(), maparray.data() + nEdgeIntCoeffs);
        }
    }
    else
    {
        ASSERTL0(false, "Mapping not defined for this type of basis");
    }
}

///////////////////////
// Wrapper Functions //
///////////////////////

DNekMatSharedPtr StdQuadExp::v_GenMatrix(const StdMatrixKey &mkey)
{
    int i, j;
    int order0       = GetBasisNumModes(0);
    int order1       = GetBasisNumModes(1);
    MatrixType mtype = mkey.GetMatrixType();

    DNekMatSharedPtr Mat;

    switch (mtype)
    {
        case ePhysInterpToEquiSpaced:
        {
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            int nq;

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = max(nq0, nq1);
            }

            int neq =
                LibUtilities::StdQuadData::getNumberOfCoefficients(nq, nq);
            Array<OneD, Array<OneD, NekDouble>> coords(neq);
            Array<OneD, NekDouble> coll(2);
            Array<OneD, DNekMatSharedPtr> I(2);
            Array<OneD, NekDouble> tmp(nq0);

            Mat     = MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1);
            int cnt = 0;

            for (i = 0; i < nq; ++i)
            {
                for (j = 0; j < nq; ++j, ++cnt)
                {
                    coords[cnt]    = Array<OneD, NekDouble>(2);
                    coords[cnt][0] = -1.0 + 2 * j / (NekDouble)(nq - 1);
                    coords[cnt][1] = -1.0 + 2 * i / (NekDouble)(nq - 1);
                }
            }

            for (i = 0; i < neq; ++i)
            {
                LocCoordToLocCollapsed(coords[i], coll);

                I[0] = m_base[0]->GetI(coll);
                I[1] = m_base[1]->GetI(coll + 1);

                // interpolate first coordinate direction
                for (j = 0; j < nq1; ++j)
                {
                    NekDouble fac = (I[1]->GetPtr())[j];
                    Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                    Vmath::Vcopy(nq0, &tmp[0], 1,
                                 Mat->GetRawPtr() + j * nq0 * neq + i, neq);
                }
            }
            break;
        }
        case ePhysInterpToGLL:
        {
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            int nq;

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = max(nq0, nq1);
            }

            int neq =
                LibUtilities::StdQuadData::getNumberOfCoefficients(nq, nq);
            Array<OneD, Array<OneD, NekDouble>> coords(neq);
            Array<OneD, NekDouble> coll(2);
            Array<OneD, DNekMatSharedPtr> I(2);
            Array<OneD, NekDouble> tmp(nq0);

            Mat     = MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0 * nq1);
            int cnt = 0;

            const LibUtilities::PointsKey key(
                nq, LibUtilities::eGaussLobattoLegendre);

            Array<OneD, const NekDouble> z;
            LibUtilities::PointsManager()[key]->GetPoints(z);

            for (int i = 0; i < nq; ++i)
            {
                for (int j = 0; j < nq; ++j, ++cnt)
                {
                    coords[cnt]    = Array<OneD, NekDouble>(2);
                    coords[cnt][0] = z[j];
                    coords[cnt][1] = z[i];
                }
            }

            for (int i = 0; i < neq; ++i)
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

            StdQuadExp Exp2D(ba, bb);
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
        case eMass:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
            // For Fourier basis set the imaginary component of mean mode
            // to have a unit diagonal component in mass matrix
            if (m_base[0]->GetBasisType() == LibUtilities::eFourier)
            {
                for (i = 0; i < order1; ++i)
                {
                    (*Mat)(order0 *i + 1, i * order0 + 1) = 1.0;
                }
            }

            if (m_base[1]->GetBasisType() == LibUtilities::eFourier)
            {
                for (i = 0; i < order0; ++i)
                {
                    (*Mat)(order0 + i, order0 + i) = 1.0;
                }
            }
            break;
        }
        case eFwdTrans:
        {
            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
            StdMatrixKey iprodkey(eIProductWRTBase, DetShapeType(), *this);
            DNekMat &Iprod = *GetStdMatrix(iprodkey);
            StdMatrixKey imasskey(eInvMass, DetShapeType(), *this);
            DNekMat &Imass = *GetStdMatrix(imasskey);

            (*Mat) = Imass * Iprod;
            break;
        }
        case eGaussDG:
        {
            ConstFactorMap factors = mkey.GetConstFactors();

            int edge    = (int)factors[StdRegions::eFactorGaussEdge];
            int dir     = (edge + 1) % 2;
            int nCoeffs = m_base[dir]->GetNumModes();

            const LibUtilities::PointsKey BS_p(
                nCoeffs, LibUtilities::eGaussGaussLegendre);
            const LibUtilities::BasisKey BS_k(LibUtilities::eGauss_Lagrange,
                                              nCoeffs, BS_p);

            Array<OneD, NekDouble> coords(1, 0.0);
            coords[0] = (edge == 0 || edge == 3) ? -1.0 : 1.0;

            LibUtilities::BasisSharedPtr basis =
                LibUtilities::BasisManager()[BS_k];
            DNekMatSharedPtr m_Ix = basis->GetI(coords);

            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(1.0, nCoeffs);
            Vmath::Vcopy(nCoeffs, m_Ix->GetPtr(), 1, Mat->GetPtr(), 1);
            break;
        }
        default:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
            break;
        }
    }

    return Mat;
}

DNekMatSharedPtr StdQuadExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return GenMatrix(mkey);
}

///////////////////////////////////
// Operator evaluation functions //
///////////////////////////////////

void StdQuadExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                      const StdMatrixKey &mkey)
{
    // Generate an orthonogal expansion
    int qa       = m_base[0]->GetNumPoints();
    int qb       = m_base[1]->GetNumPoints();
    int nmodes_a = m_base[0]->GetNumModes();
    int nmodes_b = m_base[1]->GetNumModes();
    int nmodes   = min(nmodes_a, nmodes_b);
    // Declare orthogonal basis.
    LibUtilities::PointsKey pa(qa, m_base[0]->GetPointsType());
    LibUtilities::PointsKey pb(qb, m_base[1]->GetPointsType());

    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, nmodes_a, pa);
    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, nmodes_b, pb);
    StdQuadExp OrthoExp(Ba, Bb);

    // SVV parameters loaded from the .xml case file
    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());

    // project onto modal  space.
    OrthoExp.FwdTrans(array, orthocoeffs);

    if (mkey.ConstFactorExists(
            eFactorSVVPowerKerDiffCoeff)) // Rodrigo's power kernel
    {
        NekDouble cutoff = mkey.GetConstFactor(eFactorSVVCutoffRatio);
        NekDouble SvvDiffCoeff =
            mkey.GetConstFactor(eFactorSVVPowerKerDiffCoeff) *
            mkey.GetConstFactor(eFactorSVVDiffCoeff);

        for (int j = 0; j < nmodes_a; ++j)
        {
            for (int k = 0; k < nmodes_b; ++k)
            {
                // linear space but makes high modes very negative
                NekDouble fac = std::max(
                    pow((1.0 * j) / (nmodes_a - 1), cutoff * nmodes_a),
                    pow((1.0 * k) / (nmodes_b - 1), cutoff * nmodes_b));

                orthocoeffs[j * nmodes_b + k] *= SvvDiffCoeff * fac;
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

        for (int j = 0; j < nmodes_a; ++j)
        {
            for (int k = 0; k < nmodes_b; ++k)
            {
                int maxjk = max(j, k);
                maxjk     = min(maxjk, kSVVDGFiltermodesmax - 1);

                orthocoeffs[j * nmodes_b + k] *=
                    SvvDiffCoeff * kSVVDGFilter[max_ab][maxjk];
            }
        }
    }
    else
    {
        NekDouble SvvDiffCoeff = mkey.GetConstFactor(eFactorSVVDiffCoeff);
        // Exponential Kernel implementation
        int cutoff = (int)(mkey.GetConstFactor(eFactorSVVCutoffRatio) *
                           min(nmodes_a, nmodes_b));

        //------"New" Version August 22nd '13--------------------
        for (int j = 0; j < nmodes_a; ++j)
        {
            for (int k = 0; k < nmodes_b; ++k)
            {
                if (j + k >= cutoff) // to filter out only the "high-modes"
                {
                    orthocoeffs[j * nmodes_b + k] *=
                        (SvvDiffCoeff *
                         exp(-(j + k - nmodes) * (j + k - nmodes) /
                             ((NekDouble)((j + k - cutoff + 1) *
                                          (j + k - cutoff + 1)))));
                }
                else
                {
                    orthocoeffs[j * nmodes_b + k] *= 0.0;
                }
            }
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdQuadExp::v_ExponentialFilter(Array<OneD, NekDouble> &array,
                                     const NekDouble alpha,
                                     const NekDouble exponent,
                                     const NekDouble cutoff)
{
    // Generate an orthogonal expansion
    int qa      = m_base[0]->GetNumPoints();
    int qb      = m_base[1]->GetNumPoints();
    int nmodesA = m_base[0]->GetNumModes();
    int nmodesB = m_base[1]->GetNumModes();
    int P       = nmodesA - 1;
    int Q       = nmodesB - 1;

    // Declare orthogonal basis.
    LibUtilities::PointsKey pa(qa, m_base[0]->GetPointsType());
    LibUtilities::PointsKey pb(qb, m_base[1]->GetPointsType());

    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, nmodesA, pa);
    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, nmodesB, pb);
    StdQuadExp OrthoExp(Ba, Bb);

    // Cutoff
    int Pcut = cutoff * P;
    int Qcut = cutoff * Q;

    // Project onto orthogonal space.
    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
    OrthoExp.FwdTrans(array, orthocoeffs);

    //
    NekDouble fac, fac1, fac2;
    for (int i = 0; i < nmodesA; ++i)
    {
        for (int j = 0; j < nmodesB; ++j)
        {
            // to filter out only the "high-modes"
            if (i > Pcut || j > Qcut)
            {
                fac1 = (NekDouble)(i - Pcut) / ((NekDouble)(P - Pcut));
                fac2 = (NekDouble)(j - Qcut) / ((NekDouble)(Q - Qcut));
                fac  = max(fac1, fac2);
                fac  = pow(fac, exponent);
                orthocoeffs[i * nmodesB + j] *= exp(-alpha * fac);
            }
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdQuadExp::v_ReduceOrderCoeffs(
    int numMin, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int n_coeffs = inarray.size();

    Array<OneD, NekDouble> coeff(n_coeffs);
    Array<OneD, NekDouble> coeff_tmp(n_coeffs, 0.0);
    Array<OneD, NekDouble> tmp;
    Array<OneD, NekDouble> tmp2;

    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();
    int numMax  = nmodes0;

    Vmath::Vcopy(n_coeffs, inarray, 1, coeff_tmp, 1);

    const LibUtilities::PointsKey Pkey0(nmodes0,
                                        LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey Pkey1(nmodes1,
                                        LibUtilities::eGaussLobattoLegendre);

    LibUtilities::BasisKey b0(m_base[0]->GetBasisType(), nmodes0, Pkey0);
    LibUtilities::BasisKey b1(m_base[1]->GetBasisType(), nmodes1, Pkey1);

    LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A, nmodes0, Pkey0);
    LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_A, nmodes1, Pkey1);

    LibUtilities::InterpCoeff2D(b0, b1, coeff_tmp, bortho0, bortho1, coeff);

    Vmath::Zero(n_coeffs, coeff_tmp, 1);

    int cnt = 0;
    for (int i = 0; i < numMin + 1; ++i)
    {
        Vmath::Vcopy(numMin, tmp = coeff + cnt, 1, tmp2 = coeff_tmp + cnt, 1);

        cnt = i * numMax;
    }

    LibUtilities::InterpCoeff2D(bortho0, bortho1, coeff_tmp, b0, b1, outarray);
}

void StdQuadExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray,
                                const StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdQuadExp::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    StdQuadExp::v_LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdQuadExp::v_LaplacianMatrixOp(
    const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void StdQuadExp::v_WeakDerivMatrixOp(
    const int i, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void StdQuadExp::v_HelmholtzMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    StdQuadExp::v_HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdQuadExp::v_GetSimplexEquiSpacedConnectivity(
    Array<OneD, int> &conn, [[maybe_unused]] bool standard)
{
    int np1 = m_base[0]->GetNumPoints();
    int np2 = m_base[1]->GetNumPoints();
    int np  = max(np1, np2);

    conn = Array<OneD, int>(6 * (np - 1) * (np - 1));

    int row   = 0;
    int rowp1 = 0;
    int cnt   = 0;
    for (int i = 0; i < np - 1; ++i)
    {
        rowp1 += np;
        for (int j = 0; j < np - 1; ++j)
        {
            conn[cnt++] = row + j;
            conn[cnt++] = row + j + 1;
            conn[cnt++] = rowp1 + j;

            conn[cnt++] = rowp1 + j + 1;
            conn[cnt++] = rowp1 + j;
            conn[cnt++] = row + j + 1;
        }
        row += np;
    }
}

} // namespace Nektar::StdRegions
