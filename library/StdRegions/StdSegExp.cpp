///////////////////////////////////////////////////////////////////////////////
//
// File: StdSegExp.cpp
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
// Description: Routines within Standard Segment Expansions
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LibUtilities/Foundations/ManagerAccess.h> // for PointsManager, etc
#include <StdRegions/StdSegExp.h>

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

/** \brief Constructor using BasisKey class for quadrature points and
 *  order definition
 *
 *  \param Ba BasisKey class definition containing order and quadrature
 *  points.
 */
StdSegExp::StdSegExp(const LibUtilities::BasisKey &Ba)
    : StdExpansion(Ba.GetNumModes(), 1, Ba),
      StdExpansion1D(Ba.GetNumModes(), Ba)
{
    // cache integration weights for future use
    m_weights.push_back(m_base[0]->GetW());
}

/** \brief Return Shape of region, using  ShapeType enum list.
 *  i.e. Segment
 */
LibUtilities::ShapeType StdSegExp::v_DetShapeType() const
{
    return LibUtilities::eSegment;
}

bool StdSegExp::v_IsBoundaryInteriorExpansion() const
{

    bool returnval = false;

    if (m_base[0]->GetBasisType() == LibUtilities::eModified_A)
    {
        returnval = true;
    }

    if (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
    {
        returnval = true;
    }

    return returnval;
}

//---------------------------------------------------------------------
// Differentiation Methods
//---------------------------------------------------------------------
/** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the physical
 *  quadrature points given by \a inarray and return in \a outarray.
 *
 *  This is a wrapper around StdExpansion1D::Tensor_Deriv
 *  \param inarray array of a function evaluated at the quadrature points
 *  \param  outarray the resulting array of the derivative \f$
 *  du/d_{\xi_1}|_{\xi_{1i}} \f$ will be stored in the array \a outarra
 */

void StdSegExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &out_d0,
                               [[maybe_unused]] Array<OneD, NekDouble> &out_d1,
                               [[maybe_unused]] Array<OneD, NekDouble> &out_d2)
{
    PhysTensorDeriv(inarray, out_d0);
}

//---------------------------------------------------------------------
// Transforms
//---------------------------------------------------------------------

/** \brief Backward transform from coefficient space given
 *  in \a inarray and evaluate at the physical quadrature
 *  points \a outarray
 *
 *  Operation can be evaluated as \f$ u(\xi_{1i}) =
 *  \sum_{p=0}^{order-1} \hat{u}_p \phi_p(\xi_{1i}) \f$ or equivalently
 *  \f$ {\bf u} = {\bf B}^T {\bf \hat{u}} \f$ where
 *  \f${\bf B}[i][j] = \phi_i(\xi_{1j}), \mbox{\_coeffs}[p] = {\bf
 *  \hat{u}}[p] \f$
 *
 *  The function takes the coefficient array \a inarray as
 *  input for the transformation
 *
 *  \param inarray: the coeffficients of the expansion
 *
 *  \param outarray: the resulting array of the values of the function at
 *  the physical quadrature points will be stored in the array \a outarray
 */

void StdSegExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();

    if (m_base[0]->Collocation())
    {
        std::memcpy(outarray.data(), inarray.data(),
                    nquad0 * sizeof(NekDouble));
    }
    else
    {
        const Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();

        int nmodes0 = m_base[0]->GetNumModes();

        // Switch statment using boost_pp and macros. This unfolls intwo a
        // nested swtich statement where the outer swtich statement runs
        // from SMIN to SMAX for modal order and the inner switch
        // statemets run from the outer value of the case to 2*SMAX for
        // the quadrature order. If you want to see it unwrapped compile
        // in verbose mode and add --preprocess to the c++ command.
        // Default case
#undef BWDTRANS_DEF
#define BWDTRANS_DEF                                                           \
    BwdTransSegKernel(nmodes0, nquad0, (const vec_t *)base0.data(),            \
                      (const vec_t *)inarray.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef BWDTRANS_Q
#define BWDTRANS_Q(r, i)                                                       \
    case NQ(i):                                                                \
        BwdTransSegKernel(NM(i), NQ(i), (const vec_t *)base0.data(),           \
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
        switch (nmodes0)
        {
            BOOST_PP_FOR((SMIN, 0, SMAX), STDLEV2TEST, STDLEV2UPDATE,
                         BWDTRANS_M)
            default:
                BWDTRANS_DEF;
                break;
        }
    }
}

void StdSegExp::v_ReduceOrderCoeffs(int numMin,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    int n_coeffs = inarray.size();

    Array<OneD, NekDouble> coeff(n_coeffs);
    Array<OneD, NekDouble> coeff_tmp(n_coeffs, 0.0);
    Array<OneD, NekDouble> tmp;
    Array<OneD, NekDouble> tmp2;

    int nmodes0 = m_base[0]->GetNumModes();

    Vmath::Vcopy(n_coeffs, inarray, 1, coeff_tmp, 1);

    const LibUtilities::PointsKey Pkey0(nmodes0,
                                        LibUtilities::eGaussLobattoLegendre);

    LibUtilities::BasisKey b0(m_base[0]->GetBasisType(), nmodes0, Pkey0);

    LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A, nmodes0, Pkey0);

    LibUtilities::InterpCoeff1D(b0, coeff_tmp, bortho0, coeff);

    Vmath::Zero(n_coeffs, coeff_tmp, 1);

    Vmath::Vcopy(numMin, tmp = coeff, 1, tmp2 = coeff_tmp, 1);

    LibUtilities::InterpCoeff1D(bortho0, coeff_tmp, b0, outarray);
}

void StdSegExp::v_FwdTransBndConstrained(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    if (m_base[0]->Collocation())
    {
        Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
    }
    else
    {
        int nInteriorDofs = m_ncoeffs - 2;
        int offset        = 0;

        switch (m_base[0]->GetBasisType())
        {
            case LibUtilities::eGLL_Lagrange:
            {
                offset = 1;
            }
            break;
            case LibUtilities::eGauss_Lagrange:
            {
                nInteriorDofs = m_ncoeffs;
                offset        = 0;
            }
            break;
            case LibUtilities::eModified_A:
            case LibUtilities::eModified_B:
            {
                offset = 2;
            }
            break;
            default:
                ASSERTL0(false, "This type of FwdTrans is not defined for this "
                                "expansion type");
        }

        fill(outarray.data(), outarray.data() + m_ncoeffs, 0.0);

        if (m_base[0]->GetBasisType() != LibUtilities::eGauss_Lagrange)
        {
            outarray[GetVertexMap(0)] = inarray[0];
            outarray[GetVertexMap(1)] = inarray[m_base[0]->GetNumPoints() - 1];

            if (m_ncoeffs > 2)
            {
                // ideally, we would like to have tmp0 to be replaced by
                // outarray (currently MassMatrixOp does not allow aliasing)
                Array<OneD, NekDouble> tmp0(m_ncoeffs);
                Array<OneD, NekDouble> tmp1(m_ncoeffs);

                StdMatrixKey masskey(eMass, v_DetShapeType(), *this);
                MassMatrixOp(outarray, tmp0, masskey);
                v_IProductWRTBase(inarray, tmp1);

                Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                // get Mass matrix inverse (only of interior DOF)
                DNekMatSharedPtr matsys =
                    (m_stdStaticCondMatrixManager[masskey])->GetBlock(1, 1);

                Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, 1.0,
                            &(matsys->GetPtr())[0], nInteriorDofs,
                            tmp1.data() + offset, 1, 0.0,
                            outarray.data() + offset, 1);
            }
        }
        else
        {
            v_FwdTrans(inarray, outarray);
        }
    }
}

//---------------------------------------------------------------------
// Inner product functions
//---------------------------------------------------------------------
/** \brief Inner product of \a inarray over region with respect to the
 *  expansion basis (this)->m_base[0] and return in \a outarray
 *
 *  @param base0 - An array containing the values of the basis in the
 *  0-direction at the quarature poitns
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
void StdSegExp::v_IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, const NekDouble> &jac,
    const bool Deformed)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int order0 = m_base[0]->GetNumModes();

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
    IProductSegKernel<false, false, true>(                                     \
        order0, nquad0, (const vec_t *)inarray.data(),                         \
        (const vec_t *)base0.data(), (const vec_t *)m_weights[0].data(),       \
        (const vec_t *)jac.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductSegKernel<false, false, true>(                                 \
            NM(i), NQ(i), (const vec_t *)inarray.data(),                       \
            (const vec_t *)base0.data(), (const vec_t *)m_weights[0].data(),   \
            (const vec_t *)jac.data(), (vec_t *)outarray.data());              \
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
        switch (order0)
        {
            BOOST_PP_FOR((SMIN, 0, SMAX), STDLEV2TEST, STDLEV2UPDATE,
                         IPRODUCTWRTBASE_M)
            default:
                IPRODUCTWRTBASE_DEF;
                break;
        }
    }
    else // non-deformed case
    {
        // Default case
#undef IPRODUCTWRTBASE_DEF
#define IPRODUCTWRTBASE_DEF                                                    \
    IProductSegKernel<false, false, false>(                                    \
        order0, nquad0, (const vec_t *)inarray.data(),                         \
        (const vec_t *)base0.data(), (const vec_t *)m_weights[0].data(),       \
        (const vec_t *)jac.data(), (vec_t *)outarray.data())

        // Inner loop case over quarature points
#undef IPRODUCTWRTBASE_Q
#define IPRODUCTWRTBASE_Q(r, i)                                                \
    case NQ(i):                                                                \
        IProductSegKernel<false, false, false>(                                \
            NM(i), NQ(i), (const vec_t *)inarray.data(),                       \
            (const vec_t *)base0.data(), (const vec_t *)m_weights[0].data(),   \
            (const vec_t *)jac.data(), (vec_t *)outarray.data());              \
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
        switch (order0)
        {
            BOOST_PP_FOR((SMIN, 0, SMAX), STDLEV2TEST, STDLEV2UPDATE,
                         IPRODUCTWRTBASE_M)
            default:
                IPRODUCTWRTBASE_DEF;
                break;
        }
    }
}

void StdSegExp::v_IProductWRTDerivBase(
    [[maybe_unused]] const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    ASSERTL1(dir == 0, "input dir is out of range");
    const Array<OneD, const NekDouble> one(1, 1.0);
    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), inarray, outarray, one,
                            false);
}

//----------------------------
// Evaluation
//----------------------------
void StdSegExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble> &xi,
                                         Array<OneD, NekDouble> &eta)
{
    eta[0] = xi[0];
}

void StdSegExp::v_LocCollapsedToLocCoord(
    const Array<OneD, const NekDouble> &eta, Array<OneD, NekDouble> &xi)
{
    xi[0] = eta[0];
}

void StdSegExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
{
    int nquad             = m_base[0]->GetNumPoints();
    const NekDouble *base = m_base[0]->GetBdata().data();

    ASSERTL2(mode <= m_ncoeffs,
             "calling argument mode is larger than total expansion order");

    Vmath::Vcopy(nquad, (NekDouble *)base + mode * nquad, 1, &outarray[0], 1);
}

NekDouble StdSegExp::v_PhysEvaluateBasis(
    const Array<OneD, const NekDouble> &coords, int mode)
{
    return StdExpansion::BaryEvaluateBasis<0>(coords[0], mode);
}

NekDouble StdSegExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    return StdExpansion1D::BaryTensorDeriv(coord, inarray, firstOrderDerivs);
}

NekDouble StdSegExp::v_PhysEvalFirstSecondDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs,
    std::array<NekDouble, 6> &secondOrderDerivs)
{
    return StdExpansion1D::BaryTensorDeriv(coord, inarray, firstOrderDerivs,
                                           secondOrderDerivs);
}
void StdSegExp::v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    [[maybe_unused]] const StdMatrixKey &mkey)
{
    int nquad = m_base[0]->GetNumPoints();

    Array<OneD, NekDouble> physValues(nquad);
    Array<OneD, NekDouble> dPhysValuesdx(nquad);

    v_BwdTrans(inarray, physValues);

    // Laplacian matrix operation
    v_PhysDeriv(physValues, dPhysValuesdx, NullNekDouble1DArray,
                NullNekDouble1DArray);
    v_IProductWRTBase(dPhysValuesdx, outarray);
}

void StdSegExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void StdSegExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    const StdMatrixKey &mkey)
{
    int nquad = m_base[0]->GetNumPoints();

    Array<OneD, NekDouble> physValues(nquad);
    Array<OneD, NekDouble> dPhysValuesdx(nquad);
    Array<OneD, NekDouble> wsp(m_ncoeffs);

    v_BwdTrans(inarray, physValues);

    // mass matrix operation
    v_IProductWRTBase(physValues, wsp);

    // Laplacian matrix operation
    const Array<OneD, const NekDouble> one(1, 1.0);
    v_PhysDeriv(physValues, dPhysValuesdx, NullNekDouble1DArray,
                NullNekDouble1DArray);
    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), dPhysValuesdx, outarray,
                            one, false);
    Blas::Daxpy(m_ncoeffs, mkey.GetConstFactor(eFactorLambda), wsp.data(), 1,
                outarray.data(), 1);
}

void StdSegExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                     const StdMatrixKey &mkey)
{
    // Generate an orthogonal expansion
    int nq     = m_base[0]->GetNumPoints();
    int nmodes = m_base[0]->GetNumModes();
    // Declare orthogonal basis.
    LibUtilities::PointsKey pKey(nq, m_base[0]->GetPointsType());

    LibUtilities::BasisKey B(LibUtilities::eOrtho_A, nmodes, pKey);
    StdSegExp OrthoExp(B);

    // SVV parameters loaded from the .xml case file
    NekDouble SvvDiffCoeff = mkey.GetConstFactor(eFactorSVVDiffCoeff);
    int cutoff = (int)(mkey.GetConstFactor(eFactorSVVCutoffRatio)) * nmodes;

    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());

    // project onto modal  space.
    OrthoExp.FwdTrans(array, orthocoeffs);

    //
    for (int j = 0; j < nmodes; ++j)
    {
        if (j >= cutoff) // to filter out only the "high-modes"
        {
            orthocoeffs[j] *=
                (SvvDiffCoeff *
                 exp(-(j - nmodes) * (j - nmodes) /
                     ((NekDouble)((j - cutoff + 1) * (j - cutoff + 1)))));
        }
        else
        {
            orthocoeffs[j] *= 0.0;
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdSegExp::v_ExponentialFilter(Array<OneD, NekDouble> &array,
                                    const NekDouble alpha,
                                    const NekDouble exponent,
                                    const NekDouble cutoff)
{
    // Generate an orthogonal expansion
    int nq     = m_base[0]->GetNumPoints();
    int nmodes = m_base[0]->GetNumModes();
    int P      = nmodes - 1;
    // Declare orthogonal basis.
    LibUtilities::PointsKey pKey(nq, m_base[0]->GetPointsType());

    LibUtilities::BasisKey B(LibUtilities::eOrtho_A, nmodes, pKey);
    StdSegExp OrthoExp(B);

    // Cutoff
    int Pcut = cutoff * P;

    // Project onto orthogonal space.
    Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
    OrthoExp.FwdTrans(array, orthocoeffs);

    //
    NekDouble fac;
    for (int j = 0; j < nmodes; ++j)
    {
        // to filter out only the "high-modes"
        if (j > Pcut)
        {
            fac = (NekDouble)(j - Pcut) / ((NekDouble)(P - Pcut));
            fac = pow(fac, exponent);
            orthocoeffs[j] *= exp(-alpha * fac);
        }
    }

    // backward transform to physical space
    OrthoExp.BwdTrans(orthocoeffs, array);
}

void StdSegExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                            [[maybe_unused]] Array<OneD, NekDouble> &coords_1,
                            [[maybe_unused]] Array<OneD, NekDouble> &coords_2)
{
    Vmath::Vcopy(GetNumPoints(0), (m_base[0]->GetZ()).data(), 1, &coords_0[0],
                 1);
}

//---------------------------------------------------------------------
// Helper functions
//---------------------------------------------------------------------

int StdSegExp::v_GetNverts() const
{
    return 2;
}

int StdSegExp::v_GetNtraces() const
{
    return 2;
}

int StdSegExp::v_GetTraceNcoeffs([[maybe_unused]] const int i) const
{
    return 1;
}

int StdSegExp::v_GetTraceIntNcoeffs([[maybe_unused]] const int i) const
{
    return 0;
}

int StdSegExp::v_GetTraceNumPoints([[maybe_unused]] const int i) const
{
    return 1;
}

int StdSegExp::v_NumBndryCoeffs() const
{
    return 2;
}

int StdSegExp::v_NumDGBndryCoeffs() const
{
    return 2;
}

int StdSegExp::v_CalcNumberOfCoefficients(
    const std::vector<unsigned int> &nummodes, int &modes_offset)
{
    int nmodes = nummodes[modes_offset];
    modes_offset += 1;

    return nmodes;
}

//---------------------------------------------------------------------
// Wrapper functions
//---------------------------------------------------------------------

DNekMatSharedPtr StdSegExp::v_GenMatrix(const StdMatrixKey &mkey)
{
    DNekMatSharedPtr Mat;
    MatrixType mattype;

    switch (mattype = mkey.GetMatrixType())
    {
        case ePhysInterpToEquiSpaced:
        {
            int nq0 = m_base[0]->GetNumPoints();
            int nq;

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = nq0;
            }

            int neq = LibUtilities::StdSegData::getNumberOfCoefficients(nq);
            Array<OneD, NekDouble> coords(1);
            DNekMatSharedPtr I;
            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0);

            for (int i = 0; i < neq; ++i)
            {
                coords[0] = -1.0 + 2 * i / (NekDouble)(neq - 1);
                I         = m_base[0]->GetI(coords);
                Vmath::Vcopy(nq0, I->GetRawPtr(), 1, Mat->GetRawPtr() + i, neq);
            }
        }
        break;
        case ePhysInterpToGLL:
        {
            int nq0 = m_base[0]->GetNumPoints();
            int nq;

            // take definition from key
            if (mkey.ConstFactorExists(eFactorConst))
            {
                nq = (int)mkey.GetConstFactor(eFactorConst);
            }
            else
            {
                nq = nq0;
            }

            int neq = LibUtilities::StdSegData::getNumberOfCoefficients(nq);
            Array<OneD, NekDouble> coords(1);
            DNekMatSharedPtr I;
            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(neq, nq0);

            const LibUtilities::PointsKey key(
                neq, LibUtilities::eGaussLobattoLegendre);

            Array<OneD, const NekDouble> z;
            LibUtilities::PointsManager()[key]->GetPoints(z);

            for (int i = 0; i < neq; ++i)
            {
                coords[0] = z[i];
                I         = m_base[0]->GetI(coords);
                Vmath::Vcopy(nq0, I->GetRawPtr(), 1, Mat->GetRawPtr() + i, neq);
            }
        }
        break;
        case eFwdTrans:
        {
            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
            StdMatrixKey iprodkey(eIProductWRTBase, v_DetShapeType(), *this);
            DNekMat &Iprod = *GetStdMatrix(iprodkey);
            StdMatrixKey imasskey(eInvMass, v_DetShapeType(), *this);
            DNekMat &Imass = *GetStdMatrix(imasskey);

            (*Mat) = Imass * Iprod;
        }
        break;
        default:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);

            if (mattype == eMass)
            {
                // For Fourier basis set the imaginary component
                // of mean mode to have a unit diagonal component
                // in mass matrix
                if (m_base[0]->GetBasisType() == LibUtilities::eFourier)
                {
                    (*Mat)(1, 1) = 1.0;
                }
            }
        }
        break;
    }

    return Mat;
}

DNekMatSharedPtr StdSegExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

//---------------------------------------------------------------------
// Mappings
//---------------------------------------------------------------------

void StdSegExp::v_GetBoundaryMap(Array<OneD, unsigned int> &outarray)
{
    if (outarray.size() != NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
    }
    const LibUtilities::BasisType Btype = GetBasisType(0);
    int nummodes                        = m_base[0]->GetNumModes();

    outarray[0] = 0;

    switch (Btype)
    {
        case LibUtilities::eGLL_Lagrange:
        case LibUtilities::eGauss_Lagrange:
        case LibUtilities::eChebyshev:
        case LibUtilities::eFourier:
            outarray[1] = nummodes - 1;
            break;
        case LibUtilities::eModified_A:
        case LibUtilities::eModified_B:
            outarray[1] = 1;
            break;
        default:
            ASSERTL0(0, "Mapping array is not defined for this expansion");
            break;
    }
}

void StdSegExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
{
    int i;
    if (outarray.size() != GetNcoeffs() - NumBndryCoeffs())
    {
        outarray = Array<OneD, unsigned int>(GetNcoeffs() - NumBndryCoeffs());
    }
    const LibUtilities::BasisType Btype = GetBasisType(0);

    switch (Btype)
    {
        case LibUtilities::eGLL_Lagrange:
        case LibUtilities::eGauss_Lagrange:
        case LibUtilities::eChebyshev:
        case LibUtilities::eFourier:
            for (i = 0; i < GetNcoeffs() - 2; i++)
            {
                outarray[i] = i + 1;
            }
            break;
        case LibUtilities::eModified_A:
        case LibUtilities::eModified_B:
            for (i = 0; i < GetNcoeffs() - 2; i++)
            {
                outarray[i] = i + 2;
            }
            break;
        default:
            ASSERTL0(0, "Mapping array is not defined for this expansion");
            break;
    }
}

int StdSegExp::v_GetVertexMap(int localVertexId,
                              [[maybe_unused]] bool useCoeffPacking)
{
    ASSERTL0((localVertexId == 0) || (localVertexId == 1),
             "local vertex id"
             "must be between 0 or 1");

    int localDOF = localVertexId;

    if ((m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange) &&
        (localVertexId == 1))
    {
        localDOF = m_base[0]->GetNumModes() - 1;
    }
    return localDOF;
}

void StdSegExp::v_GetSimplexEquiSpacedConnectivity(
    Array<OneD, int> &conn, [[maybe_unused]] bool standard)
{
    int np = m_base[0]->GetNumPoints();

    conn    = Array<OneD, int>(2 * (np - 1));
    int cnt = 0;
    for (int i = 0; i < np - 1; ++i)
    {
        conn[cnt++] = i;
        conn[cnt++] = i + 1;
    }
}

/**  \brief Get the map of the coefficient location to teh
 *    local trace coefficients
 */

void StdSegExp::v_GetTraceCoeffMap(const unsigned int traceid,
                                   Array<OneD, unsigned int> &maparray)
{
    int order0 = m_base[0]->GetNumModes();

    ASSERTL0(traceid < 2, "eid must be between 0 and 1");

    if (maparray.size() != 1)
    {
        maparray = Array<OneD, unsigned int>(1);
    }

    const LibUtilities::BasisType bType = GetBasisType(0);

    if (bType == LibUtilities::eModified_A)
    {
        maparray[0] = (traceid == 0) ? 0 : 1;
    }
    else if (bType == LibUtilities::eGLL_Lagrange ||
             bType == LibUtilities::eGauss_Lagrange)
    {
        maparray[0] = (traceid == 0) ? 0 : order0 - 1;
    }
    else
    {
        ASSERTL0(false, "Unknown Basis");
    }
}

void StdSegExp::v_GetTraceToElementMap(const int tid,
                                       Array<OneD, unsigned int> &maparray,
                                       Array<OneD, int> &signarray,
                                       [[maybe_unused]] Orientation orient,
                                       [[maybe_unused]] int P,
                                       [[maybe_unused]] int Q)
{
    v_GetTraceCoeffMap(tid, maparray);

    if (signarray.size() != 1)
    {
        signarray = Array<OneD, int>(1, 1);
    }
    else
    {
        signarray[0] = 1;
    }
}

void StdSegExp::v_GetElmtTraceToTraceMap(
    [[maybe_unused]] const unsigned int eid,
    Array<OneD, unsigned int> &maparray, Array<OneD, int> &signarray,
    [[maybe_unused]] Orientation orient, [[maybe_unused]] int P,
    [[maybe_unused]] int Q)
{
    // parameters for higher dimnesion traces
    if (maparray.size() != 1)
    {
        maparray = Array<OneD, unsigned int>(1);
    }

    maparray[0] = 0;

    if (signarray.size() != 1)
    {
        signarray = Array<OneD, int>(1, 1);
    }
    else
    {
        signarray[0] = 1;
    }
}

} // namespace Nektar::StdRegions
