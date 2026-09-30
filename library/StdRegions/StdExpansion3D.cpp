///////////////////////////////////////////////////////////////////////////////
//
// File: StdExpansion3D.cpp
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
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion3D.h>

#include <LibUtilities/BasicUtils/NekInline.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <StdRegions/Operators/SwitchLevel1.h>
#include <StdRegions/Operators/SwitchLevel2.h>

#ifdef max
#undef max
#endif

namespace Nektar::StdRegions
{
// Declaretion of scalar routine
using vec_t = tinysimd::scalarT<double>;
#include <StdRegions/Operators/PhysDerivSumFacStdKernels.hpp>

StdExpansion3D::StdExpansion3D(
    [[maybe_unused]] int numcoeffs,
    [[maybe_unused]] const LibUtilities::BasisKey &Ba,
    [[maybe_unused]] const LibUtilities::BasisKey &Bb,
    [[maybe_unused]] const LibUtilities::BasisKey &Bc)
{
}

void StdExpansion3D::IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &base1,
    const Array<OneD, const NekDouble> &base2,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, NekDouble> &jac,
    const bool Deformed, [[maybe_unused]] bool CollDir0,
    [[maybe_unused]] bool CollDir1, [[maybe_unused]] bool CollDir2)
{
    v_IProductWRTBaseKernel(base0, base1, base2, inarray, outarray, jac,
                            Deformed, CollDir0, CollDir1, CollDir2);
}

void StdExpansion3D::v_GenStdMatBwdDeriv(const int dir, DNekMatSharedPtr &mat)
{
    ASSERTL1((dir == 0) || (dir == 1) || (dir == 2), "Invalid direction.");

    const int nq0 = m_base[0]->GetNumPoints();
    const int nq1 = m_base[1]->GetNumPoints();
    const int nq2 = m_base[2]->GetNumPoints();
    const int nq  = nq0 * nq1 * nq2;

    const bool CollDir0 = m_base[0]->Collocation();
    const bool CollDir1 = m_base[1]->Collocation();
    const bool CollDir2 = m_base[2]->Collocation();

    Array<OneD, NekDouble> in(nq, 0.0);
    Array<OneD, NekDouble> out(m_ncoeffs);
    Array<OneD, NekDouble> one(1, 1.0);

    for (int i = 0; i < nq; i++)
    {
        int l = i % nq0;
        int m = (i / nq0) % nq1;
        int n = i / (nq0 * nq1);

        // initialise with inverse of weights t
        in[i] = 1.0 / (m_weights[0][l] * m_weights[1][m] * m_weights[2][n]);

        // do standard iproduct
        if (dir == 0)
        {
            v_IProductWRTBaseKernel(m_base[0]->GetDbdata(),
                                    m_base[1]->GetBdata(),
                                    m_base[2]->GetBdata(), in, out, one, false,
                                    false, CollDir1, CollDir2);
        }
        else if (dir == 1)
        {
            v_IProductWRTBaseKernel(m_base[0]->GetBdata(),
                                    m_base[1]->GetDbdata(),
                                    m_base[2]->GetBdata(), in, out, one, false,
                                    CollDir0, false, CollDir2);
        }
        else // dir == 2
        {
            v_IProductWRTBaseKernel(m_base[0]->GetBdata(),
                                    m_base[1]->GetBdata(),
                                    m_base[2]->GetDbdata(), in, out, one, false,
                                    CollDir0, CollDir1, false);
        }
        in[i] = 0.0;

        for (int j = 0; j < m_ncoeffs; j++)
        {
            (*mat)(j, i) = out[j];
        }
    }
}

void StdExpansion3D::PhysTensorDeriv(
    const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &out_d0,
    Array<OneD, NekDouble> &out_d1, Array<OneD, NekDouble> &out_d2)
{
    const int nquad0 = m_base[0]->GetNumPoints();
    const int nquad1 = m_base[1]->GetNumPoints();
    const int nquad2 = m_base[2]->GetNumPoints();

    bool Deriv0         = (out_d0.size() > 0);
    bool Deriv1         = (out_d1.size() > 0);
    bool Deriv2         = (out_d2.size() > 0);
    const NekDouble *D0 = m_base[0]->GetD()->GetRawPtr();
    const NekDouble *D1 = m_base[1]->GetD()->GetRawPtr();
    const NekDouble *D2 = m_base[2]->GetD()->GetRawPtr();

    Array<OneD, const NekDouble> intmp;
    // copy inarray data if inarray and outarray are the same.
    if ((inarray.data() == out_d0.data()) ||
        (inarray.data() == out_d1.data()) || (inarray.data() == out_d2.data()))
    {
        Array<OneD, NekDouble> wsp(nquad0 * nquad1 * nquad2);
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
#undef PHYSDERIV_DEF
#define PHYSDERIV_DEF                                                          \
    PhysDerivTensor3DKernel(nquad0, nquad1, nquad2,                            \
                            (const vec_t *)intmp.data(), (const vec_t *)D0,    \
                            (const vec_t *)D1, (const vec_t *)D2,              \
                            (vec_t *)out_d0.data(), (vec_t *)out_d1.data(),    \
                            (vec_t *)out_d2.data(), Deriv0, Deriv1, Deriv2)

    // Loop case over quarature points
#undef PHYSDERIV_Q
#define PHYSDERIV_Q(r, i)                                                      \
    case NQ1(i):                                                               \
        PhysDerivTensor3DKernel(                                               \
            NQ1(i), NQ1(i), NQ1(i), (const vec_t *)intmp.data(),               \
            (const vec_t *)D0, (const vec_t *)D1, (const vec_t *)D2,           \
            (vec_t *)out_d0.data(), (vec_t *)out_d1.data(),                    \
            (vec_t *)out_d2.data(), Deriv0, Deriv1, Deriv2);                   \
        break;

    // templated cases on  standard quadrature
    // usage where quad order goes from SMIN to SMAX
    if ((nquad0 == nquad1) && (nquad1 == nquad2))
    {
        switch (nquad0)
        {
            BOOST_PP_FOR((SMIN, SMAX), STDLEV1TEST, STDLEV1UPDATE, PHYSDERIV_Q);
            default:
                PHYSDERIV_DEF;
                break;
        }
    }
    else
    {
        PHYSDERIV_DEF;
    }
}

void StdExpansion3D::v_PhysDeriv(const int dir,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray)
{
    switch (dir)
    {
        case 0:
        {
            v_PhysDeriv(inarray, outarray, NullNekDouble1DArray,
                        NullNekDouble1DArray);
            break;
        }

        case 1:
        {
            v_PhysDeriv(inarray, NullNekDouble1DArray, outarray,
                        NullNekDouble1DArray);
            break;
        }

        case 2:
        {
            v_PhysDeriv(inarray, NullNekDouble1DArray, NullNekDouble1DArray,
                        outarray);
            break;
        }

        default:
        {
            ASSERTL1(false, "input dir is out of range");
        }
        break;
    }
}

NekDouble StdExpansion3D::v_StdPhysEvaluate(
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &physvals)
{
    Array<OneD, NekDouble> eta(3);

    WARNINGL2(coords[0] >= -1 - NekConstants::kNekZeroTol, "coord[0] < -1");
    WARNINGL2(coords[0] <= 1 + NekConstants::kNekZeroTol, "coord[0] >  1");
    WARNINGL2(coords[1] >= -1 - NekConstants::kNekZeroTol, "coord[1] < -1");
    WARNINGL2(coords[1] <= 1 + NekConstants::kNekZeroTol, "coord[1] >  1");
    WARNINGL2(coords[2] >= -1 - NekConstants::kNekZeroTol, "coord[2] < -1");
    WARNINGL2(coords[2] <= 1 + NekConstants::kNekZeroTol, "coord[2] >  1");

    // Obtain local collapsed coordinate from Cartesian coordinate.
    LocCoordToLocCollapsed(coords, eta);

    const int nq0 = m_base[0]->GetNumPoints();
    const int nq1 = m_base[1]->GetNumPoints();
    const int nq2 = m_base[2]->GetNumPoints();

    Array<OneD, NekDouble> wsp1(nq1 * nq2), wsp2(nq2);

    // Construct the 2D square...
    const NekDouble *ptr = &physvals[0];
    for (int i = 0; i < nq1 * nq2; ++i, ptr += nq0)
    {
        wsp1[i] = StdExpansion::BaryEvaluate<0>(eta[0], ptr);
    }

    for (int i = 0; i < nq2; ++i)
    {
        wsp2[i] = StdExpansion::BaryEvaluate<1>(eta[1], &wsp1[i * nq1]);
    }

    return StdExpansion::BaryEvaluate<2>(eta[2], &wsp2[0]);
}

/**
 * \f$
 * \begin{array}{rcl}
 * I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
 * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
 * \psi_{p}^{a}(\xi_{1i}) \psi_{q}^{a}(\xi_{2j}) \psi_{r}^{a}(\xi_{3k})
 * w_i w_j w_k u(\xi_{1,i} \xi_{2,j} \xi_{3,k})
 *
 * J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\xi_{1,i})
 *                   \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j})
 *                   \sum_{k=0}^{nq_2} \psi_{r}^a
 *                   u(\xi_{1i},\xi_{2j},\xi_{3k}) J_{i,j,k}
 * \end{array} \f$ \n
 * where
 * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3)
 *  = \psi_p^a( \xi_1) \psi_{q}^a(\xi_2) \psi_{r}^a(\xi_3) \f$ \n
 * which can be implemented as \n
 * \f$f_{r} (\xi_{3k})
 *  = \sum_{k=0}^{nq_3} \psi_{r}^a u(\xi_{1i},\xi_{2j}, \xi_{3k})
 * J_{i,j,k} = {\bf B_3 U}   \f$ \n
 * \f$ g_{q} (\xi_{3k})
 *  = \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2j}) f_{r}(\xi_{3k})
 *  = {\bf B_2 F}  \f$ \n
 * \f$ (\phi_{pqr}, u)_{\delta}
 *  = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k})  g_{q} (\xi_{3k})
 *  = {\bf B_1 G} \f$
 *
 * @param   inarray     Physical space function definition
 * @param   outarray    Inner product with respect to basis
 *
 *
 * This is a wrapper function around \a IProductWRTBaseKernel()
 */
void StdExpansion3D::v_IProductWRTBase(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    const bool CollDir0 = m_base[0]->Collocation();
    const bool CollDir1 = m_base[1]->Collocation();
    const bool CollDir2 = m_base[2]->Collocation();

    if (CollDir0 && CollDir1 && CollDir2)
    {
        MultiplyByStdQuadratureMetric(inarray, outarray);
    }
    else
    {
        const Array<OneD, const NekDouble> one(1, 1.0);
        v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                                m_base[2]->GetBdata(), inarray, outarray, one,
                                false, CollDir0, CollDir1, CollDir2);
    }
}

NekDouble StdExpansion3D::v_PhysEvaluateInterp(
    const Array<OneD, DNekMatSharedPtr> &I,
    const Array<OneD, const NekDouble> &physvals)
{
    NekDouble value;

    int Qx = m_base[0]->GetNumPoints();
    int Qy = m_base[1]->GetNumPoints();
    int Qz = m_base[2]->GetNumPoints();

    Array<OneD, NekDouble> sumFactorization_qr =
        Array<OneD, NekDouble>(Qy * Qz);
    Array<OneD, NekDouble> sumFactorization_r = Array<OneD, NekDouble>(Qz);

    // Lagrangian interpolation matrix
    NekDouble *interpolatingNodes = nullptr;

    // Interpolate first coordinate direction
    interpolatingNodes = &I[0]->GetPtr()[0];

    Blas::Dgemv('T', Qx, Qy * Qz, 1.0, &physvals[0], Qx, &interpolatingNodes[0],
                1, 0.0, &sumFactorization_qr[0], 1);

    // Interpolate in second coordinate direction
    interpolatingNodes = &I[1]->GetPtr()[0];

    Blas::Dgemv('T', Qy, Qz, 1.0, &sumFactorization_qr[0], Qy,
                &interpolatingNodes[0], 1, 0.0, &sumFactorization_r[0], 1);

    // Interpolate in third coordinate direction
    interpolatingNodes = &I[2]->GetPtr()[0];
    value = Vmath::Dot(Qz, interpolatingNodes, 1, &sumFactorization_r[0], 1);

    return value;
}

void StdExpansion3D::v_MultiplyByStdQuadratureMetric(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    int cnt = 0;
    for (int i = 0; i < nquad2; ++i)
    {
        NekDouble w2 = m_weights[2][i];
        for (int j = 0; j < nquad1; ++j)
        {
            NekDouble w1w2 = m_weights[1][j] * w2;
            for (int k = 0; k < nquad0; ++k, ++cnt)
            {
                outarray[cnt] = inarray[cnt] * m_weights[0][k] * w1w2;
            }
        }
    }
}

/**
 * @param   inarray     Input coefficients.
 * @param   output      Output coefficients.
 * @param   mkey        Matrix key
 */
void StdExpansion3D::v_LaplacianMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    if (mkey.GetNVarCoeff() == 0 &&
        !mkey.ConstFactorExists(StdRegions::eFactorCoeffD00) &&
        !mkey.ConstFactorExists(eFactorSVVCutoffRatio))
    {
        // This implementation is only valid when there are no
        // coefficients associated to the Laplacian operator
        int nqtot = GetTotPoints();

        // Allocate temporary storage
        Array<OneD, NekDouble> wsp0(7 * nqtot);
        Array<OneD, NekDouble> wsp1(wsp0 + nqtot);

        if (!(m_base[0]->Collocation() && m_base[1]->Collocation() &&
              m_base[2]->Collocation()))
        {
            // LAPLACIAN MATRIX OPERATION
            // wsp0 = u       = B   * u_hat
            // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
            // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
            BwdTrans(inarray, wsp0);
            LaplacianMatrixOp_MatFree_Kernel(wsp0, outarray, wsp1);
        }
        else
        {
            LaplacianMatrixOp_MatFree_Kernel(inarray, outarray, wsp1);
        }
    }
    else
    {
        StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                            mkey);
    }
}

void StdExpansion3D::v_HelmholtzMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    if (mkey.GetNVarCoeff() == 0 &&
        !mkey.ConstFactorExists(StdRegions::eFactorCoeffD00))
    {
        using std::max;

        int nquad0  = m_base[0]->GetNumPoints();
        int nquad1  = m_base[1]->GetNumPoints();
        int nquad2  = m_base[2]->GetNumPoints();
        int nmodes0 = m_base[0]->GetNumModes();
        int nmodes1 = m_base[1]->GetNumModes();
        int nmodes2 = m_base[2]->GetNumModes();
        int wspsize = max(nquad0 * nmodes2 * (nmodes1 + nquad1),
                          nquad0 * nquad1 * (nquad2 + nmodes0) +
                              nmodes0 * nmodes1 * nquad2);

        NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);

        Array<OneD, NekDouble> wsp0(8 * wspsize);
        Array<OneD, NekDouble> wsp1(wsp0 + 1 * wspsize);
        Array<OneD, NekDouble> wsp2(wsp0 + 2 * wspsize);

        if (!(m_base[0]->Collocation() && m_base[1]->Collocation() &&
              m_base[2]->Collocation()))
        {
            // MASS MATRIX OPERATION
            // The following is being calculated:
            // wsp0     = B   * u_hat = u
            // wsp1     = W   * wsp0
            // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
            BwdTrans(inarray, wsp0);
            IProductWRTBase(wsp0, outarray);
            LaplacianMatrixOp_MatFree_Kernel(wsp0, wsp1, wsp2);
        }
        else
        {
            // specialised implementation for the classical spectral
            // element method
            MultiplyByQuadratureMetric(inarray, outarray);
            LaplacianMatrixOp_MatFree_Kernel(inarray, wsp1, wsp2);
        }

        // outarray = lambda * outarray + wsp1
        //          = (lambda * M + L ) * u_hat
        Vmath::Svtvp(m_ncoeffs, lambda, &outarray[0], 1, &wsp1[0], 1,
                     &outarray[0], 1);
    }
    else
    {
        StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                            mkey);
    }
}

int StdExpansion3D::v_GetNedges(void) const
{
    NEKERROR(ErrorUtil::efatal, "This function is not valid or not defined");
    return 0;
}

int StdExpansion3D::v_GetEdgeNcoeffs([[maybe_unused]] const int i) const
{
    NEKERROR(ErrorUtil::efatal, "This function is not valid or not defined");
    return 0;
}

void StdExpansion3D::v_GetEdgeInteriorToElementMap(
    [[maybe_unused]] const int tid,
    [[maybe_unused]] Array<OneD, unsigned int> &maparray,
    [[maybe_unused]] Array<OneD, int> &signarray,
    [[maybe_unused]] Orientation traceOrient)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion3D::v_GetTraceToElementMap(const int tid,
                                            Array<OneD, unsigned int> &maparray,
                                            Array<OneD, int> &signarray,
                                            Orientation traceOrient, int P,
                                            int Q)
{
    Array<OneD, unsigned int> map1, map2;
    GetTraceCoeffMap(tid, map1);
    GetElmtTraceToTraceMap(tid, map2, signarray, traceOrient, P, Q);

    if (maparray.size() != map2.size())
    {
        maparray = Array<OneD, unsigned int>(map2.size());
    }

    for (int i = 0; i < map2.size(); ++i)
    {
        maparray[i] = map1[map2[i]];
    }
}

/**
 * @brief This method produces a mapping @param idmap which
 * reorientates face data according to the input parameter @param
 * Orient. The sign convention is assumed to take the element local
 * face to a global trace face and this is denoted by the boolean
 * Forwards, i..e globaltrace[i] = localtrace[idmap[i]]. If the
 * boolean is set to Forwards == false then a mapping is produced
 * which maps the gloabl trace back to the local elemental trace such that
 * localtrace[i] = globaltrace[idmap[i]].
 */
void StdExpansion3D::v_ReOrientTracePhysMap(
    const StdRegions::Orientation orient, Array<OneD, int> &idmap,
    const int nq0, const int nq1, bool Forwards)
{
    if (idmap.size() != nq0 * nq1)
    {
        idmap = Array<OneD, int>(nq0 * nq1);
    }

    switch (orient)
    {
        case StdRegions::eDir1FwdDir1_Dir2FwdDir2: // Used in Tri & Quad faces
            // eseentially straight copy
            for (int i = 0; i < nq0 * nq1; ++i)
            {
                idmap[i] = i;
            }
            break;
        case StdRegions::eDir1BwdDir1_Dir2FwdDir2: // Used in Tri & Quad faces
        {
            // Direction A negative and B positive
            for (int j = 0; j < nq1; j++)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    idmap[j * nq0 + i] = nq0 - 1 - i + j * nq0;
                }
            }
        }
        break;
        case StdRegions::eDir1FwdDir1_Dir2BwdDir2:
        {
            // Direction A positive and B negative
            for (int j = 0; j < nq1; j++)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    idmap[j * nq0 + i] = nq0 * (nq1 - 1 - j) + i;
                }
            }
        }
        break;
        case StdRegions::eDir1BwdDir1_Dir2BwdDir2:
        {
            // Direction A negative and B negative
            for (int j = 0; j < nq1; j++)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    idmap[j * nq0 + i] = nq0 * nq1 - 1 - j * nq0 - i;
                }
            }
        }
        break;
        case StdRegions::eDir1FwdDir2_Dir2FwdDir1:
        {
            // Transposed, Direction A and B positive
            if (Forwards)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    for (int j = 0; j < nq1; ++j)
                    {
                        idmap[i * nq1 + j] = i + j * nq0;
                    }
                }
            }
            else // inverse case - different if nq0 != nq1
            {
                for (int j = 0; j < nq1; ++j)
                {
                    for (int i = 0; i < nq0; ++i)
                    {
                        idmap[j * nq0 + i] = i * nq1 + j;
                    }
                }
            }
        }
        break;
        case StdRegions::eDir1FwdDir2_Dir2BwdDir1:
        {
            // Transposed, Direction A positive with mapped direction
            // B and direction B negative with mapped direction A
            if (Forwards)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    for (int j = 0; j < nq1; ++j)
                    {
                        idmap[i * nq1 + j] = i + nq0 * (nq1 - 1) - j * nq0;
                    }
                }
            }
            else
            {
                for (int j = 0; j < nq1; ++j)
                {
                    for (int i = 0; i < nq0; ++i)
                    {
                        idmap[j * nq0 + i] = nq1 - 1 - j + i * nq1;
                    }
                }
            }
        }
        break;
        case StdRegions::eDir1BwdDir2_Dir2FwdDir1:
        {
            // Transposed, Direction A negative with mapped directon B and
            // direction B positive with mapped direction A
            if (Forwards)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    for (int j = 0; j < nq1; ++j)
                    {
                        idmap[i * nq1 + j] = nq0 - 1 - i + j * nq0;
                    }
                }
            }
            else
            {
                for (int j = 0; j < nq1; ++j)
                {
                    for (int i = 0; i < nq0; ++i)
                    {
                        idmap[j * nq0 + i] = nq1 * (nq0 - 1) - i * nq1 + j;
                    }
                }
            }
        }
        break;
        case StdRegions::eDir1BwdDir2_Dir2BwdDir1:
        {
            // Transposed, Direction A and B negative
            if (Forwards)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    for (int j = 0; j < nq1; ++j)
                    {
                        idmap[i * nq1 + j] = nq0 * nq1 - 1 - i - j * nq0;
                    }
                }
            }
            else
            {
                for (int j = 0; j < nq1; ++j)
                {
                    for (int i = 0; i < nq0; ++i)
                    {
                        idmap[j * nq0 + i] = nq0 * nq1 - 1 - j - i * nq1;
                    }
                }
            }
        }
        break;
        default:
            ASSERTL0(false, "Unknow orientation");
            break;
    }
}

LibUtilities::BasisKey EvaluateQuadFaceBasisKey(
    [[maybe_unused]] const int facedir,
    const LibUtilities::BasisSharedPtr &faceDirBasis)
{
    auto faceDirBasisType = faceDirBasis->GetBasisType();
    auto pointsType       = faceDirBasis->GetPointsType();
    auto nummodes         = faceDirBasis->GetNumModes();
    auto numpoints        = faceDirBasis->GetNumPoints();

    switch (faceDirBasisType)
    {
        case LibUtilities::eModified_A:
        case LibUtilities::eModified_B:
        case LibUtilities::eModified_C:
        {
            LibUtilities::PointsType pType;
            switch (pointsType)
            {
                case LibUtilities::eGaussRadauMLegendre:
                case LibUtilities::eGaussRadauMAlpha2Beta0:
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                {
                    numpoints = numpoints + 1;
                    pType     = LibUtilities::eGaussLobattoLegendre;
                }
                break;
                case LibUtilities::eGaussLegendreWithM:
                {
                    numpoints = numpoints + 1;
                    pType     = LibUtilities::eGaussLegendreWithMP;
                }
                break;
                default: // do not change points
                {
                    pType = faceDirBasis->GetPointsType();
                }
            }
            const LibUtilities::PointsKey pkey(numpoints, pType);
            return LibUtilities::BasisKey(LibUtilities::eModified_A, nummodes,
                                          pkey);
        }
        case LibUtilities::eGLL_Lagrange:
        {
            const LibUtilities::PointsKey pkey(
                numpoints, LibUtilities::eGaussLobattoLegendre);
            return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange, nummodes,
                                          pkey);
        }
        case LibUtilities::eOrtho_A:
        case LibUtilities::eOrtho_B:
        case LibUtilities::eOrtho_C:
        {
            LibUtilities::PointsType pType;
            switch (pointsType)
            {
                case LibUtilities::eGaussRadauMLegendre:
                case LibUtilities::eGaussRadauMAlpha2Beta0:
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                {
                    numpoints = numpoints + 1;
                    pType     = LibUtilities::eGaussLobattoLegendre;
                }
                break;
                case LibUtilities::eGaussLegendreWithM:
                {
                    numpoints = numpoints + 1;
                    pType     = LibUtilities::eGaussLegendreWithMP;
                }
                break;
                default: // do not change points
                {
                    pType = faceDirBasis->GetPointsType();
                    break;
                }
            }
            const LibUtilities::PointsKey pkey(numpoints, pType);
            return LibUtilities::BasisKey(LibUtilities::eOrtho_A, nummodes,
                                          pkey);
        }
        default:
        {
            NEKERROR(ErrorUtil::efatal, "expansion type unknown");
            break;
        }
    }

    // Keep things happy by returning a value.
    return LibUtilities::NullBasisKey;
}

LibUtilities::BasisKey EvaluateTriFaceBasisKey(
    const int facedir, const LibUtilities::BasisSharedPtr &faceDirBasis,
    bool UseGLL)
{
    auto faceDirBasisType = faceDirBasis->GetBasisType();
    auto pointsType       = faceDirBasis->GetPointsType();
    auto nummodes         = faceDirBasis->GetNumModes();
    auto numpoints        = faceDirBasis->GetNumPoints();

    switch (faceDirBasisType)
    {
        case LibUtilities::eModified_A:
        case LibUtilities::eModified_B:
        case LibUtilities::eModified_C:
        case LibUtilities::eModifiedPyr_C:
        {
            LibUtilities::BasisType bType = LibUtilities::eNoBasisType;
            LibUtilities::PointsKey pkey  = LibUtilities::NullPointsKey;
            switch (facedir) // determine the basis type
            {
                case 0:
                {
                    bType = LibUtilities::eModified_A;

                    switch (pointsType) // determine the points type
                    {
                        case LibUtilities::eGaussRadauMLegendre:
                            // case LibUtilities::eGaussRadauMAlpha2Beta0: don't
                            // this this is posssible
                        case LibUtilities::eGaussRadauMAlpha1Beta0:
                        {
                            pkey = LibUtilities::PointsKey(
                                numpoints + 1,
                                LibUtilities::eGaussLobattoLegendre);
                        }
                        break;
                        case LibUtilities::eGaussLegendreWithM:
                        {
                            pkey = LibUtilities::PointsKey(
                                numpoints + 1,
                                LibUtilities::eGaussLegendreWithMP);
                        }
                        break;
                        default: // For other points type, just return the
                                 // points key
                        {
                            pkey = faceDirBasis->GetPointsKey();
                        }
                    }
                }
                break;
                case 1: // this never appears together with Modified_A
                {
                    bType = LibUtilities::eModified_B;

                    switch (pointsType) // determine the points type
                    {
                        case LibUtilities::eGaussRadauMLegendre:
                        case LibUtilities::eGaussRadauMAlpha2Beta0:
                        case LibUtilities::eGaussRadauMAlpha1Beta0:
                        {
                            if (UseGLL) //  force to use GLL
                            {
                                pkey = LibUtilities::PointsKey(
                                    numpoints + 1,
                                    LibUtilities::eGaussLobattoLegendre);
                            }
                            else
                            {
                                pkey = LibUtilities::PointsKey(
                                    numpoints,
                                    LibUtilities::eGaussRadauMAlpha1Beta0);
                            }
                        }
                        break;
                        default: // For other points type, just return the
                                 // points key
                        {
                            pkey = faceDirBasis->GetPointsKey();
                        }
                    }
                }
                break;
                default:
                {
                    NEKERROR(ErrorUtil::efatal, "invalid value to flag");
                    break;
                }
            }
            return LibUtilities::BasisKey(bType, nummodes, pkey);
        }

        case LibUtilities::eGLL_Lagrange:
        {
            switch (facedir)
            {
                case 0:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints, LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eOrtho_A,
                                                  nummodes, pkey);
                }
                break;
                case 1:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints, LibUtilities::eGaussRadauMAlpha1Beta0);
                    return LibUtilities::BasisKey(LibUtilities::eOrtho_B,
                                                  nummodes, pkey);
                }
                break;
                default:
                {
                    NEKERROR(ErrorUtil::efatal, "invalid value to flag");
                    break;
                }
            }
            break;
        }

        case LibUtilities::eOrtho_A:
        case LibUtilities::eOrtho_B:
        case LibUtilities::eOrtho_C:
        case LibUtilities::eOrthoPyr_C:
        {
            LibUtilities::BasisType bType = LibUtilities::eNoBasisType;
            LibUtilities::PointsKey pkey  = LibUtilities::NullPointsKey;
            switch (facedir) // determine the basis type
            {
                case 0:
                {
                    bType = LibUtilities::eOrtho_A;

                    switch (pointsType) // determine the points type
                    {
                        case LibUtilities::eGaussRadauMLegendre:
                        case LibUtilities::eGaussRadauMAlpha2Beta0:
                        case LibUtilities::eGaussRadauMAlpha1Beta0:
                        {
                            pkey = LibUtilities::PointsKey(
                                numpoints + 1,
                                LibUtilities::eGaussLobattoLegendre);
                        }
                        break;
                        case LibUtilities::eGaussLegendreWithM:
                        {
                            pkey = LibUtilities::PointsKey(
                                numpoints + 1,
                                LibUtilities::eGaussLegendreWithMP);
                        }
                        break;
                        default: // For other points type, just return the
                                 // points key
                        {
                            pkey = faceDirBasis->GetPointsKey();
                        }
                    }
                }
                break;
                case 1: // this never appears together with Ortho_A
                {
                    bType = LibUtilities::eOrtho_B;

                    switch (pointsType) // determine the points type
                    {
                        case LibUtilities::eGaussRadauMLegendre:
                        case LibUtilities::eGaussRadauMAlpha2Beta0:
                        {
                            pkey = LibUtilities::PointsKey(
                                numpoints,
                                LibUtilities::eGaussRadauMAlpha1Beta0);
                        }
                        break;
                        default: // For other points type, just return the
                                 // points key
                        {
                            pkey = faceDirBasis->GetPointsKey();
                        }
                    }
                }
                break;
                default:
                {
                    NEKERROR(ErrorUtil::efatal, "invalid value to flag");
                    break;
                }
            }
            return LibUtilities::BasisKey(bType, nummodes, pkey);
        }
        default:
        {
            NEKERROR(ErrorUtil::efatal, "expansion type unknown");
            break;
        }
    }

    // Keep things happy by returning a value.
    return LibUtilities::NullBasisKey;
}

void StdExpansion3D::v_PhysInterp(std::shared_ptr<StdExpansion> fromExp,
                                  const Array<OneD, const NekDouble> &fromData,
                                  Array<OneD, NekDouble> &toData,
                                  [[maybe_unused]] bool Transpose)
{

    LibUtilities::Interp3D(fromExp->GetBasis(0)->GetPointsKey(),
                           fromExp->GetBasis(1)->GetPointsKey(),
                           fromExp->GetBasis(2)->GetPointsKey(), fromData,
                           m_base[0]->GetPointsKey(), m_base[1]->GetPointsKey(),
                           m_base[2]->GetPointsKey(), toData);
}

} // namespace Nektar::StdRegions
