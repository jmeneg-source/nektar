///////////////////////////////////////////////////////////////////////////////
//
// File: StdExpansion2D.cpp
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
// which are common to 2D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion2D.h>

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

StdExpansion2D::StdExpansion2D(
    [[maybe_unused]] int numcoeffs,
    [[maybe_unused]] const LibUtilities::BasisKey &Ba,
    [[maybe_unused]] const LibUtilities::BasisKey &Bb)
{
}

//----------------------------
// Differentiation Methods
//----------------------------
/**
 *   Calculate the derivative along the tenosr directions. This function was
 *  originally in StdEpxansion2D but due to the boost_pp switch statement is
 *  currently shape dependent
 */
void StdExpansion2D::PhysTensorDeriv(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray_d0, Array<OneD, NekDouble> &outarray_d1)
{
    int nquad0          = m_base[0]->GetNumPoints();
    int nquad1          = m_base[1]->GetNumPoints();
    bool Deriv0         = (outarray_d0.size() > 0);
    bool Deriv1         = (outarray_d1.size() > 0);
    const NekDouble *D0 = m_base[0]->GetD()->GetRawPtr();
    const NekDouble *D1 = m_base[1]->GetD()->GetRawPtr();

    Array<OneD, const NekDouble> intmp;
    // copy inarray data if inarray and outarray are the same.
    if ((inarray.data() == outarray_d0.data()) ||
        (inarray.data() == outarray_d1.data()))
    {
        Array<OneD, NekDouble> wsp(nquad0 * nquad1);
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
    PhysDerivTensor2DKernel(nquad0, nquad1, (const vec_t *)intmp.data(),       \
                            (const vec_t *)D0, (const vec_t *)D1,              \
                            (vec_t *)outarray_d0.data(),                       \
                            (vec_t *)outarray_d1.data(), Deriv0, Deriv1)

    // Loop case over quarature points
#undef PHYSDERIV_Q
#define PHYSDERIV_Q(r, i)                                                      \
    case NQ1(i):                                                               \
        PhysDerivTensor2DKernel(NQ1(i), NQ1(i), (const vec_t *)intmp.data(),   \
                                (const vec_t *)D0, (const vec_t *)D1,          \
                                (vec_t *)outarray_d0.data(),                   \
                                (vec_t *)outarray_d1.data(), Deriv0, Deriv1);  \
        break;

    // templated cases on  standard quadrature
    // usage where quad order goes from SMIN to SMAX
    if (nquad0 == nquad1)
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

void StdExpansion2D::v_PhysDeriv(const int dir,
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
        default:
        {
            ASSERTL1(false, "input dir is out of range");
            break;
        }
    }
}

NekDouble StdExpansion2D::v_StdPhysEvaluate(
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &physvals)
{
    ASSERTL2(coords[0] > -1 - NekConstants::kNekZeroTol, "coord[0] < -1");
    ASSERTL2(coords[0] < 1 + NekConstants::kNekZeroTol, "coord[0] >  1");
    ASSERTL2(coords[1] > -1 - NekConstants::kNekZeroTol, "coord[1] < -1");
    ASSERTL2(coords[1] < 1 + NekConstants::kNekZeroTol, "coord[1] >  1");

    Array<OneD, NekDouble> coll(2);
    LocCoordToLocCollapsed(coords, coll);

    const int nq0 = m_base[0]->GetNumPoints();
    const int nq1 = m_base[1]->GetNumPoints();

    Array<OneD, NekDouble> wsp(nq1);
    for (int i = 0; i < nq1; ++i)
    {
        wsp[i] = StdExpansion::BaryEvaluate<0>(coll[0], &physvals[0] + i * nq0);
    }

    return StdExpansion::BaryEvaluate<1>(coll[1], &wsp[0]);
}

NekDouble StdExpansion2D::v_PhysEvaluateInterp(
    const Array<OneD, DNekMatSharedPtr> &I,
    const Array<OneD, const NekDouble> &physvals)
{
    NekDouble val;
    int i;
    int nq0 = m_base[0]->GetNumPoints();
    int nq1 = m_base[1]->GetNumPoints();
    Array<OneD, NekDouble> wsp1(nq1);

    // interpolate first coordinate direction
    for (i = 0; i < nq1; ++i)
    {
        wsp1[i] =
            Vmath::Dot(nq0, &(I[0]->GetPtr())[0], 1, &physvals[i * nq0], 1);
    }

    // interpolate in second coordinate direction
    val = Vmath::Dot(nq1, I[1]->GetPtr(), 1, wsp1, 1);

    return val;
}

/** \brief Calculate the inner product of inarray with respect to
 *  the basis B=base0*base1 and put into outarray
 *
 *  \f$
 *  \begin{array}{rcl}
 *  I_{pq} = (\phi_p \phi_q, u) & = & \sum_{i=0}^{nq_0}
 *  \sum_{j=0}^{nq_1}
 *  \phi_p(\xi_{0,i}) \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i}
 *  \xi_{1,j}) \\
 *  & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
 *  \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}
 *  \end{array}
 *  \f$
 *
 *  where
 *
 *  \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
 *
 *  which can be implemented as
 *
 *  \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j})
 *  \tilde{u}_{i,j} = {\bf B_1 U}  \f$
 *  \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} =
 *  {\bf B_0 F}  \f$
 *
 * This is a wrapper function around \a IProductWRTBaseKernel()
 */
void StdExpansion2D::v_IProductWRTBase(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    const bool CollDir0 = m_base[0]->Collocation();
    const bool CollDir1 = m_base[1]->Collocation();

    if (CollDir0 && CollDir1)
    {
        v_MultiplyByStdQuadratureMetric(inarray, outarray);
    }
    else
    {
        const Array<OneD, const NekDouble> one(1, 1.0);
        v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                                inarray, outarray, one, false, CollDir0,
                                CollDir1);
    }
}

void StdExpansion2D::IProductWRTBaseKernel(
    const Array<OneD, const NekDouble> &base0,
    const Array<OneD, const NekDouble> &base1,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const Array<OneD, NekDouble> &jac,
    const bool Deformed, const bool CollDir0, const bool CollDir1)
{
    v_IProductWRTBaseKernel(base0, base1, inarray, outarray, jac, Deformed,
                            CollDir0, CollDir1);
}

void StdExpansion2D::v_GenStdMatBwdDeriv(const int dir, DNekMatSharedPtr &mat)
{
    ASSERTL1((dir == 0) || (dir == 1), "Invalid direction.");

    const int nq0 = m_base[0]->GetNumPoints();
    const int nq1 = m_base[1]->GetNumPoints();
    const int nq  = nq0 * nq1;

    Array<OneD, NekDouble> in(nq, 0.0);
    Array<OneD, NekDouble> out(m_ncoeffs);
    Array<OneD, NekDouble> one(1, 1.0);

    for (int i = 0; i < nq; i++)
    {
        int l = i % nq0;
        int m = (i / nq0);

        // initialise with inverse of weights t
        in[i] = 1.0 / (m_weights[0][l] * m_weights[1][m]);

        // do standard iproduct
        if (dir == 0)
        {
            v_IProductWRTBaseKernel(m_base[0]->GetDbdata(),
                                    m_base[1]->GetBdata(), in, out, one, false,
                                    false, m_base[1]->Collocation());
        }
        else if (dir == 1)
        {
            v_IProductWRTBaseKernel(m_base[0]->GetBdata(),
                                    m_base[1]->GetDbdata(), in, out, one, false,
                                    m_base[0]->Collocation(), false);
        }
        in[i] = 0.0;

        for (int j = 0; j < m_ncoeffs; j++)
        {
            (*mat)(j, i) = out[j];
        }
    }
}

void StdExpansion2D::v_MultiplyByStdQuadratureMetric(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();

    int cnt = 0;
    for (int i = 0; i < nquad1; ++i)
    {
        for (int j = 0; j < nquad0; ++j, ++cnt)
        {
            outarray[cnt] = inarray[cnt] * m_weights[0][j] * m_weights[1][i];
        }
    }
}

void StdExpansion2D::v_LaplacianMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    if (mkey.GetNVarCoeff() == 0 &&
        !mkey.ConstFactorExists(StdRegions::eFactorCoeffD00) &&
        !mkey.ConstFactorExists(StdRegions::eFactorSVVCutoffRatio))
    {
        using std::max;

        // This implementation is only valid when there are no
        // coefficients associated to the Laplacian operator
        int nquad0  = m_base[0]->GetNumPoints();
        int nquad1  = m_base[1]->GetNumPoints();
        int nqtot   = nquad0 * nquad1;
        int nmodes0 = m_base[0]->GetNumModes();
        int nmodes1 = m_base[1]->GetNumModes();
        int wspsize =
            max(max(max(nqtot, m_ncoeffs), nquad1 * nmodes0), nquad0 * nmodes1);

        // Allocate temporary storage
        Array<OneD, NekDouble> wsp0(4 * wspsize);    // size wspsize
        Array<OneD, NekDouble> wsp1(wsp0 + wspsize); // size 3*wspsize

        if (!(m_base[0]->Collocation() && m_base[1]->Collocation()))
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

void StdExpansion2D::v_HelmholtzMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    if (mkey.GetNVarCoeff() == 0 &&
        !mkey.ConstFactorExists(StdRegions::eFactorCoeffD00) &&
        !mkey.ConstFactorExists(StdRegions::eFactorSVVCutoffRatio))
    {
        using std::max;

        int nquad0  = m_base[0]->GetNumPoints();
        int nquad1  = m_base[1]->GetNumPoints();
        int nqtot   = nquad0 * nquad1;
        int nmodes0 = m_base[0]->GetNumModes();
        int nmodes1 = m_base[1]->GetNumModes();
        int wspsize =
            max(max(max(nqtot, m_ncoeffs), nquad1 * nmodes0), nquad0 * nmodes1);
        NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);

        // Allocate temporary storage
        Array<OneD, NekDouble> wsp0(5 * wspsize);        // size wspsize
        Array<OneD, NekDouble> wsp1(wsp0 + wspsize);     // size wspsize
        Array<OneD, NekDouble> wsp2(wsp0 + 2 * wspsize); // size 3*wspsize

        if (!(m_base[0]->Collocation() && m_base[1]->Collocation()))
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

void StdExpansion2D::v_GetTraceCoeffMap(
    [[maybe_unused]] const unsigned int traceid,
    [[maybe_unused]] Array<OneD, unsigned int> &maparray)
{
    ASSERTL0(false,
             "This method must be defined at the individual shape level");
}

/** \brief Determine the mapping to re-orientate the
    coefficients along the element trace (assumed to align
    with the standard element) into the orientation of the
    local trace given by edgeOrient.
 */
void StdExpansion2D::v_GetElmtTraceToTraceMap(
    const unsigned int eid, Array<OneD, unsigned int> &maparray,
    Array<OneD, int> &signarray, Orientation edgeOrient, int P,
    [[maybe_unused]] int Q)
{
    unsigned int i;

    int dir;
    // determine basis direction for edge.
    if (DetShapeType() == LibUtilities::eTriangle)
    {
        dir = (eid == 0) ? 0 : 1;
    }
    else
    {
        dir = eid % 2;
    }

    int numModes = m_base[dir]->GetNumModes();

    // P is the desired length of the map
    P = (P == -1) ? numModes : P;

    // decalare maparray
    if (maparray.size() != P)
    {
        maparray = Array<OneD, unsigned int>(P);
    }

    // fill default mapping as increasing index
    for (i = 0; i < P; ++i)
    {
        maparray[i] = i;
    }

    if (signarray.size() != P)
    {
        signarray = Array<OneD, int>(P, 1);
    }
    else
    {
        std::fill(signarray.data(), signarray.data() + P, 1);
    }

    // Zero signmap and set maparray to zero if
    // elemental modes are not as large as trace modes
    for (i = numModes; i < P; ++i)
    {
        signarray[i] = 0.0;
        maparray[i]  = maparray[0];
    }

    if (edgeOrient == eBackwards)
    {
        const LibUtilities::BasisType bType = GetBasisType(dir);

        if ((bType == LibUtilities::eModified_A) ||
            (bType == LibUtilities::eModified_B))
        {
            std::swap(maparray[0], maparray[1]);

            for (i = 3; i < std::min(P, numModes); i += 2)
            {
                signarray[i] *= -1;
            }
        }
        else if (bType == LibUtilities::eGLL_Lagrange ||
                 bType == LibUtilities::eGauss_Lagrange)
        {
            ASSERTL1(P == numModes, "Different trace space edge dimension "
                                    "and element edge dimension not currently "
                                    "possible for GLL-Lagrange bases");

            std::reverse(maparray.data(), maparray.data() + P);
        }
        else
        {
            ASSERTL0(false, "Mapping not defined for this type of basis");
        }
    }
}

void StdExpansion2D::v_GetTraceToElementMap(const int eid,
                                            Array<OneD, unsigned int> &maparray,
                                            Array<OneD, int> &signarray,
                                            Orientation edgeOrient, int P,
                                            int Q)
{
    Array<OneD, unsigned int> map1, map2;
    v_GetTraceCoeffMap(eid, map1);
    v_GetElmtTraceToTraceMap(eid, map2, signarray, edgeOrient, P, Q);

    if (maparray.size() != map2.size())
    {
        maparray = Array<OneD, unsigned int>(map2.size());
    }

    for (int i = 0; i < map2.size(); ++i)
    {
        maparray[i] = map1[map2[i]];
    }
}

void StdExpansion2D::v_PhysInterp(std::shared_ptr<StdExpansion> fromExp,
                                  const Array<OneD, const NekDouble> &fromData,
                                  Array<OneD, NekDouble> &toData,
                                  bool Transpose)
{
    if (Transpose)
    {
        LibUtilities::Interp2D(fromExp->GetBasis(0)->GetPointsKey(),
                               fromExp->GetBasis(1)->GetPointsKey(), fromData,
                               m_base[1]->GetPointsKey(),
                               m_base[0]->GetPointsKey(), toData);
    }
    else
    {
        LibUtilities::Interp2D(fromExp->GetBasis(0)->GetPointsKey(),
                               fromExp->GetBasis(1)->GetPointsKey(), fromData,
                               m_base[0]->GetPointsKey(),
                               m_base[1]->GetPointsKey(), toData);
    }
}

void StdExpansion2D::v_ReOrientTracePhysMap(
    const StdRegions::Orientation orient, Array<OneD, int> &idmap,
    const int nq0, [[maybe_unused]] const int nq1,
    [[maybe_unused]] bool Forwards)
{
    if (idmap.size() != nq0)
    {
        idmap = Array<OneD, int>(nq0);
    }
    switch (orient)
    {
        case StdRegions::eForwards:
            // Fwd
            for (int i = 0; i < nq0; ++i)
            {
                idmap[i] = i;
            }
            break;
        case StdRegions::eBackwards:
        {
            // Bwd
            for (int i = 0; i < nq0; ++i)
            {
                idmap[i] = nq0 - 1 - i;
            }
        }
        break;
        default:
            ASSERTL0(false, "Unknown orientation");
            break;
    }
}

} // namespace Nektar::StdRegions
