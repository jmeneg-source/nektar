///////////////////////////////////////////////////////////////////////////////
//
// File: TetExp.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LocalRegions/TetExp.h>
#include <SpatialDomains/SegGeom.h>

using namespace std;

namespace Nektar::LocalRegions
{
/**
 * @class TetExp
 * Defines a Tetrahedral local expansion.
 */

/**
 * \brief Constructor using BasisKey class for quadrature points and
 * order definition
 *
 * @param   Ba          Basis key for first coordinate.
 * @param   Bb          Basis key for second coordinate.
 * @param   Bc          Basis key for third coordinate.
 */
TetExp::TetExp(const LibUtilities::BasisKey &Ba,
               const LibUtilities::BasisKey &Bb,
               const LibUtilities::BasisKey &Bc,
               SpatialDomains::Geometry3D *geom)
    : StdExpansion(LibUtilities::StdTetData::getNumberOfCoefficients(
                       Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                   3, Ba, Bb, Bc),
      StdExpansion3D(LibUtilities::StdTetData::getNumberOfCoefficients(
                         Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                     Ba, Bb, Bc),
      StdTetExp(Ba, Bb, Bc), Expansion(geom), Expansion3D(geom),
      m_matrixManager(
          std::bind(&Expansion3D::CreateMatrix, this, std::placeholders::_1)),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1))
{
}

/**
 * \brief Copy Constructor
 */
TetExp::TetExp(const TetExp &T)
    : StdRegions::StdExpansion(T), StdRegions::StdExpansion3D(T),
      StdRegions::StdTetExp(T), Expansion(T), Expansion3D(T),
      m_matrixManager(T.m_matrixManager),
      m_staticCondMatrixManager(T.m_staticCondMatrixManager)
{
}

//-----------------------------
// Inner product functions
//-----------------------------
/**
 * @brief Calculates the inner product \f$ I_{pqr} = (u,
 * \partial_{x_i} \phi_{pqr}) \f$.
 *
 * The derivative of the basis functions is performed using the chain
 * rule in order to incorporate the geometric factors. Assuming that
 * the basis functions are a tensor product
 * \f$\phi_{pqr}(\eta_1,\eta_2,\eta_3) =
 * \phi_1(\eta_1)\phi_2(\eta_2)\phi_3(\eta_3)\f$, this yields the
 * result
 *
 * \f[
 * I_{pqr} = \sum_{j=1}^3 \left(u, \frac{\partial u}{\partial \eta_j}
 * \frac{\partial \eta_j}{\partial x_i}\right)
 * \f]
 *
 * In the prismatic element, we must also incorporate a second set of
 * geometric factors which incorporate the collapsed co-ordinate
 * system, so that
 *
 * \f[ \frac{\partial\eta_j}{\partial x_i} = \sum_{k=1}^3
 * \frac{\partial\eta_j}{\partial\xi_k}\frac{\partial\xi_k}{\partial
 * x_i} \f]
 *
 * These derivatives can be found on p152 of Sherwin & Karniadakis.
 *
 * @param dir       Direction in which to take the derivative.
 * @param inarray   The function \f$ u \f$.
 * @param outarray  Value of the inner product.
 */
void TetExp::v_IProductWRTDerivBase(const int dir,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    const int nquad0 = m_base[0]->GetNumPoints();
    const int nquad1 = m_base[1]->GetNumPoints();
    const int nquad2 = m_base[2]->GetNumPoints();
    const int nqtot  = nquad0 * nquad1 * nquad2;

    Array<OneD, NekDouble> tmp2(nqtot);
    Array<OneD, NekDouble> tmp3(nqtot);
    Array<OneD, NekDouble> tmp4(nqtot);
    Array<OneD, NekDouble> tmp6(m_ncoeffs);

    Array<OneD, Array<OneD, NekDouble>> tmp2D{3};
    tmp2D[0] = tmp2;
    tmp2D[1] = tmp3;
    tmp2D[2] = tmp4;

    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    bool Deformed = (m_geomFactors->GetGtype() == SpatialDomains::eDeformed);

    v_AlignVectorToCollapsedDir(dir, inarray, tmp2D);

    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                            m_base[2]->GetBdata(), tmp2, outarray, jac,
                            Deformed);

    v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                            m_base[2]->GetBdata(), tmp3, tmp6, jac, Deformed);
    Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);

    v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                            m_base[2]->GetDbdata(), tmp4, tmp6, jac, Deformed);
    Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);
}

void TetExp::v_AlignVectorToCollapsedDir(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int i, j;

    const int nquad0 = m_base[0]->GetNumPoints();
    const int nquad1 = m_base[1]->GetNumPoints();
    const int nquad2 = m_base[2]->GetNumPoints();
    const int nqtot  = nquad0 * nquad1 * nquad2;

    const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
    const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();
    const Array<OneD, const NekDouble> &z2 = m_base[2]->GetZ();

    Array<OneD, NekDouble> tmp2(nqtot);
    Array<OneD, NekDouble> tmp3(nqtot);

    const Array<TwoD, const NekDouble> &df = m_geomFactors->GetDerivFactors();

    if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vmul(nqtot, &df[3 * dir][0], 1, inarray.data(), 1, tmp2.data(),
                    1);
        Vmath::Vmul(nqtot, &df[3 * dir + 1][0], 1, inarray.data(), 1,
                    tmp3.data(), 1);
        Vmath::Vmul(nqtot, &df[3 * dir + 2][0], 1, inarray.data(), 1,
                    outarray[2].data(), 1);
    }
    else
    {
        Vmath::Smul(nqtot, df[3 * dir][0], inarray.data(), 1, tmp2.data(), 1);
        Vmath::Smul(nqtot, df[3 * dir + 1][0], inarray.data(), 1, tmp3.data(),
                    1);
        Vmath::Smul(nqtot, df[3 * dir + 2][0], inarray.data(), 1,
                    outarray[2].data(), 1);
    }

    NekDouble g0, g1, g1a, g2, g3;
    int k, cnt;

    for (cnt = 0, k = 0; k < nquad2; ++k)
    {
        g2 = 2.0 / (1.0 - z2[k]);
        for (j = 0; j < nquad1; ++j)
        {
            g1 = g2 / (1.0 - z1[j]);
            g0 = 2.0 * g1;
            g3 = (1.0 + z1[j]) * g2 * 0.5;

            for (i = 0; i < nquad0; ++i, ++cnt)
            {
                g1a = g1 * (1 + z0[i]);

                outarray[0][cnt] =
                    g0 * tmp2[cnt] + g1a * (tmp3[cnt] + outarray[2][cnt]);

                outarray[1][cnt] = g2 * tmp3[cnt] + g3 * outarray[2][cnt];
            }
        }
    }
}

//-----------------------------
// Evaluation functions
//-----------------------------
NekDouble TetExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    Array<OneD, NekDouble> Lcoord(3);
    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);
    return StdTetExp::v_PhysEvalFirstDeriv(Lcoord, inarray, firstOrderDerivs);
}

/**
 * \brief Get the coordinates "coords" at the local coordinates "Lcoords"
 */
void TetExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                        Array<OneD, NekDouble> &coords)
{
    int i;

    ASSERTL1(Lcoords[0] <= -1.0 && Lcoords[0] >= 1.0 && Lcoords[1] <= -1.0 &&
                 Lcoords[1] >= 1.0 && Lcoords[2] <= -1.0 && Lcoords[2] >= 1.0,
             "Local coordinates are not in region [-1,1]");

    // m_geom->FillGeom(); // TODO: implement FillGeom()

    for (i = 0; i < m_geom->GetCoordim(); ++i)
    {
        coords[i] = m_geom->GetCoord(i, Lcoords);
    }
}

void TetExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                         Array<OneD, NekDouble> &coords_1,
                         Array<OneD, NekDouble> &coords_2)
{
    Expansion::v_GetCoords(coords_0, coords_1, coords_2);
}

//-----------------------------
// Helper functions
//-----------------------------

StdRegions::StdExpansionSharedPtr TetExp::v_GetStdExp(void) const
{
    return MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey(), m_base[1]->GetBasisKey(),
        m_base[2]->GetBasisKey());
}

StdRegions::StdExpansionSharedPtr TetExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());
    LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(), 2,
                                 m_base[1]->GetPointsKey());
    LibUtilities::BasisKey bkey2(m_base[2]->GetBasisType(), 2,
                                 m_base[2]->GetPointsKey());

    return MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(bkey0, bkey1,
                                                                   bkey2);
}

void TetExp::v_ExtractDataToCoeffs(
    const NekDouble *data, const std::vector<unsigned int> &nummodes,
    const int mode_offset, NekDouble *coeffs,
    [[maybe_unused]] std::vector<LibUtilities::BasisType> &fromType)
{
    int data_order0 = nummodes[mode_offset];
    int fillorder0  = min(m_base[0]->GetNumModes(), data_order0);
    int data_order1 = nummodes[mode_offset + 1];
    int order1      = m_base[1]->GetNumModes();
    int fillorder1  = min(order1, data_order1);
    int data_order2 = nummodes[mode_offset + 2];
    int order2      = m_base[2]->GetNumModes();
    int fillorder2  = min(order2, data_order2);

    switch (m_base[0]->GetBasisType())
    {
        case LibUtilities::eModified_A:
        {
            int i, j;
            int cnt  = 0;
            int cnt1 = 0;

            ASSERTL1(m_base[1]->GetBasisType() == LibUtilities::eModified_B,
                     "Extraction routine not set up for this basis");
            ASSERTL1(m_base[2]->GetBasisType() == LibUtilities::eModified_C,
                     "Extraction routine not set up for this basis");

            Vmath::Zero(m_ncoeffs, coeffs, 1);
            for (j = 0; j < fillorder0; ++j)
            {
                for (i = 0; i < fillorder1 - j; ++i)
                {
                    Vmath::Vcopy(fillorder2 - j - i, &data[cnt], 1,
                                 &coeffs[cnt1], 1);
                    cnt += data_order2 - j - i;
                    cnt1 += order2 - j - i;
                }

                // count out data for j iteration
                for (i = fillorder1 - j; i < data_order1 - j; ++i)
                {
                    cnt += data_order2 - j - i;
                }

                for (i = fillorder1 - j; i < order1 - j; ++i)
                {
                    cnt1 += order2 - j - i;
                }
            }
        }
        break;
        default:
            ASSERTL0(false, "basis is either not set up or not "
                            "hierarchicial");
    }
}

/**
 * \brief Returns the physical values at the quadrature points of a face
 */
void TetExp::v_GetTracePhysMap(const int face, Array<OneD, int> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    int nq0 = 0;
    int nq1 = 0;

    // get forward aligned faces.
    switch (face)
    {
        case 0:
        {
            nq0 = nquad0;
            nq1 = nquad1;
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            for (int i = 0; i < nquad0 * nquad1; ++i)
            {
                outarray[i] = i;
            }

            break;
        }
        case 1:
        {
            nq0 = nquad0;
            nq1 = nquad2;
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            // Direction A and B positive
            for (int k = 0; k < nquad2; k++)
            {
                for (int i = 0; i < nquad0; ++i)
                {
                    outarray[k * nquad0 + i] = (nquad0 * nquad1 * k) + i;
                }
            }
            break;
        }
        case 2:
        {
            nq0 = nquad1;
            nq1 = nquad2;
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            // Directions A and B positive
            for (int j = 0; j < nquad1 * nquad2; ++j)
            {
                outarray[j] = nquad0 - 1 + j * nquad0;
            }
            break;
        }
        case 3:
        {
            nq0 = nquad1;
            nq1 = nquad2;
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            // Directions A and B positive
            for (int j = 0; j < nquad1 * nquad2; ++j)
            {
                outarray[j] = j * nquad0;
            }
        }
        break;
        default:
            ASSERTL0(false, "face value (> 3) is out of range");
            break;
    }
}

/**
 * \brief Compute the normal of a triangular face
 */
void TetExp::v_ComputeTraceNormal(const int face)
{
    int i;
    LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
    for (int i = 0; i < ptsKeys.size(); ++i)
    {
        // Need at least 2 points for computing normals
        if (ptsKeys[i].GetNumPoints() == 1)
        {
            LibUtilities::PointsKey pKey(2, ptsKeys[i].GetPointsType());
            ptsKeys[i] = pKey;
        }
    }

    SpatialDomains::GeomType type = m_geomFactors->GetGtype();
    const Array<TwoD, const NekDouble> &df =
        m_geomFactors->ComputeDerivFactors(ptsKeys);
    const Array<OneD, const NekDouble> &jac =
        m_geomFactors->ComputeJac(ptsKeys);

    LibUtilities::BasisKey tobasis0 = GetTraceBasisKey(face, 0);
    LibUtilities::BasisKey tobasis1 = GetTraceBasisKey(face, 1);

    // number of face quadrature points
    int nq_face = tobasis0.GetNumPoints() * tobasis1.GetNumPoints();

    int vCoordDim = GetCoordim();

    m_traceNormals[face] = Array<OneD, Array<OneD, NekDouble>>(vCoordDim);
    Array<OneD, Array<OneD, NekDouble>> &normal = m_traceNormals[face];
    for (i = 0; i < vCoordDim; ++i)
    {
        normal[i] = Array<OneD, NekDouble>(nq_face);
    }

    size_t nqb                     = nq_face;
    size_t nbnd                    = face;
    m_elmtBndNormDirElmtLen[nbnd]  = Array<OneD, NekDouble>{nqb, 0.0};
    Array<OneD, NekDouble> &length = m_elmtBndNormDirElmtLen[nbnd];

    // Regular geometry case
    if (type == SpatialDomains::eRegular ||
        type == SpatialDomains::eMovingRegular)
    {
        NekDouble fac;

        // Set up normals
        switch (face)
        {
            case 0:
            {
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = -df[3 * i + 2][0];
                }

                break;
            }
            case 1:
            {
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = -df[3 * i + 1][0];
                }

                break;
            }
            case 2:
            {
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] =
                        df[3 * i][0] + df[3 * i + 1][0] + df[3 * i + 2][0];
                }

                break;
            }
            case 3:
            {
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = -df[3 * i][0];
                }
                break;
            }
            default:
                ASSERTL0(false, "face is out of range (edge < 3)");
        }

        // normalise
        fac = 0.0;
        for (i = 0; i < vCoordDim; ++i)
        {
            fac += normal[i][0] * normal[i][0];
        }
        fac = 1.0 / sqrt(fac);
        Vmath::Fill(nqb, fac, length, 1);

        for (i = 0; i < vCoordDim; ++i)
        {
            Vmath::Fill(nq_face, fac * normal[i][0], normal[i], 1);
        }
    }
    else
    {
        // Set up deformed normals
        int j, k;

        int nq0 = ptsKeys[0].GetNumPoints();
        int nq1 = ptsKeys[1].GetNumPoints();
        int nq2 = ptsKeys[2].GetNumPoints();
        int nqtot;
        int nq01 = nq0 * nq1;

        // number of elemental quad points
        if (face == 0)
        {
            nqtot = nq01;
        }
        else if (face == 1)
        {
            nqtot = nq0 * nq2;
        }
        else
        {
            nqtot = nq1 * nq2;
        }

        LibUtilities::PointsKey points0;
        LibUtilities::PointsKey points1;

        Array<OneD, NekDouble> faceJac(nqtot);
        Array<OneD, NekDouble> normals(vCoordDim * nqtot, 0.0);

        // Extract Jacobian along face and recover local derivates
        // (dx/dr) for polynomial interpolation by multiplying m_gmat by
        // jacobian
        switch (face)
        {
            case 0:
            {
                for (j = 0; j < nq01; ++j)
                {
                    normals[j]             = -df[2][j] * jac[j];
                    normals[nqtot + j]     = -df[5][j] * jac[j];
                    normals[2 * nqtot + j] = -df[8][j] * jac[j];
                    faceJac[j]             = jac[j];
                }

                points0 = ptsKeys[0];
                points1 = ptsKeys[1];
                break;
            }

            case 1:
            {
                for (j = 0; j < nq0; ++j)
                {
                    for (k = 0; k < nq2; ++k)
                    {
                        int tmp                      = j + nq01 * k;
                        normals[j + k * nq0]         = -df[1][tmp] * jac[tmp];
                        normals[nqtot + j + k * nq0] = -df[4][tmp] * jac[tmp];
                        normals[2 * nqtot + j + k * nq0] =
                            -df[7][tmp] * jac[tmp];
                        faceJac[j + k * nq0] = jac[tmp];
                    }
                }

                points0 = ptsKeys[0];
                points1 = ptsKeys[2];
                break;
            }

            case 2:
            {
                for (j = 0; j < nq1; ++j)
                {
                    for (k = 0; k < nq2; ++k)
                    {
                        int tmp = nq0 - 1 + nq0 * j + nq01 * k;
                        normals[j + k * nq1] =
                            (df[0][tmp] + df[1][tmp] + df[2][tmp]) * jac[tmp];
                        normals[nqtot + j + k * nq1] =
                            (df[3][tmp] + df[4][tmp] + df[5][tmp]) * jac[tmp];
                        normals[2 * nqtot + j + k * nq1] =
                            (df[6][tmp] + df[7][tmp] + df[8][tmp]) * jac[tmp];
                        faceJac[j + k * nq1] = jac[tmp];
                    }
                }

                points0 = ptsKeys[1];
                points1 = ptsKeys[2];
                break;
            }

            case 3:
            {
                for (j = 0; j < nq1; ++j)
                {
                    for (k = 0; k < nq2; ++k)
                    {
                        int tmp                      = j * nq0 + nq01 * k;
                        normals[j + k * nq1]         = -df[0][tmp] * jac[tmp];
                        normals[nqtot + j + k * nq1] = -df[3][tmp] * jac[tmp];
                        normals[2 * nqtot + j + k * nq1] =
                            -df[6][tmp] * jac[tmp];
                        faceJac[j + k * nq1] = jac[tmp];
                    }
                }

                points0 = ptsKeys[1];
                points1 = ptsKeys[2];
                break;
            }

            default:
                ASSERTL0(false, "face is out of range (face < 3)");
        }

        Array<OneD, NekDouble> work(nq_face, 0.0);
        // Interpolate Jacobian and invert
        LibUtilities::Interp2D(points0, points1, faceJac,
                               tobasis0.GetPointsKey(), tobasis1.GetPointsKey(),
                               work);
        Vmath::Sdiv(nq_face, 1.0, &work[0], 1, &work[0], 1);

        // Interpolate normal and multiply by inverse Jacobian.
        for (i = 0; i < vCoordDim; ++i)
        {
            LibUtilities::Interp2D(points0, points1, &normals[i * nqtot],
                                   tobasis0.GetPointsKey(),
                                   tobasis1.GetPointsKey(), &normal[i][0]);
            Vmath::Vmul(nq_face, work, 1, normal[i], 1, normal[i], 1);
        }

        // Normalise to obtain unit normals.
        Vmath::Zero(nq_face, work, 1);
        for (i = 0; i < GetCoordim(); ++i)
        {
            Vmath::Vvtvp(nq_face, normal[i], 1, normal[i], 1, work, 1, work, 1);
        }

        Vmath::Vsqrt(nq_face, work, 1, work, 1);
        Vmath::Sdiv(nq_face, 1.0, work, 1, work, 1);

        Vmath::Vcopy(nqb, work, 1, length, 1);

        for (i = 0; i < GetCoordim(); ++i)
        {
            Vmath::Vmul(nq_face, normal[i], 1, work, 1, normal[i], 1);
        }
    }
}

//-----------------------------
// Operator creation functions
//-----------------------------
void TetExp::v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    TetExp::v_LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
}

void TetExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void TetExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                  const StdRegions::StdMatrixKey &mkey)
{
    int nq = GetTotPoints();

    // Calculate sqrt of the Jacobian
    Array<OneD, const NekDouble> jac = m_geomFactors->GetJac();
    Array<OneD, NekDouble> sqrt_jac(nq);
    if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vsqrt(nq, jac, 1, sqrt_jac, 1);
    }
    else
    {
        Vmath::Fill(nq, sqrt(jac[0]), sqrt_jac, 1);
    }

    // Multiply array by sqrt(Jac)
    Vmath::Vmul(nq, sqrt_jac, 1, array, 1, array, 1);

    // Apply std region filter
    StdTetExp::v_SVVLaplacianFilter(array, mkey);

    // Divide by sqrt(Jac)
    Vmath::Vdiv(nq, array, 1, sqrt_jac, 1, array, 1);
}

//-----------------------------
// Matrix creation functions
//-----------------------------
DNekMatSharedPtr TetExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
{
    DNekMatSharedPtr returnval;

    switch (mkey.GetMatrixType())
    {
        case StdRegions::eHybridDGHelmholtz:
        case StdRegions::eHybridDGLamToU:
        case StdRegions::eHybridDGLamToQ0:
        case StdRegions::eHybridDGLamToQ1:
        case StdRegions::eHybridDGLamToQ2:
        case StdRegions::eHybridDGHelmBndLam:
        case StdRegions::eInvLaplacianWithUnityMean:
            returnval = Expansion3D::v_GenMatrix(mkey);
            break;
        default:
            returnval = StdTetExp::v_GenMatrix(mkey);
    }

    return returnval;
}

DNekMatSharedPtr TetExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
    LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
    LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
    StdRegions::StdTetExpSharedPtr tmp =
        MemoryManager<StdTetExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr TetExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}

void TetExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

DNekScalBlkMatSharedPtr TetExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}

void TetExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
{
    m_staticCondMatrixManager.DeleteObject(mkey);
}

void TetExp::GeneralMatrixOp_MatOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
{
    DNekScalMatSharedPtr mat = GetLocMatrix(mkey);

    if (inarray.data() == outarray.data())
    {
        Array<OneD, NekDouble> tmp(m_ncoeffs);
        Vmath::Vcopy(m_ncoeffs, inarray.data(), 1, tmp.data(), 1);

        Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, mat->Scale(),
                    (mat->GetOwnedMatrix())->GetPtr().data(), m_ncoeffs,
                    tmp.data(), 1, 0.0, outarray.data(), 1);
    }
    else
    {
        Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, mat->Scale(),
                    (mat->GetOwnedMatrix())->GetPtr().data(), m_ncoeffs,
                    inarray.data(), 1, 0.0, outarray.data(), 1);
    }
}

void TetExp::v_LaplacianMatrixOp_MatFree_Kernel(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, Array<OneD, NekDouble> &wsp)
{
    // This implementation is only valid when there are no
    // coefficients associated to the Laplacian operator
    if (m_metrics.count(eMetricLaplacian00) == 0)
    {
        ComputeLaplacianMetric();
    }

    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();
    int nqtot  = nquad0 * nquad1 * nquad2;

    ASSERTL1(wsp.size() >= 6 * nqtot, "Insufficient workspace size.");
    ASSERTL1(m_ncoeffs <= nqtot, "Workspace not set up for ncoeffs > nqtot");

    const Array<OneD, const NekDouble> &base0  = m_base[0]->GetBdata();
    const Array<OneD, const NekDouble> &base1  = m_base[1]->GetBdata();
    const Array<OneD, const NekDouble> &base2  = m_base[2]->GetBdata();
    const Array<OneD, const NekDouble> &dbase0 = m_base[0]->GetDbdata();
    const Array<OneD, const NekDouble> &dbase1 = m_base[1]->GetDbdata();
    const Array<OneD, const NekDouble> &dbase2 = m_base[2]->GetDbdata();
    const Array<OneD, const NekDouble> &metric00 =
        m_metrics[eMetricLaplacian00];
    const Array<OneD, const NekDouble> &metric01 =
        m_metrics[eMetricLaplacian01];
    const Array<OneD, const NekDouble> &metric02 =
        m_metrics[eMetricLaplacian02];
    const Array<OneD, const NekDouble> &metric11 =
        m_metrics[eMetricLaplacian11];
    const Array<OneD, const NekDouble> &metric12 =
        m_metrics[eMetricLaplacian12];
    const Array<OneD, const NekDouble> &metric22 =
        m_metrics[eMetricLaplacian22];

    // Allocate temporary storage
    Array<OneD, NekDouble> wsp0(2 * nqtot, wsp);
    Array<OneD, NekDouble> wsp1(nqtot, wsp + 1 * nqtot);
    Array<OneD, NekDouble> wsp2(nqtot, wsp + 2 * nqtot);
    Array<OneD, NekDouble> wsp3(nqtot, wsp + 3 * nqtot);
    Array<OneD, NekDouble> wsp4(nqtot, wsp + 4 * nqtot);
    Array<OneD, NekDouble> wsp5(nqtot, wsp + 5 * nqtot);

    // LAPLACIAN MATRIX OPERATION
    // wsp1 = du_dxi1 = D_xi1 * inarray = D_xi1 * u
    // wsp2 = du_dxi2 = D_xi2 * inarray = D_xi2 * u
    // wsp2 = du_dxi3 = D_xi3 * inarray = D_xi3 * u
    PhysTensorDeriv(inarray, wsp0, wsp1, wsp2);

    // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
    // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
    // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
    // especially for this purpose
    Vmath::Vvtvvtp(nqtot, &metric00[0], 1, &wsp0[0], 1, &metric01[0], 1,
                   &wsp1[0], 1, &wsp3[0], 1);
    Vmath::Vvtvp(nqtot, &metric02[0], 1, &wsp2[0], 1, &wsp3[0], 1, &wsp3[0], 1);
    Vmath::Vvtvvtp(nqtot, &metric01[0], 1, &wsp0[0], 1, &metric11[0], 1,
                   &wsp1[0], 1, &wsp4[0], 1);
    Vmath::Vvtvp(nqtot, &metric12[0], 1, &wsp2[0], 1, &wsp4[0], 1, &wsp4[0], 1);
    Vmath::Vvtvvtp(nqtot, &metric02[0], 1, &wsp0[0], 1, &metric12[0], 1,
                   &wsp1[0], 1, &wsp5[0], 1);
    Vmath::Vvtvp(nqtot, &metric22[0], 1, &wsp2[0], 1, &wsp5[0], 1, &wsp5[0], 1);

    // outarray = m = (D_xi1 * B)^T * k
    // wsp1     = n = (D_xi2 * B)^T * l
    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    bool Deformed = (m_geomFactors->GetGtype() == SpatialDomains::eDeformed);

    v_IProductWRTBaseKernel(dbase0, base1, base2, wsp3, outarray, jac,
                            Deformed);
    v_IProductWRTBaseKernel(base0, dbase1, base2, wsp4, wsp2, jac, Deformed);
    Vmath::Vadd(m_ncoeffs, wsp2.data(), 1, outarray.data(), 1, outarray.data(),
                1);
    v_IProductWRTBaseKernel(base0, base1, dbase2, wsp5, wsp2, jac, Deformed);
    Vmath::Vadd(m_ncoeffs, wsp2.data(), 1, outarray.data(), 1, outarray.data(),
                1);
}

void TetExp::v_ComputeLaplacianMetric()
{
    int i, j;
    const unsigned int nqtot = GetTotPoints();
    const unsigned int dim   = 3;
    const MetricType m[3][3] = {
        {eMetricLaplacian00, eMetricLaplacian01, eMetricLaplacian02},
        {eMetricLaplacian01, eMetricLaplacian11, eMetricLaplacian12},
        {eMetricLaplacian02, eMetricLaplacian12, eMetricLaplacian22}};

    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int j = i; j < dim; ++j)
        {
            m_metrics[m[i][j]] = Array<OneD, NekDouble>(nqtot);
        }
    }

    // Define shorthand synonyms for m_metrics storage
    Array<OneD, NekDouble> g0(m_metrics[m[0][0]]);
    Array<OneD, NekDouble> g1(m_metrics[m[1][1]]);
    Array<OneD, NekDouble> g2(m_metrics[m[2][2]]);
    Array<OneD, NekDouble> g3(m_metrics[m[0][1]]);
    Array<OneD, NekDouble> g4(m_metrics[m[0][2]]);
    Array<OneD, NekDouble> g5(m_metrics[m[1][2]]);

    // Allocate temporary storage
    Array<OneD, NekDouble> alloc(7 * nqtot, 0.0);
    Array<OneD, NekDouble> h0(alloc);               // h0
    Array<OneD, NekDouble> h1(alloc + 1 * nqtot);   // h1
    Array<OneD, NekDouble> h2(alloc + 2 * nqtot);   // h2
    Array<OneD, NekDouble> h3(alloc + 3 * nqtot);   // h3
    Array<OneD, NekDouble> wsp4(alloc + 4 * nqtot); // wsp4
    Array<OneD, NekDouble> wsp5(alloc + 5 * nqtot); // wsp5
    Array<OneD, NekDouble> wsp6(alloc + 6 * nqtot); // wsp6
    // Reuse some of the storage as workspace
    Array<OneD, NekDouble> wsp7(alloc);             // wsp7
    Array<OneD, NekDouble> wsp8(alloc + 1 * nqtot); // wsp8
    Array<OneD, NekDouble> wsp9(alloc + 2 * nqtot); // wsp9

    const Array<TwoD, const NekDouble> &df = m_geomFactors->GetDerivFactors();
    const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
    const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();
    const Array<OneD, const NekDouble> &z2 = m_base[2]->GetZ();
    const unsigned int nquad0              = m_base[0]->GetNumPoints();
    const unsigned int nquad1              = m_base[1]->GetNumPoints();
    const unsigned int nquad2              = m_base[2]->GetNumPoints();

    for (j = 0; j < nquad2; ++j)
    {
        for (i = 0; i < nquad1; ++i)
        {
            Vmath::Fill(nquad0, 4.0 / (1.0 - z1[i]) / (1.0 - z2[j]),
                        &h0[0] + i * nquad0 + j * nquad0 * nquad1, 1);
            Vmath::Fill(nquad0, 2.0 / (1.0 - z1[i]) / (1.0 - z2[j]),
                        &h1[0] + i * nquad0 + j * nquad0 * nquad1, 1);
            Vmath::Fill(nquad0, 2.0 / (1.0 - z2[j]),
                        &h2[0] + i * nquad0 + j * nquad0 * nquad1, 1);
            Vmath::Fill(nquad0, (1.0 + z1[i]) / (1.0 - z2[j]),
                        &h3[0] + i * nquad0 + j * nquad0 * nquad1, 1);
        }
    }
    for (i = 0; i < nquad0; i++)
    {
        Blas::Dscal(nquad1 * nquad2, 1 + z0[i], &h1[0] + i, nquad0);
    }

    // Step 3. Construct combined metric terms for physical space to
    // collapsed coordinate system.
    // Order of construction optimised to minimise temporary storage
    if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
    {
        // wsp4
        Vmath::Vadd(nqtot, &df[1][0], 1, &df[2][0], 1, &wsp4[0], 1);
        Vmath::Vvtvvtp(nqtot, &df[0][0], 1, &h0[0], 1, &wsp4[0], 1, &h1[0], 1,
                       &wsp4[0], 1);
        // wsp5
        Vmath::Vadd(nqtot, &df[4][0], 1, &df[5][0], 1, &wsp5[0], 1);
        Vmath::Vvtvvtp(nqtot, &df[3][0], 1, &h0[0], 1, &wsp5[0], 1, &h1[0], 1,
                       &wsp5[0], 1);
        // wsp6
        Vmath::Vadd(nqtot, &df[7][0], 1, &df[8][0], 1, &wsp6[0], 1);
        Vmath::Vvtvvtp(nqtot, &df[6][0], 1, &h0[0], 1, &wsp6[0], 1, &h1[0], 1,
                       &wsp6[0], 1);

        // g0
        Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0],
                       1, &g0[0], 1);
        Vmath::Vvtvp(nqtot, &wsp6[0], 1, &wsp6[0], 1, &g0[0], 1, &g0[0], 1);

        // g4
        Vmath::Vvtvvtp(nqtot, &df[2][0], 1, &wsp4[0], 1, &df[5][0], 1, &wsp5[0],
                       1, &g4[0], 1);
        Vmath::Vvtvp(nqtot, &df[8][0], 1, &wsp6[0], 1, &g4[0], 1, &g4[0], 1);

        // overwrite h0, h1, h2
        // wsp7 (h2f1 + h3f2)
        Vmath::Vvtvvtp(nqtot, &df[1][0], 1, &h2[0], 1, &df[2][0], 1, &h3[0], 1,
                       &wsp7[0], 1);
        // wsp8 (h2f4 + h3f5)
        Vmath::Vvtvvtp(nqtot, &df[4][0], 1, &h2[0], 1, &df[5][0], 1, &h3[0], 1,
                       &wsp8[0], 1);
        // wsp9 (h2f7 + h3f8)
        Vmath::Vvtvvtp(nqtot, &df[7][0], 1, &h2[0], 1, &df[8][0], 1, &h3[0], 1,
                       &wsp9[0], 1);

        // g3
        Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp7[0], 1, &wsp5[0], 1, &wsp8[0],
                       1, &g3[0], 1);
        Vmath::Vvtvp(nqtot, &wsp6[0], 1, &wsp9[0], 1, &g3[0], 1, &g3[0], 1);

        // overwrite wsp4, wsp5, wsp6
        // g1
        Vmath::Vvtvvtp(nqtot, &wsp7[0], 1, &wsp7[0], 1, &wsp8[0], 1, &wsp8[0],
                       1, &g1[0], 1);
        Vmath::Vvtvp(nqtot, &wsp9[0], 1, &wsp9[0], 1, &g1[0], 1, &g1[0], 1);

        // g5
        Vmath::Vvtvvtp(nqtot, &df[2][0], 1, &wsp7[0], 1, &df[5][0], 1, &wsp8[0],
                       1, &g5[0], 1);
        Vmath::Vvtvp(nqtot, &df[8][0], 1, &wsp9[0], 1, &g5[0], 1, &g5[0], 1);

        // g2
        Vmath::Vvtvvtp(nqtot, &df[2][0], 1, &df[2][0], 1, &df[5][0], 1,
                       &df[5][0], 1, &g2[0], 1);
        Vmath::Vvtvp(nqtot, &df[8][0], 1, &df[8][0], 1, &g2[0], 1, &g2[0], 1);
    }
    else
    {
        // wsp4
        Vmath::Svtsvtp(nqtot, df[0][0], &h0[0], 1, df[1][0] + df[2][0], &h1[0],
                       1, &wsp4[0], 1);
        // wsp5
        Vmath::Svtsvtp(nqtot, df[3][0], &h0[0], 1, df[4][0] + df[5][0], &h1[0],
                       1, &wsp5[0], 1);
        // wsp6
        Vmath::Svtsvtp(nqtot, df[6][0], &h0[0], 1, df[7][0] + df[8][0], &h1[0],
                       1, &wsp6[0], 1);

        // g0
        Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0],
                       1, &g0[0], 1);
        Vmath::Vvtvp(nqtot, &wsp6[0], 1, &wsp6[0], 1, &g0[0], 1, &g0[0], 1);

        // g4
        Vmath::Svtsvtp(nqtot, df[2][0], &wsp4[0], 1, df[5][0], &wsp5[0], 1,
                       &g4[0], 1);
        Vmath::Svtvp(nqtot, df[8][0], &wsp6[0], 1, &g4[0], 1, &g4[0], 1);

        // overwrite h0, h1, h2
        // wsp7 (h2f1 + h3f2)
        Vmath::Svtsvtp(nqtot, df[1][0], &h2[0], 1, df[2][0], &h3[0], 1,
                       &wsp7[0], 1);
        // wsp8 (h2f4 + h3f5)
        Vmath::Svtsvtp(nqtot, df[4][0], &h2[0], 1, df[5][0], &h3[0], 1,
                       &wsp8[0], 1);
        // wsp9 (h2f7 + h3f8)
        Vmath::Svtsvtp(nqtot, df[7][0], &h2[0], 1, df[8][0], &h3[0], 1,
                       &wsp9[0], 1);

        // g3
        Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp7[0], 1, &wsp5[0], 1, &wsp8[0],
                       1, &g3[0], 1);
        Vmath::Vvtvp(nqtot, &wsp6[0], 1, &wsp9[0], 1, &g3[0], 1, &g3[0], 1);

        // overwrite wsp4, wsp5, wsp6
        // g1
        Vmath::Vvtvvtp(nqtot, &wsp7[0], 1, &wsp7[0], 1, &wsp8[0], 1, &wsp8[0],
                       1, &g1[0], 1);
        Vmath::Vvtvp(nqtot, &wsp9[0], 1, &wsp9[0], 1, &g1[0], 1, &g1[0], 1);

        // g5
        Vmath::Svtsvtp(nqtot, df[2][0], &wsp7[0], 1, df[5][0], &wsp8[0], 1,
                       &g5[0], 1);
        Vmath::Svtvp(nqtot, df[8][0], &wsp9[0], 1, &g5[0], 1, &g5[0], 1);

        // g2
        Vmath::Fill(nqtot,
                    df[2][0] * df[2][0] + df[5][0] * df[5][0] +
                        df[8][0] * df[8][0],
                    &g2[0], 1);
    }
}
} // namespace Nektar::LocalRegions
