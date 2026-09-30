///////////////////////////////////////////////////////////////////////////////
//
// File: HexExp.cpp
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
// Description: Methods for Hex expansion in local regoins
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LocalRegions/HexExp.h>
#include <SpatialDomains/HexGeom.h>

using namespace std;

namespace Nektar::LocalRegions
{
/**
 * @class HexExp
 * Defines a hexahedral local expansion.
 */

/**
 * \brief Constructor using BasisKey class for quadrature points and
 * order definition
 *
 * @param   Ba          Basis key for first coordinate.
 * @param   Bb          Basis key for second coordinate.
 * @param   Bc          Basis key for third coordinate.
 */
HexExp::HexExp(const LibUtilities::BasisKey &Ba,
               const LibUtilities::BasisKey &Bb,
               const LibUtilities::BasisKey &Bc,
               SpatialDomains::Geometry3D *geom)
    : StdExpansion(Ba.GetNumModes() * Bb.GetNumModes() * Bc.GetNumModes(), 3,
                   Ba, Bb, Bc),
      StdExpansion3D(Ba.GetNumModes() * Bb.GetNumModes() * Bc.GetNumModes(), Ba,
                     Bb, Bc),
      StdHexExp(Ba, Bb, Bc), Expansion(geom), Expansion3D(geom),
      m_matrixManager(
          std::bind(&Expansion3D::CreateMatrix, this, std::placeholders::_1)),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1))
{
}

/**
 * \brief Copy Constructor
 *
 * @param   T           HexExp to copy.
 */
HexExp::HexExp(const HexExp &T)
    : StdExpansion(T), StdExpansion3D(T), StdHexExp(T), Expansion(T),
      Expansion3D(T), m_matrixManager(T.m_matrixManager),
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
 * \f$\phi_{pqr}(\xi_1,\xi_2,\xi_3) =
 * \phi_1(\xi_1)\phi_2(\xi_2)\phi_3(\xi_3)\f$, in the hexahedral
 * element, this is straightforward and yields the result
 *
 * \f[
 * I_{pqr} = \sum_{k=1}^3 \left(u, \frac{\partial u}{\partial \xi_k}
 * \frac{\partial \xi_k}{\partial x_i}\right)
 * \f]
 *
 * @param dir       Direction in which to take the derivative.
 * @param inarray   The function \f$ u \f$.
 * @param outarray  Value of the inner product.
 */
void HexExp::v_IProductWRTDerivBase(const int dir,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    ASSERTL1((dir == 0) || (dir == 1) || (dir == 2), "Invalid direction.");

    const int nq0 = m_base[0]->GetNumPoints();
    const int nq1 = m_base[1]->GetNumPoints();
    const int nq2 = m_base[2]->GetNumPoints();
    const int nq  = nq0 * nq1 * nq2;

    Array<OneD, NekDouble> tmp2(nq);        // Dir1 metric
    Array<OneD, NekDouble> tmp3(nq);        // Dir2 metric
    Array<OneD, NekDouble> tmp4(nq);        // Dir3 metric
    Array<OneD, NekDouble> tmp5(m_ncoeffs); // iprod tmp

    Array<OneD, Array<OneD, NekDouble>> tmp2D{3};
    tmp2D[0] = tmp2;
    tmp2D[1] = tmp3;
    tmp2D[2] = tmp4;

    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    bool Deformed = (m_geomFactors->GetGtype() == SpatialDomains::eDeformed);

    const bool CollDir0 = m_base[0]->Collocation();
    const bool CollDir1 = m_base[1]->Collocation();
    const bool CollDir2 = m_base[2]->Collocation();

    HexExp::v_AlignVectorToCollapsedDir(dir, inarray, tmp2D);

    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                            m_base[2]->GetBdata(), tmp2, outarray, jac,
                            Deformed, false, CollDir1, CollDir2);

    v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                            m_base[2]->GetBdata(), tmp3, tmp5, jac, Deformed,
                            CollDir0, false, CollDir2);
    Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);

    v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                            m_base[2]->GetDbdata(), tmp4, tmp5, jac, Deformed,
                            CollDir0, CollDir1, false);
    Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);
}

void HexExp::v_AlignVectorToCollapsedDir(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    ASSERTL1((dir == 0) || (dir == 1) || (dir == 2), "Invalid direction.");

    const int nq0 = m_base[0]->GetNumPoints();
    const int nq1 = m_base[1]->GetNumPoints();
    const int nq2 = m_base[2]->GetNumPoints();
    const int nq  = nq0 * nq1 * nq2;

    const Array<TwoD, const NekDouble> &df = m_geomFactors->GetDerivFactors();

    Array<OneD, NekDouble> tmp1(nq); // Quad metric

    Array<OneD, NekDouble> tmp2 = outarray[0]; // Dir1 metric
    Array<OneD, NekDouble> tmp3 = outarray[1]; // Dir2 metric
    Array<OneD, NekDouble> tmp4 = outarray[2];

    Vmath::Vcopy(nq, inarray, 1, tmp1, 1); // Dir3 metric

    if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vmul(nq, &df[3 * dir][0], 1, tmp1.data(), 1, tmp2.data(), 1);
        Vmath::Vmul(nq, &df[3 * dir + 1][0], 1, tmp1.data(), 1, tmp3.data(), 1);
        Vmath::Vmul(nq, &df[3 * dir + 2][0], 1, tmp1.data(), 1, tmp4.data(), 1);
    }
    else
    {
        Vmath::Smul(nq, df[3 * dir][0], tmp1.data(), 1, tmp2.data(), 1);
        Vmath::Smul(nq, df[3 * dir + 1][0], tmp1.data(), 1, tmp3.data(), 1);
        Vmath::Smul(nq, df[3 * dir + 2][0], tmp1.data(), 1, tmp4.data(), 1);
    }
}

/**
 *
 * @param dir       Vector direction in which to take the derivative.
 * @param inarray   The function \f$ u \f$.
 * @param outarray  Value of the inner product.
 */
void HexExp::v_IProductWRTDirectionalDerivBase(
    const Array<OneD, const NekDouble> &direction,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int shapedim  = 3;
    const int nq0 = m_base[0]->GetNumPoints();
    const int nq1 = m_base[1]->GetNumPoints();
    const int nq2 = m_base[2]->GetNumPoints();
    const int nq  = nq0 * nq1 * nq2;

    Array<OneD, NekDouble> tmp2(nq);        // Dir1 metric
    Array<OneD, NekDouble> tmp3(nq);        // Dir2 metric
    Array<OneD, NekDouble> tmp4(nq);        // Dir3 metric
    Array<OneD, NekDouble> tmp5(m_ncoeffs); // iprod tmp

    Array<OneD, Array<OneD, NekDouble>> tmp2D{3};
    tmp2D[0] = tmp2;
    tmp2D[1] = tmp3;
    tmp2D[2] = tmp4;

    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    bool Deformed = (m_geomFactors->GetGtype() == SpatialDomains::eDeformed);

    const bool CollDir0 = m_base[0]->Collocation();
    const bool CollDir1 = m_base[1]->Collocation();
    const bool CollDir2 = m_base[2]->Collocation();

    const Array<TwoD, const NekDouble> &df = m_geomFactors->GetDerivFactors();

    Array<OneD, Array<OneD, NekDouble>> dfdir(shapedim);
    Expansion::ComputeGmatcdotMF(df, direction, dfdir);

    Vmath::Vmul(nq, &dfdir[0][0], 1, inarray.data(), 1, tmp2.data(), 1);
    Vmath::Vmul(nq, &dfdir[1][0], 1, inarray.data(), 1, tmp3.data(), 1);
    Vmath::Vmul(nq, &dfdir[2][0], 1, inarray.data(), 1, tmp4.data(), 1);

    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                            m_base[2]->GetBdata(), tmp2, outarray, jac,
                            Deformed, false, CollDir1, CollDir2);

    v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                            m_base[2]->GetBdata(), tmp3, tmp5, jac, Deformed,
                            CollDir0, false, CollDir2);
    Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);

    v_IProductWRTBaseKernel(m_base[0]->GetBdata(), m_base[1]->GetBdata(),
                            m_base[2]->GetDbdata(), tmp4, tmp5, jac, Deformed,
                            CollDir0, CollDir1, false);
    Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);
}

//-----------------------------
// Evaluation functions
//-----------------------------
NekDouble HexExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    Array<OneD, NekDouble> Lcoord(3);
    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);
    return StdHexExp::v_PhysEvalFirstDeriv(Lcoord, inarray, firstOrderDerivs);
}

StdRegions::StdExpansionSharedPtr HexExp::v_GetStdExp(void) const
{
    return MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey(), m_base[1]->GetBasisKey(),
        m_base[2]->GetBasisKey());
}

StdRegions::StdExpansionSharedPtr HexExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());
    LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(), 2,
                                 m_base[1]->GetPointsKey());
    LibUtilities::BasisKey bkey2(m_base[2]->GetBasisType(), 2,
                                 m_base[2]->GetPointsKey());

    return MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(bkey0, bkey1,
                                                                   bkey2);
}

/**
 * \brief Retrieves the physical coordinates of a given set of
 * reference coordinates.
 *
 * @param   Lcoords     Local coordinates in reference space.
 * @param   coords      Corresponding coordinates in physical space.
 */
void HexExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                        Array<OneD, NekDouble> &coords)
{
    int i;

    ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[0] <= 1.0 && Lcoords[1] >= -1.0 &&
                 Lcoords[1] <= 1.0 && Lcoords[2] >= -1.0 && Lcoords[2] <= 1.0,
             "Local coordinates are not in region [-1,1]");

    m_geom->FillGeom();

    for (i = 0; i < m_geom->GetCoordim(); ++i)
    {
        coords[i] = m_geom->GetCoord(i, Lcoords);
    }
}

void HexExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                         Array<OneD, NekDouble> &coords_1,
                         Array<OneD, NekDouble> &coords_2)
{
    Expansion::v_GetCoords(coords_0, coords_1, coords_2);
}

//-----------------------------
// Helper functions
//-----------------------------
void HexExp::v_ExtractDataToCoeffs(
    const NekDouble *data, const std::vector<unsigned int> &nummodes,
    const int mode_offset, NekDouble *coeffs,
    std::vector<LibUtilities::BasisType> &fromType)
{
    int data_order0 = nummodes[mode_offset];
    int fillorder0  = min(m_base[0]->GetNumModes(), data_order0);
    int data_order1 = nummodes[mode_offset + 1];
    int order1      = m_base[1]->GetNumModes();
    int fillorder1  = min(order1, data_order1);
    int data_order2 = nummodes[mode_offset + 2];
    int order2      = m_base[2]->GetNumModes();
    int fillorder2  = min(order2, data_order2);

    // Check if same basis
    if (fromType[0] != m_base[0]->GetBasisType() ||
        fromType[1] != m_base[1]->GetBasisType() ||
        fromType[2] != m_base[2]->GetBasisType())
    {
        // Construct a hex with the appropriate basis type at our
        // quadrature points, and one more to do a forwards
        // transform. We can then copy the output to coeffs.
        StdRegions::StdHexExp tmpHex(
            LibUtilities::BasisKey(fromType[0], data_order0,
                                   m_base[0]->GetPointsKey()),
            LibUtilities::BasisKey(fromType[1], data_order1,
                                   m_base[1]->GetPointsKey()),
            LibUtilities::BasisKey(fromType[2], data_order2,
                                   m_base[2]->GetPointsKey()));
        StdRegions::StdHexExp tmpHex2(m_base[0]->GetBasisKey(),
                                      m_base[1]->GetBasisKey(),
                                      m_base[2]->GetBasisKey());

        Array<OneD, const NekDouble> tmpData(tmpHex.GetNcoeffs(), data);
        Array<OneD, NekDouble> tmpBwd(tmpHex2.GetTotPoints());
        Array<OneD, NekDouble> tmpOut(tmpHex2.GetNcoeffs());

        tmpHex.BwdTrans(tmpData, tmpBwd);
        tmpHex2.FwdTrans(tmpBwd, tmpOut);
        Vmath::Vcopy(tmpOut.size(), &tmpOut[0], 1, coeffs, 1);

        return;
    }

    switch (m_base[0]->GetBasisType())
    {
        case LibUtilities::eModified_A:
        {
            int i, j;
            int cnt  = 0;
            int cnt1 = 0;

            ASSERTL1(m_base[1]->GetBasisType() == LibUtilities::eModified_A,
                     "Extraction routine not set up for this basis");
            ASSERTL1(m_base[2]->GetBasisType() == LibUtilities::eModified_A,
                     "Extraction routine not set up for this basis");

            Vmath::Zero(m_ncoeffs, coeffs, 1);
            for (j = 0; j < fillorder0; ++j)
            {
                for (i = 0; i < fillorder1; ++i)
                {
                    Vmath::Vcopy(fillorder2, &data[cnt], 1, &coeffs[cnt1], 1);
                    cnt += data_order2;
                    cnt1 += order2;
                }

                // count out data for j iteration
                for (i = fillorder1; i < data_order1; ++i)
                {
                    cnt += data_order2;
                }

                for (i = fillorder1; i < order1; ++i)
                {
                    cnt1 += order2;
                }
            }
            break;
        }
        case LibUtilities::eGLL_Lagrange:
        {
            LibUtilities::PointsKey p0(nummodes[0],
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::PointsKey p1(nummodes[1],
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::PointsKey p2(nummodes[2],
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::PointsKey t0(m_base[0]->GetNumModes(),
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::PointsKey t1(m_base[1]->GetNumModes(),
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::PointsKey t2(m_base[2]->GetNumModes(),
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::Interp3D(p0, p1, p2, data, t0, t1, t2, coeffs);
        }
        break;
        default:
            ASSERTL0(false, "basis is either not set up or not "
                            "hierarchicial");
    }
}

void HexExp::v_GetTracePhysMap(const int face, Array<OneD, int> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nquad2 = m_base[2]->GetNumPoints();

    int nq0 = 0;
    int nq1 = 0;

    switch (face)
    {
        case 0:
            nq0 = nquad0;
            nq1 = nquad1;

            // Directions A and B positive
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            for (int i = 0; i < nquad0 * nquad1; ++i)
            {
                outarray[i] = i;
            }

            break;
        case 1:
            nq0 = nquad0;
            nq1 = nquad2;
            // Direction A and B positive
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            // Direction A and B positive
            for (int k = 0; k < nquad2; k++)
            {
                for (int i = 0; i < nquad0; ++i)
                {
                    outarray[k * nquad0 + i] = nquad0 * nquad1 * k + i;
                }
            }
            break;
        case 2:
            nq0 = nquad1;
            nq1 = nquad2;

            // Direction A and B positive
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            for (int i = 0; i < nquad1 * nquad2; i++)
            {
                outarray[i] = nquad0 - 1 + i * nquad0;
            }
            break;
        case 3:
            nq0 = nquad0;
            nq1 = nquad2;

            // Direction A and B positive
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            for (int k = 0; k < nquad2; k++)
            {
                for (int i = 0; i < nquad0; i++)
                {
                    outarray[k * nquad0 + i] =
                        (nquad0 * (nquad1 - 1)) + (k * nquad0 * nquad1) + i;
                }
            }
            break;
        case 4:
            nq0 = nquad1;
            nq1 = nquad2;

            // Direction A and B positive
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            for (int i = 0; i < nquad1 * nquad2; i++)
            {
                outarray[i] = i * nquad0;
            }
            break;
        case 5:
            nq0 = nquad0;
            nq1 = nquad1;
            // Directions A and B positive
            if (outarray.size() != nq0 * nq1)
            {
                outarray = Array<OneD, int>(nq0 * nq1);
            }

            for (int i = 0; i < nquad0 * nquad1; i++)
            {
                outarray[i] = nquad0 * nquad1 * (nquad2 - 1) + i;
            }

            break;
        default:
            ASSERTL0(false, "face value (> 5) is out of range");
            break;
    }
}

void HexExp::v_ComputeTraceNormal(const int face)
{
    int i;
    SpatialDomains::GeomType type = m_geomFactors->GetGtype();

    LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
    for (i = 0; i < ptsKeys.size(); ++i)
    {
        // Need at least 2 points for computing normals
        if (ptsKeys[i].GetNumPoints() == 1)
        {
            LibUtilities::PointsKey pKey(2, ptsKeys[i].GetPointsType());
            ptsKeys[i] = pKey;
        }
    }

    const Array<TwoD, const NekDouble> &df =
        m_geomFactors->ComputeDerivFactors(ptsKeys);
    const Array<OneD, const NekDouble> &jac =
        m_geomFactors->ComputeJac(ptsKeys);

    LibUtilities::BasisKey tobasis0 = GetTraceBasisKey(face, 0);
    LibUtilities::BasisKey tobasis1 = GetTraceBasisKey(face, 1);

    // Number of quadrature points in face expansion.
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
    if ((type == SpatialDomains::eRegular) ||
        (type == SpatialDomains::eMovingRegular))
    {
        NekDouble fac;
        // Set up normals
        switch (face)
        {
            case 0:
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = -df[3 * i + 2][0];
                }
                break;
            case 1:
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = -df[3 * i + 1][0];
                }
                break;
            case 2:
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = df[3 * i][0];
                }
                break;
            case 3:
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = df[3 * i + 1][0];
                }
                break;
            case 4:
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = -df[3 * i][0];
                }
                break;
            case 5:
                for (i = 0; i < vCoordDim; ++i)
                {
                    normal[i][0] = df[3 * i + 2][0];
                }
                break;
            default:
                ASSERTL0(false, "face is out of range (edge < 5)");
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
    else // Set up deformed normals
    {
        int j, k;

        int nqe0  = ptsKeys[0].GetNumPoints();
        int nqe1  = ptsKeys[1].GetNumPoints();
        int nqe2  = ptsKeys[2].GetNumPoints();
        int nqe01 = nqe0 * nqe1;
        int nqe02 = nqe0 * nqe2;
        int nqe12 = nqe1 * nqe2;

        int nqe;
        if (face == 0 || face == 5)
        {
            nqe = nqe01;
        }
        else if (face == 1 || face == 3)
        {
            nqe = nqe02;
        }
        else
        {
            nqe = nqe12;
        }

        LibUtilities::PointsKey points0;
        LibUtilities::PointsKey points1;

        Array<OneD, NekDouble> faceJac(nqe);
        Array<OneD, NekDouble> normals(vCoordDim * nqe, 0.0);

        // Extract Jacobian along face and recover local
        // derivates (dx/dr) for polynomial interpolation by
        // multiplying m_gmat by jacobian
        switch (face)
        {
            case 0:
                for (j = 0; j < nqe; ++j)
                {
                    normals[j]           = -df[2][j] * jac[j];
                    normals[nqe + j]     = -df[5][j] * jac[j];
                    normals[2 * nqe + j] = -df[8][j] * jac[j];
                    faceJac[j]           = jac[j];
                }

                points0 = ptsKeys[0];
                points1 = ptsKeys[1];
                break;
            case 1:
                for (j = 0; j < nqe0; ++j)
                {
                    for (k = 0; k < nqe2; ++k)
                    {
                        int idx                     = j + nqe01 * k;
                        normals[j + k * nqe0]       = -df[1][idx] * jac[idx];
                        normals[nqe + j + k * nqe0] = -df[4][idx] * jac[idx];
                        normals[2 * nqe + j + k * nqe0] =
                            -df[7][idx] * jac[idx];
                        faceJac[j + k * nqe0] = jac[idx];
                    }
                }
                points0 = ptsKeys[0];
                points1 = ptsKeys[2];
                break;
            case 2:
                for (j = 0; j < nqe1; ++j)
                {
                    for (k = 0; k < nqe2; ++k)
                    {
                        int idx               = nqe0 - 1 + nqe0 * j + nqe01 * k;
                        normals[j + k * nqe1] = df[0][idx] * jac[idx];
                        normals[nqe + j + k * nqe1]     = df[3][idx] * jac[idx];
                        normals[2 * nqe + j + k * nqe1] = df[6][idx] * jac[idx];
                        faceJac[j + k * nqe1]           = jac[idx];
                    }
                }
                points0 = ptsKeys[1];
                points1 = ptsKeys[2];
                break;
            case 3:
                for (j = 0; j < nqe0; ++j)
                {
                    for (k = 0; k < nqe2; ++k)
                    {
                        int idx = nqe0 * (nqe1 - 1) + j + nqe01 * k;
                        normals[j + k * nqe0]           = df[1][idx] * jac[idx];
                        normals[nqe + j + k * nqe0]     = df[4][idx] * jac[idx];
                        normals[2 * nqe + j + k * nqe0] = df[7][idx] * jac[idx];
                        faceJac[j + k * nqe0]           = jac[idx];
                    }
                }
                points0 = ptsKeys[0];
                points1 = ptsKeys[2];
                break;
            case 4:
                for (j = 0; j < nqe1; ++j)
                {
                    for (k = 0; k < nqe2; ++k)
                    {
                        int idx                     = j * nqe0 + nqe01 * k;
                        normals[j + k * nqe1]       = -df[0][idx] * jac[idx];
                        normals[nqe + j + k * nqe1] = -df[3][idx] * jac[idx];
                        normals[2 * nqe + j + k * nqe1] =
                            -df[6][idx] * jac[idx];
                        faceJac[j + k * nqe1] = jac[idx];
                    }
                }
                points0 = ptsKeys[1];
                points1 = ptsKeys[2];
                break;
            case 5:
                for (j = 0; j < nqe01; ++j)
                {
                    int idx              = j + nqe01 * (nqe2 - 1);
                    normals[j]           = df[2][idx] * jac[idx];
                    normals[nqe + j]     = df[5][idx] * jac[idx];
                    normals[2 * nqe + j] = df[8][idx] * jac[idx];
                    faceJac[j]           = jac[idx];
                }
                points0 = ptsKeys[0];
                points1 = ptsKeys[1];
                break;
            default:
                ASSERTL0(false, "face is out of range (face < 5)");
        }

        Array<OneD, NekDouble> work(nq_face, 0.0);
        // Interpolate Jacobian and invert
        LibUtilities::Interp2D(points0, points1, faceJac,
                               tobasis0.GetPointsKey(), tobasis1.GetPointsKey(),
                               work);

        Vmath::Sdiv(nq_face, 1.0, &work[0], 1, &work[0], 1);

        // interpolate
        for (i = 0; i < GetCoordim(); ++i)
        {
            LibUtilities::Interp2D(points0, points1, &normals[i * nqe],
                                   tobasis0.GetPointsKey(),
                                   tobasis1.GetPointsKey(), &normal[i][0]);
            Vmath::Vmul(nq_face, work, 1, normal[i], 1, normal[i], 1);
        }

        // normalise normal vectors
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
void HexExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void HexExp::v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    HexExp::v_LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
}

void HexExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void HexExp::v_WeakDerivMatrixOp(const int i,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void HexExp::v_WeakDirectionalDerivMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(inarray, outarray, mkey);
}

void HexExp::v_MassLevelCurvatureMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::MassLevelCurvatureMatrixOp_MatFree(inarray, outarray, mkey);
}

void HexExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    HexExp::v_HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
}

/**
 * This function is used to compute exactly the advective numerical flux
 * on the interface of two elements with different expansions, hence an
 * appropriate number of Gauss points has to be used. The number of
 * Gauss points has to be equal to the number used by the highest
 * polynomial degree of the two adjacent elements
 *
 * @param   numMin     Is the reduced polynomial order
 * @param   inarray    Input array of coefficients
 * @param   dumpVar    Output array of reduced coefficients.
 */
void HexExp::v_ReduceOrderCoeffs(int numMin,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray)
{
    int n_coeffs = inarray.size();
    int nmodes0  = m_base[0]->GetNumModes();
    int nmodes1  = m_base[1]->GetNumModes();
    int nmodes2  = m_base[2]->GetNumModes();
    int numMax   = nmodes0;

    Array<OneD, NekDouble> coeff(n_coeffs);
    Array<OneD, NekDouble> coeff_tmp1(nmodes0 * nmodes1, 0.0);
    Array<OneD, NekDouble> coeff_tmp2(n_coeffs, 0.0);
    Array<OneD, NekDouble> tmp, tmp2, tmp3, tmp4;

    Vmath::Vcopy(n_coeffs, inarray, 1, coeff_tmp2, 1);

    const LibUtilities::PointsKey Pkey0(nmodes0,
                                        LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey Pkey1(nmodes1,
                                        LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey Pkey2(nmodes2,
                                        LibUtilities::eGaussLobattoLegendre);

    LibUtilities::BasisKey b0(m_base[0]->GetBasisType(), nmodes0, Pkey0);
    LibUtilities::BasisKey b1(m_base[1]->GetBasisType(), nmodes1, Pkey1);
    LibUtilities::BasisKey b2(m_base[2]->GetBasisType(), nmodes2, Pkey2);
    LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A, nmodes0, Pkey0);
    LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_A, nmodes1, Pkey1);
    LibUtilities::BasisKey bortho2(LibUtilities::eOrtho_A, nmodes2, Pkey2);

    LibUtilities::InterpCoeff3D(b0, b1, b2, coeff_tmp2, bortho0, bortho1,
                                bortho2, coeff);

    Vmath::Zero(n_coeffs, coeff_tmp2, 1);

    int cnt = 0, cnt2 = 0;

    for (int u = 0; u < numMin + 1; ++u)
    {
        for (int i = 0; i < numMin; ++i)
        {
            Vmath::Vcopy(numMin, tmp = coeff + cnt + cnt2, 1,
                         tmp2 = coeff_tmp1 + cnt, 1);

            cnt = i * numMax;
        }

        Vmath::Vcopy(nmodes0 * nmodes1, tmp3 = coeff_tmp1, 1,
                     tmp4 = coeff_tmp2 + cnt2, 1);

        cnt2 = u * nmodes0 * nmodes1;
    }

    LibUtilities::InterpCoeff3D(bortho0, bortho1, bortho2, coeff_tmp2, b0, b1,
                                b2, outarray);
}

void HexExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
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
    StdHexExp::v_SVVLaplacianFilter(array, mkey);

    // Divide by sqrt(Jac)
    Vmath::Vdiv(nq, array, 1, sqrt_jac, 1, array, 1);
}

//-----------------------------
// Matrix creation functions
//-----------------------------
DNekMatSharedPtr HexExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
            returnval = StdHexExp::v_GenMatrix(mkey);
    }

    return returnval;
}

DNekMatSharedPtr HexExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
    LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
    LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();

    StdRegions::StdHexExpSharedPtr tmp =
        MemoryManager<StdHexExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr HexExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}

void HexExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

DNekScalBlkMatSharedPtr HexExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}

void HexExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
{
    m_staticCondMatrixManager.DeleteObject(mkey);
}

void HexExp::v_LaplacianMatrixOp_MatFree_Kernel(
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
    Array<OneD, NekDouble> wsp0(wsp);
    Array<OneD, NekDouble> wsp1(wsp + 1 * nqtot);
    Array<OneD, NekDouble> wsp2(wsp + 2 * nqtot);
    Array<OneD, NekDouble> wsp3(wsp + 3 * nqtot);
    Array<OneD, NekDouble> wsp4(wsp + 4 * nqtot);
    Array<OneD, NekDouble> wsp5(wsp + 5 * nqtot);

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

    const bool CollDir0 = m_base[0]->Collocation();
    const bool CollDir1 = m_base[1]->Collocation();
    const bool CollDir2 = m_base[2]->Collocation();

    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    bool Deformed = (m_geomFactors->GetGtype() == SpatialDomains::eDeformed);

    // outarray = m = (D_xi1 * B)^T * k
    // wsp1     = n = (D_xi2 * B)^T * l
    v_IProductWRTBaseKernel(dbase0, base1, base2, wsp3, outarray, jac, Deformed,
                            false, CollDir1, CollDir2);
    v_IProductWRTBaseKernel(base0, dbase1, base2, wsp4, wsp2, jac, Deformed,
                            CollDir0, false, CollDir2);
    Vmath::Vadd(m_ncoeffs, wsp2.data(), 1, outarray.data(), 1, outarray.data(),
                1);
    v_IProductWRTBaseKernel(base0, base1, dbase2, wsp5, wsp2, jac, Deformed,
                            CollDir0, CollDir1, false);
    Vmath::Vadd(m_ncoeffs, wsp2.data(), 1, outarray.data(), 1, outarray.data(),
                1);
}

void HexExp::v_ComputeLaplacianMetric()
{
    const SpatialDomains::GeomType type = m_geomFactors->GetGtype();
    const unsigned int nqtot            = GetTotPoints();
    const unsigned int dim              = 3;
    const MetricType m[3][3]            = {
        {eMetricLaplacian00, eMetricLaplacian01, eMetricLaplacian02},
        {eMetricLaplacian01, eMetricLaplacian11, eMetricLaplacian12},
        {eMetricLaplacian02, eMetricLaplacian12, eMetricLaplacian22}};

    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int j = i; j < dim; ++j)
        {
            m_metrics[m[i][j]] = Array<OneD, NekDouble>(nqtot);
            const Array<TwoD, const NekDouble> &gmat =
                m_geomFactors->GetGmat(GetPointsKeys());
            if (type == SpatialDomains::eDeformed)
            {
                Vmath::Vcopy(nqtot, &gmat[i * dim + j][0], 1,
                             &m_metrics[m[i][j]][0], 1);
            }
            else
            {
                Vmath::Fill(nqtot, gmat[i * dim + j][0], &m_metrics[m[i][j]][0],
                            1);
            }
        }
    }
}

} // namespace Nektar::LocalRegions
