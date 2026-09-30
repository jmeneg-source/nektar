///////////////////////////////////////////////////////////////////////////////
//
// File: SegExp.cpp
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
// Description: SegExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Interp.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/SegExp.h>

using namespace std;

namespace Nektar::LocalRegions
{

/**
 * @class SegExp
 * Defines a Segment local expansion.
 */

/// Constructor using BasisKey class for quadrature points and
/// order definition.
/**
 * @param   Ba          Basis key of segment expansion.
 * @param   geom        Description of geometry.
 */
SegExp::SegExp(const LibUtilities::BasisKey &Ba,
               SpatialDomains::Geometry1D *geom)
    : StdExpansion(Ba.GetNumModes(), 1, Ba),
      StdExpansion1D(Ba.GetNumModes(), Ba), StdRegions::StdSegExp(Ba),
      Expansion(geom), Expansion1D(geom),
      m_matrixManager(
          std::bind(&SegExp::CreateMatrix, this, std::placeholders::_1)),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1))
{
}

/// Copy Constructor
/**
 * @param   S           Existing segment to duplicate.
 */
SegExp::SegExp(const SegExp &S)
    : StdExpansion(S), StdExpansion1D(S), StdRegions::StdSegExp(S),
      Expansion(S), Expansion1D(S), m_matrixManager(S.m_matrixManager),
      m_staticCondMatrixManager(S.m_staticCondMatrixManager)
{
}

//-----------------------------
// Transforms
//-----------------------------
void SegExp::v_FwdTransBndConstrained(
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
        bool hasEndPoints = true;
        bool hasEndModes  = true;

        switch (m_base[0]->GetBasisType())
        {
            case LibUtilities::eGLL_Lagrange:
            {
                nInteriorDofs = m_ncoeffs - 2;
                offset        = 1;
                hasEndModes   = true;
            }
            break;
            case LibUtilities::eOrtho_A:
            case LibUtilities::eGauss_Lagrange:
            {
                nInteriorDofs = m_ncoeffs;
                offset        = 0;
                hasEndModes   = false;
            }
            break;
            case LibUtilities::eModified_A:
            case LibUtilities::eModified_B:
            {
                nInteriorDofs = m_ncoeffs - 2;
                offset        = 2;
                hasEndModes   = true;
            }
            break;
            default:
                ASSERTL0(false, "This type of FwdTrans is not defined"
                                "for this expansion type");
        }

        switch (m_base[0]->GetPointsType())
        {
            case LibUtilities::eGaussGaussLegendre:
            case LibUtilities::eGaussGaussChebyshev:
            case LibUtilities::eGaussKronrodLegendre:
            {
                hasEndPoints = false;
            }
            break;
            case LibUtilities::eGaussLegendreWithMP:
            case LibUtilities::eGaussLobattoLegendre:
            case LibUtilities::eGaussLobattoChebyshev:
            case LibUtilities::eGaussLobattoKronrodLegendre:
            case LibUtilities::ePolyEvenlySpaced:
            case LibUtilities::eFourierEvenlySpaced:
            {
                hasEndPoints = true;
            }
            break;
            default:
                ASSERTL0(false, "FwdTransBndConstrained cannot be used "
                                "with this point type");
        }

        fill(outarray.data(), outarray.data() + m_ncoeffs, 0.0);

        if (hasEndPoints && hasEndModes)
        {

            outarray[GetVertexMap(0)] = inarray[0];
            outarray[GetVertexMap(1)] = inarray[m_base[0]->GetNumPoints() - 1];

            if (m_ncoeffs > 2)
            {
                // ideally, we would like to have tmp0 to be replaced
                // by outarray (currently MassMatrixOp does not allow
                // aliasing)
                Array<OneD, NekDouble> tmp0(m_ncoeffs);
                Array<OneD, NekDouble> tmp1(m_ncoeffs);

                StdRegions::StdMatrixKey stdmasskey(StdRegions::eMass,
                                                    DetShapeType(), *this);
                MassMatrixOp(outarray, tmp0, stdmasskey);
                v_IProductWRTBase(inarray, tmp1);

                Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                // get Mass matrix inverse (only of interior DOF)
                MatrixKey masskey(StdRegions::eMass, DetShapeType(), *this);
                DNekScalMatSharedPtr matsys =
                    (m_staticCondMatrixManager[masskey])->GetBlock(1, 1);

                Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, matsys->Scale(),
                            &((matsys->GetOwnedMatrix())->GetPtr())[0],
                            nInteriorDofs, tmp1.data() + offset, 1, 0.0,
                            outarray.data() + offset, 1);
            }
        }
        else
        {
            v_FwdTrans(inarray, outarray);
        }
    }
}

//-----------------------------
// Inner product functions
//-----------------------------
void SegExp::v_IProductWRTDerivBase(const int dir,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    ASSERTL1(dir < 3, "input dir is out of range");
    ASSERTL1((dir == 2) ? m_geom->GetCoordim() == 3 : true,
             "input dir is out of range");

    int nquad                                = m_base[0]->GetNumPoints();
    const Array<TwoD, const NekDouble> &gmat = m_geomFactors->GetDerivFactors();

    Array<OneD, NekDouble> tmp1(nquad);
    const bool Deformed =
        m_geomFactors->GetGtype() == SpatialDomains::eDeformed;

    if (Deformed)
    {
        Vmath::Vmul(nquad, gmat[dir], 1, inarray, 1, tmp1, 1);
    }
    else
    {
        Vmath::Smul(nquad, gmat[dir][0], inarray, 1, tmp1, 1);
    }

    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), tmp1, outarray, jac,
                            Deformed);
}

void SegExp::v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx,
                                         const Array<OneD, const NekDouble> &Fy,
                                         Array<OneD, NekDouble> &outarray)
{
    int nq = m_base[0]->GetNumPoints();
    Array<OneD, NekDouble> Fn(nq);

    // @TODO: This routine no longer makes sense as a normal is not unique to an
    // edge
    const Array<OneD, const Array<OneD, NekDouble>> &normals =
        GetLeftAdjacentElementExp()->GetTraceNormal(
            GetLeftAdjacentElementTrace());
    Vmath::Vmul(nq, &Fx[0], 1, &normals[0][0], 1, &Fn[0], 1);
    Vmath::Vvtvp(nq, &Fy[0], 1, &normals[1][0], 1, &Fn[0], 1, &Fn[0], 1);

    v_IProductWRTBase(Fn, outarray);
}

void SegExp::v_NormVectorIProductWRTBase(
    const Array<OneD, const Array<OneD, NekDouble>> &Fvec,
    Array<OneD, NekDouble> &outarray)
{
    NormVectorIProductWRTBase(Fvec[0], Fvec[1], outarray);
}

//-----------------------------
// Evaluation functions
//-----------------------------
NekDouble SegExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    Array<OneD, NekDouble> Lcoord(1);
    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);
    return StdSegExp::v_PhysEvalFirstDeriv(Lcoord, inarray, firstOrderDerivs);
}

NekDouble SegExp::v_PhysEvalFirstSecondDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs,
    std::array<NekDouble, 6> &secondOrderDerivs)
{
    Array<OneD, NekDouble> Lcoord(1);
    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);
    return StdSegExp::v_PhysEvalFirstSecondDeriv(
        Lcoord, inarray, firstOrderDerivs, secondOrderDerivs);
}

void SegExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                        Array<OneD, NekDouble> &coords)
{
    int i;

    ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[0] <= 1.0,
             "Local coordinates are not in region [-1,1]");

    m_geom->FillGeom();
    for (i = 0; i < m_geom->GetCoordim(); ++i)
    {
        coords[i] = m_geom->GetCoord(i, Lcoords);
    }
}

void SegExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                         Array<OneD, NekDouble> &coords_1,
                         Array<OneD, NekDouble> &coords_2)
{
    Expansion::v_GetCoords(coords_0, coords_1, coords_2);
}

// Get vertex value from the 1D Phys space.
void SegExp::v_GetVertexPhysVals(const int vertex,
                                 const Array<OneD, const NekDouble> &inarray,
                                 NekDouble &outarray)
{
    int nquad = m_base[0]->GetNumPoints();

    if (m_base[0]->GetPointsType() != LibUtilities::eGaussGaussLegendre)
    {
        switch (vertex)
        {
            case 0:
                outarray = inarray[0];
                break;
            case 1:
                outarray = inarray[nquad - 1];
                break;
        }
    }
    else
    {
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorGaussVertex] = vertex;

        StdRegions::StdMatrixKey key(StdRegions::eInterpGauss, DetShapeType(),
                                     *this, factors);

        DNekScalMatSharedPtr mat_gauss = m_matrixManager[key];

        outarray =
            Vmath::Dot(nquad, mat_gauss->GetOwnedMatrix()->GetPtr().data(), 1,
                       &inarray[0], 1);
    }
}

// Get vertex value from the 1D Phys space.
void SegExp::v_GetTracePhysVals(
    const int edge,
    [[maybe_unused]] const StdRegions::StdExpansionSharedPtr &EdgeExp,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] StdRegions::Orientation orient)
{
    NekDouble result;
    v_GetVertexPhysVals(edge, inarray, result);
    outarray[0] = result;
}

// Get vertex map from the 1D Phys space.
void SegExp::v_GetTracePhysMap(const int vertex, Array<OneD, int> &map)
{
    int nquad = m_base[0]->GetNumPoints();

    ASSERTL1(vertex == 0 || vertex == 1, "Vertex value should be 0 or 1");

    map = Array<OneD, int>(1);

    map[0] = vertex == 0 ? 0 : nquad - 1;
}

//-----------------------------
// Helper functions
//-----------------------------

void SegExp::v_SetCoeffsToOrientation(StdRegions::Orientation dir,
                                      Array<OneD, const NekDouble> &inarray,
                                      Array<OneD, NekDouble> &outarray)
{

    if (dir == StdRegions::eBackwards)
    {
        if (&inarray[0] != &outarray[0])
        {
            Array<OneD, NekDouble> intmp(inarray);
            ReverseCoeffsAndSign(intmp, outarray);
        }
        else
        {
            ReverseCoeffsAndSign(inarray, outarray);
        }
    }
}

StdRegions::StdExpansionSharedPtr SegExp::v_GetStdExp() const
{
    return MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey());
}

StdRegions::StdExpansionSharedPtr SegExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());

    return MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bkey0);
}

const Array<OneD, const NekDouble> &SegExp::v_GetPhysNormals(void)
{
    NEKERROR(ErrorUtil::efatal, "Got to SegExp");
    return NullNekDouble1DArray;
}

int SegExp::v_NumBndryCoeffs() const
{
    return 2;
}

int SegExp::v_NumDGBndryCoeffs() const
{
    return 2;
}

/// Unpack data from input file assuming it comes from
// the same expansion type
void SegExp::v_ExtractDataToCoeffs(
    const NekDouble *data, const std::vector<unsigned int> &nummodes,
    const int mode_offset, NekDouble *coeffs,
    [[maybe_unused]] std::vector<LibUtilities::BasisType> &fromType)
{
    switch (m_base[0]->GetBasisType())
    {
        case LibUtilities::eModified_A:
        {
            int fillorder = min((int)nummodes[mode_offset], m_ncoeffs);

            Vmath::Zero(m_ncoeffs, coeffs, 1);
            Vmath::Vcopy(fillorder, &data[0], 1, &coeffs[0], 1);
        }
        break;
        case LibUtilities::eGLL_Lagrange:
        {
            // Assume that input is also Gll_Lagrange
            // but no way to check;
            LibUtilities::PointsKey f0(nummodes[mode_offset],
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::PointsKey t0(m_base[0]->GetNumModes(),
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::Interp1D(f0, data, t0, coeffs);
        }
        break;
        case LibUtilities::eGauss_Lagrange:
        {
            // Assume that input is also Gauss_Lagrange
            // but no way to check;
            LibUtilities::PointsKey f0(nummodes[mode_offset],
                                       LibUtilities::eGaussGaussLegendre);
            LibUtilities::PointsKey t0(m_base[0]->GetNumModes(),
                                       LibUtilities::eGaussGaussLegendre);
            LibUtilities::Interp1D(f0, data, t0, coeffs);
        }
        break;
        default:
            ASSERTL0(false, "basis is either not set up or not hierarchicial");
    }
}

void SegExp::v_ComputeTraceNormal(const int vertex)
{
    int i;
    SpatialDomains::GeomType type            = m_geomFactors->GetGtype();
    const Array<TwoD, const NekDouble> &gmat = m_geomFactors->GetDerivFactors();
    int nqe                                  = 1;
    int vCoordDim                            = GetCoordim();

    m_traceNormals[vertex] = Array<OneD, Array<OneD, NekDouble>>(vCoordDim);
    Array<OneD, Array<OneD, NekDouble>> &normal = m_traceNormals[vertex];
    for (i = 0; i < vCoordDim; ++i)
    {
        normal[i] = Array<OneD, NekDouble>(nqe);
    }

    size_t nqb                     = nqe;
    size_t nbnd                    = vertex;
    m_elmtBndNormDirElmtLen[nbnd]  = Array<OneD, NekDouble>{nqb, 0.0};
    Array<OneD, NekDouble> &length = m_elmtBndNormDirElmtLen[nbnd];

    // Regular geometry case
    if ((type == SpatialDomains::eRegular) ||
        (type == SpatialDomains::eMovingRegular))
    {
        NekDouble vert;
        // Set up normals
        switch (vertex)
        {
            case 0:
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Fill(nqe, -gmat[i][0], normal[i], 1);
                }
                break;
            case 1:
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Fill(nqe, gmat[i][0], normal[i], 1);
                }
                break;
            default:
                ASSERTL0(false, "point is out of range (point < 2)");
        }

        // normalise
        vert = 0.0;
        for (i = 0; i < vCoordDim; ++i)
        {
            vert += normal[i][0] * normal[i][0];
        }
        vert = 1.0 / sqrt(vert);

        Vmath::Fill(nqb, vert, length, 1);

        for (i = 0; i < vCoordDim; ++i)
        {
            Vmath::Smul(nqe, vert, normal[i], 1, normal[i], 1);
        }
    }
}

//-----------------------------
// Operator creation functions
//-----------------------------

void SegExp::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] const StdRegions::StdMatrixKey &mkey)
{
    int nquad                                = m_base[0]->GetNumPoints();
    const Array<TwoD, const NekDouble> &gmat = m_geomFactors->GetDerivFactors();

    Array<OneD, NekDouble> physValues(nquad);
    Array<OneD, NekDouble> dPhysValuesdx(nquad);

    BwdTrans(inarray, physValues);

    // Laplacian matrix operation
    switch (m_geom->GetCoordim())
    {
        case 1:
        {
            v_PhysDeriv(physValues, dPhysValuesdx, NullNekDouble1DArray,
                        NullNekDouble1DArray);

            // multiply with the proper geometric factors
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 2:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy);

            // multiply with the proper geometric factors
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 3:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);
            Array<OneD, NekDouble> dPhysValuesdz(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy, dPhysValuesdz);

            // multiply with the proper geometric factors
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[2][0], 1, dPhysValuesdz.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        default:
            ASSERTL0(false, "Wrong number of dimensions");
            break;
    }

    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    const bool Deformed =
        m_geomFactors->GetGtype() == SpatialDomains::eDeformed;
    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), dPhysValuesdx, outarray,
                            jac, Deformed);
}

void SegExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void SegExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    int nquad                                = m_base[0]->GetNumPoints();
    const Array<TwoD, const NekDouble> &gmat = m_geomFactors->GetDerivFactors();
    const NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);

    Array<OneD, NekDouble> physValues(nquad);
    Array<OneD, NekDouble> dPhysValuesdx(nquad);
    Array<OneD, NekDouble> wsp(m_ncoeffs);

    BwdTrans(inarray, physValues);

    // mass matrix operation
    v_IProductWRTBase(physValues, wsp);

    // Laplacian matrix operation
    switch (m_geom->GetCoordim())
    {
        case 1:
        {
            PhysDeriv(physValues, dPhysValuesdx);

            // multiply with the proper geometric factors
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 2:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy);

            // multiply with the proper geometric factors
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 3:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);
            Array<OneD, NekDouble> dPhysValuesdz(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy, dPhysValuesdz);

            // multiply with the proper geometric factors
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[2][0], 1, dPhysValuesdz.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        default:
            ASSERTL0(false, "Wrong number of dimensions");
            break;
    }

    const Array<OneD, const NekDouble> &jac = m_geomFactors->GetJac();
    const bool Deformed =
        m_geomFactors->GetGtype() == SpatialDomains::eDeformed;
    v_IProductWRTBaseKernel(m_base[0]->GetDbdata(), dPhysValuesdx, outarray,
                            jac, Deformed);
    Blas::Daxpy(m_ncoeffs, lambda, wsp.data(), 1, outarray.data(), 1);
}

//-----------------------------
// Matrix creation functions
//-----------------------------

DNekScalBlkMatSharedPtr SegExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}

void SegExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
{
    m_staticCondMatrixManager.DeleteObject(mkey);
}

DNekScalMatSharedPtr SegExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}

void SegExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

DNekMatSharedPtr SegExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey = m_base[0]->GetBasisKey();
    StdRegions::StdSegExpSharedPtr tmp =
        MemoryManager<StdSegExp>::AllocateSharedPtr(bkey);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr SegExp::CreateMatrix(const MatrixKey &mkey)
{
    DNekScalMatSharedPtr returnval;
    NekDouble fac;
    LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();

    ASSERTL2(m_geomFactors->GetGtype() != SpatialDomains::eNoGeomType,
             "Geometric information is not set up");

    switch (mkey.GetMatrixType())
    {
        case StdRegions::eMass:
        {
            if ((m_geomFactors->GetGtype() == SpatialDomains::eDeformed) ||
                (mkey.GetNVarCoeff()))
            {
                fac = 1.0;
                goto UseLocRegionsMatrix;
            }
            else
            {
                fac = (m_geomFactors->GetJac())[0];
                goto UseStdRegionsMatrix;
            }
        }
        break;
        case StdRegions::eInvMass:
        {
            if ((m_geomFactors->GetGtype() == SpatialDomains::eDeformed) ||
                (mkey.GetNVarCoeff()))
            {
                NekDouble one = 1.0;
                StdRegions::StdMatrixKey masskey(StdRegions::eMass,
                                                 DetShapeType(), *this);
                DNekMatSharedPtr mat = GenMatrix(masskey);
                mat->Invert();

                returnval =
                    MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
            }
            else
            {
                fac = 1.0 / (m_geomFactors->GetJac())[0];
                goto UseStdRegionsMatrix;
            }
        }
        break;
        case StdRegions::eWeakDeriv0:
        case StdRegions::eWeakDeriv1:
        case StdRegions::eWeakDeriv2:
        {
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed ||
                mkey.GetNVarCoeff())
            {
                fac = 1.0;
                goto UseLocRegionsMatrix;
            }
            else
            {
                int dir = 0;
                switch (mkey.GetMatrixType())
                {
                    case StdRegions::eWeakDeriv0:
                        dir = 0;
                        break;
                    case StdRegions::eWeakDeriv1:
                        ASSERTL1(m_geom->GetCoordim() >= 2,
                                 "Cannot call eWeakDeriv2 in a "
                                 "coordinate system which is not at "
                                 "least two-dimensional");
                        dir = 1;
                        break;
                    case StdRegions::eWeakDeriv2:
                        ASSERTL1(m_geom->GetCoordim() == 3,
                                 "Cannot call eWeakDeriv2 in a "
                                 "coordinate system which is not "
                                 "three-dimensional");
                        dir = 2;
                        break;
                    default:
                        break;
                }

                MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                    mkey.GetShapeType(), *this);

                DNekMatSharedPtr WeakDerivStd = GetStdMatrix(deriv0key);
                fac = m_geomFactors->GetDerivFactors()[dir][0] *
                      m_geomFactors->GetJac()[0];

                returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                    fac, WeakDerivStd);
            }
        }
        break;
        case StdRegions::eLaplacian:
        {
            if (m_geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                fac = 1.0;
                goto UseLocRegionsMatrix;
            }
            else
            {
                int coordim = m_geom->GetCoordim();
                fac         = 0.0;
                for (int i = 0; i < coordim; ++i)
                {
                    fac += m_geomFactors->GetDerivFactors()[i][0] *
                           m_geomFactors->GetDerivFactors()[i][0];
                }
                fac *= m_geomFactors->GetJac()[0];
                goto UseStdRegionsMatrix;
            }
        }
        break;
        case StdRegions::eLinearAdvection:
        {
            if ((m_geomFactors->GetGtype() == SpatialDomains::eDeformed) ||
                (mkey.GetNVarCoeff()))
            {
                fac = 1.0;
                goto UseLocRegionsMatrix;
            }
            else
            {
                fac = (m_geomFactors->GetJac())[0];
                goto UseStdRegionsMatrix;
            }
        }
        break;
        case StdRegions::eHelmholtz:
        {
            NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
            MatrixKey masskey(StdRegions::eMass, mkey.GetShapeType(), *this);
            DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
            MatrixKey lapkey(StdRegions::eLaplacian, mkey.GetShapeType(), *this,
                             mkey.GetConstFactors(), mkey.GetVarCoeffs());
            DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

            int rows = LapMat.GetRows();
            int cols = LapMat.GetColumns();

            DNekMatSharedPtr helm =
                MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

            NekDouble one = 1.0;
            (*helm)       = LapMat + factor * MassMat;

            returnval =
                MemoryManager<DNekScalMat>::AllocateSharedPtr(one, helm);
        }
        break;
        case StdRegions::eLinearAdvectionReaction:
        {
            NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);

            // Construct mass matrix
            // Check for mass-specific varcoeffs to avoid unncessary
            // re-computation of the elemental matrix every time step
            StdRegions::VarCoeffMap massVarcoeffs = StdRegions::NullVarCoeffMap;
            if (mkey.HasVarCoeff(StdRegions::eVarCoeffMass))
            {
                massVarcoeffs[StdRegions::eVarCoeffMass] =
                    mkey.GetVarCoeff(StdRegions::eVarCoeffMass);
            }
            MatrixKey masskey(StdRegions::eMass, mkey.GetShapeType(), *this,
                              mkey.GetConstFactors(), massVarcoeffs);
            DNekScalMat &MassMat = *GetLocMatrix(masskey);

            // Construct advection matrix
            // Check for varcoeffs not required;
            // assume advection velocity is always time-dependent
            MatrixKey advkey(mkey, StdRegions::eLinearAdvection);
            DNekScalMat &AdvMat = *GetLocMatrix(advkey);

            int rows = MassMat.GetRows();
            int cols = MassMat.GetColumns();

            DNekMatSharedPtr adr =
                MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

            NekDouble one = 1.0;
            (*adr)        = -lambda * MassMat + AdvMat;

            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, adr);

            // Clear memory for time-dependent matrices
            DropLocMatrix(advkey);
            if (!massVarcoeffs.empty())
            {
                DropLocMatrix(masskey);
            }
        }
        break;
        case StdRegions::eLinearAdvectionDiffusionReaction:
        {
            NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);

            // Construct mass matrix
            // Check for mass-specific varcoeffs to avoid unncessary
            // re-computation of the elemental matrix every time step
            StdRegions::VarCoeffMap massVarcoeffs = StdRegions::NullVarCoeffMap;
            if (mkey.HasVarCoeff(StdRegions::eVarCoeffMass))
            {
                massVarcoeffs[StdRegions::eVarCoeffMass] =
                    mkey.GetVarCoeff(StdRegions::eVarCoeffMass);
            }
            MatrixKey masskey(StdRegions::eMass, mkey.GetShapeType(), *this,
                              mkey.GetConstFactors(), massVarcoeffs);
            DNekScalMat &MassMat = *GetLocMatrix(masskey);

            // Construct laplacian matrix (Check for varcoeffs)
            // Take all varcoeffs if one or more are detected
            // TODO We might want to have a map
            // from MatrixType to Vector of Varcoeffs and vice-versa
            StdRegions::VarCoeffMap lapVarcoeffs = StdRegions::NullVarCoeffMap;
            if ((mkey.HasVarCoeff(StdRegions::eVarCoeffLaplacian)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD00)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD01)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD10)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD02)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD20)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD11)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD12)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD21)) ||
                (mkey.HasVarCoeff(StdRegions::eVarCoeffD22)))
            {
                lapVarcoeffs = mkey.GetVarCoeffs();
            }
            MatrixKey lapkey(StdRegions::eLaplacian, mkey.GetShapeType(), *this,
                             mkey.GetConstFactors(), lapVarcoeffs);
            DNekScalMat &LapMat = *GetLocMatrix(lapkey);

            // Construct advection matrix
            // Check for varcoeffs not required;
            // assume advection velocity is always time-dependent
            MatrixKey advkey(mkey, StdRegions::eLinearAdvection);
            DNekScalMat &AdvMat = *GetLocMatrix(advkey);

            int rows = LapMat.GetRows();
            int cols = LapMat.GetColumns();

            DNekMatSharedPtr adr =
                MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

            NekDouble one = 1.0;
            (*adr)        = LapMat - lambda * MassMat + AdvMat;

            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, adr);

            // Clear memory for time-dependent matrices
            DropLocMatrix(advkey);
            if (!massVarcoeffs.empty())
            {
                DropLocMatrix(masskey);
            }
            if (!lapVarcoeffs.empty())
            {
                DropLocMatrix(lapkey);
            }
        }
        break;
        case StdRegions::eHybridDGHelmholtz:
        case StdRegions::eHybridDGLamToU:
        case StdRegions::eHybridDGLamToQ0:
        case StdRegions::eHybridDGHelmBndLam:
        {
            NekDouble one = 1.0;

            DNekMatSharedPtr mat = GenMatrix(mkey);
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
        }
        break;
        case StdRegions::eInvHybridDGHelmholtz:
        {
            NekDouble one = 1.0;

            //                  StdRegions::StdMatrixKey
            //                  hkey(StdRegions::eHybridDGHelmholtz,
            //                                                DetShapeType(),*this,
            //                                                mkey.GetConstant(0),
            //                                                mkey.GetConstant(1));
            MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetShapeType(),
                           *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
            DNekMatSharedPtr mat = GenMatrix(hkey);

            mat->Invert();
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
        }
        break;
        case StdRegions::eInterpGauss:
        {
            DNekMatSharedPtr m_Ix;
            Array<OneD, NekDouble> coords(1, 0.0);
            StdRegions::ConstFactorMap factors = mkey.GetConstFactors();
            int vertex = (int)factors[StdRegions::eFactorGaussVertex];

            coords[0] = (vertex == 0) ? -1.0 : 1.0;

            m_Ix = m_base[0]->GetI(coords);
            returnval =
                MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, m_Ix);
        }
        break;

        UseLocRegionsMatrix:
        {
            DNekMatSharedPtr mat = GenMatrix(mkey);
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac, mat);
        }
        break;
        UseStdRegionsMatrix:
        {
            DNekMatSharedPtr mat = GetStdMatrix(mkey);
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac, mat);
        }
        break;
        default:
        {
            NekDouble one        = 1.0;
            DNekMatSharedPtr mat = GenMatrix(mkey);

            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
        }
        break;
    }

    return returnval;
}

DNekMatSharedPtr SegExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
            returnval = Expansion1D::v_GenMatrix(mkey);
            break;
        default:
            returnval = StdSegExp::v_GenMatrix(mkey);
            break;
    }

    return returnval;
}

//-----------------------------
// Private methods
//-----------------------------

///  Reverse the coefficients in a boundary interior expansion
///  this routine is of use when we need the segment
///  coefficients corresponding to a expansion in the reverse
///  coordinate direction
void SegExp::ReverseCoeffsAndSign(const Array<OneD, NekDouble> &inarray,
                                  Array<OneD, NekDouble> &outarray)
{

    int m;
    NekDouble sgn = 1;

    ASSERTL1(&inarray[0] != &outarray[0],
             "inarray and outarray can not be the same");
    switch (GetBasisType(0))
    {
        case LibUtilities::eModified_A:
            // Swap vertices
            outarray[0] = inarray[1];
            outarray[1] = inarray[0];
            // negate odd modes
            for (m = 2; m < m_ncoeffs; ++m)
            {
                outarray[m] = sgn * inarray[m];
                sgn         = -sgn;
            }
            break;
        case LibUtilities::eGLL_Lagrange:
        case LibUtilities::eGauss_Lagrange:
            for (m = 0; m < m_ncoeffs; ++m)
            {
                outarray[m_ncoeffs - 1 - m] = inarray[m];
            }
            break;
        default:
            ASSERTL0(false, "This basis is not allowed in this method");
            break;
    }
}

/* \brief Mass inversion product from \a inarray to \a outarray
 *
 *     Multiply by the inverse of the mass matrix
 *     \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$
 *
 **/
void SegExp::MultiplyByElmtInvMass(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray)
{
    // get Mass matrix inverse
    MatrixKey masskey(StdRegions::eInvMass, DetShapeType(), *this);
    DNekScalMatSharedPtr matsys = m_matrixManager[masskey];

    NekVector<NekDouble> in(m_ncoeffs, inarray, eCopy);
    NekVector<NekDouble> out(m_ncoeffs, outarray, eWrapper);

    out = (*matsys) * in;
}

} // namespace Nektar::LocalRegions
