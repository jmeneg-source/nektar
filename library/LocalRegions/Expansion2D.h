///////////////////////////////////////////////////////////////////////////////
//
// File: Expansion2D.h
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
// Description: Header file for Expansion2D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION2D_H
#define EXPANSION2D_H

#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/MatrixKey.h>
#include <SpatialDomains/Geometry2D.h>
#include <StdRegions/StdExpansion2D.h>

namespace Nektar::LocalRegions
{
class Expansion3D;
typedef std::shared_ptr<Expansion3D> Expansion3DSharedPtr;
typedef std::weak_ptr<Expansion3D> Expansion3DWeakPtr;

class Expansion2D;
typedef std::shared_ptr<Expansion2D> Expansion2DSharedPtr;
typedef std::weak_ptr<Expansion2D> Expansion2DWeakPtr;
typedef std::vector<Expansion2DSharedPtr> Expansion2DVector;

class Expansion2D : virtual public Expansion,
                    virtual public StdRegions::StdExpansion2D
{
public:
    LOCAL_REGIONS_EXPORT Expansion2D(SpatialDomains::Geometry2D *pGeom);

    LOCAL_REGIONS_EXPORT ~Expansion2D() override = default;

    LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr
    CreateMatrix(const MatrixKey &mkey);

    LOCAL_REGIONS_EXPORT void SetTraceToGeomOrientation(
        Array<OneD, ExpansionSharedPtr> &EdgeExp,
        Array<OneD, NekDouble> &inout);

    LOCAL_REGIONS_EXPORT Array<OneD, unsigned int> GetTraceInverseBoundaryMap(
        int eid);

    inline void AddNormTraceInt(const int dir,
                                Array<OneD, ExpansionSharedPtr> &EdgeExp,
                                Array<OneD, Array<OneD, NekDouble>> &edgeCoeffs,
                                Array<OneD, NekDouble> &outarray);

    inline void AddNormTraceInt(const int dir,
                                Array<OneD, const NekDouble> &inarray,
                                Array<OneD, ExpansionSharedPtr> &EdgeExp,
                                Array<OneD, NekDouble> &outarray,
                                const StdRegions::VarCoeffMap &varcoeffs);

    inline void AddEdgeBoundaryInt(
        const int edge, ExpansionSharedPtr &EdgeExp,
        Array<OneD, NekDouble> &edgePhys, Array<OneD, NekDouble> &outarray,
        const StdRegions::VarCoeffMap &varcoeffs = StdRegions::NullVarCoeffMap);

    inline void AddHDGHelmholtzEdgeTerms(
        const NekDouble tau, const int edge,
        Array<OneD, ExpansionSharedPtr> &EdgeExp,
        Array<OneD, NekDouble> &edgePhys,
        const StdRegions::VarCoeffMap &dirForcing,
        Array<OneD, NekDouble> &outarray);

    inline void AddHDGHelmholtzTraceTerms(
        const NekDouble tau, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, ExpansionSharedPtr> &EdgeExp,
        const StdRegions::VarCoeffMap &dirForcing,
        Array<OneD, NekDouble> &outarray);

    inline SpatialDomains::Geometry2D *GetGeom2D() const;

    LOCAL_REGIONS_EXPORT void ReOrientEdgePhysMap(
        const int nvert, const StdRegions::Orientation orient, const int nq0,
        Array<OneD, int> &idmap);

    LOCAL_REGIONS_EXPORT void v_NormalTraceDerivFactors(
        Array<OneD, Array<OneD, NekDouble>> &factors,
        Array<OneD, Array<OneD, NekDouble>> &d0factors,
        Array<OneD, Array<OneD, NekDouble>> &d1factors) override;

protected:
    std::vector<bool> m_requireNeg;

    //----------------------------
    // Differentiation Methods
    //----------------------------
    LOCAL_REGIONS_EXPORT void v_PhysDeriv(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_PhysDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d0, Array<OneD, NekDouble> &out_d1,
        Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray) override;
    LOCAL_REGIONS_EXPORT void v_PhysDirectionalDeriv(
        const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const NekDouble> &direction,
        Array<OneD, NekDouble> &out) override;

    LOCAL_REGIONS_EXPORT void v_IProductWRTBase(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    // Hybridized DG routines
    void v_DGDeriv(const int dir, const Array<OneD, const NekDouble> &incoeffs,
                   Array<OneD, ExpansionSharedPtr> &EdgeExp,
                   Array<OneD, Array<OneD, NekDouble>> &edgeCoeffs,
                   Array<OneD, NekDouble> &out_d) override;

    DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey) override;

    void v_GenTraceExp(const int traceid, ExpansionSharedPtr &exp) override;

    void v_AddEdgeNormBoundaryInt(const int edge,
                                  const ExpansionSharedPtr &EdgeExp,
                                  const Array<OneD, const NekDouble> &Fx,
                                  const Array<OneD, const NekDouble> &Fy,
                                  Array<OneD, NekDouble> &outarray) override;

    void v_AddEdgeNormBoundaryInt(const int edge,
                                  const ExpansionSharedPtr &EdgeExp,
                                  const Array<OneD, const NekDouble> &Fn,
                                  Array<OneD, NekDouble> &outarray) override;

    void v_AddRobinMassMatrix(const int edgeid,
                              const Array<OneD, const NekDouble> &primCoeffs,
                              DNekMatSharedPtr &inoutmat) override;

    void v_AddRobinTraceContribution(
        const int traceid, const Array<OneD, const NekDouble> &primCoeffs,
        const Array<OneD, NekDouble> &incoeffs,
        Array<OneD, NekDouble> &coeffs) override;

    DNekMatSharedPtr v_BuildVertexMatrix(
        const DNekScalMatSharedPtr &r_bnd) override;

    void v_ReOrientTracePhysVals(const StdRegions::Orientation orient,
                                 const Array<OneD, const NekDouble> &in,
                                 Array<OneD, NekDouble> &out, const int nq0,
                                 const int nq1, bool Forwards) override;

    void v_SetUpPhysNormals(const int edge) override;
    NekDouble v_VectorFlux(
        const Array<OneD, Array<OneD, NekDouble>> &vec) override;
    void v_TraceNormLen(const int traceid, NekDouble &h, NekDouble &p) override;

private:
    void GetPhysEdgeVarCoeffsFromElement(
        const int edge, ExpansionSharedPtr &EdgeExp,
        const Array<OneD, const NekDouble> &varcoeff,
        Array<OneD, NekDouble> &outarray);

    Array<OneD, NekDouble> GetnEdgecdotMF(
        const int dir, const int edge, ExpansionSharedPtr &EdgeExp_e,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const StdRegions::VarCoeffMap &varcoeffs);
};

inline SpatialDomains::Geometry2D *Expansion2D::GetGeom2D() const
{
    return static_cast<SpatialDomains::Geometry2D *>(m_geom);
}
} // namespace Nektar::LocalRegions

#endif
