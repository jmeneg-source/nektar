////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurveOCE.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: CAD object curve methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMesh/CADSystem/OCE/CADCurveOCE.h>
#include <array>

using namespace std;

namespace Nektar::NekMesh
{

std::string CADCurveOCE::key = GetCADCurveFactory().RegisterCreatorFunction(
    "oce", CADCurveOCE::create, "CADCurveOCE");

void CADCurveOCE::Initialise(int i, TopoDS_Shape in)
{
    m_occEdge = TopoDS::Edge(in);

    GProp_GProps System;
    BRepGProp::LinearProperties(m_occEdge, System);
    m_length = System.Mass() / 1000.0;

    m_c  = BRep_Tool::Curve(TopoDS::Edge(in), m_b[0], m_b[1]);
    m_id = i;
}

NekDouble CADCurveOCE::tAtArcLength(NekDouble s)
{
    GeomAdaptor_Curve c(m_c);
    GCPnts_AbscissaPoint ap(c, s * 1000.0, m_b[0]);
    return ap.Parameter();
}

NekDouble CADCurveOCE::Length(NekDouble ti, NekDouble tf)
{
    Handle(Geom_Curve) NewCurve = new Geom_TrimmedCurve(m_c, ti, tf);
    TopoDS_Edge NewEdge         = BRepBuilderAPI_MakeEdge(NewCurve);
    GProp_GProps System;
    BRepGProp::LinearProperties(NewEdge, System);
    return System.Mass() / 1000.0;
}

NekDouble CADCurveOCE::GetMinDistance(std::array<NekDouble, 3> &xyz)
{
    gp_Pnt loc(xyz[0] * 1000.0, xyz[1] * 1000.0, xyz[2] * 1000.0);
    GeomAPI_ProjectPointOnCurve proj(loc, m_c, m_b[0], m_b[1]);
    if (proj.NbPoints())
    {
        return proj.LowerDistance() / 1000.0;
    }
    else
    {
        return std::min(loc.Distance(m_c->Value(m_b[0])),
                        loc.Distance(m_c->Value(m_b[1]))) /
               1000.0;
    }
}

NekDouble CADCurveOCE::loct(std::array<NekDouble, 3> xyz, NekDouble &t)
{
    t = 0.0;
    gp_Pnt loc(xyz[0] * 1000.0, xyz[1] * 1000.0, xyz[2] * 1000.0);
    ShapeAnalysis_Curve sac;
    gp_Pnt p;
    sac.Project(m_c, loc, Precision::Confusion(), p, t);
    return p.Distance(loc) / 1000.0;
}

// This function updates the parametric coordinate t of a cartesian point xyz
// within the bounds of the CADCurve cf and cl ; Fixes the bug with respect to
// Periodic CADCurves like Circles
NekDouble CADCurveOCE::loct(std::array<NekDouble, 3> xyz, NekDouble &t,
                            NekDouble cf, NekDouble cl)
{
    t = 0.0;

    gp_Pnt loc(xyz[0] * 1000.0, xyz[1] * 1000.0, xyz[2] * 1000.0);

    ShapeAnalysis_Curve sac;
    gp_Pnt p;
    sac.Project(m_c, loc, Precision::Confusion(), p, t, cf, cl);

    return p.Distance(loc) / 1000.0;
}

std::array<NekDouble, 3> CADCurveOCE::P(NekDouble t)
{
    gp_Pnt loc = m_c->Value(t);
    return {loc.X() / 1000.0, loc.Y() / 1000.0, loc.Z() / 1000.0};
}

void CADCurveOCE::P(NekDouble t, NekDouble &x, NekDouble &y, NekDouble &z)
{
    gp_Pnt loc = m_c->Value(t);
    x          = loc.X() / 1000.0;
    y          = loc.Y() / 1000.0;
    z          = loc.Z() / 1000.0;
}

std::array<NekDouble, 9> CADCurveOCE::D2(NekDouble t)
{
    gp_Pnt loc;
    gp_Vec d1, d2;
    m_c->D2(t, loc, d1, d2);

    return {loc.X() / 1000.0, loc.Y() / 1000.0, loc.Z() / 1000.0,
            d1.X() / 1000.0,  d1.Y() / 1000.0,  d1.Z() / 1000.0,
            d2.X() / 1000.0,  d2.Y() / 1000.0,  d2.Z() / 1000.0};
}

std::array<NekDouble, 3> CADCurveOCE::N(NekDouble t)
{
    GeomLProp_CLProps d(m_c, 2, Precision::Confusion());
    d.SetParameter(t + 1e-8);

    gp_Vec d2 = d.D2();
    if (d2.Magnitude() < 1e-8)
    {
        // no normal, stright line
        return {0.0, 0.0, 0.0};
    }

    gp_Dir n;
    d.Normal(n);

    return {n.X(), n.Y(), n.Z()};
}

NekDouble CADCurveOCE::Curvature(NekDouble t)
{
    GeomLProp_CLProps d(m_c, 2, Precision::Confusion());
    d.SetParameter(t);

    return d.Curvature() * 1000.0;
}

std::array<NekDouble, 2> CADCurveOCE::GetBounds()
{
    return m_b;
}

void CADCurveOCE::GetBounds(NekDouble &tmin, NekDouble &tmax)
{
    tmin = m_b[0];
    tmax = m_b[1];
}

std::array<NekDouble, 6> CADCurveOCE::GetMinMax()
{
    gp_Pnt start =
        BRep_Tool::Pnt(TopExp::FirstVertex(m_occEdge, Standard_True));
    gp_Pnt end = BRep_Tool::Pnt(TopExp::LastVertex(m_occEdge, Standard_True));

    return {start.X() / 1000.0, start.Y() / 1000.0, start.Z() / 1000.0,
            end.X() / 1000.0,   end.Y() / 1000.0,   end.Z() / 1000.0};
}

std::array<NekDouble, 6> CADCurveOCE::BoundingBox(NekDouble scale)
{
    BRepMesh_IncrementalMesh brmsh(m_occEdge, 0.005 * scale);

    Bnd_Box B;
    BRepBndLib::Add(m_occEdge, B);
    NekDouble e = sqrt(B.SquareExtent()) * 0.01 * scale;
    e           = min(e, 5e-3 * scale);
    B.Enlarge(e);

    std::array<NekDouble, 6> ret;
    B.Get(ret[0], ret[1], ret[2], ret[3], ret[4], ret[5]);
    ret[0] /= 1000.0;
    ret[1] /= 1000.0;
    ret[2] /= 1000.0;
    ret[3] /= 1000.0;
    ret[4] /= 1000.0;
    ret[5] /= 1000.0;

    return ret;
}

} // namespace Nektar::NekMesh
