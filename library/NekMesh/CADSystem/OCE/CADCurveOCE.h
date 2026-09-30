////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurveOCE.h
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
//  Description: CAD object curve.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESH_CADSYSTEM_OCE_CADCURVEOCE
#define NEKMESH_CADSYSTEM_OCE_CADCURVEOCE

#include <NekMesh/CADSystem/CADCurve.h>
#include <NekMesh/CADSystem/OCE/OpenCascade.h>

namespace Nektar::NekMesh
{

class CADCurveOCE : public CADCurve
{
public:
    static CADCurveSharedPtr create()
    {
        return MemoryManager<CADCurveOCE>::AllocateSharedPtr();
    }

    static std::string key;

    CADCurveOCE()
    {
    }

    ~CADCurveOCE() override
    {
    }

    std::array<NekDouble, 2> GetBounds() override;
    void GetBounds(NekDouble &tmin, NekDouble &tmax) override;
    NekDouble Length(NekDouble ti, NekDouble tf) override;
    std::array<NekDouble, 3> P(NekDouble t) override;
    void P(NekDouble t, NekDouble &x, NekDouble &y, NekDouble &z) override;
    std::array<NekDouble, 9> D2(NekDouble t) override;
    NekDouble tAtArcLength(NekDouble s) override;
    std::array<NekDouble, 6> GetMinMax() override;
    NekDouble loct(std::array<NekDouble, 3> xyz, NekDouble &t) override;
    NekDouble GetMinDistance(std::array<NekDouble, 3> &xyz) override;
    NekDouble Curvature(NekDouble t) override;
    std::array<NekDouble, 3> N(NekDouble t) override;
    std::array<NekDouble, 6> BoundingBox(NekDouble scale) override;

    void Initialise(int i, TopoDS_Shape in);

    NekDouble loct(std::array<NekDouble, 3> xyz, NekDouble &t, NekDouble cf,
                   NekDouble cl) override;

private:
    /// OpenCascade edge
    TopoDS_Edge m_occEdge;
    /// object used for reverse lookups
    Handle(Geom_Curve) m_c;
    /// store the parametric bounds of the curve
    std::array<NekDouble, 2> m_b;
};
} // namespace Nektar::NekMesh

#endif
