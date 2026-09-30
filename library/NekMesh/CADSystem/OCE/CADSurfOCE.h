////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurfOCE.h
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
//  Description: CAD object surface.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMesh_CADSYSTEM_OCE_CADSURFOCE
#define NekMesh_CADSYSTEM_OCE_CADSURFOCE

#include <NekMesh/CADSystem/CADSurf.h>
#include <NekMesh/CADSystem/OCE/OpenCascade.h>

namespace Nektar::NekMesh
{

class CADSurfOCE : public CADSurf
{
public:
    static CADSurfSharedPtr create()
    {
        return MemoryManager<CADSurfOCE>::AllocateSharedPtr();
    }

    static std::string key;

    CADSurfOCE()
    {
    }

    ~CADSurfOCE() override
    {
    }

    void Initialise(int i, TopoDS_Shape in);

    std::array<NekDouble, 4> GetBounds() override;
    void GetBounds(NekDouble &umin, NekDouble &umax, NekDouble &vmin,
                   NekDouble &vmax) override;
    std::array<NekDouble, 3> N(std::array<NekDouble, 2> uv) override;
    std::array<NekDouble, 9> D1(std::array<NekDouble, 2> uv) override;
    std::array<NekDouble, 18> D2(std::array<NekDouble, 2> uv) override;
    std::array<NekDouble, 3> P(std::array<NekDouble, 2> uv) override;
    void P(std::array<NekDouble, 2> uv, NekDouble &x, NekDouble &y,
           NekDouble &z) override;
    std::array<NekDouble, 2> locuv(std::array<NekDouble, 3> p,
                                   NekDouble &dist) override;
    NekDouble Curvature(std::array<NekDouble, 2> uv) override;
    std::array<NekDouble, 6> BoundingBox() override;
    std::array<NekDouble, 6> BoundingBox(NekDouble scale) override;
    bool IsPlanar() override;
    std::array<NekDouble, 2> locuv(std::array<NekDouble, 3> p, NekDouble &dist,
                                   NekDouble Umin, NekDouble Usup,
                                   NekDouble Vmin, NekDouble Vsup) override;

private:
    /// Function which tests the the value of uv used is within the surface
    void Test(std::array<NekDouble, 2> uv) override;
    /// OpenCascade object for surface.
    Handle(Geom_Surface) m_s;
    /// parametric bounds
    std::array<NekDouble, 4> m_bounds;
    /// locuv object (stored because it gets faster with stored information)
    ShapeAnalysis_Surface *m_sas;
    /// original shape
    TopoDS_Shape m_shape;
    ///
    BRepTopAdaptor_FClass2d *m_2Dclass;
    /// True if we're a transfinite surface (used for Geo)
    bool m_isTransfiniteSurf;
};

} // namespace Nektar::NekMesh

#endif
