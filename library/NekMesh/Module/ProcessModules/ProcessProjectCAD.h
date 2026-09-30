////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessProjectCAD.h
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
//  Description:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_PROCESSPROJECTCAD
#define UTILITIES_NEKMESH_PROCESSPROJECTCAD

#include <NekMesh/Module/Module.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>

#include <boost/geometry/index/rtree.hpp>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, unsigned> boxI;

namespace Nektar::NekMesh
{

class ProcessProjectCAD : public NekMesh::ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(NekMesh::MeshSharedPtr m)
    {
        return MemoryManager<ProcessProjectCAD>::AllocateSharedPtr(m);
    }
    static NekMesh::ModuleKey className;

    ProcessProjectCAD(NekMesh::MeshSharedPtr m);
    ~ProcessProjectCAD() override;

    /// Write mesh to output file.
    void Process() override;

    std::string GetModuleName() override
    {
        return "ProcessProjectCAD";
    }

private:
    // Min edge length in an element
    std::map<NodeSharedPtr, NekDouble> minConEdge; // should be unordered?!
    // Vertices on the surface
    NodeSet surfNodes;
    // Edges on the surface
    EdgeSet surfEdges;
    // Surface vertices to 3D elements
    std::map<NodeSharedPtr, std::vector<ElementSharedPtr>> surfNodeToEl;
    // this is a set of nodes which have a CAD failure
    // if touched in the HO stage they should be ignored and linearised
    NodeSet lockedNodes;

    void LoadCAD(std::string filename);
    void CreateBoundingBoxes(bgi::rtree<boxI, bgi::quadratic<16>> &rtree,
                             bgi::rtree<boxI, bgi::quadratic<16>> &rtreeCurve,
                             bgi::rtree<boxI, bgi::quadratic<16>> &rtreeNode,
                             NekDouble tolv1, NekDouble scale);
    bool IsNotValid(std::vector<NekMesh::ElementSharedPtr> &els);
    void CalculateMinEdgeLength();
    void Auxilaries();
    void LinkVertexToCAD(NekMesh::MeshSharedPtr &m_mesh, bool CADCurve,
                         NodeSet &lockedNodes, NekDouble tolv1, NekDouble tolv2,
                         bgi::rtree<boxI, bgi::quadratic<16>> &rtree,
                         bgi::rtree<boxI, bgi::quadratic<16>> &rtreeCurve,
                         bgi::rtree<boxI, bgi::quadratic<16>> &rtreeNode);

    bool FindAndProject(bgi::rtree<boxI, bgi::quadratic<16>> &rtree,
                        std::array<NekDouble, 3> &in, int &surf);
    void LinkEdgeToCAD(EdgeSet &surfEdges, NekDouble tolv1);
    void LinkFaceToCAD(NekDouble tolv1);

    // for CASE 3 elements between two surfaces
    void ProjectEdges(EdgeSet &surfEdges, int order,
                      bgi::rtree<boxI, bgi::quadratic<16>> &rtree);
    std::vector<int> IntersectCADSurf(std::vector<CADSurfSharedPtr> v1_CADs,
                                      std::vector<CADSurfSharedPtr> v2_CADs);

    std::vector<int> IntersectCADCurve(std::vector<CADCurveSharedPtr> v1_CADs,
                                       std::vector<CADCurveSharedPtr> v2_CADs);
    void LinkHOtoCAD(EdgeSet &surfEdges, NekDouble tolv1);

    void Diagnostics();
    void ExportCAD();
};
} // namespace Nektar::NekMesh

#endif
