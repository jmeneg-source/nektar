////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessProjectCAD.cpp
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

#include "ProcessProjectCAD.h"
#include <NekMesh/MeshElements/Element.h>

#include <NekMesh/CADSystem/CADCurve.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar::NekMesh
{

const NekDouble prismU1[6] = {-1.0, 1.0, 1.0, -1.0, -1.0, -1.0};
const NekDouble prismV1[6] = {-1.0, -1.0, 1.0, 1.0, -1.0, 1.0};
const NekDouble prismW1[6] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0};

const int pyramidV0[4] = {0, 1, 2, 3};
const int pyramidV1[4] = {1, 2, 3, 0};
const int pyramidV2[4] = {3, 0, 1, 2};
const int pyramidV3[4] = {4, 4, 4, 4};

ModuleKey ProcessProjectCAD::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "projectcad"), ProcessProjectCAD::create,
        "Projects mesh to CAD");

ProcessProjectCAD::ProcessProjectCAD(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["file"]  = ConfigOption(false, "", "CAD file");
    m_config["order"] = ConfigOption(false, "4", "Enforce a polynomial order");
    m_config["surfopti"] = ConfigOption(true, "1", "Run HO-Surface Module");
    m_config["varopti"] =
        ConfigOption(false, "0", "Run the Variational Optmiser");
    m_config["cLength"] =
        ConfigOption(false, "1", "Characteristic Length CAD-Reconstruction");
    m_config["tolv1"] =
        ConfigOption(false, "1e-6",
                     "(Optional) min distance of initial Vertex to CADSurface");
    m_config["tolv2"] = ConfigOption(
        false, "1e-5",
        " (Optional) max distance of initial Vertex to CADSurface");
    m_config["ho"] =
        ConfigOption(false, "", "Pass when the input is already HO-mesh.");
    m_config["extract"] = ConfigOption(false, "", "Export the CAD to cad.txt.");
}

ProcessProjectCAD::~ProcessProjectCAD()
{
}

bool ProcessProjectCAD::FindAndProject(
    bgi::rtree<boxI, bgi::quadratic<16>> &rtree, std::array<NekDouble, 3> &in,
    [[maybe_unused]] int &surf)
{
    point q(in[0], in[1], in[2]);
    vector<boxI> result;
    rtree.query(bgi::intersects(q), back_inserter(result));

    if (result.size() == 0)
    {
        // along a projecting edge the node is too far from any surface boxes
        // this is hardly surprising but rare, return false and linearise
        m_log(VERBOSE) << "FindAndProject() - Edge node not in any of the "
                          "Bounding boxes - No projection! "
                       << endl;
        return false;
    }

    int minsurf       = 0;
    NekDouble minDist = numeric_limits<double>::max();
    NekDouble dist;

    for (int j = 0; j < result.size(); j++)
    {
        m_mesh->m_cad->GetSurf(result[j].second)->locuv(in, dist);

        if (dist < minDist)
        {
            minDist = dist;
            minsurf = result[j].second;
        }
    }

    auto uv = m_mesh->m_cad->GetSurf(minsurf)->locuv(in, dist);

    in = m_mesh->m_cad->GetSurf(minsurf)->P(uv);

    return true;
}

bool ProcessProjectCAD::IsNotValid(vector<ElementSharedPtr> &els)
{
    // short algebraic method to figure out the vailidy of elements
    // test the volume of tetrahedrons constituted by three sibling edges
    for (int i = 0; i < els.size(); i++)
    {
        if (els[i]->GetShapeType() == LibUtilities::ePrism)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            for (int j = 0; j < 6; j++)
            {
                NekDouble a2 = 0.5 * (1 + prismU1[j]);
                NekDouble b1 = 0.5 * (1 - prismV1[j]);
                NekDouble b2 = 0.5 * (1 + prismV1[j]);
                NekDouble c2 = 0.5 * (1 + prismW1[j]);
                NekDouble d  = 0.5 * (prismU1[j] + prismW1[j]);

                std::array<NekDouble, 9> jac;

                jac[0] = -0.5 * b1 * ns[0]->m_x + 0.5 * b1 * ns[1]->m_x +
                         0.5 * b2 * ns[2]->m_x - 0.5 * b2 * ns[3]->m_x;
                jac[1] = -0.5 * b1 * ns[0]->m_y + 0.5 * b1 * ns[1]->m_y +
                         0.5 * b2 * ns[2]->m_y - 0.5 * b2 * ns[3]->m_y;
                jac[2] = -0.5 * b1 * ns[0]->m_z + 0.5 * b1 * ns[1]->m_z +
                         0.5 * b2 * ns[2]->m_z - 0.5 * b2 * ns[3]->m_z;

                jac[3] = 0.5 * d * ns[0]->m_x - 0.5 * a2 * ns[1]->m_x +
                         0.5 * a2 * ns[2]->m_x - 0.5 * d * ns[3]->m_x -
                         0.5 * c2 * ns[4]->m_x + 0.5 * c2 * ns[5]->m_x;
                jac[4] = 0.5 * d * ns[0]->m_y - 0.5 * a2 * ns[1]->m_y +
                         0.5 * a2 * ns[2]->m_y - 0.5 * d * ns[3]->m_y -
                         0.5 * c2 * ns[4]->m_y + 0.5 * c2 * ns[5]->m_y;
                jac[5] = 0.5 * d * ns[0]->m_z - 0.5 * a2 * ns[1]->m_z +
                         0.5 * a2 * ns[2]->m_z - 0.5 * d * ns[3]->m_z -
                         0.5 * c2 * ns[4]->m_z + 0.5 * c2 * ns[5]->m_z;

                jac[6] = -0.5 * b1 * ns[0]->m_x - 0.5 * b2 * ns[3]->m_x +
                         0.5 * b1 * ns[4]->m_x + 0.5 * b2 * ns[5]->m_x;
                jac[7] = -0.5 * b1 * ns[0]->m_y - 0.5 * b2 * ns[3]->m_y +
                         0.5 * b1 * ns[4]->m_y + 0.5 * b2 * ns[5]->m_y;
                jac[8] = -0.5 * b1 * ns[0]->m_z - 0.5 * b2 * ns[3]->m_z +
                         0.5 * b1 * ns[4]->m_z + 0.5 * b2 * ns[5]->m_z;

                NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                               jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                               jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

                if (jc < NekConstants::kNekZeroTol)
                {
                    return true;
                }
            }
        }
        else if (els[i]->GetShapeType() == LibUtilities::ePyramid)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            for (int j = 0; j < 4; j++)
            {
                std::array<NekDouble, 9> jac;

                jac[0] = 0.5 * (ns[pyramidV1[j]]->m_x - ns[pyramidV0[j]]->m_x);
                jac[1] = 0.5 * (ns[pyramidV1[j]]->m_y - ns[pyramidV0[j]]->m_y);
                jac[2] = 0.5 * (ns[pyramidV1[j]]->m_z - ns[pyramidV0[j]]->m_z);
                jac[3] = 0.5 * (ns[pyramidV2[j]]->m_x - ns[pyramidV0[j]]->m_x);
                jac[4] = 0.5 * (ns[pyramidV2[j]]->m_y - ns[pyramidV0[j]]->m_y);
                jac[5] = 0.5 * (ns[pyramidV2[j]]->m_z - ns[pyramidV0[j]]->m_z);
                jac[6] = 0.5 * (ns[pyramidV3[j]]->m_x - ns[pyramidV0[j]]->m_x);
                jac[7] = 0.5 * (ns[pyramidV3[j]]->m_y - ns[pyramidV0[j]]->m_y);
                jac[8] = 0.5 * (ns[pyramidV3[j]]->m_z - ns[pyramidV0[j]]->m_z);

                NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                               jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                               jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

                if (jc < NekConstants::kNekZeroTol)
                {
                    return true;
                }
            }
        }
        else if (els[i]->GetShapeType() == LibUtilities::eTetrahedron)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            std::array<NekDouble, 9> jac;

            jac[0] = 0.5 * (ns[1]->m_x - ns[0]->m_x);
            jac[1] = 0.5 * (ns[1]->m_y - ns[0]->m_y);
            jac[2] = 0.5 * (ns[1]->m_z - ns[0]->m_z);
            jac[3] = 0.5 * (ns[2]->m_x - ns[0]->m_x);
            jac[4] = 0.5 * (ns[2]->m_y - ns[0]->m_y);
            jac[5] = 0.5 * (ns[2]->m_z - ns[0]->m_z);
            jac[6] = 0.5 * (ns[3]->m_x - ns[0]->m_x);
            jac[7] = 0.5 * (ns[3]->m_y - ns[0]->m_y);
            jac[8] = 0.5 * (ns[3]->m_z - ns[0]->m_z);

            NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                           jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                           jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

            if (jc < NekConstants::kNekZeroTol)
            {
                return true;
            }
        }
        else if (els[i]->GetShapeType() == LibUtilities::eHexahedron)
        {
            // Hexes are not checked at the moment!
            NekDouble jc = 1.0;
            if (jc < NekConstants::kNekZeroTol)
            {
                return true;
            }
        }
        else
        {
            m_log(FATAL) << "Only prisms, pyramids and tetrahedra supported."
                         << endl;
        }
    }

    return false;
}

void ProcessProjectCAD::Process()
{
    m_log(VERBOSE) << "Projecting CAD back onto linear mesh." << endl;

    m_log(WARNING) << "ProcessAssignCAD: Warning: This module is designed for "
                   << "use with Star-CCM+ meshes only; it also requires that "
                   << "the star mesh was created in a certain way." << endl;

    if (!m_config["order"].beenSet)
    {
        m_log(VERBOSE) << "Mesh order not set: will assume order 4" << endl;
    }

    // Projection Order
    int order         = m_config["order"].as<int>();
    m_mesh->m_nummode = order + 1;

    // Tolerances for vertex association
    NekDouble tolv1, tolv2;
    tolv1 = m_config["tolv1"].as<NekDouble>();
    tolv2 = m_config["tolv2"].as<NekDouble>();

    // Characteristic Length for CAD Reconstruction (auto calculation of tolv1
    // and tolv2)
    NekDouble cLength = 1.0;
    if (m_config["cLength"].beenSet)
    {
        cLength = m_config["cLength"].as<NekDouble>();
        tolv1 *= cLength;
        tolv2 *= cLength;
    }
    m_log(WARNING) << "Tolerances = " << tolv1 << " " << tolv2 << endl;
    // 1. Load CAD instance of the CAD model
    std::string filename = m_config["file"].as<string>();
    LoadCAD(filename);

    // 2. Create Bounding boxes of the CAD surfaces into a k-d tree
    NekDouble scale = cLength;
    bgi::rtree<boxI, bgi::quadratic<16>> rtree;
    bgi::rtree<boxI, bgi::quadratic<16>> rtreeCurve;
    bgi::rtree<boxI, bgi::quadratic<16>> rtreeNode;
    CreateBoundingBoxes(rtree, rtreeCurve, rtreeNode, tolv2, scale);

    m_log(VERBOSE) << "Bounding Boxes Surf/Curv/Vertex= " << rtree.size() << " "
                   << rtreeCurve.size() << " " << rtreeNode.size() << endl;
    // 3. Auxilaries ( can be moved to Module.cpp)
    // SurfNodes , surfNodeToEl, minConEdge
    Auxilaries();

    // 4.  Link Surface Vertices to CAD and Project them to the closest CAD
    LinkVertexToCAD(m_mesh, true, lockedNodes, tolv1, tolv2, rtree, rtreeCurve,
                    rtreeNode);

    // 5. clear the associations with CAD surfaces
    // necessary since the projection of the edges will change the surface uv
    // and the association will be wrong to some surfaces that were closed
    // beforehand
    for (auto vertex = surfNodes.begin(); vertex != surfNodes.end(); vertex++)
    {
        (*vertex)->ClearCADSurfs();
        (*vertex)->ClearCADCurves();
        // (*vertex)->ClearCADVert();
    }

    // 6. Update the secondary tolerances on already projected nodes and do the
    //  final Linking Vertex - CAD Surface / Curve
    tolv1 = tolv1 * 0.1, tolv2 = tolv2 * 0.1;
    LinkVertexToCAD(m_mesh, true, lockedNodes, tolv1, tolv2, rtree, rtreeCurve,
                    rtreeNode);

    // 7. Associate Edges to CAD
    LinkEdgeToCAD(surfEdges, tolv1);

    // 8. Associate Faces to CAD
    LinkFaceToCAD(tolv1);

    // Project the Edges to CAD that
    int orderCASE3 = 2;
    ProjectEdges(surfEdges, orderCASE3, rtree);

    if (m_config["ho"].beenSet)
    {
        LinkHOtoCAD(surfEdges, tolv1 * 10.0);
    }

    ////**** HOSurface ****////
    int m_surfopti         = m_config["surfopti"].as<bool>();
    ModuleSharedPtr module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "hosurface"), m_mesh);
    module->SetLogger(m_log);
    if (m_surfopti == 1)
    {
        //        module->RegisterConfig("opti","");
        module->RegisterConfig("third_party", "");

        try
        {
            module->SetDefaults();
            module->Process();
        }
        catch (runtime_error &e)
        {
            m_log(WARNING)
                << "High-order surface meshing has failed with message:"
                << endl;
            m_log(WARNING) << e.what() << endl;
            m_log(WARNING) << "The mesh will be written as normal but the "
                           << "incomplete surface will remain faceted" << endl;
            return;
        }
    }

    Diagnostics();
    if (m_config["extract"].beenSet)
    {
        m_log(VERBOSE) << "Extract CAD " << endl;
        ExtractCAD();
    }

    m_log(VERBOSE) << "HO-Surface CAD complete." << endl;
}

void ProcessProjectCAD::LoadCAD(std::string filename)
{
    m_log(VERBOSE) << "Start Loading the CAD file " << endl;
    ModuleSharedPtr module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh);
    module->RegisterConfig("filename", filename);
    module->SetDefaults();
    module->Process();
    m_log(VERBOSE) << "CAD loaded succesfully!" << endl;
}

void ProcessProjectCAD::Auxilaries()
{
    // find nodes on the surface
    // find unique nodes on the surface
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el      = m_mesh->m_element[2][i];
        vector<NodeSharedPtr> ns = el->GetVertexList();
        for (int j = 0; j < ns.size(); j++)
        {
            surfNodes.insert(ns[j]);
        }
    }

    // link surface nodes to their 3D element
    for (int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        if (m_mesh->m_element[3][i]->HasBoundaryLinks())
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[3][i]->GetVertexList();
            for (int j = 0; j < ns.size(); j++)
            {
                if (surfNodes.count(ns[j]) > 0)
                {
                    surfNodeToEl[ns[j]].push_back(m_mesh->m_element[3][i]);
                }
            }
        }
    }

    // Calculate the min edge length of the elements
    // Calculate MinConEdge
    CalculateMinEdgeLength();

    // make edges of surface mesh unique
    ClearElementLinks();
    // EdgeSet surfEdges;
    vector<ElementSharedPtr> &elmt = m_mesh->m_element[2];
    map<int, int> surfIdToLoc;
    for (int i = 0; i < elmt.size(); i++)
    {
        surfIdToLoc.insert(make_pair(elmt[i]->GetId(), i));
        for (int j = 0; j < elmt[i]->GetEdgeCount(); ++j)
        {
            pair<EdgeSet::iterator, bool> testIns;
            EdgeSharedPtr ed = elmt[i]->GetEdge(j);
            testIns          = surfEdges.insert(ed);

            if (testIns.second)
            {
                EdgeSharedPtr ed2 = *testIns.first;
                ed2->m_elLink.push_back(
                    pair<ElementSharedPtr, int>(elmt[i], j));
            }
            else
            {
                EdgeSharedPtr e2 = *(testIns.first);
                elmt[i]->SetEdge(j, e2);

                // Update edge to element map.
                e2->m_elLink.push_back(pair<ElementSharedPtr, int>(elmt[i], j));
            }
        }
    }
}

void ProcessProjectCAD::CreateBoundingBoxes(
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeSurf,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeCurve,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeNode, NekDouble tolv2,
    NekDouble scale)
{
    // CAD Surfs
    vector<boxI> boxes;
    for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
    {
        // m_log(VERBOSE).Progress(i, m_mesh->m_cad->GetNumSurf(),
        //                         "building surface bboxes", i - 1);
        auto bx = m_mesh->m_cad->GetSurf(i)->BoundingBox(scale);
        boxes.push_back(make_pair(
            box(point(bx[0], bx[1], bx[2]), point(bx[3], bx[4], bx[5])), i));

        m_log(VERBOSE) << " Bounding box = " << i << endl;
        m_log(VERBOSE) << bx[0] << " " << bx[1] << " " << bx[2] << endl;
        m_log(VERBOSE) << bx[3] << " " << bx[4] << " " << bx[5] << endl;
    }

    m_log(VERBOSE).Newline();
    m_log(VERBOSE) << "Building Surf admin data structures." << endl;
    rtreeSurf.insert(boxes.begin(), boxes.end());
    m_log(VERBOSE) << "Bounding Box ." << endl;

    // CAD Curves
    boxes.clear();
    for (int i = 1; i <= m_mesh->m_cad->GetNumCurve(); i++)
    {
        auto bx = m_mesh->m_cad->GetCurve(i)->BoundingBox(scale);

        boxes.push_back(make_pair(
            box(point(bx[0], bx[1], bx[2]), point(bx[3], bx[4], bx[5])), i));
    }
    rtreeCurve.insert(boxes.begin(), boxes.end());

    // CAD Vertices
    boxes.clear();
    NekDouble tol = tolv2; // 1e-8 * scale
    for (int i = 1; i <= m_mesh->m_cad->GetNumVerts(); i++)
    {
        auto vert                     = m_mesh->m_cad->GetVert(i);
        std::array<NekDouble, 3> locT = vert->GetLoc();
        boxes.push_back(
            make_pair(box(point(locT[0] - tol, locT[1] - tol, locT[2] - tol),
                          point(locT[0] + tol, locT[1] + tol, locT[2] + tol)),
                      i));
    }
    rtreeNode.insert(boxes.begin(), boxes.end());
}

void ProcessProjectCAD::CalculateMinEdgeLength()
{
    // link the surface node to a value for the shortest connecting edge to it
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el      = m_mesh->m_element[2][i];
        vector<NodeSharedPtr> ns = el->GetVertexList();
        NekDouble l1             = ns[0]->Distance(ns[1]);
        NekDouble l2             = ns[1]->Distance(ns[2]);
        NekDouble l3             = ns[2]->Distance(ns[0]);

        if (minConEdge.count(ns[0]))
        {
            NekDouble l       = minConEdge[ns[0]];
            minConEdge[ns[0]] = min(l, min(l1, l3));
        }
        else
        {
            minConEdge.insert(make_pair(ns[0], min(l1, l3)));
        }
        if (minConEdge.count(ns[1]))
        {
            NekDouble l       = minConEdge[ns[1]];
            minConEdge[ns[1]] = min(l, min(l1, l1));
        }
        else
        {
            minConEdge.insert(make_pair(ns[1], min(l1, l1)));
        }
        if (minConEdge.count(ns[2]))
        {
            NekDouble l       = minConEdge[ns[2]];
            minConEdge[ns[2]] = min(l, min(l2, l3));
        }
        else
        {
            minConEdge.insert(make_pair(ns[2], min(l2, l3)));
        }
    }
}

void ProcessProjectCAD::LinkVertexToCAD(
    NekMesh::MeshSharedPtr &m_mesh, bool projectV, NodeSet &lockedNodes,
    NekDouble tolv1, NekDouble tolv2,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtree,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeCurve,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeNode)
{
    map<int, vector<int>> finds;

    m_log(VERBOSE) << "Searching tree." << endl;

    NekDouble maxNodeCor = 0;

    // find nodes surface and parametric location
    int ct = 0;
    for (auto i = surfNodes.begin(); i != surfNodes.end(); i++, ct++)
    {
        m_log(VERBOSE).Progress(ct, surfNodes.size(), "projecting verts",
                                ct - 1);

        point q((*i)->m_x, (*i)->m_y, (*i)->m_z);
        vector<boxI> result;
        rtree.query(bgi::intersects(q), back_inserter(result));

        if (result.size() == 0)
        {
            // Vertex is too far from any surface bounding boxes
            m_log(WARNING)
                << "Vertex  " << (*i)
                << " is not in any boundin boxes. 1. "
                   "Make sure you use the correct STEP file. 2. The problem is "
                   "likely in the linear mesh -> try to refine this region.  "
                << endl;
            continue;
        }

        // Vertex Tolerances to the CAD Surface - 0.5 * min edge length + [min
        // tol max tol]
        NekDouble tol = minConEdge[*i] * 0.5;
        tol           = min(tol, tolv1);
        tol           = max(tol, tolv2);

        vector<int> distId;
        vector<NekDouble> distList;
        // sort the surfaces by distance to the node
        for (int j = 0; j < result.size(); j++)
        {
            NekDouble dist;
            m_mesh->m_cad->GetSurf(result[j].second)
                ->locuv((*i)->GetLoc(), dist);
            distList.push_back(dist);
            distId.push_back(result[j].second);
        }

        bool repeat = true;
        while (repeat)
        {
            repeat = false;
            for (int j = 0; j < distId.size() - 1; j++)
            {
                if (distList[j + 1] < distList[j])
                {
                    repeat = true;
                    swap(distList[j + 1], distList[j]);
                    swap(distId[j + 1], distId[j]);
                }
            }
        }

        int pos = 0;
        for (int j = 0; j < distId.size(); j++)
        {
            if (distList[j] < tol)
            {
                pos++;
            }
        }

        distId.resize(pos);

        finds[pos].push_back(0);
        // if the node is not close to any surface lock it (no CAD given + no
        // projection for its edges)
        if (pos == 0)
        {
            lockedNodes.insert(*i);
            m_log(WARNING) << "surface minDist = " << distList[0] << " unknown "
                           << "(tolerance: " << tol << ")   xyz = " << (*i)
                           << endl;
        }
        else
        {
            NekDouble shift;
            bool st = false;
            for (int j = 0; j < distId.size(); j++)
            {
                if (distList[j] > tol)
                {
                    continue;
                }
                if (m_mesh->m_cad->GetSurf(distId[j])->IsPlanar())
                {
                    // continue;
                }

                shift                     = distList[j];
                NekDouble dist            = 0;
                CADSurfSharedPtr s        = m_mesh->m_cad->GetSurf(distId[j]);
                auto l                    = (*i)->GetLoc();
                [[maybe_unused]] auto uvt = s->locuv(l, dist);

                if (true)
                {
                    NekDouble tmpX = (*i)->m_x;
                    NekDouble tmpY = (*i)->m_y;
                    NekDouble tmpZ = (*i)->m_z;

                    (*i)->m_x = s->P(uvt)[0];
                    (*i)->m_y = s->P(uvt)[1];
                    (*i)->m_z = s->P(uvt)[2];

                    if (ProcessProjectCAD::IsNotValid(surfNodeToEl[*i]))
                    {
                        (*i)->m_x = tmpX;
                        (*i)->m_y = tmpY;
                        (*i)->m_z = tmpZ;

                        m_log(VERBOSE) << "Element not valid after vertex ";
                        m_log(VERBOSE)
                            << "projection reset it and lock the vertex  "
                            << tmpX << " " << tmpY << " " << tmpZ << endl;
                        break;
                    }
                }

                st = true;
                break;
            }

            if (!st)
            {
                lockedNodes.insert(*i);
                continue;
            }

            for (int j = 0; j < distId.size(); j++)
            {
                if (distList[j] > tol)
                {
                    continue;
                }
                if (m_mesh->m_cad->GetSurf(distId[j])->IsPlanar())
                {
                    // continue;
                }

                CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(distId[j]);
                NekDouble dist     = 0;
                auto loc           = (*i)->GetLoc();
                auto uv            = s->locuv(loc, dist);
                (*i)->SetCADSurf(s, uv);
            }
            maxNodeCor = max(maxNodeCor, shift);
        }
    }

    // for the vertices with multiple CAD Surfaces, check CADCurves with the
    // bounding boxes Intersect with CAD Curves rtree
    bool CADCurve = true;
    if (CADCurve)
    {
        for (auto vertex : surfNodes)
        {
            // CAD Curve
            if (vertex->GetCADSurfs().size() > 1)
            {
                point point((vertex)->m_x, (vertex)->m_y, (vertex)->m_z);
                vector<boxI> result;
                rtreeCurve.query(bgi::intersects(point), back_inserter(result));

                if (result.size() == 1)
                {
                    // Single CAD Curve
                    CADCurveSharedPtr CADCurve_t =
                        m_mesh->m_cad->GetCurve(result[0].second);
                    NekDouble t0, dist0;
                    NekDouble tmin, tmax;
                    CADCurve_t->GetBounds(tmin, tmax);

                    dist0 = CADCurve_t->loct(vertex->GetLoc(), t0, tmin, tmax);
                    if (dist0 < tolv1)
                    {
                        vertex->SetCADCurve(CADCurve_t, t0);
                        if (projectV)
                        {
                            NekDouble tmpX = vertex->m_x;
                            NekDouble tmpY = vertex->m_y;
                            NekDouble tmpZ = vertex->m_z;

                            vertex->m_x =
                                CADCurve_t->P(t0)[0]; // This gets the node far
                                                      // from CADSurf???
                            vertex->m_y = CADCurve_t->P(t0)[1];
                            vertex->m_z = CADCurve_t->P(t0)[2];

                            if (ProcessProjectCAD::IsNotValid(
                                    surfNodeToEl[vertex]))
                            {
                                vertex->m_x = tmpX;
                                vertex->m_y = tmpY;
                                vertex->m_z = tmpZ;

                                m_log(VERBOSE)
                                    << "Element not valid after vertex ";
                                m_log(VERBOSE) << "projection reset it and "
                                                  "lock the vertex  "
                                               << vertex << endl;
                            }
                        }
                    }
                    else
                    {
                        m_log(TRACE)
                            << "Vertex " << vertex
                            << " not close enough to the CAD Curve " << endl;
                    }
                }
                else if (result.size() == 0)
                {
                    m_log(WARNING) << "No CAD Curve found in the bounding box"
                                   << vertex << endl;
                    continue;
                }
                else
                {
                    m_log(TRACE) << "Multiple CAD Curves found for the "
                                 << "vertex " << vertex << endl;
                    for (auto res : result)
                    {
                        CADCurveSharedPtr CADCurve_t =
                            m_mesh->m_cad->GetCurve(res.second);
                        NekDouble t0, dist0;
                        NekDouble tmin, tmax;
                        CADCurve_t->GetBounds(tmin, tmax);
                        dist0 =
                            CADCurve_t->loct(vertex->GetLoc(), t0, tmin, tmax);
                        if (dist0 < tolv1)
                        {
                            vertex->SetCADCurve(CADCurve_t, t0);
                            if (projectV)
                            {
                                NekDouble tmpX = vertex->m_x;
                                NekDouble tmpY = vertex->m_y;
                                NekDouble tmpZ = vertex->m_z;

                                vertex->m_x =
                                    CADCurve_t->P(t0)[0]; // This gets the node
                                                          // far from CADSurf???
                                vertex->m_y = CADCurve_t->P(t0)[1];
                                vertex->m_z = CADCurve_t->P(t0)[2];

                                if (ProcessProjectCAD::IsNotValid(
                                        surfNodeToEl[vertex]))
                                {
                                    vertex->m_x = tmpX;
                                    vertex->m_y = tmpY;
                                    vertex->m_z = tmpZ;

                                    m_log(VERBOSE)
                                        << "Element not valid after vertex ";
                                    m_log(VERBOSE) << "projection reset it and "
                                                      "lock the vertex  "
                                                   << vertex << endl;
                                }
                            }
                        }
                        else
                        {
                            m_log(TRACE) << "Vertex " << vertex
                                         << " not close enough to the CAD Curve"
                                         << dist0 << " tol = " << tolv2 << endl;
                        }
                    }
                }
            }

            // CAD Vertex (for the moment no benefits )
            bool CADVertexAssociation = true;
            if (CADVertexAssociation)
            {
                if (vertex->GetCADSurfs().size() > 1)
                {
                    point point((vertex)->m_x, (vertex)->m_y, (vertex)->m_z);
                    vector<boxI> result;
                    rtreeNode.query(bgi::intersects(point),
                                    back_inserter(result));

                    if (result.size() == 1)
                    {
                        // Single CAD Vertex ( we do not project on it for now!)
                        vertex->SetCADVertex(m_mesh->m_cad->GetVert(
                            result[0]
                                .second)); // No Adj Curves Assigned to the V !
                    }
                    else if (result.size() > 1)
                    {
                        // Multiple CAD Vertices
                        m_log(WARNING) << "Multiple CAD Vertices found for the "
                                       << "vertex " << vertex
                                       << "  size= " << result.size() << endl;
                    }
                }
            }
        }
    }

    m_log(VERBOSE) << "  - max surface Vertex correction " << maxNodeCor
                   << endl;
    m_log(VERBOSE) << "  - lockedNodes N= " << lockedNodes.size() << endl;
}

void ProcessProjectCAD::LinkEdgeToCAD(EdgeSet &surfEdges, NekDouble tolv1)
{
    // Every Edge needs to have only 1 CAD Object CADCurve or CADSurf
    // This is not necessary to be perfect as it will be filled by the Associate
    // Faces However it can be used as a verification for the FACE association
    // in the future It is beneficial to associate the edge to CAD Curve due to
    // optimization and projection sliding on the CAD Curve is more rorbust than
    // to the CADSurf .

    for (auto edge : surfEdges)
    {
        NodeSharedPtr v1 = edge->m_n1;
        NodeSharedPtr v2 = edge->m_n2;

        if (lockedNodes.count(v1) || lockedNodes.count(v2))
        {
            continue;
        }

        // 1. Get CAD Curve based on the vertex CADCurves - should be 99% of
        // CADCurve edges
        if (v1->GetCADCurves().size() && v2->GetCADCurves().size())
        {
            if (v1->GetCADCurves()[0] == v2->GetCADCurves()[0])
            {
                edge->m_parentCAD = v1->GetCADCurves()[0];
                continue;
            }
            else
            {
                vector<int> cmn =
                    IntersectCADCurve(v1->GetCADCurves(), v2->GetCADCurves());
                if (cmn.size() == 1)
                {
                    edge->m_parentCAD = m_mesh->m_cad->GetCurve(cmn[0]);
                    continue;
                }
                else if (cmn.size() > 1)
                {
                    // This is often the case when you have two CAD Curves that
                    // are the same, but  topologically different and OCE CAD
                    // Sewing (sew_tolerance) has not merged them
                    m_log(WARNING)
                        << "Edge with different CAD Curves cmn.size=  "
                        << cmn.size() << " v1 = " << edge->m_n1
                        << " v2 = " << edge->m_n2 << endl;
                }
            }
        }

        // 2. Try to associate the edge to CADCurves based on the vertex CADSurf
        vector<int> cmn =
            IntersectCADSurf(v1->GetCADSurfs(), v2->GetCADSurfs());

        if (cmn.size() == 0)
        {
            // no CAD surface found for the edge (CASE3)
            m_log(TRACE)
                << "Case 3 edge association (NO-CADSurf or Curve) v1 = " << v1
                << "  = v2 " << v2 << endl;
            continue;
        }

        if (cmn.size() == 1)
        {
            // Clearly Edge is on a single CAD surface (internal)
            edge->m_parentCAD = m_mesh->m_cad->GetSurf(cmn[0]);
        }
        else if (cmn.size() == 2)
        {
            // N=2 CAD Surfaces could be CAD-curve or CAD-surface (CASE2)
            // Try to find the correct edge topologically
            vector<CADSurfSharedPtr> v1_CAD = v1->GetCADSurfs();
            vector<CADSurfSharedPtr> v2_CAD = v2->GetCADSurfs();

            // 1.Create vi1 , vi2
            vector<int> CADCurves_uv1;
            vector<int> CADCurves_uv2;

            CADSurfSharedPtr EdgeSurf1 = m_mesh->m_cad->GetSurf(cmn[0]);
            CADSurfSharedPtr EdgeSurf2 = m_mesh->m_cad->GetSurf(cmn[1]);

            // Checking overlapping CADCurves between the surfaces
            for (auto EdgeLoop : EdgeSurf1->GetEdges())
            {
                for (auto CADCurve : EdgeLoop->edges)
                {
                    CADCurves_uv1.push_back(CADCurve->GetId());
                }
            }
            for (auto EdgeLoop : EdgeSurf2->GetEdges())
            {
                for (auto CADCurve : EdgeLoop->edges)
                {
                    CADCurves_uv2.push_back(CADCurve->GetId());
                }
            }

            sort(CADCurves_uv1.begin(), CADCurves_uv1.end());
            sort(CADCurves_uv2.begin(), CADCurves_uv2.end());

            vector<int> commonCADCurves;
            set_intersection(CADCurves_uv1.begin(), CADCurves_uv1.end(),
                             CADCurves_uv2.begin(), CADCurves_uv2.end(),
                             back_inserter(commonCADCurves));

            NekDouble tolDist = tolv1; // POSSIBLE PROBLEM !!!!
            if (commonCADCurves.size() == 1)
            {
                CADCurveSharedPtr CADCurve_t =
                    m_mesh->m_cad->GetCurve(commonCADCurves[0]);

                NekDouble t0, t1, dist0, dist1;
                NekDouble tmin, tmax;
                CADCurve_t->GetBounds(tmin, tmax);

                dist0 = CADCurve_t->loct(edge->m_n1->GetLoc(), t0, tmin, tmax);
                dist1 = CADCurve_t->loct(edge->m_n2->GetLoc(), t1, tmin, tmax);

                if ((dist0 < tolDist) && (dist1 < tolDist))
                {
                    edge->m_n1->SetCADCurve(CADCurve_t, t0);
                    edge->m_n2->SetCADCurve(CADCurve_t, t1);
                    edge->m_parentCAD = CADCurve_t;
                }
                else
                {
                    // Associate to one of the CADSurf.
                    m_log(VERBOSE)
                        << " dist > distol (cmnCADCurve.size=1) = " << dist0
                        << " " << dist1 << endl;
                }
            }
            else if (commonCADCurves.size() >= 2)
            {
                // common when the CAD is not perfect (2 CADcurves that overlap
                // are two different topological objects ) in this case just
                // take the curve and assign the closest one within tolv2*0.1
                // tolerance this is a stricter due to the

                // Another test case is a CADSurf like NACA, where it fills the
                // DO WE NEED THIS ?
                edge->m_parentCAD = EdgeSurf1;
                m_log(VERBOSE) << "edge CASE commonCADCurves.size()==2 for "
                                  "comn.size()  = 2     xyz= "
                               << edge->m_n1 << endl;
            }
        }
        else
        {
            // If more than 2 common CAD Surfaces are present, we do not
            // associate CADCurve because it is too risky Closest CAD Surf
            m_log(VERBOSE) << "too many common surfaces for Edge association  "
                              "(cmn>2) will use the element to associate. "
                           << endl;
        }

        // if no CAD Curves are associated, the CAD Surf will be associated
        // through the Face Association.
    }
}

vector<int> ProcessProjectCAD::IntersectCADSurf(
    vector<CADSurfSharedPtr> v1_CADs, vector<CADSurfSharedPtr> v2_CADs)
{
    vector<int> cmn;
    if (v1_CADs.size() == 0 || v2_CADs.size() == 0)
    {
        return cmn;
    }

    vector<CADSurfSharedPtr> v1 = v1_CADs;
    vector<CADSurfSharedPtr> v2 = v2_CADs;

    vector<int> vi1, vi2;
    for (size_t j = 0; j < v1.size(); ++j)
    {
        vi1.push_back(v1[j]->GetId());
    }
    for (size_t j = 0; j < v2.size(); ++j)
    {
        vi2.push_back(v2[j]->GetId());
    }

    sort(vi1.begin(), vi1.end());
    sort(vi2.begin(), vi2.end());

    set_intersection(vi1.begin(), vi1.end(), vi2.begin(), vi2.end(),
                     back_inserter(cmn));

    return cmn;
}

vector<int> ProcessProjectCAD::IntersectCADCurve(
    vector<CADCurveSharedPtr> v1_CADs, vector<CADCurveSharedPtr> v2_CADs)
{
    vector<CADCurveSharedPtr> v1 = v1_CADs;
    vector<CADCurveSharedPtr> v2 = v2_CADs;

    vector<int> vi1, vi2;
    for (size_t j = 0; j < v1.size(); ++j)
    {
        vi1.push_back(v1[j]->GetId());
    }
    for (size_t j = 0; j < v2.size(); ++j)
    {
        vi2.push_back(v2[j]->GetId());
    }

    sort(vi1.begin(), vi1.end());
    sort(vi2.begin(), vi2.end());

    vector<int> cmn;
    set_intersection(vi1.begin(), vi1.end(), vi2.begin(), vi2.end(),
                     back_inserter(cmn));

    return cmn;
}

void ProcessProjectCAD::ProjectEdges(
    EdgeSet &surfEdges, int order, bgi::rtree<boxI, bgi::quadratic<16>> &rtree)
{
    m_log(VERBOSE) << " Projecting Edges to CAD (CASE3)" << endl;
    // Project the Edges to CAD

    LibUtilities::PointsKey ekey(order + 1,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    // make surface edges high-order
    int cnt  = 0;
    int cnt1 = 0;
    for (auto i = surfEdges.begin(); i != surfEdges.end(); i++)
    {
        // IF the edge is already associated with a CAD surface, HOSurf will do
        // the curving job
        if ((*i)->m_parentCAD)
        {
            continue;
        }
        if (lockedNodes.count((*i)->m_n1) || lockedNodes.count((*i)->m_n2))
        {
            continue;
        }
        cnt++;

        vector<CADSurfSharedPtr> v1 = (*i)->m_n1->GetCADSurfs();
        vector<CADSurfSharedPtr> v2 = (*i)->m_n2->GetCADSurfs();

        vector<int> vi1, vi2, vi1vi2;
        for (size_t j = 0; j < v1.size(); ++j)
        {
            vi1.push_back(v1[j]->GetId());
        }
        for (size_t j = 0; j < v2.size(); ++j)
        {
            vi2.push_back(v2[j]->GetId());
        }

        sort(vi1.begin(), vi1.end());
        sort(vi2.begin(), vi2.end());

        vector<int> cmn;
        set_intersection(vi1.begin(), vi1.end(), vi2.begin(), vi2.end(),
                         back_inserter(cmn));

        (*i)->m_curveType = LibUtilities::eGaussLobattoLegendre;

        if (cmn.size() == 1 || cmn.size() == 2)
        {
            for (int j = 0; j < cmn.size(); j++)
            {
                if (m_mesh->m_cad->GetSurf(cmn[j])->IsPlanar())
                {
                    // if its planar dont care
                    continue;
                }

                auto uvb = (*i)->m_n1->GetCADSurfInfo(cmn[j]);
                auto uve = (*i)->m_n2->GetCADSurfInfo(cmn[j]);

                // can compare the loction of the projection to the
                // corresponding position of the straight sided edge
                // if the two differ by more than the length of the edge
                // something has gone wrong
                NekDouble len = (*i)->m_n1->Distance((*i)->m_n2);

                for (int k = 1; k < order + 1 - 1; k++)
                {
                    std::array<NekDouble, 2> uv = {
                        uvb[0] * (1.0 - gll[k]) / 2.0 +
                            uve[0] * (1.0 + gll[k]) / 2.0,
                        uvb[1] * (1.0 - gll[k]) / 2.0 +
                            uve[1] * (1.0 + gll[k]) / 2.0};
                    auto loc = m_mesh->m_cad->GetSurf(cmn[j])->P(uv);

                    std::array<NekDouble, 3> locT;
                    locT[0] = (*i)->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
                              (*i)->m_n2->m_x * (1.0 + gll[k]) / 2.0;
                    locT[1] = (*i)->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
                              (*i)->m_n2->m_y * (1.0 + gll[k]) / 2.0;
                    locT[2] = (*i)->m_n1->m_z * (1.0 - gll[k]) / 2.0 +
                              (*i)->m_n2->m_z * (1.0 + gll[k]) / 2.0;

                    NekDouble d = sqrt((locT[0] - loc[0]) * (locT[0] - loc[0]) +
                                       (locT[1] - loc[1]) * (locT[1] - loc[1]) +
                                       (locT[2] - loc[2]) * (locT[2] - loc[2]));

                    if (d > len)
                    {
                        (*i)->m_edgeNodes.clear();
                        break;
                    }

                    NodeSharedPtr nn = std::shared_ptr<Node>(
                        new Node(0, loc[0], loc[1], loc[2]));

                    (*i)->m_edgeNodes.push_back(nn);
                }

                if ((*i)->m_edgeNodes.size() != 0)
                {
                    // it suceeded on this surface so skip the other possibility
                    break;
                }
            }
        }
        else if (cmn.size() == 0)
        {
            // projection, if the projection requires more than two surfaces
            // including the edge nodes, then,  in theory projection shouldnt be
            // used
            vi1vi2.insert(vi1vi2.end(), vi1.begin(), vi1.end());
            vi1vi2.insert(vi1vi2.end(), vi2.begin(), vi2.end());

            set<int> sused;
            for (int k = 1; k < order + 1 - 1; k++)
            {
                std::array<NekDouble, 3> locT = {
                    (*i)->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
                        (*i)->m_n2->m_x * (1.0 + gll[k]) / 2.0,
                    (*i)->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
                        (*i)->m_n2->m_y * (1.0 + gll[k]) / 2.0,
                    (*i)->m_n1->m_z * (1.0 - gll[k]) / 2.0 +
                        (*i)->m_n2->m_z * (1.0 + gll[k]) / 2.0};

                int s;
                if (!FindAndProject(rtree, locT, s))
                {
                    (*i)->m_edgeNodes.clear();
                    m_log(VERBOSE) << "failed to find CAD" << endl;
                    break;
                }

                sused.insert(s);

                if (sused.size() > 2)
                {
                    m_log(WARNING) << "found too many CAD " << endl;
                    (*i)->m_edgeNodes.clear();
                    break;
                }

                NodeSharedPtr nn = std::shared_ptr<Node>(
                    new Node(0, locT[0], locT[1], locT[2]));

                (*i)->m_edgeNodes.push_back(nn);
            }
            cnt1++;
        }
    }

    m_log(VERBOSE) << "Edges No CAD N= " << cnt << endl;
    m_log(VERBOSE) << "Edges No CAD Projected N= " << cnt1 << endl;
}

void ProcessProjectCAD::Diagnostics()
{
    m_log(VERBOSE) << endl
                   << " ---------Diagnostics---------     " << endl
                   << endl;

    if (true)
    {
        long int counterElementsNoCAD = 0;
        long int counterEdgesNoCAD    = 0;
        long int counterVerticesNoCAD = 0;

        m_log(VERBOSE) << "         CAD Stats       " << endl;
        m_log(VERBOSE) << "CAD Vertices =           "
                       << m_mesh->m_cad->GetNumVerts() << endl;
        m_log(VERBOSE) << "CAD Curves =             "
                       << m_mesh->m_cad->GetNumCurve() << endl;
        m_log(VERBOSE) << "CAD Surfaces =           "
                       << m_mesh->m_cad->GetNumSurf() << endl;

        m_log(VERBOSE) << "        Mesh Stats       " << endl;
        m_log(VERBOSE) << "Mesh Surface Vertices =  " << surfNodes.size()
                       << endl;
        m_log(VERBOSE) << "Mesh Surface Edges =     " << surfEdges.size()
                       << endl;
        m_log(VERBOSE) << "Mesh Surface Face =      "
                       << m_mesh->m_element[2].size() << endl;

        // Elements without CAD
        EdgeSet NoCADEdges;
        for (auto element : m_mesh->m_element[2])
        {
            if (element->m_parentCAD == nullptr)
            {
                counterElementsNoCAD++;
                m_log(TRACE) << "Element No CAD Vertex1 xyz= "
                             << element->GetVertex(0)->m_x << " "
                             << element->GetVertex(0)->m_y << " "
                             << element->GetVertex(0)->m_z << endl;
            }
        }

        // Edges without CAD
        int cntEdgeCurve = 0, cntEdgeSurf = 0;
        for (auto edge : surfEdges)
        {
            if (edge->m_parentCAD == nullptr)
            {
                counterEdgesNoCAD++;
            }
            else if (edge->m_parentCAD->GetType() == CADType::eCurve)
            {
                cntEdgeCurve++;
            }
            else
            {
                cntEdgeSurf++;
            }
        }

        // Vertices CAD
        int cntCAD2 = 0, cntCAD1 = 0, cntCAD3orMore = 0;
        int cntCADVertices = 0, cntCADCurves = 0;
        for (auto vertex : surfNodes)
        {
            if (vertex->GetCADSurfs().size() == 0 &&
                vertex->GetCADCurves().size() == 0)
            {
                counterVerticesNoCAD++;
            }

            if (vertex->GetCADSurfs().size() == 1)
            {
                cntCAD1++;
            }

            if (vertex->GetCADSurfs().size() == 2)
            {
                cntCAD2++;
            }

            if (vertex->GetCADSurfs().size() > 2)
            {
                cntCAD3orMore++;
            }

            if (vertex->GetCADVertex() != nullptr)
            {
                cntCADVertices++;
            }
            if (vertex->GetCADCurves().size() > 0)
            {
                cntCADCurves++;
            }
        }

        // Surface Edges without CAD

        m_log(WARNING) << "Vertices No CAD (Includes PlanarSurf) N= "
                       << counterVerticesNoCAD << endl;
        // m_log(WARNING) << "Planar Vertices =                        " <<
        // planarCnt << endl;
        m_log(WARNING) << "Edges without CADObject               N= "
                       << counterEdgesNoCAD << endl;
        m_log(WARNING) << "Faces without CADObject               N= "
                       << counterElementsNoCAD << endl;

        m_log(VERBOSE) << "Vertices 1 CADSurf                    N= " << cntCAD1
                       << endl;
        m_log(VERBOSE) << "Vertices 2 CADSurf                    N= " << cntCAD2
                       << endl;
        m_log(VERBOSE) << "Vertices CAD3 or more Surf            N= "
                       << cntCAD3orMore << endl;
        m_log(VERBOSE) << "Vertices- CADVertex (also have CADSu) N= "
                       << cntCADVertices << endl;
        m_log(VERBOSE) << "Vertices- CADCurve                    N= "
                       << cntCADCurves << endl;
        m_log(VERBOSE) << "Edges CADCurve                        N= "
                       << cntEdgeCurve << endl;
        m_log(VERBOSE) << "Edges CADSurf                         N= "
                       << cntEdgeSurf << endl;
    }
}

void ProcessProjectCAD::LinkFaceToCAD(NekDouble tolv1)
{
    for (auto element : m_mesh->m_element[2])
    {
        vector<NodeSharedPtr> vertices = element->GetVertexList();
        vector<EdgeSharedPtr> edges    = element->GetEdgeList();

        // CASE1 - all Edges internal to same CADSurf
        bool internal = true;
        for (int i = 1; i < edges.size(); i++)
        {
            if (edges[i]->m_parentCAD != edges[i - 1]->m_parentCAD)
            {
                internal = false;
                break;
            }
        }

        if (internal && edges[0]->m_parentCAD &&
            edges[0]->m_parentCAD->GetType() == 2)
        {
            element->m_parentCAD = edges[0]->m_parentCAD;
            continue;
        }

        // CASE2 and CASE3 - 2 or more CADSurfs or None (Use vertice CAD)
        vector<vector<int>> cmn;
        for (int i = 1; i < vertices.size(); i++)
        {
            vector<int> cmn_i = IntersectCADSurf(
                vertices[i]->GetCADSurfs(), vertices[i - 1]->GetCADSurfs());
            if (cmn_i.size() > 0)
            {
                cmn.push_back(cmn_i);
            }
            else
            {
                cmn.clear();
                break;
            }
        }

        if (cmn.size() == 0)
        {
            // CASE 3
            continue;
        }

        std::vector<int> commonCAD = cmn[0];
        for (auto cmn_i : cmn)
        {
            std::sort(cmn_i.begin(), cmn_i.end());
            std::vector<int> temp;
            std::set_intersection(commonCAD.begin(), commonCAD.end(),
                                  cmn_i.begin(), cmn_i.end(),
                                  std::back_inserter(temp));
            commonCAD = temp;
        }

        // commonCAD CADSurf found in the individual vertices
        if (commonCAD.size() == 1)
        {
            // Internal element based on the
            element->m_parentCAD = m_mesh->m_cad->GetSurf(commonCAD[0]);
            for (auto edge : edges)
            {
                if (edge->m_parentCAD == nullptr)
                {
                    edge->m_parentCAD = element->m_parentCAD;
                }
            }
        }
        else
        {
            // CASE 2 - 2 or more CADSurfs (max 1-2% of the faces)
            // This could be a trailing edge surface for example
            // Sliver Surface on IFW, etc
            // or very thin surface, where all vertices share >1 CAD Surf
            // Solution :: use edge nodes and face nodes to check which one is
            // closest FaceNode Effect x2, Edge

            // Get center of linear triag/quad
            std::array<NekDouble, 3> center = {0.0, 0.0, 0.0};
            for (auto vert : vertices)
            {
                center[0] += vert->m_x / vertices.size();
                center[1] += vert->m_y / vertices.size();
                center[2] += vert->m_z / vertices.size();
            }

            // Check mid of the distance of the centroids of edges and mid face
            // to every CAD
            NekDouble minDist = tolv1 * 10.0;
            int minID         = -1;
            for (int id : commonCAD)
            {
                std::array<NekDouble, 4> lim;
                CADSurfSharedPtr surf = m_mesh->m_cad->GetSurf(id);
                surf->GetBounds(lim[0], lim[1], lim[2], lim[3]);
                std::array<NekDouble, 3> loc = {0.0, 0.0, 0.0};
                NekDouble distoveral         = 0.0;
                for (auto edge : edges)
                {
                    loc[0] =
                        (edge->m_n1->GetLoc()[0] + edge->m_n2->GetLoc()[0]) *
                        0.5;
                    loc[1] =
                        (edge->m_n1->GetLoc()[1] + edge->m_n2->GetLoc()[1]) *
                        0.5;
                    loc[2] =
                        (edge->m_n1->GetLoc()[2] + edge->m_n2->GetLoc()[2]) *
                        0.5;

                    NekDouble dist = 1e7;

                    surf->locuv(loc, dist, lim[0], lim[1], lim[2], lim[3]);
                    distoveral += dist;
                }

                NekDouble dist = 1e7;
                surf->locuv(center, dist);
                distoveral += dist;

                if (distoveral < minDist)
                {
                    minID   = id;
                    minDist = distoveral;
                }
            }

            if (minID != -1)
            {
                element->m_parentCAD = m_mesh->m_cad->GetSurf(minID);
                for (auto edge : edges)
                {
                    if (!edge->m_parentCAD)
                    {
                        edge->m_parentCAD = m_mesh->m_cad->GetSurf(minID);
                    }
                }
            }

            // Choose the smallest one if within 1/10 min edge
        }
    }
}

void ProcessProjectCAD::LinkHOtoCAD(EdgeSet &surfEdges, NekDouble tolv1)
{
    for (auto edge : surfEdges)
    {
        if (edge->m_parentCAD && edge->m_edgeNodes.size() > 0)
        {
            // loop over the edges and assign CADCurve
            if (edge->m_parentCAD->GetType() == 1)
            {
                // CAD Curve
                CADCurveSharedPtr curve =
                    m_mesh->m_cad->GetCurve(edge->m_parentCAD->GetId());
                std::array<NekDouble, 2> lim;
                curve->GetBounds(lim[0], lim[1]);
                for (auto node : edge->m_edgeNodes)
                {
                    NekDouble dist               = 1e6;
                    std::array<NekDouble, 3> loc = node->GetLoc();
                    NekDouble t = curve->loct(loc, dist, lim[0], lim[1]);
                    if (dist < tolv1)
                    {
                        // Just give the node a parametric location
                        // DO NOT Project the node for the moment !
                        node->SetCADCurve(curve, t);
                    }
                }
            }
            else if (edge->m_parentCAD->GetType() == 2)
            {
                // CAD Surf
                CADSurfSharedPtr surf =
                    m_mesh->m_cad->GetSurf(edge->m_parentCAD->GetId());
                std::array<NekDouble, 4> lim;
                surf->GetBounds(lim[0], lim[1], lim[2], lim[3]);
                for (auto node : edge->m_edgeNodes)
                {
                    NekDouble dist               = 1e6;
                    std::array<NekDouble, 3> loc = node->GetLoc();
                    std::array<NekDouble, 2> uv =
                        surf->locuv(loc, dist, lim[0], lim[1], lim[2], lim[3]);
                    if (dist < tolv1)
                    {
                        // Just give the node a parametric location
                        // DO NOT Project the node for the moment !
                        node->SetCADSurf(surf, uv);
                    }
                }
            }
        }
    }

    for (auto face : m_mesh->m_element[2])
    {
        if (face->m_parentCAD)
        {
            CADSurfSharedPtr surf =
                m_mesh->m_cad->GetSurf(face->m_parentCAD->GetId());
            std::array<NekDouble, 4> lim;
            surf->GetBounds(lim[0], lim[1], lim[2], lim[3]);

            vector<NodeSharedPtr> nodelist;
            face->GetCurvedNodes(nodelist);
            for (auto node : nodelist)
            {
                NekDouble dist               = 1e6;
                std::array<NekDouble, 3> loc = node->GetLoc();
                std::array<NekDouble, 2> uv =
                    surf->locuv(loc, dist, lim[0], lim[1], lim[2], lim[3]);
                if (dist < tolv1 && node->GetCADCurves().size() == 0)
                {
                    // Just give the node a parametric location
                    // DO NOT Project the node for the moment !
                    // Prioritise the CADCurve allocation if present !
                    node->SetCADSurf(surf, uv);
                }
            }
        }
    }
}
} // namespace Nektar::NekMesh
