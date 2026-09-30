////////////////////////////////////////////////////////////////////////////////
//
//  File: Pyramid.cpp
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
//  Description: Mesh pyramid object.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMesh/MeshElements/Pyramid.h>
#include <SpatialDomains/PyrGeom.h>

using namespace std;

namespace Nektar::NekMesh
{

LibUtilities::ShapeType Pyramid::type =
    GetElementFactory().RegisterCreatorFunction(LibUtilities::ePyramid,
                                                Pyramid::create, "Pyramid");

/// Vertex IDs that make up pyramid faces.
int Pyramid::m_faceIds[5][4] = {
    {0, 1, 2, 3}, {0, 1, 4, -1}, {1, 2, 4, -1}, {3, 2, 4, -1}, {0, 3, 4, -1}};

/**
 * @brief Create a pyramidic element.
 */
Pyramid::Pyramid(ElmtConfig pConf, vector<NodeSharedPtr> pNodeList,
                 vector<int> pTagList)
    : Element(pConf, GetNumNodes(pConf), pNodeList.size())
{
    m_tag     = "P";
    m_dim     = 3;
    m_taglist = pTagList;
    int n     = m_conf.m_order - 1;

    // This edge-node map is based on Nektar++ ordering.
    map<pair<int, int>, int> edgeNodeMap;
    map<pair<int, int>, int>::iterator it;
    edgeNodeMap[pair<int, int>(1, 2)] = 6;
    edgeNodeMap[pair<int, int>(2, 3)] = 6 + n;
    edgeNodeMap[pair<int, int>(4, 3)] = 6 + 2 * n;
    edgeNodeMap[pair<int, int>(1, 4)] = 6 + 3 * n;
    edgeNodeMap[pair<int, int>(1, 5)] = 6 + 4 * n;
    edgeNodeMap[pair<int, int>(2, 5)] = 6 + 5 * n;
    edgeNodeMap[pair<int, int>(3, 5)] = 6 + 6 * n;
    edgeNodeMap[pair<int, int>(4, 5)] = 6 + 7 * n;

    // Add vertices
    for (int i = 0; i < 5; ++i)
    {
        m_vertex.push_back(pNodeList[i]);
    }

    // Create edges (with corresponding set of edge points)
    int eid = 0;
    for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
    {
        vector<NodeSharedPtr> edgeNodes;
        if (m_conf.m_order > 1)
        {
            for (int j = it->second; j < it->second + n; ++j)
            {
                edgeNodes.push_back(pNodeList[j - 1]);
            }
        }
        m_edge.push_back(EdgeSharedPtr(new Edge(
            pNodeList[it->first.first - 1], pNodeList[it->first.second - 1],
            edgeNodes, m_conf.m_edgeCurveType)));
        m_edge.back()->m_id = eid++;
    }

    // Create faces
    int face_edges[5][4];
    int faceoffset = 5 + 8 * n;
    for (int j = 0; j < 5; ++j)
    {
        vector<NodeSharedPtr> faceVertices;
        vector<EdgeSharedPtr> faceEdges;
        vector<NodeSharedPtr> faceNodes;
        int nEdge = j > 0 ? 3 : 4;

        for (int k = 0; k < nEdge; ++k)
        {
            faceVertices.push_back(m_vertex[m_faceIds[j][k]]);
            NodeSharedPtr a = m_vertex[m_faceIds[j][k]];
            NodeSharedPtr b = m_vertex[m_faceIds[j][(k + 1) % nEdge]];
            for (unsigned int i = 0; i < m_edge.size(); ++i)
            {
                if ((m_edge[i]->m_n1 == a && m_edge[i]->m_n2 == b) ||
                    (m_edge[i]->m_n1 == b && m_edge[i]->m_n2 == a))
                {
                    faceEdges.push_back(m_edge[i]);
                    face_edges[j][k] = i;
                    break;
                }
            }
        }

        if (m_conf.m_faceNodes)
        {
            int facenodes = j == 0 ? n * n : n * (n - 1) / 2;
            for (int i = 0; i < facenodes; ++i)
            {
                faceNodes.push_back(pNodeList[faceoffset + i]);
            }
            faceoffset += facenodes;
        }

        // Try to translate between common face curve types
        LibUtilities::PointsType pType = m_conf.m_faceCurveType;

        if (pType == LibUtilities::ePolyEvenlySpaced && (j > 0))
        {
            pType = LibUtilities::eNodalTriEvenlySpaced;
        }

        m_face.push_back(
            FaceSharedPtr(new Face(faceVertices, faceNodes, faceEdges, pType)));
    }

    // Reorder edges to align with Nektar++ order.
    vector<EdgeSharedPtr> tmp(8);
    tmp[0] = m_edge[face_edges[0][0]];
    tmp[1] = m_edge[face_edges[0][1]];
    tmp[2] = m_edge[face_edges[0][2]];
    tmp[3] = m_edge[face_edges[0][3]];
    tmp[4] = m_edge[face_edges[1][2]];
    tmp[5] = m_edge[face_edges[1][1]];
    tmp[6] = m_edge[face_edges[3][1]];
    tmp[7] = m_edge[face_edges[3][2]];
    m_edge = tmp;
}

SpatialDomains::Geometry *Pyramid::GetGeom(int coordDim,
                                           SpatialDomains::EntityHolder &holder)
{
    std::array<SpatialDomains::Geometry2D *, 5> faces;

    for (int i = 0; i < 5; ++i)
    {
        faces[i] = m_face[i]->GetGeom(coordDim, holder);
    }

    SpatialDomains::PyrGeomUniquePtr pyr =
        ObjPoolManager<SpatialDomains::PyrGeom>::AllocateUniquePtr(m_id, faces);
    auto ret = dynamic_cast<SpatialDomains::Geometry *>(pyr.get());
    holder.m_pyrVec.push_back(std::move(pyr));

    ret->Setup();
    return ret;
}

/**
 * @brief Return the number of nodes defining a pyramid.
 */
unsigned int Pyramid::GetNumNodes(ElmtConfig pConf)
{
    int n = pConf.m_order;

    // valid for any order pyramid
    return (5             // corners
            + 8 * (n - 1) // mid-edge
            +
            pConf.m_faceNodes * ((n - 1) * (n - 1) +
                                 4 * (n - 1) * (n - 2) /
                                     2) // square base + 4xtriangle-number faces
            + pConf.m_volumeNodes * (n - 2) * (n - 1) * (2 * n - 3) /
                  6 // square pyramidal numbers
    );
}
} // namespace Nektar::NekMesh
