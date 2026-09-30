////////////////////////////////////////////////////////////////////////////////
//
//  File: LinearMeshGraph.hpp
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
//  Description: Create a LOR representation of a high-order MeshGraph.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_LINEARMESHGRAPH_HPP
#define NEKTAR_SPATIALDOMAINS_LINEARMESHGRAPH_HPP

#include <map>
#include <vector>

#include <SpatialDomains/MeshGraph.h>

namespace Nektar::SpatialDomains
{

template <class T> class Array3D
{
public:
    Array3D(std::size_t n0, std::size_t n1, std::size_t n2)
        : m_n0(n0), m_n1(n1), m_n2(n2), m_data(n0 * n1 * n2)
    {
    }

    T &operator()(std::size_t i, std::size_t j, std::size_t k)
    {
        assert(i < m_n0 && j < m_n1 && k < m_n2);
        return m_data[(i * m_n1 + j) * m_n2 + k];
    }

    const T &operator()(std::size_t i, std::size_t j, std::size_t k) const
    {
        assert(i < m_n0 && j < m_n1 && k < m_n2);
        return m_data[(i * m_n1 + j) * m_n2 + k];
    }

    T *data()
    {
        return m_data.data();
    }
    const T *data() const
    {
        return m_data.data();
    }

    std::size_t n0() const
    {
        return m_n0;
    }
    std::size_t n1() const
    {
        return m_n1;
    }
    std::size_t n2() const
    {
        return m_n2;
    }

private:
    std::size_t m_n0, m_n1, m_n2;
    std::vector<T> m_data;
};

class LinearMeshGraph
{
public:
    LinearMeshGraph(SpatialDomains::MeshGraphSharedPtr graph)
        : m_graph(graph), m_session(m_graph->GetSession()),
          m_meshDimension(m_graph->GetMeshDimension()),
          m_spaceDimension(m_graph->GetSpaceDimension())
    {
    }

    SpatialDomains::MeshGraphSharedPtr GetLinearGraph()
    {
        return m_linMesh;
    }

    SpatialDomains::MeshGraphSharedPtr CreateLinearGraph(
        int nsplit,
        std::map<int, std::pair<int, std::vector<int>>> &LinCoeffMap,
        bool UseGLL, bool useSimplex);

private:
    void LinMeshSetUp1DGeom(int nsplit, int voffset, bool UseGLL);

    void AddSplitEdge(int nsplit, int vid0, int vid1, int maxvertid, int edgeid,
                      Array<OneD, Array<OneD, NekDouble>> &Coords);

    void LinMeshSetUp2DGeom(int nsplit, std::map<int, int> &FceEdgOffset,
                            std::map<int, std::map<int, int>> &CoeffMap,
                            int voffset, bool UseGLL, bool useSimplex);

    void AddSplitTri(int nsplit, Array<TwoD, int> &vids,
                     Array<TwoD, int> &eids_h, Array<TwoD, int> &eids_v,
                     Array<TwoD, int> &eids_d, int maxedgeid, int edgeoffset,
                     int triid, std::map<int, int> &FceEdgOffset,
                     std::map<int, std::map<int, int>> &CoeffMap);

    void LinMeshSetUpTetGeom(int nsplit, std::map<int, int> &FceEdgOffset,
                             std::map<int, std::map<int, int>> &CoeffMap,
                             bool UseGLL);

    void LinMeshSetUpPrismGeom(int nsplit, std::map<int, int> &FceEdgOffset,
                               std::map<int, std::map<int, int>> &CoeffMap,
                               bool UseGLL);

    void LinMeshSetUpCompositesDomain(
        int nsplit,
        std::map<int, std::pair<int, std::vector<int>>> &LinCoeffMap,
        std::map<int, std::map<int, int>> &CoeffMap2D,
        std::map<int, std::map<int, int>> &CoeffMap3D, bool useSimplex);

    TiXmlElement *SetupLinearExpansionType(void)
    {
        TiXmlElement *expansionTypes =
            m_session->GetElement("NEKTAR/EXPANSIONS");

        TiXmlElement *newExpType = new TiXmlElement("NEKTAR/EXPANSIONS");

        // Find the Expansions tag
        ASSERTL0(expansionTypes, "Unable to find EXPANSIONS tag in file.");

        if (expansionTypes)
        {
            TiXmlElement *expansion = expansionTypes->FirstChildElement();
            ASSERTL0(expansion, "Unable to find entries in EXPANSIONS tag in"
                                "file.");
            std::string expType           = expansion->Value();
            std::vector<std::string> vars = m_session->GetVariables();

            if (expType == "E")
            {
                while (expansion)
                {
                    TiXmlElement *exp = new TiXmlElement("E");

                    std::string expstr = expansion->Attribute("COMPOSITE");
                    exp->SetAttribute("COMPOSITE", expstr);
                    exp->SetAttribute("NUMMODES", 2);
                    const char *rType = expansion->Attribute("TYPE");
                    if (rType) // use same type if exists
                    {
                        std::string typestr = expansion->Attribute("TYPE");
                        exp->SetAttribute("TYPE", typestr);
                    }
                    else // otherwise specify as modified. Not currently
                         // using BasisType
                    {
                        exp->SetAttribute("TYPE", "MODIFIED");
                    }
                    std::string fieldstr = expansion->Attribute("FIELDS");
                    exp->SetAttribute("FIELDS", fieldstr);

                    newExpType->LinkEndChild(exp);
                    expansion = expansion->NextSiblingElement("E");
                }
            }
            else
            {
                ASSERTL0(false, "Expansion type not defined");
            }
        }

        return newExpType;
    }

    // shuffle so highest vid is on edges 1,2 and return number of
    // offset rotations.
    int TriShuffleEdges(std::array<SegGeom *, 3> &edges,
                        std::array<SegGeom *, 3> *edgesSort = nullptr)
    {
        std::array<SegGeom *, 3> save;
        int maxVid0, maxVid1, startEd;

        if (edgesSort == nullptr)
        {
            maxVid0 = std::max(edges[0]->GetVid(0), edges[0]->GetVid(1));
            maxVid1 = std::max(edges[1]->GetVid(0), edges[1]->GetVid(1));

            startEd = (maxVid0 == maxVid1) ? 2 : (maxVid0 > maxVid1) ? 1 : 0;
        }
        else
        {
            maxVid0 = std::max((*edgesSort)[0]->GetVid(0),
                               (*edgesSort)[0]->GetVid(1));
            maxVid1 = std::max((*edgesSort)[1]->GetVid(0),
                               (*edgesSort)[1]->GetVid(1));
            startEd = (maxVid0 == maxVid1) ? 2 : (maxVid0 > maxVid1) ? 1 : 0;
        }

        // // swap start and 0
        save[0] = std::move(edges[0]);
        save[1] = std::move(edges[1]);
        save[2] = std::move(edges[2]);

        edges[0] = std::move(save[startEd]);
        edges[1] = std::move(save[(startEd + 1) % 3]);
        edges[2] = std::move(save[(startEd + 2) % 3]);

        return (3 - startEd) % 3;
    }

    /* Given face 0 and 1 give the other face ordering */
    const unsigned int TetFaceOrder[4][4][2] = {
        {{0, 0}, {2, 3}, {3, 1}, {1, 2}},
        {{3, 2}, {0, 0}, {0, 3}, {2, 0}},
        {{1, 3}, {3, 0}, {0, 0}, {0, 1}},
        {{2, 1}, {0, 2}, {1, 0}, {0, 0}}};

    /** reshuffles tet to point to singular vertices assuming faces
        are ordered in conventiaional manner */
    void TetShuffleFaces(std::array<TriGeom *, 4> &faces)
    {
        std::array<int, 4> faceid = {-1, -1, -1, -1};
        std::array<TriGeom *, 4> save;

        for (int i = 0; i < 4; ++i)
        {
            ASSERTL0(faces[i] != nullptr,
                     "Null face pointer at index " + std::to_string(i));
        }

        /** We know all faces are orientated so singular vertex on vertex 2
         **/
        faceid[0] = 0;
        /* face 0  */
        for (int i = 1; i < 4; ++i)
        {
            if (faces[i]->GetVid(2) < faces[faceid[0]]->GetVid(2))
            {
                faceid[0] = i;
            }
        }

        // /* face 1  */
        int eid = faces[faceid[0]]->GetEid(0);
        int i;
        for (i = 0; i < 4; ++i)
        {
            if ((i != faceid[0]) && (faces[i]->GetEid(0) == eid))
            {
                faceid[1] = i;
                break;
            }
        }
        ASSERTL0(i != 4, "Failed to find face 1 id in shuffle_faces: faces not "
                         "set up correctly");

        // // assign last two faces according to predefined map above.
        faceid[2] = TetFaceOrder[faceid[0]][faceid[1]][0];
        faceid[3] = TetFaceOrder[faceid[0]][faceid[1]][1];

        for (int i = 0; i < 4; ++i)
        {
            save[i] = std::move(faces[i]);
        }

        for (int i = 0; i < 4; ++i)
        {
            faces[i] = std::move(save[faceid[i]]);
        }
    }

    void PrismShuffleFaces(std::array<Geometry2D *, 5> &faces)
    {
        // assume second face defined singular vertex
        int edg0 = faces[1]->GetEdge(0)->GetGlobalID();

        // find face with this same edge
        int i, j;
        for (i = 0; i < 5; i += 2)
        {
            for (j = 0; j < 4; ++j)
            {
                if (faces[i]->GetEdge(j)->GetGlobalID() == edg0)
                {
                    break;
                }
            }
            if (j != 4)
            {
                break;
            }
        }
        ASSERTL1(i != 6, "Failed to find common edge");

        /* rotate faces  */

        if (i == 2)
        {
            Geometry2D *save0 = std::move(faces[0]);
            faces[0]          = std::move(faces[2]);
            faces[2]          = std::move(faces[4]);
            faces[4]          = std::move(save0);
        }
        else if (i == 4)
        {
            Geometry2D *save0 = std::move(faces[0]);
            faces[0]          = std::move(faces[4]);
            faces[4]          = std::move(faces[2]);
            faces[2]          = std::move(save0);
        }
    }

private:
    SpatialDomains::MeshGraphSharedPtr m_graph;
    SpatialDomains::MeshGraphSharedPtr m_linMesh;
    LibUtilities::SessionReaderSharedPtr m_session;
    int m_meshDimension;
    int m_spaceDimension;
};

SpatialDomains::MeshGraphSharedPtr LinearMeshGraph::CreateLinearGraph(
    int nsplit, std::map<int, std::pair<int, std::vector<int>>> &LinCoeffMap,
    bool UseGLL, bool useSimplex)
{
    //--------------------------------------------------------
    // initial setup
    //--------------------------------------------------------
    LibUtilities::CommSharedPtr comm = m_session->GetComm();
    ASSERTL0(comm.get(), "Communication not initialised.");

    // Populate SessionReader. This should be done only on the root process
    // so that we can partition appropriately without all processes having
    // to read in the input file.
    const bool isRoot = comm->TreatAsRankZero();
    std::string geomType;

    if (isRoot)
    {
        // Get geometry type, i.e. XML (compressed/uncompressed) or HDF5.
        geomType = m_session->GetGeometryType();

        // Convert to a vector of chars so that we can broadcast.
        std::vector<char> v(geomType.c_str(),
                            geomType.c_str() + geomType.length());

        size_t length = v.size();
        comm->Bcast(length, 0);
        comm->Bcast(v, 0);
    }
    else
    {
        size_t length;
        comm->Bcast(length, 0);

        std::vector<char> v(length);
        comm->Bcast(v, 0);

        geomType = std::string(v.begin(), v.end());
    }

    m_linMesh = MemoryManager<SpatialDomains::MeshGraph>::AllocateSharedPtr();

    // copy starting info
    m_linMesh->SetPartitionNumber(m_graph->GetPartitionNumber());
    m_linMesh->SetMeshPartitioned(m_graph->GetMeshPartitioned());
    m_linMesh->SetMeshDimension(m_meshDimension);
    m_linMesh->SetSpaceDimension(m_spaceDimension);
    m_linMesh->SetSession(m_graph->GetSession());

    //--------------------------------------------------------
    // Need to offset vertices so that new face vertids are lower than
    // original vertices, edges so that new orientation are aligned
    // (particularly for prisms).
    //--------------------------------------------------------
    int voffset      = 0;
    auto &prismGeoms = m_graph->GetGeomMap<PrismGeom>();
    if (m_meshDimension == 3 && prismGeoms.size() > 0)
    {
        // check for alignment of triangular faces or edge id of bottom edges
        for (auto [id, prism] : prismGeoms)
        {
            if (prism->GetEorient(1) != prism->GetEorient(3))
            {
                NEKERROR(ErrorUtil::ewarning,
                         "triangular faces are not aligned ");
            }
        }

        ASSERTL0(useSimplex == false, "Prisms are not currently "
                                      "setup for simplex LOR scheme");
    }

    //--------------------------------------------------------
    // Set up 1D elements
    //--------------------------------------------------------
    LinMeshSetUp1DGeom(nsplit, voffset, UseGLL);

    //--------------------------------------------------------
    // Set up 2D elements
    //--------------------------------------------------------
    std::map<int, std::map<int, int>> CoeffMap2D;
    std::map<int, int> FceEdgOffset; // map of face id and offset when shuffled

    LinMeshSetUp2DGeom(nsplit, FceEdgOffset, CoeffMap2D, voffset, UseGLL,
                       useSimplex);

    //--------------------------------------------------------
    // Set up 3D elements
    //--------------------------------------------------------
    std::map<int, std::map<int, int>> CoeffMap3D;
    if (m_meshDimension == 3)
    {
        LinMeshSetUpTetGeom(nsplit, FceEdgOffset, CoeffMap3D, UseGLL);
        LinMeshSetUpPrismGeom(nsplit, FceEdgOffset, CoeffMap3D, UseGLL);
    }

    //--------------------------------------------------------
    // Reset composites for BCs and domain definitions
    //--------------------------------------------------------
    LinMeshSetUpCompositesDomain(nsplit, LinCoeffMap, CoeffMap2D, CoeffMap3D,
                                 useSimplex);

    return m_linMesh;
}

void LinearMeshGraph::LinMeshSetUp1DGeom(int nsplit, int voffset, bool UseGLL)
{
    LibUtilities::CommSharedPtr comm = m_session->GetComm();

    /** set up a deep copy of existing vertices with voffset  **/
    for (auto [id, vert] : m_graph->GetGeomMap<PointGeom>())
    {
        int vid              = voffset + vert->GetGlobalID();
        PointGeomUniquePtr v = ObjPoolManager<PointGeom>::AllocateUniquePtr(
            m_spaceDimension, vid, (*vert)[0], (*vert)[1], (*vert)[2]);

        m_linMesh->AddGeom(vid, std::move(v));
    }

    int maxvertid = m_linMesh->GetGeomMap<PointGeom>().rbegin()->first + 1;
    comm->AllReduce(maxvertid, LibUtilities::ReduceMax);

    Array<OneD, Array<OneD, NekDouble>> Coords(3);
    Coords[0] = Array<OneD, NekDouble>(nsplit + 1);
    Coords[1] = Array<OneD, NekDouble>(nsplit + 1);
    Coords[2] = Array<OneD, NekDouble>(nsplit + 1);

    //--------------------------------------------------------
    // Split edges of macro mesh into new edges and vertices
    //--------------------------------------------------------

    for (auto [id, edgegeom] : m_graph->GetGeomMap<SegGeom>())
    {
        StdRegions::StdExpansionSharedPtr Xmap = edgegeom->GetXmap();

        Array<OneD, NekDouble> tmp(Xmap->GetTotPoints());
        if (UseGLL) // get GLL points
        {
            Xmap->BwdTrans(edgegeom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToGLL(tmp, Coords[0], nsplit + 1);
            Xmap->BwdTrans(edgegeom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToGLL(tmp, Coords[1], nsplit + 1);

            if (m_meshDimension == 3)
            {
                Xmap->BwdTrans(edgegeom->GetCoeffs(2), tmp);
                Xmap->PhysInterpToGLL(tmp, Coords[2], nsplit + 1);
            }
        }
        else // get equispaced points
        {
            Xmap->BwdTrans(edgegeom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Coords[0], nsplit + 1);
            Xmap->BwdTrans(edgegeom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Coords[1], nsplit + 1);

            if (m_meshDimension == 3)
            {
                Xmap->BwdTrans(edgegeom->GetCoeffs(2), tmp);
                Xmap->PhysInterpToSimplexEquiSpaced(tmp, Coords[2], nsplit + 1);
            }
        }

        // replace with method with api EdgeAdd(nsplit,vid0,vid1,
        // maxvertid, edgeid)
        int edgeid = edgegeom->GetGlobalID();
        AddSplitEdge(nsplit, edgegeom->GetVertex(0)->GetGlobalID(),
                     edgegeom->GetVertex(1)->GetGlobalID(), maxvertid, edgeid,
                     Coords);
    }
}

void LinearMeshGraph::AddSplitEdge(int nsplit, int vid0, int vid1,
                                   int maxvertid, int edgeid,
                                   Array<OneD, Array<OneD, NekDouble>> &Coords)
{
    for (int i = 1; i < nsplit; ++i)
    {
        int vid                 = maxvertid + edgeid * (nsplit - 1) + i - 1;
        PointGeomUniquePtr vert = ObjPoolManager<PointGeom>::AllocateUniquePtr(
            m_spaceDimension, vid, Coords[0][i], Coords[1][i], Coords[2][i]);
        m_linMesh->AddGeom(vid, std::move(vert));
    }

    // add SegGeoms
    for (int i = 0; i < nsplit; ++i)
    {
        int edgid = edgeid * nsplit + i;
        std::array<PointGeom *, 2> vert;

        if (i == 0)
        {
            vert[0] = m_linMesh->GetPointGeom(vid0);

            vert[1] = (nsplit == 1)
                          ? m_linMesh->GetPointGeom(vid1)
                          : m_linMesh->GetPointGeom(maxvertid +
                                                    edgeid * (nsplit - 1) + i);
        }
        else if (i == nsplit - 1)
        {
            vert[0] = m_linMesh->GetPointGeom(maxvertid +
                                              edgeid * (nsplit - 1) + i - 1);
            vert[1] = m_linMesh->GetPointGeom(vid1);
        }
        else
        {
            vert[0] = m_linMesh->GetPointGeom(maxvertid +
                                              edgeid * (nsplit - 1) + i - 1);
            vert[1] =
                m_linMesh->GetPointGeom(maxvertid + edgeid * (nsplit - 1) + i);
        }

        SegGeomUniquePtr edge = ObjPoolManager<SegGeom>::AllocateUniquePtr(
            edgid, m_spaceDimension, vert);

        m_linMesh->AddGeom(edgid, std::move(edge));
    }
}

void LinearMeshGraph::LinMeshSetUp2DGeom(
    int nsplit, std::map<int, int> &FceEdgOffset,
    std::map<int, std::map<int, int>> &CoeffMap, int voffset, bool UseGLL,
    bool useSimplex)
{
    LibUtilities::CommSharedPtr comm = m_session->GetComm();

    Array<OneD, NekDouble> Xpts, Ypts, Zpts;
    Xpts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1));
    Ypts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1));
    Zpts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1));

    // use max vertid
    int maxvertid =
        m_graph->GetGeomMap<PointGeom>().rbegin()->first + 1 + voffset;
    comm->AllReduce(maxvertid, LibUtilities::ReduceMax);

    int maxedgeid = m_linMesh->GetGeomMap<SegGeom>().rbegin()->first + 1;
    comm->AllReduce(maxedgeid, LibUtilities::ReduceMax);

    //--------------------------------------------------------
    // Split tri faces of macro mesh into new edges and vertices
    //--------------------------------------------------------
    Array<TwoD, int> vids(nsplit + 1, nsplit + 1, 0);
    // horirzotnal edges
    Array<TwoD, int> eids_h(nsplit + 1, nsplit, 0);
    // vertical edges
    Array<TwoD, int> eids_v(nsplit, nsplit + 1, 0);
    // diagonal edges
    Array<TwoD, int> eids_d(nsplit + 1, nsplit + 1, 0);

    // number of interior edges or quad if split into simplices ~ 3
    //  npslit*nsplit so need to use that here too
    int edgeoffset = nsplit * (3 * nsplit - 2);

    // add vertices, edges for tri geoms
    for (auto [id, geom] : m_graph->GetGeomMap<TriGeom>())
    {
        StdRegions::StdExpansionSharedPtr Xmap = geom->GetXmap();

        Array<OneD, NekDouble> tmp(Xmap->GetTotPoints());

        if (UseGLL) // get GLL points
        {
            Xmap->BwdTrans(geom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToGLL(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToGLL(tmp, Ypts, nsplit + 1);

            if (m_meshDimension == 3)
            {
                Xmap->BwdTrans(geom->GetCoeffs(2), tmp);
                Xmap->PhysInterpToGLL(tmp, Zpts, nsplit + 1);
            }
        }
        else // get equispaced points
        {
            // get equispaced points
            Xmap->BwdTrans(geom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Ypts, nsplit + 1);
            if (m_meshDimension == 3)
            {
                Xmap->BwdTrans(geom->GetCoeffs(2), tmp);
                Xmap->PhysInterpToSimplexEquiSpaced(tmp, Zpts, nsplit + 1);
            }
        }

        int triid = geom->GetGlobalID();

        int edgeid0 = geom->GetEdge(0)->GetGlobalID();
        int edgeid1 = geom->GetEdge(1)->GetGlobalID();
        int edgeid2 = geom->GetEdge(2)->GetGlobalID();

        // Determine if edges are cartesian aligned:
        bool fwd0 = geom->GetEorient(0) == StdRegions::eForwards;
        bool fwd1 = geom->GetEorient(1) == StdRegions::eForwards;
        bool fwd2 = geom->GetEorient(2) == StdRegions::eForwards;

        // fill in vertex ids along boundary of element  //
        for (int i = 1; i < nsplit; ++i)
        {
            vids[0][i] =
                fwd0 ? maxvertid + edgeid0 * (nsplit - 1) + i - 1
                     : maxvertid + edgeid0 * (nsplit - 1) + nsplit - i - 1;
            vids[i][nsplit - i] =
                fwd1 ? maxvertid + edgeid1 * (nsplit - 1) + i - 1
                     : maxvertid + edgeid1 * (nsplit - 1) + nsplit - i - 1;
            vids[i][0] =
                fwd2 ? maxvertid + edgeid2 * (nsplit - 1) + i - 1
                     : maxvertid + edgeid2 * (nsplit - 1) + nsplit - i - 1;
        }

        // set up corner ids
        vids[0][0] =
            fwd0 ? m_graph->GetSegGeom(edgeid0)->GetVertex(0)->GetGlobalID()
                 : m_graph->GetSegGeom(edgeid0)->GetVertex(1)->GetGlobalID();
        vids[0][nsplit] =
            fwd0 ? m_graph->GetSegGeom(edgeid0)->GetVertex(1)->GetGlobalID()
                 : m_graph->GetSegGeom(edgeid0)->GetVertex(0)->GetGlobalID();
        vids[nsplit][0] =
            fwd2 ? m_graph->GetSegGeom(edgeid2)->GetVertex(1)->GetGlobalID()
                 : m_graph->GetSegGeom(edgeid2)->GetVertex(0)->GetGlobalID();

        // fill in edge ids  along boundary of element
        for (int i = 0; i < nsplit; ++i)
        {
            eids_h[0][i] =
                fwd0 ? edgeid0 * nsplit + i : edgeid0 * nsplit + nsplit - i - 1;
            eids_v[i][0] =
                fwd2 ? edgeid2 * nsplit + i : edgeid2 * nsplit + nsplit - i - 1;
            eids_d[i][nsplit - i] =
                fwd1 ? edgeid1 * nsplit + i : edgeid1 * nsplit + nsplit - i - 1;
        }

        // add in vertices on interior of element
        int cnt = nsplit + 1;
        int vid;
        for (int j = 1; j < nsplit; ++j)
        {
            for (int i = 1; i < nsplit - j; ++i)
            {
                vid = triid * (nsplit - 1) * (nsplit - 1) +
                      (j - 1) * (nsplit - 1) + (i - 1);
                if (voffset == 0)
                {
                    vid += maxvertid + maxedgeid * (nsplit - 1) / nsplit;
                }

                vids[j][i] = vid;
                PointGeomUniquePtr vert =
                    ObjPoolManager<PointGeom>::AllocateUniquePtr(
                        m_spaceDimension, vid, Xpts[cnt + i], Ypts[cnt + i],
                        Zpts[cnt + i]);
                m_linMesh->AddGeom(vid, std::move(vert));
            }
            cnt += nsplit + 1 - j;
        }

        AddSplitTri(nsplit, vids, eids_h, eids_v, eids_d, maxedgeid, edgeoffset,
                    triid, FceEdgOffset, CoeffMap);
    }

    //--------------------------------------------------------
    // Split quad faces of macro mesh into new edges and vertices
    //--------------------------------------------------------

    // add vertices, edges for quad geoms
    for (auto [id, quadgeom] : m_graph->GetGeomMap<QuadGeom>())
    {
        StdRegions::StdExpansionSharedPtr Xmap = quadgeom->GetXmap();

        Array<OneD, NekDouble> tmp(Xmap->GetTotPoints());

        if (UseGLL) // get GLL points
        {
            Xmap->BwdTrans(quadgeom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToGLL(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(quadgeom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToGLL(tmp, Ypts, nsplit + 1);

            if (m_meshDimension == 3)
            {
                Xmap->BwdTrans(quadgeom->GetCoeffs(2), tmp);
                Xmap->PhysInterpToGLL(tmp, Zpts, nsplit + 1);
            }
        }
        else // get equispaced points
        {
            // get equispaced points
            Xmap->BwdTrans(quadgeom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(quadgeom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Ypts, nsplit + 1);

            if (m_meshDimension == 3)
            {
                Xmap->BwdTrans(quadgeom->GetCoeffs(2), tmp);
                Xmap->PhysInterpToSimplexEquiSpaced(tmp, Zpts, nsplit + 1);
            }
        }

        int quadid = quadgeom->GetGlobalID();

        int edgeid0 = quadgeom->GetEdge(0)->GetGlobalID();
        int edgeid1 = quadgeom->GetEdge(1)->GetGlobalID();
        int edgeid2 = quadgeom->GetEdge(2)->GetGlobalID();
        int edgeid3 = quadgeom->GetEdge(3)->GetGlobalID();

        // Determine if edges are cartesian aligned:
        bool fwd0 = quadgeom->GetEorient(0) == StdRegions::eForwards;
        bool fwd1 = quadgeom->GetEorient(1) == StdRegions::eForwards;
        bool fwd2 = quadgeom->GetEorient(2) == StdRegions::eForwards;
        bool fwd3 = quadgeom->GetEorient(3) == StdRegions::eForwards;

        // fill in vertex ids along boundary of element  //
        for (int i = 1; i < nsplit; ++i)
        {
            vids[0][i] =
                fwd0 ? maxvertid + edgeid0 * (nsplit - 1) + i - 1
                     : maxvertid + edgeid0 * (nsplit - 1) + nsplit - i - 1;
            vids[i][nsplit] =
                fwd1 ? maxvertid + edgeid1 * (nsplit - 1) + i - 1
                     : maxvertid + edgeid1 * (nsplit - 1) + nsplit - i - 1;

            // Need to check if fwd is anticlockwise or cartesian -
            // set up for cartesian
            vids[nsplit][i] =
                fwd2 ? maxvertid + edgeid2 * (nsplit - 1) + i - 1
                     : maxvertid + edgeid2 * (nsplit - 1) + nsplit - i - 1;
            vids[i][0] =
                fwd3 ? maxvertid + edgeid3 * (nsplit - 1) + i - 1
                     : maxvertid + edgeid3 * (nsplit - 1) + nsplit - i - 1;
        }

        // set up corner ids
        vids[0][0] =
            fwd0 ? m_graph->GetSegGeom(edgeid0)->GetVertex(0)->GetGlobalID()
                 : m_graph->GetSegGeom(edgeid0)->GetVertex(1)->GetGlobalID();
        vids[0][nsplit] =
            fwd0 ? m_graph->GetSegGeom(edgeid0)->GetVertex(1)->GetGlobalID()
                 : m_graph->GetSegGeom(edgeid0)->GetVertex(0)->GetGlobalID();
        vids[nsplit][0] =
            fwd2 ? m_graph->GetSegGeom(edgeid2)->GetVertex(0)->GetGlobalID()
                 : m_graph->GetSegGeom(edgeid2)->GetVertex(1)->GetGlobalID();
        vids[nsplit][nsplit] =
            fwd2 ? m_graph->GetSegGeom(edgeid2)->GetVertex(1)->GetGlobalID()
                 : m_graph->GetSegGeom(edgeid2)->GetVertex(0)->GetGlobalID();

        // fill in edge ids  along boundary of element
        for (int i = 0; i < nsplit; ++i)
        {
            eids_h[0][i] =
                fwd0 ? edgeid0 * nsplit + i : edgeid0 * nsplit + nsplit - i - 1;
            eids_h[nsplit][i] =
                fwd2 ? edgeid2 * nsplit + i : edgeid2 * nsplit + nsplit - i - 1;
            eids_v[i][nsplit] =
                fwd1 ? edgeid1 * nsplit + i : edgeid1 * nsplit + nsplit - i - 1;
            eids_v[i][0] =
                fwd3 ? edgeid3 * nsplit + i : edgeid3 * nsplit + nsplit - i - 1;
        }

        // add in vertices on interior of element
        int vid;
        for (int j = 1; j < nsplit; ++j)
        {
            for (int i = 1; i < nsplit; ++i)
            {
                vid = quadid * (nsplit - 1) * (nsplit - 1) +
                      (j - 1) * (nsplit - 1) + (i - 1);

                if (voffset == 0)
                {
                    vid += maxvertid + maxedgeid * (nsplit - 1) / nsplit;
                }

                vids[j][i] = vid;
                PointGeomUniquePtr vert =
                    ObjPoolManager<PointGeom>::AllocateUniquePtr(
                        m_spaceDimension, vid, Xpts[j * (nsplit + 1) + i],
                        Ypts[j * (nsplit + 1) + i], Zpts[j * (nsplit + 1) + i]);
                m_linMesh->AddGeom(vid, std::move(vert));
            }
        }

        // Add in horizontal edges
        // KK - we might not need to create explicitly the edges now
        int cnt = 0;
        for (int j = 1; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit; ++i)
            {
                int edgid = maxedgeid + quadid * edgeoffset + cnt++;

                eids_h[j][i] = edgid;

                std::array<PointGeom *, 2> vert;
                vert[0] = m_linMesh->GetPointGeom(vids[j][i]);
                vert[1] = m_linMesh->GetPointGeom(vids[j][i + 1]);

                SegGeomUniquePtr edge =
                    ObjPoolManager<SegGeom>::AllocateUniquePtr(
                        edgid, m_spaceDimension, vert);

                m_linMesh->AddGeom(edgid, std::move(edge));
            }
        }

        // Add in vertical edges
        // KK - we might not need to create explicitly the edges now
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 1; i < nsplit; ++i)
            {
                int edgid = maxedgeid + quadid * edgeoffset + cnt++;

                eids_v[j][i] = edgid;

                std::array<PointGeom *, 2> vert;
                vert[0] = m_linMesh->GetPointGeom(vids[j][i]);
                vert[1] = m_linMesh->GetPointGeom(vids[j + 1][i]);
                SegGeomUniquePtr edge =
                    ObjPoolManager<SegGeom>::AllocateUniquePtr(
                        edgid, m_spaceDimension, vert);

                m_linMesh->AddGeom(edgid, std::move(edge));
            }
        }

        if (useSimplex)
        {
            // Add in vertical and diagonal edges
            for (int j = 0; j < nsplit; ++j)
            {
                for (int i = 0; i < nsplit; ++i)
                {
                    int edgid = maxedgeid + quadid * edgeoffset + cnt++;

                    eids_d[j][i] = edgid;

                    std::array<PointGeom *, 2> vert;
                    vert[0] = m_linMesh->GetPointGeom(vids[j][i + 1]);
                    vert[1] = m_linMesh->GetPointGeom(vids[j + 1][i]);
                    SegGeomUniquePtr edge_d =
                        ObjPoolManager<SegGeom>::AllocateUniquePtr(
                            edgid, m_spaceDimension, vert);
                    m_linMesh->AddGeom(edgid, std::move(edge_d));
                }
            }

            // add in quad elements
            for (int j = 0; j < nsplit; ++j)
            {
                for (int i = 0; i < nsplit; ++i)
                {
                    int newtriid =
                        2 * (quadid * nsplit * nsplit + j * nsplit + i);

                    std::array<SegGeom *, 3> edges;
                    edges[0]             = m_linMesh->GetSegGeom(eids_h[j][i]);
                    edges[1]             = m_linMesh->GetSegGeom(eids_d[j][i]);
                    edges[2]             = m_linMesh->GetSegGeom(eids_v[j][i]);
                    TriGeomUniquePtr tri = ObjPoolManager<
                        SpatialDomains::TriGeom>::AllocateUniquePtr(newtriid,
                                                                    edges);
                    m_linMesh->AddGeom(newtriid, std::move(tri));

                    newtriid =
                        2 * (quadid * nsplit * nsplit + j * nsplit + i) + 1;

                    edges[0] = m_linMesh->GetSegGeom(eids_h[j + 1][i]);
                    edges[1] = m_linMesh->GetSegGeom(eids_d[j][i]);
                    edges[2] = m_linMesh->GetSegGeom(eids_v[j][i + 1]);

                    TriGeomUniquePtr tri1 = ObjPoolManager<
                        SpatialDomains::TriGeom>::AllocateUniquePtr(newtriid,
                                                                    edges);

                    m_linMesh->AddGeom(newtriid, std::move(tri1));
                }
            }
        }
        else
        {
            // add in elements
            for (int j = 0; j < nsplit; ++j)
            {
                for (int i = 0; i < nsplit; ++i)
                {
                    int newquadid =
                        2 * quadid * nsplit * nsplit + j * nsplit + i;

                    std::array<SegGeom *, 4> edges;
                    edges[0] = m_linMesh->GetSegGeom(eids_h[j][i]);
                    edges[2] = m_linMesh->GetSegGeom(eids_h[j + 1][i]);
                    edges[1] = m_linMesh->GetSegGeom(eids_v[j][i + 1]);
                    edges[3] = m_linMesh->GetSegGeom(eids_v[j][i]);

                    QuadGeomUniquePtr quad = ObjPoolManager<
                        SpatialDomains::QuadGeom>::AllocateUniquePtr(newquadid,
                                                                     edges);
                    m_linMesh->AddGeom(newquadid, std::move(quad));
                }
            }
        }
    }
}

void LinearMeshGraph::AddSplitTri(int nsplit, Array<TwoD, int> &vids,
                                  Array<TwoD, int> &eids_h,
                                  Array<TwoD, int> &eids_v,
                                  Array<TwoD, int> &eids_d, int maxedgeid,
                                  int edgeoffset, int triid,
                                  std::map<int, int> &FceEdgOffset,
                                  std::map<int, std::map<int, int>> &CoeffMap)

{
    // Add in horizontal edges
    int cnt = 0;
    for (int j = 1; j < nsplit; ++j)
    {
        for (int i = 0; i < nsplit - j; ++i)
        {
            int edgid    = maxedgeid + triid * edgeoffset + cnt++;
            eids_h[j][i] = edgid;

            std::array<PointGeom *, 2> vert;
            vert[0] = m_linMesh->GetPointGeom(vids[j][i]);
            vert[1] = m_linMesh->GetPointGeom(vids[j][i + 1]);

            SegGeomUniquePtr edge = ObjPoolManager<SegGeom>::AllocateUniquePtr(
                edgid, m_spaceDimension, vert);

            m_linMesh->AddGeom(edgid, std::move(edge));
        }
    }

    // Add in vertical and diagonal edges
    for (int j = 0; j < nsplit; ++j)
    {
        for (int i = 1; i < nsplit - j; ++i)
        {
            int edgid = maxedgeid + triid * edgeoffset + cnt++;

            eids_v[j][i] = edgid;

            std::array<PointGeom *, 2> vert;
            vert[0]               = m_linMesh->GetPointGeom(vids[j][i]);
            vert[1]               = m_linMesh->GetPointGeom(vids[j + 1][i]);
            SegGeomUniquePtr edge = ObjPoolManager<SegGeom>::AllocateUniquePtr(
                edgid, m_spaceDimension, vert);

            m_linMesh->AddGeom(edgid, std::move(edge));

            edgid = maxedgeid + triid * edgeoffset + cnt++;

            eids_d[j][i] = edgid;

            vert[0] = m_linMesh->GetPointGeom(vids[j][i]);
            vert[1] = m_linMesh->GetPointGeom(vids[j + 1][i - 1]);
            SegGeomUniquePtr edge_d =
                ObjPoolManager<SegGeom>::AllocateUniquePtr(
                    edgid, m_spaceDimension, vert);

            m_linMesh->AddGeom(edgid, std::move(edge_d));
        }
    }

    // add in tri elements
    cnt = 0;
    for (int j = 0; j < nsplit; ++j)
    {
        // lower triangular elements.
        for (int i = 0; i < nsplit - j; ++i)
        {
            int newtriid = 2 * triid * nsplit * nsplit + cnt++;

            std::array<SegGeom *, 3> edges;
            edges[0] = m_linMesh->GetSegGeom(eids_h[j][i]);
            edges[1] = m_linMesh->GetSegGeom(eids_d[j][i + 1]);
            edges[2] = m_linMesh->GetSegGeom(eids_v[j][i]);
            // shuffle so lowest vid is between  edges 1,2
            FceEdgOffset[newtriid] = TriShuffleEdges(edges);

            TriGeomUniquePtr tri =
                ObjPoolManager<SpatialDomains::TriGeom>::AllocateUniquePtr(
                    newtriid, edges);

            m_linMesh->AddGeom(newtriid, std::move(tri));
        }

        // upper triangular elements.
        for (int i = 1; i < nsplit - j; ++i)
        {
            int newtriid = 2 * triid * nsplit * nsplit + cnt++;

            std::array<SegGeom *, 3> edges;
            edges[0] = m_linMesh->GetSegGeom(eids_d[j][i]);
            edges[1] = m_linMesh->GetSegGeom(eids_v[j][i]);
            edges[2] = m_linMesh->GetSegGeom(eids_h[j + 1][i - 1]);
            // shuffle so lowest vid is on edges 1,2
            FceEdgOffset[newtriid] = TriShuffleEdges(edges);
            TriGeomUniquePtr tri =
                ObjPoolManager<SpatialDomains::TriGeom>::AllocateUniquePtr(
                    newtriid, edges);
            m_linMesh->AddGeom(newtriid, std::move(tri));
        }
    }

    // collect a mapping of Vid to point vertex location
    cnt = 0;
    for (int j = 0; j < nsplit + 1; ++j)
    {
        for (int i = 0; i < nsplit + 1 - j; ++i)
        {
            CoeffMap[triid][vids[j][i]] = cnt++;
        }
    }
}

void LinearMeshGraph::LinMeshSetUpTetGeom(
    int nsplit, std::map<int, int> &FceEdgOffset,
    std::map<int, std::map<int, int>> &CoeffMap, bool UseGLL)
{
    // get vertex ids on faces of 3D shape
    Array3D<int> E_vid(nsplit + 1, nsplit + 1, nsplit + 1);
    // X-aligned edge
    Array3D<int> E_eidx(nsplit + 1, nsplit + 1, nsplit + 1);
    // Y-aligned edge
    Array3D<int> E_eidy(nsplit + 1, nsplit + 1, nsplit + 1);
    // Z-aligned edge
    Array3D<int> E_eidz(nsplit + 1, nsplit + 1, nsplit + 1);

    // XY-diagonal edge
    Array3D<int> E_eidxy(nsplit + 1, nsplit + 1, nsplit + 1);
    // ZY-diagonal edge
    Array3D<int> E_eidyz(nsplit + 1, nsplit + 1, nsplit + 1);
    // XZ-diagonal edge
    Array3D<int> E_eidxz(nsplit + 1, nsplit + 1, nsplit + 1);

    // Fill with xy face ids (lower triangular if simplex)
    Array3D<int> E_fidxy1(nsplit + 1, nsplit + 1, nsplit + 1);
    // Fill with xy face ids (upperr triangular if simplex)
    Array3D<int> E_fidxy2(nsplit + 1, nsplit + 1, nsplit + 1);

    // Fill with xz face ids (lower triangular if simplex)
    Array3D<int> E_fidxz1(nsplit + 1, nsplit + 1, nsplit + 1);
    // Fill with xz face ids (upperr triangular if simplex)
    Array3D<int> E_fidxz2(nsplit + 1, nsplit + 1, nsplit + 1);

    // Fill with yz face ids (lower triangular if simplex)
    Array3D<int> E_fidyz1(nsplit + 1, nsplit + 1, nsplit + 1);
    // Fill with yz face ids (upperr triangular if simplex)
    Array3D<int> E_fidyz2(nsplit + 1, nsplit + 1, nsplit + 1);

    // Fill with yz diagonal face ids (lower triangular if simplex)
    Array3D<int> E_fidyzd1(nsplit + 1, nsplit + 1, nsplit + 1);
    // Fill with yz diagonal face ids (upperr triangular if simplex)
    Array3D<int> E_fidyzd2(nsplit + 1, nsplit + 1, nsplit + 1);

    int maxvertid, maxedgeid;
    LibUtilities::CommSharedPtr comm = m_session->GetComm();

    Array<OneD, NekDouble> Xpts, Ypts, Zpts;
    Xpts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1) * (nsplit + 1));
    Ypts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1) * (nsplit + 1));
    Zpts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1) * (nsplit + 1));

    /* reset maxedgeid and maxvertid  */
    maxvertid = m_linMesh->GetGeomMap<PointGeom>().rbegin()->first + 1;
    comm->AllReduce(maxvertid, LibUtilities::ReduceMax);
    maxedgeid = m_linMesh->GetGeomMap<SegGeom>().rbegin()->first + 1;
    comm->AllReduce(maxedgeid, LibUtilities::ReduceMax);

    int maxfaceid = 0;

    if (m_linMesh->GetGeomMap<TriGeom>().size())
    {
        maxfaceid = std::max(
            maxfaceid, m_linMesh->GetGeomMap<TriGeom>().rbegin()->first + 1);
    }
    if (m_linMesh->GetGeomMap<QuadGeom>().size())
    {
        maxfaceid = std::max(
            maxfaceid, m_linMesh->GetGeomMap<QuadGeom>().rbegin()->first + 1);
    }

    comm->AllReduce(maxfaceid, LibUtilities::ReduceMax);
    //--------------------------------------------------------
    // Split tet elements of macro mesh into new edges and vertices
    //--------------------------------------------------------

    int edgeoffset = 3 * nsplit * (nsplit - 1) * (nsplit - 1);
    int faceoffset = 3 * nsplit * nsplit * (nsplit - 1);

    // For every original tet keep a map between the vid of the split
    // tet and the location of the points in the macro nodes to create
    // the LinCoeffMap when required. Needed since we reshuffle tet
    // orientaiton so cannot infer pattern directly
    for (auto [id, geom] : m_graph->GetGeomMap<TetGeom>())
    {
        int vcnt   = 0;
        int ecnt   = 0;
        int fcnt   = 0;
        int elmtid = geom->GetGlobalID();

        /* face 0  */
        int fid        = geom->GetFace(0)->GetGlobalID();
        int fce_offset = 2 * fid * nsplit * nsplit;
        bool FceFwd =
            geom->GetForient(0) == StdRegions::eDir1FwdDir1_Dir2FwdDir2;

        int cnt = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit - j; ++i)
            {
                int Fce_id        = FceFwd ? fce_offset + cnt + i
                                           : fce_offset + cnt + nsplit - j - 1 - i;
                E_fidxy1(0, j, i) = Fce_id;

                // due to edge shuffling earlier need to find new
                // orientation (singlar vertex at top) to fill in
                // stored edges and vertices
                int oset = FceEdgOffset[Fce_id];

                TriGeom *face = m_linMesh->GetTriGeom(Fce_id);
                // fill in vertices & edges from faces1
                if (FceFwd)
                {
                    E_vid(0, j, i)     = face->GetVid(oset);
                    E_vid(0, j, i + 1) = face->GetVid((oset + 1) % 3);

                    E_eidy(0, j, i)  = face->GetEid((oset + 2) % 3);
                    E_eidxy(0, j, i) = face->GetEid((oset + 1) % 3);
                }
                else
                {
                    E_vid(0, j, i)     = face->GetVid((oset + 1) % 3);
                    E_vid(0, j, i + 1) = face->GetVid(oset);

                    E_eidy(0, j, i)  = face->GetEid((oset + 1) % 3);
                    E_eidxy(0, j, i) = face->GetEid((oset + 2) % 3);
                }
                E_vid(0, j + 1, i) = face->GetVid((oset + 2) % 3);
                E_eidx(0, j, i)    = face->GetEid(oset);
            }
            cnt += nsplit - j;

            for (int i = 0; i < nsplit - j - 1; ++i)
            {
                int Fce_id        = FceFwd ? fce_offset + cnt + i
                                           : fce_offset + cnt + nsplit - j - 2 - i;
                E_fidxy2(0, j, i) = Fce_id;
            }
            cnt += nsplit - j - 1;
        }

        /* face 1  */
        fid        = geom->GetFace(1)->GetGlobalID();
        fce_offset = 2 * fid * nsplit * nsplit;
        FceFwd = geom->GetForient(1) == StdRegions::eDir1FwdDir1_Dir2FwdDir2;
        cnt    = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit - j; ++i)
            {
                int Fce_id        = FceFwd ? fce_offset + cnt + i
                                           : fce_offset + cnt + nsplit - j - 1 - i;
                E_fidxz1(j, 0, i) = Fce_id;

                // due to edge shuffling earlier need to identify know
                // orientation (singlar vertex at top) orientation to
                // fill in stored edges and vertices
                int oset = FceEdgOffset[Fce_id];

                TriGeom *face = m_linMesh->GetTriGeom(Fce_id);
                // fill in vertices & edges from faces
                if (FceFwd)
                {
                    E_vid(j, 0, i)     = face->GetVid(oset);
                    E_vid(j, 0, i + 1) = face->GetVid((oset + 1) % 3);

                    E_eidz(j, 0, i)  = face->GetEid((oset + 2) % 3);
                    E_eidxz(j, 0, i) = face->GetEid((oset + 1) % 3);
                }
                else
                {
                    E_vid(j, 0, i)     = face->GetVid((oset + 1) % 3);
                    E_vid(j, 0, i + 1) = face->GetVid(oset);

                    E_eidz(j, 0, i)  = face->GetEid((oset + 1) % 3);
                    E_eidxz(j, 0, i) = face->GetEid((oset + 2) % 3);
                }
                E_vid(j + 1, 0, i) = face->GetVid((oset + 2) % 3);
                E_eidx(j, 0, i)    = face->GetEid(oset);
            }
            cnt += nsplit - j;

            for (int i = 0; i < nsplit - j - 1; ++i)
            {
                int Fce_id        = FceFwd ? fce_offset + cnt + i
                                           : fce_offset + cnt + nsplit - j - 2 - i;
                E_fidxz2(j, 0, i) = Fce_id;
            }
            cnt += nsplit - j - 1;
        }

        /* face 2  */
        fid        = geom->GetFace(2)->GetGlobalID();
        fce_offset = 2 * fid * nsplit * nsplit;
        FceFwd = geom->GetForient(2) == StdRegions::eDir1FwdDir1_Dir2FwdDir2;

        cnt = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit - j; ++i)
            {
                int Fce_id = FceFwd ? fce_offset + cnt + i
                                    : fce_offset + cnt + nsplit - j - 1 - i;

                E_fidyzd1(j, i, nsplit - j - i - 1) = Fce_id;

                // due to edge shuffling earlier need to identify know
                // orientation (singlar vertex at top) orientation to
                // fill in stored edges and vertices
                int oset = FceEdgOffset[Fce_id];

                TriGeom *face = m_linMesh->GetTriGeom(Fce_id);
                // fill in vertices & edges from faces
                if (FceFwd)
                {
                    E_vid(j, i, nsplit - j - i) = face->GetVid(oset);
                    E_vid(j, i + 1, nsplit - j - i - 1) =
                        face->GetVid((oset + 1) % 3);

                    E_eidxz(j, i, nsplit - j - i - 1) =
                        face->GetEid((oset + 2) % 3);
                    E_eidyz(j, i, nsplit - j - i - 1) =
                        face->GetEid((oset + 1) % 3);
                }
                else
                {
                    E_vid(j, i, nsplit - j - i) = face->GetVid((oset + 1) % 3);
                    E_vid(j, i + 1, nsplit - j - i - 1) = face->GetVid(oset);

                    E_eidxz(j, i, nsplit - j - i - 1) =
                        face->GetEid((oset + 1) % 3);
                    E_eidyz(j, i, nsplit - j - i - 1) =
                        face->GetEid((oset + 2) % 3);
                }

                E_vid(j + 1, i, nsplit - j - i - 1) =
                    face->GetVid((oset + 2) % 3);
                E_eidxy(j, i, nsplit - j - i - 1) = face->GetEid(oset);

                // reorient face if required. Note cannto do earlier
                // since need to setup the arrays above
                // LinMesh->m_triGeoms[Fce_id] = ReOrientFace(face);
            }
            cnt += nsplit - j;

            for (int i = 0; i < nsplit - j - 1; ++i)
            {
                int Fce_id = FceFwd ? fce_offset + cnt + i
                                    : fce_offset + cnt + nsplit - j - 2 - i;
                E_fidyzd2(j, i, nsplit - j - i - 2) = Fce_id;
            }
            cnt += nsplit - j - 1;
        }

        /* face 3  */
        fid        = geom->GetFace(3)->GetGlobalID();
        fce_offset = 2 * fid * nsplit * nsplit;
        FceFwd = geom->GetForient(3) == StdRegions::eDir1FwdDir1_Dir2FwdDir2;
        cnt    = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit - j; ++i)
            {
                int Fce_id        = FceFwd ? fce_offset + cnt + i
                                           : fce_offset + cnt + nsplit - j - 1 - i;
                E_fidyz1(j, i, 0) = Fce_id;

                // due to edge shuffling earlier need to identify know
                // orientation (singlar vertex at top) orientation to
                // fill in stored edges and vertices
                int oset = FceEdgOffset[Fce_id];

                TriGeom *face = m_linMesh->GetTriGeom(Fce_id);
                // fill in vertices & edges from faces
                if (FceFwd)
                {
                    E_vid(j, i, 0)     = face->GetVid(oset);
                    E_vid(j, i + 1, 0) = face->GetVid((oset + 1) % 3);

                    E_eidz(j, i, 0)  = face->GetEid((oset + 2) % 3);
                    E_eidyz(j, i, 0) = face->GetEid((oset + 1) % 3);
                }
                else
                {
                    E_vid(j, i, 0)     = face->GetVid((oset + 1) % 3);
                    E_vid(j, i + 1, 0) = face->GetVid(oset);

                    E_eidz(j, i, 0)  = face->GetEid((oset + 1) % 3);
                    E_eidyz(j, i, 0) = face->GetEid((oset + 2) % 3);
                }
                E_vid(j + 1, i, 0) = face->GetVid((oset + 2) % 3);
                E_eidy(j, i, 0)    = face->GetEid(oset);
            }
            cnt += nsplit - j;

            for (int i = 0; i < nsplit - j - 1; ++i)
            {
                int Fce_id        = FceFwd ? fce_offset + cnt + i
                                           : fce_offset + cnt + nsplit - j - 2 - i;
                E_fidyz2(j, i, 0) = Fce_id;
            }
            cnt += nsplit - j - 1;
        }

        StdRegions::StdExpansionSharedPtr Xmap = geom->GetXmap();

        Array<OneD, NekDouble> tmp(Xmap->GetTotPoints());

        // fill in interior vids if exist
        if (UseGLL) // get GLL points
        {
            Xmap->BwdTrans(geom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToGLL(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToGLL(tmp, Ypts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(2), tmp);
            Xmap->PhysInterpToGLL(tmp, Zpts, nsplit + 1);
        }
        else // get equispaced points
        {
            Xmap->BwdTrans(geom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Ypts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(2), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Zpts, nsplit + 1);
        }

        cnt = (nsplit + 2) * (nsplit + 1) / 2;
        for (int k = 1; k < nsplit; ++k)
        {
            cnt += nsplit + 1 - k; // first row
            for (int j = 1; j < nsplit - k; ++j)
            {
                cnt++;
                for (int i = 1; i < nsplit - k - j; ++i, cnt++)
                {
                    int vid =
                        maxvertid +
                        elmtid * (nsplit - 1) * (nsplit - 1) * (nsplit - 1) +
                        vcnt++;

                    E_vid(k, j, i) = vid;

                    PointGeomUniquePtr vert =
                        ObjPoolManager<PointGeom>::AllocateUniquePtr(
                            m_spaceDimension, vid, Xpts[cnt], Ypts[cnt],
                            Zpts[cnt]);

                    m_linMesh->AddGeom(vid, std::move(vert));
                }
                cnt++;
            }
            cnt++;
        }

        // declare Tets
        cnt = 0;
        for (int k = 0; k < nsplit; ++k)
        {
            for (int j = 0; j < nsplit - k; ++j)
            {
                for (int i = 0; i < nsplit - k - j; ++i)
                {
                    if (i < nsplit - k - j - 2)
                    {
                        /* addd in 6 new edges */
                        // edge 0
                        int edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        std::array<PointGeom *, 2> vert;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k, j + 1, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i + 1));
                        SegGeomUniquePtr edge0 =
                            ObjPoolManager<SegGeom>::AllocateUniquePtr(
                                edgid, m_spaceDimension, vert);
                        m_linMesh->AddGeom(edgid, std::move(edge0));
                        E_eidz(k, j + 1, i + 1) = edgid;

                        // edge 1
                        edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i + 1));
                        SegGeomUniquePtr edge1 =
                            ObjPoolManager<SegGeom>::AllocateUniquePtr(
                                edgid, m_spaceDimension, vert);
                        m_linMesh->AddGeom(edgid, std::move(edge1));
                        E_eidx(k + 1, j + 1, i) = edgid;

                        // edge 2
                        edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i + 1));
                        SegGeomUniquePtr edge2 =
                            ObjPoolManager<SegGeom>::AllocateUniquePtr(
                                edgid, m_spaceDimension, vert);
                        m_linMesh->AddGeom(edgid, std::move(edge2));
                        E_eidy(k + 1, j, i + 1) = edgid;

                        // edge 3
                        edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k, j + 1, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j, i + 1));
                        SegGeomUniquePtr edge3 =
                            ObjPoolManager<SegGeom>::AllocateUniquePtr(
                                edgid, m_spaceDimension, vert);
                        m_linMesh->AddGeom(edgid, std::move(edge3));
                        E_eidyz(k, j, i + 1) = edgid;

                        // edge 4
                        edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k, j + 1, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i));
                        SegGeomUniquePtr edge4 =
                            ObjPoolManager<SegGeom>::AllocateUniquePtr(
                                edgid, m_spaceDimension, vert);
                        m_linMesh->AddGeom(edgid, std::move(edge4));
                        E_eidxz(k, j + 1, i) = edgid;

                        // edge 5
                        edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i));
                        SegGeomUniquePtr edge5 =
                            ObjPoolManager<SegGeom>::AllocateUniquePtr(
                                edgid, m_spaceDimension, vert);
                        m_linMesh->AddGeom(edgid, std::move(edge5));
                        E_eidxy(k + 1, j, i) = edgid;

                        // /* Add in 4 new faces */

                        /* Tri0 = bottom left Tri */
                        int newtriid = elmtid * faceoffset + maxfaceid + fcnt++;
                        std::array<SegGeom *, 3> edges;
                        edges[0] = m_linMesh->GetSegGeom(E_eidyz(k, j, i + 1));
                        edges[1] = m_linMesh->GetSegGeom(E_eidxz(k, j + 1, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidxy(k + 1, j, i));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri0 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid, edges);

                        m_linMesh->AddGeom(newtriid, std::move(tri0));
                        E_fidyzd2(k, j, i) = newtriid;

                        /* Tri1 =  */
                        newtriid = elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidyz(k, j, i + 1));
                        edges[1] =
                            m_linMesh->GetSegGeom(E_eidz(k, j + 1, i + 1));
                        edges[2] =
                            m_linMesh->GetSegGeom(E_eidy(k + 1, j, i + 1));
                        // shuffle so lowest vid is on edges 1,  2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri1 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid, edges);
                        m_linMesh->AddGeom(newtriid, std::move(tri1));
                        E_fidyz2(k, j, i + 1) = newtriid;

                        /* Tri2 =  */
                        newtriid = elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidxz(k, j + 1, i));
                        edges[1] =
                            m_linMesh->GetSegGeom(E_eidz(k, j + 1, i + 1));
                        edges[2] =
                            m_linMesh->GetSegGeom(E_eidx(k + 1, j + 1, i));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri2 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid, edges);
                        m_linMesh->AddGeom(newtriid, std::move(tri2));
                        E_fidxz2(k, j + 1, i) = newtriid;

                        /* Tri3 =  */
                        newtriid = elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidxy(k + 1, j, i));
                        edges[1] =
                            m_linMesh->GetSegGeom(E_eidy(k + 1, j, i + 1));
                        edges[2] =
                            m_linMesh->GetSegGeom(E_eidx(k + 1, j + 1, i));
                        // shuffle so lowest vid is on
                        TriShuffleEdges(edges);
                        // edges 1, 2
                        TriGeomUniquePtr tri3 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid, edges);
                        m_linMesh->AddGeom(newtriid, std::move(tri3));
                        E_fidxy2(k + 1, j, i) = newtriid;

                        /* add new top  Corner  Tet */
                        int newtetid =
                            elmtid * nsplit * nsplit * nsplit + cnt++;
                        std::array<TriGeom *, 4> faces;
                        faces[0] =
                            m_linMesh->GetTriGeom(E_fidyzd2(k, j, i)); // tri0
                        faces[1] = m_linMesh->GetTriGeom(
                            E_fidyz2(k, j, i + 1)); // tri1
                        faces[2] = m_linMesh->GetTriGeom(
                            E_fidxz2(k, j + 1, i)); // tri2
                        faces[3] = m_linMesh->GetTriGeom(
                            E_fidxy2(k + 1, j, i)); // tri3
                        TetShuffleFaces(faces);     // shuffle so max
                        // vert id at top
                        TetGeomUniquePtr tet =
                            ObjPoolManager<SpatialDomains::TetGeom>::
                                AllocateUniquePtr(newtetid, faces);
                        m_linMesh->AddGeom(newtetid, std::move(tet));
                        m_linMesh->PopulateFaceToElMap(
                            m_linMesh->GetTetGeom(newtetid), TetGeom::kNfaces);
                    }

                    if (i < nsplit - k - j - 1)
                    {
                        /* addd in new edge (off
                           diagonal) */
                        int edgid = elmtid * edgeoffset + maxedgeid + ecnt++;

                        std::array<PointGeom *, 2> vert;
                        vert[0] = m_linMesh->GetPointGeom(E_vid(k, j + 1, i));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j, i + 1));

                        SegGeomUniquePtr edge =
                            ObjPoolManager<SpatialDomains::SegGeom>::
                                AllocateUniquePtr(edgid, m_spaceDimension,
                                                  vert);
                        m_linMesh->AddGeom(edgid, std::move(edge));

                        /* Add in 8 new faces */

                        /* Tri0 = bottom left Tri */
                        int newtriid0 =
                            elmtid * faceoffset + maxfaceid + fcnt++;

                        std::array<SegGeom *, 3> edges;
                        edges[0] = m_linMesh->GetSegGeom(E_eidxy(k, j, i));
                        edges[1] = m_linMesh->GetSegGeom(E_eidyz(k, j, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidxz(k, j, i));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri0 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid0, edges);
                        m_linMesh->AddGeom(newtriid0, std::move(tri0));
                        E_fidyzd1(k, j, i) = newtriid0;

                        /* Tri1  temp Tri - yz diagonal
                         */
                        int newtriid1 =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidyz(k, j, i));
                        edges[1] = m_linMesh->GetSegGeom(edgid);
                        edges[2] = m_linMesh->GetSegGeom(E_eidx(k + 1, j, i));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri1 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid1, edges);
                        m_linMesh->AddGeom(newtriid1, std::move(tri1));

                        /* Tri2  temp Tri yz diagonal */
                        int newtriid2 =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidx(k, j + 1, i));
                        edges[1] = m_linMesh->GetSegGeom(edgid);
                        edges[2] = m_linMesh->GetSegGeom(E_eidyz(k, j, i + 1));
                        // shuffle so lowest vid is on edges 1, 2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri2 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid2, edges);
                        m_linMesh->AddGeom(newtriid2, std::move(tri2));

                        /* Tri3  temp Tri xz diagonal */
                        int newtriid3 =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidxy(k, j, i));
                        edges[1] = m_linMesh->GetSegGeom(edgid);
                        edges[2] = m_linMesh->GetSegGeom(E_eidz(k, j, i + 1));
                        // shuffle so lowest vid is on edges 1,  2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri3 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid3, edges);
                        m_linMesh->AddGeom(newtriid3, std::move(tri3));

                        /* Tri4  temp Tri xz diagonal */
                        int newtriid4 =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(edgid);
                        edges[1] = m_linMesh->GetSegGeom(E_eidxy(k + 1, j, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidz(k, j + 1, i));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri4 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid4, edges);
                        m_linMesh->AddGeom(newtriid4, std::move(tri4));

                        //     /* Tri5  Tri xz vertical */
                        int newtriid5 =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidx(k, j + 1, i));
                        edges[1] = m_linMesh->GetSegGeom(E_eidxz(k, j + 1, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidz(k, j + 1, i));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri5 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid5, edges);
                        m_linMesh->AddGeom(newtriid5, std::move(tri5));
                        E_fidxz1(k, j + 1, i) = newtriid5;

                        /* Tri6  Tri yz vertical */
                        int newtriid6 =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidy(k, j, i + 1));
                        edges[1] = m_linMesh->GetSegGeom(E_eidyz(k, j, i + 1));
                        edges[2] = m_linMesh->GetSegGeom(E_eidz(k, j, i + 1));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri6 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid6, edges);
                        m_linMesh->AddGeom(newtriid6, std::move(tri6));
                        E_fidyz1(k, j, i + 1) = newtriid6;

                        //     /* Tri7  Tri xy horizontal */
                        int newtriid7 =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidx(k + 1, j, i));
                        edges[1] = m_linMesh->GetSegGeom(E_eidy(k + 1, j, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidxy(k + 1, j, i));
                        // shuffle so lowest vid is on edges 1,2
                        TriShuffleEdges(edges);
                        TriGeomUniquePtr tri7 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid7, edges);
                        m_linMesh->AddGeom(newtriid7, std::move(tri7));
                        E_fidxy1(k + 1, j, i) = newtriid7;

                        /* add in four new Tets */
                        // Tet 0
                        int newtetid =
                            elmtid * nsplit * nsplit * nsplit + cnt++;
                        std::array<TriGeom *, 4> faces;
                        faces[0] = m_linMesh->GetTriGeom(E_fidyzd1(k, j, i));
                        faces[1] = m_linMesh->GetTriGeom(newtriid3);
                        faces[2] = m_linMesh->GetTriGeom(newtriid1);
                        faces[3] = m_linMesh->GetTriGeom(E_fidxz2(k, j, i));
                        // store location of vertid
                        TetShuffleFaces(faces); // shuffle so maxvert id at top
                        TetGeomUniquePtr tet0 =
                            ObjPoolManager<SpatialDomains::TetGeom>::
                                AllocateUniquePtr(newtetid, faces);
                        m_linMesh->AddGeom(newtetid, std::move(tet0));
                        m_linMesh->PopulateFaceToElMap(
                            m_linMesh->GetTetGeom(newtetid), 4);

                        // Tet 1
                        newtetid = elmtid * nsplit * nsplit * nsplit + cnt++;
                        faces[0] = m_linMesh->GetTriGeom(newtriid1);
                        faces[1] = m_linMesh->GetTriGeom(E_fidyz2(k, j, i));
                        faces[2] = m_linMesh->GetTriGeom(E_fidxy1(k + 1, j, i));
                        faces[3] = m_linMesh->GetTriGeom(newtriid4);
                        TetShuffleFaces(faces);
                        // shuffle so max vert id at top
                        TetGeomUniquePtr tet1 =
                            ObjPoolManager<SpatialDomains::TetGeom>::
                                AllocateUniquePtr(newtetid, faces);
                        m_linMesh->AddGeom(newtetid, std::move(tet1));
                        m_linMesh->PopulateFaceToElMap(
                            m_linMesh->GetTetGeom(newtetid), 4);

                        // Tet 2
                        newtetid = elmtid * nsplit * nsplit * nsplit + cnt++;
                        faces[0] = m_linMesh->GetTriGeom(E_fidxy2(k, j, i));
                        faces[1] = m_linMesh->GetTriGeom(newtriid3);
                        faces[2] = m_linMesh->GetTriGeom(E_fidyz1(k, j, i + 1));
                        faces[3] = m_linMesh->GetTriGeom(newtriid2);
                        TetShuffleFaces(faces);
                        // shuffle so max vert id at top
                        TetGeomUniquePtr tet2 =
                            ObjPoolManager<SpatialDomains::TetGeom>::
                                AllocateUniquePtr(newtetid, faces);
                        m_linMesh->AddGeom(newtetid, std::move(tet2));
                        m_linMesh->PopulateFaceToElMap(
                            m_linMesh->GetTetGeom(newtetid), 4);

                        // Tet 3
                        newtetid = elmtid * nsplit * nsplit * nsplit + cnt++;
                        faces[0] = m_linMesh->GetTriGeom(newtriid2);
                        faces[1] = m_linMesh->GetTriGeom(E_fidxz1(k, j + 1, i));
                        faces[2] = m_linMesh->GetTriGeom(newtriid4);
                        faces[3] = m_linMesh->GetTriGeom(E_fidyzd2(k, j, i));
                        TetShuffleFaces(faces); // shuffle so max
                        // vert id at top
                        TetGeomUniquePtr tet3 =
                            ObjPoolManager<SpatialDomains::TetGeom>::
                                AllocateUniquePtr(newtetid, faces);
                        m_linMesh->AddGeom(newtetid, std::move(tet3));
                        m_linMesh->PopulateFaceToElMap(
                            m_linMesh->GetTetGeom(newtetid), 4);
                    }

                    /* bottom Corner Tet */
                    int newtetid = elmtid * nsplit * nsplit * nsplit + cnt++;

                    std::array<TriGeom *, 4> faces;
                    faces[0] = m_linMesh->GetTriGeom(E_fidxy1(k, j, i));
                    faces[1] = m_linMesh->GetTriGeom(E_fidxz1(k, j, i));
                    faces[2] = m_linMesh->GetTriGeom(E_fidyzd1(k, j, i));
                    faces[3] = m_linMesh->GetTriGeom(E_fidyz1(k, j, i));
                    TetShuffleFaces(faces); // shuffle
                    // so max vert id at top
                    TetGeomUniquePtr tet = ObjPoolManager<
                        SpatialDomains::TetGeom>::AllocateUniquePtr(newtetid,
                                                                    faces);
                    m_linMesh->AddGeom(newtetid, std::move(tet));
                    m_linMesh->PopulateFaceToElMap(
                        m_linMesh->GetTetGeom(newtetid), 4);
                }
            }
        }
        ASSERTL1(ecnt <= edgeoffset, "This should not happen, "
                                     "edgeoffset might need to be increased");
        ASSERTL1(fcnt <= faceoffset, "This should not happen, "
                                     "faceoffset might need to be increased");

        // collect a mapping of Vid to point vertex location
        cnt = 0;
        for (int k = 0; k < nsplit + 1; ++k)
        {
            for (int j = 0; j < nsplit + 1 - k; ++j)
            {
                for (int i = 0; i < nsplit + 1 - j - k; ++i)
                {
                    CoeffMap[elmtid][E_vid(k, j, i)] = cnt++;
                }
            }
        }
    }
}

void LinearMeshGraph::LinMeshSetUpPrismGeom(
    int nsplit, std::map<int, int> &FceEdgOffset,
    std::map<int, std::map<int, int>> &CoeffMap, bool UseGLL)
{
    // get vertex ids on faces of 3D shape
    Array3D<int> E_vid(nsplit + 1, nsplit + 1, nsplit + 1);

    // X-aligned edge
    Array3D<int> E_eidx(nsplit + 1, nsplit + 1, nsplit + 1);
    // Y-aligned edge
    Array3D<int> E_eidy(nsplit + 1, nsplit + 1, nsplit + 1);
    // Z-aligned edge
    Array3D<int> E_eidz(nsplit + 1, nsplit + 1, nsplit + 1);

    // XZ-diagonal edge
    Array3D<int> E_eidxz(nsplit + 1, nsplit + 1, nsplit + 1);

    // xy face ids
    Array3D<int> E_fidxy(nsplit + 1, nsplit + 1, nsplit + 1);
    // xz face ids (lower triangular if simplex)
    Array3D<int> E_fidxz1(nsplit + 1, nsplit + 1, nsplit + 1);
    //  xz face ids (upper triangular if simplex)
    Array3D<int> E_fidxz2(nsplit + 1, nsplit + 1, nsplit + 1);
    //  yz face ids
    Array3D<int> E_fidyz(nsplit + 1, nsplit + 1, nsplit + 1);

    // Fill with yz diagonal face ids
    Array3D<int> E_fidyzd(nsplit + 1, nsplit + 1, nsplit + 1);

    int maxvertid, maxedgeid;
    LibUtilities::CommSharedPtr comm = m_session->GetComm();

    Array<OneD, NekDouble> Xpts, Ypts, Zpts;
    Xpts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1) * (nsplit + 1));
    Ypts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1) * (nsplit + 1));
    Zpts = Array<OneD, NekDouble>((nsplit + 1) * (nsplit + 1) * (nsplit + 1));

    /* reset maxedgeid and maxvertid  */
    maxvertid = m_linMesh->GetGeomMap<PointGeom>().rbegin()->first + 1;
    comm->AllReduce(maxvertid, LibUtilities::ReduceMax);
    maxedgeid = m_linMesh->GetGeomMap<SegGeom>().rbegin()->first + 1;
    comm->AllReduce(maxedgeid, LibUtilities::ReduceMax);

    int maxfaceid = 0;
    if (m_linMesh->GetGeomMap<TriGeom>().size())
    {
        maxfaceid = std::max(
            maxfaceid, m_linMesh->GetGeomMap<TriGeom>().rbegin()->first + 1);
    }
    if (m_linMesh->GetGeomMap<QuadGeom>().size())
    {
        maxfaceid = std::max(
            maxfaceid, m_linMesh->GetGeomMap<QuadGeom>().rbegin()->first + 1);
    }

    comm->AllReduce(maxfaceid, LibUtilities::ReduceMax);

    // define verte reorientaiton - might need for all square faces
    std::map<StdRegions::Orientation, std::array<int, 4>> vertMap{
        {StdRegions::eDir1FwdDir1_Dir2FwdDir2, {{0, 1, 2, 3}}},
        {StdRegions::eDir1BwdDir1_Dir2FwdDir2, {{1, 0, 3, 2}}},
        {StdRegions::eDir1FwdDir1_Dir2BwdDir2, {{3, 2, 1, 0}}},
        {StdRegions::eDir1BwdDir1_Dir2BwdDir2, {{2, 3, 0, 1}}},
        {StdRegions::eDir1FwdDir2_Dir2FwdDir1, {{0, 3, 2, 1}}},
        {StdRegions::eDir1BwdDir2_Dir2FwdDir1, {{3, 0, 1, 2}}},
        {StdRegions::eDir1FwdDir2_Dir2BwdDir1, {{1, 2, 3, 0}}},
        {StdRegions::eDir1BwdDir2_Dir2BwdDir1, {{2, 1, 0, 3}}}};

    // define edge reorientaiton - might need for all square faces
    std::map<StdRegions::Orientation, std::array<int, 4>> edgMap{
        {StdRegions::eDir1FwdDir1_Dir2FwdDir2, {{0, 1, 2, 3}}},
        {StdRegions::eDir1BwdDir1_Dir2FwdDir2, {{0, 3, 2, 1}}},
        {StdRegions::eDir1FwdDir1_Dir2BwdDir2, {{2, 1, 0, 3}}},
        {StdRegions::eDir1BwdDir1_Dir2BwdDir2, {{2, 3, 0, 1}}},
        {StdRegions::eDir1FwdDir2_Dir2FwdDir1, {{3, 2, 1, 0}}},
        {StdRegions::eDir1BwdDir2_Dir2FwdDir1, {{3, 0, 1, 2}}},
        {StdRegions::eDir1FwdDir2_Dir2BwdDir1, {{1, 2, 3, 0}}},
        {StdRegions::eDir1BwdDir2_Dir2BwdDir1, {{1, 0, 3, 2}}}};

    //--------------------------------------------------------
    // Split prism elements of macro mesh into new edges and vertices
    //--------------------------------------------------------

    int edgeoffset = 3 * nsplit * (nsplit - 1) * (nsplit - 1);
    int faceoffset = 3 * nsplit * nsplit * (nsplit - 1);

    // For every original elmt keep a map between the vid of the split
    // tet and the location of the points in the macro nodes to create
    // the LinCoeffMap when required. Needed since we reshuffle tri faces
    // orientaiton so cannot infer pattern directly
    for (auto [id, geom] : m_graph->GetGeomMap<PrismGeom>())
    {
        int vcnt   = 0;
        int ecnt   = 0;
        int fcnt   = 0;
        int elmtid = geom->GetGlobalID();

        StdRegions::StdExpansionSharedPtr Xmap = geom->GetXmap();

        /* face 0 - Quad */
        int fid        = geom->GetFace(0)->GetGlobalID();
        int fce_offset = 2 * fid * nsplit * nsplit;

        StdRegions::Orientation orient = geom->GetForient(0);
        Array<OneD, int> idmap;
        // macro split reorientation
        Xmap->ReOrientTracePhysMap(orient, idmap, nsplit, nsplit, false);

        int cnt = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit; ++i)
            {
                int Fce_id = fce_offset + idmap[i + j * nsplit];

                E_fidxy(0, j, i) = Fce_id;

                QuadGeom *face = m_linMesh->GetQuadGeom(Fce_id);

                // fill in vertices & edges from faces1
                E_vid(0, j, i)         = face->GetVid(vertMap[orient][0]);
                E_vid(0, j, i + 1)     = face->GetVid(vertMap[orient][1]);
                E_vid(0, j + 1, i)     = face->GetVid(vertMap[orient][3]);
                E_vid(0, j + 1, i + 1) = face->GetVid(vertMap[orient][2]);

                E_eidx(0, j, i)     = face->GetEid(edgMap[orient][0]);
                E_eidx(0, j + 1, i) = face->GetEid(edgMap[orient][2]);
                E_eidy(0, j, i)     = face->GetEid(edgMap[orient][3]);
                E_eidy(0, j, i + 1) = face->GetEid(edgMap[orient][1]);
            }
        }

        /* face 1 - Tri */
        fid        = geom->GetFace(1)->GetGlobalID();
        fce_offset = 2 * fid * nsplit * nsplit;
        bool FceFwd =
            geom->GetForient(1) == StdRegions::eDir1FwdDir1_Dir2FwdDir2;
        cnt = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit - j; ++i)
            {
                int Fce_id = FceFwd ? fce_offset + cnt + i
                                    : fce_offset + cnt + nsplit - j - 1 - i;

                E_fidxz1(j, 0, i) = Fce_id;

                // due to edge shuffling earlier need to identify
                // know
                // orientation (singlar vertex at top) orientation to
                // fill in stored edges and vertices
                int oset = FceEdgOffset[Fce_id];

                TriGeom *face = m_linMesh->GetTriGeom(Fce_id);
                // fill in vertices & edges from faces
                if (FceFwd)
                {
                    E_vid(j, 0, i)     = face->GetVid(oset);
                    E_vid(j, 0, i + 1) = face->GetVid((oset + 1) % 3);

                    E_eidz(j, 0, i)  = face->GetEid((oset + 2) % 3);
                    E_eidxz(j, 0, i) = face->GetEid((oset + 1) % 3);
                }
                else
                {
                    E_vid(j, 0, i)     = face->GetVid((oset + 1) % 3);
                    E_vid(j, 0, i + 1) = face->GetVid(oset);

                    E_eidz(j, 0, i)  = face->GetEid((oset + 1) % 3);
                    E_eidxz(j, 0, i) = face->GetEid((oset + 2) % 3);
                }
                E_vid(j + 1, 0, i) = face->GetVid((oset + 2) % 3);
                E_eidx(j, 0, i)    = face->GetEid(oset);
            }
            cnt += nsplit - j;

            for (int i = 0; i < nsplit - j - 1; ++i)
            {
                int Fce_id        = FceFwd ? fce_offset + cnt + i
                                           : fce_offset + cnt + nsplit - j - 2 - i;
                E_fidxz2(j, 0, i) = Fce_id;
            }
            cnt += nsplit - j - 1;
        }

        /* face 2 - Quad on diagonal */
        fid        = geom->GetFace(2)->GetGlobalID();
        fce_offset = 2 * fid * nsplit * nsplit;
        orient     = geom->GetForient(2);

        // macro split reorientation
        Xmap->ReOrientTracePhysMap(orient, idmap, nsplit, nsplit, false);

        cnt = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit; ++i)
            {
                int Fce_id = fce_offset + idmap[i + j * nsplit];

                E_fidyzd(j, i, nsplit - 1 - j) = Fce_id;

                QuadGeom *face = m_linMesh->GetQuadGeom(Fce_id);

                // fill in vertices & edges from faces1
                E_vid(j, i, nsplit - j)     = face->GetVid(vertMap[orient][0]);
                E_vid(j, i + 1, nsplit - j) = face->GetVid(vertMap[orient][1]);
                E_vid(j + 1, i, nsplit - 1 - j) =
                    face->GetVid(vertMap[orient][3]);
                E_vid(j + 1, i + 1, nsplit - 1 - j) =
                    face->GetVid(vertMap[orient][2]);

                E_eidy(j, i, nsplit - j) = face->GetEid(edgMap[orient][0]);
                E_eidy(j + 1, i, nsplit - 1 - j) =
                    face->GetEid(edgMap[orient][2]);
                E_eidxz(j, i, nsplit - 1 - j) = face->GetEid(edgMap[orient][3]);
                E_eidxz(j, i + 1, nsplit - 1 - j) =
                    face->GetEid(edgMap[orient][1]);
            }
        }
        /* face 3 - Tri */
        fid        = geom->GetFace(3)->GetGlobalID();
        fce_offset = 2 * fid * nsplit * nsplit;
        FceFwd = geom->GetForient(3) == StdRegions::eDir1FwdDir1_Dir2FwdDir2;
        cnt    = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit - j; ++i)
            {
                int Fce_id = FceFwd ? fce_offset + cnt + i
                                    : fce_offset + cnt + nsplit - j - 1 - i;

                E_fidxz1(j, nsplit, i) = Fce_id;

                // due to edge shuffling earlier need to identify
                // know
                // orientation (singlar vertex at top) orientation to
                // fill in stored edges and vertices
                int oset = FceEdgOffset[Fce_id];

                TriGeom *face = m_linMesh->GetTriGeom(Fce_id);
                // fill in vertices & edges from faces
                if (FceFwd)
                {
                    E_vid(j, nsplit, i)     = face->GetVid(oset);
                    E_vid(j, nsplit, i + 1) = face->GetVid((oset + 1) % 3);

                    E_eidz(j, nsplit, i)  = face->GetEid((oset + 2) % 3);
                    E_eidxz(j, nsplit, i) = face->GetEid((oset + 1) % 3);
                }
                else
                {
                    E_vid(j, nsplit, i)     = face->GetVid((oset + 1) % 3);
                    E_vid(j, nsplit, i + 1) = face->GetVid(oset);

                    E_eidz(j, nsplit, i)  = face->GetEid((oset + 1) % 3);
                    E_eidxz(j, nsplit, i) = face->GetEid((oset + 2) % 3);
                }
                E_vid(j + 1, nsplit, i) = face->GetVid((oset + 2) % 3);
                E_eidx(j, nsplit, i)    = face->GetEid(oset);
            }
            cnt += nsplit - j;

            for (int i = 0; i < nsplit - j - 1; ++i)
            {
                int Fce_id             = FceFwd ? fce_offset + cnt + i
                                                : fce_offset + cnt + nsplit - j - 2 - i;
                E_fidxz2(j, nsplit, i) = Fce_id;
            }
            cnt += nsplit - j - 1;
        }

        /* face 4 - Quad */
        fid        = geom->GetFace(4)->GetGlobalID();
        fce_offset = 2 * fid * nsplit * nsplit;

        orient = geom->GetForient(4);
        // macro split reorientation
        Xmap->ReOrientTracePhysMap(orient, idmap, nsplit, nsplit, false);

        cnt = 0;
        for (int j = 0; j < nsplit; ++j)
        {
            for (int i = 0; i < nsplit; ++i)
            {
                int Fce_id = fce_offset + idmap[i + j * nsplit];

                E_fidyz(j, i, 0) = Fce_id;

                QuadGeom *face = m_linMesh->GetQuadGeom(Fce_id);

                // fill in vertices & edges from faces1
                E_vid(j, i, 0)         = face->GetVid(vertMap[orient][0]);
                E_vid(j, i + 1, 0)     = face->GetVid(vertMap[orient][1]);
                E_vid(j + 1, i, 0)     = face->GetVid(vertMap[orient][3]);
                E_vid(j + 1, i + 1, 0) = face->GetVid(vertMap[orient][2]);

                E_eidy(j, i, 0)     = face->GetEid(edgMap[orient][0]);
                E_eidy(j + 1, i, 0) = face->GetEid(edgMap[orient][2]);
                E_eidz(j, i, 0)     = face->GetEid(edgMap[orient][3]);
                E_eidz(j, i + 1, 0) = face->GetEid(edgMap[orient][1]);
            }
        }

        Array<OneD, NekDouble> tmp(Xmap->GetTotPoints());

        // fill in interior vids if exist
        if (UseGLL) // get GLL points
        {
            Xmap->BwdTrans(geom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToGLL(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToGLL(tmp, Ypts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(2), tmp);
            Xmap->PhysInterpToGLL(tmp, Zpts, nsplit + 1);
        }
        else /* get equispaced points  */
        {
            Xmap->BwdTrans(geom->GetCoeffs(0), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Xpts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(1), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Ypts, nsplit + 1);
            Xmap->BwdTrans(geom->GetCoeffs(2), tmp);
            Xmap->PhysInterpToSimplexEquiSpaced(tmp, Zpts, nsplit + 1);
        }

        cnt = (nsplit + 1) * (nsplit + 1);
        for (int k = 1; k < nsplit; ++k)
        {
            cnt += nsplit + 1 - k; // first row
            for (int j = 1; j < nsplit; ++j)
            {
                cnt++;
                for (int i = 1; i < nsplit - k; ++i, cnt++)
                {
                    int vid =
                        maxvertid +
                        elmtid * (nsplit - 1) * (nsplit - 1) * (nsplit - 1) +
                        vcnt++;

                    E_vid(k, j, i) = vid;

                    PointGeomUniquePtr vert =
                        ObjPoolManager<PointGeom>::AllocateUniquePtr(
                            m_spaceDimension, vid, Xpts[cnt], Ypts[cnt],
                            Zpts[cnt]);

                    m_linMesh->AddGeom(vid, std::move(vert));
                }
                cnt++;
            }
            cnt += nsplit + 1 - k; // last row
        }

        // declare Prisms
        cnt = 0;
        for (int k = 0; k < nsplit; ++k)
        {
            for (int j = 0; j < nsplit; ++j)
            {
                for (int i = 0; i < nsplit - k; ++i)
                {

                    if ((i < nsplit - 1 - k) && (j < nsplit - 1))
                    {
                        /* add three new edges */

                        /* new diagonals */
                        int edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        std::array<PointGeom *, 2> vert;

                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k, j + 1, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i));
                        SegGeomUniquePtr edge0 =
                            ObjPoolManager<SpatialDomains::SegGeom>::
                                AllocateUniquePtr(edgid, m_spaceDimension,
                                                  vert);
                        m_linMesh->AddGeom(edgid, std::move(edge0));
                        E_eidxz(k, j + 1, i) = edgid;

                        // add in vertical edge
                        edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k, j + 1, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i + 1));
                        SegGeomUniquePtr edge1 =
                            ObjPoolManager<SpatialDomains::SegGeom>::
                                AllocateUniquePtr(edgid, m_spaceDimension,
                                                  vert);
                        m_linMesh->AddGeom(edgid, std::move(edge1));
                        E_eidz(k, j + 1, i + 1) = edgid;

                        // horizontal edge
                        edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i));
                        SegGeomUniquePtr edge2 =
                            ObjPoolManager<SpatialDomains::SegGeom>::
                                AllocateUniquePtr(edgid, m_spaceDimension,
                                                  vert);
                        m_linMesh->AddGeom(edgid, std::move(edge2));
                        E_eidx(k + 1, j + 1, i) = edgid;

                        /* Add in new faces */
                        /* Far end Tri */
                        int newtriid = elmtid * faceoffset + maxfaceid + fcnt++;
                        std::array<SegGeom *, 3>
                            edges; // BUG KK - it was 4 originally ???
                        std::array<SegGeom *, 3> edgessort;

                        edgessort[0] =
                            m_linMesh->GetSegGeom(E_eidx(k + 1, 0, i));
                        edgessort[1] = m_linMesh->GetSegGeom(E_eidxz(k, 0, i));
                        edgessort[2] =
                            m_linMesh->GetSegGeom(E_eidz(k, 0, i + 1));
                        edges[0] =
                            m_linMesh->GetSegGeom(E_eidx(k + 1, j + 1, i));
                        edges[1] = m_linMesh->GetSegGeom(E_eidxz(k, j + 1, i));
                        edges[2] =
                            m_linMesh->GetSegGeom(E_eidz(k, j + 1, i + 1));
                        // shuffle so lowest vid is
                        // on edges 1, 2
                        TriShuffleEdges(edges, &edgessort); // BUG KK ?
                        TriGeomUniquePtr tri1 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid, edges);
                        m_linMesh->AddGeom(newtriid, std::move(tri1));

                        E_fidxz2(k, j + 1, i) = newtriid;
                    }

                    if (i < (nsplit - 2 - k))
                    {
                        /* new y-parallel edge  */
                        int edgid = elmtid * edgeoffset + maxedgeid + ecnt++;
                        std::array<PointGeom *, 2> vert;

                        vert[0] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j, i + 1));
                        vert[1] =
                            m_linMesh->GetPointGeom(E_vid(k + 1, j + 1, i + 1));
                        SegGeomUniquePtr edge0 =
                            ObjPoolManager<SpatialDomains::SegGeom>::
                                AllocateUniquePtr(edgid, m_spaceDimension,
                                                  vert);
                        m_linMesh->AddGeom(edgid, std::move(edge0));
                        E_eidy(k + 1, j, i + 1) = edgid;
                    }

                    if (i < (nsplit - 1 - k))
                    {
                        std::array<SegGeom *, 4> edges;
                        /* diagonal face */
                        int newquadid =
                            elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] = m_linMesh->GetSegGeom(E_eidy(k, j, i + 1));
                        edges[1] = m_linMesh->GetSegGeom(E_eidxz(k, j + 1, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidy(k + 1, j, i));
                        edges[3] = m_linMesh->GetSegGeom(E_eidxz(k, j, i));

                        QuadGeomUniquePtr quad0 =
                            ObjPoolManager<SpatialDomains::QuadGeom>::
                                AllocateUniquePtr(newquadid, edges);
                        m_linMesh->AddGeom(newquadid, std::move(quad0));
                        E_fidyzd(k, j, i) = newquadid;

                        /* vertical face */
                        newquadid = elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0]  = m_linMesh->GetSegGeom(E_eidy(k, j, i + 1));
                        edges[1] =
                            m_linMesh->GetSegGeom(E_eidz(k, j + 1, i + 1));
                        edges[2] =
                            m_linMesh->GetSegGeom(E_eidy(k + 1, j, i + 1));
                        edges[3] = m_linMesh->GetSegGeom(E_eidz(k, j, i + 1));

                        QuadGeomUniquePtr quad1 =
                            ObjPoolManager<SpatialDomains::QuadGeom>::
                                AllocateUniquePtr(newquadid, edges);
                        m_linMesh->AddGeom(newquadid, std::move(quad1));
                        E_fidyz(k, j, i + 1) = newquadid;

                        /* horizontal face */
                        newquadid = elmtid * faceoffset + maxfaceid + fcnt++;
                        edges[0] =
                            m_linMesh->GetSegGeom(E_eidy(k + 1, j, i + 1));
                        edges[1] =
                            m_linMesh->GetSegGeom(E_eidx(k + 1, j + 1, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidy(k + 1, j, i));
                        edges[3] = m_linMesh->GetSegGeom(E_eidx(k + 1, j, i));

                        QuadGeomUniquePtr quad2 =
                            ObjPoolManager<SpatialDomains::QuadGeom>::
                                AllocateUniquePtr(newquadid, edges);
                        m_linMesh->AddGeom(newquadid, std::move(quad2));
                        E_fidxy(k + 1, j, i) = newquadid;

                        // /* top Corner Prism */
                        int newprismid =
                            elmtid * nsplit * nsplit * nsplit + cnt++;

                        std::array<Geometry2D *, 5> faces;
                        faces[0] = m_linMesh->GetQuadGeom(E_fidxy(k + 1, j, i));
                        faces[1] = m_linMesh->GetTriGeom(E_fidxz2(k, j, i));
                        faces[2] = m_linMesh->GetQuadGeom(E_fidyzd(k, j, i));
                        faces[3] = m_linMesh->GetTriGeom(E_fidxz2(k, j + 1, i));
                        faces[4] = m_linMesh->GetQuadGeom(E_fidyz(k, j, i + 1));
                        PrismShuffleFaces(faces); // ensure faces are
                        // ordered correctly
                        PrismGeomUniquePtr prism =
                            ObjPoolManager<SpatialDomains::PrismGeom>::
                                AllocateUniquePtr(newprismid, faces);
                        m_linMesh->AddGeom(newprismid, std::move(prism));
                        m_linMesh->PopulateFaceToElMap(
                            m_linMesh->GetPrismGeom(newprismid), 5);
                    }

                    if (j < nsplit - 1)
                    {
                        int newtriid = elmtid * faceoffset + maxfaceid + fcnt++;
                        std::array<SegGeom *, 3> edges;
                        std::array<SegGeom *, 3> edgessort;

                        edgessort[0] = m_linMesh->GetSegGeom(E_eidx(k, 0, i));
                        edgessort[1] = m_linMesh->GetSegGeom(E_eidxz(k, 0, i));
                        edgessort[2] = m_linMesh->GetSegGeom(E_eidz(k, 0, i));
                        edges[0] = m_linMesh->GetSegGeom(E_eidx(k, j + 1, i));
                        edges[1] = m_linMesh->GetSegGeom(E_eidxz(k, j + 1, i));
                        edges[2] = m_linMesh->GetSegGeom(E_eidz(k, j + 1, i));
                        // shuffle so lowest vid is on
                        // edges 1, 2
                        TriShuffleEdges(edges, &edgessort);
                        TriGeomUniquePtr tri0 =
                            ObjPoolManager<SpatialDomains::TriGeom>::
                                AllocateUniquePtr(newtriid, edges);
                        m_linMesh->AddGeom(newtriid, std::move(tri0));
                        E_fidxz1(k, j + 1, i) = newtriid;
                    }

                    /* bottom Corner Prism */
                    int newprismid = elmtid * nsplit * nsplit * nsplit + cnt++;

                    std::array<Geometry2D *, 5> faces;
                    faces[0] = m_linMesh->GetQuadGeom(E_fidxy(k, j, i));
                    faces[1] = m_linMesh->GetTriGeom(E_fidxz1(k, j, i));
                    faces[2] = m_linMesh->GetQuadGeom(E_fidyzd(k, j, i));
                    faces[3] = m_linMesh->GetTriGeom(E_fidxz1(k, j + 1, i));
                    faces[4] = m_linMesh->GetQuadGeom(E_fidyz(k, j, i));
                    PrismShuffleFaces(faces); // ensure faces are
                    //                         ordered correctly
                    PrismGeomUniquePtr prism =
                        ObjPoolManager<SpatialDomains::PrismGeom>::
                            AllocateUniquePtr(newprismid, faces);
                    m_linMesh->AddGeom(newprismid, std::move(prism));

                    m_linMesh->PopulateFaceToElMap(
                        m_linMesh->GetPrismGeom(newprismid), 5);
                }
            }
        }
        ASSERTL1(
            ecnt <= edgeoffset,
            "This should not happen, edgeoffset might need to be increased");
        ASSERTL1(
            fcnt <= faceoffset,
            "This should not happen, faceoffset might need to be increased");

        // collect a mapping of Vid to point vertex location
        cnt = 0;
        for (int k = 0; k < nsplit + 1; ++k)
        {
            for (int j = 0; j < nsplit + 1; ++j)
            {
                for (int i = 0; i < nsplit + 1 - k; ++i)
                {
                    CoeffMap[elmtid][E_vid(k, j, i)] = cnt++;
                }
            }
        }
    }
}

void LinearMeshGraph::LinMeshSetUpCompositesDomain(
    int nsplit, std::map<int, std::pair<int, std::vector<int>>> &LinCoeffMap,
    std::map<int, std::map<int, int>> &CoeffMap2D,
    std::map<int, std::map<int, int>> &CoeffMap3D, bool useSimplex)
{
    auto &linMeshComp = m_linMesh->GetComposites();

    for (auto &c : m_graph->GetComposites())
    {
        CompositeSharedPtr curVector =
            MemoryManager<Composite>::AllocateSharedPtr();
        linMeshComp[c.first] = curVector;

        for (auto &g : c.second->m_geomVec)
        {
            switch (g->GetShapeType())
            {
                case LibUtilities::eSegment:
                {
                    int id = g->GetGlobalID();
                    for (int i = 0; i < nsplit; ++i)
                    {
                        int newid = id * nsplit + i;

                        auto it = m_linMesh->GetGeomMap<SegGeom>().find(newid);
                        ASSERTL1(it != m_linMesh->GetGeomMap<SegGeom>().end(),
                                 "Failed to find LinMesh edge id");
                        linMeshComp[c.first]->m_geomVec.push_back(it->second);
                    }
                }
                break;
                case LibUtilities::eTriangle:
                {
                    int id = g->GetGlobalID();
                    for (int i = 0; i < nsplit * nsplit; ++i)
                    {
                        int newid = 2 * id * nsplit * nsplit + i;

                        auto it = m_linMesh->GetGeomMap<TriGeom>().find(newid);
                        ASSERTL1(it != m_linMesh->GetGeomMap<TriGeom>().end(),
                                 "Failed to find LinMesh tri id");
                        linMeshComp[c.first]->m_geomVec.push_back(it->second);
                    }
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    int id = g->GetGlobalID();

                    if (useSimplex)
                    {
                        for (int i = 0; i < nsplit * nsplit; ++i)
                        {
                            int newid = 2 * (id * nsplit * nsplit + i);

                            auto it1 =
                                m_linMesh->GetGeomMap<TriGeom>().find(newid);
                            ASSERTL1(
                                it1 != m_linMesh->GetGeomMap<TriGeom>().end(),
                                "Failed to find LinMesh tri id (from quad)");
                            linMeshComp[c.first]->m_geomVec.push_back(
                                it1->second);

                            auto it2 = m_linMesh->GetGeomMap<TriGeom>().find(
                                newid + 1);
                            ASSERTL1(
                                it2 != m_linMesh->GetGeomMap<TriGeom>().end(),
                                "Failed to find LinMesh tri id (from quad)");
                            linMeshComp[c.first]->m_geomVec.push_back(
                                it2->second);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < nsplit * nsplit; ++i)
                        {
                            int newid = 2 * id * nsplit * nsplit + i;

                            auto it =
                                m_linMesh->GetGeomMap<QuadGeom>().find(newid);
                            ASSERTL1(
                                it != m_linMesh->GetGeomMap<QuadGeom>().end(),
                                "Failed to find LinMesh quad id");
                            linMeshComp[c.first]->m_geomVec.push_back(
                                it->second);
                        }
                    }
                }
                break;
                case LibUtilities::eTetrahedron:
                {
                    int id = g->GetGlobalID();
                    for (int i = 0; i < nsplit * nsplit * nsplit; ++i)
                    {
                        int newid = id * nsplit * nsplit * nsplit + i;

                        auto it = m_linMesh->GetGeomMap<TetGeom>().find(newid);
                        ASSERTL1(it != m_linMesh->GetGeomMap<TetGeom>().end(),
                                 "Failed to find LinMesh  tet id");
                        linMeshComp[c.first]->m_geomVec.push_back(it->second);
                    }
                }
                break;
                case LibUtilities::ePrism:
                {
                    int id = g->GetGlobalID();
                    for (int i = 0; i < nsplit * nsplit * nsplit; ++i)
                    {
                        int newid = id * nsplit * nsplit * nsplit + i;

                        auto it =
                            m_linMesh->GetGeomMap<PrismGeom>().find(newid);
                        ASSERTL1(it != m_linMesh->GetGeomMap<PrismGeom>().end(),
                                 "Failed to find LinMesh prism id");
                        linMeshComp[c.first]->m_geomVec.push_back(it->second);
                    }
                }
                break;
                default:
                    NEKERROR(ErrorUtil::ewarning, "Shape  Not set up  ");
                    break;
            }
        }
    }

    // reset domain information based on new composites;
    for (auto &d : m_graph->GetDomain())
    {
        // make a list of composite ids
        std::string compstr;
        for (auto sd = d.second.begin(); sd != d.second.end();)
        {
            compstr += std::to_string(sd->first);
            compstr += (++sd != d.second.end()) ? "," : " ";
        }
        if (compstr.size()) // put in this check since in hdf5 format have empty
                            // string cases
        {
            std::map<int, CompositeSharedPtr> unrollDomain;
            m_linMesh->GetCompositeList(compstr, unrollDomain);
            m_linMesh->GetDomain()[d.first] = unrollDomain;
        }
    }

    // Setup Epxansion info as linear expansion
    m_linMesh->ReadExpansionInfo(SetupLinearExpansionType());

    // loop over first expansion entry
    for (auto &expIt : *m_graph->GetExpansionInfoMap().begin()->second)
    {
        Geometry *geom = expIt.second->m_geomPtr;
        // expIt.second->m_geomPtr
        switch (expIt.second->m_geomPtr->GetShapeType())
        {
            case LibUtilities::Tri:
            {
                int elmtid = geom->GetGlobalID();
                ASSERTL1(CoeffMap2D[elmtid].size() ==
                             (nsplit + 1) * (nsplit + 2) / 2,
                         "CoeffMap2D is not the correct length");

                for (int i = 0; i < nsplit * nsplit; ++i)
                {
                    TriGeom *tri =
                        m_linMesh->GetTriGeom(2 * elmtid * nsplit * nsplit + i);

                    std::vector<int> offset(3);

                    offset[0] = CoeffMap2D[elmtid][tri->GetVid(0)];
                    offset[1] = CoeffMap2D[elmtid][tri->GetVid(2)];
                    offset[2] = CoeffMap2D[elmtid][tri->GetVid(1)];

                    LinCoeffMap[2 * elmtid * nsplit * nsplit + i] =
                        std::make_pair(elmtid, offset);
                }
            }
            break;
            case LibUtilities::Quad:
                if (useSimplex)
                {
                    int elmtid = geom->GetGlobalID();
                    for (int i = 0, e = 0; i < nsplit; ++i)
                    {
                        for (int j = 0; j < nsplit; ++j, ++e)
                        {
                            std::vector<int> offset(3);
                            offset[0] = i * (nsplit + 1) + j;
                            offset[1] = (i + 1) * (nsplit + 1) + j;
                            offset[2] = i * (nsplit + 1) + j + 1;

                            LinCoeffMap[2 * (elmtid * nsplit * nsplit + e)] =
                                std::make_pair(elmtid, offset);

                            std::vector<int> offset1(3);
                            offset1[0] = (i + 1) * (nsplit + 1) + j + 1;
                            offset1[1] = i * (nsplit + 1) + j + 1;
                            offset1[2] = (i + 1) * (nsplit + 1) + j;

                            LinCoeffMap[2 * (elmtid * nsplit * nsplit + e) +
                                        1] = std::make_pair(elmtid, offset1);
                        }
                    }
                }
                else
                {
                    int elmtid = geom->GetGlobalID();
                    for (int i = 0, e = 0; i < nsplit; ++i)
                    {
                        for (int j = 0; j < nsplit; ++j, ++e)
                        {
                            std::vector<int> offset(4);
                            offset[0] = i * (nsplit + 1) + j;
                            offset[1] = i * (nsplit + 1) + j + 1;
                            offset[2] = (i + 1) * (nsplit + 1) + j;
                            offset[3] = (i + 1) * (nsplit + 1) + j + 1;

                            LinCoeffMap[2 * elmtid * nsplit * nsplit + e] =
                                std::make_pair(elmtid, offset);
                        }
                    }
                }
                break;
            case LibUtilities::Tet:
            {
                int elmtid = geom->GetGlobalID();
                ASSERTL1(CoeffMap3D[elmtid].size() ==
                             (nsplit + 1) * (nsplit + 2) * (nsplit + 3) / 6,
                         "CoeffMap3D (Tet) is not the correct length");

                for (int i = 0; i < nsplit * nsplit * nsplit; ++i)
                {
                    TetGeom *tet = m_linMesh->GetTetGeom(
                        elmtid * nsplit * nsplit * nsplit + i);

                    std::vector<int> offset(4);

                    offset[0] = CoeffMap3D[elmtid][tet->GetVid(0)];
                    offset[1] = CoeffMap3D[elmtid][tet->GetVid(3)];
                    offset[2] = CoeffMap3D[elmtid][tet->GetVid(2)];
                    offset[3] = CoeffMap3D[elmtid][tet->GetVid(1)];

                    LinCoeffMap[elmtid * nsplit * nsplit * nsplit + i] =
                        std::make_pair(elmtid, offset);
                }
            }
            break;
            case LibUtilities::Prism:
            {
                int elmtid = geom->GetGlobalID();
                ASSERTL1(CoeffMap3D[elmtid].size() ==
                             (nsplit + 1) * (nsplit + 1) * (nsplit + 2) / 2,
                         "CoeffMap3D (Prism) is not the correct length");

                for (int i = 0; i < nsplit * nsplit * nsplit; ++i)
                {
                    PrismGeom *prism = m_linMesh->GetPrismGeom(
                        elmtid * nsplit * nsplit * nsplit + i);

                    std::vector<int> offset(6);

                    offset[0] = CoeffMap3D[elmtid][prism->GetVid(0)];
                    offset[1] = CoeffMap3D[elmtid][prism->GetVid(4)];
                    offset[2] = CoeffMap3D[elmtid][prism->GetVid(3)];
                    offset[3] = CoeffMap3D[elmtid][prism->GetVid(5)];
                    offset[4] = CoeffMap3D[elmtid][prism->GetVid(1)];
                    offset[5] = CoeffMap3D[elmtid][prism->GetVid(2)];

                    LinCoeffMap[elmtid * nsplit * nsplit * nsplit + i] =
                        std::make_pair(elmtid, offset);
                }
            }
            break;
            default:
                NEKERROR(ErrorUtil::efatal, "Need to setup local mapping");
                break;
        }
    }
}

} // namespace Nektar::SpatialDomains

#endif
