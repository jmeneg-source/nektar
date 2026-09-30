////////////////////////////////////////////////////////////////////////////////
//
//  File: InputStar.cpp
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

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <boost/algorithm/string.hpp>

#include "InputStar.h"
#include <NekMesh/MeshElements/Element.h>

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar
{
static char const kDefaultState[] = "default";
namespace NekMesh
{

ModuleKey InputStar::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "ccm"), InputStar::create,
    "Reads mesh from Star CCM (.ccm).");

InputStar::InputStar(MeshSharedPtr m) : InputModule(m)
{
    m_config["writelabelsonly"] = ConfigOption(
        true, "0",
        "Just write out tags from star file for each surface/composite");
}

InputStar::~InputStar()
{
}

/**
 * Tecplot file Polyhedron format contains a list of nodes, a node count per
 * face, the node ids, Element ids that are on the left of each face and Element
 * ids which are on the right of each face. There are then a series of zone of
 * each surface. In the case of a surface the number of nodes is not provided
 * indicating it is a 2D zone.
 *
 * @param pFilename Filename of Tecplot file to read.
 */
void InputStar::Process()
{
    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    m_log(VERBOSE) << "Reading CCM+ file: '" << m_config["infile"].as<string>()
                   << "'" << endl;

    InitCCM();

    SetupElements();

    if (m_config["writelabelsonly"].beenSet)
    {
        return;
    }

    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    PrintSummary();
}

void InputStar::SetupElements(void)
{
    int i;
    string line, tag;
    stringstream s;
    streampos pos;
    int nComposite = 0;

    // Read in Nodes
    ReadNodes(m_mesh->m_node);

    // Get list of faces nodes and adjacents elements.
    std::unordered_map<int, vector<int>> FaceNodes;
    Array<OneD, vector<int>> ElementFaces;

    // Read interior faces and set up first part of Element
    // Faces and FaceNodes
    ReadInternalFaces(FaceNodes, ElementFaces);

    vector<vector<int>> BndElementFaces;
    vector<string> Facelabels;
    ReadBoundaryFaces(BndElementFaces, FaceNodes, ElementFaces, Facelabels);

    if (m_config["writelabelsonly"].beenSet)
    {
        nComposite = 2;
        // write boundary zones/composites
        m_log << "Element labels:" << endl;
        for (i = 0; i < BndElementFaces.size(); ++i)
        {
            m_log << " 2D Zone (composite = " << nComposite
                  << ", label = " << Facelabels[i] << ")" << endl;
            nComposite++;
        }
        return;
    }

    // 3D Zone
    // Reset node ordering so that all prism faces have
    // consistent numbering for singular vertex re-ordering
    ResetNodes(m_mesh->m_node, ElementFaces, FaceNodes);

    // create Prisms first
    int nelements = ElementFaces.size();
    m_log(VERBOSE) << nelements << " Elements" << endl;
    m_log(VERBOSE) << "Generating 3D Zones: " << endl;
    int cnt = 0;
    for (i = 0; i < nelements; ++i)
    {
        Array<OneD, int> Nodes =
            SortFaceNodes(m_mesh->m_node, ElementFaces[i], FaceNodes);
        if (ElementFaces[i].size() == 5 && Nodes.size() == 6) // if a prism
        {
            GenElement3D(m_mesh->m_node, i, ElementFaces[i], FaceNodes,
                         nComposite, true);
            ++cnt;
        }
    }
    m_log(VERBOSE) << "  - # of prisms: " << cnt << endl;

    nComposite++;

    // create Pyras second
    cnt = 0;
    for (i = 0; i < nelements; ++i)
    {
        Array<OneD, int> Nodes =
            SortFaceNodes(m_mesh->m_node, ElementFaces[i], FaceNodes);
        if (ElementFaces[i].size() == 5 && Nodes.size() == 5) // if a pyra
        {
            GenElement3D(m_mesh->m_node, i, ElementFaces[i], FaceNodes,
                         nComposite, true);
            ++cnt;
        }
    }
    m_log(VERBOSE) << "  - # of pyramids: " << cnt << endl;

    nComposite++;

    // create Tets third
    cnt = 0;
    for (i = 0; i < nelements; ++i)
    {
        if (ElementFaces[i].size() == 4) // if a tetra
        {
            GenElement3D(m_mesh->m_node, i, ElementFaces[i], FaceNodes,
                         nComposite, true);
            ++cnt;
        }
    }
    m_log(VERBOSE) << "  - # of tetrahedra: " << cnt << endl;
    nComposite++;

    // create Hexes fourth
    cnt = 0;
    for (i = 0; i < nelements; ++i)
    {
        if (ElementFaces[i].size() == 6) // if a hexa
        {
            GenElement3D(m_mesh->m_node, i, ElementFaces[i], FaceNodes,
                         nComposite, true);
            ++cnt;
        }
    }
    m_log(VERBOSE) << "  - # of hexa: " << cnt << endl;
    nComposite++;

    // Insert vertices into map.
    for (auto &node : m_mesh->m_node)
    {
        m_mesh->m_vertexSet.insert(node);
    }

    // Add boundary zones/composites
    for (i = 0; i < BndElementFaces.size(); ++i)
    {
        // m_log(VERBOSE) << "Generating 2D Zone (composite = " << nComposite
        //                << ", label = " << Facelabels[i] << ")" << endl;

        for (int j = 0; j < BndElementFaces[i].size(); ++j)
        {
            auto it = FaceNodes.find(BndElementFaces[i][j]);
            if (it != FaceNodes.end())
            {
                GenElement2D(m_mesh->m_node, j, it->second, nComposite);
            }
            else
            {
                m_log(FATAL) << "Failed to find the nodes for face "
                             << BndElementFaces[i][j] << endl;
            }
        }

        // m_mesh->m_faceLabels[nComposite] = Facelabels[i];
        nComposite++;
    }
}

static void PrismLineFaces(int prismid, map<int, int> &facelist,
                           vector<vector<int>> &FacesToPrisms,
                           vector<vector<int>> &PrismsToFaces,
                           vector<bool> &PrismDone);

/**
 * Algorithm:
 * 0. Setup
 *  - generate maps and vectors that contain information relevant for the next
 *    steps
 * 1. Pyra-prism lines
 *  - for each pyra,
 *      a) number the apex node in reverse order to ensure they have
 *         a higher global ID than the base nodes
 *      b) label prism nodes for any prism lines with a pyra on either end
 *         NOTE: this may raise an error if two pyras on opposite ends are
 *               misaligned
 * 2. Non-pyra prism lines
 *  - label prism nodes for prism lines with orientation not dictated by pyra
 * 3. Number everything else
 *  - order doesn't matter for quad faces and there is no coupling between the
 *    tri faces in tets
 */
void InputStar::ResetNodes(vector<NodeSharedPtr> &Vnodes,
                           Array<OneD, vector<int>> &ElementFaces,
                           std::unordered_map<int, vector<int>> &FaceNodes)
{
    int i, j;
    Array<OneD, int> NodeReordering(Vnodes.size(), -1);
    int face1_map[3] = {0, 1, 4};
    int face3_map[3] = {3, 2, 5};
    int nodeid       = 0;
    int revNodeid    = Vnodes.size() - 1;
    map<int, bool> FacesRenumbered;

    // Determine Prism triangular face connectivity.
    vector<vector<int>> FaceToPrisms(FaceNodes.size());
    vector<vector<int>> PrismToFaces(ElementFaces.size());
    map<int, int> Prisms;

    map<int, int> Pyras;

    // the global tri/quad face IDs for each prism/pyra
    vector<vector<int>> GlobTriFaces(ElementFaces.size());
    vector<vector<int>> GlobQuadFaces(ElementFaces.size());

    // 0.
    // generate map of prism-faces to prisms and prism to
    // triangular-faces as well as ids of each prism.
    for (i = 0; i < ElementFaces.size(); ++i)
    {
        // Find Prism (and pyramids!).
        if (ElementFaces[i].size() == 5)
        {
            vector<int> LocTriFaces;
            // Find tri and quad faces faces
            for (j = 0; j < ElementFaces[i].size(); ++j)
            {
                if (FaceNodes[ElementFaces[i][j]].size() == 3)
                {
                    LocTriFaces.push_back(j);
                    GlobTriFaces[i].push_back(ElementFaces[i][j]);
                }
                else
                {
                    GlobQuadFaces[i].push_back(ElementFaces[i][j]);
                }
            }

            if (LocTriFaces.size() == 2) // prism otherwise a pyramid
            {
                Prisms[i] = i;

                PrismToFaces[i].push_back(ElementFaces[i][LocTriFaces[0]]);
                PrismToFaces[i].push_back(ElementFaces[i][LocTriFaces[1]]);

                FaceToPrisms[ElementFaces[i][LocTriFaces[0]]].push_back(i);
                FaceToPrisms[ElementFaces[i][LocTriFaces[1]]].push_back(i);
            }
            else if (LocTriFaces.size() == 4)
            {
                Pyras[i] = i;
            }
            else
            {
                ASSERTL0(
                    false,
                    "Not set up for elements which are not Prism or Pyramid");
            }
        }
    }

    // 1.
    for (auto &PyraIt : Pyras)
    {
        // find the ID of the apex node (the node in any of the tri faces which
        // is not also in the quad (base) face)
        int apexNode          = -1;
        vector<int> baseNodes = FaceNodes[GlobQuadFaces[PyraIt.second][0]];

        // assert that none of the base nodes have already been ID'd
        for (int id : baseNodes)
        {
            ASSERTL0(NodeReordering[id] == -1,
                     "Impossible mesh: cannot enforce apex as collapsed point");
        }

        // the choice of the 0th tri face is arbitrary
        for (int id : FaceNodes[GlobTriFaces[PyraIt.second][0]])
        {
            // if a node `id` is not in the base, it must be the apex
            if (find(baseNodes.begin(), baseNodes.end(), id) == baseNodes.end())
            {
                apexNode = id;
                break;
            }
        }

        ASSERTL0(apexNode != -1, "Apex node not found in pyramid");

        // if apex node hasn't already been given an ID (because of sharing it
        // with another pyramid)
        if (NodeReordering[apexNode] == -1)
        {
            NodeReordering[apexNode] = revNodeid--;
        }

        // traverse along each prism line that exists out of the pyra
        for (int faceId : GlobTriFaces[PyraIt.second]) // for each tri face
        {
            TraversePyraPrismLine(PyraIt.second, faceId, apexNode, FaceToPrisms,
                                  GlobTriFaces, Vnodes, ElementFaces, FaceNodes,
                                  NodeReordering, revNodeid);
        }
    }

    // 2.
    // For every prism find the list of prismatic elements
    // that represent an aligned block of cells. Then renumber
    // these blocks consecutively
    vector<bool> FacesDone(FaceNodes.size(), false);
    vector<bool> PrismDone(ElementFaces.size(), false);

    for (auto &PrismIt : Prisms)
    {
        int elmtid = PrismIt.first;
        map<int, int> facelist;

        if (PrismDone[elmtid])
        {
            continue;
        }
        else
        {
            // Generate list of faces in list
            PrismLineFaces(elmtid, facelist, FaceToPrisms, PrismToFaces,
                           PrismDone);
        }
        // loop over faces and number vertices of associated prisms.
        for (auto &faceIt : facelist)
        {
            int faceid = faceIt.second;

            for (i = 0; i < FaceToPrisms[faceid].size(); ++i)
            {
                int prismid = FaceToPrisms[faceid][i];

                if ((FacesDone[PrismToFaces[prismid][0]] == true) &&
                    (FacesDone[PrismToFaces[prismid][1]] == true))
                {
                    continue;
                }

                Array<OneD, int> Nodes =
                    SortFaceNodes(Vnodes, ElementFaces[prismid], FaceNodes);

                if ((FacesDone[PrismToFaces[prismid][0]] == false) &&
                    (FacesDone[PrismToFaces[prismid][1]] == false))
                {
                    // number all nodes consecutive since
                    // already correctly re-arranged.
                    for (i = 0; i < 3; ++i)
                    {
                        if (NodeReordering[Nodes[face1_map[i]]] == -1)
                        {
                            NodeReordering[Nodes[face1_map[i]]] = nodeid++;
                        }
                    }

                    for (i = 0; i < 3; ++i)
                    {
                        if (NodeReordering[Nodes[face3_map[i]]] == -1)
                        {
                            NodeReordering[Nodes[face3_map[i]]] = nodeid++;
                        }
                    }
                }
                else if ((FacesDone[PrismToFaces[prismid][0]] == false) &&
                         (FacesDone[PrismToFaces[prismid][1]] == true))
                {
                    // find node of highest id
                    int max_id1, max_id2;

                    max_id1 = (NodeReordering[Nodes[face3_map[0]]] <
                               NodeReordering[Nodes[face3_map[1]]])
                                  ? 1
                                  : 0;
                    max_id2 = (NodeReordering[Nodes[face3_map[max_id1]]] <
                               NodeReordering[Nodes[face3_map[2]]])
                                  ? 2
                                  : max_id1;

                    // add numbering according to order of
                    int id0 = (max_id1 == 1) ? 0 : 1;

                    if (NodeReordering[Nodes[face1_map[id0]]] == -1)
                    {
                        NodeReordering[Nodes[face1_map[id0]]] = nodeid++;
                    }

                    if (NodeReordering[Nodes[face1_map[max_id1]]] == -1)
                    {
                        NodeReordering[Nodes[face1_map[max_id1]]] = nodeid++;
                    }

                    if (NodeReordering[Nodes[face1_map[max_id2]]] == -1)
                    {
                        NodeReordering[Nodes[face1_map[max_id2]]] = nodeid++;
                    }
                }
                else if ((FacesDone[PrismToFaces[prismid][0]] == true) &&
                         (FacesDone[PrismToFaces[prismid][1]] == false))
                {
                    // find node of highest id
                    int max_id1, max_id2;

                    max_id1 = (NodeReordering[Nodes[face1_map[0]]] <
                               NodeReordering[Nodes[face1_map[1]]])
                                  ? 1
                                  : 0;
                    max_id2 = (NodeReordering[Nodes[face1_map[max_id1]]] <
                               NodeReordering[Nodes[face1_map[2]]])
                                  ? 2
                                  : max_id1;

                    // add numbering according to order of
                    int id0 = (max_id1 == 1) ? 0 : 1;

                    if (NodeReordering[Nodes[face3_map[id0]]] == -1)
                    {
                        NodeReordering[Nodes[face3_map[id0]]] = nodeid++;
                    }

                    if (NodeReordering[Nodes[face3_map[max_id1]]] == -1)
                    {
                        NodeReordering[Nodes[face3_map[max_id1]]] = nodeid++;
                    }

                    if (NodeReordering[Nodes[face3_map[max_id2]]] == -1)
                    {
                        NodeReordering[Nodes[face3_map[max_id2]]] = nodeid++;
                    }
                }
            }
        }
    }

    // 3.
    // fill in any unset nodes from other shapes
    for (i = 0; i < NodeReordering.size(); ++i)
    {
        if (NodeReordering[i] == -1)
        {
            NodeReordering[i] = nodeid++;
        }
    }

    ASSERTL1(nodeid == revNodeid + 1, "Have not renumbered all nodes");

    // Renumbering successfull so reset nodes and faceNodes;
    for (auto &it : FaceNodes)
    {
        for (j = 0; j < it.second.size(); ++j)
        {
            it.second[j] = NodeReordering[it.second[j]];
        }
    }

    vector<NodeSharedPtr> save(Vnodes);
    for (i = 0; i < Vnodes.size(); ++i)
    {
        Vnodes[NodeReordering[i]] = save[i];
        Vnodes[NodeReordering[i]]->SetID(NodeReordering[i]);
    }
}

void InputStar::TraversePyraPrismLine(
    int currElemId, int currFaceId, int currApexNode,
    std::vector<std::vector<int>> FaceToPrisms,
    std::vector<std::vector<int>> GlobTriFaces,
    std::vector<NodeSharedPtr> &Vnodes,
    Array<OneD, std::vector<int>> &ElementFaces,
    std::unordered_map<int, std::vector<int>> &FaceNodes,
    Array<OneD, int> &NodeReordering, int &revNodeid)
{
    // these convert between the node ordering generated by OrderFaceNodes()
    // and the node ordering for the two triangular faces
    int face1_map[3] = {0, 1, 4};
    int face3_map[3] = {3, 2, 5};

    // iterate over the one/two prisms that connect the face
    for (int nextElemId : FaceToPrisms[currFaceId])
    {
        // only carry out the rest of the function if travelling in the right
        // direction
        if (nextElemId == currElemId)
        {
            continue;
        }

        // iterate over the two tri faces (one at each end)
        for (int nextFaceId : GlobTriFaces[nextElemId])
        {
            // ensure that we are travelling in the right direction
            if (nextFaceId == currFaceId)
            {
                continue;
            }

            // reorder the nodes
            Array<OneD, int> Nodes =
                SortFaceNodes(Vnodes, ElementFaces[nextElemId], FaceNodes);

            int nextApexNode = -1; // the ID of the opposite collapsed point
            int currFace     = -1; // either 1 or 3

            // highest ID for each face
            int highestId1 = -1;
            int highestId3 = -1;

            // find the orientation of the prism (face 1 vs 3)
            for (int i = 0; i < 3; ++i)
            {
                highestId1 = (NodeReordering[Nodes[face1_map[i]]] > highestId1)
                                 ? NodeReordering[Nodes[face1_map[i]]]
                                 : highestId1;

                highestId3 = (NodeReordering[Nodes[face3_map[i]]] > highestId3)
                                 ? NodeReordering[Nodes[face3_map[i]]]
                                 : highestId3;

                if (Nodes[face1_map[i]] == currApexNode)
                {
                    nextApexNode = Nodes[face3_map[i]];
                    currFace     = 1;
                }

                if (Nodes[face3_map[i]] == currApexNode)
                {
                    nextApexNode = Nodes[face1_map[i]];
                    currFace     = 3;
                }
            }

            if (currFace == 1)
            {
                ASSERTL0(highestId3 == NodeReordering[nextApexNode],
                         "Error meshing: incompatible collapsed point");
            }
            else if (currFace == 3)
            {
                ASSERTL0(highestId1 == NodeReordering[nextApexNode],
                         "Error meshing: incompatible collapsed point");
            }
            else
            {
                m_log(FATAL) << "Error meshing: collapsed point not found";
            }

            if (NodeReordering[nextApexNode] == -1)
            {
                NodeReordering[nextApexNode] = revNodeid--;
            }

            // continue traversing along the line
            TraversePyraPrismLine(nextElemId, nextFaceId, -1, FaceToPrisms,
                                  GlobTriFaces, Vnodes, ElementFaces, FaceNodes,
                                  NodeReordering, revNodeid);
        }
    }
}

static void PrismLineFaces(int prismid, map<int, int> &facelist,
                           vector<vector<int>> &FaceToPrisms,
                           vector<vector<int>> &PrismToFaces,
                           vector<bool> &PrismDone)
{
    if (PrismDone[prismid] == false)
    {
        PrismDone[prismid] = true;

        // Add faces0
        int face       = PrismToFaces[prismid][0];
        facelist[face] = face;
        for (int i = 0; i < FaceToPrisms[face].size(); ++i)
        {
            PrismLineFaces(FaceToPrisms[face][i], facelist, FaceToPrisms,
                           PrismToFaces, PrismDone);
        }

        // Add faces1
        face           = PrismToFaces[prismid][1];
        facelist[face] = face;
        for (int i = 0; i < FaceToPrisms[face].size(); ++i)
        {
            PrismLineFaces(FaceToPrisms[face][i], facelist, FaceToPrisms,
                           PrismToFaces, PrismDone);
        }
    }
}

void InputStar::GenElement2D(vector<NodeSharedPtr> &VertNodes,
                             [[maybe_unused]] int i, vector<int> &FaceNodes,
                             int nComposite)
{
    LibUtilities::ShapeType elType = LibUtilities::eTriangle;

    if (FaceNodes.size() == 4)
    {
        elType = LibUtilities::eQuadrilateral;
    }
    else if (FaceNodes.size() != 3)
    {
        m_log(FATAL) << "Not set up for elements which are not tets or prisms"
                     << endl;
    }

    // Create element tags
    vector<int> tags;
    tags.push_back(nComposite);

    // make unique node list
    vector<NodeSharedPtr> nodeList;
    Array<OneD, int> Nodes = SortEdgeNodes(VertNodes, FaceNodes);
    for (int j = 0; j < Nodes.size(); ++j)
    {
        nodeList.push_back(VertNodes[Nodes[j]]);
    }

    // Create element
    ElmtConfig conf(elType, 1, true, true);
    ElementSharedPtr E =
        GetElementFactory().CreateInstance(elType, conf, nodeList, tags);

    m_mesh->m_element[E->GetDim()].push_back(E);
}

void InputStar::GenElement3D(vector<NodeSharedPtr> &VertNodes,
                             [[maybe_unused]] int i, vector<int> &ElementFaces,
                             std::unordered_map<int, vector<int>> &FaceNodes,
                             int nComposite, bool DoOrient)
{
    LibUtilities::ShapeType elType = LibUtilities::eTetrahedron;

    // set up Node list
    Array<OneD, int> Nodes = SortFaceNodes(VertNodes, ElementFaces, FaceNodes);
    int nnodes             = Nodes.size();
    map<LibUtilities::ShapeType, int> domainComposite;

    // element type
    if (nnodes == 5)
    {
        elType = LibUtilities::ePyramid;
    }
    else if (nnodes == 6)
    {
        elType = LibUtilities::ePrism;
    }
    else if (nnodes == 8)
    {
        elType = LibUtilities::eHexahedron;
    }
    else if (nnodes != 4)
    {
        m_log(FATAL) << "Unknown element type" << endl;
    }

    // Create element tags
    vector<int> tags;
    tags.push_back(nComposite);

    // make unique node list
    vector<NodeSharedPtr> nodeList;
    for (int j = 0; j < Nodes.size(); ++j)
    {
        nodeList.push_back(VertNodes[Nodes[j]]);
    }

    ElmtConfig conf(elType, 1, true, true, DoOrient);
    ElementSharedPtr E =
        GetElementFactory().CreateInstance(elType, conf, nodeList, tags);

    m_mesh->m_element[E->GetDim()].push_back(E);
}

Array<OneD, int> InputStar::SortEdgeNodes(vector<NodeSharedPtr> &Vnodes,
                                          vector<int> &FaceNodes)
{
    Array<OneD, int> returnval;

    if (FaceNodes.size() == 3) // Triangle
    {
        returnval = Array<OneD, int>(3);

        returnval[0] = FaceNodes[0];
        returnval[1] = FaceNodes[1];
        returnval[2] = FaceNodes[2];
    }
    else if (FaceNodes.size() == 4) // quadrilateral
    {
        returnval = Array<OneD, int>(4);

        int indx0 = FaceNodes[0];
        int indx1 = FaceNodes[1];
        int indx2 = FaceNodes[2];
        int indx3 = FaceNodes[3];

        // calculate 0-1,
        Node a = *(Vnodes[indx1]) - *(Vnodes[indx0]);
        // calculate 0-2,
        Node b      = *(Vnodes[indx2]) - *(Vnodes[indx0]);
        Node acurlb = a.curl(b);

        // calculate 2-1,
        Node c = *(Vnodes[indx1]) - *(Vnodes[indx2]);
        // calculate 3-2,
        Node d      = *(Vnodes[indx3]) - *(Vnodes[indx2]);
        Node acurld = a.curl(d);

        NekDouble acurlb_dot_acurld = acurlb.dot(acurld);
        if (acurlb_dot_acurld > 0.0)
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[2] = indx2;
            returnval[3] = indx3;
        }
        else
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[2] = indx3;
            returnval[3] = indx2;
        }
    }

    return returnval;
}

Array<OneD, int> InputStar::SortFaceNodes(
    vector<NodeSharedPtr> &Vnodes, vector<int> &ElementFaces,
    std::unordered_map<int, vector<int>> &FaceNodes)
{
    int i, j;
    Array<OneD, int> returnval;

    if (ElementFaces.size() == 4) // Tetrahedron
    {
        ASSERTL1(FaceNodes[ElementFaces[0]].size() == 3,
                 "Face is not triangular");

        returnval = Array<OneD, int>(4);

        auto it   = FaceNodes.find(ElementFaces[0]);
        int indx0 = it->second[0];
        int indx1 = it->second[1];
        int indx2 = it->second[2];
        int indx3 = -1;

        // calculate 0-1,
        Node a = *(Vnodes[indx1]) - *(Vnodes[indx0]);
        // calculate 0-2,
        Node b = *(Vnodes[indx2]) - *(Vnodes[indx0]);

        // Find fourth node index;
        ASSERTL1(FaceNodes[ElementFaces[1]].size() == 3,
                 "Face is not triangular");

        auto it2 = FaceNodes.find(ElementFaces[1]);
        for (i = 0; i < 3; ++i)
        {
            if ((it2->second[i] != indx0) && (it2->second[i] != indx1) &&
                (it2->second[i] != indx2))
            {
                indx3 = it2->second[i];
                break;
            }
        }

        // calculate 0-3,
        Node c      = *(Vnodes[indx3]) - *(Vnodes[indx0]);
        Node acurlb = a.curl(b);

        NekDouble acurlb_dotc = acurlb.dot(c);
        if (acurlb_dotc < 0.0)
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[2] = indx2;
            returnval[3] = indx3;
        }
        else
        {
            returnval[0] = indx1;
            returnval[1] = indx0;
            returnval[2] = indx2;
            returnval[3] = indx3;
        }
    }
    else if (ElementFaces.size() == 5) // prism or pyramid
    {
        int triface0, triface1, triface2, triface3;
        int quadface0, quadface1, quadface2;
        bool isPrism = true;

        // find ids of tri faces and first quad face
        triface0 = triface1 = triface2 = triface3 = -1;
        quadface0 = quadface1 = quadface2 = -1;
        for (i = 0; i < 5; ++i)
        {
            auto it = FaceNodes.find(ElementFaces[i]);
            if (it->second.size() == 3)
            {
                if (triface0 == -1)
                {
                    triface0 = i;
                }
                else if (triface1 == -1)
                {
                    triface1 = i;
                }
                else if (triface2 == -1)
                {
                    triface2 = i;
                    isPrism  = false;
                }
                else if (triface3 == -1)
                {
                    triface3 = i;
                }
            }

            if (it->second.size() == 4)
            {
                if (quadface0 == -1)
                {
                    quadface0 = i;
                }
                else if (quadface1 == -1)
                {
                    quadface1 = i;
                }
                else if (quadface2 == -1)
                {
                    quadface2 = i;
                }
            }
        }

        if (isPrism) // Prism
        {
            returnval = Array<OneD, int>(6);
            ASSERTL1(quadface0 != -1, "Quad face 0 not found");
            ASSERTL1(quadface1 != -1, "Quad face 1 not found");
            ASSERTL1(quadface2 != -1, "Quad face 2 not found");
            ASSERTL1(triface0 != -1, "Tri face 0 not found");
            ASSERTL1(triface1 != -1, "Tri face 1 not found");
        }
        else // Pyramid
        {
            returnval = Array<OneD, int>(5);
            ASSERTL1(quadface0 != -1, "Quad face 0 not found");
            ASSERTL1(triface0 != -1, "Tri face 0 not found");
            ASSERTL1(triface1 != -1, "Tri face 1 not found");
            ASSERTL1(triface2 != -1, "Tri face 2 not found");
            ASSERTL1(triface3 != -1, "Tri face 3 not found");
        }

        // find matching nodes between triface0 and triquad0
        int indx0, indx1, indx2, indx3, indx4;

        indx0 = indx1 = indx2 = indx3 = indx4 = -1;
        // Loop over all quad nodes and if they match any
        // triangular nodes If they do set these to indx0 and
        // indx1 and if not set it to indx2, indx3

        auto &triface0_vec  = FaceNodes.find(ElementFaces[triface0])->second;
        auto &quadface0_vec = FaceNodes.find(ElementFaces[quadface0])->second;
        for (i = 0; i < 4; ++i)
        {
            for (j = 0; j < 3; ++j)
            {
                if (triface0_vec[j] == quadface0_vec[i])
                {
                    break; // same node break
                }
            }

            if (j == 3) // Vertex not in quad face
            {
                if (indx2 == -1)
                {
                    indx2 = quadface0_vec[i];
                }
                else if (indx3 == -1)
                {
                    indx3 = quadface0_vec[i];
                }
                else
                {
                    m_log(FATAL) << "More than two vertices do not match "
                                 << "triangular face" << endl;
                }
            }
            else // if found match then set indx0,indx1;
            {
                if (indx0 == -1)
                {
                    indx0 = quadface0_vec[i];
                }
                else
                {
                    indx1 = quadface0_vec[i];
                }
            }
        }

        // Finally check for top vertex
        for (int i = 0; i < 3; ++i)
        {
            if (triface0_vec[i] != indx0 && triface0_vec[i] != indx1 &&
                triface0_vec[i] != indx2)
            {
                indx4 = triface0_vec[i];
                break;
            }
        }

        // calculate 0-1,
        Node a = *(Vnodes[indx1]) - *(Vnodes[indx0]);
        // calculate 0-4,
        Node b = *(Vnodes[indx4]) - *(Vnodes[indx0]);
        // calculate 0-2,
        Node c      = *(Vnodes[indx2]) - *(Vnodes[indx0]);
        Node acurlb = a.curl(b);

        NekDouble acurlb_dotc = acurlb.dot(c);
        if (acurlb_dotc < 0.0)
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[4] = indx4;
        }
        else
        {
            returnval[0] = indx1;
            returnval[1] = indx0;
            returnval[4] = indx4;
        }

        // check to see if two vertices are shared between one of the other
        // faces
        // to define which is indx2 and indx3
        if (isPrism == true)
        {
            auto &quadface1_vec =
                FaceNodes.find(ElementFaces[quadface1])->second;
            auto &quadface2_vec =
                FaceNodes.find(ElementFaces[quadface2])->second;
            int cnt = 0;
            for (int i = 0; i < 4; ++i)
            {
                if (quadface1_vec[i] == returnval[1] ||
                    quadface1_vec[i] == indx2)
                {
                    cnt++;
                }
            }

            if (cnt == 2) // have two matching vertices
            {
                returnval[2] = indx2;
                returnval[3] = indx3;
            }
            else
            {
                cnt = 0;
                for (int i = 0; i < 4; ++i)
                {
                    if (quadface2_vec[i] == returnval[1] ||
                        quadface2_vec[i] == indx2)
                    {
                        cnt++;
                    }
                }

                if (cnt != 2) // neither of the other faces has two matching
                              // nodes so reverse
                {
                    returnval[2] = indx3;
                    returnval[3] = indx2;
                }
                else // have two matching vertices
                {
                    returnval[2] = indx2;
                    returnval[3] = indx3;
                }
            }
            // finally need to find last vertex from second triangular face.
            auto &triface1_vec = FaceNodes.find(ElementFaces[triface1])->second;
            for (int i = 0; i < 3; ++i)
            {
                if (triface1_vec[i] != indx2 && triface1_vec[i] != indx3)
                {
                    returnval[5] = triface1_vec[i];
                    break;
                }
            }
        }
        else
        {
            vector<int> trifaceid;
            trifaceid.push_back(triface1);
            trifaceid.push_back(triface2);
            trifaceid.push_back(triface3);
            int cnt = 0;
            for (auto id : trifaceid)
            {
                auto &triface_vec = FaceNodes.find(ElementFaces[id])->second;
                cnt               = 0;
                for (int i = 0; i < 3; ++i)
                {
                    if (triface_vec[i] == returnval[1] ||
                        triface_vec[i] == indx2)
                    {
                        cnt++;
                    }
                }

                if (cnt == 2) // have two matching vertices
                {
                    returnval[2] = indx2;
                    returnval[3] = indx3;
                    break;
                }
            }
            if (cnt != 2)
            {
                returnval[2] = indx3;
                returnval[3] = indx2;
            }
        }
    }
    else if (ElementFaces.size() == 6) // hexahedron
    {
        ASSERTL1(FaceNodes[ElementFaces[0]].size() == 4 &&
                     FaceNodes[ElementFaces[1]].size() == 4 &&
                     FaceNodes[ElementFaces[2]].size() == 4 &&
                     FaceNodes[ElementFaces[3]].size() == 4 &&
                     FaceNodes[ElementFaces[4]].size() == 4 &&
                     FaceNodes[ElementFaces[5]].size() == 4,
                 "Shape not recognised");

        returnval = Array<OneD, int>(8);

        // neighbours map contains a set of each neighbouring node for every
        // node
        map<int, std::unordered_set<int>> neighbours;

        // populate neighbours
        for (int faceId : ElementFaces)
        {
            vector<int> nodes = FaceNodes[faceId];
            neighbours[nodes[0]].insert({nodes[1], nodes[3]});
            neighbours[nodes[1]].insert({nodes[2], nodes[0]});
            neighbours[nodes[2]].insert({nodes[3], nodes[1]});
            neighbours[nodes[3]].insert({nodes[0], nodes[2]});
        }

        // the bottom face is chosen arbitrarily as the 0th indexed one
        vector<int> bottomFace = FaceNodes[ElementFaces[0]];
        std::unordered_set<int> bottomFaceSet(bottomFace.begin(),
                                              bottomFace.end());

        // the top face is the face opposite, sharing no edges or vertices
        vector<int> topFace;

        // find the top face
        for (int i = 1; i < 6; i++)
        {
            bool sharedEdge = false;
            for (int id : FaceNodes[ElementFaces[i]])
            {
                if (bottomFaceSet.count(id)) // if node @ id is in the
                                             // bottom face
                {
                    sharedEdge = true;
                    break;
                }
            }
            if (sharedEdge)
            {
                continue;
            }
            topFace = FaceNodes[ElementFaces[i]];
        }

        ASSERTL0(topFace.size() == 4, "Top face of hex not found");

        // insert the bottom face nodes into returnval[0:4]
        // insert the top face nodes into returnval[4:]
        // returnval[i+4] is the node opposite returnval[i]
        for (int i = 0; i < 4; i++)
        {
            returnval[i]     = bottomFace[i];
            bool oppositeSet = false;

            for (int id : neighbours[bottomFace[i]])
            {
                if (!bottomFaceSet.count(id))
                {
                    returnval[i + 4] = id;
                    oppositeSet      = true;
                    break;
                }
            }
            ASSERTL0(oppositeSet, "Node in the top face of hex not numbered");
        }
    }
    else
    {
        m_log(FATAL) << "SortFaceNodes not set up for this number of faces"
                     << endl;
    }

    return returnval;
}

// initialise and read ccm file to ccm structure
void InputStar::InitCCM(void)
{
    // Open ccm file for reading.
    CCMIOID root;
    // Open the file.  Because we did not initialize 'err' we
    // need to pass in NULL (which always means kCCMIONoErr)
    // and then assign the return value to 'err'.).
    string fname = m_config["infile"].as<string>();
    m_ccmErr     = CCMIOOpenFile(nullptr, fname.c_str(), kCCMIORead, &root);

    if (m_ccmErr != kCCMIONoErr)
    {
        m_log(FATAL) << "Error opening file '" << fname << "'" << endl;
    }

    int i = 0;
    CCMIOID state, problem;

    // We are going to assume that we have a state with a
    // known name.  We could instead use CCMIONextEntity() to
    // walk through all the states in the file and present the
    // list to the user for selection.
    CCMIOGetState(&m_ccmErr, root, kDefaultState, &problem, &state);
    if (m_ccmErr != kCCMIONoErr)
    {
        m_log(FATAL) << "No state named '" << kDefaultState << "'" << endl;
    }

    // Find the first processor (i has previously been
    // initialized to 0) and read the mesh and solution
    // information.
    CCMIONextEntity(&m_ccmErr, state, kCCMIOProcessor, &i, &m_ccmProcessor);

    if (m_ccmErr != kCCMIONoErr)
    {
        m_log(FATAL) << "CCM error: failed to find next entity" << endl;
    }
}

void InputStar::ReadNodes(std::vector<NodeSharedPtr> &Nodes)
{
    CCMIOID mapID, vertices;
    CCMIOSize nVertices;
    int dims = 1;

    CCMIOReadProcessor(&m_ccmErr, m_ccmProcessor, &vertices, &m_ccmTopology,
                       nullptr, nullptr);

    if (m_ccmErr != kCCMIONoErr)
    {
        m_log(FATAL) << "CCM error: error reading processor" << endl;
    }

    CCMIOEntitySize(&m_ccmErr, vertices, &nVertices, nullptr);
    if (m_ccmErr != kCCMIONoErr)
    {
        m_log(FATAL) << "CCM error: error reading NextEntitySize in ReadNodes"
                     << endl;
    }

    // Read the vertices.  This involves reading both the vertex data and
    // the map, which maps the index into the data array with the ID number.
    // As we process the vertices we need to be sure to scale them by the
    // appropriate scaling factor.  The offset is just to show you can read
    // any chunk.  Normally this would be in a for loop.
    float scale;
    int nvert = nVertices;
    vector<int> mapData;
    mapData.resize(nvert);
    vector<float> verts;
    verts.resize(3 * nvert);

    for (int k = 0; k < nvert; ++k)
    {
        verts[3 * k] = verts[3 * k + 1] = verts[3 * k + 2] = 0.0;
        mapData[k]                                         = 0;
    }

    CCMIOReadVerticesf(&m_ccmErr, vertices, &dims, &scale, &mapID, &verts[0], 0,
                       nVertices);

    if (m_ccmErr != kCCMIONoErr)
    {
        m_log(FATAL) << "Error reading vertices from CCM file in ReadNodes"
                     << endl;
    }

    CCMIOReadMap(&m_ccmErr, mapID, &mapData[0], 0, nVertices);

    if (m_ccmErr != kCCMIONoErr)
    {
        m_log(FATAL) << "Error reading node map in ReadNodes" << endl;
    }

    for (int i = 0; i < nVertices; ++i)
    {
        Nodes.push_back(std::make_shared<Node>(
            i, verts[3 * i], verts[3 * i + 1], verts[3 * i + 2]));
    }
}

void InputStar::ReadInternalFaces(
    std::unordered_map<int, vector<int>> &FacesNodes,
    Array<OneD, vector<int>> &ElementFaces)
{

    CCMIOID mapID, id;
    CCMIOSize nFaces, size;
    vector<int> faces, faceCells, mapData;

    // Read the internal faces.
    CCMIOGetEntity(&m_ccmErr, m_ccmTopology, kCCMIOInternalFaces, 0, &id);
    CCMIOEntitySize(&m_ccmErr, id, &nFaces, nullptr);

    int nf = nFaces;
    mapData.resize(nf);
    faceCells.resize(2 * nf);

    CCMIOReadFaces(&m_ccmErr, id, kCCMIOInternalFaces, nullptr, &size, nullptr,
                   kCCMIOStart, kCCMIOEnd);
    faces.resize((size_t)size);
    CCMIOReadFaces(&m_ccmErr, id, kCCMIOInternalFaces, &mapID, nullptr,
                   &faces[0], kCCMIOStart, kCCMIOEnd);
    CCMIOReadFaceCells(&m_ccmErr, id, kCCMIOInternalFaces, &faceCells[0],
                       kCCMIOStart, kCCMIOEnd);
    CCMIOReadMap(&m_ccmErr, mapID, &mapData[0], kCCMIOStart, kCCMIOEnd);

    // Add face nodes
    int cnt = 0;
    for (int i = 0; i < nf; ++i)
    {
        vector<int> Fnodes;
        int j;
        if (cnt < faces.size())
        {
            int nv = faces[cnt];
            if (nv > 4)
            {
                m_log(FATAL) << "Can only handle meshes with up to four nodes "
                             << "per face" << endl;
            }

            for (j = 0; j < nv; ++j)
            {
                if (cnt + 1 + j < faces.size())
                {
                    Fnodes.push_back(faces[cnt + 1 + j] - 1);
                }
            }
            cnt += nv + 1;
        }
        FacesNodes[mapData[i] - 1] = Fnodes;
    }

    // find number of elements;
    int nelmt = 0;
    for (int i = 0; i < faceCells.size(); ++i)
    {
        nelmt = max(nelmt, faceCells[i]);
    }

    ElementFaces = Array<OneD, vector<int>>(nelmt);
    for (int i = 0; i < nf; ++i)
    {
        // left element
        if (faceCells[2 * i])
        {
            ElementFaces[faceCells[2 * i] - 1].push_back(mapData[i] - 1);
        }

        // right element
        if (faceCells[2 * i + 1])
        {
            ElementFaces[faceCells[2 * i + 1] - 1].push_back(mapData[i] - 1);
        }
    }
}

void InputStar::ReadBoundaryFaces(
    vector<vector<int>> &BndElementFaces,
    std::unordered_map<int, vector<int>> &FacesNodes,
    Array<OneD, vector<int>> &ElementFaces, vector<string> &Facelabels)
{
    // Read the boundary faces.
    int index = 0;
    CCMIOID mapID, id;
    CCMIOSize nFaces, size;
    vector<int> faces, faceCells, mapData;
    vector<string> facelabel;

    while (CCMIONextEntity(nullptr, m_ccmTopology, kCCMIOBoundaryFaces, &index,
                           &id) == kCCMIONoErr)
    {
        int boundaryVal;

        CCMIOEntitySize(&m_ccmErr, id, &nFaces, nullptr);
        CCMIOSize nf = nFaces;
        mapData.resize(nf);
        faceCells.resize(nf);
        CCMIOReadFaces(&m_ccmErr, id, kCCMIOBoundaryFaces, nullptr, &size,
                       nullptr, kCCMIOStart, kCCMIOEnd);

        faces.resize((size_t)size);
        CCMIOReadFaces(&m_ccmErr, id, kCCMIOBoundaryFaces, &mapID, nullptr,
                       &faces[0], kCCMIOStart, kCCMIOEnd);
        CCMIOReadFaceCells(&m_ccmErr, id, kCCMIOBoundaryFaces, &faceCells[0],
                           kCCMIOStart, kCCMIOEnd);
        CCMIOReadMap(&m_ccmErr, mapID, &mapData[0], kCCMIOStart, kCCMIOEnd);

        CCMIOGetEntityIndex(&m_ccmErr, id, &boundaryVal);

        // check to see if we have a label for this boundary faces
        int size;
        char *name;
        if (CCMIOReadOptstr(nullptr, id, "Label", &size, nullptr) ==
            kCCMIONoErr)
        {
            name = new char[size + 1];
            CCMIOReadOptstr(nullptr, id, "Label", nullptr, name);
            Facelabels.push_back(string(name));
        }
        else
        {
            Facelabels.push_back("Not known");
        }

        // Add face nodes
        int cnt = 0;
        for (int i = 0; i < nf; ++i)
        {
            vector<int> Fnodes;
            int j;
            if (cnt < faces.size())
            {
                int nv = faces[cnt];

                if (nv > 4)
                {
                    m_log(FATAL) << "Can only handle meshes with up to four "
                                 << "nodes per face" << endl;
                }

                for (j = 0; j < nv; ++j)
                {
                    if (cnt + 1 + j < faces.size())
                    {
                        Fnodes.push_back(faces[cnt + 1 + j] - 1);
                    }
                }
                cnt += nv + 1;
            }
            FacesNodes[mapData[i] - 1] = Fnodes;
        }

        vector<int> BndFaces;
        for (int i = 0; i < nf; ++i)
        {
            if (faceCells[i])
            {
                ElementFaces[faceCells[i] - 1].push_back(mapData[i] - 1);
            }
            BndFaces.push_back(mapData[i] - 1);
        }
        BndElementFaces.push_back(BndFaces);
    }
}
} // namespace NekMesh
} // namespace Nektar
