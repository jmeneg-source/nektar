////////////////////////////////////////////////////////////////////////////////
//
//  File: InputCGNS.h
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
//  Description: CGNS converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_INPUTCGNS
#define UTILITIES_NEKMESH_INPUTCGNS

#include <NekMesh/Module/Module.h>
#include <cgnslib.h>

namespace Nektar::NekMesh
{

/**
 * Converter for CGNS files.
 */

class InputCGNS : public NekMesh::InputModule
{
public:
    InputCGNS(NekMesh::MeshSharedPtr m);
    ~InputCGNS() override;
    void Process() override;

    /// Creates an instance of this class
    static NekMesh::ModuleSharedPtr create(NekMesh::MeshSharedPtr m)
    {
        return MemoryManager<InputCGNS>::AllocateSharedPtr(m);
    }
    /// ModuleKey for class.
    static NekMesh::ModuleKey className;

    string GetModuleName() override
    {
        return "InputCGNS";
    }

private:
    void SetupElements();
    void ResetNodes(vector<NodeSharedPtr> &Vnodes,
                    Array<OneD, vector<int>> &ElementFaces,
                    std::unordered_map<int, vector<int>> &FaceNodes,
                    vector<pair<ElementType_t, vector<int>>> &elemInfo,
                    vector<int> &VolumeElems);

    Array<OneD, int> SortFaceNodes(
        vector<NodeSharedPtr> &Vnodes, vector<int> &ElementFaces,
        std::unordered_map<int, vector<int>> &FaceNodes);

    void GenElement2D(vector<NodeSharedPtr> &VertNodes, vector<int> elemNodes,
                      ElementType_t elemType, int nComposite, bool DoOrient);

    void GenElement3D(vector<NodeSharedPtr> &VertNodes, vector<int> elemNodes,
                      vector<int> &ElementFaces,
                      std::unordered_map<int, vector<int>> &FaceNodes,
                      ElementType_t elemType, int nComposite, bool DoOrient);

    vector<int> CGNSReordering(LibUtilities::ShapeType shapeType, int order);

    Array<OneD, int> SortEdgeNodes(vector<NodeSharedPtr> &Vnodes,
                                   vector<int> &FaceNodes);

    void ReadFaces(vector<pair<ElementType_t, vector<int>>> &elemInfo,
                   std::unordered_map<int, vector<int>> &FaceNodes,
                   Array<OneD, vector<int>> &ElementFaces,
                   vector<int> &BoundaryElems, vector<int> &VolumeElems);

    void ExtractSeparatedElemInfo(
        cgsize_t elemDataSize, int elemSize, ElementType_t elemType,
        int sectionInd, vector<pair<ElementType_t, vector<int>>> &elemInfo);

    void ExtractMixedElemInfo(
        cgsize_t elemDataSize, int elemSize, int sectionInd,
        vector<pair<ElementType_t, vector<int>>> &elemInfo);

    void SaveNode(int id, NekDouble x = 0, NekDouble y = 0, NekDouble z = 0);

    void PyramidShielding(vector<NodeSharedPtr> &Vnodes,
                          Array<OneD, vector<int>> &ElementFaces,
                          std::unordered_map<int, vector<int>> &FaceNodes,
                          vector<pair<ElementType_t, vector<int>>> &elemInfo,
                          vector<int> &VolumeElems, vector<int> &NodeReordering,
                          int pyraElemIdx);

    void TraversePyraPrismLine(int currElemId, int currFaceId, int currApexNode,
                               vector<vector<int>> FaceToPrisms,
                               vector<vector<int>> GlobTriFaces,
                               vector<NodeSharedPtr> &Vnodes,
                               Array<OneD, vector<int>> &ElementFaces,
                               std::unordered_map<int, vector<int>> &FaceNodes,
                               vector<int> &NodeReordering, int &revNodeid);

    int m_fileIndex;
    int m_baseIndex = 1;
    int m_zoneIndex = 1;
    vector<double> m_x; // x-coordinates from CGNS file
    vector<double> m_y; // y-coordinates from CGNS file
    vector<double> m_z; // z-coordinates from CGNS file

    std::unordered_map<ElementType_t, vector<int>> m_orderingMap;

    std::unordered_map<LibUtilities::ShapeType, vector<vector<int>>>
        m_shapeType2LocElemNodes{
            // Gives the ordering of the nodes for each face for each shape
            // type.
            // This ordering is to be compatible with the code in from
            // InputStar.ccp
            {LibUtilities::eTetrahedron,
             {{0, 2, 1}, {3, 1, 2}, {3, 0, 1}, {3, 2, 0}}},
            {LibUtilities::ePyramid,
             {{1, 4, 0}, {2, 4, 1}, {3, 4, 2}, {0, 4, 3}, {3, 2, 1, 0}}},
            {LibUtilities::ePrism,
             {{2, 0, 3, 5}, {2, 1, 0}, {0, 1, 4, 3}, {3, 4, 5}, {1, 2, 5, 4}}},
            {LibUtilities::eHexahedron,
             {{0, 1, 2, 3},
              {0, 4, 5, 1},
              {4, 0, 3, 7},
              {7, 3, 2, 6},
              {2, 1, 5, 6},
              {4, 7, 6, 5}}},
        };

    map<int, ElementType_t> m_ind2ElemType{
        // converts from the element type index to the type
        // {CG_Null, ElementTypeNull},     {CG_UserDefined,
        // ElementTypeUserDefined},
        {2, NODE},      {3, BAR_2},       {4, BAR_3},       {5, TRI_3},
        {6, TRI_6},     {7, QUAD_4},      {8, QUAD_8},      {9, QUAD_9},
        {10, TETRA_4},  {11, TETRA_10},   {12, PYRA_5},     {13, PYRA_14},
        {14, PENTA_6},  {15, PENTA_15},   {16, PENTA_18},   {17, HEXA_8},
        {18, HEXA_20},  {19, HEXA_27},    {20, MIXED},      {21, PYRA_13},
        {22, NGON_n},   {23, NFACE_n},    {24, BAR_4},      {25, TRI_9},
        {26, TRI_10},   {27, QUAD_12},    {28, QUAD_16},    {29, TETRA_16},
        {30, TETRA_20}, {31, PYRA_21},    {32, PYRA_29},    {33, PYRA_30},
        {34, PENTA_24}, {35, PENTA_38},   {36, PENTA_40},   {37, HEXA_32},
        {38, HEXA_56},  {39, HEXA_64},    {40, BAR_5},      {41, TRI_12},
        {42, TRI_15},   {43, QUAD_P4_16}, {44, QUAD_25},    {45, TETRA_22},
        {46, TETRA_34}, {47, TETRA_35},   {48, PYRA_P4_29}, {49, PYRA_50},
        {50, PYRA_55},  {51, PENTA_33},   {52, PENTA_66},   {53, PENTA_75},
        {54, HEXA_44},  {55, HEXA_98},    {56, HEXA_125}};

    map<ElementType_t, LibUtilities::ShapeType> m_elemType2ShapeType{
        {NODE, LibUtilities::ePoint},
        {BAR_2, LibUtilities::eSegment},
        {BAR_3, LibUtilities::eSegment},
        {BAR_4, LibUtilities::eSegment},
        {BAR_5, LibUtilities::eSegment},
        {TRI_3, LibUtilities::eTriangle},
        {TRI_6, LibUtilities::eTriangle},
        {TRI_9, LibUtilities::eTriangle},
        {TRI_10, LibUtilities::eTriangle},
        {TRI_12, LibUtilities::eTriangle},
        {TRI_15, LibUtilities::eTriangle},
        {QUAD_4, LibUtilities::eQuadrilateral},
        {QUAD_8, LibUtilities::eQuadrilateral},
        {QUAD_9, LibUtilities::eQuadrilateral},
        {QUAD_12, LibUtilities::eQuadrilateral},
        {QUAD_16, LibUtilities::eQuadrilateral},
        {QUAD_P4_16, LibUtilities::eQuadrilateral},
        {QUAD_25, LibUtilities::eQuadrilateral},
        {TETRA_4, LibUtilities::eTetrahedron},
        {TETRA_10, LibUtilities::eTetrahedron},
        {TETRA_16, LibUtilities::eTetrahedron},
        {TETRA_20, LibUtilities::eTetrahedron},
        {TETRA_22, LibUtilities::eTetrahedron},
        {TETRA_34, LibUtilities::eTetrahedron},
        {TETRA_35, LibUtilities::eTetrahedron},
        {PYRA_5, LibUtilities::ePyramid},
        {PYRA_13, LibUtilities::ePyramid},
        {PYRA_14, LibUtilities::ePyramid},
        {PYRA_21, LibUtilities::ePyramid},
        {PYRA_29, LibUtilities::ePyramid},
        {PYRA_30, LibUtilities::ePyramid},
        {PYRA_P4_29, LibUtilities::ePyramid},
        {PYRA_50, LibUtilities::ePyramid},
        {PYRA_55, LibUtilities::ePyramid},
        {PENTA_6, LibUtilities::ePrism},
        {PENTA_15, LibUtilities::ePrism},
        {PENTA_18, LibUtilities::ePrism},
        {PENTA_24, LibUtilities::ePrism},
        {PENTA_38, LibUtilities::ePrism},
        {PENTA_40, LibUtilities::ePrism},
        {PENTA_33, LibUtilities::ePrism},
        {PENTA_66, LibUtilities::ePrism},
        {PENTA_75, LibUtilities::ePrism},
        {HEXA_8, LibUtilities::eHexahedron},
        {HEXA_20, LibUtilities::eHexahedron},
        {HEXA_27, LibUtilities::eHexahedron},
        {HEXA_32, LibUtilities::eHexahedron},
        {HEXA_56, LibUtilities::eHexahedron},
        {HEXA_64, LibUtilities::eHexahedron},
        {HEXA_44, LibUtilities::eHexahedron},
        {HEXA_98, LibUtilities::eHexahedron},
        {HEXA_125, LibUtilities::eHexahedron},
    };

    map<LibUtilities::ShapeType, uint> m_shapeType2ExpDim{
        {LibUtilities::ePoint, 0},       {LibUtilities::eSegment, 1},
        {LibUtilities::eTriangle, 2},    {LibUtilities::eQuadrilateral, 2},
        {LibUtilities::eTetrahedron, 3}, {LibUtilities::ePyramid, 3},
        {LibUtilities::ePrism, 3},       {LibUtilities::eHexahedron, 3}};

    map<ElementType_t, int> m_elemType2Order{
        {NODE, 0},       {BAR_2, 1},    {BAR_3, 2},    {BAR_4, 3},
        {BAR_5, 4},      {TRI_3, 1},    {TRI_6, 2},    {TRI_9, 3},
        {TRI_10, 3},     {TRI_12, 4},   {TRI_15, 4},   {QUAD_4, 1},
        {QUAD_8, 2},     {QUAD_9, 2},   {QUAD_12, 3},  {QUAD_16, 3},
        {QUAD_P4_16, 4}, {QUAD_25, 4},  {TETRA_4, 1},  {TETRA_10, 2},
        {TETRA_16, 3},   {TETRA_20, 3}, {TETRA_22, 4}, {TETRA_34, 4},
        {TETRA_35, 4},   {PYRA_5, 1},   {PYRA_13, 2},  {PYRA_14, 2},
        {PYRA_21, 3},    {PYRA_29, 3},  {PYRA_30, 3},  {PYRA_P4_29, 4},
        {PYRA_50, 4},    {PYRA_55, 4},  {PENTA_6, 1},  {PENTA_15, 2},
        {PENTA_18, 2},   {PENTA_24, 3}, {PENTA_38, 3}, {PENTA_40, 3},
        {PENTA_33, 4},   {PENTA_66, 4}, {PENTA_75, 4}, {HEXA_8, 1},
        {HEXA_20, 2},    {HEXA_27, 2},  {HEXA_32, 3},  {HEXA_56, 3},
        {HEXA_64, 3},    {HEXA_44, 4},  {HEXA_98, 4},  {HEXA_125, 4},
    };

    // summary: true for (QUAD_9, PYRA_14, PENTA_18, HEXA_27, TRI_10, QUAD_16,
    //                    TETRA_20, PYRA_29, PYRA_30, PENTA_38, PENTA_40,
    //                    HEXA_56, HEXA_64, TRI_15, QUAD_25, TETRA_34,
    //                    TETRA_35, PYRA_50, PYRA_55, PENTA_66, PENTA_75,
    //                    HEXA_98, HEXA_125)
    //          false for the rest
    map<ElementType_t, int> m_elemType2FaceNodes{
        // {ElementTypeNull, },   {ElementTypeUserDefined, },
        {NODE, false},    {BAR_2, false},      {BAR_3, false},
        {BAR_4, false},   {BAR_5, false},      {TRI_3, false},
        {TRI_6, false},   {TRI_9, false},      {TRI_10, true},
        {TRI_12, false},  {TRI_15, true},      {QUAD_4, false},
        {QUAD_8, false},  {QUAD_9, true},      {QUAD_12, false},
        {QUAD_16, true},  {QUAD_P4_16, false}, {QUAD_25, true},
        {TETRA_4, false}, {TETRA_10, false},   {TETRA_16, false},
        {TETRA_20, true}, {TETRA_22, false},   {TETRA_34, true},
        {TETRA_35, true}, {PYRA_5, false},     {PYRA_13, false},
        {PYRA_14, true},  {PYRA_21, false},    {PYRA_29, true},
        {PYRA_30, true},  {PYRA_P4_29, false}, {PYRA_50, true},
        {PYRA_55, true},  {PENTA_6, false},    {PENTA_15, false},
        {PENTA_18, true}, {PENTA_24, false},   {PENTA_38, true},
        {PENTA_40, true}, {PENTA_33, false},   {PENTA_66, true},
        {PENTA_75, true}, {HEXA_8, false},     {HEXA_20, false},
        {HEXA_27, true},  {HEXA_32, false},    {HEXA_56, true},
        {HEXA_64, true},  {HEXA_44, false},    {HEXA_98, true},
        {HEXA_125, true},
    };

    // summary: true for (HEXA_27, PYRA_30, PENTA_40, HEXA_64, TETRA_35,
    //                    PYRA_55, PENTA_75, HEXA_125)
    //          false otherwise
    map<ElementType_t, int> m_elemType2VolNodes{
        {NODE, false},     {BAR_2, false},      {BAR_3, false},
        {BAR_4, false},    {BAR_5, false},      {TRI_3, false},
        {TRI_6, false},    {TRI_9, false},      {TRI_12, false},
        {TRI_15, false},   {TRI_10, false},     {QUAD_4, false},
        {QUAD_8, false},   {QUAD_9, false},     {QUAD_12, false},
        {QUAD_16, false},  {QUAD_P4_16, false}, {QUAD_25, false},
        {TETRA_4, false},  {TETRA_10, false},   {TETRA_16, false},
        {TETRA_20, false}, {TETRA_22, false},   {TETRA_34, false},
        {TETRA_35, true},  {PYRA_5, false},     {PYRA_13, false},
        {PYRA_14, false},  {PYRA_21, false},    {PYRA_29, false},
        {PYRA_30, true},   {PYRA_P4_29, false}, {PYRA_50, false},
        {PYRA_55, true},   {PENTA_6, false},    {PENTA_15, false},
        {PENTA_18, false}, {PENTA_24, false},   {PENTA_38, false},
        {PENTA_40, true},  {PENTA_33, false},   {PENTA_66, false},
        {PENTA_75, true},  {HEXA_8, false},     {HEXA_20, false},
        {HEXA_27, true},   {HEXA_32, false},    {HEXA_56, false},
        {HEXA_64, true},   {HEXA_44, false},    {HEXA_98, false},
        {HEXA_125, true},
    };

    map<ElementType_t, ElementType_t> m_pyrShieldingTetType{
        {PYRA_5, TETRA_4},      {PYRA_13, TETRA_10}, {PYRA_14, TETRA_10},
        {PYRA_21, TETRA_16},    {PYRA_29, TETRA_20}, {PYRA_30, TETRA_20},
        {PYRA_P4_29, TETRA_22}, {PYRA_50, TETRA_34}, {PYRA_55, TETRA_35}};
};
} // namespace Nektar::NekMesh

#endif
