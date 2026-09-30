///////////////////////////////////////////////////////////////////////////////
//
// File: Element.cpp
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
// Description: Python wrapper for Element.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <NekMesh/MeshElements/Element.h>

// For the inclusion of opaque types.
#include <NekMesh/Python/NekMesh.h>

using namespace Nektar;
using namespace Nektar::NekMesh;

ElmtConfig *ElmtConfig_Init(LibUtilities::ShapeType shapeType,
                            unsigned int order, bool faceNodes, bool volNodes,
                            bool reorient,
                            LibUtilities::PointsType edgeNodeType,
                            LibUtilities::PointsType faceNodeType)
{
    // Create a new ElmtConfig object. We return a pointer deliberately, since
    // boost::python has to store either a pointer or a shared_ptr to this
    // object.
    return new ElmtConfig(shapeType, order, faceNodes, volNodes, reorient,
                          edgeNodeType, faceNodeType);
}

std::shared_ptr<Element> Element_Create(ElmtConfig *conf, py::list &nodes,
                                        py::list &tags)
{
    int numNodes = py::len(nodes), numTags = py::len(tags);

    std::vector<NodeSharedPtr> nodes_(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
        nodes_[i] = py::cast<NodeSharedPtr>(nodes[i]);
    }

    std::vector<int> tags_(numTags);
    for (int i = 0; i < numTags; ++i)
    {
        tags_[i] = py::cast<int>(tags[i]);
    }

    return GetElementFactory().CreateInstance(conf->m_e, *conf, nodes_, tags_);
}

void export_Element(py::module &m)
{
    // Export element configuration struct
    py::class_<ElmtConfig>(m, "ElmtConfig")
        .def(py::init<>(&ElmtConfig_Init), py::arg("shapeType"),
             py::arg("order"), py::arg("faceNodes"), py::arg("volNodes"),
             py::arg("reorient")     = true,
             py::arg("edgeNodeType") = LibUtilities::eNoPointsType,
             py::arg("faceNodeType") = LibUtilities::eNoPointsType);

    // Export element base class
    py::class_<Element, std::shared_ptr<Element>>(m, "Element")
        .def("GetId", &Element::GetId)
        .def("GetDim", &Element::GetDim)
        .def("GetShapeType", &Element::GetShapeType)
        .def("GetTag", &Element::GetTag)

        // Factory methods
        .def_static("Create", &Element_Create);

    py::bind_vector<std::vector<ElementSharedPtr>>(m, "ElementVector");
}
