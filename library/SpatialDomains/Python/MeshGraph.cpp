////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph.cpp
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
//  Description: Python wrapper for MeshGraph.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <SpatialDomains/Movement/Movement.h>
#include <SpatialDomains/Python/SpatialDomains.h>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

/*
 * @brief Simple wrapper to build Composite objects from lists of
 * Geometry objects.
 */
CompositeSharedPtr Composite_Init(py::list geometries)
{
    CompositeSharedPtr composite = std::make_shared<Composite>();
    composite->m_geomVec.clear();
    for (int i = 0; i < py::len(geometries); i++)
    {
        composite->m_geomVec.emplace_back(py::cast<Geometry *>(geometries[i]));
    }
    return composite;
}

template <typename T>
void MeshGraph_AddGeom(MeshGraphSharedPtr graph, unique_ptr_objpool<T> geom)
{
    auto id = geom->GetGlobalID();
    graph->AddGeom(id, std::move(geom));
}

template <typename T>
void MeshGraph_GeomMapView(py::module_ &m, const std::string &name)
{
    using MapView = GeomMapView<T>;

    py::class_<MapView>(m, name.c_str())
        .def("__getitem__", &MapView::at, py::return_value_policy::reference)
        .def("__len__", &MapView::size)
        .def(
            "__iter__",
            [](const MapView &self) {
                return py::make_key_iterator(self.begin(), self.end());
            },
            py::keep_alive<0, 1>()) // keep object alive while iterating
        .def(
            "items",
            [](const MapView &self) {
                return py::make_iterator(self.begin(), self.end());
            },
            py::keep_alive<0, 1>())
        .def("values",
             [](const MapView &self) {
                 py::list result;
                 for (const auto &[_, ptr] : self)
                 {
                     result.append(ptr);
                 }
                 return result;
             })
        .def("keys",
             [](const MapView &self) {
                 py::list result;
                 for (const auto &[id, _] : self)
                 {
                     result.append(id);
                 }
                 return result;
             })
        .def(
            "get",
            [](const MapView &self, int id) -> T * {
                auto it = self.find(id);
                return (*it).second;
            },
            py::return_value_policy::reference);
}

/**
 * @brief MeshGraph exports.
 */
void export_MeshGraph(py::module &m)
{
    py::bind_map<LibUtilities::FieldMetaDataMap>(m, "FieldMetaDataMap");
    py::bind_vector<std::vector<Geometry *>>(m, "GeometryList");

    py::class_<Composite, std::shared_ptr<Composite>>(m, "Composite")
        .def(py::init<>())
        .def(py::init<>(&Composite_Init))
        .def_readwrite("geometries", &Composite::m_geomVec);

    py::bind_map<CurveMap>(m, "CurveMap");
    py::bind_map<CompositeMap>(m, "CompositeMap");
    py::bind_map<std::map<int, CompositeMap>>(m, "DomainMap");

    MeshGraph_GeomMapView<PointGeom>(m, "PointGeomView");
    MeshGraph_GeomMapView<SegGeom>(m, "SegGeomView");
    MeshGraph_GeomMapView<TriGeom>(m, "TriGeomView");
    MeshGraph_GeomMapView<QuadGeom>(m, "QuadGeomView");
    MeshGraph_GeomMapView<TetGeom>(m, "TetGeomView");
    MeshGraph_GeomMapView<PyrGeom>(m, "PyrGeomView");
    MeshGraph_GeomMapView<PrismGeom>(m, "PrismGeomView");
    MeshGraph_GeomMapView<HexGeom>(m, "HexGeomView");

    py::class_<MeshGraph, std::shared_ptr<MeshGraph>>(m, "MeshGraph")
        .def(py::init<>())

        .def("Empty", &MeshGraph::Empty)

        .def("GetMeshDimension", &MeshGraph::GetMeshDimension)
        .def("GetSpaceDimension", &MeshGraph::GetSpaceDimension)
        .def("SetMeshDimension", &MeshGraph::SetMeshDimension)
        .def("SetSpaceDimension", &MeshGraph::SetSpaceDimension)

        .def_property_readonly("points", &MeshGraph::GetGeomMap<PointGeom>,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("segments", &MeshGraph::GetGeomMap<SegGeom>,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("quads", &MeshGraph::GetGeomMap<QuadGeom>,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("tris", &MeshGraph::GetGeomMap<TriGeom>,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("tets", &MeshGraph::GetGeomMap<TetGeom>,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("pyrs", &MeshGraph::GetGeomMap<PyrGeom>,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("prisms", &MeshGraph::GetGeomMap<PrismGeom>,
                               py::return_value_policy::reference_internal)
        .def_property_readonly("hexes", &MeshGraph::GetGeomMap<HexGeom>,
                               py::return_value_policy::reference_internal)

        .def("AddGeom", &MeshGraph_AddGeom<PointGeom>)
        .def("AddGeom", &MeshGraph_AddGeom<SegGeom>)
        .def("AddGeom", &MeshGraph_AddGeom<TriGeom>)
        .def("AddGeom", &MeshGraph_AddGeom<QuadGeom>)

        .def("GetVertex", &MeshGraph::GetPointGeom,
             py::return_value_policy::reference_internal)
        .def("GetPointGeom", &MeshGraph::GetPointGeom,
             py::return_value_policy::reference_internal)
        .def("GetSegGeom", &MeshGraph::GetSegGeom,
             py::return_value_policy::reference_internal)
        .def("GetTriGeom", &MeshGraph::GetTriGeom,
             py::return_value_policy::reference_internal)
        .def("GetQuadGeom", &MeshGraph::GetQuadGeom,
             py::return_value_policy::reference_internal)
        .def("GetHexGeom", &MeshGraph::GetHexGeom,
             py::return_value_policy::reference_internal)
        .def("GetPrismGeom", &MeshGraph::GetPrismGeom,
             py::return_value_policy::reference_internal)
        .def("GetTetGeom", &MeshGraph::GetTetGeom,
             py::return_value_policy::reference_internal)
        .def("GetPyrGeom", &MeshGraph::GetPyrGeom,
             py::return_value_policy::reference_internal)

        //.def("AddGeom", &MeshGraph::AddGeom<PointGeom>)
        .def("GetCurvedEdges", &MeshGraph::GetCurvedEdges,
             py::return_value_policy::reference_internal)
        .def("GetCurvedFaces", &MeshGraph::GetCurvedFaces,
             py::return_value_policy::reference_internal)
        .def("GetComposites", &MeshGraph::GetComposites,
             py::return_value_policy::reference_internal)
        .def<std::map<int, CompositeMap> &(MeshGraph::*)()>(
            "GetDomain", &MeshGraph::GetDomain,
            py::return_value_policy::reference_internal)

        .def("GetMovement", &MeshGraph::GetMovement,
             py::return_value_policy::reference_internal)

        .def("GetNumElements", &MeshGraph::GetNumElements)

        .def("SetExpansionInfosToEvenlySpacedPoints",
             &MeshGraph::SetExpansionInfoToEvenlySpacedPoints)
        .def("SetExpansionInfosToPolyOrder",
             &MeshGraph::SetExpansionInfoToNumModes)
        .def("SetExpansionInfosToPointOrder",
             &MeshGraph::SetExpansionInfoToPointOrder);
}
