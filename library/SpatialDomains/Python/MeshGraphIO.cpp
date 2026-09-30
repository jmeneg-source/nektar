////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraphIO.cpp
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
//  Description: Python wrapper for MeshGraphIO.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SpatialDomains/Python/SpatialDomains.h>

#include <SpatialDomains/MeshGraphIO.h>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

/*
 * @brief Lightweight wrapper around MeshGraph::Read to avoid wrapping
 * DomainRange struct.
 */
MeshGraphSharedPtr MeshGraphIO_Read(
    const LibUtilities::SessionReaderSharedPtr &session)
{
    return MeshGraphIO::Read(session);
}

MeshGraphIOSharedPtr MeshGraphIO_Create(std::string ioType)
{
    return GetMeshGraphIOFactory().CreateInstance(ioType);
}

void export_MeshGraphIO(py::module &m)
{
    py::class_<MeshGraphIO, std::shared_ptr<MeshGraphIO>>(m, "MeshGraphIO")

        .def("Write", &MeshGraphIO::WriteGeometry, py::arg("outfile"),
             py::arg("defaultExp") = false,
             py::arg("metadata")   = LibUtilities::NullFieldMetaDataMap)

        .def("SetMeshGraph", &MeshGraphIO::SetMeshGraph)

        .def_static("Read", MeshGraphIO_Read)
        .def_static("Create", MeshGraphIO_Create);
}
