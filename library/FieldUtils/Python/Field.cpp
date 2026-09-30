///////////////////////////////////////////////////////////////////////////////
//
// File: Field.cpp
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
// Description: Python wrapper for Field class.
//
///////////////////////////////////////////////////////////////////////////////

#include <FieldUtils/Field.hpp>
#include <FieldUtils/FieldConvertComm.hpp>
#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::FieldUtils;

// Called at the start of every loop over nparts.
// Args: FieldSharedPtr "f",
// Python's sys.argv stored as a Python list called "py_argv",
// and the partition number stored as an integer called "part".
// Function: clears data stored in f, defines f->m_partComm.
void NewPartition(FieldSharedPtr f, py::list &py_argv, int part)
{
    std::cout << std::endl << "Processing partition: " << part << std::endl;
    f->ClearField();

    PyCppCommandLine cpp(py_argv);

    f->m_partComm = std::shared_ptr<FieldConvertComm>(
        new FieldConvertComm(cpp.GetArgc(), cpp.GetArgv(), f->m_nParts, part));
}

// Wrapper around the Field constructor
// Args: FieldConvert command line arguments.
// Function: constructs a FieldSharedPtr "f", defines f->m_defComm if
// nparts specified, and stores a line arguments in f->m_vm.
// Returns: f
FieldSharedPtr Field_Init(py::list &argv, int nparts = 0, int output_points = 0,
                          int output_points_hom_z = 0, bool error = false,
                          bool force_output = false, bool no_equispaced = false,
                          int npz = 0, std::string onlyshape = "",
                          int part_only = 0, int part_only_overlapping = 0,
                          bool useSessionVariables = false,
                          bool useSessionExpansion = false,
                          bool verbose = false, std::string domain = "")
{
    // Construct shared pointer to a Field object.
    std::shared_ptr<Field> f = MemoryManager<Field>::AllocateSharedPtr();

    if (py::len(argv) > 0)
    {
        // Get argc and argv from the Python command line.
        PyCppCommandLine cpp(argv);

        // Define the MPI Communicator.
        f->m_comm = LibUtilities::GetCommFactory().CreateInstance(
            "Serial", cpp.GetArgc(), cpp.GetArgv());

        if (nparts)
        {
            f->m_vm.insert(
                std::make_pair("nparts", po::variable_value(nparts, false)));
            if (nparts > 1)
            {
                f->m_nParts  = nparts;
                f->m_defComm = f->m_comm;
            }
        }
    }

    // Populate m_vm.
    if (output_points)
    {
        f->m_vm.insert(std::make_pair(
            "output-points", po::variable_value(output_points, false)));
    }

    if (output_points_hom_z)
    {
        f->m_vm.insert(
            std::make_pair("output-points-hom-z",
                           po::variable_value(output_points_hom_z, false)));
    }

    if (error)
    {
        f->m_vm.insert(std::make_pair("error", po::variable_value()));
    }

    if (force_output)
    {
        f->m_vm.insert(std::make_pair("force-output", po::variable_value()));
    }

    if (domain.size())
    {
        NEKERROR(ErrorUtil::efatal,
                 "The doamin option in field is deprecated. Please use "
                 "the xml option range=\"xmax,xmin,ymax,ymin\", \n\t i.e."
                 "InputModule.Create(\"xml\",  field, \"myfile.xml\", "
                 "range=\"-1,1,-1,1\").Run()");
    }

    if (no_equispaced)
    {
        f->m_vm.insert(std::make_pair("no-equispaced", po::variable_value()));
    }

    if (npz)
    {
        f->m_vm.insert(std::make_pair("npz", po::variable_value(npz, false)));
    }

    if (onlyshape.size())
    {
        f->m_vm.insert(
            std::make_pair("onlyshape", po::variable_value(onlyshape, false)));
    }

    if (part_only)
    {
        f->m_vm.insert(
            std::make_pair("part-only", po::variable_value(part_only, false)));
    }

    if (part_only_overlapping)
    {
        f->m_vm.insert(
            std::make_pair("part-only-overlapping",
                           po::variable_value(part_only_overlapping, false)));
    }

    if (useSessionVariables)
    {
        f->m_vm.insert(
            std::make_pair("use_session_variables", po::variable_value()));
    }

    if (useSessionExpansion)
    {
        f->m_vm.insert(
            std::make_pair("use_session_expansion", po::variable_value()));
    }

    if (verbose)
    {
        f->m_vm.insert(std::make_pair("verbose", po::variable_value()));
    }

    return f;
}

// Get the i-th Pts field
const Array<OneD, const NekDouble> Field_GetPts(FieldSharedPtr f, const int i)
{
    return f->m_fieldPts->GetPts(i);
}

// Set the i-th Pts field
// TODO: The reference for inarray is unavailable in Python if it's not a const
// [work] Array<OneD, NekDouble> inarray
// [fail] Array<OneD, NekDouble> &inarray
// [work] const Array<OneD, const NekDouble> &inarray
void Field_SetPts(FieldSharedPtr f, const int i,
                  const Array<OneD, const NekDouble> &inarray)
{
    f->m_fieldPts->SetPts(i, inarray);
}

void export_Field(py::module &m)
{
    py::class_<Field, std::shared_ptr<Field>>(m, "Field")
        .def(py::init<>(&Field_Init), py::arg("argv") = py::list(),
             py::arg("nparts") = 0, py::arg("output_points") = 0,
             py::arg("output_points_hom_z") = 0, py::arg("error") = false,
             py::arg("force_output") = false, py::arg("no_equispaced") = false,
             py::arg("npz") = 0, py::arg("onlyshape") = "",
             py::arg("part_only") = 0, py::arg("part_only_overlapping") = 0,
             py::arg("use_session_variables") = false,
             py::arg("use_session_expansion") = false,
             py::arg("verbose") = false, py::arg("domain") = "")

        .def("GetPts", &Field_GetPts)
        .def("SetPts", &Field_SetPts)
        .def("ClearField", &Field::ClearField)
        .def("NewPartition", &NewPartition)
        .def("ReadFieldDefs", &Field::ReadFieldDefs)

        .def("SetupFromExpList", &Field::SetupFromExpList)

        .def_readwrite("graph", &Field::m_graph)
        .def_readwrite("session", &Field::m_session)
        .def_readwrite("verbose", &Field::m_verbose)
        .def_readwrite("comm", &Field::m_comm);
}
