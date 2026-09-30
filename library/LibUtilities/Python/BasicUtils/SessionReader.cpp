///////////////////////////////////////////////////////////////////////////////
//
// File: SessionReader.cpp
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
// Description: Python wrapper for SessionReader.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <pybind11/stl.h>

using namespace Nektar::LibUtilities;

#ifdef NEKTAR_USE_MPI
#include <LibUtilities/Communication/CommMpi.h>
CommSharedPtr MPICOMM = CommSharedPtr();
#endif

/*
 * @brief Thin wrapper around SessionReader to provide a nicer Pythonic
 * interface.
 *
 * This allows us to do, for example
 *
 *     session = SessionReader.CreateInstance(sys.argv)
 *
 * which is more natural in Python.
 */

SessionReaderSharedPtr SessionReader_CreateInstance(py::list &ns)
{
    PyCppCommandLine cpp(ns);

    int argc    = cpp.GetArgc();
    char **argv = cpp.GetArgv();

#ifdef NEKTAR_USE_MPI
    // In the case we're using MPI, it may already have been initialised. So to
    // handle this, we'll construct our own CommMpi object and pass through to
    // the SessionReader. This will persist indefinitely, or at least until the
    // library is unloaded by Python.

    if (!MPICOMM)
    {
        MPICOMM = GetCommFactory().CreateInstance("ParallelMPI", argc, argv);
    }

    std::vector<std::string> filenames(argc - 1);
    for (int i = 1; i < argc; ++i)
    {
        filenames[i - 1] = std::string(argv[i]);
    }

    // Create session reader.
    SessionReaderSharedPtr sr =
        SessionReader::CreateInstance(argc, argv, filenames, MPICOMM);
#else
    // Create session reader.
    SessionReaderSharedPtr sr = SessionReader::CreateInstance(argc, argv);
#endif

    return sr;
}

/**
 * @brief SessionReader exports.
 */

void export_SessionReader(py::module &m)
{
    py::class_<SessionReader, std::shared_ptr<SessionReader>>(m,
                                                              "SessionReader")

        .def_static("CreateInstance", SessionReader_CreateInstance)

        .def("GetSessionName", &SessionReader::GetSessionName,
             py::return_value_policy::copy)

        .def("InitSession", &SessionReader::InitSession,
             py::arg("filenames") = py::list())

        .def("Finalise", &SessionReader::Finalise)

        .def("DefinesParameter", &SessionReader::DefinesParameter)
        .def("GetParameter", &SessionReader::GetParameter,
             py::return_value_policy::copy)
        .def("GetParameters", &SessionReader::GetParameters)

        .def("SetParameter", py::overload_cast<const std::string &, int &>(
                                 &SessionReader::SetParameter))
        .def("SetParameter", py::overload_cast<const std::string &, double &>(
                                 &SessionReader::SetParameter))
        .def("SetParameter", py::overload_cast<const std::string &, size_t &>(
                                 &SessionReader::SetParameter))

        .def("DefinesSolverInfo", &SessionReader::DefinesSolverInfo)
        .def("GetSolverInfo", &SessionReader::GetSolverInfo,
             py::return_value_policy::copy)
        .def("SetSolverInfo", &SessionReader::SetSolverInfo)

        .def("GetVariable", &SessionReader::GetVariable,
             py::return_value_policy::copy)
        .def("GetVariables", &SessionReader::GetVariables)

        .def("GetFunction",
             py::overload_cast<const std::string &, const std::string &,
                               const int>(&SessionReader::GetFunction,
                                          py::const_),
             py::arg("name"), py::arg("var"), py::arg("domain") = 0)
        .def("GetFunction",
             py::overload_cast<const std::string &, const unsigned int &,
                               const int>(&SessionReader::GetFunction,
                                          py::const_),
             py::arg("name"), py::arg("var"), py::arg("domain") = 0)

        .def("GetComm", &SessionReader::GetComm)

        .def("GetSharedFilesystem", &SessionReader::GetSharedFilesystem);
}
