///////////////////////////////////////////////////////////////////////////////
//
// File: NekPyConfig.hpp
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
// Description: NekPy configuration to include boost headers and define
// commonly-used macros.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_LIBUTILITIES_PYTHON_NEKPYCONFIG_HPP
#define NEKTAR_LIBRARY_LIBUTILITIES_PYTHON_NEKPYCONFIG_HPP

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <memory>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <LibUtilities/BasicUtils/CppCommandLine.hpp>

namespace py = pybind11;

// Define some common STL opaque types
PYBIND11_MAKE_OPAQUE(std::vector<unsigned int>)

#define SIZENAME(s) SIZE_##s
#define NEKPY_WRAP_ENUM(MOD, ENUMNAME, MAPNAME)                                \
    {                                                                          \
        py::enum_<ENUMNAME> tmp(MOD, #ENUMNAME);                               \
        for (int a = 0; a < (int)SIZENAME(ENUMNAME); ++a)                      \
        {                                                                      \
            tmp.value(MAPNAME[a], (ENUMNAME)a);                                \
        }                                                                      \
        tmp.export_values();                                                   \
    }
#define NEKPY_WRAP_ENUM_STRING(MOD, ENUMNAME, MAPNAME)                         \
    {                                                                          \
        py::enum_<ENUMNAME> tmp(MOD, #ENUMNAME);                               \
        for (int a = 0; a < (int)SIZENAME(ENUMNAME); ++a)                      \
        {                                                                      \
            tmp.value(MAPNAME[a].c_str(), (ENUMNAME)a);                        \
        }                                                                      \
        tmp.export_values();                                                   \
    }
#define NEKPY_WRAP_ENUM_STRING_DOCS(MOD, ENUMNAME, MAPNAME, DOCSTRING)         \
    {                                                                          \
        py::enum_<ENUMNAME> tmp(MOD, #ENUMNAME);                               \
        for (int a = 0; a < (int)SIZENAME(ENUMNAME); ++a)                      \
        {                                                                      \
            tmp.value(MAPNAME[a].c_str(), (ENUMNAME)a);                        \
        }                                                                      \
        tmp.export_values();                                                   \
        PyTypeObject *pto = reinterpret_cast<PyTypeObject *>(tmp.ptr());       \
        PyDict_SetItemString(pto->tp_dict, "__doc__",                          \
                             PyUnicode_FromString(DOCSTRING));                 \
    }

/**
 * @brief Helper structure to construct C++ command line `argc` and `argv`
 * variables from a Python list.
 */
struct PyCppCommandLine : public Nektar::LibUtilities::CppCommandLine
{
    /**
     * @brief Constructor.
     *
     * @param py_argv   List of command line arguments from Python.
     */
    PyCppCommandLine(py::list &py_argv) : Nektar::LibUtilities::CppCommandLine()
    {
        int argc = py::len(py_argv);
        std::vector<std::string> argv(argc);

        for (int i = 0; i < argc; ++i)
        {
            argv[i] = py::cast<std::string>(py_argv[i]);
        }

        Setup(argv);
    }
};

#endif
