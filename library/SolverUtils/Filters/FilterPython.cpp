///////////////////////////////////////////////////////////////////////////////
//
// File: FilterPython.cpp
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
// Description: Run a Python script during runtime through Filters.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterPython.h>

#include <fstream>
#include <streambuf>
#include <string>

namespace Nektar::SolverUtils
{

/// Temporarily stolen from boost examples
std::string parse_python_exception()
{
    PyObject *type_ptr = nullptr, *value_ptr = nullptr,
             *traceback_ptr = nullptr;

    // Fetch the exception info from the Python C API
    PyErr_Fetch(&type_ptr, &value_ptr, &traceback_ptr);

    // Fallback error
    std::string ret("Unfetchable Python error");

    // If the fetch got a type pointer, parse the type into the exception string
    if (type_ptr != nullptr)
    {
        py::handle h_type(type_ptr);
        ret = py::str(h_type);
    }

    // Do the same for the exception value (the stringification of the
    // exception)
    if (value_ptr != nullptr)
    {
        py::handle h_val(value_ptr);
        ret += ": " + std::string(py::str(h_val));
    }

    // Parse lines from the traceback using the Python traceback module
    if (traceback_ptr != nullptr)
    {
        py::handle h_tb(traceback_ptr);

        // Load the traceback module and the format_tb function
        py::object tb(py::module_::import("traceback"));
        py::object fmt_tb(tb.attr("format_tb"));

        // Call format_tb to get a list of traceback strings
        py::object tb_list(fmt_tb(h_tb));

        // Join the traceback strings into a single string
        py::object tb_str(py::str("\n").attr("join")(tb_list));

        // Extract the string, check the extraction, and fallback in necessary
        ret += ": " + std::string(py::str(tb_str));
    }
    return ret;
}

// Converts a OneD array of ExpLists to a Python list.
inline py::list ArrayOneDToPyList(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    py::list expLists;

    for (int i = 0; i < pFields.size(); ++i)
    {
        expLists.append(py::cast(pFields[i]));
    }

    return expLists;
}

std::string FilterPython::className =
    GetFilterFactory().RegisterCreatorFunction("Python", FilterPython::create);

FilterPython::FilterPython(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const std::shared_ptr<EquationSystem> &pEquation,
                           const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    if (!Py_IsInitialized())
    {
        py::initialize_interpreter();
    }

    auto it = pParams.find("PythonFile");
    ASSERTL0(it != pParams.end(), "Empty parameter 'PythonFile'.");
    std::string pythonFile = it->second, pythonCode;

    // Set the global namespace.
    m_global = py::module_::import("__main__").attr("__dict__");

    try
    {
        std::ifstream t(pythonFile);
        std::string str((std::istreambuf_iterator<char>(t)),
                        std::istreambuf_iterator<char>());
        pythonCode = str;
    }
    catch (...)
    {
        ASSERTL0(false, "Error reading Python file: " + pythonFile);
    }

    // Print out the loaded session file if we're running in verbose mode.
    if (pSession->DefinesCmdLineArgument("verbose"))
    {
        std::cout << "-----------------" << std::endl;
        std::cout << "BEGIN PYTHON CODE" << std::endl;
        std::cout << "-----------------" << std::endl;
        std::cout << pythonCode << std::endl;
        std::cout << "-----------------" << std::endl;
        std::cout << "END PYTHON CODE" << std::endl;
        std::cout << "-----------------" << std::endl;
    }

    try
    {
        // Import NekPy libraries. This ensures we have all of our boost python
        // wrappers. I guess this could also be done by the calling script...
        auto nekpy    = py::module_::import("NekPy");
        auto multireg = py::module_::import("NekPy.MultiRegions");

        // Eval the python code. We can then grab functions etc from the global
        // namespace.
        py::exec(pythonCode.c_str(), m_global, m_global);
    }
    catch (py::error_already_set const &)
    {
        // Parse and output the exception
        ASSERTL0(false,
                 "Failed to load Python code: " + parse_python_exception());
    }

    // Give the option as well to just create a filter that's defined via
    // Python.
    it = pParams.find("FilterName");
    if (it != pParams.end())
    {
        std::string filterName = it->second;

        // Check filter factory to make sure we have this class.
        ASSERTL0(GetFilterFactory().ModuleExists(filterName),
                 "Unable to locate filter named '" + filterName + "'");

        try
        {
            m_pyFilter = GetFilterFactory().CreateInstance(filterName, pSession,
                                                           pEquation, pParams);
        }
        catch (py::error_already_set const &)
        {
            // Parse and output the exception
            ASSERTL0(false,
                     "Failed to load Python code: " + parse_python_exception());
        }
    }
}

void FilterPython::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    try
    {
        if (m_pyFilter)
        {
            m_pyFilter->Initialise(pFields, time);
        }
        else
        {
            m_global["filter_initialise"](ArrayOneDToPyList(pFields), time);
        }
    }
    catch (py::error_already_set const &)
    {
        std::cout << "Error in Python: " << parse_python_exception()
                  << std::endl;
    }
}

void FilterPython::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    try
    {
        if (m_pyFilter)
        {
            m_pyFilter->Update(pFields, time);
        }
        else
        {
            m_global["filter_update"](ArrayOneDToPyList(pFields), time);
        }
    }
    catch (py::error_already_set const &)
    {
        std::cout << "Error in Python: " << parse_python_exception()
                  << std::endl;
    }
}

void FilterPython::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    try
    {
        if (m_pyFilter)
        {
            m_pyFilter->Finalise(pFields, time);
        }
        else
        {
            m_global["filter_finalise"](ArrayOneDToPyList(pFields), time);
        }
    }
    catch (py::error_already_set const &)
    {
        std::cout << "Error in Python: " << parse_python_exception()
                  << std::endl;
    }
}

bool FilterPython::v_IsTimeDependent()
{
    if (m_pyFilter)
    {
        return m_pyFilter->IsTimeDependent();
    }

    return true;
}

} // namespace Nektar::SolverUtils
