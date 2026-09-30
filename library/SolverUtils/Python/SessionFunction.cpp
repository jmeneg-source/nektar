//////////////////////////////////////////////////////////////////////////////
//
// File: SessionFunction.cpp
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
// Description: Python wrapper for SessionFunction.
//
///////////////////////////////////////////////////////////////////////////////

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SolverUtils/Core/SessionFunction.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

std::shared_ptr<SessionFunction> SessionFunction_Init(
    const LibUtilities::SessionReaderSharedPtr session,
    const MultiRegions::ExpListSharedPtr field, std::string functionName)
{
    return MemoryManager<SessionFunction>::AllocateSharedPtr(session, field,
                                                             functionName);
}

std::string SessionFunction_Describe(std::shared_ptr<SessionFunction> func,
                                     std::string fieldName)
{
    return func->Describe(fieldName);
}

Array<OneD, MultiRegions::ExpListSharedPtr> SessionFunction_Evaluate_ExpList(
    std::shared_ptr<SessionFunction> func, std::vector<std::string> fieldNames,
    Array<OneD, MultiRegions::ExpListSharedPtr> explists, NekDouble time)
{
    func->Evaluate(fieldNames, explists, time);
    return explists;
}

Array<OneD, MultiRegions::ExpListSharedPtr> SessionFunction_Evaluate(
    std::shared_ptr<SessionFunction> func, py::list fieldNames,
    Array<OneD, MultiRegions::ExpListSharedPtr> explists, NekDouble time)
{
    std::vector<std::string> fields(py::len(fieldNames));

    for (std::size_t i = 0; i < py::len(fieldNames); ++i)
    {
        fields[i] = py::cast<std::string>(fieldNames[i]);
    }

    if (explists.size() == 0)
    {
        throw new ErrorUtil::NekError("List of fields is empty");
    }

    // Check first entry in list to guess which version of this function needs
    // to be called.
    return SessionFunction_Evaluate_ExpList(func, fields, explists, time);
}

Array<OneD, MultiRegions::ExpListSharedPtr> SessionFunction_Evaluate1(
    std::shared_ptr<SessionFunction> func, py::list fieldNames,
    Array<OneD, MultiRegions::ExpListSharedPtr> explists)
{
    return SessionFunction_Evaluate(func, fieldNames, explists, 0.0);
}

void export_SessionFunction(py::module &m)
{
    py::class_<SessionFunction, std::shared_ptr<SessionFunction>>(
        m, "SessionFunction")
        .def(py::init<>(&SessionFunction_Init))
        .def("Describe", &SessionFunction_Describe)
        .def("Evaluate", &SessionFunction_Evaluate)
        .def("Evaluate", &SessionFunction_Evaluate1);
}
