///////////////////////////////////////////////////////////////////////////////
//
// File: FilterPython.h
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

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERPYTHON_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERPYTHON_H

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SolverUtils/Filters/Filter.h>

#include <pybind11/embed.h>

namespace Nektar::SolverUtils
{

#pragma GCC visibility push(hidden)
class FilterPython : public Filter
{
public:
    friend class MemoryManager<FilterPython>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<EquationSystem> &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterPython>::AllocateSharedPtr(
            pSession, pEquation, pParams);
        return p;
    }

    /// Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterPython(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<EquationSystem> &pEquation,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT ~FilterPython() override = default;

protected:
    void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    bool v_IsTimeDependent() override;

private:
    /// Container for Python global variables.
    py::object m_global;
    /// Python filter class to run (optionally)
    std::shared_ptr<Filter> m_pyFilter = nullptr;
};
#pragma GCC visibility pop

} // namespace Nektar::SolverUtils

#endif
