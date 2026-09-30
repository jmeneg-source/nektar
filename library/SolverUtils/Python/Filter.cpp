///////////////////////////////////////////////////////////////////////////////
//
// File: Filter.cpp
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
// Description: Python wrapper for Filter.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <MultiRegions/ContField.h>

#include <SolverUtils/Filters/Filter.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

/**
 * @brief Filter wrapper to handle virtual function calls in @c Filter and its
 * subclasses.
 */
#pragma GCC visibility push(hidden)
struct FilterWrap : public Filter, public py::trampoline_self_life_support
{
    /**
     * @brief Constructor, which is identical to Filter.
     *
     * @param field  Input field.
     */
    FilterWrap(LibUtilities::SessionReaderSharedPtr session,
               std::shared_ptr<EquationSystem> eqn)
        : Filter(session, eqn)
    {
    }

    void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override
    {
        PYBIND11_OVERRIDE_PURE_NAME(void, Filter, "Initialise", v_Initialise,
                                    pFields, time);
    }

    void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override
    {
        PYBIND11_OVERRIDE_PURE_NAME(void, Filter, "Update", v_Update, pFields,
                                    time);
    }

    void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override
    {
        PYBIND11_OVERRIDE_PURE_NAME(void, Filter, "Finalise", v_Finalise,
                                    pFields, time);
    }

    bool v_IsTimeDependent() override
    {
        PYBIND11_OVERRIDE_PURE_NAME(bool, Filter, "IsTimeDependent",
                                    v_IsTimeDependent, );
    }
};
#pragma GCC visibility pop

struct FilterPublic : public Filter
{
    using Filter::v_Finalise;
    using Filter::v_Initialise;
    using Filter::v_IsTimeDependent;
    using Filter::v_Update;
};

/**
 * @brief Lightweight wrapper for Filter factory creation function.
 */
FilterSharedPtr Filter_Create(py::args args, const py::kwargs &kwargs)
{
    using NekError = ErrorUtil::NekError;

    if (py::len(args) != 3)
    {
        throw NekError("Filter.Create() requires three arguments: "
                       "filter name, a SessionReader object and an "
                       "EquationSystem object.");
    }

    std::string filterName = py::cast<std::string>(args[0]);

    LibUtilities::SessionReaderSharedPtr session;
    EquationSystemSharedPtr eqsys;

    try
    {
        session = py::cast<LibUtilities::SessionReaderSharedPtr>(args[1]);
    }
    catch (...)
    {
        throw NekError("Second argument to Create() should be a SessionReader "
                       "object.");
    }

    try
    {
        eqsys = py::cast<EquationSystemSharedPtr>(args[2]);
    }
    catch (...)
    {
        throw NekError("Second argument to Create() should be a EquationSystem "
                       "object.");
    }

    // Process keyword arguments.
    Filter::ParamMap params;

    for (auto &kwarg : kwargs)
    {
        std::string arg = py::cast<std::string>(kwarg.first);
        std::string val = py::cast<std::string>(py::str(kwarg.second));
        params[arg]     = val;
    }

    std::shared_ptr<Filter> filter =
        GetFilterFactory().CreateInstance(filterName, session, eqsys, params);

    return filter;
}

void export_Filter(py::module &m)
{
    static NekFactory_Register<FilterFactory> fac(GetFilterFactory());

    // Wrapper for the Filter class. Note that since Filter contains a pure
    // virtual function, we need the FilterWrap helper class to handle this for
    // us. In the lightweight wrappers above, we therefore need to ensure we're
    // passing std::shared_ptr<Filter> as the first argument, otherwise they
    // won't accept objects constructed from Python.
    py::classh<Filter, FilterWrap>(m, "Filter")
        .def(py::init<LibUtilities::SessionReaderSharedPtr,
                      EquationSystemSharedPtr>())

        .def("Initialise", &FilterPublic::v_Initialise)
        .def("Update", &FilterPublic::v_Update)
        .def("Finalise", &FilterPublic::v_Finalise)
        .def("IsTimeDependent", &FilterPublic::v_IsTimeDependent)

        // Factory functions.
        .def_static("Create", &Filter_Create)

        .def_static("Register", [](std::string &filterName, py::object &obj) {
            fac(filterName, obj, filterName);
        });
}
