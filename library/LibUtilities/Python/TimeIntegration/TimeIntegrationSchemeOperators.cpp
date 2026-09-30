//////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeOperators.cpp
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
// Description: Python wrapper for TimeIntegrationScheme.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>

using namespace Nektar;
using namespace Nektar::LibUtilities;

/**
 * @brief Helper class for holding a reference to a Python function to act as a
 * wrapper for TimeIntegrationScheme::FunctorType1.
 *
 * This wrapper is used for TimeIntegrationSchemeOperators::DefineOdeRhs and
 * TimeIntegrationSchemeOperators::DefineProjection.
 */
#pragma GCC visibility push(hidden)
struct CallbackHolderT1
{
    /// Default constructor
    CallbackHolderT1(py::object cb) : m_cb(cb)
    {
    }

    /**
     * @brief C++ callback function to invoke Python function stored in #m_cb.
     */
    void call(TimeIntegrationSchemeOperators::InArrayType &in,
              TimeIntegrationSchemeOperators::OutArrayType &out,
              const NekDouble time)
    {
        py::object ret = m_cb(in, time);

        py::list outList = py::cast<py::list>(ret);

        for (std::size_t i = 0; i < py::len(outList); ++i)
        {
            out[i] = py::cast<Array<OneD, NekDouble>>(outList[i]);
        }
    }

private:
    /// Callback defined in Python code.
    py::object m_cb;
};
#pragma GCC visibility pop

/**
 * @brief Helper class for holding a reference to a Python function to act as a
 * wrapper for TimeIntegrationScheme::FunctorType1.
 *
 * This wrapper is used for TimeIntegrationSchemeOperators::DefineImplicitSolve.
 */
#pragma GCC visibility push(hidden)
struct CallbackHolderT2
{
    /// Default constructor
    CallbackHolderT2(py::object cb) : m_cb(cb)
    {
    }

    /**
     * @brief C++ callback function to invoke Python function stored in #m_cb.
     */
    void call(TimeIntegrationSchemeOperators::InArrayType &in,
              TimeIntegrationSchemeOperators::OutArrayType &out,
              const NekDouble time, const NekDouble lambda)
    {
        py::object ret = m_cb(in, time, lambda);

        py::list outList = py::cast<py::list>(ret);

        for (std::size_t i = 0; i < py::len(outList); ++i)
        {
            out[i] = py::cast<Array<OneD, NekDouble>>(outList[i]);
        }
    }

private:
    /// Callback defined in Python code.
    py::object m_cb;
};
#pragma GCC visibility pop

void TimeIntegrationSchemeOperators_DefineOdeRhs(
    TimeIntegrationSchemeOperatorsSharedPtr op, py::object callback)
{
    CallbackHolderT1 *cb = new CallbackHolderT1(callback);

    op->DefineOdeRhs(&CallbackHolderT1::call, cb);
}

void TimeIntegrationSchemeOperators_DefineProjection(
    TimeIntegrationSchemeOperatorsSharedPtr op, py::object callback)
{
    CallbackHolderT1 *cb = new CallbackHolderT1(callback);

    op->DefineProjection(&CallbackHolderT1::call, cb);
}

void TimeIntegrationSchemeOperators_DefineImplicitSolve(
    TimeIntegrationSchemeOperatorsSharedPtr op, py::object callback)
{
    CallbackHolderT2 *cb = new CallbackHolderT2(callback);

    op->DefineImplicitSolve(&CallbackHolderT2::call, cb);
}

void export_TimeIntegrationSchemeOperators(py::module &m)
{
    py::class_<TimeIntegrationSchemeOperators,
               std::shared_ptr<TimeIntegrationSchemeOperators>>(
        m, "TimeIntegrationSchemeOperators")
        .def("DefineOdeRhs", &TimeIntegrationSchemeOperators_DefineOdeRhs)
        .def("DefineProjection",
             &TimeIntegrationSchemeOperators_DefineProjection)
        .def("DefineImplicitSolve",
             &TimeIntegrationSchemeOperators_DefineImplicitSolve);
}
