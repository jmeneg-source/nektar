//////////////////////////////////////////////////////////////////////////////
//
// File: EquationSystem.h
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
// Description: Python wrapper for EquationSystem.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_SOLVERUTILS_PYTHON_EQUATIONSYSTEM_H
#define NEKTAR_LIBRARY_SOLVERUTILS_PYTHON_EQUATIONSYSTEM_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SpatialDomains/MeshGraph.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

/**
 * @brief EquationSystem wrapper to handle virtual function calls in @c
 * EquationSystem and its subclasses.
 */
#pragma GCC visibility push(hidden)
template <typename T>
struct EquationSystemWrap : public T, public py::trampoline_self_life_support
{
    /**
     * @brief Constructor, which is identical to Filter.
     *
     * @param field  Input field.
     */
    EquationSystemWrap(LibUtilities::SessionReaderSharedPtr session,
                       SpatialDomains::MeshGraphSharedPtr graph)
        : T(session, graph)
    {
    }

    void v_InitObject(bool declareField) override
    {
        PYBIND11_OVERRIDE_NAME(void, T, "InitObject", v_InitObject,
                               declareField);
    }

    void v_DoInitialise(bool dumpInitialConditions) override
    {
        PYBIND11_OVERRIDE_NAME(void, T, "DoInitialize", v_DoInitialise,
                               dumpInitialConditions);
    }

    void v_DoSolve() override
    {
        PYBIND11_OVERRIDE_NAME(void, T, "DoSolve", v_DoSolve, );
    }

    void v_SetInitialConditions(NekDouble initialtime,
                                bool dumpInitialConditions,
                                const int domain) override
    {
        PYBIND11_OVERRIDE_NAME(void, T, "SetInitialConditions",
                               v_SetInitialConditions, initialtime,
                               dumpInitialConditions, domain);
    }

    static Array<OneD, NekDouble> EvaluateExactSolution(
        std::shared_ptr<EquationSystemWrap> eqsys, unsigned int field,
        const NekDouble time)
    {
        py::gil_scoped_acquire gil;
        py::function override =
            py::get_override(eqsys.get(), "EvaluateExactSolution");

        if (override)
        {
            auto obj = override(field, time);
            return py::cast<Array<OneD, NekDouble>>(obj);
        }

        Array<OneD, NekDouble> outfield(
            eqsys->UpdateFields()[field]->GetNpoints());
        eqsys->EquationSystemWrap<T>::v_EvaluateExactSolution(field, outfield,
                                                              time);
        return outfield;
    }

    static NekDouble LinfError(std::shared_ptr<EquationSystemWrap> eqsys,
                               unsigned int field)
    {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(eqsys.get(), "LinfError");

        if (override)
        {
            return py::cast<NekDouble>(override(field));
        }

        return eqsys->EquationSystemWrap<T>::v_LinfError(field);
    }

    static NekDouble L2Error(std::shared_ptr<EquationSystemWrap> eqsys,
                             unsigned int field)
    {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(eqsys.get(), "L2Error");

        if (override)
        {
            return py::cast<NekDouble>(override(field));
        }

        return eqsys->EquationSystemWrap<T>::v_L2Error(field);
    }

    using T::v_EvaluateExactSolution;
    using T::v_L2Error;
    using T::v_LinfError;
};
#pragma GCC visibility pop

template <class T> struct EquationSystemPublic : public T
{
public:
    using T::v_DoInitialise;
    using T::v_DoSolve;
    using T::v_EvaluateExactSolution;
    using T::v_InitObject;
    using T::v_SetInitialConditions;
};

#endif
