///////////////////////////////////////////////////////////////////////////////
//
// File: EquationSystem.cpp
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
// Description: Python wrapper for the EquationSystem.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Python/EquationSystem.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

/**
 * @brief Dummy equation system that can be used for Python testing.
 *
 * This class simply evaluates a known function as the implementation as part of
 * the DoInitialise method.
 */
class DummyEquationSystem : public EquationSystem
{
public:
    friend class MemoryManager<DummyEquationSystem>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        EquationSystemSharedPtr p =
            MemoryManager<DummyEquationSystem>::AllocateSharedPtr(pSession,
                                                                  pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    ~DummyEquationSystem() override
    {
    }

protected:
    DummyEquationSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : EquationSystem(pSession, pGraph)
    {
    }

    void v_InitObject(bool DeclareField) override
    {
        EquationSystem::v_InitObject(DeclareField);

        for (int i = 0; i < m_fields.size(); ++i)
        {
            int nq = m_fields[0]->GetNpoints();

            Array<OneD, NekDouble> x(nq), y(nq), z(nq);
            m_fields[i]->GetCoords(x, y, z);

            for (int j = 0; j < nq; ++j)
            {
                m_fields[i]->UpdatePhys()[j] = sin(x[j]) * sin(y[j]);
            }

            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                  m_fields[i]->UpdateCoeffs());
        }
    }
};

std::string DummyEquationSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "Dummy", DummyEquationSystem::create, "Dummy equation system.");

EquationSystemSharedPtr EquationSystem_Create(
    std::string eqnSysName, LibUtilities::SessionReaderSharedPtr session,
    SpatialDomains::MeshGraphSharedPtr mesh)
{
    EquationSystemSharedPtr ret =
        GetEquationSystemFactory().CreateInstance(eqnSysName, session, mesh);
    return ret;
}

Array<OneD, MultiRegions::ExpListSharedPtr> EquationSystem_GetFields(
    EquationSystemSharedPtr eqSys)
{
    return eqSys->UpdateFields();
}

void EquationSystem_PrintSummary(EquationSystemSharedPtr eqSys)
{
    eqSys->PrintSummary(std::cout);
}

std::shared_ptr<SessionFunction> EquationSystem_GetFunction1(
    EquationSystemSharedPtr eqSys, std::string name)
{
    return eqSys->GetFunction(name);
}

std::shared_ptr<SessionFunction> EquationSystem_GetFunction2(
    EquationSystemSharedPtr eqSys, std::string name,
    MultiRegions::ExpListSharedPtr field)
{
    return eqSys->GetFunction(name, field);
}

void EquationSystem_WriteFld(EquationSystemSharedPtr eqSys, std::string name)
{
    eqSys->WriteFld(name);
}

void EquationSystem_Checkpoint_Output(EquationSystemSharedPtr eqSys, int n)
{
    eqSys->Checkpoint_Output(n);
}

void export_EquationSystem(py::module &m)
{
    using EqSysWrap = EquationSystemWrap<EquationSystem>;
    using EqSysPub  = EquationSystemPublic<EquationSystem>;

    static NekFactory_Register<EquationSystemFactory> fac(
        GetEquationSystemFactory());

    py::classh<EquationSystem, EqSysWrap>(m, "EquationSystem")
        .def(py::init<LibUtilities::SessionReaderSharedPtr,
                      SpatialDomains::MeshGraphSharedPtr>())

        // Virtual functions that can be optionally overridden
        .def("InitObject", &EqSysPub::v_InitObject)
        .def("DoInitialise", &EqSysPub::v_DoInitialise)
        .def("DoSolve", &EqSysPub::v_DoSolve)
        .def("SetInitialConditions", &EqSysPub::v_SetInitialConditions)
        .def("EvaluateExactSolution", &EqSysWrap::EvaluateExactSolution)
        .def("LinfError", &EqSysWrap::LinfError)
        .def("L2Error", &EqSysWrap::L2Error)

        // Fields accessors (read-only)
        .def("GetFields", &EquationSystem_GetFields)
        .def_property_readonly("fields", &EquationSystem_GetFields)

        // Various utility functions
        .def("GetNvariables", &EquationSystem::GetNvariables)
        .def("GetVariable", &EquationSystem::GetVariable)
        .def("GetNpoints", &EquationSystem::GetNpoints)
        .def("SetInitialStep", &EquationSystem::SetInitialStep)
        .def("ZeroPhysFields", &EquationSystem::ZeroPhysFields)

        // Time accessors/properties
        .def("GetTime", &EquationSystem::GetTime)
        .def("SetTime", &EquationSystem::SetTime)
        .def_property("time", &EquationSystem::GetTime,
                      &EquationSystem::SetTime)

        // Timestep accessors/properties
        .def("GetTimeStep", &EquationSystem::GetTimeStep)
        .def("SetTimeStep", &EquationSystem::SetTimeStep)
        .def_property("timestep", &EquationSystem::GetTimeStep,
                      &EquationSystem::SetTimeStep)

        // Steps accessors/properties
        .def("GetSteps", &EquationSystem::GetSteps)
        .def("SetSteps", &EquationSystem::SetSteps)
        .def_property("steps", &EquationSystem::GetSteps,
                      &EquationSystem::SetSteps)

        // Print a summary
        .def("PrintSummary", &EquationSystem_PrintSummary)

        // Access functions from the session file
        .def("GetFunction", &EquationSystem_GetFunction1)
        .def("GetFunction", &EquationSystem_GetFunction2)

        // I/O utility functions
        .def("WriteFld", &EquationSystem_WriteFld)
        .def("Checkpoint_Output", &EquationSystem_Checkpoint_Output)

        // Factory functions.
        .def_static("Create", &EquationSystem_Create)
        .def_static("Register", [](std::string &filterName, py::object &obj) {
            fac(filterName, obj, filterName);
        });
}
