///////////////////////////////////////////////////////////////////////////////
//
// File: VCSFSI.cpp
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
// Description: Velocity Correction Scheme for fluid-structure interaction
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VCSFSI.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

namespace Nektar
{
std::string VCSFSI::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "VCSFSI", VCSFSI::create);

std::string VCSFSI::solverTypeLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("SolverType", "VCSFSI",
                                                   eVCSFSI);

/**
 * Constructor. Creates ...
 *
 * \param
 * \param
 */
VCSFSI::VCSFSI(const LibUtilities::SessionReaderSharedPtr &pSession,
               const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph),
      VelocityCorrectionScheme(pSession, pGraph)
{
}

void VCSFSI::v_InitObject(bool DeclareField)
{
    VelocityCorrectionScheme::v_InitObject(DeclareField);
    Array<OneD, NekDouble> tmp = m_movingFrameData + 18;
    m_rigidSolver.InitObject(m_session, m_fields[0], tmp);
    m_rigidSolver.SetMovableDoFs(m_movableDoFs);
}

/**
 * Destructor
 */
VCSFSI::~VCSFSI(void)
{
}

void VCSFSI::v_DoInitialise(bool dumpInitialConditions)
{
    m_rigidSolver.SetInitialConditions(m_session, m_movingFrameData);
    VelocityCorrectionScheme::v_DoInitialise(dumpInitialConditions);
    Array<OneD, NekDouble> AddedMass;
    m_rigidSolver.SetNewmarkBetaSolver(AddedMass);
    Array<OneD, NekDouble> aeroforce(12, 0.);
    InitialiseFilter(aeroforce);
    m_rigidSolver.SetOldFvis(aeroforce);
}

void VCSFSI::InitialiseFilter(Array<OneD, NekDouble> aeroforce)
{
    if (!m_rigidSolver.HasFreeDoFs())
    {
        return;
    }
    std::map<std::string, std::string> vParams;
    m_rigidSolver.GetFilterInfo(m_session, vParams);
    vParams["OutputFile"] = ".dummyMRFForceFile";
    m_aeroforceFilter =
        MemoryManager<SolverUtils::FilterAeroForces>::AllocateSharedPtr(
            m_session, shared_from_this(), vParams);
    m_aeroforceFilter->Initialise(m_fields, 0.0);
    m_aeroforceFilter->GetForces(m_fields, NullNekDouble1DArray, 0.);
    GetAeroForce(aeroforce);
}

void VCSFSI::v_SolveSolid(NekDouble time)
{
    // call rigid solver
    Array<OneD, NekDouble> aeroforce(12, 0.);
    if (m_rigidSolver.HasFreeDoFs())
    {
        m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());
        m_aeroforceFilter->GetForces(m_fields, NullNekDouble1DArray, time);
        GetAeroForce(aeroforce);
    }
    // 0-5 pressure force at n+1; 6-11 viscous force at n
    m_rigidSolver.UpdateFrameVelocity(aeroforce, time, m_movingFrameData);
    // update velocity boundary condition
    UpdateVelocityBCs(time);
}

} // namespace Nektar
