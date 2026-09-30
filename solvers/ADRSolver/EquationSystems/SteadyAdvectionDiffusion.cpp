/////////////////////////////////////////////////////////////////////////////////
//
// File: SteadyAdvectionDiffusion.cpp
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
// Description: Steady advection-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/SteadyAdvectionDiffusion.h>

namespace Nektar
{

std::string SteadyAdvectionDiffusion::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "SteadyAdvectionDiffusion", SteadyAdvectionDiffusion::create);

/**
 * @class SteadyAdvectionDiffusion
 * This is a solver class for solving the  problems.
 * - SteadyAdvectionDiffusion:
 *   \f$ c \cdot \nabla u -\nabla \cdot (\nabla u)  = f(x)\f$
 */

SteadyAdvectionDiffusion::SteadyAdvectionDiffusion(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph), m_lambda(0.0)
{
}

void SteadyAdvectionDiffusion::v_InitObject(bool DeclareFields)
{
    EquationSystem::v_InitObject(DeclareFields);

    m_session->LoadParameter("epsilon", m_epsilon, 1.0);

    std::vector<std::string> vel;
    vel.push_back("Vx");
    vel.push_back("Vy");
    vel.push_back("Vz");

    // Resize the advection velocities vector to dimension of the problem
    vel.resize(m_spacedim);

    // Store in the global variable m_velocity the advection velocities
    m_velocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    GetFunction("BaseFlow")->Evaluate(vel, m_velocity);
}

void SteadyAdvectionDiffusion::v_GenerateSummary(
    [[maybe_unused]] SolverUtils::SummaryList &s)
{
}

void SteadyAdvectionDiffusion::v_DoInitialise(
    [[maybe_unused]] bool dumpInitialConditions)
{
    // Set initial forcing from session file
    GetFunction("Forcing")->Evaluate(m_session->GetVariables(), m_fields);
}

void SteadyAdvectionDiffusion::v_DoSolve()
{
    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;

    factors[StdRegions::eFactorLambda] = m_lambda / m_epsilon;

    // Set advection velocities
    StdRegions::VarCoeffType varcoefftypes[] = {StdRegions::eVarCoeffVelX,
                                                StdRegions::eVarCoeffVelY,
                                                StdRegions::eVarCoeffVelZ};
    for (int i = 0; i < m_spacedim; i++)
    {
        // Scale advection velocities by diffusion coefficient
        Vmath::Smul(m_velocity[i].size(), 1.0 / m_epsilon, m_velocity[i], 1,
                    m_velocity[i], 1);
        varcoeffs[varcoefftypes[i]] = m_velocity[i];
    }

    // Solve for velocity
    for (int i = 0; i < m_fields.size(); ++i)
    {
        // Zero initial guess
        Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(), 1);
        // Scale forcing term by diffusion coefficient
        Vmath::Smul(m_fields[i]->GetTotPoints(), 1.0 / m_epsilon,
                    m_fields[i]->GetPhys(), 1, m_fields[i]->UpdatePhys(), 1);
        // Solve system
        m_fields[i]->LinearAdvectionDiffusionReactionSolve(
            m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(), factors,
            varcoeffs);
        m_fields[i]->SetPhysState(false);
    }
}

} // namespace Nektar
