///////////////////////////////////////////////////////////////////////////////
//
// File: UnsteadyInviscidBurgers.cpp
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
// Description: Unsteady inviscid Burgers solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/UnsteadyInviscidBurgers.h>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{

std::string UnsteadyInviscidBurgers::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "UnsteadyInviscidBurgers", UnsteadyInviscidBurgers::create,
        "Inviscid Burgers equation");

UnsteadyInviscidBurgers::UnsteadyInviscidBurgers(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph)
{
}

/**
 * @brief Initialisation object for the inviscid Burgers equation.
 */
void UnsteadyInviscidBurgers::v_InitObject(bool DeclareFields)
{
    AdvectionSystem::v_InitObject(DeclareFields);

    // Type of advection class to be used
    switch (m_projectionType)
    {
        // Discontinuous field
        case MultiRegions::eDiscontinuous:
        {
            // Do not forwards transform initial condition
            m_homoInitialFwd = false;

            // Define the normal velocity fields
            if (m_fields[0]->GetTrace())
            {
                m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
            }

            // Advection term
            std::string advName;
            std::string riemName;
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
            m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(
                advName, advName);
            m_advObject->SetFluxVector(&UnsteadyInviscidBurgers::GetFluxVector,
                                       this);
            m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
            m_riemannSolver =
                SolverUtils::GetRiemannSolverFactory().CreateInstance(
                    riemName, m_session);
            m_riemannSolver->SetScalar(
                "Vn", &UnsteadyInviscidBurgers::GetNormalVelocity, this);
            m_advObject->SetRiemannSolver(m_riemannSolver);
            m_advObject->InitObject(m_session, m_fields);
            break;
        }
        // Continuous field
        case MultiRegions::eGalerkin:
        {
            std::string advName;
            m_session->LoadSolverInfo("AdvectionType", advName,
                                      "NonConservative");
            m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(
                advName, advName);
            m_advObject->SetFluxVector(&UnsteadyInviscidBurgers::GetFluxVector,
                                       this);
            break;
        }
        default:
        {
            ASSERTL0(false, "Unsupported projection type.");
            break;
        }
    }

    // Forcing terms
    m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                           m_fields, m_fields.size());

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&UnsteadyInviscidBurgers::DoOdeRhs, this);
        m_ode.DefineProjection(&UnsteadyInviscidBurgers::DoOdeProjection, this);
    }
    else
    {
        ASSERTL0(false, "Implicit unsteady Advection not set up.");
    }
}

/**
 * @brief Compute the right-hand side for the inviscid Burgers equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyInviscidBurgers::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Number of fields (variables of the problem)
    int nVariables = inarray.size();

    // Number of solution points
    int nSolutionPts = GetNpoints();

    LibUtilities::Timer timer;
    timer.Start();
    // RHS computation using the new advection base class
    m_advObject->Advect(nVariables, m_fields, inarray, inarray, outarray, time);
    timer.Stop();
    // Elapsed time
    timer.AccumulateRegion("Advect");

    // Negate the RHS
    for (int i = 0; i < nVariables; ++i)
    {
        Vmath::Neg(nSolutionPts, outarray[i], 1);
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        // set up non-linear terms
        x->Apply(m_fields, inarray, outarray, time);
    }
}

/**
 * @brief Compute the projection for the inviscid Burgers equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyInviscidBurgers::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Number of variables of the problem
    int nVariables = inarray.size();

    // Set the boundary conditions
    SetBoundaryConditions(time);

    // Switch on the projection type (Discontinuous or Continuous)
    switch (m_projectionType)
    {
        // Discontinuous projection
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            if (inarray != outarray)
            {
                int npoints = GetNpoints();

                for (int i = 0; i < nVariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
            }
            break;
        }
        // Continuous projection
        case MultiRegions::eGalerkin:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for (int i = 0; i < nVariables; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans(coeffs, outarray[i]);
            }
            break;
        }
        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
            break;
        }
    }
}

/**
 * @brief Get the normal velocity for the inviscid Burgers equation.
 */
Array<OneD, NekDouble> &UnsteadyInviscidBurgers::GetNormalVelocity()
{
    // Number of trace (interface) points
    int nTracePts = GetTraceNpoints();

    // Number of solution points
    int nSolutionPts = GetNpoints();

    // Auxiliary variables to compute the normal velocity
    Array<OneD, NekDouble> Fwd(nTracePts);
    Array<OneD, NekDouble> Bwd(nTracePts);
    Array<OneD, NekDouble> physfield(nSolutionPts);

    // Reset the normal velocity
    Vmath::Zero(nTracePts, m_traceVn, 1);

    for (int i = 0; i < m_spacedim; ++i)
    {
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), physfield);
        m_fields[i]->GetFwdBwdTracePhys(physfield, Fwd, Bwd, true);
        Vmath::Vadd(nTracePts, Fwd, 1, Bwd, 1, Fwd, 1);
        Vmath::Smul(nTracePts, 0.5, Fwd, 1, Fwd, 1);
        Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Fwd, 1, m_traceVn, 1,
                     m_traceVn, 1);
    }
    Vmath::Smul(nTracePts, 0.5, m_traceVn, 1, m_traceVn, 1);

    return m_traceVn;
}

/**
 * @brief Return the flux vector for the inviscid Burgers equation.
 *
 * @param physfield   Fields.
 * @param flux        Resulting flux.
 */
void UnsteadyInviscidBurgers::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    const int nq = GetNpoints();

    for (int i = 0; i < flux.size(); ++i)
    {
        for (int j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(nq, physfield[i], 1, physfield[i], 1, flux[i][j], 1);
            Vmath::Smul(nq, 0.5, flux[i][j], 1, flux[i][j], 1);
        }
    }
}

void UnsteadyInviscidBurgers::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    AdvectionSystem::v_GenerateSummary(s);
}

} // namespace Nektar
