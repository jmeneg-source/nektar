///////////////////////////////////////////////////////////////////////////////
//
// File: LinearSWE.cpp
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
// Description: Linear Shallow water equations in primitive variables
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/LinearSWE.h>

namespace Nektar
{

std::string LinearSWE::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "LinearSWE", LinearSWE::create,
        "Linear shallow water equation in primitive variables.");

LinearSWE::LinearSWE(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : ShallowWaterSystem(pSession, pGraph)
{
}

void LinearSWE::v_InitObject(bool DeclareFields)
{
    ShallowWaterSystem::v_InitObject(DeclareFields);

    // Type of advection class to be used
    switch (m_projectionType)
    {
        // Continuous field
        case MultiRegions::eGalerkin:
        {
            //  Do nothing
            break;
        }
        // Discontinuous field
        case MultiRegions::eDiscontinuous:
        {
            std::string advName;
            std::string diffName;
            std::string riemName;

            //---------------------------------------------------------------
            // Setting up advection and diffusion operators
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
            m_advection = SolverUtils::GetAdvectionFactory().CreateInstance(
                advName, advName);
            m_advection->SetFluxVector(&LinearSWE::GetFluxVector, this);

            // Setting up Riemann solver for advection operator
            m_session->LoadSolverInfo("UpwindType", riemName, "NoSolver");
            if ((riemName == "LinearHLL") && (m_constantDepth != true))
            {
                ASSERTL0(false, "LinearHLL only valid for constant depth");
            }
            m_riemannSolver =
                SolverUtils::GetRiemannSolverFactory().CreateInstance(
                    riemName, m_session);

            // Setting up parameters for advection operator Riemann solver
            m_riemannSolver->SetParam("gravity", &LinearSWE::GetGravity, this);
            m_riemannSolver->SetAuxVec("vecLocs", &LinearSWE::GetVecLocs, this);
            m_riemannSolver->SetVector("N", &LinearSWE::GetNormals, this);

            // The numerical flux for linear SWE requires depth information
            int nTracePointsTot = m_fields[0]->GetTrace()->GetTotPoints();
            m_dFwd              = Array<OneD, NekDouble>(nTracePointsTot);
            m_dBwd              = Array<OneD, NekDouble>(nTracePointsTot);
            m_fields[0]->GetFwdBwdTracePhys(m_depth, m_dFwd, m_dBwd);
            CopyBoundaryTrace(m_dFwd, m_dBwd);
            m_riemannSolver->SetScalar("depthFwd", &LinearSWE::GetDepthFwd,
                                       this);
            m_riemannSolver->SetScalar("depthBwd", &LinearSWE::GetDepthBwd,
                                       this);

            // Concluding initialisation of advection operators
            m_advection->SetRiemannSolver(m_riemannSolver);
            m_advection->InitObject(m_session, m_fields);
            break;
        }
        default:
        {
            ASSERTL0(false, "Unsupported projection type.");
            break;
        }
    }

    m_ode.DefineOdeRhs(&LinearSWE::v_DoOdeRhs, this);
    m_ode.DefineProjection(&LinearSWE::DoOdeProjection, this);
    m_ode.DefineImplicitSolve(&LinearSWE::DoImplicitSolve, this);
}

void LinearSWE::v_DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int ndim       = m_spacedim;
    int nvariables = inarray.size();
    int nq         = GetTotPoints();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            //-------------------------------------------------------
            // Compute the DG advection including the numerical flux
            // by using SolverUtils/Advection
            // Input and output in physical space
            m_advection->Advect(nvariables, m_fields, NullNekDoubleArrayOfArray,
                                inarray, outarray, time);
            //-------------------------------------------------------

            //-------------------------------------------------------
            // negate the outarray since moving terms to the rhs
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(nq, outarray[i], 1);
            }
            //-------------------------------------------------------

            //-------------------------------------------------
            // Add "source terms"
            // Input and output in physical space

            // Coriolis forcing
            if (m_coriolis.size() != 0)
            {
                AddCoriolis(inarray, outarray);
            }
            //-------------------------------------------------
        }
        break;
        case MultiRegions::eGalerkin:
        {
            //-------------------------------------------------------
            // Compute the fluxvector in physical space
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>> fluxvector(
                nvariables);

            for (int i = 0; i < nvariables; ++i)
            {
                fluxvector[i] = Array<OneD, Array<OneD, NekDouble>>(ndim);
                for (int j = 0; j < ndim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nq);
                }
            }

            LinearSWE::GetFluxVector(inarray, fluxvector);
            //-------------------------------------------------------

            //-------------------------------------------------------
            // Take the derivative of the flux terms
            // and negate the outarray since moving terms to the rhs
            Array<OneD, NekDouble> tmp0(nq);
            Array<OneD, NekDouble> tmp1(nq);

            for (int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[0],
                                       fluxvector[i][0], tmp0);
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[1],
                                       fluxvector[i][1], tmp1);
                Vmath::Vadd(nq, tmp0, 1, tmp1, 1, outarray[i], 1);
                Vmath::Neg(nq, outarray[i], 1);
            }

            //-------------------------------------------------
            // Add "source terms"
            // Input and output in physical space

            // Coriolis forcing
            if (m_coriolis.size() != 0)
            {
                AddCoriolis(inarray, outarray);
            }
            //-------------------------------------------------
        }
        break;
        default:
            ASSERTL0(false, "Unknown projection scheme for the LinearSWE");
            break;
    }
}

void LinearSWE::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    ShallowWaterSystem::v_GenerateSummary(s);
    if (m_session->DefinesSolverInfo("UpwindType"))
    {
        std::string UpwindType;
        UpwindType = m_session->GetSolverInfo("UpwindType");
        if (UpwindType == "LinearAverage")
        {
            SolverUtils::AddSummaryItem(s, "Riemann Solver", "Linear Average");
        }
        else if (UpwindType == "LinearHLL")
        {
            SolverUtils::AddSummaryItem(s, "Riemann Solver", "Linear HLL");
        }
    }
    SolverUtils::AddSummaryItem(s, "Variables", "eta  should be in field[0]");
    SolverUtils::AddSummaryItem(s, "", "u    should be in field[1]");
    SolverUtils::AddSummaryItem(s, "", "v    should be in field[2]");
}

// Physfield in primitive Form
void LinearSWE::GetFluxVector(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    int nq = m_fields[0]->GetTotPoints();

    NekDouble g = m_g;

    // Flux vector for the mass equation
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vmul(nq, m_depth, 1, physfield[i + 1], 1, flux[0][i], 1);
    }

    // Flux vector for the momentum equations
    for (int i = 0; i < m_spacedim; ++i)
    {
        for (int j = 0; j < m_spacedim; ++j)
        {
            // must zero fluxes as not initialised to zero in AdvectionWeakDG
            Vmath::Zero(nq, flux[i + 1][j], 1);
        }

        // Add (g eta) to appropriate field
        Vmath::Smul(nq, g, physfield[0], 1, flux[i + 1][i], 1);
    }
}

/**
 * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
 * \f$ h\mathbf{v} \f$.
 *
 * @param physfield  Velocity field.
 * @param velocity   Velocity field.
 */
void LinearSWE::GetVelocityVector(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    const int npts = physfield[0].size();

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vcopy(npts, physfield[1 + i], 1, velocity[i], 1);
    }
}

void LinearSWE::CopyBoundaryTrace(const Array<OneD, const NekDouble> &Fwd,
                                  Array<OneD, NekDouble> &Bwd)
{
    int cnt = 0;

    // loop over Boundary Regions
    for (int bcRegion = 0; bcRegion < m_fields[0]->GetBndConditions().size();
         ++bcRegion)
    {
        if (m_fields[0]
                ->GetBndConditions()[bcRegion]
                ->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Copy the forward trace of the field to the backward trace
        int id2, npts;

        for (int e = 0;
             e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
             ++e)
        {
            npts = m_fields[0]
                       ->GetBndCondExpansions()[bcRegion]
                       ->GetExp(e)
                       ->GetTotPoints();
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt +
                                                                        e));

            Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
        }

        cnt += m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
    }
}

} // namespace Nektar
