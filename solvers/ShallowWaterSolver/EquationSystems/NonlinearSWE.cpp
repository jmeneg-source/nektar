///////////////////////////////////////////////////////////////////////////////
//
// File: NonlinearSWE.cpp
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
// Description: Nonlinear Shallow water equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/NonlinearSWE.h>

namespace Nektar
{

std::string NonlinearSWE::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "NonlinearSWE", NonlinearSWE::create,
        "Nonlinear shallow water equation in conservative variables.");

NonlinearSWE::NonlinearSWE(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : ShallowWaterSystem(pSession, pGraph)
{
}

void NonlinearSWE::v_InitObject(bool DeclareFields)
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
            m_advection->SetFluxVector(&NonlinearSWE::GetFluxVector, this);

            // Setting up Riemann solver for advection operator
            m_session->LoadSolverInfo("UpwindType", riemName, "Average");
            m_riemannSolver =
                SolverUtils::GetRiemannSolverFactory().CreateInstance(
                    riemName, m_session);

            // Setting up parameters for advection operator Riemann solver
            m_riemannSolver->SetParam("gravity", &NonlinearSWE::GetGravity,
                                      this);
            m_riemannSolver->SetAuxVec("vecLocs", &NonlinearSWE::GetVecLocs,
                                       this);
            m_riemannSolver->SetVector("N", &NonlinearSWE::GetNormals, this);
            m_riemannSolver->SetScalar("depth", &NonlinearSWE::GetDepth, this);

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

    m_ode.DefineOdeRhs(&NonlinearSWE::v_DoOdeRhs, this);
    m_ode.DefineProjection(&NonlinearSWE::DoOdeProjection, this);
    m_ode.DefineImplicitSolve(&NonlinearSWE::DoImplicitSolve, this);
}

void NonlinearSWE::v_DoOdeRhs(
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

            // Variable Depth
            if (m_constantDepth != true)
            {
                AddVariableDepth(inarray, outarray);
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

            NonlinearSWE::GetFluxVector(inarray, fluxvector);
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

            // Variable Depth
            if (m_constantDepth != true)
            {
                AddVariableDepth(inarray, outarray);
            }
            //-------------------------------------------------
        }
        break;
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }
}

void NonlinearSWE::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    ShallowWaterSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Variables", "h  should be in field[0]");
    SolverUtils::AddSummaryItem(s, "", "hu should be in field[1]");
    SolverUtils::AddSummaryItem(s, "", "hv should be in field[2]");
}

// Physfield in conservative Form
void NonlinearSWE::GetFluxVector(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    int nq = m_fields[0]->GetTotPoints();

    NekDouble g = m_g;
    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);

    // Flux vector for the mass equation
    for (int i = 0; i < m_spacedim; ++i)
    {
        velocity[i] = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, physfield[i + 1], 1, flux[0][i], 1);
    }

    GetVelocityVector(physfield, velocity);

    // Put (0.5 g h h) in tmp
    Array<OneD, NekDouble> tmp(nq);
    Vmath::Vmul(nq, physfield[0], 1, physfield[0], 1, tmp, 1);
    Vmath::Smul(nq, 0.5 * g, tmp, 1, tmp, 1);

    // Flux vector for the momentum equations
    for (int i = 0; i < m_spacedim; ++i)
    {
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, physfield[i + 1], 1, flux[i + 1][j],
                        1);
        }

        // Add (0.5 g h h) to appropriate field
        Vmath::Vadd(nq, flux[i + 1][i], 1, tmp, 1, flux[i + 1][i], 1);
    }
}

/**
 * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
 * \f$ h\mathbf{v} \f$.
 *
 * @param physfield  Momentum field.
 * @param velocity   Velocity field.
 */
void NonlinearSWE::GetVelocityVector(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    const int npts = physfield[0].size();

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vdiv(npts, physfield[1 + i], 1, physfield[0], 1, velocity[i], 1);
    }
}

// physarray contains the conservative variables
void NonlinearSWE::AddVariableDepth(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int ncoeffs = GetNcoeffs();
    int nq      = GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mod(ncoeffs);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vmul(nq, m_bottomSlope[i], 1, physarray[0], 1, tmp, 1);
                Vmath::Smul(nq, m_g, tmp, 1, tmp, 1);
                m_fields[0]->IProductWRTBase(tmp, mod);
                m_fields[0]->MultiplyByElmtInvMass(mod, mod);
                m_fields[0]->BwdTrans(mod, tmp);
                Vmath::Vadd(nq, tmp, 1, outarray[i + 1], 1, outarray[i + 1], 1);
            }
        }
        break;
        case MultiRegions::eGalerkin:
        {
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vmul(nq, m_bottomSlope[i], 1, physarray[0], 1, tmp, 1);
                Vmath::Smul(nq, m_g, tmp, 1, tmp, 1);
                Vmath::Vadd(nq, tmp, 1, outarray[i + 1], 1, outarray[i + 1], 1);
            }
        }
        break;
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }
}

} // namespace Nektar
