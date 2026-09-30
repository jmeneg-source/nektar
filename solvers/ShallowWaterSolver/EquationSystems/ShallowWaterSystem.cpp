///////////////////////////////////////////////////////////////////////////////
//
// File: ShallowWaterSystem.cpp
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
// Description: Generic timestepping for shallow water solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

namespace Nektar
{
/**
 * @class ShallowWaterSystem
 *
 * Provides the underlying timestepping framework for shallow water flow solvers
 * including the general timestepping routines. This class is not intended
 * to be directly instantiated, but rather is a base class on which to
 * define shallow water solvers, e.g. SWE, Boussinesq, linear and nonlinear
 * versions.
 *
 * For details on implementing unsteady solvers see
 * \ref sectionADRSolverModuleImplementation
 */

/**
 * Processes SolverInfo parameters from the session file and sets up
 * timestepping-specific code.
 * @param   pSession        Session object to read parameters from.
 */

std::string ShallowWaterSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "ShallowWaterSystem", ShallowWaterSystem::create,
        "Auxiliary functions for the shallow water system.");

ShallowWaterSystem::ShallowWaterSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

void ShallowWaterSystem::v_InitObject(bool DeclareFields)
{
    UnsteadySystem::v_InitObject(DeclareFields);

    // Set up locations of velocity vector.
    m_vecLocs    = Array<OneD, Array<OneD, NekDouble>>(1);
    m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_vecLocs[0][i] = 1 + i;
    }

    // Load acceleration of gravity
    m_session->LoadParameter("Gravity", m_g, 9.81);

    EvaluateWaterDepth();

    m_constantDepth = true;
    NekDouble depth = m_depth[0];
    for (int i = 0; i < GetTotPoints(); ++i)
    {
        if (m_depth[i] != depth)
        {
            m_constantDepth = false;
            break;
        }
    }

    // Compute the bottom slopes
    int nq = GetTotPoints();
    if (m_constantDepth != true)
    {
        m_bottomSlope = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_bottomSlope[i] = Array<OneD, NekDouble>(nq);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], m_depth,
                                   m_bottomSlope[i]);
            Vmath::Neg(nq, m_bottomSlope[i], 1);
        }
    }

    EvaluateCoriolis();

    if (!m_explicitAdvection)
    {
        InitialiseNonlinSysSolver();
    }
}

void ShallowWaterSystem::v_DoOdeRhs(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time)
{
}

void ShallowWaterSystem::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    UnsteadySystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Depth",
                                m_constantDepth ? "constant" : "variable");
}

void ShallowWaterSystem::InitialiseNonlinSysSolver()
{
    unsigned int nvariables = m_fields.size();
    int ntotal              = nvariables * m_fields[0]->GetNpoints();

    // Create the key to hold settings for nonlin solver
    LibUtilities::NekSysKey key = LibUtilities::NekSysKey();

    // Load required LinSys parameters:
    m_session->LoadParameter("NekLinSysMaxIterations",
                             key.m_NekLinSysMaxIterations, 30);
    m_session->LoadParameter("LinSysMaxStorage", key.m_LinSysMaxStorage, 30);
    m_session->LoadParameter("LinSysRelativeTolInNonlin",
                             key.m_NekLinSysTolerance, 5.0E-2);
    m_session->LoadParameter("GMRESMaxHessMatBand", key.m_KrylovMaxHessMatBand,
                             31);

    // Load required NonLinSys parameters:
    m_session->LoadParameter("JacobiFreeEps", m_jacobiFreeEps, 5.0E-8);
    m_session->LoadParameter("NekNonlinSysMaxIterations",
                             key.m_NekNonlinSysMaxIterations, 10);
    m_session->LoadParameter("NewtonRelativeIteTol",
                             key.m_NekNonLinSysTolerance, 1.0E-12);
    WARNINGL0(!m_session->DefinesParameter("NewtonAbsoluteIteTol"),
              "Please specify NewtonRelativeIteTol instead of "
              "NewtonAbsoluteIteTol in XML session file");
    m_session->LoadParameter("NonlinIterTolRelativeL2",
                             key.m_NonlinIterTolRelativeL2, 1.0E-3);
    m_session->LoadSolverInfo("LinSysIterSolverTypeInNonlin",
                              key.m_LinSysIterSolverTypeInNonlin, "GMRES");

    LibUtilities::NekSysOperators nekSysOp;
    nekSysOp.DefineNekSysResEval(&ShallowWaterSystem::NonlinSysEvaluator1D,
                                 this);
    nekSysOp.DefineNekSysLhsEval(&ShallowWaterSystem::MatrixMultiplyMatrixFree,
                                 this);
    nekSysOp.DefineNekSysPrecon(&ShallowWaterSystem::DoNullPrecon, this);

    // Initialize non-linear system
    m_nonlinsol = LibUtilities::GetNekNonlinSysIterFactory().CreateInstance(
        "Newton", m_session, m_comm->GetRowComm(), ntotal, key);
    m_nonlinsol->SetSysOperators(nekSysOp);
}

void ShallowWaterSystem::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
    Array<OneD, Array<OneD, NekDouble>> &outpnt, const NekDouble time,
    const NekDouble lambda)
{
    m_TimeIntegLambda       = lambda;
    m_bndEvaluateTime       = time;
    unsigned int npoints    = m_fields[0]->GetNpoints();
    unsigned int nvariables = m_fields.size();

    Array<OneD, NekDouble> inarray(nvariables * npoints);
    Array<OneD, NekDouble> outarray(nvariables * npoints);
    Array<OneD, NekDouble> tmp;

    for (int i = 0; i < nvariables; ++i)
    {
        int noffset = i * npoints;
        Vmath::Vcopy(npoints, inpnts[i], 1, tmp = inarray + noffset, 1);
    }

    DoImplicitSolve1D(inarray, outarray);

    for (int i = 0; i < nvariables; ++i)
    {
        int noffset = i * npoints;
        Vmath::Vcopy(npoints, outarray + noffset, 1, outpnt[i], 1);
    }
}

void ShallowWaterSystem::DoImplicitSolve1D(
    const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &out)
{
    CalcRefValues(inarray);

    m_nonlinsol->SetRhsMagnitude(m_inArrayNorm);

    m_TotNewtonIts += m_nonlinsol->SolveSystem(inarray.size(), inarray, out, 0);

    m_TotLinIts += m_nonlinsol->GetNtotLinSysIts();

    m_TotImpStages++;
}

void ShallowWaterSystem::CalcRefValues(
    const Array<OneD, const NekDouble> &inarray)
{
    unsigned int npoints = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> magnitdEstimat(3, 0.0);

    for (int i = 0; i < 3; ++i)
    {
        int offset = i * npoints;
        magnitdEstimat[i] =
            Vmath::Dot(npoints, inarray + offset, inarray + offset);
    }
    m_comm->GetSpaceComm()->AllReduce(magnitdEstimat,
                                      Nektar::LibUtilities::ReduceSum);

    m_inArrayNorm = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        m_inArrayNorm += magnitdEstimat[i];
    }
}

void ShallowWaterSystem::NonlinSysEvaluator1D(
    const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &out,
    [[maybe_unused]] const bool &flag)
{
    unsigned int npoints    = m_fields[0]->GetNpoints();
    unsigned int nvariables = m_fields.size();
    Array<OneD, Array<OneD, NekDouble>> in2D(nvariables);
    Array<OneD, Array<OneD, NekDouble>> out2D(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        int offset = i * npoints;
        in2D[i]    = inarray + offset;
        out2D[i]   = out + offset;
    }
    NonlinSysEvaluator(in2D, out2D);
}

void ShallowWaterSystem::NonlinSysEvaluator(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &out)
{
    unsigned int npoints    = m_fields[0]->GetNpoints();
    unsigned int nvariables = m_fields.size();
    Array<OneD, Array<OneD, NekDouble>> inpnts(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        inpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    DoOdeProjection(inarray, inpnts, m_bndEvaluateTime);
    v_DoOdeRhs(inpnts, out, m_bndEvaluateTime);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Svtvp(npoints, -m_TimeIntegLambda, out[i], 1, inarray[i], 1,
                     out[i], 1);
        Vmath::Vsub(npoints, out[i], 1,
                    m_nonlinsol->GetRefSourceVec() + i * npoints, 1, out[i], 1);
    }
}

void ShallowWaterSystem::MatrixMultiplyMatrixFree(
    const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &out,
    [[maybe_unused]] const bool &flag)
{
    const Array<OneD, const NekDouble> solref = m_nonlinsol->GetRefSolution();
    const Array<OneD, const NekDouble> resref = m_nonlinsol->GetRefResidual();

    unsigned int ntotal   = inarray.size();
    NekDouble magninarray = Vmath::Dot(ntotal, inarray, inarray);
    m_comm->GetSpaceComm()->AllReduce(magninarray,
                                      Nektar::LibUtilities::ReduceSum);
    NekDouble eps =
        m_jacobiFreeEps * sqrt((sqrt(m_inArrayNorm) + 1.0) / magninarray);

    Array<OneD, NekDouble> solplus{ntotal};
    Array<OneD, NekDouble> resplus{ntotal};

    Vmath::Svtvp(ntotal, eps, inarray, 1, solref, 1, solplus, 1);
    NonlinSysEvaluator1D(solplus, resplus, flag);
    Vmath::Vsub(ntotal, resplus, 1, resref, 1, out, 1);
    Vmath::Smul(ntotal, 1.0 / eps, out, 1, out, 1);
}

void ShallowWaterSystem::DoNullPrecon(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, [[maybe_unused]] const bool &flag)
{
    Vmath::Vcopy(inarray.size(), inarray, 1, outarray, 1);
}

void ShallowWaterSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            if (inarray != outarray)
            {
                int npoints = GetNpoints();

                for (int i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
            }

            SetBoundaryConditions(outarray, time);
            break;
        }
        case MultiRegions::eGalerkin:
        {
            EquationSystem::SetBoundaryConditions(time);
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(), 0.0);

            for (int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans(coeffs, outarray[i]);
            }
            break;
        }
        default:
            ASSERTL0(false, "Unknown projection scheme");
            break;
    }
}

void ShallowWaterSystem::SetBoundaryConditions(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray, NekDouble time)
{
    std::string varName;
    int nvariables = 3;
    int cnt        = 0;
    int nTracePts  = GetTraceTotPoints();

    // Extract trace for boundaries. Needs to be done on all processors to avoid
    // deadlock.
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }

    // Loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),
                           "Wall"))
        {
            WallBoundary2D(n, cnt, Fwd);
        }

        // Time Dependent Boundary Condition (specified in meshfile)
        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
        {
            for (int i = 0; i < nvariables; ++i)
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

void ShallowWaterSystem::WallBoundary2D(
    int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &Fwd)
{
    int nvariables = 3;

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int id1, id2, npts;

    for (int e = 0;
         e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
    {
        npts = m_fields[0]
                   ->GetBndCondExpansions()[bcRegion]
                   ->GetExp(e)
                   ->GetNumPoints(0);
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        switch (m_expdim)
        {
            case 1:
            {
                // negate the forward flux
                Vmath::Neg(npts, &Fwd[1][id2], 1);
            }
            break;
            case 2:
            {
                Array<OneD, NekDouble> tmp_n(npts);
                Array<OneD, NekDouble> tmp_t(npts);

                Vmath::Vmul(npts, &Fwd[1][id2], 1, &m_traceNormals[0][id2], 1,
                            &tmp_n[0], 1);
                Vmath::Vvtvp(npts, &Fwd[2][id2], 1, &m_traceNormals[1][id2], 1,
                             &tmp_n[0], 1, &tmp_n[0], 1);

                Vmath::Vmul(npts, &Fwd[1][id2], 1, &m_traceNormals[1][id2], 1,
                            &tmp_t[0], 1);
                Vmath::Vvtvm(npts, &Fwd[2][id2], 1, &m_traceNormals[0][id2], 1,
                             &tmp_t[0], 1, &tmp_t[0], 1);

                // negate the normal flux
                Vmath::Neg(npts, tmp_n, 1);

                // rotate back to Cartesian
                Vmath::Vmul(npts, &tmp_t[0], 1, &m_traceNormals[1][id2], 1,
                            &Fwd[1][id2], 1);
                Vmath::Vvtvm(npts, &tmp_n[0], 1, &m_traceNormals[0][id2], 1,
                             &Fwd[1][id2], 1, &Fwd[1][id2], 1);

                Vmath::Vmul(npts, &tmp_t[0], 1, &m_traceNormals[0][id2], 1,
                            &Fwd[2][id2], 1);
                Vmath::Vvtvp(npts, &tmp_n[0], 1, &m_traceNormals[1][id2], 1,
                             &Fwd[2][id2], 1, &Fwd[2][id2], 1);
            }
            break;
            case 3:
                ASSERTL0(false,
                         "3D not implemented for Shallow Water Equations");
                break;
            default:
                ASSERTL0(false, "Illegal expansion dimension");
        }

        // copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npts, &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

// physarray contains the conservative variables
void ShallowWaterSystem::AddCoriolis(
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
            // add to u equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[2], 1, tmp, 1);
            m_fields[0]->IProductWRTBase(tmp, mod);
            m_fields[0]->MultiplyByElmtInvMass(mod, mod);
            m_fields[0]->BwdTrans(mod, tmp);
            Vmath::Vadd(nq, tmp, 1, outarray[1], 1, outarray[1], 1);

            // add to v equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
            Vmath::Neg(nq, tmp, 1);
            m_fields[0]->IProductWRTBase(tmp, mod);
            m_fields[0]->MultiplyByElmtInvMass(mod, mod);
            m_fields[0]->BwdTrans(mod, tmp);
            Vmath::Vadd(nq, tmp, 1, outarray[2], 1, outarray[2], 1);
        }
        break;
        case MultiRegions::eGalerkin:
        {
            // add to u equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[2], 1, tmp, 1);
            Vmath::Vadd(nq, tmp, 1, outarray[1], 1, outarray[1], 1);

            // add to v equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
            Vmath::Neg(nq, tmp, 1);
            Vmath::Vadd(nq, tmp, 1, outarray[2], 1, outarray[2], 1);
        }
        break;
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }
}

void ShallowWaterSystem::ConservativeToPrimitive()
{
    int nq = GetTotPoints();

    // \eta = h - d
    Vmath::Vsub(nq, m_fields[0]->GetPhys(), 1, m_depth, 1,
                m_fields[0]->UpdatePhys(), 1);

    // u = hu / h
    Vmath::Vdiv(nq, m_fields[1]->GetPhys(), 1, m_fields[0]->GetPhys(), 1,
                m_fields[1]->UpdatePhys(), 1);

    // v = hv / v
    Vmath::Vdiv(nq, m_fields[2]->GetPhys(), 1, m_fields[0]->GetPhys(), 1,
                m_fields[2]->UpdatePhys(), 1);
}

void ShallowWaterSystem::PrimitiveToConservative()
{
    int nq = GetTotPoints();

    // h = \eta + d
    Vmath::Vadd(nq, m_fields[0]->GetPhys(), 1, m_depth, 1,
                m_fields[0]->UpdatePhys(), 1);

    // hu = h * u
    Vmath::Vmul(nq, m_fields[0]->GetPhys(), 1, m_fields[1]->GetPhys(), 1,
                m_fields[1]->UpdatePhys(), 1);

    // hv = h * v
    Vmath::Vmul(nq, m_fields[0]->GetPhys(), 1, m_fields[2]->GetPhys(), 1,
                m_fields[2]->UpdatePhys(), 1);
}

void ShallowWaterSystem::EvaluateWaterDepth(void)
{
    GetFunction("WaterDepth")->Evaluate("d", m_depth);
}

void ShallowWaterSystem::EvaluateCoriolis(void)
{
    GetFunction("Coriolis")->Evaluate("f", m_coriolis);
}

} // namespace Nektar
