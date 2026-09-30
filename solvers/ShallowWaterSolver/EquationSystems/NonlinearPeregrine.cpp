///////////////////////////////////////////////////////////////////////////////
//
// File: NonlinearPeregrine.cpp
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
// Description: Nonlinear Boussinesq equations of Peregrine in
//              conservative variables (constant depth case)
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/NonlinearPeregrine.h>

namespace Nektar
{

std::string NonlinearPeregrine::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "NonlinearPeregrine", NonlinearPeregrine::create,
        "Nonlinear Peregrine equations in conservative variables.");

NonlinearPeregrine::NonlinearPeregrine(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : NonlinearSWE(pSession, pGraph), m_factors()
{
    m_factors[StdRegions::eFactorLambda] = 0.0;
    m_factors[StdRegions::eFactorTau]    = 1000000.0;
    // note: eFactorTau = 1.0 becomes unstable...
    // we need to investigate the behaviuor w.r.t. tau
}

void NonlinearPeregrine::v_InitObject(bool DeclareFields)
{
    NonlinearSWE::v_InitObject(DeclareFields);

    if (m_session->DefinesSolverInfo("PROBLEMTYPE"))
    {
        std::string ProblemTypeStr = m_session->GetSolverInfo("PROBLEMTYPE");
        for (int i = 0; i < (int)SIZE_ProblemType; ++i)
        {
            if (boost::iequals(ProblemTypeMap[i], ProblemTypeStr))
            {
                m_problemType = (ProblemType)i;
                break;
            }
        }
    }
    else
    {
        m_problemType = (ProblemType)0;
    }

    // NB! At the moment only the constant depth case is
    // supported for the Peregrine eq.
    if (m_session->DefinesParameter("ConstDepth"))
    {
        m_const_depth = m_session->GetParameter("ConstDepth");
    }
    else
    {
        ASSERTL0(false, "Constant Depth not specified");
    }

    ASSERTL0(m_projectionType != MultiRegions::eGalerkin,
             "Continuous projection type not supported for Peregrine.");

    m_ode.DefineOdeRhs(&NonlinearPeregrine::v_DoOdeRhs, this);
    m_ode.DefineProjection(&NonlinearPeregrine::DoOdeProjection, this);
    m_ode.DefineImplicitSolve(&NonlinearPeregrine::DoImplicitSolve, this);
}

/**
 * @brief Set the initial conditions.
 */
void NonlinearPeregrine::v_SetInitialConditions(
    NekDouble initialtime, bool dumpInitialConditions,
    [[maybe_unused]] const int domain)
{
    switch (m_problemType)
    {
        case eSolitaryWave:
        {
            LaitoneSolitaryWave(0.1, m_const_depth, 0.0, 0.0);
            break;
        }
        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
            m_nchk--; // Note: m_nchk has been incremented in EquationSystem.
            break;
        }
    }

    // Update time in field info if required
    if (m_fieldMetaDataMap.find("Time") != m_fieldMetaDataMap.end())
    {
        m_fieldMetaDataMap["Time"] = boost::lexical_cast<std::string>(m_time);
    }

    if (dumpInitialConditions && m_checksteps && m_nchk == 0 &&
        !m_comm->IsParallelInTime())
    {
        Checkpoint_Output(m_nchk);
    }
    else if (dumpInitialConditions && m_nchk == 0 && m_comm->IsParallelInTime())
    {
        std::string newdir = m_sessionName + ".pit";
        if (!fs::is_directory(newdir))
        {
            fs::create_directory(newdir);
        }
        if (m_comm->GetTimeComm()->GetRank() == 0)
        {
            WriteFld(newdir + "/" + m_sessionName + "_0.fld");
        }
    }
    m_nchk++;
}

void NonlinearPeregrine::v_DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size() - 1;
    int ncoeffs    = GetNcoeffs();
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
                ASSERTL0(false,
                         "Variable depth not supported for the Peregrine "
                         "equations");
            }

            //-------------------------------------------------

            //---------------------------------------
            // As no more terms is required for the
            // continuity equation and we have aleady evaluated
            //  the values for h_t we are done for h
            //---------------------------------------

            //---------------------------------------------

            //-------------------------------------------------
            // create tmp fields to be used during
            // the dispersive section

            int nTraceNumPoints = GetTraceTotPoints();
            Array<OneD, NekDouble> upwindX(nTraceNumPoints);
            Array<OneD, NekDouble> upwindY(nTraceNumPoints);

            Array<OneD, Array<OneD, NekDouble>> coeffsfield(2);
            Array<OneD, Array<OneD, NekDouble>> physfield(2);
            for (int i = 0; i < 2; ++i)
            {
                coeffsfield[i] = Array<OneD, NekDouble>(ncoeffs);
                physfield[i]   = Array<OneD, NekDouble>(nq);
            }

            Vmath::Vcopy(nq, outarray[1], 1, physfield[0], 1);
            Vmath::Vcopy(nq, outarray[2], 1, physfield[1], 1);

            //---------------------------------------
            // Start for solve of mixed dispersive terms using the 'WCE method'
            // (Eskilsson & Sherwin, JCP 2006)

            // constant depth case
            // \nabla \cdot (\nabla z) - invgamma z
            //    = - invgamma (\nabla \cdot {\bf f}_(2,3)

            NekDouble gamma    = (m_const_depth * m_const_depth) / 3.0;
            NekDouble invgamma = 1.0 / gamma;
            //--------------------------------------------

            //--------------------------------------------
            // Compute the forcing function for the wave continuity equation (eq
            // 26a)

            // Set boundary condidtions for z
            SetBoundaryConditionsForcing(physfield, time);

            // \nabla \phi \cdot f_{2,3}
            m_fields[0]->IProductWRTDerivBase(0, physfield[0], coeffsfield[0]);
            m_fields[0]->IProductWRTDerivBase(1, physfield[1], coeffsfield[1]);
            Vmath::Vadd(ncoeffs, coeffsfield[0], 1, coeffsfield[1], 1,
                        coeffsfield[0], 1);
            Vmath::Neg(ncoeffs, coeffsfield[0], 1);

            // Evaluate  upwind numerical flux (physical space)
            NumericalFluxForcing(physfield, upwindX, upwindY);
            Array<OneD, NekDouble> normflux(nTraceNumPoints);
            Vmath::Vvtvvtp(nTraceNumPoints, upwindX, 1, m_traceNormals[0], 1,
                           upwindY, 1, m_traceNormals[1], 1, normflux, 1);
            m_fields[0]->AddTraceIntegral(normflux, coeffsfield[0]);
            m_fields[0]->MultiplyByElmtInvMass(coeffsfield[0], coeffsfield[0]);
            m_fields[0]->BwdTrans(coeffsfield[0], physfield[0]);

            Vmath::Smul(nq, -invgamma, physfield[0], 1, physfield[0], 1);

            //--------------------------------------

            //--------------------------------------
            // Solve the Helmhotz-type equation for the wave continuity equation
            // (eq. 26b)

            WCESolve(physfield[0], invgamma);

            Vmath::Vcopy(nq, physfield[0], 1, outarray[3], 1); // store z

            //------------------------------------

            //------------------------------------
            // Return to the primary variables (eq. 26c)

            // Compute gamma \nabla z
            Vmath::Smul(nq, gamma, physfield[0], 1, physfield[0], 1);

            m_fields[0]->IProductWRTDerivBase(0, physfield[0], coeffsfield[0]);
            m_fields[0]->IProductWRTDerivBase(1, physfield[0], coeffsfield[1]);

            Vmath::Neg(ncoeffs, coeffsfield[0], 1);
            Vmath::Neg(ncoeffs, coeffsfield[1], 1);

            // Set boundary conditions
            SetBoundaryConditionsContVariables(physfield[0], time);

            // Evaluate upwind numerical flux (physical space)
            NumericalFluxConsVariables(physfield[0], upwindX, upwindY);

            Vmath::Vmul(nTraceNumPoints, upwindX, 1, m_traceNormals[0], 1,
                        normflux, 1);
            m_fields[0]->AddTraceIntegral(normflux, coeffsfield[0]);
            Vmath::Vmul(nTraceNumPoints, upwindY, 1, m_traceNormals[1], 1,
                        normflux, 1);
            m_fields[0]->AddTraceIntegral(normflux, coeffsfield[1]);

            // Solve the remaining block-diagonal systems
            Array<OneD, Array<OneD, NekDouble>> modarray(2);
            for (int i = 0; i < 2; ++i)
            {
                modarray[i] = Array<OneD, NekDouble>(ncoeffs);
            }
            m_fields[0]->IProductWRTBase(outarray[1], modarray[0]);
            m_fields[0]->IProductWRTBase(outarray[2], modarray[1]);

            Vmath::Vadd(ncoeffs, modarray[0], 1, coeffsfield[0], 1, modarray[0],
                        1);
            Vmath::Vadd(ncoeffs, modarray[1], 1, coeffsfield[1], 1, modarray[1],
                        1);

            m_fields[0]->MultiplyByElmtInvMass(modarray[0], modarray[0]);
            m_fields[0]->MultiplyByElmtInvMass(modarray[1], modarray[1]);

            m_fields[0]->BwdTrans(modarray[0], outarray[1]);
            m_fields[0]->BwdTrans(modarray[1], outarray[2]);

            //------------------------------------

            break;
        }
        case MultiRegions::eGalerkin:
            ASSERTL0(false, "Unknown projection scheme for the Peregrine");
            break;
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }
}

void NonlinearPeregrine::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    NonlinearSWE::v_GenerateSummary(s);
}

// initial condition Laitone's first order solitary wave
void NonlinearPeregrine::LaitoneSolitaryWave(NekDouble amp, NekDouble d,
                                             NekDouble time, NekDouble x_offset)
{
    int nq = GetTotPoints();

    NekDouble A = 1.0;
    NekDouble C = sqrt(m_g * d) * (1.0 + 0.5 * (amp / d));

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> zeros(nq, 0.0);

    // get the coordinates (assuming all fields have the same
    // discretisation)
    m_fields[0]->GetCoords(x0, x1);

    for (int i = 0; i < nq; i++)
    {
        (m_fields[0]->UpdatePhys())[i] =
            amp * pow((1.0 / cosh(sqrt(0.75 * (amp / (d * d * d))) *
                                  (A * (x0[i] + x_offset) - C * time))),
                      2.0);
        (m_fields[1]->UpdatePhys())[i] =
            (amp / d) *
            pow((1.0 / cosh(sqrt(0.75 * (amp / (d * d * d))) *
                            (A * (x0[i] + x_offset) - C * time))),
                2.0) *
            sqrt(m_g * d);
    }

    Vmath::Sadd(nq, d, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);
    Vmath::Vmul(nq, m_fields[0]->GetPhys(), 1, m_fields[1]->GetPhys(), 1,
                m_fields[1]->UpdatePhys(), 1);
    Vmath::Vcopy(nq, zeros, 1, m_fields[2]->UpdatePhys(), 1);
    Vmath::Vcopy(nq, zeros, 1, m_fields[3]->UpdatePhys(), 1);

    // Forward transform to fill the coefficient space
    for (int i = 0; i < 4; ++i)
    {
        m_fields[i]->SetPhysState(true);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }
}

void NonlinearPeregrine::WCESolve(Array<OneD, NekDouble> &fce, NekDouble lambda)
{
    // As of now we need not to specify any BC routine for the WCE: periodic
    // and zero Neumann (for walls)

    // Note: this is just valid for the constant depth case:

    m_factors[StdRegions::eFactorLambda] = lambda;

    m_fields[3]->HelmSolve(fce, m_fields[3]->UpdateCoeffs(), m_factors);

    m_fields[3]->BwdTrans(m_fields[3]->GetCoeffs(), m_fields[3]->UpdatePhys());

    Vmath::Vcopy(fce.size(), m_fields[3]->GetPhys(), 1, fce, 1);
}

void NonlinearPeregrine::SetBoundaryConditionsForcing(
    Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] NekDouble time)
{
    int cnt = 0;

    // Loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),
                           "Wall"))
        {
            WallBoundaryForcing(n, cnt, inarray);
        }

        // Timedependent Boundary Condition
        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
        {
            ASSERTL0(false, "time-dependent BC not implemented for Boussinesq");
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

void NonlinearPeregrine::WallBoundaryForcing(
    int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &inarray)
{
    int nTraceNumPoints = GetTraceTotPoints();

    // Get physical values of f1 and f2 for the forward trace
    Array<OneD, Array<OneD, NekDouble>> Fwd(2);
    for (int i = 0; i < 2; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int id1, id2, npts;
    MultiRegions::ExpListSharedPtr bcexp =
        m_fields[0]->GetBndCondExpansions()[bcRegion];

    for (int e = 0; e < bcexp->GetExpSize(); ++e)
    {
        npts = bcexp->GetExp(e)->GetTotPoints();
        id1  = bcexp->GetPhys_Offset(e);
        id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
            m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        switch (m_expdim)
        {
            case 1:
            {
                ASSERTL0(false, "1D not yet implemented for Boussinesq");
                break;
            }
            case 2:
            {
                Array<OneD, NekDouble> tmp_n(npts);
                Array<OneD, NekDouble> tmp_t(npts);

                Vmath::Vmul(npts, &Fwd[0][id2], 1, &m_traceNormals[0][id2], 1,
                            &tmp_n[0], 1);
                Vmath::Vvtvp(npts, &Fwd[1][id2], 1, &m_traceNormals[1][id2], 1,
                             &tmp_n[0], 1, &tmp_n[0], 1);

                Vmath::Vmul(npts, &Fwd[0][id2], 1, &m_traceNormals[1][id2], 1,
                            &tmp_t[0], 1);
                Vmath::Vvtvm(npts, &Fwd[1][id2], 1, &m_traceNormals[0][id2], 1,
                             &tmp_t[0], 1, &tmp_t[0], 1);

                // Negate the normal flux
                Vmath::Neg(npts, tmp_n, 1);

                // Rotate back to Cartesian
                Vmath::Vmul(npts, &tmp_t[0], 1, &m_traceNormals[1][id2], 1,
                            &Fwd[0][id2], 1);
                Vmath::Vvtvm(npts, &tmp_n[0], 1, &m_traceNormals[0][id2], 1,
                             &Fwd[0][id2], 1, &Fwd[0][id2], 1);

                Vmath::Vmul(npts, &tmp_t[0], 1, &m_traceNormals[0][id2], 1,
                            &Fwd[1][id2], 1);
                Vmath::Vvtvp(npts, &tmp_n[0], 1, &m_traceNormals[1][id2], 1,
                             &Fwd[1][id2], 1, &Fwd[1][id2], 1);
                break;
            }
            case 3:
                ASSERTL0(false, "3D not implemented for Boussinesq equations");
                break;
            default:
                ASSERTL0(false, "Illegal expansion dimension");
        }

        // Copy boundary adjusted values into the boundary expansion
        bcexp = m_fields[1]->GetBndCondExpansions()[bcRegion];
        Vmath::Vcopy(npts, &Fwd[0][id2], 1, &(bcexp->UpdatePhys())[id1], 1);

        bcexp = m_fields[2]->GetBndCondExpansions()[bcRegion];
        Vmath::Vcopy(npts, &Fwd[1][id2], 1, &(bcexp->UpdatePhys())[id1], 1);
    }
}

void NonlinearPeregrine::NumericalFluxForcing(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, NekDouble> &numfluxX, Array<OneD, NekDouble> &numfluxY)
{
    int nTraceNumPoints = GetTraceTotPoints();

    //-----------------------------------------------------
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(2);
    Array<OneD, Array<OneD, NekDouble>> Bwd(2);

    for (int i = 0; i < 2; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
    }
    //-----------------------------------------------------

    //-----------------------------------------------------
    // get the physical values at the trace
    // (any time-dependent BC previously put in fields[1] and [2]

    m_fields[1]->GetFwdBwdTracePhys(inarray[0], Fwd[0], Bwd[0]);
    m_fields[2]->GetFwdBwdTracePhys(inarray[1], Fwd[1], Bwd[1]);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // use centred fluxes for the numerical flux
    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        numfluxX[i] = 0.5 * (Fwd[0][i] + Bwd[0][i]);
        numfluxY[i] = 0.5 * (Fwd[1][i] + Bwd[1][i]);
    }
    //-----------------------------------------------------
}

void NonlinearPeregrine::SetBoundaryConditionsContVariables(
    const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] NekDouble time)
{
    int cnt = 0;

    // Loop over Boundary Regions
    for (int n = 0; n < m_fields[3]->GetBndConditions().size(); ++n)
    {
        // Wall Boundary Condition
        if (boost::iequals(m_fields[3]->GetBndConditions()[n]->GetUserDefined(),
                           "Wall") ||
            m_fields[3]->GetBndConditions()[n]->IsTimeDependent())
        {
            WallBoundaryContVariables(n, cnt, inarray);
        }

        cnt += m_fields[3]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

void NonlinearPeregrine::WallBoundaryContVariables(
    int bcRegion, int cnt, const Array<OneD, const NekDouble> &inarray)
{
    int nTraceNumPoints = GetTraceTotPoints();

    // Get physical values of z for the forward trace
    Array<OneD, NekDouble> z(nTraceNumPoints);
    m_fields[3]->ExtractTracePhys(inarray, z);

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int id1, id2, npts;
    MultiRegions::ExpListSharedPtr bcexp =
        m_fields[3]->GetBndCondExpansions()[bcRegion];

    for (int e = 0; e < bcexp->GetExpSize(); ++e)
    {
        npts = bcexp->GetExp(e)->GetTotPoints();
        id1  = bcexp->GetPhys_Offset(e);
        id2  = m_fields[3]->GetTrace()->GetPhys_Offset(
            m_fields[3]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        // Copy boundary adjusted values into the boundary expansion
        // field[3]
        bcexp = m_fields[3]->GetBndCondExpansions()[bcRegion];
        Vmath::Vcopy(npts, &z[id2], 1, &(bcexp->UpdatePhys())[id1], 1);
    }
}

void NonlinearPeregrine::NumericalFluxConsVariables(
    const Array<OneD, const NekDouble> &physfield, Array<OneD, NekDouble> &outX,
    Array<OneD, NekDouble> &outY)
{
    int nTraceNumPoints = GetTraceTotPoints();

    //-----------------------------------------------------
    // get temporary arrays
    Array<OneD, NekDouble> Fwd(nTraceNumPoints);
    Array<OneD, NekDouble> Bwd(nTraceNumPoints);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // get the physical values at the trace
    // (we have put any time-dependent BC in field[3])

    m_fields[3]->GetFwdBwdTracePhys(physfield, Fwd, Bwd);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // use centred fluxes for the numerical flux
    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        outX[i] = 0.5 * (Fwd[i] + Bwd[i]);
        outY[i] = 0.5 * (Fwd[i] + Bwd[i]);
    }
    //-----------------------------------------------------
}

} // namespace Nektar
