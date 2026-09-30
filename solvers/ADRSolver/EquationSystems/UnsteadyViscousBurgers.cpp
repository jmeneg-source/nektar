///////////////////////////////////////////////////////////////////////////////
//
// File: UnsteadyViscousBurgers.cpp
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
// Description: Unsteady viscous Burgers solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/UnsteadyViscousBurgers.h>

namespace Nektar
{

std::string UnsteadyViscousBurgers::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "UnsteadyViscousBurgers", UnsteadyViscousBurgers::create,
        "Viscous Burgers equation");

UnsteadyViscousBurgers::UnsteadyViscousBurgers(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph),
      UnsteadyInviscidBurgers(pSession, pGraph)
{
}

/**
 * @brief Initialisation object for the viscous Burgers equation.
 */
void UnsteadyViscousBurgers::v_InitObject(bool DeclareFields)
{
    UnsteadyInviscidBurgers::v_InitObject(DeclareFields);

    m_session->LoadParameter("epsilon", m_epsilon, 1.0);
    m_session->MatchSolverInfo("SpectralVanishingViscosity", "True",
                               m_useSpecVanVisc, false);
    m_session->MatchSolverInfo("SpectralVanishingViscosity", "VarDiff",
                               m_useSpecVanViscVarDiff, false);

    if (m_useSpecVanViscVarDiff)
    {
        m_useSpecVanVisc = true;
    }

    if (m_useSpecVanVisc)
    {
        m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
        m_session->LoadParameter("SVVDiffCoeff", m_sVVDiffCoeff, 0.1);
    }

    // Type of diffusion classes to be used
    switch (m_projectionType)
    {
        // Discontinuous field
        case MultiRegions::eDiscontinuous:
        {
            if (m_explicitDiffusion)
            {
                std::string diffName;
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
                    diffName, diffName);
                m_diffusion->SetFluxVector(
                    &UnsteadyViscousBurgers::GetFluxVectorDiff, this);
                m_diffusion->InitObject(m_session, m_fields);
            }
            break;
        }
        // Continuous field
        case MultiRegions::eGalerkin:
        {
            std::string advName;
            m_session->LoadSolverInfo("AdvectionType", advName,
                                      "NonConservative");
            if (advName.compare("WeakDG") == 0)
            {
                // Define the normal velocity fields
                if (m_fields[0]->GetTrace())
                {
                    m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
                }

                std::string riemName;
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver =
                    SolverUtils::GetRiemannSolverFactory().CreateInstance(
                        riemName, m_session);
                m_riemannSolver->SetScalar(
                    "Vn", &UnsteadyViscousBurgers::GetNormalVelocity, this);
                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject(m_session, m_fields);
            }

            // In case of Galerkin explicit diffusion gives an error
            if (m_explicitDiffusion)
            {
                ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
            }
            break;
        }
        default:
        {
            ASSERTL0(false, "Unsupported projection type.");
            break;
        }
    }

    m_ode.DefineImplicitSolve(&UnsteadyViscousBurgers::DoImplicitSolve, this);
    m_ode.DefineOdeRhs(&UnsteadyViscousBurgers::DoOdeRhs, this);
}

/**
 * @brief Compute the right-hand side for the viscous Burgers equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyViscousBurgers::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Number of fields (variables of the problem)
    int nVariables = inarray.size();

    // Number of solution points
    int nSolutionPts = GetNpoints();

    UnsteadyInviscidBurgers::DoOdeRhs(inarray, outarray, time);

    if (m_explicitDiffusion)
    {
        Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nVariables);

        for (int i = 0; i < nVariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
        }

        m_diffusion->Diffuse(nVariables, m_fields, inarray, outarrayDiff);

        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vadd(nSolutionPts, &outarrayDiff[i][0], 1, &outarray[i][0],
                        1, &outarray[i][0], 1);
        }
    }
}

/* @brief Compute the diffusion term implicitly.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 * @param lambda     Diffusion coefficient.
 */
void UnsteadyViscousBurgers::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    int nvariables = inarray.size();
    int nq         = m_fields[0]->GetNpoints();

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;

    if (m_useSpecVanVisc)
    {
        factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
        factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff / m_epsilon;
    }

    if (m_projectionType == MultiRegions::eDiscontinuous)
    {
        factors[StdRegions::eFactorTau] = 1.0;
    }

    Array<OneD, Array<OneD, NekDouble>> F(nvariables);
    for (int n = 0; n < nvariables; ++n)
    {
        F[n] = Array<OneD, NekDouble>(nq);
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: nabla^2 \hat{Y}
    // where \hat = modal coeffs
    for (int i = 0; i < nvariables; ++i)
    {
        // Multiply 1.0/timestep/lambda
        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1,
                    F[i], 1);
    }

    // Setting boundary conditions
    SetBoundaryConditions(time);

    if (m_useSpecVanViscVarDiff)
    {
        static int cnt = 0;

        if (cnt % 10 == 0)
        {
            Array<OneD, Array<OneD, NekDouble>> vel(m_fields.size());
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->ClearGlobalLinSysManager();
                vel[i] = m_fields[i]->UpdatePhys();
            }
            SVVVarDiffCoeff(vel, m_varCoeffLap);
        }
        ++cnt;
    }

    for (int i = 0; i < nvariables; ++i)
    {
        // Solve a system of equations with Helmholtz solver
        m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), factors,
                               m_varCoeffLap);

        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
    }
}

/**
 * @brief Return the flux vector for the diffusion part.
 *
 */
void UnsteadyViscousBurgers::GetFluxVectorDiff(
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor)
{
    unsigned int nDim              = qfield.size();
    unsigned int nConvectiveFields = qfield[0].size();
    unsigned int nPts              = qfield[0][0].size();
    for (unsigned int j = 0; j < nDim; ++j)
    {
        for (unsigned int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Smul(nPts, m_epsilon, qfield[j][i], 1, viscousTensor[j][i],
                        1);
        }
    }
}

void UnsteadyViscousBurgers::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    UnsteadyInviscidBurgers::v_GenerateSummary(s);
    if (m_useSpecVanVisc)
    {
        std::stringstream ss;
        ss << "SVV (cut off = " << m_sVVCutoffRatio
           << ", coeff = " << m_sVVDiffCoeff << ")";
        SolverUtils::AddSummaryItem(s, "Smoothing", ss.str());
    }
}

} // namespace Nektar
