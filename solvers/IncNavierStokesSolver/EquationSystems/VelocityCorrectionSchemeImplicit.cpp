///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionSchemeImplicit.cpp
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
// Description: Implicit Velocity Correction Scheme for the Incompressible
// Navier Stokes equations similar to Dong and Shen 2010
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionSchemeImplicit.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/ContField.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
using namespace MultiRegions;

string VCSImplicit::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "VCSImplicit", VCSImplicit::create);

string VCSImplicit::solverTypeLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("SolverType", "VCSImplicit",
                                                   eVCSImplicit);

/**
 * Constructor. No special instructions here.
 * v_DoInitialise sets the scheme up.
 *
 * \param
 * \param
 */
VCSImplicit::VCSImplicit(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph),
      VelocityCorrectionScheme(pSession, pGraph)
{
}

/**
 *
 */
void VCSImplicit::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    AdvectionSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Splitting Scheme",
                                "Implicit velocity correction");

    if (m_extrapolation->GetSubStepName().size())
    {
        SolverUtils::AddSummaryItem(s, "Substepping",
                                    m_extrapolation->GetSubStepName());
    }

    string dealias = m_homogen_dealiasing ? "Homogeneous1D" : "";
    if (m_specHP_dealiasing)
    {
        dealias += (dealias == "" ? "" : " + ") + string("spectral/hp");
    }
    if (dealias != "")
    {
        SolverUtils::AddSummaryItem(s, "Dealiasing", dealias);
    }

    string smoothing = m_useSpecVanVisc ? "spectral/hp" : "";
    if (smoothing != "")
    {
        if (m_svvVarDiffCoeff == NullNekDouble1DArray)
        {
            SolverUtils::AddSummaryItem(
                s, "Smoothing-SpecHP",
                "SVV (" + smoothing + " Exp Kernel(cut-off = " +
                    boost::lexical_cast<string>(m_sVVCutoffRatio) +
                    ", diff coeff = " +
                    boost::lexical_cast<string>(m_sVVDiffCoeff) + "))");
        }
        else
        {
            if (m_IsSVVPowerKernel)
            {
                SolverUtils::AddSummaryItem(
                    s, "Smoothing-SpecHP",
                    "SVV (" + smoothing + " Power Kernel (Power ratio =" +
                        boost::lexical_cast<string>(m_sVVCutoffRatio) +
                        ", diff coeff = " +
                        boost::lexical_cast<string>(m_sVVDiffCoeff) +
                        "*Uh/p))");
            }
            else
            {
                SolverUtils::AddSummaryItem(
                    s, "Smoothing-SpecHP",
                    "SVV (" + smoothing + " DG Kernel (diff coeff = " +
                        boost::lexical_cast<string>(m_sVVDiffCoeff) +
                        "*Uh/p))");
            }
        }
    }

    if (m_useHomo1DSpecVanVisc && (m_HomogeneousType == eHomogeneous1D))
    {
        SolverUtils::AddSummaryItem(
            s, "Smoothing-Homo1D",
            "SVV (Homogeneous1D - Exp Kernel(cut-off = " +
                boost::lexical_cast<string>(m_sVVCutoffRatioHomo1D) +
                ", diff coeff = " +
                boost::lexical_cast<string>(m_sVVDiffCoeffHomo1D) + "))");
    }

    if (m_useGJPStabilisation)
    {
        SolverUtils::AddSummaryItem(
            s, "GJP Stab. Impl.    ",
            m_session->GetSolverInfo("GJPStabilisation"));
        SolverUtils::AddSummaryItem(s, "GJP Stab. JumpScale", m_GJPJumpScale);

        if (boost::iequals(m_session->GetSolverInfo("GJPStabilisation"),
                           "Explicit"))
        {
            SolverUtils::AddSummaryItem(
                s, "GJP Normal Velocity",
                m_session->GetSolverInfo("GJPNormalVelocity"));
        }
    }
}

/**
 *
 */
void VCSImplicit::v_DoInitialise(bool dumpInitialConditions)
{
    VelocityCorrectionScheme::v_DoInitialise(dumpInitialConditions);

    // Get integration order for extrapolation
    m_intOrder = m_intScheme->GetOrder();

    // Initialise Advection velocity
    m_AdvVel = Array<OneD, Array<OneD, NekDouble>>(m_velocity.size());
    for (int i = 0; i < m_velocity.size(); i++)
    {
        m_AdvVel[i] = Array<OneD, NekDouble>(m_fields[i]->GetTotPoints(), 0.0);
    }

    // Initialise storage for explicit advection terms
    m_advection = Array<OneD, Array<OneD, NekDouble>>(m_nConvectiveFields);
    for (int i = 0; i < m_nConvectiveFields; i++)
    {
        m_advection[i] =
            Array<OneD, NekDouble>(m_fields[i]->GetTotPoints(), 0.0);
    }

    // Check advection velocity: Extrapolated or Updated scheme
    // Default to Extrapolated scheme
    m_session->MatchSolverInfo("AdvectionVelocity", "Extrapolated",
                               m_advectionVelocity, true);

    // Extrapolated velocity (for 2nd order scheme)
    // Only initiliase for Extrapolated scheme
    if (m_advectionVelocity)
    {
        m_extVel =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_velocity.size());
        for (int i = 0; i < m_velocity.size(); i++)
        {
            m_extVel[i] = Array<OneD, Array<OneD, NekDouble>>(m_intOrder);
            for (int j = 0; j < m_intOrder; j++)
            {
                m_extVel[i][j] =
                    Array<OneD, NekDouble>(m_fields[i]->GetTotPoints(), 0.0);
            }
        }
    }

    // Check whether multiple or single GlobalLinSys are
    // required and when to delete them
    // Initialise to always reset GlobalLinSys
    m_unsetGlobalLinSys = Array<OneD, int32_t>(m_velocity.size(), 1);
    CheckUnsetGlobalLinSys(m_unsetGlobalLinSys);

    // Check skew-symmetric advection operator
    // Default to Convective operator = false
    m_session->MatchSolverInfo("AdvectionOperator", "SkewSymmetric",
                               m_implicitSkewSymAdvection, false);
}

/**
 * Computes the forcing term for the pressure solve in
 * VCSImplicit::v_SolvePressure(). It uses the weak pressure forcing similar to
 * VCSWeakPressure (coefficient space).
 *
 * @param fields    Holds the BDF formula of previous solutions \f$\sum_q
 * \frac{\alpha_q}{\gamma} [u^{n-q}, v^{n-q}, \ldots]\f$.
 * @param Forcing   Array for the forcing term of the pressure solve. May
 * contain temporary values on input. On output, Forcing = \f$ \int_\Omega
 * \nabla \phi \cdot \left( \frac{1}{\Delta t} \sum_q \alpha_q \mathbf{u}^{n-q}
 * + \mathbf{f}^{n+1} - N(\mathbf{u})^n \right) \f$.
 *
 */
void VCSImplicit::v_SetUpPressureForcing(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    Array<OneD, Array<OneD, NekDouble>> &Forcing, const NekDouble aii_Dt)
{
    int ncoeffs = m_pressure->GetNcoeffs();
    int phystot = m_fields[0]->GetTotPoints();
    int nvel    = m_velocity.size();

    // Add 1/dt \sum_q \alpha_q u^{n-q} - N(u^n)
    for (int i = 0; i < nvel; i++)
    {
        Vmath::Svtvp(phystot, 1.0 / aii_Dt, fields[i], 1, m_advection[i], 1,
                     Forcing[i], 1);
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, Forcing, Forcing, m_time + m_timestep);
    }

    // Do Innerproduct with derivative base
    m_pressure->IProductWRTDerivBase(Forcing, Forcing[0]);

    // Negate to reverse negation in HelmSolve
    Vmath::Neg(ncoeffs, Forcing[0], 1);
}

/**
 * Computes the forcing term and advection velocity for ADR solves.
 *
 * @param inarray   Holds the BDF formula of previous solutions \f$\sum_q
 \frac{\alpha_q}{\gamma} [u^{n-q}, v^{n-q}, \ldots]\f$.
 * @param Forcing   Shared forcing array m_F. Used in pressure and velocity
 solve. On input, holds temporary values from pressure solve. On output: Forcing
 = \f$ - \frac{1}{\nu} \left( \frac{1}{\Delta t} \sum_q \alpha_q
 \mathbf{u}^{n-q} - \nabla p^{n+1} + \mathbf{f}^{n+1} \right) \f$.

 * @param aii_Dt    \f$\frac{\Delta t}{\gamma}\f$.
 *
 * Additionally, the advection velocity \f$ \tilde{\mathbf{u}}^{n+1} \f$ is
 computed and handed to the velocity solve via the array m_AdvVel. It is
 computed as \f$ \tilde{\mathbf{u}}^{n+1} = \sum_q \frac{\alpha_q}{\gamma}
 \mathbf{u}^{n-q} - \frac{\Delta t}{\gamma} \left[ \nabla p^{n+1} -
 N(\mathbf{u})^n - \nu \nabla \times \omega^n + \mathbf{f}^{n+1} \right] \f$.
 *
 */
void VCSImplicit::v_SetUpViscousForcing(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &Forcing, const NekDouble aii_Dt)
{
    int phystot = m_fields[0]->GetTotPoints();
    int nvel    = m_velocity.size();

    // Update (PhysSpace) pressure to n+1
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());

    // Compute gradient \nabla p^{n+1}
    if (nvel == 2)
    {
        m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[m_velocity[0]],
                              Forcing[m_velocity[1]]);
    }
    else
    {
        m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[m_velocity[0]],
                              Forcing[m_velocity[1]], Forcing[m_velocity[2]]);
    }

    // Negate pressure gradient: -1 * \nabla p
    for (int i = 0; i < nvel; ++i)
    {
        Vmath::Neg(phystot, Forcing[i], 1);
    }

    // Zero convective fields
    for (int i = nvel; i < m_nConvectiveFields; ++i)
    {
        Vmath::Zero(phystot, Forcing[i], 1);
    }

    // Add forcing terms: += f^{n+1}
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, Forcing, Forcing, m_time + m_timestep);
    }

    // Build Advection Velocity
    for (int i = 0; i < nvel; ++i)
    {
        // Extrapolated Advection Velocity = u^*n+1
        if (m_advectionVelocity)
        {
            Vmath::Smul(phystot, 1.0 / m_diffCoeff[i],
                        m_extVel[i][m_intOrder - 1], 1, m_AdvVel[i], 1);
        }
        // Updated Advection Velocity = u^*n+1 - (...)
        else
        {
            // Get velocity fields [u^n, v^n, ..] for curl of vorticity
            Array<OneD, Array<OneD, NekDouble>> velocity(nvel);
            for (int i = 0; i < nvel; i++)
            {
                velocity[i] = Array<OneD, NekDouble>(m_fields[i]->GetPhys());
            }
            // Curl of vorticity \nabla \times \nabla \times [u^n, v^n, ..]
            m_fields[0]->CurlCurl(velocity, velocity);

            // \frac{\Delta t}{\gamma}(-\nabla p + f) + inarray
            Vmath::Svtvp(phystot, aii_Dt, Forcing[i], 1, inarray[i], 1,
                         m_AdvVel[i], 1);
            // += \frac{\Delta t}{\gamma} -N(u)
            Vmath::Svtvp(phystot, aii_Dt, m_advection[i], 1, m_AdvVel[i], 1,
                         m_AdvVel[i], 1);
            // += -\frac{\Delta t}{\gamma} CurlCurl(u)
            Vmath::Svtvp(phystot, -aii_Dt * m_diffCoeff[i], velocity[i], 1,
                         m_AdvVel[i], 1, m_AdvVel[i], 1);
            // *= 1/\nu
            Vmath::Smul(phystot, 1.0 / m_diffCoeff[i], m_AdvVel[i], 1,
                        m_AdvVel[i], 1);
        }
    }

    // Build Forcing term
    for (int i = 0; i < nvel; ++i)
    {
        // \frac{\gamma}{\Delta t} * inarray - \nabla p + f
        Vmath::Svtvp(phystot, 1.0 / aii_Dt, inarray[i], 1, Forcing[i], 1,
                     Forcing[i], 1);

        // *= -1/kinvis
        Vmath::Smul(phystot, -1.0 / m_diffCoeff[i], Forcing[i], 1, Forcing[i],
                    1);
    }
}

/**
 * Solve pressure system via a Poisson problem with ContField::HelmSolve().
 * Uses coefficient space forcing and hence sets the
 * argument PhysSpaceForcing=false in ContFieldHelmSolve().
 *
 * @param Forcing See output description in
 * VCSImplicit::v_SetUpPressureForcing().
 *
 */
void VCSImplicit::v_SolvePressure(const Array<OneD, NekDouble> &Forcing)
{
    StdRegions::ConstFactorMap factors;
    // Setup coefficient for equation, Lambda = 0 ensures Poisson instead of
    // Helmholtz problem
    factors[StdRegions::eFactorLambda] = 0.0;

    // Solve Pressure Poisson Equation (with Weak Forcing)
    m_pressure->HelmSolve(Forcing, m_pressure->UpdateCoeffs(), factors,
                          StdRegions::NullVarCoeffMap,
                          StdRegions::NullVarFactorsMap, NullNekDouble1DArray,
                          false);

    // Add presure to outflow bc if using convective like BCs
    m_extrapolation->AddPressureToOutflowBCs(m_kinvis);
}

/**
 * Solve velocity system via the
 * ContField::LinearAdvectionDiffusionReactionSolve().
 *
 * @param Forcing Holds the forcing term for each velocity component. See output
 * description in VCSImplicit::v_SetUpViscousForcing().
 * @param inarray  Unused in this routine. Still, holds \f$\sum_q
 * \frac{\alpha_q}{\gamma} [u^{n-q}, v^{n-q}, \ldots] \f$.
 * @param outarray Holds vector of previous solution \f$ [u^n, v^n, \ldots] \f$.
 * Overwritten upon output with new solution [u^{n+1}, v^{n+1}, \ldots].
 *
 * Additionally, the advection velocity \f$ \tilde{\mathbf{u}}^{n+1} \f$ is
 * handed to the velocity solve via the array m_AdvVel. It is computed in
 * VCSImplicit::v_SetUpViscousForcing(). It is defined as \f$
 * \tilde{\mathbf{u}}^{n+1} = \sum_q \frac{\alpha_q}{\gamma} \mathbf{u}^{n-q} -
 * \frac{\Delta t}{\gamma} \left[ \nabla p^{n+1} - N(\mathbf{u})^n - \nu \nabla
 * \times \omega^n + \mathbf{f}^{n+1} \right] \f$.
 *
 */
void VCSImplicit::v_SolveViscous(
    const Array<OneD, const Array<OneD, NekDouble>> &Forcing,
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble aii_Dt)
{
    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;
    StdRegions::VarFactorsMap varfactors = StdRegions::NullVarFactorsMap;

    AppendSVVFactors(factors, varfactors);
    ComputeGJPNormalVelocity(inarray, varcoeffs);

    // Compute mass varcoeff for implicit Skew-symmetric operator
    AddImplicitSkewSymAdvection(varcoeffs, aii_Dt);

    // Set advection velocities
    StdRegions::VarCoeffType varcoefftypes[] = {StdRegions::eVarCoeffVelX,
                                                StdRegions::eVarCoeffVelY,
                                                StdRegions::eVarCoeffVelZ};
    for (int i = 0; i < m_velocity.size(); i++)
    {
        varcoeffs[varcoefftypes[i]] = m_AdvVel[i];
    }

    // Solve Advection-Diffusion-Reaction system
    for (int i = 0; i < m_nConvectiveFields; ++i)
    {
        // Add diffusion coefficient to GJP matrix operator (Implicit part)
        if (m_useGJPStabilisation)
        {
            factors[StdRegions::eFactorGJP] = m_GJPJumpScale / m_diffCoeff[i];
        }

        // \lambda = - /frac{\gamma}{\Delta t \nu}
        factors[StdRegions::eFactorLambda] = -1.0 / aii_Dt / m_diffCoeff[i];

        auto gkey = m_fields[i]->LinearAdvectionDiffusionReactionSolve(
            Forcing[i], m_fields[i]->UpdateCoeffs(), factors, varcoeffs,
            varfactors);

        // Nuke GlobalLinSys to avoid memory leak
        // The condition nukes, if:
        // - Remove after last velocity solve (w-velocity in 3D,
        // assumes same matrix for each velocity)
        // - Catches Null return from 3DH1D solve by checking matrix type
        // is a LinearADR matrix (variant)
        if (m_unsetGlobalLinSys[i] &&
            (gkey.GetMatrixType() ==
                 StdRegions::eLinearAdvectionDiffusionReaction ||
             gkey.GetMatrixType() ==
                 StdRegions::eLinearAdvectionDiffusionReactionGJP))
        {
            m_fields[i]->UnsetGlobalLinSys(gkey, true);
        }

        // Transform solution to PhysSpace
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
    }
}

/**
 * Explicit part of the method - HOPBCs based on VCSWeakPressure.
 * For the implicit scheme, we do not extrapolate the advection and rotational
 * diffusion terms, see also S. Dong and J. Shen (2010).
 *
 * @param inarray Holds the vector of previous solutions \f$ [u^n, v^n, \ldots]
 * \f$.
 * @param outarray In the implicit scheme, currently an empty array. TODO Don't
 * waste memory here.
 * @param time Holds the new time \f$ t^{n+1} = (n+1) \Delta t \f$.
 *
 */
void VCSImplicit::v_EvaluateAdvection_SetPressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Zero RHS derivatives, avoids undefined values
    for (int i = 0; i < m_velocity.size(); ++i)
    {
        Vmath::Zero(outarray[i].size(), outarray[i], 1);
    }

    // Calculate High-Order pressure boundary conditions
    LibUtilities::Timer timer;
    timer.Start();
    m_extrapolation->EvaluatePressureBCs(inarray, outarray, m_kinvis);
    timer.Stop();
    timer.AccumulateRegion("Pressure BCs");

    // Extrapolate advection velocity with StifflyStableBeta coefficients
    if (m_advectionVelocity)
    {
        for (int i = 0; i < m_velocity.size(); ++i)
        {
            Vmath::Vcopy(inarray[i].size(), inarray[i], 1,
                         m_extVel[i][m_intOrder - 1], 1);
            m_extrapolation->ExtrapolateArray(m_extVel[i]);
        }
    }

    // Compute explicit advection terms
    EvaluateAdvectionTerms(inarray, m_advection, time);
}

void VCSImplicit::AddImplicitSkewSymAdvection(StdRegions::VarCoeffMap varcoeffs,
                                              NekDouble aii_Dt)
{
    if (!m_implicitSkewSymAdvection)
    {
        return;
    }

    // Setup temporary storages
    int nvel    = m_velocity.size();
    int phystot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> divAdvVel(phystot, 0.0);
    Array<OneD, NekDouble> tmpstore(phystot);
    Array<OneD, NekDouble> ones(phystot, 1.0);

    // Compute divergence of advection velocity: 1/\nu \tilde{u}
    m_fields[0]->PhysDeriv(eX, m_AdvVel[0], divAdvVel);
    for (int i = 1; i < nvel; ++i)
    {
        m_fields[i]->PhysDeriv(DirCartesianMap[i], m_AdvVel[i], tmpstore);
        Vmath::Vadd(phystot, divAdvVel, 1, tmpstore, 1, divAdvVel, 1);
    }
    // Smul: 1/(2*lambda) * divAdvVel
    NekDouble lambda = 1.0 / aii_Dt / m_diffCoeff[0];
    Vmath::Smul(phystot, 1.0 / 2.0 / lambda, divAdvVel, 1, divAdvVel, 1);

    // Vadd: 1 + 1/(2*lambda) * divAdvVel
    Vmath::Vadd(phystot, ones, 1, divAdvVel, 1, divAdvVel, 1);

    // Assign VarcoeffMass
    varcoeffs[StdRegions::eVarCoeffMass] = divAdvVel;
}

void VCSImplicit::CheckUnsetGlobalLinSys(
    Array<OneD, int32_t> &unsetGlobalLinSys)
{
    // Loop over all convected Fields (excludes pressure)
    for (int i = 0; i < m_nConvectiveFields; ++i)
    {
        // Always unset after last solve
        if (i == m_nConvectiveFields - 1)
        {
            unsetGlobalLinSys[i] = 1;
            continue;
        }

        auto contField =
            std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[i]);
        // Check against all solves after this one, whether same
        // Global Linear system is required
        for (int j = i + 1; j < m_nConvectiveFields; ++j)
        {
            // If same BCs (ie GlobalLinSys) in j-solves, save the GlobalLinSys
            // for re-use
            if (std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[j])
                    ->SameTypeOfBoundaryConditions(*contField))
            {
                unsetGlobalLinSys[i] = 0;
            }
            // If no same BCs (ie GlobalLinSys) in j-solves, UnsetGlobalLinSys()
            // for i-solve
            else
            {
                unsetGlobalLinSys[i] = 1;
            }
        }
    }
}

} // namespace Nektar
