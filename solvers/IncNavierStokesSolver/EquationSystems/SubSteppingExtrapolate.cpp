///////////////////////////////////////////////////////////////////////////////
//
// File: SubSteppingExtrapolate.cpp
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
// Description: Abstract base class for SubSteppingExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/SubSteppingExtrapolate.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/DisContField.h>

namespace Nektar
{

using namespace LibUtilities;

/**
 * Registers the class with the Factory.
 */
std::string SubSteppingExtrapolate::className =
    GetExtrapolateFactory().RegisterCreatorFunction(
        "SubStepping", SubSteppingExtrapolate::create, "SubStepping");

SubSteppingExtrapolate::SubSteppingExtrapolate(
    const LibUtilities::SessionReaderSharedPtr pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    MultiRegions::ExpListSharedPtr pPressure, const Array<OneD, int> pVel,
    const SolverUtils::AdvectionSharedPtr advObject)
    : Extrapolate(pSession, pFields, pPressure, pVel, advObject)
{
    m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
    m_session->LoadParameter("SubStepCFL", m_cflSafetyFactor, 0.5);
    m_session->LoadParameter("MinSubSteps", m_minsubsteps, 1);
    m_session->LoadParameter("MaxSubSteps", m_maxsubsteps, 100);

    size_t dim     = m_fields[0]->GetCoordim(0);
    m_traceNormals = Array<OneD, Array<OneD, NekDouble>>(dim);

    size_t nTracePts = m_fields[0]->GetTrace()->GetNpoints();
    for (size_t i = 0; i < dim; ++i)
    {
        m_traceNormals[i] = Array<OneD, NekDouble>(nTracePts);
    }
    m_fields[0]->GetTrace()->GetNormals(m_traceNormals);
}

void SubSteppingExtrapolate::v_EvaluatePressureBCs(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &fields,
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &N,
    [[maybe_unused]] NekDouble kinvis)
{
    ASSERTL0(false, "This method should not be called by Substepping routine");
}

void SubSteppingExtrapolate::v_SubSteppingTimeIntegration(
    const LibUtilities::TimeIntegrationSchemeSharedPtr &IntegrationScheme)
{
    Nektar::LibUtilities::Timer timer;
    m_intScheme = IntegrationScheme;

    timer.Start();
    size_t order = IntegrationScheme->GetOrder();

    // Set to 1 for first step and it will then be increased in
    // time advance routines
    if ((IntegrationScheme->GetName() == "Euler" &&
         IntegrationScheme->GetVariant() == "Backward") ||

        (IntegrationScheme->GetName() == "BDFImplicit" &&
         (order == 1 || order == 2)))
    {
        // Note RK first order SSP is just Forward Euler.
        std::string vSubStepIntScheme        = "RungeKutta";
        std::string vSubStepIntSchemeVariant = "SSP";
        int vSubStepIntSchemeOrder           = order;

        if (m_session->DefinesSolverInfo("SubStepIntScheme"))
        {
            vSubStepIntScheme = m_session->GetSolverInfo("SubStepIntScheme");
            vSubStepIntSchemeVariant = "";
            vSubStepIntSchemeOrder   = order;
        }

        m_subStepIntegrationScheme =
            LibUtilities::GetTimeIntegrationSchemeFactory().CreateInstance(
                vSubStepIntScheme, vSubStepIntSchemeVariant,
                vSubStepIntSchemeOrder, std::vector<NekDouble>());

        size_t nvel = m_velocity.size();
        size_t ndim = order + 1;

        // Fields for linear/quadratic interpolation
        m_previousVelFields = Array<OneD, Array<OneD, NekDouble>>(ndim * nvel);
        int ntotpts         = m_fields[0]->GetTotPoints();
        m_previousVelFields[0] = Array<OneD, NekDouble>(ndim * nvel * ntotpts);

        for (size_t i = 1; i < ndim * nvel; ++i)
        {
            m_previousVelFields[i] = m_previousVelFields[i - 1] + ntotpts;
        }
        // Vn fields
        m_previousVnFields    = Array<OneD, Array<OneD, NekDouble>>(ndim);
        ntotpts               = m_fields[0]->GetTrace()->GetTotPoints();
        m_previousVnFields[0] = Array<OneD, NekDouble>(ndim * ntotpts);
        for (size_t i = 1; i < ndim; ++i)
        {
            m_previousVnFields[i] = m_previousVnFields[i - 1] + ntotpts;
        }
    }
    else
    {
        ASSERTL0(0, "Integration method not suitable: Options include "
                    "BackwardEuler or BDFImplicitOrder{1,2}");
    }

    m_intSteps = IntegrationScheme->GetNumIntegrationPhases();

    // set explicit time-integration class operators
    m_subStepIntegrationOps.DefineOdeRhs(
        &SubSteppingExtrapolate::SubStepAdvection, this);
    m_subStepIntegrationOps.DefineProjection(
        &SubSteppingExtrapolate::SubStepProjection, this);
    timer.Stop();
    timer.AccumulateRegion(
        "SubSteppingExtrapolate:v_SubSteppingTimeIntegration", 10);
}

/**
 * Explicit Advection terms used by SubStepAdvance time integration
 */
void SubSteppingExtrapolate::SubStepAdvection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    Nektar::LibUtilities::Timer timer, timer1;
    size_t i;
    size_t nVariables     = inarray.size();
    size_t nQuadraturePts = inarray[0].size();

    timer.Start();
    /// Get the number of coefficients
    size_t ncoeffs = m_fields[0]->GetNcoeffs();

    /// Define an auxiliary variable to compute the RHS
    Array<OneD, Array<OneD, NekDouble>> WeakAdv(nVariables);
    WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs * nVariables);
    for (i = 1; i < nVariables; ++i)
    {
        WeakAdv[i] = WeakAdv[i - 1] + ncoeffs;
    }

    Array<OneD, Array<OneD, NekDouble>> Velfields(m_velocity.size());

    Velfields[0] = Array<OneD, NekDouble>(nQuadraturePts * m_velocity.size());

    for (i = 1; i < m_velocity.size(); ++i)
    {
        Velfields[i] = Velfields[i - 1] + nQuadraturePts;
    }

    Array<OneD, NekDouble> Vn(m_fields[0]->GetTrace()->GetTotPoints());

    SubStepExtrapolateField(fmod(time, m_timestep), Velfields, Vn);

    for (auto &x : m_forcing)
    {
        x->PreApply(m_fields, Velfields, Velfields, time);
    }
    timer1.Start();
    m_advObject->Advect(m_velocity.size(), m_fields, Velfields, inarray,
                        outarray, time);
    timer1.Stop();
    timer1.AccumulateRegion("SubStepAdvection:Advect", 10);

    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, outarray, outarray, time);
    }

    for (i = 0; i < nVariables; ++i)
    {
        m_fields[i]->IProductWRTBase(outarray[i], WeakAdv[i]);
        // negation requried due to sign of DoAdvection term to be consistent
        Vmath::Neg(ncoeffs, WeakAdv[i], 1);
    }

    AddAdvectionPenaltyFlux(Vn, inarray, WeakAdv);

    /// Operations to compute the RHS
    for (i = 0; i < nVariables; ++i)
    {
        // Negate the RHS
        Vmath::Neg(ncoeffs, WeakAdv[i], 1);

        /// Multiply the flux by the inverse of the mass matrix
        m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);

        /// Store in outarray the physical values of the RHS
        m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
    }
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:SubStepAdvection", 10);
}

/**
 * Projection used by SubStepAdvance time integration
 */
void SubSteppingExtrapolate::SubStepProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time)
{
    Nektar::LibUtilities::Timer timer;
    ASSERTL1(inarray.size() == outarray.size(),
             "Inarray and outarray of different sizes ");

    timer.Start();
    for (size_t i = 0; i < inarray.size(); ++i)
    {
        Vmath::Vcopy(inarray[i].size(), inarray[i], 1, outarray[i], 1);
    }
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:SubStepProjection", 10);
}

/**
 *
 */
void SubSteppingExtrapolate::v_SubStepSetPressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] const NekDouble Aii_Dt, NekDouble kinvis)
{
    Nektar::LibUtilities::Timer timer;
    // int nConvectiveFields =m_fields.size()-1;
    Array<OneD, Array<OneD, NekDouble>> nullvelfields;

    timer.Start();
    m_pressureCalls++;

    // Calculate viscous BCs at current level and
    // put in m_pressureHBCs[0]
    CalcNeumannPressureBCs(m_previousVelFields, nullvelfields, kinvis);

    // Extrapolate to m_pressureHBCs to n+1
    ExtrapolateArray(m_pressureHBCs);

    // Add (phi,Du/Dt) term to m_presureHBC
    AddDuDt();

    // Copy m_pressureHBCs to m_PbndExp
    CopyPressureHBCsToPbndExp();

    // Evaluate High order outflow conditiosn if required.
    CalcOutflowBCs(inarray, kinvis);
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:SubStepSetPressureBCs", 10);
}

/**
 *    At the start, the newest value is stored in array[nlevels-1]
 *        and the previous values in the first positions
 *    At the end, the acceleration from BDF is stored in array[nlevels-1]
 *        and the storage has been updated to included the new value
 */
void SubSteppingExtrapolate::v_AccelerationBDF(
    Array<OneD, Array<OneD, NekDouble>> &array)
{
    Nektar::LibUtilities::Timer timer;
    int nlevels = array.size();
    int nPts    = array[0].size();

    timer.Start();
    if (nPts)
    {
        // Update array
        RollOver(array);

        // Calculate acceleration using Backward Differentiation Formula
        Array<OneD, NekDouble> accelerationTerm(nPts, 0.0);

        int acc_order = std::min(m_pressureCalls, m_intSteps);
        Vmath::Smul(nPts, StifflyStable_Gamma0_Coeffs[acc_order - 1], array[0],
                    1, accelerationTerm, 1);

        for (int i = 0; i < acc_order; i++)
        {
            Vmath::Svtvp(
                nPts, -1 * StifflyStable_Alpha_Coeffs[acc_order - 1][i],
                array[i + 1], 1, accelerationTerm, 1, accelerationTerm, 1);
        }
        array[nlevels - 1] = accelerationTerm;
    }
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:v_AccelerationBDF", 10);
}

/**
 * Save the current fields in \a m_previousVelFields and cycle out older
 * previous fields so that it can be extrapolated to the new substeps
 * of the next time step. Also extract the normal trace of the
 * velocity field into \a m_previousVnFields along the trace so that this
 * can also be extrapolated.
 */
void SubSteppingExtrapolate::v_SubStepSaveFields(const int nstep)
{
    Nektar::LibUtilities::Timer timer;
    size_t i, n;
    timer.Start();
    size_t nvel      = m_velocity.size();
    size_t npts      = m_fields[0]->GetTotPoints();
    size_t ntracepts = m_fields[0]->GetTrace()->GetTotPoints();

    // rotate fields
    size_t nblocks = m_previousVelFields.size() / nvel;
    Array<OneD, NekDouble> save;

    // rotate storage space
    for (n = 0; n < nvel; ++n)
    {
        save = m_previousVelFields[(nblocks - 1) * nvel + n];

        for (i = nblocks - 1; i > 0; --i)
        {
            m_previousVelFields[i * nvel + n] =
                m_previousVelFields[(i - 1) * nvel + n];
        }

        m_previousVelFields[n] = save;
    }

    save = m_previousVnFields[nblocks - 1];
    for (i = nblocks - 1; i > 0; --i)
    {
        m_previousVnFields[i] = m_previousVnFields[i - 1];
    }
    m_previousVnFields[0] = save;

    // Put previous field
    for (i = 0; i < nvel; ++i)
    {
        m_fields[m_velocity[i]]->BwdTrans(
            m_fields[m_velocity[i]]->GetCoeffs(),
            m_fields[m_velocity[i]]->UpdatePhys());
        Vmath::Vcopy(npts, m_fields[m_velocity[i]]->GetPhys(), 1,
                     m_previousVelFields[i], 1);
    }

    Array<OneD, NekDouble> Fwd(ntracepts);

    // calculate Vn
    m_fields[0]->ExtractTracePhys(m_previousVelFields[0], Fwd);
    Vmath::Vmul(ntracepts, m_traceNormals[0], 1, Fwd, 1, m_previousVnFields[0],
                1);
    for (i = 1; i < m_bnd_dim; ++i)
    {
        m_fields[0]->ExtractTracePhys(m_previousVelFields[i], Fwd);
        Vmath::Vvtvp(ntracepts, m_traceNormals[i], 1, Fwd, 1,
                     m_previousVnFields[0], 1, m_previousVnFields[0], 1);
    }

    if (nstep == 0) // initialise all levels with first field
    {
        for (n = 0; n < nvel; ++n)
        {
            for (i = 1; i < nblocks; ++i)
            {
                Vmath::Vcopy(npts, m_fields[m_velocity[n]]->GetPhys(), 1,
                             m_previousVelFields[i * nvel + n], 1);
            }
        }

        for (i = 1; i < nblocks; ++i)
        {
            Vmath::Vcopy(ntracepts, m_previousVnFields[0], 1,
                         m_previousVnFields[i], 1);
        }
    }
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:v_SubStepSaveFields", 10);
}

/**
 *
 */
void SubSteppingExtrapolate::v_SubStepAdvance(int nstep, NekDouble time)
{
    int n;
    int nsubsteps;
    Nektar::LibUtilities::Timer timer;

    NekDouble dt;

    Array<OneD, Array<OneD, NekDouble>> fields;

    timer.Start();
    static int ncalls = 1;
    size_t nint       = std::min(ncalls++, m_intSteps);

    // this needs to change
    m_comm = m_fields[0]->GetComm()->GetRowComm();

    // Get the proper time step with CFL control
    dt = GetSubstepTimeStep();

    nsubsteps = (m_timestep > dt) ? ((int)(m_timestep / dt) + 1) : 1;
    nsubsteps = std::max(m_minsubsteps, nsubsteps);

    ASSERTL0(nsubsteps < m_maxsubsteps,
             "Number of substeps has exceeded maximum");

    dt = m_timestep / nsubsteps;

    if (m_infosteps && !((nstep + 1) % m_infosteps) && m_comm->GetRank() == 0)
    {
        std::cout << "Sub-integrating using " << nsubsteps
                  << " steps over Dt = " << m_timestep
                  << " (SubStep CFL=" << m_cflSafetyFactor << ")" << std::endl;
    }

    const TripleArray &solutionVector = m_intScheme->GetSolutionVector();

    for (size_t m = 0; m < nint; ++m)
    {
        // We need to update the fields held by the m_intScheme
        fields = solutionVector[m];

        // Initialise NS solver which is set up to use a GLM method
        // with calls to EvaluateAdvection_SetPressureBCs and
        // SolveUnsteadyStokesSystem
        m_subStepIntegrationScheme->InitializeScheme(dt, fields, time,
                                                     m_subStepIntegrationOps);

        for (n = 0; n < nsubsteps; ++n)
        {
            fields = m_subStepIntegrationScheme->TimeIntegrate(n, dt);
        }

        // set up HBC m_acceleration field for Pressure BCs
        IProductNormVelocityOnHBC(fields, m_iprodnormvel[m]);

        // Reset time integrated solution in m_intScheme
        m_intScheme->SetSolutionVector(m, fields);
    }
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:v_SubStepAdvance", 10);
}

/**
 *
 */
NekDouble SubSteppingExtrapolate::GetSubstepTimeStep()
{
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    size_t n_element = m_fields[0]->GetExpSize();

    const Array<OneD, int> ExpOrder = m_fields[0]->EvalBasisNumModesMaxPerExp();

    const NekDouble cLambda = 0.2; // Spencer book pag. 317

    Array<OneD, NekDouble> tstep(n_element, 0.0);
    Array<OneD, NekDouble> stdVelocity(n_element, 0.0);
    Array<OneD, Array<OneD, NekDouble>> velfields(m_velocity.size());

    for (size_t i = 0; i < m_velocity.size(); ++i)
    {
        velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
    }
    stdVelocity = GetMaxStdVelocity(velfields);

    for (size_t el = 0; el < n_element; ++el)
    {
        tstep[el] =
            m_cflSafetyFactor / (stdVelocity[el] * cLambda *
                                 (ExpOrder[el] - 1) * (ExpOrder[el] - 1));
    }

    NekDouble TimeStep = Vmath::Vmin(n_element, tstep, 1);
    m_comm->AllReduce(TimeStep, LibUtilities::ReduceMin);
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:GetSubStepTimeStep", 10);

    return TimeStep;
}

/**
 * Add the advection penalty term \f$ (\hat{u} - u^e)V_n \f$ given the
 * normal velocity \a Vn at this time level and the \a physfield values
 * containing the velocity field at this time level
 */

void SubSteppingExtrapolate::AddAdvectionPenaltyFlux(
    const Array<OneD, NekDouble> &Vn,
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &Outarray)
{
    ASSERTL1(physfield.size() == Outarray.size(),
             "Physfield and outarray are of different dimensions");

    size_t i;
    Nektar::LibUtilities::Timer timer;
    timer.Start();

    /// Number of trace points
    size_t nTracePts = m_fields[0]->GetTrace()->GetNpoints();

    /// Forward state array
    Array<OneD, NekDouble> Fwd(3 * nTracePts);

    /// Backward state array
    Array<OneD, NekDouble> Bwd = Fwd + nTracePts;

    /// upwind numerical flux state array
    Array<OneD, NekDouble> numflux = Bwd + nTracePts;

    for (i = 0; i < physfield.size(); ++i)
    {
        /// Extract forwards/backwards trace spaces
        /// Note it is important to use the zeroth field but with the
        /// specialised method to use boudnary conditions from other
        /// fields since trace spaces may not be the same if there are
        /// mixed boundary conditions
        std::dynamic_pointer_cast<MultiRegions::DisContField>(m_fields[0])
            ->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd,
                                 m_fields[i]->GetBndConditions(),
                                 m_fields[i]->GetBndCondExpansions());

        /// Upwind between elements
        m_fields[0]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux);

        /// Construct difference between numflux and Fwd,Bwd
        Vmath::Vsub(nTracePts, numflux, 1, Fwd, 1, Fwd, 1);
        Vmath::Vsub(nTracePts, numflux, 1, Bwd, 1, Bwd, 1);

        /// Calculate the numerical fluxes multipling Fwd, Bwd and
        /// numflux by the normal advection velocity
        Vmath::Vmul(nTracePts, Fwd, 1, Vn, 1, Fwd, 1);
        Vmath::Vmul(nTracePts, Bwd, 1, Vn, 1, Bwd, 1);

        m_fields[0]->AddFwdBwdTraceIntegral(Fwd, Bwd, Outarray[i]);
    }
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:AddAdvectionPenaltyFlux",
                           10);
}

/**
 * Extrapolate field using equally time spaced field un,un-1,un-2, (at
 * dt intervals) to time for substep n+t at order \a m_intSteps. Also
 * extrapolate the normal velocity along the trace at the same order
 */
void SubSteppingExtrapolate::SubStepExtrapolateField(
    NekDouble toff, Array<OneD, Array<OneD, NekDouble>> &ExtVel,
    Array<OneD, NekDouble> &ExtVn)
{
    size_t npts      = m_fields[0]->GetTotPoints();
    size_t ntracepts = m_fields[0]->GetTrace()->GetTotPoints();
    size_t nvel      = m_velocity.size();
    size_t i, j;
    Array<OneD, NekDouble> l(4);
    Nektar::LibUtilities::Timer timer;
    timer.Start();

    size_t ord = m_intSteps;

    // calculate Lagrange interpolants
    Vmath::Fill(4, 1.0, l, 1);

    for (i = 0; i <= ord; ++i)
    {
        for (j = 0; j <= ord; ++j)
        {
            if (i != j)
            {
                l[i] *= (j * m_timestep + toff);
                l[i] /= (j * m_timestep - i * m_timestep);
            }
        }
    }

    for (i = 0; i < nvel; ++i)
    {
        Vmath::Smul(npts, l[0], m_previousVelFields[i], 1, ExtVel[i], 1);

        for (j = 1; j <= ord; ++j)
        {
            Blas::Daxpy(npts, l[j], m_previousVelFields[j * nvel + i], 1,
                        ExtVel[i], 1);
        }
    }

    Vmath::Smul(ntracepts, l[0], m_previousVnFields[0], 1, ExtVn, 1);
    for (j = 1; j <= ord; ++j)
    {
        Blas::Daxpy(ntracepts, l[j], m_previousVnFields[j], 1, ExtVn, 1);
    }
    timer.Stop();
    timer.AccumulateRegion("SubSteppingExtrapolate:SubStepExtrapolateFields",
                           10);
}

/**
 *
 */
void SubSteppingExtrapolate::v_MountHOPBCs(
    int HBCdata, NekDouble kinvis, Array<OneD, NekDouble> &Q,
    [[maybe_unused]] Array<OneD, const NekDouble> &Advection)
{
    Vmath::Smul(HBCdata, -kinvis, Q, 1, Q, 1);
}

std::string SubSteppingExtrapolate::v_GetSubStepName(void)
{
    return m_subStepIntegrationScheme->GetName();
}

} // namespace Nektar
