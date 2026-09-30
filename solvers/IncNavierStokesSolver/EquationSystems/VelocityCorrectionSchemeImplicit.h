///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionSchemeImplicit.h
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
// Description: Velocity Correction Scheme header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEMEIMPLICIT_H
#define NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEMEIMPLICIT_H

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>

namespace Nektar
{
class VCSImplicit : public VelocityCorrectionScheme
{
public:
    friend class MemoryManager<VCSImplicit>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<VCSImplicit>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

protected:
    /// 2D Array for Advection Velocities [dir][dof]
    Array<OneD, Array<OneD, NekDouble>> m_advection, m_AdvVel;
    /// 3D Array for extrapolated Advection Velocities [dir][time-levle][dof]
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_extVel;
    /// bool to identify implicit scheme ie advection velocity
    bool m_advectionVelocity;
    /// bool to identify advection operator
    bool m_implicitSkewSymAdvection;
    /// integer for advection velocity
    int m_intOrder;
    std::string m_convectiveType;
    /// Array checking whether GlobalLinSys needs to be unset
    Array<OneD, int32_t> m_unsetGlobalLinSys;

    VCSImplicit(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph);

    ~VCSImplicit() override = default;

    // Virtual functions
    void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    void v_SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble>> &fields,
        Array<OneD, Array<OneD, NekDouble>> &Forcing,
        const NekDouble aii_Dt) override;

    void v_SetUpViscousForcing(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &Forcing,
        const NekDouble aii_Dt) override;

    void v_SolvePressure(const Array<OneD, NekDouble> &Forcing) override;

    void v_SolveViscous(
        const Array<OneD, const Array<OneD, NekDouble>> &Forcing,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble aii_Dt) override;

    void v_DoInitialise(bool dumpInitialConditions) override;

    std::string v_GetExtrapolateStr(void) override
    {
        return "Implicit"; // Use ImplicitExtrapolate.cpp
    }

    void v_EvaluateAdvection_SetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble time) override;

    static std::string solverTypeLookupId;

    void AddImplicitSkewSymAdvection(StdRegions::VarCoeffMap varcoeffs,
                                     NekDouble aii_Dt);

    // Function to check whether UnsetGlobalLinSys has to be called
    // at i-th call of LinearAdvectionDiffusionReactionSolve in SolveViscous
    void CheckUnsetGlobalLinSys(Array<OneD, int32_t> &unsetGlobalLinSys);

private:
};

typedef std::shared_ptr<VCSImplicit> VCSImplicitSharedPtr;

} // namespace Nektar

#endif // NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEMEIMPLICIT_H
