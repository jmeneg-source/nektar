///////////////////////////////////////////////////////////////////////////////
//
// File: UnsteadyInviscidBurgers.h
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

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYINVISCIDBURGERS_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYINVISCIDBURGERS_H

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Forcing/Forcing.h>

namespace Nektar
{

class UnsteadyInviscidBurgers : public SolverUtils::AdvectionSystem
{
public:
    friend class MemoryManager<UnsteadyInviscidBurgers>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<UnsteadyInviscidBurgers>::AllocateSharedPtr(pSession,
                                                                      pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

protected:
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;
    Array<OneD, NekDouble> m_traceVn;
    /// Forcing terms
    std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

    UnsteadyInviscidBurgers(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph);

    ~UnsteadyInviscidBurgers() override = default;

    /// Initialise the object
    void v_InitObject(bool DeclareFields = true) override;

    /// Compute the RHS
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Compute the projection
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    /// Get the normal velocity
    Array<OneD, NekDouble> &GetNormalVelocity();

    /// Evaluate the flux at each solution point
    void GetFluxVector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                       Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

    /// Print Summary
    void v_GenerateSummary(SolverUtils::SummaryList &s) override;

private:
};

} // namespace Nektar

#endif
