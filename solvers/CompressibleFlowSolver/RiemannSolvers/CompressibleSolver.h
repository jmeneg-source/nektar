///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleSolver.h
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
// Description: Compressible Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_COMPRESSIBLESOLVER
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_COMPRESSIBLESOLVER

#include <CompressibleFlowSolver/Misc/EquationOfState.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class CompressibleSolver : public RiemannSolver
{
protected:
    bool m_pointSolve;
    EquationOfStateSharedPtr m_eos;
    bool m_idealGas;

    /// Session ctor
    CompressibleSolver(const LibUtilities::SessionReaderSharedPtr &pSession);

    /// Programmatic ctor
    CompressibleSolver();

    using ND = NekDouble;

    void v_Solve(const int nDim, const Array<OneD, const Array<OneD, ND>> &Fwd,
                 const Array<OneD, const Array<OneD, ND>> &Bwd,
                 Array<OneD, Array<OneD, ND>> &flux) override;

    virtual void v_ArraySolve(
        [[maybe_unused]] const Array<OneD, const Array<OneD, ND>> &Fwd,
        [[maybe_unused]] const Array<OneD, const Array<OneD, ND>> &Bwd,
        [[maybe_unused]] Array<OneD, Array<OneD, ND>> &flux)
    {
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }

    virtual void v_PointSolve(
        [[maybe_unused]] ND rhoL, [[maybe_unused]] ND rhouL,
        [[maybe_unused]] ND rhovL, [[maybe_unused]] ND rhowL,
        [[maybe_unused]] ND EL, [[maybe_unused]] ND rhoR,
        [[maybe_unused]] ND rhouR, [[maybe_unused]] ND rhovR,
        [[maybe_unused]] ND rhowR, [[maybe_unused]] ND ER,
        [[maybe_unused]] ND &rhof, [[maybe_unused]] ND &rhouf,
        [[maybe_unused]] ND &rhovf, [[maybe_unused]] ND &rhowf,
        [[maybe_unused]] ND &Ef)
    {
        NEKERROR(ErrorUtil::efatal,
                 "This function should be defined by subclasses.");
    }

    ND GetRoeSoundSpeed(ND rhoL, ND pL, ND eL, ND HL, ND srL, ND rhoR, ND pR,
                        ND eR, ND HR, ND srR, ND HRoe, ND URoe2, ND srLR);
};

} // namespace Nektar

#endif
