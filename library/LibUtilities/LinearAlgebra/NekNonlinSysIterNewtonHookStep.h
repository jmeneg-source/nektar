///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewtonHookStep.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: NekNonlinSysIterNewtonHookStep definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_HOOKSTEP_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_HOOKSTEP_H

#include <LibUtilities/LinearAlgebra/NekLinSysIterGMRES.h>
#include <LibUtilities/LinearAlgebra/NekNonlinSysIterNewton.h>

namespace Nektar::LibUtilities
{

/**
 * Newton hook-step solver uses GMRES data to compute a trust-region
 * constrained step. During the GMRES solution the Arnoldi process constructs
 * a Krylov basis and an upper Hessenberg matrix such that the linear residual
 * is represented by the reduced model:
 *
 *     || beta e_1 - H y ||_2.
 *
 * The hook-step solver reuses this reduced problem to solve the constrained
 * trust-region problem
 *
 *     min || beta e_1 - H y ||_2,  subject to ||y||_2 <= Delta,
 *
 * and reconstructs the full-space correction from the Krylov basis.
 *
 * This implementation requires GMRES and uses GMRESHookData from the recently
 * completed GMRES solve.
 * Left-preconditioned GMRES is not supported.
 * The reconstruction assumes the exported GMRES basis represents
 * the same search space used by the linear solve.
 */
class NekNonlinSysIterNewtonHookStep : public NekNonlinSysIterNewton
{
public:
    friend class MemoryManager<NekNonlinSysIterNewtonHookStep>;

    LIB_UTILITIES_EXPORT static NekNonlinSysIterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekNonlinSysIterSharedPtr p =
            MemoryManager<NekNonlinSysIterNewtonHookStep>::AllocateSharedPtr(
                pSession, vComm, nDimen, pKey);
        p->InitObject();
        return p;
    }

    static std::string className;

    LIB_UTILITIES_EXPORT NekNonlinSysIterNewtonHookStep(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey);

    LIB_UTILITIES_EXPORT ~NekNonlinSysIterNewtonHookStep() override = default;

protected:
    // Initial, current and min and max trust radi
    NekDouble m_trustRadiusInit = 1.0e-2; // can be set via parameters
    NekDouble m_trustRadius     = 1.0e-2;
    NekDouble m_trustRadiusMin  = 1.0e-8; // can be set via parameters
    NekDouble m_trustRadiusMax  = 1.0e2;  // can be set via parameters

    // Acceptance measure rho
    // minimum rho to accept the trial update
    NekDouble m_rhoAccept = 1.0e-1;
    // If rho < m_rhoShrink will result in trust radius shrink
    NekDouble m_rhoShrink = 2.5e-1;
    // If rho > m_rhoGrow will result in trust radius growth
    NekDouble m_rhoGrow = 7.5e-1;

    // If rho < m_rhoShrink, shrink by m_shrinkFactor
    NekDouble m_shrinkFactor = 0.5;
    // If rho > m_rhoGrow, grow by m_growFactor
    NekDouble m_growFactor = 2.0;

    // in fallback mode, used to test if the residual grcrease is sufficient
    NekDouble m_meaningfulAredRel = 1.0e-3;
    NekDouble m_meaningfulAredAbs = 1.0e-12;

    // In fallback mode, the maximum number of shrink steps
    int m_maxTrustRegionSteps = 4;
    // In oportunistic mode, the maximum number of growth steps
    int m_maxTrustRegionGrowthSteps = 3;
    // Limit the residuum growth befor failing
    NekDouble m_maxFallbackResidualGrowth = 10.0;

    // temporary trial vectors
    Array<OneD, NekDouble> m_solutionTrial;
    Array<OneD, NekDouble> m_residualTrial;
    Array<OneD, NekDouble> m_deltaTrial;

    // temporary exploration vecotrs
    Array<OneD, NekDouble> m_solutionBest;
    Array<OneD, NekDouble> m_residualBest;

    void v_InitObject() override;

    bool v_ApplyNewtonUpdate(const int ntotal,
                             const NekDouble oldResNorm) override;

private:
    NekDouble CalcL2Norm(const int ntotal,
                         const Array<OneD, const NekDouble> &x) const;

    std::vector<NekDouble> SolveDenseLinearSystem(
        std::vector<std::vector<NekDouble>> A, std::vector<NekDouble> b) const;

    Array<OneD, NekDouble> SolveReducedHookStep(
        const NekLinSysIterGMRES::GMRESHookData &data,
        const NekDouble Delta) const;

    void BuildHookStepFromKrylovBasis(
        const NekLinSysIterGMRES::GMRESHookData &data,
        const Array<OneD, const NekDouble> &y, Array<OneD, NekDouble> &delta);

    NekDouble ComputeReducedResidualNorm(
        const NekLinSysIterGMRES::GMRESHookData &data,
        const Array<OneD, const NekDouble> &y) const;

    bool HasMeaningfulDecrease(const NekDouble oldResNorm,
                               const NekDouble newResNorm) const;
};

} // namespace Nektar::LibUtilities

#endif
