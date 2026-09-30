///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewtonHookStep.cpp
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

#include <LibUtilities/LinearAlgebra/NekNonlinSysIterNewtonHookStep.h>
#include <algorithm>
#include <cmath>
#include <vector>

/*
<SOLVERINFO>
    <I PROPERTY="NonlinearSolver" VALUE="NewtonHookStep" />
    <I PROPERTY="LinSysIterSolver" VALUE="GMRES" />
</SOLVERINFO>

<PARAMETERS>
    <P> TrustRadiusInit          = 1.0e-2 </P>
    <P> TrustRadiusMin           = 1.0e-8 </P>
    <P> TrustRadiusMax           = 1.0e0  </P>

    <P> TrustRegionRhoAccept     = 1.0e-1 </P>
    <P> TrustRegionRhoShrink     = 2.5e-1 </P>
    <P> TrustRegionRhoGrow       = 7.5e-1 </P>
    <P> TrustRegionShrinkFactor  = 0.5    </P>
    <P> TrustRegionGrowFactor    = 2.0    </P>

    <P> MaxTrustRegionSteps      = 4      </P>
    <P> MaxTrustRegionGrowthSteps = 3 </P>
    <P> TrustRegionMaxFallbackResidualGrowth = 10.0 </P>

    <P> LinSysMaxStorage         = 1000   </P>

    <P> MeaningfulAredRel = 1.0e-3 </P>
    <P> MeaningfulAredAbs = 1.0e-12 </P>
</PARAMETERS>
*/

using namespace std;

namespace Nektar::LibUtilities
{
string NekNonlinSysIterNewtonHookStep::className =
    LibUtilities::GetNekNonlinSysIterFactory().RegisterCreatorFunction(
        "NewtonHookStep", NekNonlinSysIterNewtonHookStep::create,
        "Newton solver with GMRES-based hook step.");

NekNonlinSysIterNewtonHookStep::NekNonlinSysIterNewtonHookStep(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm, const int nDimen,
    const NekSysKey &pKey)
    : NekNonlinSysIterNewton(pSession, vComm, nDimen, pKey)
{
    pSession->LoadParameter("TrustRadiusInit", m_trustRadiusInit, 1.0e-2);
    pSession->LoadParameter("TrustRadiusMin", m_trustRadiusMin, 1.0e-8);
    pSession->LoadParameter("TrustRadiusMax", m_trustRadiusMax, 1.0e2);

    pSession->LoadParameter("TrustRegionRhoAccept", m_rhoAccept, 1.0e-1);
    pSession->LoadParameter("TrustRegionRhoShrink", m_rhoShrink, 2.5e-1);
    pSession->LoadParameter("TrustRegionRhoGrow", m_rhoGrow, 7.5e-1);
    pSession->LoadParameter("TrustRegionShrinkFactor", m_shrinkFactor, 0.5);
    pSession->LoadParameter("TrustRegionGrowFactor", m_growFactor, 2.0);

    pSession->LoadParameter("MaxTrustRegionSteps", m_maxTrustRegionSteps, 4);
    pSession->LoadParameter("MaxTrustRegionGrowthSteps",
                            m_maxTrustRegionGrowthSteps, 3);
    pSession->LoadParameter("TrustRegionMaxFallbackResidualGrowth",
                            m_maxFallbackResidualGrowth, 10.0);
    pSession->LoadParameter("MeaningfulAredRel", m_meaningfulAredRel, 1.0e-3);
    pSession->LoadParameter("MeaningfulAredAbs", m_meaningfulAredAbs, 1.0e-12);
}

void NekNonlinSysIterNewtonHookStep::v_InitObject()
{
    NekNonlinSysIterNewton::v_InitObject();

    m_solutionTrial = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_residualTrial = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_deltaTrial    = Array<OneD, NekDouble>(m_SysDimen, 0.0);

    m_solutionBest = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_residualBest = Array<OneD, NekDouble>(m_SysDimen, 0.0);

    m_trustRadius = m_trustRadiusInit;

    auto gmres = std::dynamic_pointer_cast<NekLinSysIterGMRES>(m_linsol);
    ASSERTL0(gmres, "NewtonHookStep requires LinSysIterSolver=GMRES.");
    gmres->EnableHookStepExport(true);
}

NekDouble NekNonlinSysIterNewtonHookStep::CalcL2Norm(
    const int ntotal, const Array<OneD, const NekDouble> &x) const
{
    NekDouble nrm2 = Vmath::Dot(ntotal, x, x);
    m_rowComm->AllReduce(nrm2, LibUtilities::ReduceSum);
    return sqrt(nrm2);
}

/*
    Solve a small dense linear system Ax=b for x
    using Gaussian elimination with pivoting
    The system to be solved is
    $$(H^T H + μI)y=\beta H^t e_1​$$, so a symmetric one
    could use Cholesky factorisation
*/
std::vector<NekDouble> NekNonlinSysIterNewtonHookStep::SolveDenseLinearSystem(
    std::vector<std::vector<NekDouble>> A, std::vector<NekDouble> b) const
{
    int n = static_cast<int>(b.size());

    for (int k = 0; k < n; ++k)
    {
        int piv        = k;
        NekDouble amax = std::abs(A[k][k]);
        for (int i = k + 1; i < n; ++i)
        {
            if (std::abs(A[i][k]) > amax)
            {
                amax = std::abs(A[i][k]);
                piv  = i;
            }
        }

        ASSERTL0(amax > 1.0e-14, "Hook-step dense solve breakdown.");

        if (piv != k)
        {
            std::swap(A[piv], A[k]);
            std::swap(b[piv], b[k]);
        }

        NekDouble diag = A[k][k];
        for (int j = k; j < n; ++j)
        {
            A[k][j] /= diag;
        }
        b[k] /= diag;

        for (int i = k + 1; i < n; ++i)
        {
            NekDouble fac = A[i][k];
            for (int j = k; j < n; ++j)
            {
                A[i][j] -= fac * A[k][j];
            }
            b[i] -= fac * b[k];
        }
    }

    std::vector<NekDouble> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j)
        {
            x[i] -= A[i][j] * x[j];
        }
    }

    return x;
}

/*
    Using GMRES data look for $\mu$ such thatm_TrustRadiusInit
    $(H^T H + μI)y=\beta H^t e_1​$ and $||y(\mu)<\Delta||$
*/
Array<OneD, NekDouble> NekNonlinSysIterNewtonHookStep::SolveReducedHookStep(
    const NekLinSysIterGMRES::GMRESHookData &data, const NekDouble Delta) const
{
    int m = data.nswp;
    Array<OneD, NekDouble> y(m, 0.0);

    std::vector<std::vector<NekDouble>> H(m + 1,
                                          std::vector<NekDouble>(m, 0.0));

    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < m + 1; ++i)
        {
            H[i][j] = data.H[j][i];
        }
    }

    std::vector<std::vector<NekDouble>> AtA(m, std::vector<NekDouble>(m, 0.0));
    std::vector<NekDouble> rhs(m, 0.0);

    for (int i = 0; i < m; ++i)
    {
        rhs[i] = H[0][i] * data.beta;
        for (int j = 0; j < m; ++j)
        {
            NekDouble sum = 0.0;
            for (int k = 0; k < m + 1; ++k)
            {
                sum += H[k][i] * H[k][j];
            }
            AtA[i][j] = sum;
        }
    }

    auto computeY = [&](NekDouble mu) {
        std::vector<std::vector<NekDouble>> A = AtA;
        for (int i = 0; i < m; ++i)
        {
            A[i][i] += mu;
        }
        return SolveDenseLinearSystem(A, rhs);
    };

    auto normY = [&](const std::vector<NekDouble> &yv) {
        NekDouble s = 0.0;
        for (auto v : yv)
        {
            s += v * v;
        }
        return sqrt(s);
    };

    std::vector<NekDouble> y0 = computeY(0.0);
    if (normY(y0) <= Delta)
    {
        for (int i = 0; i < m; ++i)
        {
            y[i] = y0[i];
        }
        return y;
    }

    NekDouble muL = 0.0;
    NekDouble muR = 1.0;

    while (normY(computeY(muR)) > Delta)
    {
        muR *= 10.0;
        ASSERTL0(muR < 1.0e20, "Unable to bracket hook-step multiplier.");
    }

    // For now just do 40 iterations, consider a stop cryterion tolerance
    for (int it = 0; it < 40; ++it)
    {
        NekDouble mu              = 0.5 * (muL + muR);
        std::vector<NekDouble> ym = computeY(mu);
        if (normY(ym) > Delta)
        {
            muL = mu;
        }
        else
        {
            muR = mu;
        }
    }

    std::vector<NekDouble> yv = computeY(muR);
    for (int i = 0; i < m; ++i)
    {
        y[i] = yv[i];
    }

    return y;
}

/*
Take computed y spanning the Krylov space and convert to the actual step delta.
*/
void NekNonlinSysIterNewtonHookStep::BuildHookStepFromKrylovBasis(
    const NekLinSysIterGMRES::GMRESHookData &data,
    const Array<OneD, const NekDouble> &y, Array<OneD, NekDouble> &delta)
{
    int nNonDir = data.nGlobal - data.nDir;
    Vmath::Zero(nNonDir, delta, 1);

    for (int i = 0; i < data.nswp; ++i)
    {
        // Can not mix Array and T* ...
        Vmath::Svtvp(nNonDir, y[i], &data.V[i][0] + data.nDir, 1, delta.data(),
                     1, delta.data(), 1);
    }

    auto gmres = std::dynamic_pointer_cast<NekLinSysIterGMRES>(m_linsol);
    if (gmres && gmres->UsingRightPrecon())
    {
        Array<OneD, NekDouble> tmp(nNonDir, 0.0);
        gmres->ApplyRightPrecon(delta, tmp);
        Vmath::Vcopy(nNonDir, tmp, 1, delta, 1);
    }
}

/*
Hook step is based on the quadratic model of the problem
This computes the predicted linearised residual norm in the reduced
GMRES/Krylov space
$$
\sqrt{2m(y)} = || \beta e_1 - H y ||_2
$$.
*/
NekDouble NekNonlinSysIterNewtonHookStep::ComputeReducedResidualNorm(
    const NekLinSysIterGMRES::GMRESHookData &data,
    const Array<OneD, const NekDouble> &y) const
{
    int m = data.nswp;
    Array<OneD, NekDouble> rred(m + 1, 0.0);
    rred[0] = data.beta;

    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < m + 1; ++i)
        {
            rred[i] -= data.H[j][i] * y[j];
        }
    }

    NekDouble nrm2 = 0.0;
    for (int i = 0; i < m + 1; ++i)
    {
        nrm2 += rred[i] * rred[i];
    }
    return sqrt(nrm2);
}

bool NekNonlinSysIterNewtonHookStep::v_ApplyNewtonUpdate(
    const int ntotal, const NekDouble oldResNorm)
{
    auto gmres = std::dynamic_pointer_cast<NekLinSysIterGMRES>(m_linsol);
    ASSERTL0(gmres, "NewtonHookStep requires GMRES.");

    ASSERTL0(!gmres->UsingLeftPrecon(),
             "Initial NewtonHookStep implementation does not support "
             "GMRESLeftPrecon.");

    const auto &data = gmres->GetHookData();
    ASSERTL0(data.valid, "GMRES hook-step data not available.");

    NekDouble trustRadiusTrial = m_trustRadius;
    cout << " Attempting to apply Newton update with Hookstep "
         << " Trust radius = " << trustRadiusTrial
         << " GMRES beta = " << data.beta << " oldResNorm= " << oldResNorm
         << endl;

    bool haveAcceptedCandidate = false;
    bool hadShrink             = false;

    NekDouble bestResNorm     = 0.0;
    NekDouble bestRho         = 0.0;
    NekDouble bestTrustRadius = 0.0;
    bool bestBoundaryStep     = false;

    // ---------------------------------------------------------------------
    // Phase 1: recovery by shrinking radius until the first acceptable step.
    // ---------------------------------------------------------------------
    for (int it = 0; it < m_maxTrustRegionSteps; ++it)
    {
        Array<OneD, NekDouble> yHook =
            SolveReducedHookStep(data, trustRadiusTrial);

        Vmath::Zero(ntotal, m_deltaTrial, 1);
        BuildHookStepFromKrylovBasis(data, yHook, m_deltaTrial);

        // predicted norm of the residual
        NekDouble predResNorm = ComputeReducedResidualNorm(data, yHook);

        Vmath::Vcopy(ntotal, m_Solution, 1, m_solutionTrial, 1);
        Vmath::Svtvp(ntotal, -1.0, m_deltaTrial, 1, m_solutionTrial, 1,
                     m_solutionTrial, 1);

        m_operator.DoNekSysResEval(m_solutionTrial, m_residualTrial);

        NekDouble newResNormSq =
            Vmath::Dot(ntotal, m_residualTrial, m_residualTrial);
        m_rowComm->AllReduce(newResNormSq, LibUtilities::ReduceSum);
        NekDouble newResNorm = sqrt(newResNormSq);

        NekDouble ared = oldResNorm - newResNorm;
        NekDouble pred = oldResNorm - predResNorm;
        // Protect from small pred
        NekDouble rho = (pred > 1.0e-14 * oldResNorm) ? ared / pred : -1.0;

        NekDouble yNorm   = CalcL2Norm(data.nswp, yHook);
        bool boundaryStep = std::abs(yNorm - trustRadiusTrial) <=
                            1.0e-2 * (trustRadiusTrial + 1.0e-14);

        cout << "  HookStep TrustRadius trial = " << trustRadiusTrial
             << " yNorm = " << yNorm << " oldRes = " << oldResNorm
             << " newRes = " << newResNorm << " predRes = " << predResNorm
             << " ared = " << ared << " pred = " << pred << " rho = " << rho
             << endl;

        bool meaningfulDecrease = HasMeaningfulDecrease(oldResNorm, newResNorm);

        if (rho > m_rhoAccept || meaningfulDecrease)
        {
            cout << "Accepting HookStep with rho = " << rho << " > "
                 << " m_RhoAccept = " << m_rhoAccept << endl;

            Vmath::Vcopy(ntotal, m_solutionTrial, 1, m_solutionBest, 1);
            Vmath::Vcopy(ntotal, m_residualTrial, 1, m_residualBest, 1);

            bestResNorm           = newResNorm;
            bestRho               = rho;
            bestTrustRadius       = trustRadiusTrial;
            bestBoundaryStep      = boundaryStep;
            haveAcceptedCandidate = true;

            break;
        }

        hadShrink = true;

        trustRadiusTrial =
            std::max(m_trustRadiusMin, m_shrinkFactor * trustRadiusTrial);

        cout << " HookStep shrink trial radius to " << trustRadiusTrial
             << " because current trial was not acceptable." << endl;

        if (trustRadiusTrial <= m_trustRadiusMin + 1.0e-16)
        {
            WARNINGL0(false, "Hook-step globalization failed at TrustRadiusMin "
                             "without meaningful residual decrease.");
            m_trustRadius = m_trustRadiusMin;
            return false;
        }
    }

    cout << "HookStep finished the shrink phase with: "
         << " hadShrink=" << hadShrink << " bestRho=" << bestRho
         << " bestBoundaryStep=" << bestBoundaryStep << endl;

    // ---------------------------------------------------------------------
    // No acceptable step found in Phase 1: fallback at TrustRadiusMin.
    // ---------------------------------------------------------------------
    if (!haveAcceptedCandidate)
    {
        Array<OneD, NekDouble> yFallback =
            SolveReducedHookStep(data, m_trustRadiusMin);

        Vmath::Zero(ntotal, m_deltaTrial, 1);
        BuildHookStepFromKrylovBasis(data, yFallback, m_deltaTrial);

        Vmath::Vcopy(ntotal, m_Solution, 1, m_solutionTrial, 1);
        Vmath::Svtvp(ntotal, -1.0, m_deltaTrial, 1, m_solutionTrial, 1,
                     m_solutionTrial, 1);

        m_operator.DoNekSysResEval(m_solutionTrial, m_residualTrial);

        NekDouble fallbackResNormSq =
            Vmath::Dot(ntotal, m_residualTrial, m_residualTrial);
        m_rowComm->AllReduce(fallbackResNormSq, LibUtilities::ReduceSum);
        NekDouble fallbackResNorm = sqrt(fallbackResNormSq);

        cout << "  No acceptable hook-step found."
             << " Taking fallback TrustRadiusMin = " << m_trustRadiusMin
             << " residual = " << fallbackResNorm << endl;

        if (!std::isfinite(fallbackResNorm) ||
            fallbackResNorm > m_maxFallbackResidualGrowth * oldResNorm)
        {
            WARNINGL0(false,
                      "Fallback hook-step caused non-finite or excessive "
                      "residual growth.");
            m_trustRadius = m_trustRadiusMin;
            return false;
        }

        Vmath::Vcopy(ntotal, m_solutionTrial, 1, m_Solution, 1);
        Vmath::Vcopy(ntotal, m_residualTrial, 1, m_Residual, 1);
        m_haveUpdatedResidual = true;

        m_trustRadius = m_trustRadiusMin;
        return true;
    }

    // ---------------------------------------------------------------------
    // Phase 2: exploration by enlarging radius, but only if Phase 1
    // succeeded without shrink and the accepted step is a strong
    // boundary-limited candidate.
    // ---------------------------------------------------------------------
    if (!hadShrink && bestRho > m_rhoGrow && bestBoundaryStep)
    {
        cout << " HookStep entering exploration phase from TrustRadius = "
             << bestTrustRadius << " because bestRho = " << bestRho
             << " > m_rhoGrow = " << m_rhoGrow << endl;

        NekDouble exploreRadius = bestTrustRadius;

        for (int ig = 0; ig < m_maxTrustRegionGrowthSteps; ++ig)
        {
            NekDouble enlargedRadius =
                std::min(m_trustRadiusMax, m_growFactor * exploreRadius);

            if (enlargedRadius <= exploreRadius + 1.0e-16)
            {
                break;
            }

            Array<OneD, NekDouble> yHook =
                SolveReducedHookStep(data, enlargedRadius);

            Vmath::Zero(ntotal, m_deltaTrial, 1);
            BuildHookStepFromKrylovBasis(data, yHook, m_deltaTrial);

            NekDouble predResNorm = ComputeReducedResidualNorm(data, yHook);

            Vmath::Vcopy(ntotal, m_Solution, 1, m_solutionTrial, 1);
            Vmath::Svtvp(ntotal, -1.0, m_deltaTrial, 1, m_solutionTrial, 1,
                         m_solutionTrial, 1);

            m_operator.DoNekSysResEval(m_solutionTrial, m_residualTrial);

            NekDouble newResNormSq =
                Vmath::Dot(ntotal, m_residualTrial, m_residualTrial);
            m_rowComm->AllReduce(newResNormSq, LibUtilities::ReduceSum);
            NekDouble newResNorm = sqrt(newResNormSq);

            NekDouble ared = oldResNorm - newResNorm;
            NekDouble pred = oldResNorm - predResNorm;
            NekDouble rho  = (pred > 1.0e-14 * oldResNorm) ? ared / pred : -1.0;

            NekDouble yNorm   = CalcL2Norm(data.nswp, yHook);
            bool boundaryStep = std::abs(yNorm - enlargedRadius) <=
                                1.0e-2 * (enlargedRadius + 1.0e-14);

            cout << "  HookStep exploration trial = " << enlargedRadius
                 << " yNorm = " << yNorm << " oldRes = " << oldResNorm
                 << " newRes = " << newResNorm << " predRes = " << predResNorm
                 << " ared = " << ared << " pred = " << pred << " rho = " << rho
                 << endl;

            bool meaningfulDecrease =
                HasMeaningfulDecrease(oldResNorm, newResNorm);

            // Candidate must still be admissible.
            if (!(rho > m_rhoAccept || meaningfulDecrease))
            {
                cout << " Exploration step rejected, stop growing radius."
                     << endl;
                break;
            }

            // Keep only if actual residual improves on the best accepted one.
            if (newResNorm < bestResNorm)
            {
                cout << " Exploration improved residual from " << bestResNorm
                     << " to " << newResNorm << ", keeping enlarged step."
                     << endl;

                Vmath::Vcopy(ntotal, m_solutionTrial, 1, m_solutionBest, 1);
                Vmath::Vcopy(ntotal, m_residualTrial, 1, m_residualBest, 1);

                bestResNorm      = newResNorm;
                bestRho          = rho;
                bestTrustRadius  = enlargedRadius;
                bestBoundaryStep = boundaryStep;

                exploreRadius = enlargedRadius;

                if (!(bestRho > m_rhoGrow && bestBoundaryStep))
                {
                    cout << " Exploration stopping because enlarged step is "
                         << "no longer a strong boundary-limited candidate."
                         << endl;
                    break;
                }
            }
            else
            {
                cout << " Exploration did not improve residual (best = "
                     << bestResNorm << ", trial = " << newResNorm
                     << "), stopping." << endl;
                break;
            }
        }
    }

    // ---------------------------------------------------------------------
    // Finalize with the best accepted candidate found.
    // ---------------------------------------------------------------------
    Vmath::Vcopy(ntotal, m_solutionBest, 1, m_Solution, 1);
    Vmath::Vcopy(ntotal, m_residualBest, 1, m_Residual, 1);
    m_haveUpdatedResidual = true;

    m_trustRadius = bestTrustRadius;

    if (bestRho < m_rhoShrink)
    {
        cout << " Trust Radius shrink because rho < m_rhoShrink = "
             << m_rhoShrink << " m_shrinkFactor = " << m_shrinkFactor << endl;
        m_trustRadius =
            std::max(m_trustRadiusMin, m_shrinkFactor * m_trustRadius);
        cout << " Next Trust Radius = " << m_trustRadius << endl;
    }
    else if (bestRho > m_rhoGrow && bestBoundaryStep)
    {
        cout << " Trust Radius grow because rho > m_rhoGrow = " << m_rhoGrow
             << " m_growFactor = " << m_growFactor << endl;
        m_trustRadius =
            std::min(m_trustRadiusMax, m_growFactor * m_trustRadius);
        cout << " Next Trust Radius = " << m_trustRadius << endl;
    }

    return true;
}

bool NekNonlinSysIterNewtonHookStep::HasMeaningfulDecrease(
    const NekDouble oldResNorm, const NekDouble newResNorm) const
{
    NekDouble ared = oldResNorm - newResNorm;
    NekDouble thresh =
        std::max(m_meaningfulAredAbs, m_meaningfulAredRel * oldResNorm);
    return ared > thresh;
}
} // namespace Nektar::LibUtilities
