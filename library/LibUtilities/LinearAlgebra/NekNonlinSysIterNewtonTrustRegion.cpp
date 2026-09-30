///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewtonTrustRegion.cpp
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
// Description: NekNonlinSysIterNewtonTrustRegion definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSysIterNewtonTrustRegion.h>
#include <algorithm>
#include <cmath>

/* -- XML Parameter section --
<SOLVERINFO>
    <I PROPERTY="NonlinearSolver" VALUE="NewtonTrustRegion" />
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
    <P> TrustRegionMaxFallbackResidualGrowth  = 10.0   </P>
</PARAMETERS>
*/

using namespace std;
namespace Nektar::LibUtilities
{

string NekNonlinSysIterNewtonTrustRegion::className =
    LibUtilities::GetNekNonlinSysIterFactory().RegisterCreatorFunction(
        "NewtonTrustRegion", NekNonlinSysIterNewtonTrustRegion::create,
        "Newton solver with trust-region clipped step.");

NekNonlinSysIterNewtonTrustRegion::NekNonlinSysIterNewtonTrustRegion(
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
    pSession->LoadParameter("TrustRegionMaxFallbackResidualGrowth",
                            m_maxFallbackResidualGrowth, 10.0);
}

void NekNonlinSysIterNewtonTrustRegion::v_InitObject()
{
    NekNonlinSysIterNewton::v_InitObject();

    m_solutionTrial = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_residualTrial = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_deltaTrial    = Array<OneD, NekDouble>(m_SysDimen, 0.0);

    m_trustRadius = m_trustRadiusInit;
}

NekDouble NekNonlinSysIterNewtonTrustRegion::CalcL2Norm(
    const int ntotal, const Array<OneD, const NekDouble> &x) const
{
    NekDouble nrm2 = Vmath::Dot(ntotal, x, x);
    m_rowComm->AllReduce(nrm2, LibUtilities::ReduceSum);
    return sqrt(nrm2);
}

bool NekNonlinSysIterNewtonTrustRegion::v_ApplyNewtonUpdate(
    const int ntotal, const NekDouble oldResNorm)
{
    NekDouble deltaNorm        = CalcL2Norm(ntotal, m_DeltSltn);
    NekDouble trustRadiusTrial = m_trustRadius;
    if (m_verbose)
    {
        cout << " Attempting to apply Newton update Trust Region "
             << " Trust radius = " << trustRadiusTrial
             << " deltaNorm = " << deltaNorm << " oldResNorm= " << oldResNorm
             << endl;
    }
    for (int it = 0; it < m_maxTrustRegionSteps; ++it)
    {
        if (m_verbose)
        {
            cout << " Trust region test " << it << " of "
                 << m_maxTrustRegionSteps << endl;
        }
        Vmath::Vcopy(ntotal, m_DeltSltn, 1, m_deltaTrial, 1);

        bool boundaryStep = false;
        if (deltaNorm > trustRadiusTrial)
        {
            NekDouble scale = trustRadiusTrial / deltaNorm;
            Vmath::Smul(ntotal, scale, m_deltaTrial, 1, m_deltaTrial, 1);
            boundaryStep = true;
            if (m_verbose)
            {
                cout << " deltaNorm > trustRadiusTrial - Boundary Step" << endl;
            }
        }
        NekDouble trialNorm = CalcL2Norm(ntotal, m_deltaTrial);
        NekDouble alphaPred =
            (deltaNorm > 1.0e-14) ? (trialNorm / deltaNorm) : 1.0;

        Vmath::Vcopy(ntotal, m_Solution, 1, m_solutionTrial, 1);
        Vmath::Svtvp(ntotal, -1.0, m_deltaTrial, 1, m_solutionTrial, 1,
                     m_solutionTrial, 1);

        m_operator.DoNekSysResEval(m_solutionTrial, m_residualTrial);

        NekDouble newResNormSq =
            Vmath::Dot(ntotal, m_residualTrial, m_residualTrial);
        m_rowComm->AllReduce(newResNormSq, LibUtilities::ReduceSum);
        NekDouble newResNorm = sqrt(newResNormSq);

        NekDouble pred =
            std::max(1.0e-14 * oldResNorm, 0.5 * alphaPred * oldResNorm);
        NekDouble ared = oldResNorm - newResNorm;
        NekDouble rho  = ared / pred;

        if (m_verbose)
        {
            cout << "  TrustRadius trial = " << trustRadiusTrial
                 << " stepNorm = " << trialNorm << " oldRes = " << oldResNorm
                 << " newRes = " << newResNorm << " ared = " << ared
                 << " pred = " << pred << " rho = " << rho << endl;
        }

        if (rho > m_rhoAccept)
        {
            cout << " Accepting step with rho = " << rho
                 << "> m_rhoAccept = " << m_rhoAccept << endl;
            Vmath::Vcopy(ntotal, m_solutionTrial, 1, m_Solution, 1);
            Vmath::Vcopy(ntotal, m_residualTrial, 1, m_Residual, 1);
            m_haveUpdatedResidual = true;

            m_trustRadius = trustRadiusTrial;

            if (rho < m_rhoShrink)
            {
                cout << "Trust radius shrink with m_rhoShrink=" << m_rhoShrink
                     << endl;
                m_trustRadius =
                    std::max(m_trustRadiusMin, m_shrinkFactor * m_trustRadius);
            }
            else if (rho > m_rhoGrow && boundaryStep)
            {
                cout << "Trust radius grow with m_rhoGrow=" << m_rhoGrow
                     << endl;
                m_trustRadius =
                    std::min(m_trustRadiusMax, m_growFactor * m_trustRadius);
            }

            return true;
        }

        trustRadiusTrial =
            std::max(m_trustRadiusMin, m_shrinkFactor * trustRadiusTrial);

        if (trustRadiusTrial <= m_trustRadiusMin + 1.0e-16)
        {
            break;
        }
    }

    cout << "Failed to find a correct step "
         << " Fallback: force a minimum-radius step and continue." << endl;
    Vmath::Vcopy(ntotal, m_DeltSltn, 1, m_deltaTrial, 1);

    if (deltaNorm > m_trustRadiusMin)
    {
        NekDouble scale = m_trustRadiusMin / deltaNorm;
        Vmath::Smul(ntotal, scale, m_deltaTrial, 1, m_deltaTrial, 1);
    }

    Vmath::Vcopy(ntotal, m_Solution, 1, m_solutionTrial, 1);
    Vmath::Svtvp(ntotal, -1.0, m_deltaTrial, 1, m_solutionTrial, 1,
                 m_solutionTrial, 1);

    m_operator.DoNekSysResEval(m_solutionTrial, m_residualTrial);

    NekDouble fallbackResNormSq =
        Vmath::Dot(ntotal, m_residualTrial, m_residualTrial);
    m_rowComm->AllReduce(fallbackResNormSq, LibUtilities::ReduceSum);
    NekDouble fallbackResNorm = sqrt(fallbackResNormSq);

    cout << " No acceptable trust-region step found."
         << " Taking fallback TrustRadiusMin = " << m_trustRadiusMin
         << " residual = " << fallbackResNorm << endl;

    if (!std::isfinite(fallbackResNorm) ||
        fallbackResNorm > m_maxFallbackResidualGrowth * oldResNorm)
    {
        WARNINGL0(false,
                  "Fallback trust-region step caused non-finite or excessive "
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

} // namespace Nektar::LibUtilities
