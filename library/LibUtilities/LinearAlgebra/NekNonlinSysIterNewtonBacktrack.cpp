///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewtonBacktrack.cpp
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
// Description: NekNonlinSysIterNewtonBacktrack definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSysIterNewtonBacktrack.h>
#include <algorithm>

/* -- XML Parameter section --
<SOLVERINFO>
    <I PROPERTY="NonlinearSolver" VALUE="NewtonBacktrack" />
</SOLVERINFO>

<PARAMETERS>
    <P> NewtonScaleInit        = 1.0e-2 </P>
    <P> NewtonScaleMax         = 1.0e-1 </P>
    <P> NewtonScaleMin         = 1.0e-6 </P>
    <P> NewtonScaleGrowFactor  = 2.0    </P>
    <P> BacktrackFactor        = 0.5    </P>
    <P> MaxBacktrackSteps      = 4      </P>
    <P> AcceptReductionEta     = 1.0e-4 </P>
</PARAMETERS>
*/

using namespace std;

namespace Nektar::LibUtilities
{

string NekNonlinSysIterNewtonBacktrack::className =
    LibUtilities::GetNekNonlinSysIterFactory().RegisterCreatorFunction(
        "NewtonBacktrack", NekNonlinSysIterNewtonBacktrack::create,
        "Newton solver with residual-based backtracking.");

NekNonlinSysIterNewtonBacktrack::NekNonlinSysIterNewtonBacktrack(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm, const int nDimen,
    const NekSysKey &pKey)
    : NekNonlinSysIterNewton(pSession, vComm, nDimen, pKey)
{
    pSession->LoadParameter("NewtonScaleInit", m_newtonScaleInit, 1.0e-2);
    pSession->LoadParameter("NewtonScaleMax", m_newtonScaleMax, 1.0e-1);
    pSession->LoadParameter("NewtonScaleMin", m_newtonScaleMin, 1.0e-6);
    pSession->LoadParameter("NewtonScaleGrowFactor", m_newtonScaleGrowFactor,
                            2.0);
    pSession->LoadParameter("BacktrackFactor", m_backtrackFactor, 0.5);
    pSession->LoadParameter("MaxBacktrackSteps", m_maxBacktrackSteps, 4);
    pSession->LoadParameter("AcceptReductionEta", m_acceptReductionEta, 1.0e-4);
}

void NekNonlinSysIterNewtonBacktrack::v_InitObject()
{
    NekNonlinSysIterNewton::v_InitObject();

    m_solutionTrial = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_residualTrial = Array<OneD, NekDouble>(m_SysDimen, 0.0);

    // In the backtracking solver, m_NewtonScale is the adaptive starting
    // trial scale for the current Newton step.
    m_NewtonScale = m_newtonScaleInit;
}

bool NekNonlinSysIterNewtonBacktrack::v_ApplyNewtonUpdate(
    const int ntotal, const NekDouble oldResNorm)
{
    NekDouble alphaStart = std::min(m_NewtonScale, m_newtonScaleMax);
    alphaStart           = std::max(alphaStart, m_newtonScaleMin);

    if (m_verbose)
    {
        cout << " Attempting to apply Newton update with NewtonBacktrack"
             << endl;
    }
    NekDouble alpha = alphaStart;

    for (int bt = 0; bt < m_maxBacktrackSteps; ++bt)
    {
        Vmath::Vcopy(ntotal, m_Solution, 1, m_solutionTrial, 1);
        Vmath::Svtvp(ntotal, -alpha, m_DeltSltn, 1, m_solutionTrial, 1,
                     m_solutionTrial, 1);

        m_operator.DoNekSysResEval(m_solutionTrial, m_residualTrial);

        NekDouble newResNormSq =
            Vmath::Dot(ntotal, m_residualTrial, m_residualTrial);
        m_rowComm->AllReduce(newResNormSq, LibUtilities::ReduceSum);
        NekDouble newResNorm = sqrt(newResNormSq);

        if (m_verbose)
        {
            cout << "Bactracking trial " << bt + 1 << " of "
                 << m_maxBacktrackSteps + 1 << " with scale = " << alpha
                 << "\n new residual = " << newResNorm
                 << "\n old residual = " << oldResNorm
                 << "\n newResNorm / oldResNorm=" << newResNorm / oldResNorm
                 << endl;
        }

        if (newResNorm <= (1.0 - m_acceptReductionEta * alpha) * oldResNorm)
        {
            cout << "newResNorm <= (1.0 - m_acceptReductionEta * alpha) * "
                    "oldResNorm"
                 << endl;
            cout << "newResNorm = " << newResNorm
                 << " m_acceptReductionEta = " << m_acceptReductionEta
                 << " alpha = " << alpha << endl;
            Vmath::Vcopy(ntotal, m_solutionTrial, 1, m_Solution, 1);
            Vmath::Vcopy(ntotal, m_residualTrial, 1, m_Residual, 1);
            m_haveUpdatedResidual = true;

            if (bt == 0)
            {
                cout << "Good on the first try, grow alpha" << endl;
                m_NewtonScale = std::min(m_newtonScaleMax,
                                         m_newtonScaleGrowFactor * alphaStart);
            }
            else
            {
                cout << "Good step, conserve alpha" << endl;
                m_NewtonScale = alpha;
            }

            cout << "  Accepted Newton scale = " << alpha
                 << " next starting scale = " << m_NewtonScale << endl;
            return true;
        }

        alpha *= m_backtrackFactor;
        if (alpha < m_newtonScaleMin)
        {
            break;
        }
    }

    cout << "Fallback: take minimum allowed step and continue." << endl;
    NekDouble alphaFallback = m_newtonScaleMin;

    Vmath::Vcopy(ntotal, m_Solution, 1, m_solutionTrial, 1);
    Vmath::Svtvp(ntotal, -alphaFallback, m_DeltSltn, 1, m_solutionTrial, 1,
                 m_solutionTrial, 1);

    m_operator.DoNekSysResEval(m_solutionTrial, m_residualTrial);

    NekDouble fallbackResNormSq =
        Vmath::Dot(ntotal, m_residualTrial, m_residualTrial);
    m_rowComm->AllReduce(fallbackResNormSq, LibUtilities::ReduceSum);
    NekDouble fallbackResNorm = sqrt(fallbackResNormSq);

    cout << "  No acceptable backtracking step found."
         << " Taking fallback Newton scale = " << alphaFallback
         << "\n new residual = " << fallbackResNorm
         << "\n old residual = " << oldResNorm << endl;
    cout << "newResNorm / oldResNorm=" << fallbackResNorm / oldResNorm << endl;

    Vmath::Vcopy(ntotal, m_solutionTrial, 1, m_Solution, 1);
    Vmath::Vcopy(ntotal, m_residualTrial, 1, m_Residual, 1);
    m_haveUpdatedResidual = true;

    // Keep next iteration conservative.
    m_NewtonScale = m_newtonScaleMin;

    return true;
}

} // namespace Nektar::LibUtilities