///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewton.cpp
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
// Description: NekNonlinSysIterNewton definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSysIterNewton.h>

using namespace std;

namespace Nektar::LibUtilities
{
/**
 * @class  NekNonlinSysIterNewton
 *
 * Solves a nonlinear system using iterative methods.
 */
string NekNonlinSysIterNewton::className =
    LibUtilities::GetNekNonlinSysIterFactory().RegisterCreatorFunction(
        "Newton", NekNonlinSysIterNewton::create,
        "NekNonlinSysIterNewton solver.");

NekNonlinSysIterNewton::NekNonlinSysIterNewton(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
    const NekSysKey &pKey)
    : NekNonlinSysIter(pSession, vRowComm, nDimen, pKey)
{
    int inexactNewtonForcing = 0;
    pSession->LoadParameter("InexactNewtonForcing", inexactNewtonForcing, 0);
    m_inexactNewtonForcing = (bool)inexactNewtonForcing;
    pSession->LoadParameter("ForcingEtaInit", m_forcingEtaInit, 1.0e-2);
    pSession->LoadParameter("ForcingEtaMin", m_forcingEtaMin,
                            m_NekLinSysTolerance);
    pSession->LoadParameter("ForcingEtaMax", m_forcingEtaMax, 5.0e-2);
    pSession->LoadParameter("ForcingGamma", m_forcingGamma, 0.9);
    pSession->LoadParameter("ForcingAlpha", m_forcingAlpha, 1.5);
    pSession->LoadParameter("NewtonScale", m_NewtonScale, 1.0);
}

void NekNonlinSysIterNewton::v_InitObject()
{
    NekNonlinSysIter::v_InitObject();
}

/**
 *
 **/
int NekNonlinSysIterNewton::v_SolveSystem(
    const int nGlobal, const Array<OneD, const NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, [[maybe_unused]] const int nDir)
{
    ASSERTL0(nGlobal == m_SysDimen, "nGlobal != m_SysDimen");

    m_SourceVec     = pInput;
    m_Solution      = pOutput;
    m_NtotLinSysIts = 0;
    m_converged     = false;

    Vmath::Vcopy(nGlobal, pInput, 1, m_Solution, 1);

    NekDouble resnormOld = 0.0;
    int NttlNonlinIte    = 0;
    for (; NttlNonlinIte < m_NekNonlinSysMaxIterations; ++NttlNonlinIte)
    {
        if (!m_haveUpdatedResidual)
        {
            m_operator.DoNekSysResEval(m_Solution, m_Residual, true);
        }
        else
        {
            m_haveUpdatedResidual = false;
        }

        ConvergenceCheck(NttlNonlinIte, m_Residual);
        if (m_converged)
        {
            break;
        }

        NekDouble LinSysRelativeIteTol =
            CalcInexactNewtonForcing(NttlNonlinIte, resnormOld, m_SysResNorm);
        resnormOld = m_SysResNorm;
        m_linsol->SetRhsMagnitude(m_SysResNorm);
        m_linsol->SetNekLinSysTolerance(LinSysRelativeIteTol);
        if (m_verbose)
        {
            cout << "Newton Non-It=" << NttlNonlinIte
                 << " m_SysResNorm0=" << m_SysResNorm0
                 << " m_SysResNorm=" << m_SysResNorm << " "
                 << "(RES=" << sqrt(m_SysResNorm)
                 << " Res/Res0= " << sqrt(m_SysResNorm / m_SysResNorm0)
                 << " Res/DtRHS= " << sqrt(m_SysResNorm / m_rhs_magnitude)
                 << " LinSysRelativeIteTol=" << LinSysRelativeIteTol << endl;
            cout << "Will apply Newton correction with Newton Scale = "
                 << m_NewtonScale << endl;
        }
        int ntmpLinSysIts =
            m_linsol->SolveSystem(nGlobal, m_Residual, m_DeltSltn, 0);
        m_NtotLinSysIts += ntmpLinSysIts;

        NekDouble oldResNorm = sqrt(m_SysResNorm);
        bool accepted        = v_ApplyNewtonUpdate(nGlobal, oldResNorm);
        if (!accepted)
        {
            WARNINGL0(false, "Newton step rejected.");
            break;
        }
    }

    if ((!m_converged || m_verbose) && m_root && m_FlagWarnings)
    {
        int nwidthcolm = 11;

        WARNINGL0(m_converged,
                  "     # Nonlinear solver not converge in DoImplicitSolve");
        cout << right << scientific << setw(nwidthcolm)
             << setprecision(nwidthcolm - 6)
             << "     * Newton-Its converged (RES=" << sqrt(m_SysResNorm)
             << " Res/Res0= " << sqrt(m_SysResNorm / m_SysResNorm0)
             << " Res/DtRHS= " << sqrt(m_SysResNorm / m_rhs_magnitude)
             << " with " << setw(3) << NttlNonlinIte << " Non-Its)" << endl;
    }

    return NttlNonlinIte;
}

bool NekNonlinSysIterNewton::v_ApplyNewtonUpdate(
    const int ntotal, [[maybe_unused]] const NekDouble oldResNorm)
{
    Vmath::Svtvp(ntotal, -1.0 * m_NewtonScale, m_DeltSltn, 1, m_Solution, 1,
                 m_Solution, 1);
    return true;
}

NekDouble NekNonlinSysIterNewton::CalcInexactNewtonForcing(
    const int &nIteration, const NekDouble &resnormOld,
    const NekDouble &resnorm)
{
    if (!m_inexactNewtonForcing)
    {
        return m_NekLinSysTolerance;
    }

    if (nIteration == 0 || resnormOld <= 0.0)
    {
        return m_forcingEtaInit;
    }

    // resnorm and resnormOld are stored as squared norms in this class
    NekDouble rk   = sqrt(resnorm);
    NekDouble rkm1 = sqrt(resnormOld);

    NekDouble eta = m_forcingGamma * pow(rk / rkm1, m_forcingAlpha);

    eta = std::max(m_forcingEtaMin, std::min(m_forcingEtaMax, eta));

    return eta;
}

} // namespace Nektar::LibUtilities
