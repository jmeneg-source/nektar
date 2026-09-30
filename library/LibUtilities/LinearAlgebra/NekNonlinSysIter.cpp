///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIter.cpp
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
// Description: NekNonlinSysIter definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSysIter.h>

namespace Nektar::LibUtilities
{
/**
 * @class  NekNonlinSysIter
 *
 * Solves a nonlinear system using iterative methods.
 */
NekNonlinSysIterFactory &GetNekNonlinSysIterFactory()
{
    static NekNonlinSysIterFactory instance;
    return instance;
}

NekNonlinSysIter::NekNonlinSysIter(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
    const NekSysKey &pKey)
    : NekSys(pSession, vRowComm, nDimen, pKey)
{
    m_NekLinSysTolerance        = pKey.m_NekLinSysTolerance;
    m_NekNonLinSysTolerance     = pKey.m_NekNonLinSysTolerance;
    m_NonlinIterTolRelativeL2   = pKey.m_NonlinIterTolRelativeL2;
    m_NekNonlinSysMaxIterations = pKey.m_NekNonlinSysMaxIterations;
    m_LinSysIterSolverType      = pKey.m_LinSysIterSolverTypeInNonlin;

    ASSERTL0(LibUtilities::GetNekLinSysIterFactory().ModuleExists(
                 m_LinSysIterSolverType),
             "NekLinSysIter '" + m_LinSysIterSolverType +
                 "' is not defined.\n");

    m_linsol = LibUtilities::GetNekLinSysIterFactory().CreateInstance(
        m_LinSysIterSolverType, pSession, m_rowComm, m_SysDimen, pKey);
    m_linsol->SetFlagWarnings(false);
}

void NekNonlinSysIter::v_InitObject()
{
    NekSys::v_InitObject();
    m_Residual = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_DeltSltn = Array<OneD, NekDouble>(m_SysDimen, 0.0);
}

void NekNonlinSysIter::v_SetSysOperators(const NekSysOperators &in)
{
    NekSys::v_SetSysOperators(in);
    m_linsol->SetSysOperators(in);
}

void NekNonlinSysIter::ConvergenceCheck(
    const int nIteration, const Array<OneD, const NekDouble> &Residual)
{
    m_SysResNorm = Vmath::Dot(Residual.size(), Residual, Residual);
    m_rowComm->AllReduce(m_SysResNorm, Nektar::LibUtilities::ReduceSum);

    if (nIteration == 0)
    {
        m_SysResNorm0 = m_SysResNorm;
    }

    NekDouble resratio = m_SysResNorm / m_SysResNorm0;
    NekDouble tol      = m_NekNonLinSysTolerance;
    NekDouble restol   = m_NonlinIterTolRelativeL2;

    if (m_rhs_magnitude == NekConstants::kNekUnsetDouble)
    {
        m_rhs_magnitude = 1.0;
    }

    m_converged = resratio < restol * restol ||
                  m_SysResNorm < tol * tol * m_rhs_magnitude;
}

} // namespace Nektar::LibUtilities
