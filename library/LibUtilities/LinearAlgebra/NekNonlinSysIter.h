///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIter.h
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
// Description: NekNonlinSysIter header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H

#include <LibUtilities/LinearAlgebra/NekLinSysIter.h>
#include <LibUtilities/LinearAlgebra/NekSys.h>

namespace Nektar::LibUtilities
{
class NekNonlinSysIter;

typedef std::shared_ptr<NekNonlinSysIter> NekNonlinSysIterSharedPtr;

typedef LibUtilities::NekFactory<
    std::string, NekNonlinSysIter, const LibUtilities::SessionReaderSharedPtr &,
    const LibUtilities::CommSharedPtr &, const int, const NekSysKey &>
    NekNonlinSysIterFactory;
LIB_UTILITIES_EXPORT NekNonlinSysIterFactory &GetNekNonlinSysIterFactory();

class NekNonlinSysIter : public NekSys
{
public:
    LIB_UTILITIES_EXPORT NekNonlinSysIter(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
        const NekSysKey &pKey);
    LIB_UTILITIES_EXPORT ~NekNonlinSysIter() override = default;

    LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> &GetRefSolution()
        const
    {
        return m_Solution;
    }

    LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> &GetRefResidual()
        const
    {
        return m_Residual;
    }

    LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> &GetRefSourceVec()
        const
    {
        return m_SourceVec;
    }

    LIB_UTILITIES_EXPORT const NekLinSysIterSharedPtr &GetLinSys() const
    {
        return m_linsol;
    }

    LIB_UTILITIES_EXPORT int GetNtotLinSysIts() const
    {
        return m_NtotLinSysIts;
    }

    LIB_UTILITIES_EXPORT void SetNekLinSysTolerance(const NekDouble in)
    {
        m_NekLinSysTolerance = in;
    }

    LIB_UTILITIES_EXPORT void SetNekNonlinSysTolerance(const NekDouble in)
    {
        m_NekNonLinSysTolerance = in;
    }

    LIB_UTILITIES_EXPORT void SetNonlinIterTolRelativeL2(const NekDouble in)
    {
        m_NonlinIterTolRelativeL2 = in;
    }

    LIB_UTILITIES_EXPORT void SetNekNonlinSysMaxIterations(const int in)
    {
        m_NekNonlinSysMaxIterations = in;
    }

protected:
    NekDouble m_SysResNorm0;
    NekDouble m_SysResNorm;

    NekLinSysIterSharedPtr m_linsol;

    NekDouble m_NekLinSysTolerance;
    NekDouble m_NekNonLinSysTolerance;
    NekDouble m_NonlinIterTolRelativeL2;

    int m_NekNonlinSysMaxIterations;
    int m_totalIterations = 0;
    int m_NtotLinSysIts   = 0;

    std::string m_LinSysIterSolverType;

    Array<OneD, NekDouble> m_Solution;
    Array<OneD, NekDouble> m_Residual;
    Array<OneD, NekDouble> m_DeltSltn;
    Array<OneD, NekDouble> m_SourceVec;

    void v_InitObject() override;

    void v_SetSysOperators(const NekSysOperators &in) override;

    void ConvergenceCheck(const int nIteration,
                          const Array<OneD, const NekDouble> &Residual);

private:
};
} // namespace Nektar::LibUtilities
#endif
