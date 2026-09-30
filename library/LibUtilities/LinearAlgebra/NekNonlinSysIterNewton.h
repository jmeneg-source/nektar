///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewton.h
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
// Description: NekNonlinSysIterNewton header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_H

#include <LibUtilities/LinearAlgebra/NekNonlinSysIter.h>

namespace Nektar::LibUtilities
{

class NekNonlinSysIterNewton : public NekNonlinSysIter
{
public:
    /// Constructor for full direct matrix solve.
    friend class MemoryManager<NekNonlinSysIterNewton>;

    LIB_UTILITIES_EXPORT static NekNonlinSysIterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekNonlinSysIterSharedPtr p =
            MemoryManager<NekNonlinSysIterNewton>::AllocateSharedPtr(
                pSession, vRowComm, nDimen, pKey);
        p->InitObject();
        return p;
    }

    static std::string className;

    LIB_UTILITIES_EXPORT NekNonlinSysIterNewton(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
        const NekSysKey &pKey);
    LIB_UTILITIES_EXPORT ~NekNonlinSysIterNewton() override = default;

protected:
    NekDouble m_NewtonScale;

    bool m_inexactNewtonForcing = false;
    // Classic bounded Eisenstat-Walker style forcing
    NekDouble m_forcingEtaInit = 1.0e-2;
    NekDouble m_forcingEtaMin  = 1.0e-6;
    NekDouble m_forcingEtaMax  = 5.0e-2;
    NekDouble m_forcingGamma   = 0.9;
    NekDouble m_forcingAlpha   = 1.5;

    bool m_haveUpdatedResidual = false;

    void v_InitObject() override;

    int v_SolveSystem(const int nGlobal,
                      const Array<OneD, const NekDouble> &pInput,
                      Array<OneD, NekDouble> &pOutput, const int nDir) override;

    virtual bool v_ApplyNewtonUpdate(const int ntotal,
                                     const NekDouble oldResNorm);

private:
    NekDouble CalcInexactNewtonForcing(const int &nIteration,
                                       const NekDouble &resnormOld,
                                       const NekDouble &resnorm);
};

} // namespace Nektar::LibUtilities

#endif
