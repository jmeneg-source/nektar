///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewtonTrustRegion.h
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
// Description: NekNonlinSysIterNewtonTrustRegion header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_TRUST_REGION_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_TRUST_REGION_H

#include <LibUtilities/LinearAlgebra/NekNonlinSysIterNewton.h>

namespace Nektar::LibUtilities
{

class NekNonlinSysIterNewtonTrustRegion : public NekNonlinSysIterNewton
{
public:
    friend class MemoryManager<NekNonlinSysIterNewtonTrustRegion>;

    LIB_UTILITIES_EXPORT static NekNonlinSysIterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekNonlinSysIterSharedPtr p =
            MemoryManager<NekNonlinSysIterNewtonTrustRegion>::AllocateSharedPtr(
                pSession, vComm, nDimen, pKey);
        p->InitObject();
        return p;
    }

    static std::string className;

    LIB_UTILITIES_EXPORT NekNonlinSysIterNewtonTrustRegion(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey);

    LIB_UTILITIES_EXPORT ~NekNonlinSysIterNewtonTrustRegion() override =
        default;

protected:
    NekDouble m_trustRadiusInit = 1.0e-2;
    NekDouble m_trustRadius     = 1.0e-2;
    NekDouble m_trustRadiusMin  = 1.0e-8;
    NekDouble m_trustRadiusMax  = 1.0e2;

    NekDouble m_rhoAccept    = 1.0e-1;
    NekDouble m_rhoShrink    = 2.5e-1;
    NekDouble m_rhoGrow      = 7.5e-1;
    NekDouble m_shrinkFactor = 0.5;
    NekDouble m_growFactor   = 2.0;

    int m_maxTrustRegionSteps             = 4;
    NekDouble m_maxFallbackResidualGrowth = 10.0;

    Array<OneD, NekDouble> m_solutionTrial;
    Array<OneD, NekDouble> m_residualTrial;
    Array<OneD, NekDouble> m_deltaTrial;

    void v_InitObject() override;

    bool v_ApplyNewtonUpdate(const int ntotal,
                             const NekDouble oldResNorm) override;

private:
    NekDouble CalcL2Norm(const int ntotal,
                         const Array<OneD, const NekDouble> &x) const;
};

} // namespace Nektar::LibUtilities

#endif
