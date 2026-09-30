///////////////////////////////////////////////////////////////////////////////
//
// File: NekSys.cpp
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
// Description: NekSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekSys.h>

namespace Nektar::LibUtilities
{
/**
 * @class  NekSys
 *
 * Solves a nonlinear or linear system.
 */

NekSys::NekSys(const LibUtilities::SessionReaderSharedPtr &pSession,
               const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
               [[maybe_unused]] const NekSysKey &pKey)
{
    m_rowComm      = vRowComm;
    m_root         = m_rowComm->GetRank() == 0;
    m_FlagWarnings = true;
    m_verbose      = pSession->DefinesCmdLineArgument("verbose");
    m_converged    = false;
    m_SysDimen     = nDimen;
}

} // namespace Nektar::LibUtilities
