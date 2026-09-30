///////////////////////////////////////////////////////////////////////////////
//
// File: NekLinSysIterEvsDirect.cpp
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
// Description: NekLinSysIterEvsDirect definition
//
///////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/LinearAlgebra/NekLinSysIterEvsDirect.h>

namespace Nektar::LibUtilities
{
/**
 * @class  NekLinSysIterEvsDirect
 *
 * Solves a linear system using iterative methods.
 */
std::string NekLinSysIterEvsDirect::className =
    LibUtilities::GetNekLinSysIterFactory().RegisterCreatorFunction(
        "EvsDirect", NekLinSysIterEvsDirect::create,
        "NekLinSysIterEvsDirect solver.");

NekLinSysIterEvsDirect::NekLinSysIterEvsDirect(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm, const int nDimen,
    const NekSysKey &pKey)
    : NekLinSysIter(pSession, vComm, nDimen, pKey)
{
}

void NekLinSysIterEvsDirect::v_InitObject()
{
    NekLinSysIter::v_InitObject();
}

NekLinSysIterEvsDirect::~NekLinSysIterEvsDirect()
{
}

/**
 *
 */
int NekLinSysIterEvsDirect::v_SolveSystem(
    const int nGlobal,
    [[maybe_unused]] const Array<OneD, const NekDouble> &pInput,
    [[maybe_unused]] Array<OneD, NekDouble> &pOutput, const int nDir)
{
    int nNonDir = nGlobal - nDir;

    Array<OneD, NekDouble> Mat(nNonDir * nNonDir);
    Array<OneD, NekDouble> wk1(nGlobal), wk2(nGlobal), tmp;

    for (int i = 0; i < nNonDir; ++i)
    {
        // need to zero global array since in LhsEval have a global to
        // local call
        Vmath::Zero(nGlobal, wk1, 1);
        wk1[nDir + i] = 1.0;

        m_operator.DoNekSysLhsEval(wk1, wk2);
        m_operator.DoNekSysPrecon(wk2 + nDir, tmp = wk1 + nDir);

        Vmath::Vcopy(nNonDir, wk1 + nDir, 1, tmp = Mat + i, nNonDir);
    }

    ////////////////////////////////////////////////////////
    // Print Matrix
    std::ofstream mFile("Matrix.txt");

    for (int j = 0; j < nNonDir; j++)
    {
        for (int k = 0; k < nNonDir; k++)
        {
            if (fabs(Mat[j * nNonDir + k]) < 1e-15)
            {
                Mat[j * nNonDir + k] = 0.0;
            }
            mFile << Mat[j * nNonDir + k];
        }
        mFile << std::endl;
    }
    mFile.close();

    char jobvl = 'N', jobvr = 'N';
    int info = 0, lwork = 3 * nNonDir;
    NekDouble dum;

    Array<OneD, NekDouble> EIG_R(nNonDir), EIG_I(nNonDir), work(lwork);

    Lapack::Dgeev(jobvl, jobvr, nNonDir, Mat.data(), nNonDir, EIG_R.data(),
                  EIG_I.data(), &dum, 1, &dum, 1, &work[0], lwork, info);

    ////////////////////////////////////////////////////////
    // Output of the EigenValues
    std::ofstream pFile("Eigenvalues.txt");

    for (int j = 0; j < nNonDir; j++)
    {
        pFile << EIG_R[j] << " " << EIG_I[j] << std::endl;
    }
    pFile.close();

    std::cout << std::endl << "Eigenvalues : " << std::endl;
    for (int j = 0; j < nNonDir; j++)
    {
        std::cout << EIG_R[j] << "\t" << EIG_I[j] << std::endl;
    }
    std::cout << std::endl;

    return 0;
}
} // namespace Nektar::LibUtilities
