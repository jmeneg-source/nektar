///////////////////////////////////////////////////////////////////////////////
//
// File: NonlinSysIterDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Demo for testing functionality of NekNonlinSys classes
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekNonlinSysIter.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>
using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;

class NonlinSysDemo
{
    typedef const Array<OneD, const NekDouble> InArrayType;
    typedef Array<OneD, NekDouble> OutArrayType;

public:
    NonlinSysDemo(const LibUtilities::SessionReaderSharedPtr &pSession,
                  const LibUtilities::CommSharedPtr &pComm)
        : m_session(pSession), m_comm(pComm)
    {
        AllocateInitMatrix();

        std::string SolverType = "Newton";
        if (pSession->DefinesSolverInfo("NonlinearSolver"))
        {
            SolverType = pSession->GetSolverInfo("NonlinearSolver");
        }
        ASSERTL0(
            LibUtilities::GetNekNonlinSysIterFactory().ModuleExists(SolverType),
            "NekNonlinSys '" + SolverType + "' is not defined.\n");

        // Create the key to hold solver settings
        auto sysKey = LibUtilities::NekSysKey();
        // Load required LinSys parameters:
        m_session->LoadParameter("NekLinSysMaxIterations",
                                 sysKey.m_NekLinSysMaxIterations, 5000);
        m_session->LoadParameter("LinSysMaxStorage", sysKey.m_LinSysMaxStorage,
                                 100);
        m_session->LoadParameter("LinSysRelativeTolInNonlin",
                                 sysKey.m_NekLinSysTolerance, 1.0E-2);
        // Load required NonLinSys parameters:
        m_session->LoadParameter("NekNonlinSysMaxIterations",
                                 sysKey.m_NekNonlinSysMaxIterations, 100);
        m_session->LoadParameter("NekNonLinSysTolerance",
                                 sysKey.m_NekNonLinSysTolerance, 1.0E-09);
        m_session->LoadParameter("NonlinIterTolRelativeL2",
                                 sysKey.m_NonlinIterTolRelativeL2, 1.0E-6);

        m_nonlinsol = LibUtilities::GetNekNonlinSysIterFactory().CreateInstance(
            SolverType, m_session, m_comm, m_matDim, sysKey);

        m_NekSysOp.DefineNekSysResEval(&NonlinSysDemo::DoRhs, this);
        m_NekSysOp.DefineNekSysLhsEval(&NonlinSysDemo::DoLhs, this);
        m_nonlinsol->SetSysOperators(m_NekSysOp);
    }
    ~NonlinSysDemo()
    {
    }

    void DoSolve()
    {
        Array<OneD, NekDouble> pOutput(m_matDim, 0.9);

        int ntmpIts = m_nonlinsol->SolveSystem(m_matDim, pOutput, pOutput);

        int ndigits    = 6;
        int nothers    = 10;
        int nwidthcolm = nothers + ndigits - 1;
        cout << "ntmpIts = " << ntmpIts << endl
             << std::scientific << std::setw(nwidthcolm)
             << std::setprecision(ndigits - 1);

        string vars = "uvwx";
        for (int i = 0; i < m_matDim; ++i)
        {
            cout << "L 2 error (variable " << vars[i] << ") : " << pOutput[i]
                 << endl;
        }
    }

    void AllocateInitMatrix()
    {
        m_matDim = 2;
    }

    void DoLhs(InArrayType &inarray, OutArrayType &outarray,
               [[maybe_unused]] const bool &flag = false)
    {
        const Array<OneD, const NekDouble> refsol =
            m_nonlinsol->GetRefSolution();

        NekDouble x = refsol[0];
        NekDouble y = refsol[1];

        NekDouble f1 = 3.0 * x * x * inarray[0] + inarray[1];
        NekDouble f2 = 3.0 * y * y * inarray[1] - inarray[0];

        outarray[0] = f1;
        outarray[1] = f2;
    }

    void DoRhs(InArrayType &inarray, OutArrayType &outarray,
               [[maybe_unused]] const bool &flag = false)
    {
        ASSERTL1(m_matDim == inarray.size(),
                 "CoeffMat dim not equal to NekSys dim in DoRhs");
        NekDouble x  = inarray[0];
        NekDouble y  = inarray[1];
        NekDouble f1 = x * x * x + y - 1.0;
        NekDouble f2 = -x + y * y * y + 1.0;
        outarray[0]  = f1;
        outarray[1]  = f2;
    }

protected:
    int m_matDim;
    DNekMatSharedPtr m_matrix;
    Array<OneD, NekDouble> m_matDat;
    Array<OneD, NekDouble> m_SysRhs;
    NekNonlinSysIterSharedPtr m_nonlinsol;
    LibUtilities::NekSysOperators m_NekSysOp;
    LibUtilities::SessionReaderSharedPtr m_session;
    LibUtilities::CommSharedPtr m_comm;
    Array<OneD, int> m_map;
};

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;

    session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    session->InitSession();

    NonlinSysDemo nonlinsys(session, session->GetComm());

    nonlinsys.DoSolve();

    // Finalise communications
    session->Finalise();
    return 0;
}
