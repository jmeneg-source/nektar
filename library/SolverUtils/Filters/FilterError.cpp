///////////////////////////////////////////////////////////////////////////////
//
// File: FilterError.cpp
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
// Description: Outputs errors during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterError.h>

namespace Nektar::SolverUtils
{
std::string FilterError::className =
    GetFilterFactory().RegisterCreatorFunction("Error", FilterError::create);

FilterError::FilterError(const LibUtilities::SessionReaderSharedPtr &pSession,
                         const std::shared_ptr<EquationSystem> &pEquation,
                         const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    std::string outName;
    std::string ext = ".err";
    outName         = Filter::SetupOutput(ext, pParams);

    // Lock equation system pointer
    auto equationSys = m_equ.lock();
    ASSERTL0(equationSys, "Weak pointer expired");

    m_numVariables = equationSys->GetNvariables();

    // Check if we have homogeneous expansion
    // If yes, disable H1 error norm as it is not implemented for homogeneous
    // expansions
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isHomogeneous1D, false);

    m_comm = pSession->GetComm();
    if (m_comm->GetRank() == 0)
    {
        m_outFile.open(outName);
        ASSERTL0(m_outFile.good(), "Unable to open: '" + outName + "'");
        m_outFile.setf(std::ios::scientific, std::ios::floatfield);

        m_outFile << "Time";
        for (size_t i = 0; i < m_numVariables; ++i)
        {
            std::string varName = equationSys->GetVariable(i);
            m_outFile << " " + varName + "_L2" << " " + varName + "_Linf";

            if (!m_isHomogeneous1D)
            {
                m_outFile << " " + varName + "_H1";
            }
        }

        m_outFile << std::endl;
    }

    // OutputFrequency
    auto it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Empty parameter 'OutputFrequency'.");
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // ConsoleOutput
    it = pParams.find("ConsoleOutput");
    if (it == pParams.end())
    {
        m_consoleOutput = false;
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Empty parameter 'ConsoleOutput'.");
        ASSERTL0(it->second == "0" || it->second == "1",
                 "Parameter 'ConsoleOutput' can only be '0' or '1'.");
        m_consoleOutput = boost::lexical_cast<bool>(it->second);
    }
}

void FilterError::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    v_Update(pFields, time);
}

void FilterError::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Check whether output frequency matches
    if (m_index++ % m_outputFrequency > 0)
    {
        return;
    }

    // Output to file L2 and Linf error for each variable
    if (m_comm->GetRank() == 0)
    {
        m_outFile << time;
    }

    // Lock equation system pointer
    auto equationSys = m_equ.lock();
    ASSERTL0(equationSys, "Weak pointer expired");

    for (size_t i = 0; i < m_numVariables; ++i)
    {
        // Evaluate "ExactSolution" function, or zero array
        Array<OneD, NekDouble> exactsoln(pFields[i]->GetTotPoints(), 0.0);
        equationSys->EvaluateExactSolution(i, exactsoln, time);

        NekDouble vL2Error   = equationSys->L2Error(i, exactsoln);
        NekDouble vLinfError = equationSys->LinfError(i, exactsoln);

        if (m_comm->GetRank() == 0)
        {
            m_outFile << " " << vL2Error << " " << vLinfError;

            if (m_consoleOutput)
            {
                std::cout << "L 2 error (variable "
                          << equationSys->GetVariable(i) << ") : " << vL2Error
                          << std::endl;
                std::cout << "L inf error (variable "
                          << equationSys->GetVariable(i) << ") : " << vLinfError
                          << std::endl;
            }
        }

        // Only evaluate H1 for non homogeneous case
        if (!m_isHomogeneous1D)
        {
            NekDouble vH1Error = equationSys->H1Error(i, exactsoln);
            if (m_comm->GetRank() == 0)
            {
                m_outFile << " " << vH1Error;

                if (m_consoleOutput)
                {
                    std::cout << "H 1 error (variable "
                              << equationSys->GetVariable(i)
                              << ") : " << vH1Error << std::endl;
                }
            }
        }
    }

    if (m_comm->GetRank() == 0)
    {
        m_outFile << std::endl;
    }
}

void FilterError::v_Finalise(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    [[maybe_unused]] const NekDouble &time)
{
    if (m_comm->GetRank() == 0)
    {
        m_outFile.close();
    }
}

bool FilterError::v_IsTimeDependent()
{
    return true;
}
} // namespace Nektar::SolverUtils
