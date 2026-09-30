///////////////////////////////////////////////////////////////////////////////
//
// File: FilterAverageFields.cpp
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
// Description: Average solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterAverageFields.h>

namespace Nektar::SolverUtils
{
std::string FilterAverageFields::cmdSetStartFilterFileNum =
    LibUtilities::SessionReader::RegisterCmdLineArgument(
        "set-filter-averagefield-start-number", "",
        "Set the starting number of the average field filter file number.");

std::string FilterAverageFields::className =
    GetFilterFactory().RegisterCreatorFunction("AverageFields",
                                               FilterAverageFields::create);

FilterAverageFields::FilterAverageFields(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::shared_ptr<EquationSystem> &pEquation, const ParamMap &pParams)
    : FilterFieldConvert(pSession, pEquation, pParams)
{
    // Load sampling frequency
    auto it = pParams.find("SampleFrequency");
    if (it == pParams.end())
    {
        m_sampleFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_sampleFrequency = round(equ.Evaluate());
    }
}

void FilterAverageFields::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Initialise output arrays
    FilterFieldConvert::v_Initialise(pFields, time);

    if (m_session->DefinesCmdLineArgument(
            "set-filter-averagefield-start-number"))
    {
        m_outputIndex = std::stoi(m_session->GetCmdLineArgument<std::string>(
            "set-filter-averagefield-start-number"));
    }
}

void FilterAverageFields::v_ProcessSample(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    [[maybe_unused]] const NekDouble &time)
{
    for (int n = 0; n < m_outFields.size(); ++n)
    {
        Vmath::Vadd(m_outFields[n].size(), fieldcoeffs[n], 1, m_outFields[n], 1,
                    m_outFields[n], 1);
    }
}

void FilterAverageFields::v_PrepareOutput(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    [[maybe_unused]] const NekDouble &time)
{
    m_fieldMetaData["NumberOfFieldDumps"] = std::to_string(m_numSamples);
}

NekDouble FilterAverageFields::v_GetScale()
{
    return 1.0 / m_numSamples;
}

} // namespace Nektar::SolverUtils
