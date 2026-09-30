///////////////////////////////////////////////////////////////////////////////
//
// File: FilterCheckpoint.cpp
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
// Description: Outputs solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterCheckpoint.h>

namespace Nektar::SolverUtils
{
std::string FilterCheckpoint::cmdSetStartFilterFileNum =
    LibUtilities::SessionReader::RegisterCmdLineArgument(
        "set-filter-checkpoint-start-number", "",
        "Set the starting number of the check point filter file number.");

std::string FilterCheckpoint::className =
    GetFilterFactory().RegisterCreatorFunction("Checkpoint",
                                               FilterCheckpoint::create);

FilterCheckpoint::FilterCheckpoint(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::shared_ptr<EquationSystem> &pEquation, const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    // OutputFile
    std::string ext = ".chk";
    m_outputFile    = Filter::SetupOutput(ext, pParams);

    // OutputFrequency
    auto it = pParams.find("OutputFrequency");
    ASSERTL0(it != pParams.end(), "Missing parameter 'OutputFrequency'.");
    LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
    m_outputFrequency = round(equ.Evaluate());

    // Time after which we need to write checkfiles
    it = pParams.find("OutputStartTime");
    if (it == pParams.end())
    {
        m_outputStartTime = 0;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_outputStartTime = equ.Evaluate();
    }

    m_fld = LibUtilities::FieldIO::CreateDefault(pSession);
}

void FilterCheckpoint::v_Initialise(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    [[maybe_unused]] const NekDouble &time)
{
    m_index       = 0;
    m_outputIndex = 0;

    // Load information for outputIndex
    auto equationSys = m_equ.lock();
    if (equationSys->GetFieldMetaDataMap().count("FilterFileNum"))
    {
        m_outputIndex =
            atoi(equationSys->GetFieldMetaDataMap()["FilterFileNum"].c_str());
    }

    if (m_session->DefinesCmdLineArgument("set-filter-checkpoint-start-number"))
    {
        m_outputIndex = std::stoi(m_session->GetCmdLineArgument<std::string>(
            "set-filter-checkpoint-start-number"));
    }

    if (m_updateOnInitialise)
    {
        v_Update(pFields, time);
    }
}

void FilterCheckpoint::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{

    if (m_index++ % m_outputFrequency > 0 ||
        (time - m_outputStartTime) < -1.0e-07)
    {
        return;
    }

    // get the extension of file name
    std::string ext = fs::path(m_outputFile).extension().string();

    std::stringstream vTmpFilename;
    std::string vOutputFilename;
    vTmpFilename << fs::path(m_outputFile).replace_extension("").string() << "_"
                 << m_outputIndex << ".chk";
    vOutputFilename = Filter::SetupOutput(ext, vTmpFilename.str());

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        pFields[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());

    // copy Data into FieldData and set variable
    for (int j = 0; j < pFields.size(); ++j)
    {
        for (int i = 0; i < FieldDef.size(); ++i)
        {
            // Could do a search here to find correct variable
            FieldDef[i]->m_fields.push_back(m_session->GetVariable(j));
            pFields[0]->AppendFieldData(FieldDef[i], FieldData[i],
                                        pFields[j]->UpdateCoeffs());
        }
    }
    auto equationSys = m_equ.lock();
    LibUtilities::FieldMetaDataMap fieldMetaDataMap(
        equationSys->GetFieldMetaDataMap());
    fieldMetaDataMap["FilterFileNum"] = std::to_string(m_outputIndex);
    m_fld->Write(vOutputFilename, FieldDef, FieldData, fieldMetaDataMap);
    m_outputIndex++;
}

void FilterCheckpoint::v_Finalise(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    [[maybe_unused]] const NekDouble &time)
{
}

bool FilterCheckpoint::v_IsTimeDependent()
{
    return true;
}
} // namespace Nektar::SolverUtils
