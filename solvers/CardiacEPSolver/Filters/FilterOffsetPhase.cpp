///////////////////////////////////////////////////////////////////////////////
//
// File FilterOffsetPhase.cpp
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
// Description: Outputs phase fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <CardiacEPSolver/Filters/FilterOffsetPhase.h>

namespace Nektar
{
std::string FilterOffsetPhase::className =
    SolverUtils::GetFilterFactory().RegisterCreatorFunction(
        "OffsetPhase", FilterOffsetPhase::create);

FilterOffsetPhase::FilterOffsetPhase(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::shared_ptr<SolverUtils::EquationSystem> &pEquation,
    const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    std::string ext = "";
    m_outputFile    = Filter::SetupOutput(ext, pParams);

    ASSERTL0(pParams.find("OutputFrequency") != pParams.end(),
             "Missing parameter 'OutputFrequency'.");
    LibUtilities::Equation equ(m_session->GetInterpreter(),
                               pParams.find("OutputFrequency")->second);
    m_outputFrequency = floor(equ.Evaluate());

    if (pParams.find("MeanV") != pParams.end())
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(),
                                   pParams.find("MeanV")->second);
        meanV = equ.Evaluate();
    }

    m_outputIndex = 0;
    m_index       = 0;
    m_fld         = LibUtilities::FieldIO::CreateDefault(pSession);
}

FilterOffsetPhase::~FilterOffsetPhase()
{
}

void FilterOffsetPhase::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const NekDouble &time)
{
    m_index       = 0;
    m_outputIndex = 0;

    oldV = pFields[0]->GetPhys(); // Allocate storage for past data

    v_Update(pFields, 0.0);
}

void FilterOffsetPhase::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    if (m_index++ % m_outputFrequency > 0)
    {
        return;
    }

    std::stringstream vTmpFilename;
    std::string vOutputFilename;
    std::string ext = ".chk";

    vTmpFilename << fs::path(m_outputFile).replace_extension("").string() << "_"
                 << m_outputIndex << ext;
    vOutputFilename = Filter::SetupOutput(ext, vTmpFilename.str());

    // Get FieldDefinitions
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        pFields[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());

    Array<OneD, NekDouble> currentV = pFields[0]->GetPhys();
    Array<OneD, NekDouble> phase(currentV.size());
    Array<OneD, NekDouble> phase_coeff(pFields[0]->GetNcoeffs());

    // Calculate phase angle in physical space:
    for (int i = 0; i < currentV.size(); ++i)
    {
        phase[i] = atan2(currentV[i] - meanV,
                         oldV[i] - meanV); // subtract mean -57.2986
    }

    // Convert to coefficient space:
    pFields[0]->FwdTransLocalElmt(phase, phase_coeff);

    // for each item in FieldDef
    //        push "phase" to FieldDef[i]->m_fields
    //        call AppendFieldData(FieldDef[i], FieldData[i], data)
    for (int i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("phase");
        pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], phase_coeff);
    }

    // Save currentV
    oldV = currentV;

    // Update time in field info if required
    LibUtilities::FieldMetaDataMap fieldMetaDataMap;
    fieldMetaDataMap["Time"] = boost::lexical_cast<std::string>(time);

    m_fld->Write(vOutputFilename, FieldDef, FieldData, fieldMetaDataMap);
    m_outputIndex++;
}

void FilterOffsetPhase::v_Finalise(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    [[maybe_unused]] const NekDouble &time)
{
}

bool FilterOffsetPhase::v_IsTimeDependent()
{
    return true;
}
} // namespace Nektar
