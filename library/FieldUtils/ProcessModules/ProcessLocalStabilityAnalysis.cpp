////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLocalStabilityAnalysis.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Local linear stability analysis of compressible flow.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

#include "ProcessLocalStabilityAnalysis.h"

#define _DEBUG_

using namespace std;

namespace Nektar::FieldUtils
{

ModuleKey ProcessLocalStabilityAnalysis::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "localStabilityAnalysis"),
        ProcessLocalStabilityAnalysis::create,
        "Perform local (linear) stability analysis (LST).");

ProcessLocalStabilityAnalysis::ProcessLocalStabilityAnalysis(FieldSharedPtr f)
    : ProcessBoundaryExtract(f)
{
    f->m_writeBndFld = false; // turned on in upstream ProcessBoundaryExtract

    m_config["opt"] =
        ConfigOption(false, "6", "Option for quantity to compute.\n \
                              6 - scan frequency\n \
                              7 - scan spanwise wavenumber\n \
                              8 - scan both frequency and spanwise wavenumber");
    m_config["finVal1"] =
        ConfigOption(false, "0.1", "Final value for the 1st quantity.");
    m_config["numStep1"] =
        ConfigOption(false, "10", "Number of steps to scan the 1st range.");
    m_config["finVal2"] = ConfigOption(
        false, "0.1", "Final value for the 2nd quantity (optional).");
    m_config["numStep2"] = ConfigOption(
        false, "10", "Number of steps to scan the 2nd range (optional).");
}

ProcessLocalStabilityAnalysis::~ProcessLocalStabilityAnalysis()
{
}

void ProcessLocalStabilityAnalysis::GenerateInputFiles()
{
    // Update settins from the session file
    if (m_f->m_session->DefinesParameter("Pr"))
    {
        m_PRANDTL = m_f->m_session->GetParameter("Pr");
    }

    if (m_f->m_session->DefinesParameter("Gamma"))
    {
        m_GAMMA = m_f->m_session->GetParameter("Gamma");
    }

    if (m_f->m_session->DefinesParameter("RGAS"))
    {
        m_RGAS = m_f->m_session->GetParameter("RGAS");
    }

    if (m_f->m_session->DefinesParameter("ITOLRAY"))
    {
        m_ITOLRAY = m_f->m_session->GetParameter("ITOLRAY");
    }

    if (m_f->m_session->DefinesParameter("NGLOBAL"))
    {
        m_NGLOBAL = m_f->m_session->GetParameter("NGLOBAL");
    }

    if (m_f->m_session->DefinesParameter("NLOCAL"))
    {
        m_NLOCAL = m_f->m_session->GetParameter("NLOCAL");
    }

    if (m_f->m_session->DefinesParameter("BETALIN"))
    {
        m_BETALIN = m_f->m_session->GetParameter("BETALIN");
    }

    ASSERTL0(m_f->m_session->DefinesParameter("FREQLIN"),
             "The disturbance frequency for local linear stability analysis "
             "must be defined.");
    m_FREQLIN = m_f->m_session->GetParameter("FREQLIN");

    if (m_f->m_session->DefinesParameter("ITSE"))
    {
        m_ITSE = m_f->m_session->GetParameter("ITSE");
    }

    if (m_f->m_session->DefinesParameter("INEUTRL"))
    {
        m_INEUTRL = m_f->m_session->GetParameter("INEUTRL");
    }

    ofstream ofs;
    ofs.open("PSEREAD.NAM", ios::out);

    // Write $PSEUSER tag
    ofs << "$PSEUSER\n"
        << "NSTART=" << m_NSTART << "\n"
        << "ITSE=" << m_ITSE << "\n"
        << "INEUTRL=" << m_INEUTRL << "\n"
        << "BETALIN=" << m_BETALIN << "\n"
        << "FREQLIN=" << m_FREQLIN << "\n"
        << "$END\n"
        << endl;
    // Write $PSEUSER tag
    ofs << "$PSENONL\n$END\n" << endl;

    // Weite $PSECNTL tag
    ofs << "$PSECNTL\n"
        << "NGLOBAL=" << m_NGLOBAL << "\n"
        << "NLOCAL=" << m_NLOCAL << "\n"
        << "ITOLRAY=" << m_ITOLRAY << "\n"
        << "$END\n"
        << endl;

    // Write $PSEOUT tag
    ofs << "$PSEOUT\n"
        << "LABEL=" << m_LABEL << "\n"
        << "ITSPRES=" << m_ITSPRES << "\n"
        << "IEIGL=" << m_IEIGL << "\n"
        << "$END\n"
        << endl;

    // Write $PSEGAS tag
    ofs << "$PSEGAS\n"
        << "RGAS=" << m_RGAS << "\n"
        << "GAMMA=" << m_GAMMA << "\n"
        << "PRANDTL=" << m_PRANDTL << "\n"
        << "$END\n"
        << endl;

    ofs.close();
}

void ProcessLocalStabilityAnalysis::v_Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::v_Process(vm);

    // Initialize sampling parameters
    const long int option       = m_config["opt"].as<long int>(); // 6,7,8
    const NekDouble finalValue1 = m_config["finVal1"].as<NekDouble>();
    const long int numStep1     = m_config["numStep1"].as<long int>();
    const NekDouble finalValue2 = m_config["finVal2"].as<NekDouble>();
    const long int numStep2     = m_config["numStep2"].as<long int>();

    // Create PSEREAD.NAMgit d
    GenerateInputFiles();

    // call_hello(); // test calling
    call_lst(option, &finalValue1, numStep1, &finalValue2, numStep2);
}
} // namespace Nektar::FieldUtils
