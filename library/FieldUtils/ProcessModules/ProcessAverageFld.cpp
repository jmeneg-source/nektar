////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessAverageFld.cpp
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
//  Description: Take averaging of input fields
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <boost/format.hpp>

#include "ProcessAverageFld.h"

namespace Nektar::FieldUtils
{

ModuleKey ProcessAverageFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "averagefld"), ProcessAverageFld::create,
        "averaging of several field files. Must specify inputfld and range of "
        "file numbers.");

ProcessAverageFld::ProcessAverageFld(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["inputfld"] =
        ConfigOption(false, "NotSet", "Fld file form which to average");
    m_config["range"] = ConfigOption(false, "NotSet", "range of file numbers");
}

ProcessAverageFld::~ProcessAverageFld()
{
}

void ProcessAverageFld::v_Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);
    ASSERTL0(m_config["inputfld"].as<string>().compare("NotSet") != 0,
             "Need to specify inputfld=file_%.fld ");
    ASSERTL0(m_config["range"].as<string>().compare("NotSet") != 0,
             "Need to specify range=ns,dn,ne[include]");

    std::string infilename = m_config["inputfld"].as<string>();
    std::vector<int> range;
    ASSERTL0(ParseUtils::GenerateVector(m_config["range"].as<string>(), range),
             "Failed to interpret range string");
    ASSERTL0(
        range.size() == 3,
        "range string should contain 3 values nstart, deltan, nend[include].");

    vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    vector<vector<double>> FieldData;
    NekDouble scale = 0.;
    for (int i = range[0]; i <= range[2]; i += range[1])
    {
        boost::format filenameformat(infilename);
        filenameformat % i;
        string inputfldname = filenameformat.str();
        vector<LibUtilities::FieldDefinitionsSharedPtr> tempFieldDef;
        vector<vector<double>> tempFieldData;

        if (m_f->m_graph)
        {
            const SpatialDomains::ExpansionInfoMap &expansions =
                m_f->m_graph->GetExpansionInfo();

            // if Range has been speficied it is possible to have a
            // partition which is empty so check this and return if
            // no elements present.

            if (!expansions.size())
            {
                return;
            }

            Array<OneD, int> ElementGIDs(expansions.size());

            int i = 0;
            for (auto &expIt : expansions)
            {
                ElementGIDs[i++] = expIt.second->m_geomPtr->GetGlobalID();
            }
            m_f->FieldIOForFile(inputfldname)
                ->Import(inputfldname, tempFieldDef, tempFieldData,
                         LibUtilities::NullFieldMetaDataMap, ElementGIDs);
        }
        else
        {
            m_f->FieldIOForFile(inputfldname)
                ->Import(inputfldname, tempFieldDef, tempFieldData,
                         LibUtilities::NullFieldMetaDataMap);
        }

        if (FieldDef.size() == 0)
        {
            FieldDef = tempFieldDef;
        }

        if (FieldData.size() < tempFieldData.size())
        {
            FieldData.resize(tempFieldData.size());
        }

        for (size_t i = 0; i < tempFieldData.size(); ++i)
        {
            if (FieldData[i].size() < tempFieldData[i].size())
            {
                FieldData[i].resize(tempFieldData[i].size(), 0.);
            }
            Vmath::Vadd(FieldData[i].size(), tempFieldData[i].data(), 1,
                        FieldData[i].data(), 1, FieldData[i].data(), 1);
        }
        if (m_f->m_comm->GetSpaceComm()->GetRank() == 0 && vm.count("verbose"))
        {
            std::cout << "File " << inputfldname << " processed." << std::endl;
        }
        scale += 1.;
    }
    if (scale < 0.5)
    {
        return; // no input file
    }
    scale = 1. / scale;

    // Filling expansion
    int numHomoDir = m_f->m_numHomogeneousDir;
    int nfields    = FieldDef[0]->m_fields.size();
    m_f->m_exp.resize(nfields);
    for (int i = 1; i < nfields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(numHomoDir);
    }
    m_f->m_variables = FieldDef[0]->m_fields;

    for (int j = 0; j < nfields; ++j)
    {
        // load new field
        for (int i = 0; i < FieldDef.size(); ++i)
        {
            m_f->m_exp[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                               m_f->m_variables[j],
                                               m_f->m_exp[j]->UpdateCoeffs());
        }
        Vmath::Smul(m_f->m_exp[j]->GetNcoeffs(), scale,
                    m_f->m_exp[j]->UpdateCoeffs(), 1,
                    m_f->m_exp[j]->UpdateCoeffs(), 1);
        m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(),
                                m_f->m_exp[j]->UpdatePhys());
    }
}
} // namespace Nektar::FieldUtils
