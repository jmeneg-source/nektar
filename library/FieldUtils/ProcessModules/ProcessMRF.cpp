////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessMRF.cpp
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
//  Description: Add transformed coordinates to field
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessMRF.h"
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar::FieldUtils
{
ModuleKey ProcessMRF::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "MRF"), ProcessMRF::create,
    "Add transformed coordinates to output file.");

ProcessMRF::ProcessMRF(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["vectors"] = ConfigOption(false, "NotSet", "Select variables");
}

ProcessMRF::~ProcessMRF()
{
}

void ProcessMRF::v_Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);
    ReadMRFData();

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // Determine dimensions
    m_spacedim = m_f->m_graph->GetMeshDimension() + m_f->m_numHomogeneousDir;

    // transform coordinates
    int nfields          = m_f->m_variables.size();
    int addfields        = m_spacedim;
    int npoints          = m_f->m_exp[0]->GetNpoints();
    string fieldNames[3] = {"xCoord", "yCoord", "zCoord"};
    for (int i = 0; i < addfields; ++i)
    {
        m_f->m_variables.push_back(fieldNames[i]);
    }
    m_f->m_exp.resize(nfields + addfields);
    vector<Array<OneD, NekDouble>> coords(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(npoints);
    }
    if (m_spacedim == 1)
    {
        m_f->m_exp[0]->GetCoords(coords[0]);
    }
    else if (m_spacedim == 2)
    {
        m_f->m_exp[0]->GetCoords(coords[0], coords[1]);
    }
    else
    {
        m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);
    }
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Sadd(npoints, -m_pivot[i], coords[i], 1, coords[i], 1);
    }
    TransformVector(coords);
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Sadd(npoints, m_pivot[i] + m_origin[i], coords[i], 1, coords[i],
                    1);
    }
    // Add new information to m_f
    for (int i = 0; i < addfields; ++i)
    {
        m_f->m_exp[nfields + i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, coords[i], 1,
                     m_f->m_exp[nfields + i]->UpdatePhys(), 1);
        m_f->m_exp[nfields + i]->FwdTransLocalElmt(
            coords[i], m_f->m_exp[nfields + i]->UpdateCoeffs());
    }

    // tranform vectors
    vector<Array<OneD, NekDouble>> data(m_spacedim);
    vector<int> vars;
    if (m_config["vectors"].as<string>().compare("NotSet"))
    {
        ParseUtils::GenerateVariableVector(m_config["vectors"].as<string>(),
                                           m_f->m_variables, vars);
    }
    ASSERTL0(vars.size() % m_spacedim == 0,
             "The number of vector variables is not divisible by the space "
             "dimension.");
    for (int i = 0; i < vars.size() / m_spacedim; ++i)
    {
        for (int d = 0; d < m_spacedim; ++d)
        {
            data[d] = m_f->m_exp[vars[m_spacedim * i + d]]->UpdatePhys();
        }
        TransformVector(data);
        for (int d = 0; d < m_spacedim; ++d)
        {
            int nv = vars[m_spacedim * i + d];
            m_f->m_exp[nv]->FwdTransLocalElmt(data[d],
                                              m_f->m_exp[nv]->UpdateCoeffs());
        }
    }
}

void ProcessMRF::ReadMRFData()
{
    vector<string> strOrigin = {"X", "Y", "Z"};
    vector<string> strTheta  = {"Theta_x", "Theta_y", "Theta_z"};
    vector<string> strPivot  = {"X0", "Y0", "Z0"};
    m_origin.resize(3, 0.);
    m_pivot.resize(3, 0.);
    for (size_t i = 0; i < strOrigin.size(); ++i)
    {
        if (m_f->m_fieldMetaDataMap.count(strOrigin[i]))
        {
            m_origin[i] = stod(m_f->m_fieldMetaDataMap[strOrigin[i]]);
        }
    }
    for (size_t i = 0; i < strTheta.size(); ++i)
    {
        if (m_f->m_fieldMetaDataMap.count(strTheta[i]))
        {
            if (m_theta.size() < 3)
            {
                m_theta.resize(3, 0.);
            }
            m_theta[i] = stod(m_f->m_fieldMetaDataMap[strTheta[i]]);
        }
    }
    for (size_t i = 0; i < strPivot.size(); ++i)
    {
        if (m_f->m_fieldMetaDataMap.count(strPivot[i]))
        {
            m_pivot[i] = stod(m_f->m_fieldMetaDataMap[strPivot[i]]);
        }
    }
}

void ProcessMRF::TransformVector(vector<Array<OneD, NekDouble>> &data)
{
    if (m_theta.size() < 3)
    {
        return;
    }
    if (data.size() < m_spacedim)
    {
        NEKERROR(ErrorUtil::efatal,
                 "data size is small than the space dimension.");
    }

    int npoint   = data[0].size();
    int dim      = 2;
    NekDouble sz = sin(m_theta[2]), cz = cos(m_theta[2]);
    Array<OneD, Array<OneD, NekDouble>> tmp(dim);
    for (int i = 0; i < dim; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>(npoint);
        Vmath::Vcopy(npoint, data[i], 1, tmp[i], 1);
    }
    Vmath::Svtsvtp(npoint, cz, tmp[0], 1, -sz, tmp[1], 1, data[0], 1);
    Vmath::Svtsvtp(npoint, sz, tmp[0], 1, cz, tmp[1], 1, data[1], 1);
}
} // namespace Nektar::FieldUtils
