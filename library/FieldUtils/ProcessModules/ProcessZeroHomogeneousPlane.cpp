////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessZeroHomogeneousPlane.cpp
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
//  Description: Zero the homogeneous plane in wavespace of a 3DH1D field to
//  visualise the fluctuation.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessZeroHomogeneousPlane.h"

namespace Nektar::FieldUtils
{

ModuleKey ProcessZeroHomogeneousPlane::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "zero-homo-plane"),
        ProcessZeroHomogeneousPlane::create,
        "Extracts a plane from a 3DH1D expansion, requires planeid to be "
        "defined.");

ProcessZeroHomogeneousPlane::ProcessZeroHomogeneousPlane(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["planeid"] = ConfigOption(
        false, "0",
        "plane id in wavespace to be set zero. Default planeid is zero");
}

ProcessZeroHomogeneousPlane::~ProcessZeroHomogeneousPlane()
{
}

void ProcessZeroHomogeneousPlane::v_Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    if (m_f->m_numHomogeneousDir != 1)
    {
        NEKERROR(ErrorUtil::efatal,
                 "ProcessZeroHomogeneousPlane only works for Homogeneous1D.");
    }

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int planeid = m_config["planeid"].as<int>();
    int nfields = m_f->m_variables.size();

    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    // Look for correct plane (because of parallel case)
    int plane = -1;
    for (int i = 0; i < m_f->m_exp[0]->GetZIDs().size(); ++i)
    {
        if (m_f->m_exp[0]->GetZIDs()[i] == planeid)
        {
            plane = i;
        }
    }

    if (plane != -1)
    {
        for (int s = 0; s < nstrips; ++s)
        {
            for (int i = 0; i < nfields; ++i)
            {
                int n = s * nfields + i;

                // Find the memory location, set zeros
                int Ncoeff = m_f->m_exp[n]->GetPlane(plane)->GetNcoeffs();
                Vmath::Zero(Ncoeff,
                            m_f->m_exp[n]->GetPlane(plane)->UpdateCoeffs(), 1);

                m_f->m_exp[n]->BwdTrans(m_f->m_exp[n]->GetCoeffs(),
                                        m_f->m_exp[n]->UpdatePhys());
            }
        }
    }
    else
    {
        NEKERROR(ErrorUtil::efatal, "Plane not found.");
    }
}
} // namespace Nektar::FieldUtils
