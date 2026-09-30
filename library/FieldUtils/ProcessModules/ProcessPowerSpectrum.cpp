////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessPowerSpectrum.cpp
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
//  Description: Find the spanwise energy distribution of a 3DH1D field.
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessPowerSpectrum.h"
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <iostream>
#include <string>

using namespace std;

namespace Nektar::FieldUtils
{

ModuleKey ProcessPowerSpectrum::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "powerspectrum"),
        ProcessPowerSpectrum::create,
        "Output power spectrum at given regions from a 3DH1D expansion.");

ProcessPowerSpectrum::ProcessPowerSpectrum(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["box"] =
        ConfigOption(false, "NotSet",
                     "Specify a rectangular box by box=xmin,xmax,ymin,ymax");
    m_config["vars"] = ConfigOption(false, "NotSet", "signal variable");
}

ProcessPowerSpectrum::~ProcessPowerSpectrum()
{
}

void ProcessPowerSpectrum::v_Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    ASSERTL0(m_f->m_numHomogeneousDir == 1,
             "ProcessPowerSpectrum only works for Homogeneous1D fields.");
    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    std::set<int> variables;
    if (m_config["vars"].as<string>().compare("NotSet") != 0)
    {
        ParseUtils::GenerateVariableSet(m_config["vars"].as<string>(),
                                        m_f->m_variables, variables);
    }
    ASSERTL0(variables.size() >= 1,
             "At least one variable should be selected.");
    if (variables.empty())
    {
        return;
    }
    std::vector<NekDouble> box;
    if (m_config["box"].as<string>().compare("NotSet") != 0)
    {
        ASSERTL0(ParseUtils::GenerateVector(m_config["box"].as<string>(), box),
                 "Failed to interpret box string");
        ASSERTL0(
            box.size() >= 4,
            "box string should contain at least 4 values xmin,xmax,ymin,ymax");
    }
    box.resize(4, 0.);
    if (m_f->m_comm->GetSpaceComm()->GetRank() == 0)
    {

        cout << "Processing frequency spectra of variables (";
        for (auto nsignal : variables)
        {
            cout << m_f->m_variables[nsignal] << ", ";
        }
        cout << ") in box [" << box[0] << ", " << box[1] << "]X[" << box[2]
             << ", " << box[3] << "].\n";
    }

    for (auto nsignal : variables)
    {
        m_f->m_exp[nsignal]->SetWaveSpace(true);
        m_f->m_exp[nsignal]->BwdTrans(m_f->m_exp[nsignal]->GetCoeffs(),
                                      m_f->m_exp[nsignal]->UpdatePhys());
    }

    int nplanes                           = m_f->m_exp[0]->GetZIDs().size();
    NekDouble lz                          = m_f->m_exp[0]->GetHomoLen();
    MultiRegions::ExpListSharedPtr plane0 = m_f->m_exp[0]->GetPlane(0);
    int ntot                              = plane0->GetNpoints();
    Array<OneD, NekDouble> masks(ntot, 0.);
    for (size_t i = 0; i < plane0->GetExpSize(); ++i)
    {
        LocalRegions::ExpansionSharedPtr exp = plane0->GetExp(i);
        SpatialDomains::Geometry *geom       = exp->GetGeom();
        int nv                               = geom->GetNumVerts();
        NekDouble gc[3]                      = {0., 0., 0.};
        NekDouble gct[3]                     = {0., 0., 0.};
        for (size_t j = 0; j < nv; ++j)
        {
            SpatialDomains::PointGeom *vertex = geom->GetVertex(j);
            vertex->GetCoords(gct[0], gct[1], gct[2]);
            gc[0] += gct[0] / NekDouble(nv);
            gc[1] += gct[1] / NekDouble(nv);
        }
        if (box[0] <= gc[0] && gc[0] <= box[1] && box[2] <= gc[1] &&
            gc[1] <= box[3])
        {
            Vmath::Fill(exp->GetTotPoints(), 1.,
                        &masks[plane0->GetPhys_Offset(i)], 1);
        }
    }
    NekDouble area = m_f->m_exp[0]->GetPlane(0)->Integral(masks);
    m_f->m_comm->GetSpaceComm()->AllReduce(area, LibUtilities::ReduceMax);
    if (area <= NekConstants::kNekZeroTol)
    {
        area = 1.;
    }

    Array<OneD, NekDouble> value(nplanes / 2, 0.);
    for (size_t m = 0; m < value.size(); ++m)
    {
        for (auto nsignal : variables)
        {
            Array<OneD, NekDouble> wavecoef(ntot, 0.);
            Vmath::Vvtvp(ntot, m_f->m_exp[nsignal]->GetPlane(2 * m)->GetPhys(),
                         1, m_f->m_exp[nsignal]->GetPlane(2 * m)->GetPhys(), 1,
                         wavecoef, 1, wavecoef, 1);
            Vmath::Vvtvp(ntot,
                         m_f->m_exp[nsignal]->GetPlane(2 * m + 1)->GetPhys(), 1,
                         m_f->m_exp[nsignal]->GetPlane(2 * m + 1)->GetPhys(), 1,
                         wavecoef, 1, wavecoef, 1);
            Vmath::Vmul(ntot, wavecoef, 1, masks, 1, wavecoef, 1);
            if (m > 0)
            {
                Vmath::Smul(ntot, 0.25, wavecoef, 1, wavecoef, 1);
            }
            value[m] += m_f->m_exp[0]->GetPlane(2 * m)->Integral(wavecoef);
        }
    }
    m_f->m_comm->GetSpaceComm()->AllReduce(value, LibUtilities::ReduceMax);
    if (m_f->m_comm->GetSpaceComm()->GetRank() == 0)
    {
        NekDouble energy = 0.;
        ofstream powersp("spectra.dat");
        powersp << "variables = lambda beta energy" << endl;
        for (size_t m = 1; m < value.size(); ++m)
        {
            powersp << (lz / m) << " " << (2. * M_PI * m / lz) << " "
                    << value[m] << endl;
            energy += value[m] * value[m];
        }
        powersp.close();
        if (vm.count("error"))
        {
            cout << "L 2 error (variable E) : " << energy << endl;
        }
    }

    int nfields = m_f->m_variables.size();
    m_f->m_exp.resize(nfields + 1);
    m_f->m_variables.push_back("regions");
    m_f->m_exp[nfields] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    for (size_t m = 0; m < nplanes; ++m)
    {
        Vmath::Vcopy(ntot, masks, 1,
                     m_f->m_exp[nfields]->GetPlane(m)->UpdatePhys(), 1);
    }
    m_f->m_exp[nfields]->FwdTrans(m_f->m_exp[nfields]->GetPhys(),
                                  m_f->m_exp[nfields]->UpdateCoeffs());
}
} // namespace Nektar::FieldUtils
