////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVortexInducedVelocity.cpp
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
//  Description: Computes velocity induced by vortex filaments.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessVortexInducedVelocity.h"
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar::FieldUtils
{

ModuleKey ProcessVortexInducedVelocity::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "vortexinducedvelocity"),
        ProcessVortexInducedVelocity::create,
        "Computes velocity induced by vortex filaments.");

ProcessVortexInducedVelocity::ProcessVortexInducedVelocity(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["vortex"] = ConfigOption(false, "vortexfilament.dat",
                                      "vortex filaments infomation");
    // format of vortexfilament.dat
    // uniform background flow [3]: u, v, w
    // vortex filament format [10] :
    // start point (x, y, z, infinity 1/finite 0);
    // end point (x, y, z, infinity 1/finite 0); sigma, circulation
    // 0., 0.4,1.,1; 3., 0.4,1.,0; 0.1,  1.
    // 3., 0.4,1.,0; 3.,-0.4,1.,0; 0.1,  1.
    // 0.,-0.4,1.,1; 3.,-0.4,1.,0; 0.1, -1.
    // 0, 0, 0
}

ProcessVortexInducedVelocity::~ProcessVortexInducedVelocity()
{
}

void ProcessVortexInducedVelocity::ParserDouble(const char *cstr, int n,
                                                std::vector<NekDouble> &value)
{
    value.clear();
    std::vector<int> digs;
    std::vector<int> dige;
    int i    = 0;
    int flag = 0; // digit chunk
    while (i < n)
    {
        if ((cstr[i] >= '0' && cstr[i] <= '9') || cstr[i] == '.' ||
            cstr[i] == 'e' || cstr[i] == 'E' || cstr[i] == '+' ||
            cstr[i] == '-')
        {
            if (flag == 0)
            {
                digs.push_back(i);
            }
            flag = 1;
        }
        else
        {
            if (flag == 1)
            {
                dige.push_back(i);
            }
            flag = 0;
        }
        if (cstr[i] == 0)
        {
            break;
        }
        ++i;
    }
    for (int j = 0; j < digs.size(); ++j)
    {
        std::string cuts(cstr + digs[j], dige[j] - digs[j]);
        value.push_back(std::stod(cuts));
    }
}

NekDouble ProcessVortexInducedVelocity::Smooths(NekDouble s, NekDouble sigma)
{
    if (s < NekConstants::kNekZeroTol)
    {
        return 0.;
    }
    else if (sigma < NekConstants::kNekZeroTol)
    {
        return 1. / (s * s);
    }
    else if (s < 1E-3 * sigma)
    {
        return 1. / (2. * sigma * sigma);
    }
    else
    {
        s = s * s;
        return (1. - exp(-s / (2. * sigma * sigma))) / s;
    }
}

void ProcessVortexInducedVelocity::v_Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }
    // read vortex filaments
    string vorfile = m_config["vortex"].as<string>();
    ifstream filament(vorfile.c_str());
    ASSERTL0(filament.is_open(), "Unable to open file " + vorfile);
    std::vector<NekDouble> Uniformflow(3, 0.);
    char buffer[1000];
    int count  = 0;
    int Nparam = 10;
    while (!filament.eof())
    {
        filament.getline(buffer, sizeof(buffer));
        std::vector<double> value;
        ParserDouble(buffer, sizeof(buffer), value);
        if (value.size() >= Nparam)
        {
            m_vortex.resize(count + 1);
            m_vortex[count].resize(Nparam);
            for (int i = 0; i < Nparam; ++i)
            {
                m_vortex[count][i] = value[i];
            }
            ++count;
        }
        else if (value.size() == 3)
        {
            Uniformflow = value;
        }
    }
    filament.close();
    // get coordinates
    int Np = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble>> coord(3);
    Array<OneD, Array<OneD, NekDouble>> vel(3);
    for (int i = 0; i < 3; ++i)
    {
        coord[i] = Array<OneD, NekDouble>(Np);
        vel[i]   = Array<OneD, NekDouble>(Np, Uniformflow[i]);
    }
    m_f->m_exp[0]->GetCoords(coord[0], coord[1], coord[2]);
    for (int v = 0; v < m_vortex.size(); ++v)
    {
        SpatialDomains::PointGeom pstart(3, 0, m_vortex[v][0], m_vortex[v][1],
                                         m_vortex[v][2]);
        SpatialDomains::PointGeom pend(3, 0, m_vortex[v][4], m_vortex[v][5],
                                       m_vortex[v][6]);
        SpatialDomains::PointGeom e;
        e.Sub(pend, pstart);
        NekDouble dist = sqrt(e.dot(e));
        if (dist < NekConstants::kNekZeroTol * NekConstants::kNekZeroTol)
        {
            continue;
        }
        for (int i = 0; i < 3; ++i)
        {
            e[i] /= dist;
        }
        NekDouble sigma = std::fabs(m_vortex[v][8]);
        NekDouble coeff = m_vortex[v][9] / (4. * M_PI);
        for (int p = 0; p < Np; ++p)
        {
            SpatialDomains::PointGeom r(3, 0, coord[0][p], coord[1][p],
                                        coord[2][p]);
            SpatialDomains::PointGeom e_r, s_r, cross;
            e_r.Sub(pend, r);
            s_r.Sub(pstart, r);
            cross.Mult(s_r, e);
            NekDouble s = sqrt(cross.dot(cross));
            if (s < NekConstants::kNekZeroTol)
            {
                continue;
            }
            NekDouble ls = -1., le = 1.;
            if (m_vortex[v][3] < 0.5)
            {
                ls = e.dot(s_r) / s;
                ls = ls / sqrt(ls * ls + 1.);
            }
            if (m_vortex[v][7] < 0.5)
            {
                le = e.dot(e_r) / s;
                le = le / sqrt(le * le + 1.);
            }
            NekDouble velcoeff = coeff * Smooths(s, sigma) * (le - ls);
            for (int i = 0; i < 3; ++i)
            {
                vel[i][p] += velcoeff * cross[i];
            }
        }
    }
    // add field
    m_f->m_variables.resize(3);
    vector<string> varsname = {"u", "v", "w"};
    MultiRegions::ExpListSharedPtr Exp;
    for (int n = 0; n < 3; ++n)
    {
        m_f->m_variables[n] = varsname[n];
        if (n < m_f->m_exp.size())
        {
            Vmath::Vcopy(Np, vel[n], 1, m_f->m_exp[n]->UpdatePhys(), 1);
        }
        else
        {
            Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
            Vmath::Vcopy(Np, vel[n], 1, Exp->UpdatePhys(), 1);
            m_f->m_exp.insert(m_f->m_exp.begin() + n, Exp);
        }
        m_f->m_exp[n]->FwdTransLocalElmt(vel[n], m_f->m_exp[n]->UpdateCoeffs());
    }
}
} // namespace Nektar::FieldUtils
