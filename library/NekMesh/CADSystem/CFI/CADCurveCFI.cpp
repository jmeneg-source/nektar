////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurveCFI.cpp
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
//  Description: CAD object curve methods.
//
////////////////////////////////////////////////////////////////////////////////

#include "CADCurveCFI.h"

using namespace std;

namespace Nektar
{
namespace NekMesh
{

std::string CADCurveCFI::key = GetCADCurveFactory().RegisterCreatorFunction(
    "cfi", CADCurveCFI::create, "CADCurveCFI");

void CADCurveCFI::Initialise(int i, cfi::Line *in, NekDouble s)
{
    m_cfiEdge = in;
    m_scal    = s;
    m_length  = m_cfiEdge->calcLength() * m_scal;

    m_id = i;
}

NekDouble CADCurveCFI::tAtArcLength(NekDouble s)
{
    s /= m_scal;
    auto bds     = GetBounds();
    NekDouble dt = (bds[1] - bds[0]) / 1000;

    NekDouble t   = bds[0];
    NekDouble len = 0.0;

    while (len <= s)
    {
        auto drdt1 = D2(t);
        t += dt;
        auto drdt2 = D2(t);

        NekDouble mag1 = sqrt(drdt1[3] * drdt1[3] + drdt1[4] * drdt1[4] +
                              drdt1[5] * drdt1[5]);
        NekDouble mag2 = sqrt(drdt2[3] * drdt2[3] + drdt2[4] * drdt2[4] +
                              drdt2[5] * drdt2[5]);

        len += (mag1 + mag2) / 2.0 * dt;
    }

    return t - dt;
}

NekDouble CADCurveCFI::loct(std::array<NekDouble, 3> xyz, NekDouble > xyz,
                            NekDouble &t)
{
    cfi::Position p;
    p.x = xyz[0] / m_scal;
    p.y = xyz[1] / m_scal;
    p.z = xyz[2] / m_scal;

    boost::optional<cfi::Projected<double>> pj = m_cfiEdge->calcTFromXYZ(p, -1);

    t = pj.value().parameters;

    return pj.value().distance * m_scal;
}

NekDouble CADCurveCFI::Length(NekDouble ti, NekDouble tf)
{
    auto bds     = GetBounds();
    NekDouble dt = (bds[1] - bds[0]) / 1000;

    NekDouble t   = ti;
    NekDouble len = 0.0;

    while (t <= tf)
    {
        auto drdt1 = D2(t);
        t += dt;
        auto drdt2 = D2(t);

        NekDouble mag1 = sqrt(drdt1[3] * drdt1[3] + drdt1[4] * drdt1[4] +
                              drdt1[5] * drdt1[5]);
        NekDouble mag2 = sqrt(drdt2[3] * drdt2[3] + drdt2[4] * drdt2[4] +
                              drdt2[5] * drdt2[5]);

        len += (mag1 + mag2) / 2.0 * dt;
    }

    return len * m_scal;
}

std::array<NekDouble, 3> CADCurveCFI::P(NekDouble t)
{
    cfi::Position p = m_cfiEdge->calcXYZAtT(t);

    return {p.x * m_scal, p.y * m_scal, p.z * m_scal};
}

void CADCurveCFI::P(NekDouble t, NekDouble &x, NekDouble &y, NekDouble &z)
{
    cfi::Position p = m_cfiEdge->calcXYZAtT(t);

    x = p.x * m_scal;
    y = p.y * m_scal;
    z = p.z * m_scal;
}

std::array<NekDouble, 9> CADCurveCFI::D2(NekDouble t)
{
    vector<cfi::DerivativeList> *d = m_cfiEdge->calcDerivAtT(t);
    cfi::Position p                = m_cfiEdge->calcXYZAtT(t);

    cfi::DerivativeList d1 = d->at(0);
    cfi::DerivativeList d2 = d->at(1);

    return {p.x * m_scal,
            p.y * m_scal,
            p.z * m_scal,
            d1.getDeriv(0) * m_scal,
            d1.getDeriv(1) * m_scal,
            d1.getDeriv(2) * m_scal,
            d2.getDeriv(0) * m_scal,
            d2.getDeriv(0) * m_scal,
            d2.getDeriv(0) * m_scal};
}

std::array<NekDouble, 2> CADCurveCFI::GetBounds()
{
    std::array<NekDouble, 2> t;
    t[0] = 0.0;
    t[1] = 1.0;

    return t;
}

void CADCurveCFI::GetBounds(NekDouble &tmin, NekDouble &tmax)
{
    tmin = 0.0;
    tmax = 1.0;
}

std::array<NekDouble, 6> CADCurveCFI::GetMinMax()
{
    auto bds = GetBounds();

    cfi::Position x1 = m_cfiEdge->calcXYZAtT(bds[0]);
    cfi::Position x2 = m_cfiEdge->calcXYZAtT(bds[1]);

    return {x1.x * m_scal, x1.y * m_scal, x1.z * m_scal,
            x2.x * m_scal, x2.y * m_scal, x2.z * m_scal};
}
} // namespace NekMesh
} // namespace Nektar
