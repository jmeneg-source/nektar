///////////////////////////////////////////////////////////////////////////////
//
// File: IncBoundaryConditions.cpp
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
// Description: Base class for boundary conditions.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{

std::set<std::string> IncBoundaryConditions::m_BndType = {
    "MRFFar", "StaticWall", "TransMovingWall", "MRFWall", "MRFWallPressDecomp"};

IncBoundaryConditions::IncBoundaryConditions()
{
}

/**
 *  Initialize the partial slip boundary conditions
 *  The total points, unit normal vector, and the slip length are assigned
 */
void IncBoundaryConditions::Initialize(
    const LibUtilities::SessionReaderSharedPtr pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    switch (pFields[0]->GetExpType())
    {
        case MultiRegions::e2D:
        {
            m_spacedim = 2;
            m_bnd_dim  = 2;
        }
        break;
        case MultiRegions::e3DH1D:
        {
            m_spacedim = 3;
            m_bnd_dim  = 2;
        }
        break;
        case MultiRegions::e3DH2D:
        {
            m_spacedim = 3;
            m_bnd_dim  = 1;
        }
        break;
        case MultiRegions::e3D:
        {
            m_spacedim = 3;
            m_bnd_dim  = 3;
        }
        break;
        default:
            ASSERTL0(0, "Dimension not supported");
            break;
    }

    Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr>>
        BndConds(m_spacedim + 1);
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr>> BndExp(m_spacedim +
                                                                    1);

    for (int i = 0; i < m_spacedim; ++i)
    {
        BndConds[i] = pFields[i]->GetBndConditions();
        BndExp[i]   = pFields[i]->GetBndCondExpansions();
    }
    int npress           = pFields.size() - 1;
    BndConds[m_spacedim] = pFields[npress]->GetBndConditions();
    BndExp[m_spacedim]   = pFields[npress]->GetBndCondExpansions();

    for (size_t n = 0; n < BndExp[0].size(); ++n)
    {
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond(
            BndConds.size());
        Array<OneD, MultiRegions::ExpListSharedPtr> exp(BndConds.size());
        std::string bndtype;
        for (int k = 0; k < BndConds.size(); ++k)
        {
            cond[k] = BndConds[k][n];
            exp[k]  = BndExp[k][n];
            if (bndtype.size() == 0 &&
                m_BndType.find(cond[k]->GetUserDefined()) != m_BndType.end())
            {
                bndtype = cond[k]->GetUserDefined();
            }
        }
        if (bndtype.size())
        {
            m_bounds[n] = GetIncBCFactory().CreateInstance(
                bndtype, pSession, pFields, cond, exp, n, m_spacedim,
                m_bnd_dim);
        }
    }
}

void IncBoundaryConditions::Update(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &Adv,
    std::map<std::string, NekDouble> &params)
{
    for (auto &it : m_bounds)
    {
        it.second->Update(fields, Adv, params);
    }
}

} // namespace Nektar
