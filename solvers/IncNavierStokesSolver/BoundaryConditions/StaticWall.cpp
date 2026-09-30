///////////////////////////////////////////////////////////////////////////////
//
// File: StaticWall.cpp
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
// Description: Stationary wall boundary condition.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/StaticWall.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{

std::string StaticWall::className = GetIncBCFactory().RegisterCreatorFunction(
    "StaticWall", StaticWall::create, "Stationary wall boundary condition");

StaticWall::StaticWall(
    [[maybe_unused]] const LibUtilities::SessionReaderSharedPtr pSession,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    [[maybe_unused]] Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> exp,
    [[maybe_unused]] int nbnd, [[maybe_unused]] int spacedim,
    [[maybe_unused]] int bnddim)
    : IncBaseCondition(pSession, pFields, cond, exp, nbnd, spacedim, bnddim)
{
    classname  = "StaticWall";
    m_pressure = cond.size() - 1;
    if (cond[m_pressure]->GetUserDefined() == classname)
    {
        m_BndExp[m_pressure] = exp[m_pressure];
    }
}

void StaticWall::v_Initialise(
    const LibUtilities::SessionReaderSharedPtr &pSession)
{
    if (m_BndExp.empty())
    {
        return;
    }
    IncBaseCondition::v_Initialise(pSession);
    m_field->GetBndElmtExpansion(m_nbnd, m_bndElmtExps, false);
    m_bndElmtExps->SetWaveSpace(m_field->GetWaveSpace());
    m_extrapArray =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_intSteps);
    for (int n = 0; n < m_intSteps; ++n)
    {
        m_extrapArray[n] = Array<OneD, Array<OneD, NekDouble>>(m_bnddim);
        for (int i = 0; i < m_bnddim; ++i)
        {
            m_extrapArray[n][i] = Array<OneD, NekDouble>(m_npoints, 0.0);
        }
    }
}

void StaticWall::v_Update(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &fields,
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Adv,
    std::map<std::string, NekDouble> &params)
{
    if (m_BndExp.empty() || fields.size() == 0 ||
        params.find("pressure") == params.end())
    {
        return;
    }
    ++m_numCalls;
    // pressure
    Array<OneD, Array<OneD, NekDouble>> rhs(m_bnddim);
    for (int i = 0; i < m_bnddim; ++i)
    {
        rhs[i] = Array<OneD, NekDouble>(m_npoints, 0.);
    }
    AddExtrapVisPressureBCs(fields, rhs, params);
    m_BndExp[m_pressure]->NormVectorIProductWRTBase(
        rhs, m_BndExp[m_pressure]->UpdateCoeffs());
}

} // namespace Nektar
