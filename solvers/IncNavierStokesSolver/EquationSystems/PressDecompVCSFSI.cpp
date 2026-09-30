///////////////////////////////////////////////////////////////////////////////
//
// File: PressDecompVCSFSI.cpp
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
// Description: Velocity Correction Scheme for fluid-structure interaction with
// pressure decomposition
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/PressDecompVCSFSI.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

namespace Nektar
{
std::string PressDecompVCSFSI::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "PressDecompVCSFSI", PressDecompVCSFSI::create);

std::string PressDecompVCSFSI::solverTypeLookupId =
    LibUtilities::SessionReader::RegisterEnumValue(
        "SolverType", "PressDecompVCSFSI", ePressDecompVCSFSI);

/**
 * Constructor. Creates ...
 *
 * \param
 * \param
 */
PressDecompVCSFSI::PressDecompVCSFSI(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), VCSFSI(pSession, pGraph)
{
}

void PressDecompVCSFSI::v_InitObject(bool DeclareField)
{
    VelocityCorrectionScheme::v_InitObject(DeclareField);
    Array<OneD, NekDouble> tmp = m_movingFrameData + 18;
    m_rigidSolver.InitObject(m_session, m_fields[0], tmp);
    m_rigidSolver.SetMovableDoFs(m_movableDoFs);
    m_MRFABCname = "MRFWallPressDecomp";
}

/**
 * Destructor
 */
PressDecompVCSFSI::~PressDecompVCSFSI(void)
{
}

void PressDecompVCSFSI::v_DoInitialise(bool dumpInitialConditions)
{
    m_rigidSolver.SetInitialConditions(m_session, m_movingFrameData);
    VelocityCorrectionScheme::v_DoInitialise(dumpInitialConditions);
    std::set<int> dofs; // 0,1,2;3,4,5 six dofs
    GetMovableDoFs(dofs);
    SolvePotentials(dofs);
    OutputPotentials();
    m_rigidSolver.SetNewmarkBetaSolver(m_addedMass);
    Array<OneD, NekDouble> aeroforce(12, 0.);
    InitialiseFilter(aeroforce);
    m_rigidSolver.SetOldFvis(aeroforce);
}

void PressDecompVCSFSI::SolvePotentials(std::set<int> &dofs)
{
    std::map<int, Array<OneD, NekDouble>> boundValues;
    std::map<int, Array<OneD, NekDouble>> pPhys;
    CalculateBCs(dofs, boundValues);
    for (auto i : dofs)
    {
        pPhys[i] = Array<OneD, NekDouble>(m_pressure->GetTotPoints(), 0.);
        SolvePa(i, boundValues[i], pPhys[i]);
    }
    CalculateAddedMass(boundValues, pPhys);
}

void PressDecompVCSFSI::OutputPotentials()
{
    if (m_pCoef.size() == 0)
    {
        return;
    }
    std::string varname("phi0");
    std::vector<std::string> vars;
    std::vector<Array<OneD, NekDouble>> fields;
    for (auto it : m_pCoef)
    {
        varname[3] = '0' + it.first;
        vars.push_back(varname);
        fields.push_back(it.second);
    }
    WriteFld("potential.fld", m_pressure, fields, vars);
}

void PressDecompVCSFSI::CalculateBCs(std::set<int> &dofs,
                                     std::map<int, Array<OneD, NekDouble>> &bcs)
{
    int ndim = m_fields[0]->GetCoordim(0);
    // number of physics points for each field
    int numpts = 0;
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds =
        m_pressure->GetBndConditions();
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp =
        m_pressure->GetBndCondExpansions();
    for (int i = 0; i < BndConds.size(); ++i)
    {
        if (boost::iequals(BndConds[i]->GetUserDefined(), m_MRFABCname))
        {
            numpts += BndExp[i]->GetTotPoints();
        }
    }
    // load pivot point
    Array<OneD, NekDouble> pivot(3, 0.);
    if (m_movingFrameData.size() >= 21)
    {
        pivot[0] = m_movingFrameData[18];
        pivot[1] = m_movingFrameData[19];
        pivot[2] = m_movingFrameData[20];
    }
    // allocate boundary condition storage
    for (auto i : dofs)
    {
        bcs[i] = Array<OneD, NekDouble>(numpts, 0.);
    }
    // calculate values
    int offset = 0;
    for (int i = 0; i < BndConds.size(); ++i)
    {
        if (boost::iequals(BndConds[i]->GetUserDefined(), m_MRFABCname))
        {
            int npts = BndExp[i]->GetTotPoints();
            Array<OneD, Array<OneD, NekDouble>> n(3);
            Array<OneD, Array<OneD, NekDouble>> x(3);
            for (int j = 0; j < 3; ++j)
            {
                n[j] = Array<OneD, NekDouble>(npts, 0.);
                x[j] = Array<OneD, NekDouble>(npts, 0.);
            }
            BndExp[i]->GetNormals(n);
            BndExp[i]->GetCoords(x[0], x[1], x[2]);
            for (int j = 0; j < 3; ++j)
            {
                Vmath::Sadd(npts, -pivot[j], x[j], 1, x[j], 1);
            }
            Array<OneD, NekDouble> atmp;
            for (int j = 0; j < ndim; ++j)
            {
                if (dofs.find(j) != dofs.end())
                {
                    Vmath::Smul(npts, -1., n[j], 1, atmp = bcs[j] + offset, 1);
                }
            }
            if (dofs.find(5) != dofs.end())
            {
                atmp = bcs[5] + offset;
                Vmath::Vvtvvtm(npts, &x[1][0], 1, &n[0][0], 1, &x[0][0], 1,
                               &n[1][0], 1, &atmp[0], 1);
            }
            offset += npts;
        }
    }
}

void PressDecompVCSFSI::CalculateAddedMass(
    std::map<int, Array<OneD, NekDouble>> &boundValues,
    std::map<int, Array<OneD, NekDouble>> &pPhys)
{
    bool homo1D = false;
    m_session->MatchSolverInfo("Homogeneous", "1D", homo1D, false);
    NekDouble lz = 1.;
    if (homo1D)
    {
        lz = m_fields[0]->GetHomoLen();
    }
    // number of physics points for each field
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds =
        m_pressure->GetBndConditions();
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp =
        m_pressure->GetBndCondExpansions();

    int offset = 0;
    int nfld   = pPhys.size();
    Array<OneD, NekDouble> value(nfld * nfld,
                                 std::numeric_limits<NekDouble>::lowest());
    for (int i = 0; i < BndConds.size(); ++i)
    {
        if (boost::iequals(BndConds[i]->GetUserDefined(), m_MRFABCname))
        {
            MultiRegions::ExpListSharedPtr edgeExplist = BndExp[i];
            Array<OneD, Array<OneD, NekDouble>> phi(nfld);
            MultiRegions::ExpListSharedPtr BndElmtExp;
            m_pressure->GetBndElmtExpansion(i, BndElmtExp, false);
            Array<OneD, NekDouble> phiElm(BndElmtExp->GetTotPoints(), 0.);
            int j = 0;
            for (auto it : pPhys)
            {
                m_pressure->ExtractPhysToBndElmt(i, it.second, phiElm);
                m_pressure->ExtractElmtToBndPhys(i, phiElm, phi[j]);
                ++j;
            }
            ////phi_0 n_0, phi_0 n_1, phi_1 n_0, phi_1 n_1
            Array<OneD, NekDouble> mularray(edgeExplist->GetTotPoints(), 0.);
            for (int j = 0; j < nfld; ++j)
            {
                int k = 0;
                for (const auto &it : boundValues)
                {
                    Vmath::Vmul(edgeExplist->GetTotPoints(), phi[j], 1,
                                it.second + offset, 1, mularray, 1);
                    NekDouble temp = edgeExplist->Integral(mularray) / lz;
                    if (std::numeric_limits<NekDouble>::lowest() ==
                        value[k + j * nfld])
                    {
                        value[k + j * nfld] = temp;
                    }
                    else
                    {
                        value[k + j * nfld] += temp;
                    }
                    ++k;
                }
            }
            offset += edgeExplist->GetTotPoints();
        }
    }
    m_session->GetComm()->AllReduce(value, LibUtilities::ReduceMax);
    for (int j = 0; j < value.size(); ++j)
    {
        if (std::numeric_limits<NekDouble>::lowest() == value[j])
        {
            value[j] = 0;
        }
    }
    int NumDofs = m_spacedim + 1;
    m_addedMass = Array<OneD, NekDouble>(NumDofs * NumDofs, 0.);
    int i       = 0;
    for (auto it : pPhys)
    {
        int i1 = std::min(it.first, m_spacedim);
        int j  = 0;
        for (auto jt : pPhys)
        {
            int j1                         = std::min(jt.first, m_spacedim);
            m_addedMass[i1 + j1 * NumDofs] = value[i + nfld * j];
            ++j;
        }
        ++i;
    }
    if (m_session->GetComm()->GetRank() == 0)
    {
        for (int j = 0; j < NumDofs; ++j)
        {
            for (int k = 0; k < NumDofs; ++k)
            {
                std::cout << "value[" << j << ", " << k
                          << "] = " << std::scientific << std::setprecision(7)
                          << m_addedMass[k + j * NumDofs] << std::endl;
            }
        }
    }
}

void PressDecompVCSFSI::SolvePa(int i, Array<OneD, NekDouble> bc,
                                Array<OneD, NekDouble> pPhys)
{
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds =
        m_pressure->GetBndConditions();
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp =
        m_pressure->GetBndCondExpansions();
    int offset = 0;
    for (int i = 0; i < BndConds.size(); ++i)
    {
        if (boost::iequals(BndConds[i]->GetUserDefined(), m_MRFABCname))
        {
            BndExp[i]->IProductWRTBase(bc + offset, BndExp[i]->UpdateCoeffs());
            offset += BndExp[i]->GetTotPoints();
        }
    }
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;
    Array<OneD, NekDouble> forcing(m_pressure->GetTotPoints(), 0.);
    m_pCoef[i] = Array<OneD, NekDouble>(m_pressure->GetNcoeffs(), 0.);
    m_pressure->HelmSolve(forcing, m_pCoef[i], factors);
    m_pressure->BwdTrans(m_pCoef[i], pPhys);
}

/**
 * Correct pressure by adding potential terms
 */
void PressDecompVCSFSI::CorrectPressure()
{
    for (const auto &it : m_pCoef)
    {
        NekDouble acceleration = m_movingFrameData[it.first + 12];
        if (std::fabs(acceleration) != 0.)
        {
            Vmath::Svtvp(m_pressure->GetNcoeffs(), acceleration, it.second, 1,
                         m_pressure->UpdateCoeffs(), 1,
                         m_pressure->UpdateCoeffs(), 1);
        }
    }
}

void PressDecompVCSFSI::v_SolveSolid(NekDouble time)
{
    VCSFSI::v_SolveSolid(time);
    CorrectPressure();
}

} // namespace Nektar
