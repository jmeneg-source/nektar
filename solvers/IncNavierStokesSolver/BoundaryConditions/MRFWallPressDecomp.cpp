///////////////////////////////////////////////////////////////////////////////
//
// File: MRFWallPressDecomp.cpp
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
// Description: Wall boundary condition of moving reference frame with pressure
// decomposition.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/MRFWallPressDecomp.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
std::string MRFWallPressDecomp::className =
    GetIncBCFactory().RegisterCreatorFunction(
        "MRFWallPressDecomp", MRFWallPressDecomp::create,
        "Wall boundary condition of moving reference frame with pressure "
        "decomposition");

MRFWallPressDecomp::MRFWallPressDecomp(
    [[maybe_unused]] const LibUtilities::SessionReaderSharedPtr pSession,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    [[maybe_unused]] Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> exp,
    [[maybe_unused]] int nbnd, [[maybe_unused]] int spacedim,
    [[maybe_unused]] int bnddim)
    : StaticWall(pSession, pFields, cond, exp, nbnd, spacedim, bnddim)
{
    classname = "MRFWallPressDecomp";
    m_hasVels = false;
    for (size_t i = 0; i < m_bnddim; ++i)
    {
        if (cond[i]->GetUserDefined() == classname)
        {
            m_hasVels   = true;
            m_BndExp[i] = exp[i];
        }
    }
    m_pressure    = cond.size() - 1;
    m_hasPressure = cond[m_pressure]->GetUserDefined() == classname;
    if (m_hasPressure)
    {
        m_BndExp[m_pressure] = exp[m_pressure];
    }
}

MRFWallPressDecomp::~MRFWallPressDecomp()
{
}

void MRFWallPressDecomp::v_Initialise(
    const LibUtilities::SessionReaderSharedPtr &pSession)
{
    IncBaseCondition::v_Initialise(pSession);
    if (!(pSession->GetSolverInfo("SolverType") == "PressDecompVCSFSI"))
    {
        ASSERTL0(false,
                 "The boundary condition MRFWallPressDecomp is only "
                 "supported for use in the SolverType 'PressDecompVCSFSI'");
    }
    m_field->GetBndElmtExpansion(m_nbnd, m_bndElmtExps, false);
    if (m_hasPressure)
    {
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
}

/// @brief v_Update set correct BCs (in wavespace)
/// @param fields
/// @param Adv is in wavespace for 3DH1D
/// @param params
void MRFWallPressDecomp::v_Update(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &fields,
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Adv,
    std::map<std::string, NekDouble> &params)
{
    int nptsPlane0 = 0;
    SetNumPointsOnPlane0(nptsPlane0);
    if (params.find("Omega_x") != params.end() ||
        params.find("Omega_y") != params.end() ||
        params.find("Omega_z") != params.end())
    {
        InitialiseCoords(params);
    }
    // pressure, do if not define velocity
    if (params.find("pressure") != params.end() && m_hasPressure &&
        fields.size() > 0)
    {
        ++m_numCalls;

        Array<OneD, Array<OneD, NekDouble>> rhs(m_bnddim);
        for (int i = 0; i < m_bnddim; ++i)
        {
            rhs[i] = Array<OneD, NekDouble>(m_npoints, 0.);
        }
        // add viscous term and centripetal acceleration
        AddExtrapCentVisPressureBCs(fields, rhs, params, nptsPlane0);
        m_BndExp[m_pressure]->NormVectorIProductWRTBase(
            rhs, m_BndExp[m_pressure]->UpdateCoeffs());
    }
    // velocity, do if define velocity
    if (params.find("velocity") != params.end() && m_hasVels && nptsPlane0)
    {
        Array<OneD, Array<OneD, NekDouble>> velocities(m_bnddim);
        for (size_t k = 0; k < m_bnddim; ++k)
        {
            if (m_BndExp.find(k) != m_BndExp.end())
            {
                velocities[k] = Array<OneD, NekDouble>(nptsPlane0, 0.0);
            }
        }
        RigidBodyVelocity(velocities, params, nptsPlane0);
        for (int k = 0; k < m_bnddim; ++k)
        {
            if (m_BndExp.find(k) != m_BndExp.end())
            {
                if (m_BndExp[k]->GetExpType() == MultiRegions::e2DH1D)
                {
                    m_BndExp[k]->GetPlane(0)->FwdTransBndConstrained(
                        velocities[k],
                        m_BndExp[k]->GetPlane(0)->UpdateCoeffs());
                }
                else
                {
                    m_BndExp[k]->FwdTransBndConstrained(
                        velocities[k], m_BndExp[k]->UpdateCoeffs());
                }
            }
        }
    }
}

void MRFWallPressDecomp::AddCentripetalAcc(
    Array<OneD, Array<OneD, NekDouble>> &N,
    std::map<std::string, NekDouble> &params, int npts0)
{
    if (npts0 == 0 || params.find("Omega_z") == params.end())
    {
        return;
    }
    // add centripetal acceleration
    NekDouble Omega = params["Omega_z"];
    NekDouble Wz2   = Omega * Omega;
    Vmath::Svtvp(npts0, Wz2, m_coords[0], 1, N[0], 1, N[0], 1);
    Vmath::Svtvp(npts0, Wz2, m_coords[1], 1, N[1], 1, N[1], 1);
    if (params.find("U") != params.end())
    {
        NekDouble mOmegaU0 = -Omega * params["U"];
        Vmath::Sadd(npts0, mOmegaU0, N[1], 1, N[1], 1);
    }
    if (params.find("V") != params.end())
    {
        NekDouble OmegaV0 = Omega * params["V"];
        Vmath::Sadd(npts0, OmegaV0, N[0], 1, N[0], 1);
    }
}

void MRFWallPressDecomp::AddExtrapCentVisPressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    Array<OneD, Array<OneD, NekDouble>> &N,
    std::map<std::string, NekDouble> &params, int npts0)
{
    for (int i = 0; i < m_bnddim; ++i)
    {
        Vmath::Zero(m_npoints, m_extrapArray[m_intSteps - 1][i], 1);
    }
    AddVisPressureBCs(fields, m_extrapArray[m_intSteps - 1], params);
    AddCentripetalAcc(m_extrapArray[m_intSteps - 1], params, npts0);
    ExtrapolateArray(m_numCalls, m_extrapArray);
    for (int i = 0; i < m_bnddim; i++)
    {
        Vmath::Vadd(m_npoints, m_extrapArray[m_intSteps - 1][i], 1, N[i], 1,
                    N[i], 1);
    }
}

} // namespace Nektar
