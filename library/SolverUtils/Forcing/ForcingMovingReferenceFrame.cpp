///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.cpp
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
// Description: Solving the absolute flow in a moving body frame,
// by adding (U0 + Omega X (x - x0)) . grad u - Omega X u
// as the body force.
// U0 is the translational velocity of the body frame.
// Omega is the angular velocity.
// x0 is the rotation pivot in the body frame.
// All vectors use the basis of the body frame.
// Translational motion is allowed for all dimensions.
// Rotation is not allowed for 1D, 2DH1D, 3DH2D.
// Rotation in z direction is allowed for 2D and 3DH1D.
// Rotation in 3 directions are allowed for 3D.
// TODO: add suport for 3D rotation using Quaternion
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <SolverUtils/Forcing/ForcingMovingReferenceFrame.h>

namespace Nektar::SolverUtils
{

std::string ForcingMovingReferenceFrame::classNameBody =
    GetForcingFactory().RegisterCreatorFunction(
        "MovingReferenceFrame", ForcingMovingReferenceFrame::create,
        "Moving Frame");

/**
 * @brief
 * @param pSession
 * @param pEquation
 */
ForcingMovingReferenceFrame::ForcingMovingReferenceFrame(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

ForcingMovingReferenceFrame::~ForcingMovingReferenceFrame(void)
{
}

/**
 * @brief Initialise the forcing module
 * @param pFields
 * @param pNumForcingFields
 * @param pForce
 */
void ForcingMovingReferenceFrame::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const unsigned int &pNumForcingFields,
    [[maybe_unused]] const TiXmlElement *pForce)
{
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isH1d, false);
    m_session->MatchSolverInfo("Homogeneous", "2D", m_isH2d, false);
    bool singleMode, halfMode;
    m_session->MatchSolverInfo("ModeType", "SingleMode", singleMode, false);
    m_session->MatchSolverInfo("ModeType", "HalfMode", halfMode, false);
    if (singleMode || halfMode)
    {
        m_isH1d = false;
    }
    int expdim  = m_isH2d ? 1 : pFields[0]->GetGraph()->GetMeshDimension();
    m_spacedim  = expdim + (m_isH1d ? 1 : 0) + (m_isH2d ? 2 : 0);
    m_hasPlane0 = true;
    if (m_isH1d)
    {
        m_hasPlane0 = pFields[0]->GetZIDs()[0] == 0;
    }
}

/**
 * @brief Adds the body force, -Omega X u.
 * @param fields
 * @param inarray
 * @param outarray
 * @param time
 */
void ForcingMovingReferenceFrame::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble &time)
{
    // If there is no rotation, body force is zero,
    // nothing needs to be done here.
    if (m_hasRotation)
    {
        int npoints = fields[0]->GetNpoints();
        addRotation(npoints, outarray, -1., inarray, outarray);
    }
}

/**
 * @brief outarray = inarray0 + angVelScale Omega x inarray1
 */
void ForcingMovingReferenceFrame::addRotation(
    int nPnts, // number of points
    const Array<OneD, Array<OneD, NekDouble>> &inarray0, NekDouble angVelScale,
    const Array<OneD, Array<OneD, NekDouble>> &inarray1,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    ASSERTL0(&inarray1 != &outarray, "inarray1 and outarray "
                                     "should not be the same.");

    // TODO: In case of having support for all three components of Omega,
    // they should be transformed into the rotating frame first!

    // In case that the inarray0 and outarry are different, to avoid using
    // un-initialized array, copy the array first
    if (&inarray0 != &outarray)
    {
        ASSERTL0(inarray0.size() == outarray.size(),
                 "inarray0 and outarray must have same dimentions");
        for (int i = 0; i < inarray0.size(); ++i)
        {
            Vmath::Vcopy(nPnts, inarray0[i], 1, outarray[i], 1);
        }
    }

    if (m_spacedim >= 2 && m_hasOmega[2])
    {
        NekDouble cp = m_omegaxyz[2] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cm, inarray1[1], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cp, inarray1[0], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[0])
    {
        NekDouble cp = m_omegaxyz[0] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[1], 1, outarray[2], 1, outarray[2], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[2], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[1])
    {
        NekDouble cp = m_omegaxyz[1] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[2], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[0], 1, outarray[2], 1, outarray[2], 1);
    }
}

/**
 * @brief Compute the moving frame velocity at given time
 */
void ForcingMovingReferenceFrame::v_PreApply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble &time)
{
    UpdateMRFStatus(fields[0]);
    int npoints = fields[0]->GetNpoints();
    if (m_isH2d && fields[0]->GetWaveSpace())
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Array<OneD, NekDouble> tmpphys(npoints, -m_velxyz[i]);
                Array<OneD, NekDouble> tmpcoef(npoints);
                fields[0]->HomogeneousFwdTrans(npoints, tmpphys, tmpcoef);
                Vmath::Vadd(npoints, tmpcoef, 1, inarray[i], 1, outarray[i], 1);
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
    }
    else
    {
        int npoints0 = npoints;
        if (m_isH1d && fields[0]->GetWaveSpace())
        {
            npoints0 = m_hasPlane0 ? fields[0]->GetPlane(0)->GetNpoints() : 0;
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Vmath::Sadd(npoints0, -m_velxyz[i], inarray[i], 1, outarray[i],
                            1);
                if (&inarray != &outarray && npoints != npoints0)
                {
                    Array<OneD, NekDouble> tmp = outarray[i] + npoints0;
                    Vmath::Vcopy(npoints - npoints0, inarray[i] + npoints0, 1,
                                 tmp, 1);
                }
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
        if (m_hasRotation)
        {
            addRotation(npoints0, outarray, -1., m_coords, outarray);
        }
    }
}

void ForcingMovingReferenceFrame::UpdateMRFStatus(
    MultiRegions::ExpListSharedPtr field)
{
    if (m_velxyz.size() == 0)
    {
        // initialize variables
        m_velxyz      = Array<OneD, NekDouble>(3, 0.0);
        m_omegaxyz    = Array<OneD, NekDouble>(3, 0.0);
        m_hasVel      = Array<OneD, bool>(3, false);
        m_hasOmega    = Array<OneD, bool>(3, false);
        m_hasRotation = false;
        m_pivotPoint  = Array<OneD, NekDouble>(3, 0.0);

        // initialise pivot point for fluid interface
        auto equ = m_equ.lock();
        ASSERTL0(equ, "Weak pointer to the equation system is expired");
        m_FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
        m_FluidEq->GetMovingFramePivot(m_pivotPoint);
        std::set<int> moveDoFs;
        m_FluidEq->GetMovableDoFs(moveDoFs);
        for (int i = 0; i < 3; ++i)
        {
            if (moveDoFs.find(i) != moveDoFs.end())
            {
                m_hasVel[i] = true;
            }
            if (moveDoFs.find(i + 3) != moveDoFs.end())
            {
                m_hasOmega[i] = true;
                m_hasRotation = true;
            }
        }
        // account for the effect of rotation
        // Omega_X results in having v and w even if not defined by user
        // Omega_Y results in having u and w even if not defined by user
        // Omega_Z results in having u and v even if not defined by user
        for (int i = 0; i < 3; ++i)
        {
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            if (m_hasOmega[i])
            {
                m_hasVel[j] = true;
                m_hasVel[k] = true;
            }
        }
        if (m_hasRotation)
        {
            int npoints = field->GetNpoints();
            m_coords    = Array<OneD, Array<OneD, NekDouble>>(3);
            for (int j = 0; j < m_spacedim; ++j)
            {
                m_coords[j] = Array<OneD, NekDouble>(npoints);
            }
            field->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
            // move the origin to the pivot point
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Sadd(npoints, -m_pivotPoint[i], m_coords[i], 1,
                            m_coords[i], 1);
            }
        }
    }
    Array<OneD, NekDouble> vel(6, 0.0);
    m_FluidEq->GetMovingFrameVelocities(vel);
    Vmath::Vcopy(3, vel, 1, m_velxyz, 1);
    Vmath::Vcopy(3, vel + 3, 1, m_omegaxyz, 1);
}

} // namespace Nektar::SolverUtils
