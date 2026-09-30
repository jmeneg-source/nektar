///////////////////////////////////////////////////////////////////////////////
//
// File: RigidSolver.cpp
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

#include <IncNavierStokesSolver/EquationSystems/RigidSolver.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <boost/format.hpp>

namespace Nektar
{

void RigidSolver::InitObject(const LibUtilities::SessionReaderSharedPtr session,
                             const MultiRegions::ExpListSharedPtr &pField,
                             Array<OneD, NekDouble> pivot)
{
    bool isH1d, isH2d;
    session->MatchSolverInfo("Homogeneous", "1D", isH1d, false);
    session->MatchSolverInfo("Homogeneous", "2D", isH2d, false);
    bool singleMode, halfMode;
    session->MatchSolverInfo("ModeType", "SingleMode", singleMode, false);
    session->MatchSolverInfo("ModeType", "HalfMode", halfMode, false);
    if (singleMode || halfMode)
    {
        isH1d = false;
    }
    int expdim            = isH2d ? 1 : pField->GetGraph()->GetMeshDimension();
    m_spacedim            = expdim + (isH1d ? 1 : 0) + (isH2d ? 2 : 0);
    m_isRoot              = pField->GetComm()->TreatAsRankZero();
    m_index               = 1;
    m_currentTime         = -1.;
    m_extForceXYZ         = Array<OneD, NekDouble>(6, 0.0);
    m_oldFvis             = Array<OneD, NekDouble>(6, 0.0);
    TiXmlElement *pSolver = session->GetElement("Nektar/RIGIDSOLVER");
    LoadParameters(session, pSolver);
    InitBodySolver(session, pSolver, pivot);
    CheckParameters();
}

void RigidSolver::CheckParameters()
{
    // Count free translational DoFs
    int nFreeTrans = 0;
    bool freeX     = false;

    for (int i = 0; i < m_spacedim; ++i)
    {
        if (m_dirDoFs.find(i) == m_dirDoFs.end())
        {
            ++nFreeTrans;
            if (i == 0)
            {
                freeX = true;
            }
        }
    }
    if (!m_hasRotation)
    {
        m_solveType = 0; // inertial frame
    }
    else if (nFreeTrans == 1 && freeX)
    {
        m_solveType = 1; // with given rotation and only x free
        if (m_spacedim >= 2)
        {
            m_dirDoFs.erase(1);
        }
    }
    else if ((nFreeTrans == 2 || nFreeTrans == 0) && m_hasRotation)
    {
        m_solveType = 2; // full planar rigid solver (rotation involved)
    }
    else
    {
        ASSERTL0(false, "Unsupported rigid motion: this solver supports only "
                        "(i) translation-only motion, "
                        "(ii) x-only translation with optional rotation, or "
                        "(iii) full planar rigid-body motion.");
    }
    if (m_hasRotation)
    {
        m_inertialPosition = Array<OneD, NekDouble>(2, 0.);
        // m_K should be zero
        NekDouble sum = 0;
        for (size_t i = 0; i < m_K.size(); ++i)
        {
            sum += fabs(m_K[i]);
        }
        for (size_t i = 0; i < m_C.size(); ++i)
        {
            sum += fabs(m_C[i]);
        }
        ASSERTL0(sum < NekConstants::kNekZeroTol,
                 "C, K should be zero [in body frame solver].");
    }
}

void RigidSolver::SetMovableDoFs(std::vector<bool> &moveDoFs)
{
    for (const auto &it : m_frameVelFunction)
    {
        if (it.first < 6)
        {
            moveDoFs[it.first] = true;
        }
    }
    for (int i = 0; i < m_spacedim; ++i)
    {
        if (m_dirDoFs.find(i) == m_dirDoFs.end())
        {
            moveDoFs[i] = true;
        }
    }
    if (m_dirDoFs.find(m_spacedim) == m_dirDoFs.end())
    {
        moveDoFs[5] = true;
    }
    if (moveDoFs[5])
    {
        // since m_frameVelFunction is defined in inerital frame
        // motions in 2 directions occur if there is rotation
        if (moveDoFs[0] || moveDoFs[1])
        {
            moveDoFs[0] = true;
            moveDoFs[1] = true;
        }
    }
}

NekDouble RigidSolver::EvaluateExpression(
    const LibUtilities::SessionReaderSharedPtr session, std::string expression)
{
    NekDouble value = 0.;
    try
    {
        LibUtilities::Equation expession(session->GetInterpreter(), expression);
        value = expession.Evaluate();
    }
    catch (const std::runtime_error &)
    {
        NEKERROR(ErrorUtil::efatal, "Error evaluating expression" + expression);
    }
    return value;
}

void RigidSolver::ParserFunctionToMap(
    const bool allowzero, const LibUtilities::SessionReaderSharedPtr session,
    const std::string &FuncName, const std::string &var,
    std::map<int, LibUtilities::EquationSharedPtr> &result, const int i)
{
    if (session->DefinesFunction(FuncName, var))
    {
        LibUtilities::EquationSharedPtr equ =
            session->GetFunction(FuncName, var);
        if (allowzero || "0" != equ->GetExpression())
        {
            result[i] = equ;
        }
    }
}

void RigidSolver::LoadParameters(
    const LibUtilities::SessionReaderSharedPtr session,
    const TiXmlElement *pSolver)
{
    const TiXmlElement *funcNameElmt;
    // load frame velocity
    funcNameElmt = pSolver->FirstChildElement("FRAMEVELOCITY");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        ASSERTL0(session->DefinesFunction(FuncName),
                 "Function '" + FuncName + "' is not defined in the session.");
        // linear velocity
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string var = session->GetVariable(i);
            ParserFunctionToMap(false, session, FuncName, var,
                                m_frameVelFunction, i);
        }
        // linear displacement and acceleration
        std::vector<std::string> linearDispVar = {"X", "Y", "Z"};
        std::vector<std::string> linearAcceVar = {"A_x", "A_y", "A_z"};
        for (int i = 0; i < m_spacedim; ++i)
        {
            ParserFunctionToMap(true, session, FuncName, linearDispVar[i],
                                m_frameVelFunction, i + 6);
            ParserFunctionToMap(true, session, FuncName, linearAcceVar[i],
                                m_frameVelFunction, i + 12);
        }
        // angular velocities
        std::vector<std::string> angularVar = {"Omega_x", "Omega_y", "Omega_z"};
        m_hasRotation                       = false;
        for (int i = 0; i < 3; ++i)
        {
            std::string var = angularVar[i];
            ParserFunctionToMap(false, session, FuncName, var,
                                m_frameVelFunction, i + 3);
            m_hasRotation = m_hasRotation || m_frameVelFunction.find(i + 3) !=
                                                 m_frameVelFunction.end();
        }
        // angular displacement and acceleration
        std::vector<std::string> angularDispVar  = {"Theta_x", "Theta_y",
                                                    "Theta_z"};
        std::vector<std::string> angularAccepVar = {"DOmega_x", "DOmega_y",
                                                    "DOmega_z"};
        for (int i = 0; i < 3; ++i)
        {
            ParserFunctionToMap(true, session, FuncName, angularDispVar[i],
                                m_frameVelFunction, i + 3 + 6);
            ParserFunctionToMap(true, session, FuncName, angularAccepVar[i],
                                m_frameVelFunction, i + 3 + 12);
        }
    }

    // load external force
    funcNameElmt = pSolver->FirstChildElement("EXTERNALFORCE");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        if (session->DefinesFunction(FuncName))
        {
            std::vector<std::string> forceVar = {"Fx", "Fy", "Fz"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                std::string var = forceVar[i];
                ParserFunctionToMap(false, session, FuncName, var,
                                    m_extForceFunction, i);
            }
            std::vector<std::string> momentVar = {"Mx", "My", "Mz"};
            for (int i = 0; i < 3; ++i)
            {
                std::string var = momentVar[i];
                ParserFunctionToMap(false, session, FuncName, var,
                                    m_extForceFunction, i + 3);
            }
        }
    }

    // OutputFile
    const TiXmlElement *mssgTag = pSolver->FirstChildElement("OutputFile");
    std::string filename;
    if (mssgTag)
    {
        filename = mssgTag->GetText();
    }
    else
    {
        filename = session->GetSessionName();
    }
    if (!(filename.length() >= 4 &&
          filename.substr(filename.length() - 4) == ".mrf"))
    {
        filename += ".mrf";
    }
    if (m_isRoot)
    {
        m_outputStream.open(filename.c_str());
        if (m_spacedim == 2)
        {
            m_outputStream
                << "Variables = t, x, ux, ax, y, uy, ay, theta, omega, domega"
                << std::endl;
        }
        else if (m_spacedim == 3)
        {
            m_outputStream << "Variables = t, x, ux, ax, y, uy, ay, z, uz, az, "
                              "theta, omega, domega"
                           << std::endl;
        }
    }

    // output frequency
    m_outputFrequency = 1;
    mssgTag           = pSolver->FirstChildElement("OutputFrequency");
    if (mssgTag)
    {
        std::vector<std::string> values;
        std::string mssgStr = mssgTag->GetText();
        m_outputFrequency   = round(EvaluateExpression(session, mssgStr));
    }
    ASSERTL0(m_outputFrequency > 0,
             "OutputFrequency should be greater than zero.");
}

void RigidSolver::InitBodySolver(
    const LibUtilities::SessionReaderSharedPtr session,
    const TiXmlElement *pSolver, Array<OneD, NekDouble> pivot)
{
    int NumDof = m_spacedim + 1;
    const TiXmlElement *mssgTag;
    std::string mssgStr;
    // allocate memory and initialise
    m_vel = Array<OneD, Array<OneD, NekDouble>>(3);
    for (size_t i = 0; i < 3; ++i)
    {
        m_vel[i] = Array<OneD, NekDouble>(NumDof, 0.);
    }
    // read prescribed motion DoFs
    m_dirDoFs.clear();
    for (int i = 0; i < NumDof; ++i)
    {
        m_dirDoFs.insert(i);
    }
    mssgTag = pSolver->FirstChildElement("MOTIONPRESCRIBED");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof,
                 "MOTIONPRESCRIBED vector should be of size " +
                     std::to_string(NumDof));
        for (int i = 0; i < NumDof; ++i)
        {
            if (EvaluateExpression(session, values[i]) == 0)
            {
                m_dirDoFs.erase(i);
            }
        }
    }
    m_hasRotation =
        m_hasRotation || m_dirDoFs.find(m_spacedim) == m_dirDoFs.end();
    m_hasFreeMotion = m_dirDoFs.size() < NumDof;
    // read mass matrix
    m_M     = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag = pSolver->FirstChildElement("MASS");
    ASSERTL0(m_dirDoFs.size() == NumDof || mssgTag, "Mass matrix is required.");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == 1, "Mass matrix should be a scalar.");
        m_mass = EvaluateExpression(session, values[0]);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_M[i + i * NumDof] = m_mass;
        }
    }
    mssgTag = pSolver->FirstChildElement("ROTATIONINERTIA");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == 1, "Inertia should be a scalar.");
        m_rotaionInertia         = EvaluateExpression(session, values[0]);
        m_M[NumDof * NumDof - 1] = m_rotaionInertia;
    }
    // read damping matrix
    m_C     = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag = pSolver->FirstChildElement("DAMPING");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof * NumDof,
                 "Damping matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        int count = 0;
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                m_C[count] = EvaluateExpression(session, values[count]);
                ++count;
            }
        }
    }
    // read rigidity matrix
    m_K     = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag = pSolver->FirstChildElement("RIGIDITY");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof * NumDof,
                 "Rigidity matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        int count = 0;
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                m_K[count] = EvaluateExpression(session, values[count]);
                ++count;
            }
        }
    }
    // read pivot point
    mssgTag = pSolver->FirstChildElement("PIVOTPOINT");
    m_pivot = Array<OneD, NekDouble>(m_spacedim, 0.);
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == m_spacedim,
                 "PIVOTPOINT vector should be of size " +
                     std::to_string(m_spacedim));
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_pivot[i] = EvaluateExpression(session, values[i]);
        }
    }
    Vmath::Vcopy(m_spacedim, m_pivot, 1, pivot, 1);
    // read the distance between pivotpoint and masscenter
    mssgTag         = pSolver->FirstChildElement("PIVOTDISTANCE");
    m_pivotdistance = 0.;
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == 1, "PivotDistance should be a scalar.");
        m_pivotdistance = EvaluateExpression(session, values[0]);
    }
    // read Newmark Beta paramters
    m_timestep = session->GetParameter("TimeStep");
    m_beta     = 0.25;
    m_gamma    = 0.51;
    if (session->DefinesParameter("NewmarkBeta"))
    {
        m_beta = session->GetParameter("NewmarkBeta");
    }
    if (session->DefinesParameter("NewmarkGamma"))
    {
        m_gamma = session->GetParameter("NewmarkGamma");
    }
}

void RigidSolver::UpdatePrescribed(const NekDouble &time,
                                   std::map<int, NekDouble> &Dirs)
{
    int NumDof = m_spacedim + 1;
    for (auto it : m_frameVelFunction)
    {
        if (it.first < 3)
        {
            Dirs[it.first] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first == 5)
        {
            Dirs[m_spacedim] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first < 9)
        {
            Dirs[NumDof + it.first - 6] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first == 11)
        {
            Dirs[NumDof + m_spacedim] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first < 15)
        {
            Dirs[(NumDof << 1) + it.first - 12] =
                it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first == 17)
        {
            Dirs[(NumDof << 1) + m_spacedim] =
                it.second->Evaluate(0., 0., 0., time);
        }
    }
    for (auto i : m_dirDoFs)
    {
        if (Dirs.find(i) == Dirs.end())
        {
            Dirs[i]                 = 0.;
            Dirs[i + NumDof]        = 0.;
            Dirs[i + (NumDof << 1)] = 0.;
        }
    }
}

void RigidSolver::SetOldFvis(Array<OneD, NekDouble> force)
{
    for (int i = 0; i < 6; ++i)
    {
        m_oldFvis[i] = force[6 + i];
    }
}

/**
 * @brief Updates the forcing array with the current required forcing.
 * @param pFields
 * @param time
 */
void RigidSolver::UpdateFrameVelocity(Array<OneD, NekDouble> &aeroforce,
                                      const NekDouble &time,
                                      Array<OneD, NekDouble> &MRFData)
{
    if (m_currentTime >= time)
    {
        return;
    }

    for (int i = 0; i < 6; ++i)
    {
        aeroforce[i] +=
            2. * aeroforce[6 + i] - m_oldFvis[i]; // extrapolate viscous
        m_oldFvis[i] = aeroforce[6 + i];
    }
    m_currentTime = time;
    std::map<int, NekDouble> Dirs;
    UpdatePrescribed(time, Dirs);
    // compute the velocites whoes functions are provided in inertial frame
    if (m_hasFreeMotion)
    {
        for (auto it : m_extForceFunction)
        {
            m_extForceXYZ[it.first] = it.second->Evaluate(0., 0., 0., time);
        }
    }
    SolveBodyMotion(m_vel, aeroforce, Dirs);
    Array<OneD, Array<OneD, NekDouble>> tmpVel;
    if (m_hasRotation)
    {
        // update inertial position of the body
        tmpVel = Array<OneD, Array<OneD, NekDouble>>(m_vel.size());
        for (size_t i = 0; i < m_vel.size(); ++i)
        {
            tmpVel[i] = Array<OneD, NekDouble>(m_vel[i].size());
            Vmath::Vcopy(m_vel[i].size(), m_vel[i], 1, tmpVel[i], 1);
        }
        Array<OneD, NekDouble> angle(3, 0.);
        angle[2] = tmpVel[0][m_spacedim];
        m_frame.SetAngle(angle);
        m_frame.BodyToInerital(m_spacedim, tmpVel[1], tmpVel[1]);
        m_frame.BodyToInerital(m_spacedim, tmpVel[2], tmpVel[2]);
        tmpVel[2][0] -= tmpVel[1][m_spacedim] * tmpVel[1][1];
        tmpVel[2][1] += tmpVel[1][m_spacedim] * tmpVel[1][0];
        m_inertialPosition[0] +=
            m_timestep * (tmpVel[1][0] + 0.5 * m_timestep * tmpVel[2][0]);
        m_inertialPosition[1] +=
            m_timestep * (tmpVel[1][1] + 0.5 * m_timestep * tmpVel[2][1]);
        tmpVel[0][0] = m_inertialPosition[0];
        tmpVel[0][1] = m_inertialPosition[1];
    }
    if (m_isRoot && m_index % m_outputFrequency == 0)
    {
        if (!m_hasRotation)
        {
            tmpVel = m_vel;
        }
        m_outputStream << boost::format("%25.19e") % time << " ";
        for (size_t i = 0; i < tmpVel[0].size(); ++i)
        {
            m_outputStream << boost::format("%25.19e") % tmpVel[0][i] << " "
                           << boost::format("%25.19e") % tmpVel[1][i] << " "
                           << boost::format("%25.19e") % tmpVel[2][i] << " ";
        }
        m_outputStream << std::endl;
    }
    // set displacements at the current time step
    UpdateMRFData(MRFData);
    ++m_index;
}

void RigidSolver::UpdateMRFData(Array<OneD, NekDouble> &MRFData)
{
    /// MRFData:
    /// X, Y, Z, Theta_x, Theta_y, Theta_z, [inertial frame 0-5]
    /// U, V, W, Omega_x, Omega_y, Omega_z, [body frame 6-11]
    /// A_x, A_y, A_z, DOmega_x, DOmega_y, DOmega_z, [body frame 12-17]
    /// pivot_x, pivot_y, pivot_z, [body frame]
    for (int i = 0; i < 18; ++i)
    {
        MRFData[i] = 0.0;
    }

    if (m_hasRotation)
    {
        MRFData[0] = m_inertialPosition[0];
        MRFData[1] = m_inertialPosition[1];
    }
    else
    {
        MRFData[0] = m_vel[0][0];
        MRFData[1] = m_vel[0][1];
    }

    if (m_hasRotation)
    {
        MRFData[5] = m_vel[0][m_spacedim];
    }

    for (int i = 0; i < m_spacedim; ++i)
    {
        MRFData[i + 6]  = m_vel[1][i];
        MRFData[i + 12] = m_vel[2][i];
    }

    if (m_hasRotation)
    {
        MRFData[11] = m_vel[1][m_spacedim];
        MRFData[17] = m_vel[2][m_spacedim];
    }
}

void RigidSolver::SolveBodyMotion(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                                  const Array<OneD, NekDouble> &forcebody,
                                  std::map<int, NekDouble> &Dirs)
{
    if (0 == m_solveType)
    {
        SolveInertialFrame(bodyVel, forcebody, Dirs);
    }
    else if (1 == m_solveType)
    {
        SolveRotOneFree(bodyVel, forcebody, Dirs);
    }
    else if (2 == m_solveType)
    {
        SolveBodyFrame(bodyVel, forcebody, Dirs);
    }
    else
    {
        ASSERTL0(false, "Unsupported rigid solver type.");
    }
}

void RigidSolver::SolveInertialFrame(
    Array<OneD, Array<OneD, NekDouble>> &bodyVel,
    const Array<OneD, NekDouble> &forcebody, std::map<int, NekDouble> &Dirs)
{
    // only translational motion or prescribed tranlation with rotational
    // solve in absolute frame and transform to body frame
    m_bodySolver.SolvePrescribed(bodyVel, Dirs);
    if (m_hasFreeMotion)
    {
        Array<OneD, NekDouble> force(6, 0.);
        for (int i = 0; i < m_spacedim; ++i)
        {
            force[i] = forcebody[i] + m_extForceXYZ[i];
        }
        force[m_spacedim] = forcebody[5] + m_extForceXYZ[5];
        m_bodySolver.SolveFreeFixMat(bodyVel, force);
    }
}

// with rotational and one free tranlation
void RigidSolver::SolveRotOneFree(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                                  const Array<OneD, NekDouble> &forcebody,
                                  std::map<int, NekDouble> &Dirs)
{
    // one direction free
    NekDouble uy = 0.;
    if (Dirs.find(1) != Dirs.end())
    {
        uy = Dirs[1];
        Dirs.erase(1);
    }
    m_bodySolver.SolvePrescribed(bodyVel, Dirs);
    Array<OneD, NekDouble> force(6, 0.), angle(3, 0.);
    angle[2] = bodyVel[0][m_spacedim];
    m_frame.SetAngle(angle);
    m_frame.IneritalToBody(3, m_extForceXYZ, force);
    for (int i = 0; i < m_spacedim; ++i)
    {
        force[i] = forcebody[i] + force[i];
    }
    force[0] += m_mass * bodyVel[1][m_spacedim] * bodyVel[1][m_spacedim] *
                m_pivotdistance;
    force[1] -= m_mass * bodyVel[2][m_spacedim] * m_pivotdistance;
    m_bodySolver.SolveOneFree(bodyVel, force, angle, uy, m_mass);
}

// with rotational and all free tranlation
void RigidSolver::SolveBodyFrame(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                                 const Array<OneD, NekDouble> &forcebody,
                                 std::map<int, NekDouble> &Dirs)
{
    m_bodySolver.SolvePrescribed(bodyVel, Dirs); // at most 1 rotation
    Array<OneD, NekDouble> force(6, 0.), angle(3, 0.);
    if (!m_hasFreeMotion)
    {
        angle[2] = bodyVel[0][m_spacedim];
        m_frame.SetAngle(angle);
        m_frame.IneritalToBody(m_spacedim, bodyVel[1], bodyVel[1]);
        m_frame.IneritalToBody(m_spacedim, bodyVel[2], bodyVel[2]);
        bodyVel[2][0] += bodyVel[1][m_spacedim] * bodyVel[1][1];
        bodyVel[2][1] -= bodyVel[1][m_spacedim] * bodyVel[1][0];
        return;
    }
    if (m_dirDoFs.find(m_spacedim) != m_dirDoFs.end())
    {
        // known rotation
        angle[2] = bodyVel[0][m_spacedim];
        m_frame.SetAngle(angle);
        m_frame.IneritalToBody(m_spacedim, m_extForceXYZ, force);
        for (int i = 0; i < m_spacedim; ++i)
        {
            force[i] = forcebody[i] + force[i];
        }
        force[0] += m_mass * bodyVel[1][m_spacedim] * bodyVel[1][m_spacedim] *
                    m_pivotdistance;
        force[1] -= m_mass * bodyVel[2][m_spacedim] * m_pivotdistance;
        m_bodySolver.SolveFreeVarMat(bodyVel, force, m_mass);
    }
    else
    {
        // all free
        Array<OneD, Array<OneD, NekDouble>> tmpbodyVel(bodyVel.size());
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            tmpbodyVel[i] = Array<OneD, NekDouble>(bodyVel[i].size());
            Vmath::Vcopy(bodyVel[i].size(), bodyVel[i], 1, tmpbodyVel[i], 1);
        }
        for (int iter = 0; iter < 2; ++iter)
        {
            if (iter > 0)
            {
                for (size_t i = 0; i < bodyVel.size(); ++i)
                {
                    Vmath::Vcopy(bodyVel[i].size() - 1, bodyVel[i], 1,
                                 tmpbodyVel[i], 1);
                }
            }
            angle[2] = tmpbodyVel[0][m_spacedim];
            m_frame.SetAngle(angle);
            m_frame.IneritalToBody(3, m_extForceXYZ, force);
            for (int i = 0; i < m_spacedim; ++i)
            {
                force[i] = forcebody[i] + force[i];
            }
            force[m_spacedim] = forcebody[5] + m_extForceXYZ[5];
            m_bodySolver.SolveFreeVarMat(tmpbodyVel, force, m_mass);
        }
        // copy final results
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            Vmath::Vcopy(bodyVel[i].size(), tmpbodyVel[i], 1, bodyVel[i], 1);
        }
    }
}

void Newmark_BetaSolver::SolveFreeVarMat(Array<OneD, Array<OneD, NekDouble>> u,
                                         Array<OneD, NekDouble> force,
                                         const NekDouble mass)
{
    Array<OneD, NekDouble> bm(m_motionDofs, 0.);
    Array<OneD, NekDouble> bk(m_motionDofs, 0.);
    double *dMatrix = new double[m_motionDofs * m_motionDofs];
    int *ipiv       = new int[m_motionDofs];
    double *drhs    = new double[m_motionDofs];
    int info;

    for (int j = 0; j < m_motionDofs; ++j)
    {
        int j1 = m_index[j];
        bm[j]  = m_coeffs[0] * u[1][j1] + m_coeffs[1] * u[2][j1];
        bk[j]  = u[0][j1] + m_coeffs[3] * u[1][j1] + m_coeffs[4] * u[2][j1];
    }
    Array<OneD, NekDouble> rhs(m_motionDofs, 0.);
    for (int i = 0; i < m_motionDofs; ++i)
    {
        rhs[i] = force[m_index[i]];
        for (int j = 0; j < m_motionDofs; ++j)
        {
            rhs[i] += m_M[i][j] * bm[j] - m_K[i][j] * bk[j];
        }
        for (int j = m_motionDofs; j < m_rows; ++j)
        {
            int j1 = m_index[j];
            rhs[i] -= m_M[i][j] * u[2][j1] + m_C[i][j] * u[1][j1] +
                      m_K[i][j] * u[0][j1];
        }
    }
    for (int i = 0; i < m_motionDofs; ++i)
    {
        drhs[i] = rhs[i];
        for (int j = 0; j < m_motionDofs; ++j)
        {
            dMatrix[j * m_motionDofs + i] = m_Matrix[i][j];
        }
    }
    dMatrix[1] += mass * u[1][m_rows - 1];
    dMatrix[m_motionDofs] += -mass * u[1][m_rows - 1];

    Lapack::DoSgetrf(m_motionDofs, m_motionDofs, dMatrix, m_motionDofs, ipiv,
                     info);
    Lapack::Dgetrs('N', m_motionDofs, 1, dMatrix, m_motionDofs, ipiv, drhs,
                   m_motionDofs, info);
    for (int j = 0; j < m_motionDofs; ++j)
    {
        int j1   = m_index[j];
        u[1][j1] = drhs[j];
    }

    for (int j = 0; j < m_motionDofs; ++j)
    {
        int j1   = m_index[j];
        u[0][j1] = m_coeffs[2] * u[1][j1] + bk[j];
        u[2][j1] = m_coeffs[0] * u[1][j1] - bm[j];
    }

    delete[] dMatrix;
    delete[] drhs;
    delete[] ipiv;
}

void RigidSolver::SetNewmarkBetaSolver(Array<OneD, NekDouble> &AddedMass)
{
    int NumDof = m_spacedim + 1;
    if (AddedMass.size() >= NumDof * NumDof)
    {
        Vmath::Vadd(NumDof * NumDof, AddedMass, 1, m_M, 1, m_M, 1);
    }
    m_bodySolver.SetNewmarkBeta(m_beta, m_gamma, m_timestep, m_M, m_C, m_K,
                                m_dirDoFs, m_solveType);
}

void RigidSolver::SetInitialConditions(
    const LibUtilities::SessionReaderSharedPtr session,
    Array<OneD, NekDouble> MRFData)
{
    NekDouble time = 0.;
    std::map<std::string, std::string> fieldMetaDataMap;
    std::vector<std::string> strFrameData = {
        "X",   "Y",   "Z",   "Theta_x",  "Theta_y",  "Theta_z",
        "U",   "V",   "W",   "Omega_x",  "Omega_y",  "Omega_z",
        "A_x", "A_y", "A_z", "DOmega_x", "DOmega_y", "DOmega_z"};
    std::map<std::string, NekDouble> fileData;
    if (session->DefinesFunction("InitialConditions"))
    {
        for (int i = 0; i < session->GetVariables().size(); ++i)
        {
            if (session->GetFunctionType("InitialConditions",
                                         session->GetVariable(i)) ==
                LibUtilities::eFunctionTypeFile)
            {
                std::string filename = session->GetFunctionFilename(
                    "InitialConditions", session->GetVariable(i));
                fs::path pfilename(filename);
                // redefine path for parallel file which is in directory
                if (fs::is_directory(pfilename))
                {
                    fs::path metafile("Info.xml");
                    fs::path fullpath = pfilename / metafile;
                    filename          = LibUtilities::PortablePath(fullpath);
                }
                LibUtilities::FieldIOSharedPtr fld =
                    LibUtilities::FieldIO::CreateForFile(session, filename);
                fld->ImportFieldMetaData(filename, fieldMetaDataMap);

                // check to see if time is defined
                if (fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
                {
                    if (fieldMetaDataMap.find("Time") != fieldMetaDataMap.end())
                    {
                        time = std::stod(fieldMetaDataMap["Time"]);
                    }
                    fileData.clear();
                    for (auto &var : strFrameData)
                    {
                        if (fieldMetaDataMap.find(var) !=
                            fieldMetaDataMap.end())
                        {
                            fileData[var] = std::stod(fieldMetaDataMap[var]);
                        }
                    }
                    if (fileData.size() == strFrameData.size())
                    {
                        break;
                    }
                }
            }
        }
    }
    if (session->DefinesCmdLineArgument("set-start-time"))
    {
        time = std::stod(
            session->GetCmdLineArgument<std::string>("set-start-time"));
    }
    if (fileData.size() == strFrameData.size())
    {
        int NumDofm1 = m_vel[0].size() - 1;
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vel[0][i] = fileData[strFrameData[i]];
            m_vel[1][i] = fileData[strFrameData[i + 6]];
            m_vel[2][i] = fileData[strFrameData[i + 12]];
        }
        m_vel[0][NumDofm1] = fileData[strFrameData[5]];
        m_vel[1][NumDofm1] = fileData[strFrameData[11]];
        m_vel[2][NumDofm1] = fileData[strFrameData[17]];
        if (m_hasRotation)
        {
            m_inertialPosition[0] = m_vel[0][0];
            m_inertialPosition[1] = m_vel[0][1];
        }
    }
    std::map<int, NekDouble> Dirs;
    UpdatePrescribed(time, Dirs);
    SetInitialConditions(Dirs);
    UpdateMRFData(MRFData);
    // output initial status for rigid body
    if (m_isRoot)
    {
        Array<OneD, Array<OneD, NekDouble>> tmpVel;
        if (m_hasRotation)
        {
            tmpVel = Array<OneD, Array<OneD, NekDouble>>(m_vel.size());
            for (size_t i = 0; i < m_vel.size(); ++i)
            {
                tmpVel[i] = Array<OneD, NekDouble>(m_vel[i].size());
                Vmath::Vcopy(m_vel[i].size(), m_vel[i], 1, tmpVel[i], 1);
            }
            // transform u, du to body frame
            Array<OneD, NekDouble> angle(3, 0.);
            angle[2] = tmpVel[0][m_spacedim];
            m_frame.SetAngle(angle);
            m_frame.BodyToInerital(m_spacedim, tmpVel[1], tmpVel[1]);
            m_frame.BodyToInerital(m_spacedim, tmpVel[2], tmpVel[2]);
            tmpVel[2][0] -= tmpVel[1][m_spacedim] * tmpVel[1][1];
            tmpVel[2][1] += tmpVel[1][m_spacedim] * tmpVel[1][0];
            tmpVel[0][0] = m_inertialPosition[0];
            tmpVel[0][1] = m_inertialPosition[1];
        }
        else
        {
            tmpVel = m_vel;
        }
        m_outputStream << boost::format("%25.19e") % time << " ";
        for (size_t i = 0; i < tmpVel[0].size(); ++i)
        {
            m_outputStream << boost::format("%25.19e") % tmpVel[0][i] << " "
                           << boost::format("%25.19e") % tmpVel[1][i] << " "
                           << boost::format("%25.19e") % tmpVel[2][i] << " ";
        }
        m_outputStream << std::endl;
    }
}

void RigidSolver::GetFilterInfo(
    const LibUtilities::SessionReaderSharedPtr session,
    std::map<std::string, std::string> &vParams)
{
    TiXmlElement *pForce      = session->GetElement("Nektar/RIGIDSOLVER");
    const TiXmlElement *param = pForce->FirstChildElement("BOUNDARY");
    ASSERTL0(param, "Body surface should be assigned");

    vParams["Boundary"]           = param->GetText();
    const TiXmlElement *pivotElmt = pForce->FirstChildElement("PIVOTPOINT");
    if (pivotElmt)
    {
        std::string pstr = pivotElmt->GetText();
        std::replace(pstr.begin(), pstr.end(), ',', ' ');
        vParams["MomentPoint"] = pstr;
    }
}

void RigidSolver::SetInitialConditions(std::map<int, NekDouble> &Dirs)
{
    for (auto it : Dirs)
    {
        int NumDof  = m_vel[0].size();
        int NumDof2 = NumDof << 1;
        if (it.first < m_vel[0].size())
        {
            m_vel[1][it.first] = it.second;
        }
        else if (it.first < NumDof2)
        {
            m_vel[0][it.first - NumDof] = it.second;
        }
        else
        {
            m_vel[2][it.first - NumDof2] = it.second;
        }
    }
}

void Newmark_BetaSolver::SetNewmarkBeta(NekDouble beta, NekDouble gamma,
                                        NekDouble dt, Array<OneD, NekDouble> M,
                                        Array<OneD, NekDouble> C,
                                        Array<OneD, NekDouble> K,
                                        std::set<int> DirDoFs, int solveType)
{
    m_coeffs    = Array<OneD, NekDouble>(5, 0.);
    m_coeffs[0] = 1. / (gamma * dt);
    m_coeffs[1] = 1. / gamma - 1.;
    m_coeffs[2] = beta * dt / gamma;
    m_coeffs[3] = dt * (1. - beta / gamma);
    m_coeffs[4] = (0.5 - beta / gamma) * dt * dt;

    m_rows = sqrt(M.size());
    m_index.resize(m_rows, -1);
    m_motionDofs = 0;
    for (int i = 0; i < m_rows; ++i)
    {
        if (DirDoFs.find(i) == DirDoFs.end())
        {
            m_index[m_motionDofs++] = i;
        }
    }
    for (int i = 0, count = m_motionDofs; i < m_rows; ++i)
    {
        if (DirDoFs.find(i) != DirDoFs.end())
        {
            m_index[count++] = i;
        }
    }
    if (1 == solveType)
    {
        ASSERTL0(m_motionDofs == 2,
                 "2 Dofs if body is free only in x direction.");
    }
    if (m_motionDofs)
    {
        Array<OneD, NekDouble> temp;
        m_M = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        m_C = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        m_K = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        DNekMatSharedPtr inverseMatrix =
            MemoryManager<DNekMat>::AllocateSharedPtr(m_motionDofs,
                                                      m_motionDofs, 0.0, eFULL);
        for (int i = 0; i < m_motionDofs; ++i)
        {
            m_M[i]     = Array<OneD, NekDouble>(m_rows, 0.);
            m_C[i]     = Array<OneD, NekDouble>(m_rows, 0.);
            m_K[i]     = Array<OneD, NekDouble>(m_rows, 0.);
            int offset = m_index[i] * m_rows;
            for (int j = 0; j < m_rows; ++j)
            {
                int ind   = offset + m_index[j];
                m_M[i][j] = M[ind];
                m_C[i][j] = C[ind];
                m_K[i][j] = K[ind];
                NekDouble value =
                    m_coeffs[0] * M[ind] + C[ind] + m_coeffs[2] * K[ind];
                if (j < m_motionDofs)
                {
                    inverseMatrix->SetValue(i, j, value);
                }
            }
        }

        m_Matrix = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        if (solveType == 0)
        {
            inverseMatrix->Invert();
        }
        for (int i = 0; i < m_motionDofs; ++i)
        {
            m_Matrix[i] = Array<OneD, NekDouble>(m_motionDofs, 0.);
            for (int j = 0; j < m_rows; ++j)
            {
                if (j < m_motionDofs)
                {
                    m_Matrix[i][j] = inverseMatrix->GetValue(i, j);
                }
            }
        }
    }
}

void Newmark_BetaSolver::SolvePrescribed(
    Array<OneD, Array<OneD, NekDouble>> u,
    std::map<int, NekDouble> motionPrescribed)
{
    for (int i = 0; i < m_rows; ++i)
    {
        if (motionPrescribed.find(i) != motionPrescribed.end())
        {
            int i0       = i + m_rows;
            int i2       = i0 + m_rows;
            NekDouble bm = 0., bk = 0.;
            if (motionPrescribed.find(i2) == motionPrescribed.end())
            {
                bm = m_coeffs[0] * u[1][i] + m_coeffs[1] * u[2][i];
            }
            if (motionPrescribed.find(i0) == motionPrescribed.end())
            {
                bk = u[0][i] + m_coeffs[3] * u[1][i] + m_coeffs[4] * u[2][i];
            }

            u[1][i] = motionPrescribed[i];
            if (motionPrescribed.find(i2) == motionPrescribed.end())
            {
                u[2][i] = m_coeffs[0] * u[1][i] - bm;
            }
            else
            {
                u[2][i] = motionPrescribed[i2];
            }
            if (motionPrescribed.find(i0) == motionPrescribed.end())
            {
                u[0][i] = m_coeffs[2] * u[1][i] + bk;
            }
            else
            {
                u[0][i] = motionPrescribed[i0];
            }
        }
    }
}

void Newmark_BetaSolver::SolveFreeFixMat(Array<OneD, Array<OneD, NekDouble>> u,
                                         Array<OneD, NekDouble> force)
{
    if (m_motionDofs)
    {
        Array<OneD, NekDouble> bm(m_motionDofs, 0.);
        Array<OneD, NekDouble> bk(m_motionDofs, 0.);
        for (int j = 0; j < m_motionDofs; ++j)
        {
            int j1 = m_index[j];
            bm[j]  = m_coeffs[0] * u[1][j1] + m_coeffs[1] * u[2][j1];
            bk[j]  = u[0][j1] + m_coeffs[3] * u[1][j1] + m_coeffs[4] * u[2][j1];
        }
        Array<OneD, NekDouble> rhs(m_motionDofs, 0.);
        for (int i = 0; i < m_motionDofs; ++i)
        {
            rhs[i] = force[m_index[i]];
            for (int j = 0; j < m_motionDofs; ++j)
            {
                rhs[i] += m_M[i][j] * bm[j] - m_K[i][j] * bk[j];
            }
            for (int j = m_motionDofs; j < m_rows; ++j)
            {
                int j1 = m_index[j];
                rhs[i] -= m_M[i][j] * u[2][j1] + m_C[i][j] * u[1][j1] +
                          m_K[i][j] * u[0][j1];
            }
        }
        for (int j = 0; j < m_motionDofs; ++j)
        {
            int j1   = m_index[j];
            u[1][j1] = Vmath::Dot(m_motionDofs, m_Matrix[j], 1, rhs, 1);
            u[0][j1] = m_coeffs[2] * u[1][j1] + bk[j];
            u[2][j1] = m_coeffs[0] * u[1][j1] - bm[j];
        }
    }
}

/**
 * e_x \cdot M [du0, du1, du2]  + m (e_x \times Omega) \cdot u = F \cdot e_x
 * e_x = (c, -s, 0)
 * (c M00 - s M10, c M01 - s M11, c M02 - s M 12) [du0, du1, du2]^T +
 * (-s Omega, -c Omega, 0) [u0, u1, u2]^T =
 * c F0 - s F1
 * equation 2, e_y \cdot (u0, u1) = uy = s * u0 + c * u1
 **/
void Newmark_BetaSolver::SolveOneFree(Array<OneD, Array<OneD, NekDouble>> u,
                                      Array<OneD, NekDouble> force,
                                      const Array<OneD, NekDouble> theta,
                                      const NekDouble uy, const NekDouble mass)
{
    int iOmega  = m_rows - 1;
    NekDouble c = cos(theta[2]), s = sin(theta[2]);
    NekDouble C00, C01, C10, C11, M0, M1, M2, F0, F1;
    M0  = c * m_M[0][0] - s * m_M[1][0];
    M1  = c * m_M[0][1] - s * m_M[1][1];
    M2  = c * m_M[0][iOmega] - s * m_M[1][iOmega];
    C00 = m_coeffs[0] * M0 - s * u[1][iOmega] * mass;
    C01 = m_coeffs[0] * M1 - c * u[1][iOmega] * mass;
    C10 = s;
    C11 = c;
    F1  = uy;
    // F0
    Array<OneD, NekDouble> bm(m_motionDofs, 0.);
    Array<OneD, NekDouble> bk(m_motionDofs, 0.);
    for (int j = 0; j < m_motionDofs; ++j)
    {
        bm[j] = m_coeffs[0] * u[1][j] + m_coeffs[1] * u[2][j];
        bk[j] = u[0][j] + m_coeffs[3] * u[1][j] + m_coeffs[4] * u[2][j];
    }
    F0 = c * force[0] - s * force[1];
    F0 = F0 + M0 * bm[0] + M1 * bm[1] - M2 * u[2][iOmega];
    // solve
    NekDouble det = 1. / (C00 * C11 - C01 * C10);
    u[1][0]       = det * (C11 * F0 - C01 * F1);
    u[1][1]       = det * (-C10 * F0 + C00 * F1);
    for (int j = 0; j < m_motionDofs; ++j)
    {
        u[0][j] = m_coeffs[2] * u[1][j] + bk[j];
        u[2][j] = m_coeffs[0] * u[1][j] - bm[j];
    }
}

FrameTransform::FrameTransform()
{
    m_matrix = Array<OneD, NekDouble>(2, 0.);
}

void FrameTransform::SetAngle(const Array<OneD, NekDouble> theta)
{
    m_matrix[0] = cos(theta[2]);
    m_matrix[1] = sin(theta[2]);
}

void FrameTransform::BodyToInerital(const int dim,
                                    const Array<OneD, NekDouble> &body,
                                    Array<OneD, NekDouble> &inertial)
{
    NekDouble xi = body[0], eta = body[1];
    inertial[0] = xi * m_matrix[0] - eta * m_matrix[1];
    inertial[1] = eta * m_matrix[0] + xi * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            inertial[i] = body[i];
        }
    }
}

void FrameTransform::IneritalToBody(const int dim,
                                    const Array<OneD, NekDouble> &inertial,
                                    Array<OneD, NekDouble> &body)
{
    NekDouble x = inertial[0], y = inertial[1];
    body[0] = x * m_matrix[0] + y * m_matrix[1];
    body[1] = y * m_matrix[0] - x * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            body[i] = inertial[i];
        }
    }
}

} // namespace Nektar
