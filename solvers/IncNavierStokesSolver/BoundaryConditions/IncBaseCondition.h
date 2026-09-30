///////////////////////////////////////////////////////////////////////////////
//
// File: IncBaseCondition.h
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
// Description: Abstract base class for Extrapolate.
//    params {
//    "X",  "Y", "Z", "Theta_x", "Theta_y", "Theta_z", "X0",     "Y0",
//    "Z0", "U", "V", "W",       "Omega_x", "Omega_y", "Omega_z", "A_x", "A_y",
//    "A_z", "DOmega_x", "DOmega_y", "DOmega_z"};
// (x, y, z), (x0, y0, z0), (theta, omega, domega) vector in absolute frame
// (U, V, W), (A_x, A_y, A_z) vector in body frame
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_INCBASECONDITION_H
#define NEKTAR_SOLVERS_INCBASECONDITION_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{

class IncBaseCondition;
typedef std::shared_ptr<IncBaseCondition> IncBaseConditionSharedPtr;

typedef LibUtilities::NekFactory<
    std::string, IncBaseCondition, const LibUtilities::SessionReaderSharedPtr,
    Array<OneD, MultiRegions::ExpListSharedPtr>,
    Array<OneD, SpatialDomains::BoundaryConditionShPtr>,
    Array<OneD, MultiRegions::ExpListSharedPtr>, int, int, int>
    IncBCFactory;

SOLVER_UTILS_EXPORT IncBCFactory &GetIncBCFactory();

class IncBaseCondition
{
public:
    virtual ~IncBaseCondition() = default;

    void Initialise(const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        v_Initialise(pSession);
    }

    void Update(const Array<OneD, const Array<OneD, NekDouble>> &fields,
                const Array<OneD, const Array<OneD, NekDouble>> &Adv,
                std::map<std::string, NekDouble> &params)
    {
        v_Update(fields, Adv, params);
    }

protected:
    IncBaseCondition(const LibUtilities::SessionReaderSharedPtr pSession,
                     Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
                     Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
                     Array<OneD, MultiRegions::ExpListSharedPtr> exp, int nbnd,
                     int spacedim, int bnddim);

    virtual void v_Initialise(
        const LibUtilities::SessionReaderSharedPtr &pSession);

    virtual void v_Update(
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>>
            &fields,
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Adv,
        [[maybe_unused]] std::map<std::string, NekDouble> &params)
    {
    }

    void ExtrapolateArray(
        const int numCalls,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &array);

    void RollOver(Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &input);

    void AddVisPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble>> &fields,
        Array<OneD, Array<OneD, NekDouble>> &N,
        std::map<std::string, NekDouble> &params);

    void AddExtrapVisPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble>> &fields,
        Array<OneD, Array<OneD, NekDouble>> &N,
        std::map<std::string, NekDouble> &params);

    void AddRigidBodyAcc(Array<OneD, Array<OneD, NekDouble>> &N,
                         std::map<std::string, NekDouble> &params, int npts0);

    void RigidBodyVelocity(Array<OneD, Array<OneD, NekDouble>> &velocities,
                           std::map<std::string, NekDouble> &params, int npts0);
    void InitialiseCoords(std::map<std::string, NekDouble> &params);
    void SetNumPointsOnPlane0(int &npointsPlane0);

    int m_spacedim;
    /// bounday dimensionality
    int m_bnddim;
    int m_nbnd;

    int m_numCalls;
    int m_intSteps;
    std::map<int, SpatialDomains::BoundaryConditionShPtr> m_BndConds;
    std::map<int, MultiRegions::ExpListSharedPtr> m_BndExp;
    int m_npoints;
    Array<OneD, Array<OneD, NekDouble>> m_coords;
    MultiRegions::ExpListSharedPtr m_bndElmtExps;
    MultiRegions::ExpListSharedPtr m_field;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_extrapArray;
    int m_pressure;

    static NekDouble StifflyStable_Betaq_Coeffs[3][3];
    static NekDouble StifflyStable_Alpha_Coeffs[3][3];
    static NekDouble StifflyStable_Gamma0_Coeffs[3];
    std::string classname;
};

} // namespace Nektar

#endif
