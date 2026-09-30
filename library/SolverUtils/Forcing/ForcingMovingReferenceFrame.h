///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.h
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
// Description: Allows for a moving frame of reference, through adding c * du/dx
// to the body force, where c is the frame velocity vector
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGMOVINGREFERENCEFRAME
#define NEKTAR_SOLVERUTILS_FORCINGMOVINGREFERENCEFRAME

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <string>

namespace Nektar::SolverUtils
{

class ForcingMovingReferenceFrame : public Forcing
{

public:
    friend class MemoryManager<ForcingMovingReferenceFrame>;

    /// Creates an instance of this class
    SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
    {
        ForcingSharedPtr p =
            MemoryManager<ForcingMovingReferenceFrame>::AllocateSharedPtr(
                pSession, pEquation);
        p->InitObject(pFields, pNumForcingFields, pForce);
        return p;
    }

    /// Name of the class
    static std::string classNameBody;

protected:
    SOLVER_UTILS_EXPORT void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields,
        const TiXmlElement *pForce) override;

    SOLVER_UTILS_EXPORT void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble &time) override;

    SOLVER_UTILS_EXPORT void v_PreApply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble &time) override;

private:
    ForcingMovingReferenceFrame(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation);

    ~ForcingMovingReferenceFrame(void) override;

    void addRotation(int npoints,
                     const Array<OneD, Array<OneD, NekDouble>> &inarray0,
                     NekDouble angVelScale,
                     const Array<OneD, Array<OneD, NekDouble>> &inarray1,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);
    void UpdateMRFStatus(MultiRegions::ExpListSharedPtr field);
    // pivot point
    Array<OneD, NekDouble> m_pivotPoint;
    // a boolean switch indicating for which direction the velocities are
    // available. The available velocites are in body frame
    Array<OneD, bool> m_hasVel;
    Array<OneD, bool> m_hasOmega;
    bool m_hasRotation; // m_hasOmega[0] || m_hasOmega[1] || m_hasOmega[2]
    // frame linear velocities in local translating-rotating frame
    Array<OneD, NekDouble> m_velxyz;
    // frame angular velocities in local translating-rotating frame
    Array<OneD, NekDouble> m_omegaxyz;
    // coordinate vector
    Array<OneD, Array<OneD, NekDouble>> m_coords;
    bool m_isH1d;
    bool m_hasPlane0;
    bool m_isH2d;
    int m_spacedim;
    std::shared_ptr<Nektar::SolverUtils::FluidInterface> m_FluidEq;
};

} // namespace Nektar::SolverUtils

#endif
