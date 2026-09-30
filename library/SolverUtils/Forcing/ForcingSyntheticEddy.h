///////////////////////////////////////////////////////////////////////////////
//
// File:  ForcingSyntheticEddy.h
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
// Description: Derived base class - Synthetic turbulence generation.
//              This code implements the Synthetic Eddy Method (SEM).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGSYNTHETICEDDY
#define NEKTAR_SOLVERUTILS_FORCINGSYNTHETICEDDY

#include <SolverUtils/Forcing/Forcing.h>
#include <string>

namespace Nektar::SolverUtils
{

class ForcingSyntheticEddy : virtual public SolverUtils::Forcing
{
public:
    friend class MemoryManager<ForcingSyntheticEddy>;

    /// Creates an instance of this class
    SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
    {
        ForcingSharedPtr p =
            MemoryManager<ForcingSyntheticEddy>::AllocateSharedPtr(pSession,
                                                                   pEquation);
        p->InitObject(pFields, pNumForcingFields, pForce);
        return p;
    }

    /// Name of the class
    SOLVER_UTILS_EXPORT static std::string className;

protected:
    // Initial object
    SOLVER_UTILS_EXPORT void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields,
        const TiXmlElement *pForce) override;
    // Apply forcing term
    SOLVER_UTILS_EXPORT void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble &time) override;
    // Apply forcing term
    SOLVER_UTILS_EXPORT void v_ApplyCoeff(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble &time) override;
    // Set Cholesky decomposition of the Reynolds Stresses in the domain
    SOLVER_UTILS_EXPORT void SetCholeskyReyStresses(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);
    /// Set the Number of the Eddies in the box
    SOLVER_UTILS_EXPORT void SetNumberOfEddies();
    /// Set reference lengths
    SOLVER_UTILS_EXPORT void ComputeRefLenghts();
    /// Set Xi max
    SOLVER_UTILS_EXPORT void ComputeXiMax();
    /// Compute Random 1 or -1 value
    SOLVER_UTILS_EXPORT Array<OneD, Array<OneD, int>> GenerateRandomOneOrMinusOne();
    /// Compute Constant C
    SOLVER_UTILS_EXPORT NekDouble ComputeConstantC(int row, int col);
    /// Compute Gaussian
    SOLVER_UTILS_EXPORT NekDouble ComputeGaussian(NekDouble coord,
                                                  NekDouble xiMaxVal,
                                                  NekDouble constC = 1.0);
    /// Check if point is inside the box of eddies
    SOLVER_UTILS_EXPORT bool InsideBoxOfEddies(NekDouble coord0,
                                               NekDouble coord1,
                                               NekDouble coord2);
    /// Compute Stochastic Signal
    SOLVER_UTILS_EXPORT Array<OneD, Array<OneD, NekDouble>> ComputeStochasticSignal(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);
    /// Compute Velocity Fluctuation
    SOLVER_UTILS_EXPORT Array<OneD, Array<OneD, NekDouble>>
    ComputeVelocityFluctuation(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, Array<OneD, NekDouble>> stochasticSignal);
    /// Compute Characteristic Convective Turbulent Time
    SOLVER_UTILS_EXPORT Array<OneD, Array<OneD, NekDouble>> ComputeCharConvTurbTime(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);
    /// Compute Smoothing Factor
    SOLVER_UTILS_EXPORT Array<OneD, Array<OneD, NekDouble>> ComputeSmoothingFactor(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, Array<OneD, NekDouble>> convTurb);
    /// Set box of eddies mask
    SOLVER_UTILS_EXPORT void SetBoxOfEddiesMask(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);
    /// Initialise forcing term for each eddy
    SOLVER_UTILS_EXPORT void InitialiseForcingEddy(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);
    /// Compute the initial position of the eddies inside the box
    SOLVER_UTILS_EXPORT void ComputeInitialRandomLocationOfEddies();
    /// Update positions of the eddies inside the box
    SOLVER_UTILS_EXPORT void UpdateEddiesPositions();
    /// Remove eddy from forcing term
    SOLVER_UTILS_EXPORT void RemoveEddiesFromForcing(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);
    /// Compute the initial location of the eddies for the test case
    SOLVER_UTILS_EXPORT void ComputeInitialLocationTestCase();

    // Members
    // Expressions (functions) of the prescribed Reynolds stresses
    std::map<int, LibUtilities::EquationSharedPtr> m_R;
    /// Cholesky decomposition of the Reynolds Stresses
    /// for all points in the mesh
    Array<OneD, Array<OneD, NekDouble>> m_Cholesky;
    /// Bulk velocity
    NekDouble m_Ub;
    /// Standard deviation
    NekDouble m_sigma;
    /// Inlet length in the y-direction and z-direction
    Array<OneD, NekDouble> m_lyz;
    /// Center of the inlet plane
    Array<OneD, NekDouble> m_rc;
    /// Number of eddies in the box
    int m_N;
    /// Characteristic lenght scales
    Array<OneD, NekDouble> m_l;
    /// Reference lenghts
    Array<OneD, NekDouble> m_lref;
    /// Xi max
    Array<OneD, NekDouble> m_xiMax;
    // XiMaxMin - Value form Giangaspero et al. 2022
    NekDouble m_xiMaxMin = 0.4;
    /// Space dimension
    int m_spacedim;
    /// Ration of specific heats
    NekDouble m_gamma;
    /// Box of eddies mask
    Array<OneD, int> m_mask;
    /// Eddy position
    Array<OneD, Array<OneD, NekDouble>> m_eddyPos;
    /// Check when the forcing should be applied
    bool m_calcForcing{true};
    /// Eddies that add to the domain
    std::vector<unsigned int> m_eddiesIDForcing;
    /// Current time
    NekDouble m_currTime;
    /// Check for test case
    bool m_tCase;
    /// Forcing for each eddy
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ForcingEddy;

    SOLVER_UTILS_EXPORT ForcingSyntheticEddy(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation);
    SOLVER_UTILS_EXPORT ~ForcingSyntheticEddy(void) override = default;
};

} // namespace Nektar::SolverUtils

#endif
