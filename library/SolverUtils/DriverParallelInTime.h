///////////////////////////////////////////////////////////////////////////////
//
// File DriverParallelInTime.h
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
// Description: Driver class for the parallel-in-time solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERPARALLELINTIME_H
#define NEKTAR_SOLVERUTILS_DRIVERPARALLELINTIME_H

#include <SolverUtils/Driver.h>
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar::SolverUtils
{

/// Base class for the development of parallel-in-time solvers.
class DriverParallelInTime : public Driver
{
public:
protected:
    /// Constructor
    SOLVER_UTILS_EXPORT DriverParallelInTime(
        const LibUtilities::SessionReaderSharedPtr pSession,
        const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    SOLVER_UTILS_EXPORT ~DriverParallelInTime() override = default;

    /// Virtual function for initialisation implementation.
    SOLVER_UTILS_EXPORT void v_InitObject(
        std::ostream &out = std::cout) override;

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT void v_Execute(std::ostream &out = std::cout) override;

    void SetParallelInTimeEquationSystem(std::string AdvectiveType);

    void GetParametersFromSession(void);

    void InitialiseEqSystem(bool turnoff_output);

    void InitialiseInterpolationField(void);

    void PrintSolverInfo(std::ostream &out = std::cout);

    void PrintHeader(const std::string &title, const char c);

    void RecvFromPreviousProc(Array<OneD, Array<OneD, NekDouble>> &array,
                              int &convergence);

    void RecvFromPreviousProc(Array<OneD, NekDouble> &array);

    void SendToNextProc(Array<OneD, Array<OneD, NekDouble>> &array,
                        int &convergence);

    void SendToNextProc(Array<OneD, NekDouble> &array);

    void CopySolutionVector(const Array<OneD, const Array<OneD, NekDouble>> &in,
                            Array<OneD, Array<OneD, NekDouble>> &out);

    void CopyFromPhysField(const size_t timeLevel,
                           Array<OneD, Array<OneD, NekDouble>> &out);

    void CopyToPhysField(const size_t timeLevel,
                         const Array<OneD, const Array<OneD, NekDouble>> &in);

    void UpdateFieldCoeffs(const size_t timeLevel,
                           const Array<OneD, const Array<OneD, NekDouble>> &in =
                               NullNekDoubleArrayOfArray);

    void EvaluateExactSolution(const size_t timeLevel, const NekDouble &time);

    void SolutionConvergenceMonitoring(const size_t timeLevel,
                                       const size_t iter);

    void SolutionConvergenceSummary(const size_t timeLevel);

    void UpdateErrorNorm(const size_t timeLevel, const bool normalized);

    void PrintErrorNorm(const size_t timeLevel, const bool normalized);

    NekDouble vL2ErrorMax(void);

    NekDouble EstimateCommunicationTime(
        Array<OneD, Array<OneD, NekDouble>> &buffer1,
        Array<OneD, Array<OneD, NekDouble>> &buffer2);

    void Interpolate(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &infield,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &outfield,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    /// Total time integration interval.
    NekDouble m_totalTime;

    /// Time integration interval per chunk.
    NekDouble m_chunkTime;

    /// Local time.
    NekDouble m_time;
    NekDouble m_time0;

    /// Number of time chunks.
    size_t m_numChunks;

    /// Rank in time.
    size_t m_chunkRank;

    /// Maximum number of parallel-in-time iteration.
    size_t m_iterMaxPIT;

    // Number of windows for parallel-in-time time iteration.
    size_t m_numWindowsPIT;

    /// Using exact solution to compute error norms.
    bool m_exactSolution;

    /// ParallelInTime tolerance.
    NekDouble m_tolerPIT;

    /// Number of variables.
    size_t m_nVar;

    /// Number of time levels.
    size_t m_nTimeLevel;

    /// Number of time steps for each time level.
    Array<OneD, size_t> m_nsteps;

    /// Time step for each time level.
    Array<OneD, NekDouble> m_timestep;

    /// Number of dof for each time level.
    Array<OneD, size_t> m_npts;

    /// Equation system to solve.
    Array<OneD, std::shared_ptr<UnsteadySystem>> m_EqSys;

    /// Storage for parallel-in-time iteration.
    Array<OneD, NekDouble> m_vL2Errors;
    Array<OneD, NekDouble> m_vLinfErrors;
    Array<OneD, Array<OneD, NekDouble>> m_exactsoln;
};

} // namespace Nektar::SolverUtils

#endif // NEKTAR_SOLVERUTILS_DRIVERPARALLELINTIME_H
