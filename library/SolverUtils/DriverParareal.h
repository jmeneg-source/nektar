///////////////////////////////////////////////////////////////////////////////
//
// File DriverParareal.h
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
// Description: Driver class for the parareal solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERPARAREAL_H
#define NEKTAR_SOLVERUTILS_DRIVERPARAREAL_H

#include <SolverUtils/DriverParallelInTime.h>

namespace Nektar::SolverUtils
{

/// Base class for the development of solvers.
class DriverParareal : public DriverParallelInTime
{
public:
    friend class MemoryManager<DriverParareal>;

    /// Creates an instance of this class
    static DriverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        DriverSharedPtr p =
            MemoryManager<DriverParareal>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    /// Constructor
    SOLVER_UTILS_EXPORT DriverParareal(
        const LibUtilities::SessionReaderSharedPtr pSession,
        const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    SOLVER_UTILS_EXPORT ~DriverParareal() override = default;

    /// Virtual function for initialisation implementation.
    SOLVER_UTILS_EXPORT void v_InitObject(
        std::ostream &out = std::cout) override;

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT void v_Execute(std::ostream &out = std::cout) override;

    static std::string driverLookupId;

private:
    void AllocateMemory(void);

    void AssertParameters(void);

    void UpdateInitialConditionFromSolver(const size_t timeLevel);

    void UpdateSolverInitialCondition(const size_t timeLevel);

    void UpdateSolution(const size_t timeLevel, const NekDouble time,
                        const size_t nstep, const size_t wd, const size_t iter);

    void CorrectionWithOldCoarseSolution(void);

    void CorrectionWithNewCoarseSolution(void);

    void InterpolateCoarseSolution(void);

    void ApplyWindowing(const size_t w);

    void CopyConvergedCheckPoints(const size_t w, const size_t k);

    void WriteTimeChunkOuput(void);

    static constexpr size_t m_fineLevel   = 0;
    static constexpr size_t m_coarseLevel = 1;
    Array<OneD, Array<OneD, NekDouble>> m_initialCondition;
    Array<OneD, Array<OneD, NekDouble>> m_fineSolution;
    Array<OneD, Array<OneD, NekDouble>> m_coarseSolution;
};

} // namespace Nektar::SolverUtils

#endif // NEKTAR_SOLVERUTILS_DRIVERPARAREAL_H
