///////////////////////////////////////////////////////////////////////////////
//
// File DriverPFASST.h
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
// Description: Driver class for the PFASST solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERPFASST_H
#define NEKTAR_SOLVERUTILS_DRIVERPFASST_H

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>
#include <SolverUtils/DriverParallelInTime.h>
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar::SolverUtils
{

typedef Array<OneD, Array<OneD, Array<OneD, NekDouble>>> SDCarray;

/// Base class for the development of solvers.
class DriverPFASST : public DriverParallelInTime
{
public:
    friend class MemoryManager<DriverPFASST>;

    /// Creates an instance of this class
    static DriverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        DriverSharedPtr p =
            MemoryManager<DriverPFASST>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    /// Constructor
    SOLVER_UTILS_EXPORT DriverPFASST(
        const LibUtilities::SessionReaderSharedPtr pSession,
        const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    SOLVER_UTILS_EXPORT ~DriverPFASST() override = default;

    /// Virtual function for initialisation implementation.
    SOLVER_UTILS_EXPORT void v_InitObject(
        std::ostream &out = std::cout) override;

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT void v_Execute(std::ostream &out = std::cout) override;

    static std::string driverLookupId;

private:
    void AssertParameters(void);

    void InitialiseSDCScheme(void);

    void SetTimeInterpolator(void);

    bool IsNotInitialCondition(const size_t n);

    void PropagateQuadratureSolutionAndResidual(const size_t timeLevel,
                                                const size_t index);

    void UpdateFirstQuadrature(const size_t timeLevel);

    void RunSweep(const NekDouble time, const size_t timeLevel,
                  const bool update = false);

    void ResidualEval(const NekDouble time, const size_t timeLevel,
                      const size_t n);

    void ResidualEval(const NekDouble time, const size_t timeLevel);

    void IntegratedResidualEval(const size_t timeLevel);

    using DriverParallelInTime::Interpolate;

    void Interpolate(const size_t coarseLevel, const SDCarray &in,
                     const size_t fineLevel, SDCarray &out, bool forced);

    void InterpolateSolution(const size_t timeLevel);

    void InterpolateResidual(const size_t timeLevel);

    void Restrict(const size_t fineLevel, const SDCarray &in,
                  const size_t coarseLevel, SDCarray &out);

    void RestrictSolution(const size_t timeLevel);

    void RestrictResidual(const size_t timeLevel);

    void ComputeFASCorrection(const size_t timeLevel);

    void Correct(const size_t coarseLevel,
                 const Array<OneD, Array<OneD, NekDouble>> &in,
                 const size_t fineLevel,
                 Array<OneD, Array<OneD, NekDouble>> &out, bool forced);

    void CorrectInitialSolution(const size_t timeLevel);

    void CorrectInitialResidual(const size_t timeLevel);

    void Correct(const size_t coarseLevel, const SDCarray &rest,
                 const SDCarray &in, const size_t fineLevel, SDCarray &out,
                 bool forced);

    void CorrectSolution(const size_t timeLevel);

    void CorrectResidual(const size_t timeLevel);

    void ApplyWindowing(void);

    void EvaluateSDCResidualNorm(const size_t timeLevel);

    void WriteOutput(const size_t step, const NekDouble time);

    // Storage of PFASST
    Array<OneD, size_t> m_QuadPts;
    Array<OneD, Array<OneD, NekDouble>> m_ImatFtoC;
    Array<OneD, Array<OneD, NekDouble>> m_ImatCtoF;
    Array<OneD, SDCarray> m_solutionRest;
    Array<OneD, SDCarray> m_residualRest;
    Array<OneD, SDCarray> m_integralRest;
    Array<OneD, SDCarray> m_correction;
    Array<OneD, SDCarray> m_storage;
    Array<OneD, std::shared_ptr<LibUtilities::TimeIntegrationSchemeSDC>>
        m_SDCSolver;

    bool m_updateResidual = false;
};

} // namespace Nektar::SolverUtils

#endif // NEKTAR_SOLVERUTILS_DRIVERPFASST_H
