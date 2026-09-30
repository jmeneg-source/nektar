///////////////////////////////////////////////////////////////////////////////
//
// File DriverParareal.cpp
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

#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/DriverParareal.h>
#include <boost/format.hpp>

namespace Nektar::SolverUtils
{
std::string DriverParareal::className =
    GetDriverFactory().RegisterCreatorFunction("Parareal",
                                               DriverParareal::create);
std::string DriverParareal::driverLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("Driver", "Parareal", 0);

/**
 *
 */
DriverParareal::DriverParareal(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverParallelInTime(pSession, pGraph)
{
}

/**
 *
 */
void DriverParareal::v_InitObject(std::ostream &out)
{
    DriverParallelInTime::v_InitObject(out);
    DriverParallelInTime::GetParametersFromSession();
    DriverParallelInTime::PrintSolverInfo(out);
    DriverParallelInTime::InitialiseEqSystem(false);
    AllocateMemory();
}

/**
 *
 */
void DriverParareal::v_Execute([[maybe_unused]] std::ostream &out)
{
    // Timing.
    Nektar::LibUtilities::Timer timer;
    NekDouble totalTime = 0.0, predictorTime = 0.0, coarseSolveTime = 0.0,
              fineSolveTime = 0.0, correctionTime = 0.0;

    // Get and assert parameters from session file.
    AssertParameters();

    // Initialie time step parameters.
    m_totalTime = m_timestep[m_fineLevel] * m_nsteps[m_fineLevel];
    m_chunkTime = m_totalTime / m_numChunks / m_numWindowsPIT;
    m_nsteps[m_fineLevel] /= m_numChunks * m_numWindowsPIT;
    m_nsteps[m_coarseLevel] /= m_numChunks * m_numWindowsPIT;

    // Pre-solve for one time-step to initialize preconditioner.
    UpdateSolution(m_coarseLevel, m_time0, 1, 0, 0);

    // Start iteration windows.
    m_comm->GetTimeComm()->Block();
    UpdateInitialConditionFromSolver(m_fineLevel);
    for (size_t w = 0; w < m_numWindowsPIT; w++)
    {
        timer.Start();
        // Initialize time for the current window.
        m_time = m_time0 + (w * m_numChunks) * m_chunkTime +
                 m_chunkRank * m_chunkTime;

        // Print window number.
        PrintHeader((boost::format("WINDOWS #%1%") % (w + 1)).str(), '*');

        // Update coarse initial condition.
        UpdateSolverInitialCondition(m_coarseLevel);

        // Run predictor.
        for (size_t i = 0; i < m_nVar; ++i)
        {
            RecvFromPreviousProc(m_EqSys[m_coarseLevel]->UpdatePhysField(i));
        }
        if (m_chunkRank > 0)
        {
            UpdateInitialConditionFromSolver(m_coarseLevel);
        }
        UpdateSolution(m_coarseLevel, m_time, m_nsteps[m_coarseLevel], w, 0);
        for (size_t i = 0; i < m_nVar; ++i)
        {
            SendToNextProc(m_EqSys[m_coarseLevel]->UpdatePhysField(i));
        }

        // Interpolate coarse solution.
        InterpolateCoarseSolution();

        // Compute exact solution, if necessary.
        if (m_exactSolution)
        {
            EvaluateExactSolution(m_fineLevel, m_time + m_chunkTime);
        }
        timer.Stop();
        predictorTime += timer.Elapsed().count();
        totalTime += timer.Elapsed().count();

        // Solution convergence monitoring.
        timer.Start();
        CopyToPhysField(m_fineLevel, m_coarseSolution);
        SolutionConvergenceMonitoring(m_fineLevel, 0);
        timer.Stop();
        totalTime += timer.Elapsed().count();
        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "Total Computation Time : " << totalTime << "s"
                      << std::endl
                      << std::flush;
        }

        // Start Parareal iteration.
        size_t iter         = 1;
        int convergenceCurr = false;
        int convergencePrev = (m_chunkRank == 0);
        while (iter <= m_iterMaxPIT && !convergenceCurr)
        {
            // Use previous parareal solution as "exact solution", if necessary.
            timer.Start();
            if (!m_exactSolution)
            {
                CopyFromPhysField(m_fineLevel, m_exactsoln);
            }
            timer.Stop();
            totalTime += timer.Elapsed().count();

            // Calculate fine solution (parallel-in-time).
            timer.Start();
            UpdateSolverInitialCondition(m_fineLevel);
            UpdateSolution(m_fineLevel, m_time, m_nsteps[m_fineLevel], w, iter);
            timer.Stop();
            fineSolveTime += timer.Elapsed().count();
            totalTime += timer.Elapsed().count();

            // Compute F -> F - Gold
            timer.Start();
            CorrectionWithOldCoarseSolution();
            timer.Stop();
            correctionTime += timer.Elapsed().count();
            totalTime += timer.Elapsed().count();

            // Receive coarse solution from previous processor.
            timer.Start();
            RecvFromPreviousProc(m_initialCondition, convergencePrev);
            timer.Stop();
            totalTime += timer.Elapsed().count();

            // Calculate coarse solution (serial-in-time).
            timer.Start();
            UpdateSolverInitialCondition(m_coarseLevel);
            UpdateSolution(m_coarseLevel, m_time, m_nsteps[m_coarseLevel], w,
                           iter);
            timer.Stop();
            coarseSolveTime += timer.Elapsed().count();
            totalTime += timer.Elapsed().count();

            // Compute F -> F + Gnew
            timer.Start();
            CorrectionWithNewCoarseSolution();
            timer.Stop();
            correctionTime += timer.Elapsed().count();
            totalTime += timer.Elapsed().count();

            // Solution convergence monitoring.
            CopyToPhysField(m_fineLevel, m_fineSolution);
            SolutionConvergenceMonitoring(m_fineLevel, iter);
            if (m_chunkRank == m_numChunks - 1 &&
                m_comm->GetSpaceComm()->GetRank() == 0)
            {
                std::cout << "Total Computation Time : " << totalTime << "s"
                          << std::endl
                          << std::flush;
                std::cout << " - Predictor Time : " << predictorTime << "s"
                          << std::endl
                          << std::flush;
                std::cout << " - Coarse Solve Time : " << coarseSolveTime << "s"
                          << std::endl
                          << std::flush;
                std::cout << " - Fine Solve Time : " << fineSolveTime << "s"
                          << std::endl
                          << std::flush;
                std::cout << " - Correction Time : " << correctionTime << "s"
                          << std::endl
                          << std::flush;
            }

            // Check convergence of L2 error for each time chunk.
            convergenceCurr = (vL2ErrorMax() < m_tolerPIT && convergencePrev) ||
                              (m_chunkRank + 1 == iter);

            // Send solution to next processor.
            timer.Start();
            SendToNextProc(m_fineSolution, convergenceCurr);
            timer.Stop();
            totalTime += timer.Elapsed().count();

            // Increment iteration index.
            iter++;
        }

        // Copy converged check points.
        CopyConvergedCheckPoints(w, iter);

        // Write time chunk solution to files.
        WriteTimeChunkOuput();

        // Apply windowing.
        timer.Start();
        ApplyWindowing(w);
        timer.Stop();
        totalTime += timer.Elapsed().count();
    }

    m_comm->GetTimeComm()->Block();
    PrintHeader("SUMMARY", '*');
    EvaluateExactSolution(m_fineLevel, m_time + m_chunkTime);
    SolutionConvergenceSummary(m_fineLevel);
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "Total Computation Time : " << totalTime << "s"
                  << std::endl
                  << std::flush;
        std::cout << " - Predictor Time : " << predictorTime << "s" << std::endl
                  << std::flush;
        std::cout << " - Coarse Solve Time : " << coarseSolveTime << "s"
                  << std::endl
                  << std::flush;
        std::cout << " - Fine Solve Time : " << fineSolveTime << "s"
                  << std::endl
                  << std::flush;
        std::cout << " - Correction Time : " << correctionTime << "s"
                  << std::endl
                  << std::flush;
    }
}

/**
 *
 */
void DriverParareal::AllocateMemory(void)
{
    // Allocate storage for Parareal solver.
    m_initialCondition = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    m_fineSolution     = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    m_coarseSolution   = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        m_initialCondition[i] =
            Array<OneD, NekDouble>(m_npts[m_fineLevel], 0.0);
        m_fineSolution[i] = m_EqSys[m_fineLevel]->UpdatePhysField(i);
        m_coarseSolution[i] =
            (m_npts[m_fineLevel] == m_npts[m_coarseLevel])
                ? m_EqSys[m_coarseLevel]->UpdatePhysField(i)
                : Array<OneD, NekDouble>(m_npts[m_fineLevel], 0.0);
    }
}

/**
 *
 */
void DriverParareal::AssertParameters(void)
{
    // Assert time-stepping parameters
    ASSERTL0(
        m_nsteps[m_fineLevel] % m_numChunks == 0,
        "Total number of fine step should be divisible by number of chunks.");

    ASSERTL0(
        m_nsteps[m_coarseLevel] % m_numChunks == 0,
        "Total number of coarse step should be divisible by number of chunks.");

    ASSERTL0(m_nsteps[m_fineLevel] % (m_numChunks * m_numWindowsPIT) == 0,
             "Total number of fine step should be divisible by number of "
             "windows times number of chunks.");

    ASSERTL0(m_nsteps[m_coarseLevel] % (m_numChunks * m_numWindowsPIT) == 0,
             "Total number of coarse step should be divisible by number of "
             "windows times number of chunks.");

    ASSERTL0(fabs(m_timestep[m_coarseLevel] * m_nsteps[m_coarseLevel] -
                  m_timestep[m_fineLevel] * m_nsteps[m_fineLevel]) < 1e-12,
             "Fine and coarse total computational times do not match");

    ASSERTL0(m_EqSys[m_fineLevel]
                     ->GetTimeIntegrationScheme()
                     ->GetNumIntegrationPhases() == 1,
             "Only single step time-integration schemes currently supported "
             "for Parareal");

    ASSERTL0(m_EqSys[m_coarseLevel]
                     ->GetTimeIntegrationScheme()
                     ->GetNumIntegrationPhases() == 1,
             "Only single step time-integration schemes currently supported "
             "for Parareal");

    // Assert I/O parameters
    if (m_EqSys[m_fineLevel]->GetInfoSteps())
    {
        ASSERTL0(m_nsteps[m_fineLevel] % (m_EqSys[m_fineLevel]->GetInfoSteps() *
                                          m_numChunks * m_numWindowsPIT) ==
                     0,
                 "number of IO_InfoSteps should divide number of fine steps "
                 "per time chunk");
    }

    if (m_EqSys[m_coarseLevel]->GetInfoSteps())
    {
        ASSERTL0(m_nsteps[m_coarseLevel] %
                         (m_EqSys[m_coarseLevel]->GetInfoSteps() * m_numChunks *
                          m_numWindowsPIT) ==
                     0,
                 "number of IO_InfoSteps should divide number of coarse steps "
                 "per time chunk");
    }

    if (m_EqSys[m_fineLevel]->GetCheckpointSteps())
    {
        ASSERTL0(m_nsteps[m_fineLevel] %
                         (m_EqSys[m_fineLevel]->GetCheckpointSteps() *
                          m_numChunks * m_numWindowsPIT) ==
                     0,
                 "number of IO_CheckSteps should divide number of fine steps "
                 "per time chunk");
    }

    if (m_EqSys[m_coarseLevel]->GetCheckpointSteps())
    {
        ASSERTL0(m_nsteps[m_coarseLevel] %
                         (m_EqSys[m_coarseLevel]->GetCheckpointSteps() *
                          m_numChunks * m_numWindowsPIT) ==
                     0,
                 "number of IO_CheckSteps should divide number of coarse steps "
                 "per time chunk");
    }
}

/**
 *
 */
void DriverParareal::UpdateInitialConditionFromSolver(const size_t timeLevel)
{
    // Interpolate solution to fine field.
    Interpolate(m_EqSys[timeLevel]->UpdateFields(),
                m_EqSys[m_fineLevel]->UpdateFields(), NullNekDoubleArrayOfArray,
                m_initialCondition);
}

/**
 *
 */
void DriverParareal::UpdateSolverInitialCondition(const size_t timeLevel)
{
    // Restrict fine field to coarse solution.
    Interpolate(m_EqSys[m_fineLevel]->UpdateFields(),
                m_EqSys[timeLevel]->UpdateFields(), m_initialCondition,
                NullNekDoubleArrayOfArray);
}

/**
 *
 */
void DriverParareal::UpdateSolution(const size_t timeLevel,
                                    const NekDouble time, const size_t nstep,
                                    const size_t wd, const size_t iter)
{
    // Number of checkpoint by chunk.
    size_t nChkPts =
        m_EqSys[timeLevel]->GetCheckpointSteps()
            ? m_nsteps[timeLevel] / m_EqSys[timeLevel]->GetCheckpointSteps()
            : 1;

    // Checkpoint index.
    size_t iChkPts = (m_chunkRank + wd * m_numChunks) * nChkPts + 1;

    // Reinitialize check point number for each parallel-in-time
    // iteration.
    m_EqSys[timeLevel]->SetCheckpointNumber(iChkPts);

    // Update parallel-in-time iteration number.
    m_EqSys[timeLevel]->SetIterationNumberPIT(iter);

    // Update parallel-in-time window number.
    m_EqSys[timeLevel]->SetWindowNumberPIT(wd);

    m_EqSys[timeLevel]->SetTime(time);
    m_EqSys[timeLevel]->SetSteps(nstep);
    m_EqSys[timeLevel]->DoSolve();
}

/**
 *
 */
void DriverParareal::CorrectionWithOldCoarseSolution(void)
{
    // Correct solution F -> F - Gold.
    for (size_t i = 0; i < m_nVar; ++i)
    {
        Vmath::Vsub(m_npts[m_fineLevel], m_fineSolution[i], 1,
                    m_coarseSolution[i], 1, m_fineSolution[i], 1);
    }
}

/**
 *
 */
void DriverParareal::CorrectionWithNewCoarseSolution(void)
{
    // Interpolate coarse solution.
    InterpolateCoarseSolution();

    // Correct solution F -> F + Gnew.
    for (size_t i = 0; i < m_nVar; ++i)
    {
        Vmath::Vadd(m_npts[m_fineLevel], m_fineSolution[i], 1,
                    m_coarseSolution[i], 1, m_fineSolution[i], 1);
    }
}

/**
 *
 */
void DriverParareal::InterpolateCoarseSolution(void)
{
    if (m_npts[m_fineLevel] != m_npts[m_coarseLevel])
    {
        // Interpolate coarse solution to fine field.
        Interpolate(m_EqSys[m_coarseLevel]->UpdateFields(),
                    m_EqSys[m_fineLevel]->UpdateFields(),
                    NullNekDoubleArrayOfArray, m_coarseSolution);
    }
}

/**
 *
 */
void DriverParareal::ApplyWindowing(const size_t w)
{
    if (w == m_numWindowsPIT - 1)
    {
        // No windowing required for the last window.
        return;
    }

    // Use last chunk solution as initial condition for the next
    // window.
    if (m_chunkRank == m_numChunks - 1)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Vcopy(m_npts[m_fineLevel],
                         m_EqSys[m_fineLevel]->UpdatePhysField(i), 1,
                         m_initialCondition[i], 1);
        }
    }

    // Broadcast I.C. for windowing.
    for (size_t i = 0; i < m_nVar; ++i)
    {
        m_comm->GetTimeComm()->Bcast(m_initialCondition[i], m_numChunks - 1);
    }
}

/**
 *
 */
void DriverParareal::CopyConvergedCheckPoints(const size_t w, const size_t k)
{
    // Determine max number of iteration.
    size_t kmax = k;
    m_comm->GetTimeComm()->AllReduce(kmax, Nektar::LibUtilities::ReduceMax);

    if (m_comm->GetSpaceComm()->GetRank() == 0)
    {
        for (size_t j = k; j < kmax; j++)
        {
            // Copy converged solution files from directory corresponding to
            // iteration j - 1 to the directory corresponding to iteration j.

            auto sessionName = m_EqSys[m_fineLevel]->GetSessionName();

            // Input directory name.
            std::string indir =
                sessionName + "_" + std::to_string(j - 1) + ".pit";

            /// Output directory name.
            std::string outdir = sessionName + "_" + std::to_string(j) + ".pit";

            for (size_t timeLevel = 0; timeLevel < m_nTimeLevel; timeLevel++)
            {
                // Number of checkpoint by chunk.
                size_t nChkPts =
                    m_EqSys[timeLevel]->GetCheckpointSteps()
                        ? m_nsteps[timeLevel] /
                              m_EqSys[timeLevel]->GetCheckpointSteps()
                        : 0;

                // Checkpoint index.
                size_t iChkPts = (m_chunkRank + w * m_numChunks) * nChkPts;

                for (size_t i = 1; i <= nChkPts; i++)
                {
                    // Filename corresponding to checkpoint iChkPts.
                    std::string filename = sessionName + "_timeLevel" +
                                           std::to_string(timeLevel) + "_" +
                                           std::to_string(iChkPts + i) + ".chk";

                    // Intput full file name.
                    std::string infullname = indir + "/" + filename;

                    // Output full file name.
                    std::string outfullname = outdir + "/" + filename;

                    // Remove output file if already existing.
                    fs::remove_all(outfullname);

                    // Copy converged solution files.
                    fs::copy(infullname, outfullname);
                }
            }
        }
    }
}

/**
 *
 */
void DriverParareal::WriteTimeChunkOuput(void)
{
    PrintHeader("PRINT SOLUTION FILES", '-');

    // Update field coefficients.
    UpdateFieldCoeffs(m_fineLevel);

    // Output solution files.
    m_EqSys[m_fineLevel]->Output();
}

} // namespace Nektar::SolverUtils
