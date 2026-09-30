///////////////////////////////////////////////////////////////////////////////
//
// File FileSolution.h
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
// Description: load discrete check-point files and interpolate them into a
// continuous field
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_FILESOLUTION_H
#define NEKTAR_SOLVERS_FILESOLUTION_H

#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/FFT/NektarFFT.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/UnsteadySystem.h>
#include <complex>

namespace Nektar::SolverUtils
{

class FileFieldInterpolator;
typedef std::shared_ptr<FileFieldInterpolator> FileFieldInterpolatorSharedPtr;
class FileFieldInterpolator
{
public:
    friend class MemoryManager<FileFieldInterpolator>;

    void InitObject(
        const std::string functionName,
        LibUtilities::SessionReaderSharedPtr pSession,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> pFields,
        std::set<std::string> &variables, std::map<std::string, int> &series,
        std::map<std::string, NekDouble> &time);

    void InterpolateField(const std::string variable,
                          Array<OneD, NekDouble> &outarray, NekDouble time);

    void InterpolateField(const int v, Array<OneD, NekDouble> &outarray,
                          const NekDouble time);

    NekDouble GetStartTime();

protected:
    LibUtilities::SessionReaderSharedPtr m_session;
    /// number of slices
    int m_start;
    int m_skip;
    int m_slices;
    /// period length
    NekDouble m_timeStart;
    NekDouble m_period;
    int m_isperiodic;
    int m_interporder;
    /// interpolation vector
    std::map<int, Array<OneD, NekDouble>> m_interp;
    /// variables
    std::map<std::string, int> m_variableMap;

    DNekBlkMatSharedPtr GetFloquetBlockMatrix(int nexp);

    FileFieldInterpolator() = default;

    ~FileFieldInterpolator() = default;

    void DFT(const std::string file,
             const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
             const bool timefromfile);

    /// Import Base flow
    void ImportFldBase(
        std::string pInfile,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        int slice, std::map<std::string, NekDouble> &params);
};

/**
 * \brief This class is the base class for Navier Stokes problems
 *
 */
class FileSolution : public SolverUtils::AdvectionSystem,
                     public SolverUtils::FluidInterface
{
public:
    friend class MemoryManager<FileSolution>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<FileSolution>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

protected:
    FileSolution(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const SpatialDomains::MeshGraphSharedPtr &pGraph);

    ~FileSolution() override = default;

    /// Initialise the object
    void v_InitObject(bool DeclareField = true) override;

    void v_GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &velocity) override;

    void v_GetPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure) override;

    using SolverUtils::EquationSystem::v_GetPressure;

    void v_GetDensity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &density) override;

    bool v_HasConstantDensity() override;

    /// Compute the RHS
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Compute the projection
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, NekDouble time,
        NekDouble lambda);

    bool v_PostIntegrate(int step) override;

    bool v_RequireFwdTrans() override;

    void v_DoInitialise(bool dumpInitialConditions) override;

    void UpdateField(NekDouble time);

private:
    FileFieldInterpolatorSharedPtr m_solutionFile;
    std::set<std::string> m_variableFile;
    Array<OneD, Array<OneD, NekDouble>> m_coord;
    std::map<std::string, LibUtilities::EquationSharedPtr> m_solutionFunction;
};

} // namespace Nektar::SolverUtils

#endif // NEKTAR_SOLVERS_INCNAVIERSTOKES_H
