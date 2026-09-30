///////////////////////////////////////////////////////////////////////////////
//
// File FileSolution.cpp
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

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <tinyxml.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Communication/Comm.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <SolverUtils/Filters/Filter.h>

#include <LibUtilities/BasicUtils/Filesystem.hpp>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/Polylib/Polylib.h>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <ReviewSolution/EquationSystems/FileSolution.h>
#include <SolverUtils/Forcing/ForcingMovingReferenceFrame.h>
#include <StdRegions/StdSegExp.h>

namespace Nektar::SolverUtils
{
/**
 * Constructor. Creates ...
 *
 * \param
 * \param
 */
std::string FileSolution::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "FileSolution", FileSolution::create,
        "review a solution from check point files.");

FileSolution::FileSolution(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph)
{
    m_solutionFile = MemoryManager<FileFieldInterpolator>::AllocateSharedPtr();
}

/**
 * @brief Initialisation object for the unsteady linear advection equation.
 */
void FileSolution::v_InitObject(bool DeclareField)
{
    // Call to the initialisation object of UnsteadySystem
    AdvectionSystem::v_InitObject(DeclareField);

    // If explicit it computes RHS and PROJECTION for the time integration
    m_ode.DefineOdeRhs(&FileSolution::DoOdeRhs, this);
    m_ode.DefineProjection(&FileSolution::DoOdeProjection, this);
    m_ode.DefineImplicitSolve(&FileSolution::DoImplicitSolve, this);

    // load solution
    for (size_t i = 0; i < m_session->GetVariables().size(); ++i)
    {
        if (m_session->GetFunctionType("Solution", m_session->GetVariable(i)) ==
            LibUtilities::eFunctionTypeFile)
        {
            m_variableFile.insert(m_session->GetVariable(i));
        }
        else if (m_session->GetFunctionType("Solution",
                                            m_session->GetVariable(i)) ==
                 LibUtilities::eFunctionTypeExpression)
        {
            m_solutionFunction[m_session->GetVariable(i)] =
                m_session->GetFunction("Solution", m_session->GetVariable(i));
        }
        else
        {
            ASSERTL0(false, "solution not defined for variable " +
                                m_session->GetVariable(i));
        }
    }

    if (m_variableFile.size())
    {
        std::map<std::string, int> series; // start, skip, slices, order
        int tmp;
        if (m_session->DefinesParameter("N_slices"))
        {
            m_session->LoadParameter("N_slices", tmp, 1);
            series["slices"] = tmp;
        }
        if (m_session->DefinesParameter("N_start"))
        {
            m_session->LoadParameter("N_start", tmp, 0);
            series["start"] = tmp;
        }
        if (m_session->DefinesParameter("N_skip"))
        {
            m_session->LoadParameter("N_skip", tmp, 1);
            series["skip"] = tmp;
        }
        if (m_session->DefinesParameter("BaseFlow_interporder"))
        {
            m_session->LoadParameter("BaseFlow_interporder", tmp, 1);
            series["order"] = tmp;
        }
        if (m_session->DefinesParameter("Is_periodic"))
        {
            m_session->LoadParameter("Is_periodic", tmp, 1);
            series["isperiodic"] = tmp;
        }

        std::map<std::string, NekDouble> times;
        NekDouble dtmp;
        if (m_session->DefinesParameter("time_start"))
        {
            m_session->LoadParameter("time_start", dtmp, 0.);
            times["start"] = dtmp;
        }
        if (m_session->DefinesParameter("period"))
        {
            m_session->LoadParameter("period", dtmp, 1.);
            times["period"] = dtmp;
        }
        m_solutionFile->InitObject("Solution", m_session, m_fields,
                                   m_variableFile, series, times);
    }

    if (m_solutionFunction.size())
    {
        m_coord    = Array<OneD, Array<OneD, NekDouble>>(3);
        int nq     = m_fields[0]->GetNpoints();
        m_coord[0] = Array<OneD, NekDouble>(nq);
        m_coord[1] = Array<OneD, NekDouble>(nq);
        m_coord[2] = Array<OneD, NekDouble>(nq);
        m_fields[0]->GetCoords(m_coord[0], m_coord[1], m_coord[2]);
    }
}

void FileSolution::DoImplicitSolve(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] NekDouble time, [[maybe_unused]] NekDouble lambda)
{
}

/**
 * @brief Compute the right-hand side for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void FileSolution::DoOdeRhs(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time)
{
    int nSolutionPts = GetNpoints();
    for (size_t i = 0; i < outarray.size(); ++i)
    {
        Vmath::Zero(nSolutionPts, outarray[i], 1);
    }
}

/**
 * @brief Compute the projection for the linear advection equation.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void FileSolution::DoOdeProjection(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time)
{
}

bool FileSolution::v_PostIntegrate([[maybe_unused]] int step)
{
    UpdateField(m_time);
    return false;
}

void FileSolution::v_DoInitialise(bool dumpInitialConditions)
{
    m_time = m_solutionFile->GetStartTime();
    SetBoundaryConditions(m_time);
    UpdateField(m_time);

    // Update time in field info if required
    if (m_fieldMetaDataMap.find("Time") != m_fieldMetaDataMap.end())
    {
        m_fieldMetaDataMap["Time"] = boost::lexical_cast<std::string>(m_time);
    }

    if (dumpInitialConditions && m_checksteps && m_nchk == 0)
    {
        Checkpoint_Output(m_nchk);
    }
    ++m_nchk;
}

void FileSolution::UpdateField(NekDouble time)
{
    for (int i = 0; i < m_fields.size(); ++i)
    {
        std::string var = m_session->GetVariable(i);
        if (m_solutionFunction.count(var))
        {
            m_solutionFunction[var]->Evaluate(m_coord[0], m_coord[1],
                                              m_coord[2], time,
                                              m_fields[i]->UpdatePhys());
            m_fields[i]->SetPhysState(true);
            bool wavespace = m_fields[i]->GetWaveSpace();
            m_fields[i]->SetWaveSpace(false);
            m_fields[i]->FwdTransBndConstrained(m_fields[i]->GetPhys(),
                                                m_fields[i]->UpdateCoeffs());
            m_fields[i]->SetWaveSpace(wavespace);
        }
        else if (m_variableFile.count(var))
        {
            m_solutionFile->InterpolateField(var, m_fields[i]->UpdateCoeffs(),
                                             time);
            m_fields[i]->SetPhysState(true);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
        }
        else
        {
            ASSERTL0(false, "solution not defined for variable " + var);
        }
    }
}

bool FileSolution::v_RequireFwdTrans()
{
    return false;
}

void FileSolution::v_GetPressure(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    [[maybe_unused]] Array<OneD, NekDouble> &pressure)
{
}

void FileSolution::v_GetDensity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &density)
{
    for (size_t i = 0; i < m_session->GetVariables().size(); ++i)
    {
        if (m_session->GetVariable(i) == "rho")
        {
            int npoints = m_fields[i]->GetNpoints();
            Vmath::Vcopy(npoints, physfield[i], 1, density, 1);
        }
    }
    Vmath::Fill(m_fields[0]->GetNpoints(), 1., density, 1);
}

bool FileSolution::v_HasConstantDensity()
{
    for (size_t i = 0; i < m_session->GetVariables().size(); ++i)
    {
        if (m_session->GetVariable(i) == "rho")
        {
            return false;
        }
    }
    return true;
}

void FileSolution::v_GetVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    int npoints = m_fields[0]->GetNpoints();
    if (boost::iequals(m_session->GetVariable(0), "u"))
    {
        // IncNavierStokesSolver
        for (int i = 0; i < velocity.size(); ++i)
        {
            Vmath::Vcopy(npoints, physfield[i], 1, velocity[i], 1);
        }
    }
    else if (boost::iequals(m_session->GetVariable(0), "rho") &&
             boost::iequals(m_session->GetVariable(1), "rhou"))
    {
        // CompressibleFlowSolver
        for (int i = 0; i < velocity.size(); ++i)
        {
            Vmath::Vdiv(npoints, physfield[i], 1, physfield[0], 1, velocity[i],
                        1);
        }
    }
    else
    {
        // Unknown
        ASSERTL0(false, "Could not identify velocity for ProcessVorticity");
    }
}

void FileFieldInterpolator::InitObject(
    const std::string functionName,
    LibUtilities::SessionReaderSharedPtr pSession,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> pFields,
    std::set<std::string> &variables, std::map<std::string, int> &series,
    std::map<std::string, NekDouble> &times)
{
    m_session = pSession;
    for (size_t i = 0; i < m_session->GetVariables().size(); ++i)
    {
        std::string var = m_session->GetVariable(i);
        if (variables.count(var))
        {
            ASSERTL0(m_session->DefinesFunction(functionName, var) &&
                         (m_session->GetFunctionType(functionName, var) ==
                          LibUtilities::eFunctionTypeFile),
                     functionName + "(" + var + ") is not defined as a file.");
            m_variableMap[var] = i;
        }
    }

    if (series.count("start"))
    {
        m_start = series["start"];
    }
    else
    {
        m_start = 0;
    }

    if (series.count("skip"))
    {
        m_skip = series["skip"];
    }
    else
    {
        m_skip = 1;
    }

    if (series.count("slices"))
    {
        m_slices = series["slices"];
    }
    else
    {
        m_slices = 1;
    }

    if (series.count("order"))
    {
        m_interporder = series["order"];
    }
    else
    {
        m_interporder = 0;
    }

    if (series.count("isperiodic"))
    {
        m_isperiodic = series["isperiodic"];
    }
    else
    {
        m_isperiodic = 1;
    }

    bool timefromfile = false;
    if (times.count("period"))
    {
        m_period = times["period"];
    }
    else if (m_slices > 1)
    {
        timefromfile = true;
    }
    if (times.count("start"))
    {
        m_timeStart = times["start"];
    }
    else
    {
        m_timeStart = 0.;
    }

    if (m_session->GetComm()->GetRank() == 0)
    {
        std::cout << "baseflow info : interpolation order " << m_interporder
                  << ", period " << m_period << ", periodicity ";
        if (m_isperiodic)
        {
            std::cout << "yes\n";
        }
        else
        {
            std::cout << "no\n";
        }
        std::cout << "baseflow info : files from " << m_start << " to "
                  << (m_start + (m_slices - 1) * m_skip) << " (skip " << m_skip
                  << ") with " << (m_slices - (m_interporder > 1))
                  << " time intervals" << std::endl;
    }

    std::string file = m_session->GetFunctionFilename("Solution", 0);
    DFT(file, pFields, timefromfile);
}

/**
 * Import field from infile and load into \a m_fields. This routine will
 * also perform a \a BwdTrans to ensure data is in both the physical and
 * coefficient storage.
 * @param   pInFile          Filename to read.
 * @param   pFields          Array of expansion lists
 */
void FileFieldInterpolator::ImportFldBase(
    std::string pInfile,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    int pSlice, std::map<std::string, NekDouble> &params)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::vector<NekDouble>> FieldData;

    int numexp = pFields[0]->GetExpSize();
    Array<OneD, int> ElementGIDs(numexp);

    // Define list of global element ids
    for (int i = 0; i < numexp; ++i)
    {
        ElementGIDs[i] = pFields[0]->GetExp(i)->GetGeom()->GetGlobalID();
    }

    // Get Homogeneous
    LibUtilities::FieldIOSharedPtr fld =
        LibUtilities::FieldIO::CreateForFile(m_session, pInfile);
    fld->Import(pInfile, FieldDef, FieldData,
                LibUtilities::NullFieldMetaDataMap, ElementGIDs);

    int nFileVar = FieldDef[0]->m_fields.size();
    for (int j = 0; j < nFileVar; ++j)
    {
        std::string var = FieldDef[0]->m_fields[j];
        if (!m_variableMap.count(var))
        {
            continue;
        }
        int ncoeffs = pFields[m_variableMap[var]]->GetNcoeffs();
        Array<OneD, NekDouble> tmp_coeff(ncoeffs, 0.);
        for (int i = 0; i < FieldDef.size(); ++i)
        {
            pFields[m_variableMap[var]]->ExtractDataToCoeffs(
                FieldDef[i], FieldData[i], FieldDef[i]->m_fields[j], tmp_coeff);
        }
        Vmath::Vcopy(ncoeffs, &tmp_coeff[0], 1,
                     &m_interp[m_variableMap[var]][pSlice * ncoeffs], 1);
    }

    LibUtilities::FieldMetaDataMap fieldMetaDataMap;
    fld->ImportFieldMetaData(pInfile, fieldMetaDataMap);
    // check to see if time defined
    if (fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
    {
        auto iter = fieldMetaDataMap.find("Time");
        if (iter != fieldMetaDataMap.end())
        {
            params["time"] = std::stod(iter->second);
        }
    }
}

void FileFieldInterpolator::InterpolateField(const std::string variable,
                                             Array<OneD, NekDouble> &outarray,
                                             const NekDouble time)
{
    if (!m_variableMap.count(variable))
    {
        return;
    }
    InterpolateField(m_variableMap[variable], outarray, time);
}

void FileFieldInterpolator::InterpolateField(const int v,
                                             Array<OneD, NekDouble> &outarray,
                                             NekDouble time)
{
    // doesnot have this variable
    if (!m_interp.count(v))
    {
        return;
    }
    // one slice, steady solution
    int npoints = m_interp[v].size() / m_slices;
    if (m_slices == 1)
    {
        Vmath::Vcopy(npoints, &m_interp[v][0], 1, &outarray[0], 1);
        return;
    }
    // unsteady solution
    time -= m_timeStart;
    if (m_isperiodic && time > m_period)
    {
        time = fmod(time, m_period);
        if (time < 0.)
        {
            time += m_period;
        }
    }
    if (m_interporder < 1)
    {
        NekDouble BetaT = 2 * M_PI * fmod(time, m_period) / m_period;
        NekDouble phase;
        Array<OneD, NekDouble> auxiliary(npoints);

        Vmath::Vcopy(npoints, &m_interp[v][0], 1, &outarray[0], 1);
        Vmath::Svtvp(npoints, cos(0.5 * m_slices * BetaT),
                     &m_interp[v][npoints], 1, &outarray[0], 1, &outarray[0],
                     1);

        for (int i = 2; i < m_slices; i += 2)
        {
            phase = (i >> 1) * BetaT;

            Vmath::Svtvp(npoints, cos(phase), &m_interp[v][i * npoints], 1,
                         &outarray[0], 1, &outarray[0], 1);
            Vmath::Svtvp(npoints, -sin(phase), &m_interp[v][(i + 1) * npoints],
                         1, &outarray[0], 1, &outarray[0], 1);
        }
    }
    else
    {
        NekDouble x = time;
        x           = x / m_period * (m_slices - 1);
        int ix      = x;
        if (ix < 0)
        {
            ix = 0;
        }
        if (ix > m_slices - 2)
        {
            ix = m_slices - 2;
        }
        int padleft = std::max(0, m_interporder / 2 - 1);
        if (padleft > ix)
        {
            padleft = ix;
        }
        int padright = m_interporder - 1 - padleft;
        if (padright > m_slices - 1 - ix)
        {
            padright = m_slices - 1 - ix;
        }
        padleft = m_interporder - 1 - padright;
        Array<OneD, NekDouble> coeff(m_interporder, 1.);
        for (int i = 0; i < m_interporder; ++i)
        {
            for (int j = 0; j < m_interporder; ++j)
            {
                if (i != j)
                {
                    coeff[i] *= (x - ix + padleft - (NekDouble)j) /
                                ((NekDouble)i - (NekDouble)j);
                }
            }
        }
        Vmath::Zero(npoints, &outarray[0], 1);
        for (int i = ix - padleft; i < ix + padright + 1; ++i)
        {
            Vmath::Svtvp(npoints, coeff[i - ix + padleft],
                         &m_interp[v][i * npoints], 1, &outarray[0], 1,
                         &outarray[0], 1);
        }
    }
}

DNekBlkMatSharedPtr FileFieldInterpolator::GetFloquetBlockMatrix(int nexp)
{
    DNekMatSharedPtr loc_mat;
    DNekBlkMatSharedPtr BlkMatrix;

    Array<OneD, unsigned int> nrows(nexp);
    Array<OneD, unsigned int> ncols(nexp);

    nrows = Array<OneD, unsigned int>(nexp, m_slices);
    ncols = Array<OneD, unsigned int>(nexp, m_slices);

    MatrixStorage blkmatStorage = eDIAGONAL;
    BlkMatrix = MemoryManager<DNekBlkMat>::AllocateSharedPtr(nrows, ncols,
                                                             blkmatStorage);

    const LibUtilities::PointsKey Pkey(m_slices,
                                       LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey BK(LibUtilities::eFourier, m_slices, Pkey);
    StdRegions::StdSegExp StdSeg(BK);

    StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                    StdSeg.DetShapeType(), StdSeg);

    loc_mat = StdSeg.GetStdMatrix(matkey);

    // set up array of block matrices.
    for (int i = 0; i < nexp; ++i)
    {
        BlkMatrix->SetBlock(i, i, loc_mat);
    }

    return BlkMatrix;
}

// Discrete Fourier Transform for Floquet analysis
void FileFieldInterpolator::DFT(
    const std::string file,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const bool timefromfile)
{
    for (auto it : m_variableMap)
    {
        int ncoeffs         = pFields[it.second]->GetNcoeffs();
        m_interp[it.second] = Array<OneD, NekDouble>(ncoeffs * m_slices, 0.0);
    }

    // Import the slides into auxiliary vector
    // The base flow should be stored in the form "filename_%d.ext"
    // A subdirectory can also be included, such as "dir/filename_%d.ext"
    size_t found = file.find("%d");
    std::map<std::string, NekDouble> params;
    if (found != std::string::npos)
    {
        ASSERTL0(file.find("%d", found + 1) == std::string::npos,
                 "There are more than one '%d'.");
        int nstart = m_start;
        std::vector<NekDouble> times;
        for (int i = 0; i < m_slices; ++i)
        {
            int filen = nstart + i * m_skip;
            auto fmt  = boost::format(file) % filen;
            ImportFldBase(fmt.str(), pFields, i, params);
            if (m_session->GetComm()->GetRank() == 0)
            {
                std::cout << "read base flow file " << fmt.str() << std::endl;
            }
            if (timefromfile && params.count("time"))
            {
                times.push_back(params["time"]);
            }
        }
        if (timefromfile && times.size() == m_slices && m_slices > 1)
        {
            m_timeStart = times[0];
            if (m_interporder < 1)
            {
                m_period = m_slices * (times[m_slices - 1] - times[0]) /
                           (m_slices - 1.);
            }
            else
            {
                m_period = times[m_slices - 1] - times[0];
            }
        }
    }
    else if (m_slices == 1)
    {
        ImportFldBase(file.c_str(), pFields, 0, params);
    }
    else
    {
        ASSERTL0(
            false,
            "Since N_slices is specified, the filename provided for function "
            "'BaseFlow' must include exactly one instance of the format "
            "specifier '%d', to index the time-slices.");
    }

    if (!m_isperiodic || m_slices == 1)
    {
        return;
    }

    // Discrete Fourier Transform of the fields
    for (auto it : m_interp)
    {
        int npoints = pFields[it.first]->GetNcoeffs();
#ifdef NEKTAR_USING_FFTW

        // Discrete Fourier Transform using FFTW
        Array<OneD, NekDouble> fft_in(npoints * m_slices);
        Array<OneD, NekDouble> fft_out(npoints * m_slices);

        Array<OneD, NekDouble> m_tmpIN(m_slices);
        Array<OneD, NekDouble> m_tmpOUT(m_slices);

        // Shuffle the data
        for (int j = 0; j < m_slices; ++j)
        {
            Vmath::Vcopy(npoints, &(it.second)[j * npoints], 1, &(fft_in[j]),
                         m_slices);
        }

        m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW",
                                                                   m_slices);

        // FFT Transform
        for (int i = 0; i < npoints; i++)
        {
            m_FFT->FFTFwdTrans(m_tmpIN  = fft_in + i * m_slices,
                               m_tmpOUT = fft_out + i * m_slices);
        }

        // Reshuffle data
        for (int s = 0; s < m_slices; ++s)
        {
            Vmath::Vcopy(npoints, &fft_out[s], m_slices,
                         &(it.second)[s * npoints], 1);
        }

        Vmath::Zero(fft_in.size(), &fft_in[0], 1);
        Vmath::Zero(fft_out.size(), &fft_out[0], 1);
#else
        // Discrete Fourier Transform using MVM
        DNekBlkMatSharedPtr blkmat;
        blkmat = GetFloquetBlockMatrix(npoints);

        int nrows = blkmat->GetRows();
        int ncols = blkmat->GetColumns();

        Array<OneD, NekDouble> sortedinarray(ncols);
        Array<OneD, NekDouble> sortedoutarray(nrows);

        // Shuffle the data
        for (int j = 0; j < m_slices; ++j)
        {
            Vmath::Vcopy(npoints, &(it.second)[j * npoints], 1,
                         &(sortedinarray[j]), m_slices);
        }

        // Create NekVectors from the given data arrays
        NekVector<NekDouble> in(ncols, sortedinarray, eWrapper);
        NekVector<NekDouble> out(nrows, sortedoutarray, eWrapper);

        // Perform matrix-vector multiply.
        out = (*blkmat) * in;

        // Reshuffle data
        for (int s = 0; s < m_slices; ++s)
        {
            Vmath::Vcopy(npoints, &sortedoutarray[s], m_slices,
                         &(it.second)[s * npoints], 1);
        }

        for (int r = 0; r < sortedinarray.size(); ++r)
        {
            sortedinarray[0]  = 0;
            sortedoutarray[0] = 0;
        }

#endif
    }
}

NekDouble FileFieldInterpolator::GetStartTime()
{
    return m_timeStart;
}
} // namespace Nektar::SolverUtils
