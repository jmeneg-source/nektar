///////////////////////////////////////////////////////////////////////////////
//
// File: FilterLagrangianPoints.cpp
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
// Description: Evaluate physics values at a set of points.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
using namespace std;

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <SolverUtils/Filters/FilterLagrangianPoints.h>
#include <boost/format.hpp>

namespace Nektar::SolverUtils
{

std::string FilterLagrangianPoints::className =
    GetFilterFactory().RegisterCreatorFunction("LagrangianPoints",
                                               FilterLagrangianPoints::create);

void StatLagrangianPoints::v_GetCoords(int id, Array<OneD, NekDouble> &gcoords)
{
    for (int d = 0; d < m_dim; ++d)
    {
        gcoords[d] = m_coords[0][d][id];
    }
}

void StatLagrangianPoints::v_SetCoords(int id,
                                       const Array<OneD, NekDouble> &gcoords)
{
    for (int d = 0; d < m_dim; ++d)
    {
        m_coords[0][d][id] = gcoords[d];
    }
}

void StatLagrangianPoints::v_GetPhysics(int id, Array<OneD, NekDouble> &data)
{
    for (int d = 0; d < m_dim; ++d)
    {
        data[d] = m_velocity[0][d][id];
    }
}

void StatLagrangianPoints::v_SetPhysics(int id,
                                        const Array<OneD, NekDouble> &data)
{
    for (int d = 0; d < m_dim; ++d)
    {
        m_velocity[0][d][id] = data[d];
    }
    for (size_t d = 0; d < m_extraPhysVars.size(); ++d)
    {
        m_extraPhysics[d][id] = data[d + m_dim];
    }
}

static void RollOver(Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &data)
{
    int n = data.size();
    if (n <= 1)
    {
        return;
    }
    Array<OneD, Array<OneD, NekDouble>> res = data[n - 1];
    for (int i = n - 1; i > 0; --i)
    {
        data[i] = data[i - 1];
    }
    data[0] = res;
}

void StatLagrangianPoints::v_TimeAdvance(int order)
{
    order = min(order, m_intOrder);
    if (order == 1)
    {
        for (int d = 0; d < m_dim; ++d)
        {
            Vmath::Svtvp(m_totPts, m_dt, m_velocity[0][d], 1, m_coords[0][d], 1,
                         m_coords[0][d], 1);
        }
        RollOver(m_velocity);
    }
    else if (order == 2)
    {
        for (int d = 0; d < m_dim; ++d)
        {
            Array<OneD, NekDouble> res = m_coords[1][d];
            Vmath::Svtvm(m_totPts, 3., m_velocity[0][d], 1, m_velocity[1][d], 1,
                         res, 1);
            Vmath::Svtvp(m_totPts, 0.5 * m_dt, res, 1, m_coords[0][d], 1, res,
                         1);
        }
        RollOver(m_coords);
        RollOver(m_velocity);
    }
    else
    {
        ASSERTL0(false, "Integration order not supported.");
    }
}

void StatLagrangianPoints::v_ReSize(int Np)
{
    if (Np > m_totPts)
    {
        for (int i = 0; i < m_coords.size(); ++i)
        {
            for (int d = 0; d < m_dim; ++d)
            {
                m_coords[i][d]   = Array<OneD, NekDouble>(Np);
                m_velocity[i][d] = Array<OneD, NekDouble>(Np);
            }
        }
    }
    m_totPts = Np;
    m_localIDToGlobal.clear();
    m_globalIDToLocal.clear();
}

static int BinaryWrite(std::ofstream &ofile, std::string str)
{
    int tmp = 0;
    for (size_t i = 0; i < str.size(); ++i)
    {
        tmp = str[i];
        ofile.write((char *)&tmp, 4);
    }
    tmp = 0;
    ofile.write((char *)&tmp, 4);
    return 4;
}

static int OutputTec360_binary(const std::string filename,
                               const std::vector<std::string> &variables,
                               const std::vector<int> &rawN,
                               std::vector<Array<OneD, NekDouble>> &data,
                               int isdouble)
{
    std::vector<int> N = rawN;
    for (int i = N.size(); i < 3; ++i)
    {
        N.push_back(1);
    }
    std::ofstream odata;
    odata.open(filename, std::ios::binary);
    if (!odata.is_open())
    {
        printf("error unable to open file %s\n", filename.c_str());
        return -1;
    }
    char tecplotversion[] = "#!TDV112";
    odata.write((char *)tecplotversion, 8);
    int value1 = 1;
    odata.write((char *)&value1, 4);
    int filetype = 0;
    odata.write((char *)&filetype, 4);
    // read file title and variable names
    std::string filetitle = "";
    BinaryWrite(odata, filetitle);
    int nvar = variables.size();
    odata.write((char *)&nvar, 4); // number of variables
    std::vector<std::string> vartitle;
    for (int i = 0; i < nvar; ++i)
    {
        BinaryWrite(odata, variables[i]);
    }
    float zonemarker299I = 299.0f;
    odata.write((char *)&zonemarker299I, 4);
    // zone title
    std::string zonetitle("ZONE 0");
    BinaryWrite(odata, zonetitle);
    int parentzone = -1;
    odata.write((char *)&parentzone, 4);
    int strandid = -1;
    odata.write((char *)&strandid, 4);
    double soltime = 0.0;
    odata.write((char *)&soltime, 8);
    int unused = -1;
    odata.write((char *)&unused, 4);
    int zonetype = 0;
    odata.write((char *)&zonetype, 4);
    int zero = 0;
    odata.write((char *)&zero, 4);
    odata.write((char *)&zero, 4);
    odata.write((char *)&zero, 4);
    for (int i = 0; i < 3; ++i)
    {
        int tmp = N[i];
        odata.write((char *)&tmp, 4);
    }

    odata.write((char *)&zero, 4);
    float eohmarker357 = 357.0f;
    odata.write((char *)&eohmarker357, 4);
    float zonemarker299II = 299.0f;
    odata.write((char *)&zonemarker299II, 4);
    std::vector<int> binarydatatype(nvar, 1 + (isdouble > 0));
    odata.write((char *)binarydatatype.data(), 4 * nvar);
    odata.write((char *)&zero, 4);
    odata.write((char *)&zero, 4);
    int minus1 = -1;
    odata.write((char *)&minus1, 4);

    int datanumber, datasize;
    datanumber = N[0] * N[1] * N[2];
    datasize   = N[0] * N[1] * N[2] * 4;
    for (int i = 0; i < nvar; ++i)
    {
        double minv = 0., maxv = 1.;
        for (int p = 0; p < datanumber; ++p)
        {
            if (maxv < data[i][p])
            {
                maxv = data[i][p];
            }
            if (minv > data[i][p])
            {
                minv = data[i][p];
            }
        }
        odata.write((char *)&minv, 8);
        odata.write((char *)&maxv, 8);
    }

    std::vector<float> vardata(datanumber);
    for (int i = 0; i < nvar; ++i)
    {
        if (isdouble)
        {
            odata.write((char *)data[i].data(), datasize * 2);
        }
        else
        {
            std::vector<float> fdata(datanumber);
            for (int j = 0; j < datanumber; ++j)
            {
                fdata[j] = data[i][j];
            }
            odata.write((char *)fdata.data(), datasize);
        }
    }
    odata.close();
    return 0;
}

void StatLagrangianPoints::v_OutputData(
    std::string filename, [[maybe_unused]] bool verbose,
    std::map<std::string, NekDouble> &params)
{
    std::vector<std::string> COORDS = {"x", "y", "z"};
    std::vector<std::string> VELOCI = {"u", "v", "w"};
    std::vector<std::string> variables;
    std::vector<Array<OneD, NekDouble>> data;
    for (int d = 0; d < m_dim; ++d)
    {
        variables.push_back(COORDS[d]);
    }
    for (int d = 0; d < m_dim; ++d)
    {
        variables.push_back(VELOCI[d]);
    }
    if (params.size() > 0)
    {
        for (int d = 0; d < m_dim; ++d)
        {
            Array<OneD, NekDouble> x(m_totPts, 0.);
            Vmath::Sadd(m_totPts, params[COORDS[d]], m_coords[0][d], 1, x, 1);
            data.push_back(x);
        }
        for (int d = 0; d < m_dim; ++d)
        {
            Array<OneD, NekDouble> x(m_totPts, 0.);
            Vmath::Sadd(m_totPts, params[VELOCI[d]], m_velocity[0][d], 1, x, 1);
            data.push_back(x);
        }
    }
    else
    {
        for (int d = 0; d < m_dim; ++d)
        {
            data.push_back(m_coords[0][d]);
        }
        for (int d = 0; d < m_dim; ++d)
        {
            data.push_back(m_velocity[0][d]);
        }
    }
    for (size_t d = 0; d < m_extraPhysics.size(); ++d)
    {
        data.push_back(m_extraPhysics[d]);
        variables.push_back(m_extraPhysVars[d]);
    }
    OutputTec360_binary(filename, variables, m_N, data, 1);
    if (verbose)
    {
        cout << "Write file " << filename << endl;
    }
}

void StatLagrangianPoints::v_OutputError()
{
    int Np                          = m_coords[0][0].size();
    std::vector<std::string> COORDS = {"x", "y", "z"};
    std::vector<std::string> VELOCI = {"u", "v", "w"};
    for (int d = 0; d < m_dim; ++d)
    {
        NekDouble value = 0.;
        for (int i = 0; i < Np; ++i)
        {
            value += m_coords[0][d][i] * m_coords[0][d][i];
        }
        value = sqrt(value / Np);
        cout << "L 2 error (variable L" << COORDS[d] << ") : " << value << endl;
    }
    for (int d = 0; d < m_dim; ++d)
    {
        NekDouble value = 0.;
        for (int i = 0; i < Np; ++i)
        {
            value += m_velocity[0][d][i] * m_velocity[0][d][i];
        }
        value = sqrt(value / Np);
        cout << "L 2 error (variable L" << VELOCI[d] << ") : " << value << endl;
    }
}

void StatLagrangianPoints::v_AssignPoint(int id, int pid,
                                         const Array<OneD, NekDouble> &gcoords)
{
    StationaryPoints::v_AssignPoint(id, pid, gcoords);
    for (int d = 0; d < m_dim; ++d)
    {
        m_coords[0][d][id] = gcoords[d];
    }
}

StatLagrangianPoints::StatLagrangianPoints([[maybe_unused]] int rank, int dim,
                                           int intOrder, int idOffset,
                                           NekDouble dt,
                                           const std::vector<int> &Np,
                                           const std::vector<NekDouble> &Box,
                                           std::vector<std::string> extraVars)
{
    m_idOffset = idOffset;
    m_dim      = dim;
    if (intOrder < 0)
    {
        m_intOrder = -intOrder;
        m_dt       = -dt;
    }
    else
    {
        m_intOrder = intOrder;
        m_dt       = dt;
    }
    if (Np.size() < m_dim || Box.size() < m_dim * 2)
    {
        throw ErrorUtil::NekError("Not enough input for StationaryPoints.");
    }
    Array<OneD, NekDouble> delta(3, 1.);
    m_N.resize(3, 1);
    m_totPts = 1;
    for (int d = 0; d < m_dim; ++d)
    {
        m_N[d] = Np[d];
        m_totPts *= Np[d];
        delta[d] = Np[d] < 2 ? 1. : (Box[d * 2 + 1] - Box[d * 2]) / (Np[d] - 1);
    }
    m_coords   = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_intOrder);
    m_velocity = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_intOrder);
    for (int i = 0; i < m_intOrder; ++i)
    {
        m_coords[i]   = Array<OneD, Array<OneD, NekDouble>>(m_dim);
        m_velocity[i] = Array<OneD, Array<OneD, NekDouble>>(m_dim);
        for (int d = 0; d < m_dim; ++d)
        {
            m_coords[i][d]   = Array<OneD, NekDouble>(m_totPts, 0.);
            m_velocity[i][d] = Array<OneD, NekDouble>(m_totPts, 0.);
        }
    }
    m_extraPhysics = Array<OneD, Array<OneD, NekDouble>>(extraVars.size());
    for (size_t i = 0; i < extraVars.size(); ++i)
    {
        m_extraPhysVars.push_back(extraVars[i]);
        m_extraPhysics[i] = Array<OneD, NekDouble>(m_totPts, 0.);
    }
    // initialise points
    int count = 0;
    for (int k = 0; k < m_N[2]; ++k)
    {
        for (int j = 0; j < m_N[1]; ++j)
        {
            for (int i = 0; i < m_N[0]; ++i)
            {
                if (m_dim > 0)
                {
                    m_coords[0][0][count] = Box[0] + delta[0] * i;
                }
                if (m_dim > 1)
                {
                    m_coords[0][1][count] = Box[2] + delta[1] * j;
                }
                if (m_dim > 2)
                {
                    m_coords[0][2][count] = Box[4] + delta[2] * k;
                }
                m_localIDToGlobal[count]              = m_idOffset + count;
                m_globalIDToLocal[m_idOffset + count] = count;
                ++count;
            }
        }
    }
}

/**
 *
 */
FilterLagrangianPoints::FilterLagrangianPoints(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::shared_ptr<EquationSystem> &pEquation,
    const std::map<std::string, std::string> &pParams)
    : Filter(pSession, pEquation), v_params(pParams)
{
}

/**
 *
 */
FilterLagrangianPoints::~FilterLagrangianPoints()
{
    m_ofstreamSamplePoints.close();
}

/**
 *
 */
void FilterLagrangianPoints::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    auto it = v_params.find("DefaultValues");
    if (it != v_params.end())
    {
        ParseUtils::GenerateVector(it->second, m_defaultValues);
    }
    EvaluatePoints::Initialise(pFields, time, m_defaultValues);
    NekDouble dt;
    m_session->LoadParameter("TimeStep", dt);

    it = v_params.find("FrameVelocity");
    if (it != v_params.end())
    {
        std::vector<std::string> strs;
        ParseUtils::GenerateVector(it->second, strs);
        for (int i = strs.size(); i < m_spacedim; ++i)
        {
            strs.push_back("0");
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_frame.m_frameVelFunction[i] =
                MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
                    pFields[0]->GetSession()->GetInterpreter(), strs[i]);
        }
    }

    it = v_params.find("FrameDisplacement");
    if (it != v_params.end())
    {
        std::vector<std::string> strs;
        ParseUtils::GenerateVector(it->second, strs);
        for (int i = strs.size(); i < m_spacedim; ++i)
        {
            strs.push_back("0");
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_frame.m_frameDispFunction[i] =
                MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
                    pFields[0]->GetSession()->GetInterpreter(), strs[i]);
        }
    }

    int intOrder = 1;
    it           = v_params.find("IntOrder");
    if (it != v_params.end())
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        intOrder = round(equ.Evaluate());
    }
    if (intOrder == 0)
    {
        intOrder          = 1;
        m_isMovablePoints = false;
    }
    else
    {
        m_isMovablePoints = true;
    }

    it = v_params.find("RootOutputL2Norm");
    if (it != v_params.end())
    {
        std::string value = it->second;
        m_outputL2Norm    = m_rank == 0 && (value == "1" || value == "Yes");
    }
    else
    {
        m_outputL2Norm = false;
    }

    it = v_params.find("OutputFile");
    std::string outputFilename;
    if (it == v_params.end())
    {
        outputFilename = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        outputFilename = it->second;
    }
    m_outputFile = outputFilename + "_R%05d_T%010.6lf.plt";

    // OutputFrequency
    it = v_params.find("OutputFrequency");
    if (it == v_params.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    it = v_params.find("OutputSampleFrequency");
    if (it == v_params.end())
    {
        m_outputSampleFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_outputSampleFrequency = round(equ.Evaluate());
    }

    it = v_params.find("Box_" + std::to_string(m_rank));
    if (it != v_params.end())
    {
        std::vector<NekDouble> values;
        ParseUtils::GenerateVector(it->second, values);
        if (values.size() < m_spacedim * 3 + 1)
        {
            throw ErrorUtil::NekError(
                "Box formate error "
                "(offset,Nx,Ny,Nx,xmin,xmax,ymin,ymax,zmin,zmax): " +
                it->second);
        }
        m_frame.Update(time);
        std::vector<int> Np(m_spacedim, 1);
        m_box.resize(2 * m_spacedim, 0.);
        int idOffset = std::round(values[0]);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Np[d] = std::round(values[d + 1]);
        }
        for (int d = 0; d < 2 * m_spacedim; ++d)
        {
            m_box[d] = values[1 + d + m_spacedim] - m_frame.m_frameDisp[d / 2];
        }
        std::vector<std::string> extraVars;
        ExtraPhysicsVars(extraVars);
        m_staticPts = MemoryManager<StatLagrangianPoints>::AllocateSharedPtr(
            m_rank, m_spacedim, intOrder, idOffset, dt, Np, m_box, extraVars);
        CopyStaticPtsToMobile();

        it = v_params.find("OutputSamplePoints");
        std::vector<unsigned int> samplepts;
        if (it != v_params.end())
        {
            ParseUtils::GenerateSeqVector(it->second, samplepts);
        }
        std::set<int> sampleIds(samplepts.begin(), samplepts.end());
        it = v_params.find("OutputSamplePointsCondition");
        LibUtilities::EquationSharedPtr sampleCondition;
        if (it != v_params.end())
        {
            sampleCondition =
                MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
                    pFields[0]->GetSession()->GetInterpreter(), it->second);
        }

        for (int p = 0; p < m_staticPts->GetTotPoints(); ++p)
        {
            Array<OneD, NekDouble> data(m_spacedim, 0.);
            m_staticPts->GetCoords(p, data);
            int globalId = m_staticPts->LocalToGlobal(p);
            NekDouble x  = data[0] + m_frame.m_frameDisp[0],
                      y  = data[1] + m_frame.m_frameDisp[1];
            NekDouble z =
                m_spacedim > 2 ? data[2] : 0. + m_frame.m_frameDisp[2];
            bool issample =
                (sampleCondition == nullptr)
                    ? false
                    : (sampleCondition->Evaluate(x, y, z, time) > 0);
            if (sampleIds.find(globalId) != sampleIds.end() || issample)
            {
                m_samplePointIDs.insert(globalId);
            }
        }
    }
    if (!m_samplePointIDs.empty())
    {
        std::string sampleFilename = outputFilename + "_Sample_R%05d.dat";
        boost::format filename(sampleFilename);
        filename % m_rank;
        m_ofstreamSamplePoints.open(filename.str(), std::ofstream::out);
        std::vector<std::string> COORDS = {"x", "y", "z"};
        std::vector<std::string> VELOCI = {"u", "v", "w"};
        m_ofstreamSamplePoints << "variables = pid, time ";
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_ofstreamSamplePoints << COORDS[d] << " ";
        }
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_ofstreamSamplePoints << VELOCI[d] << " ";
        }
        m_ofstreamSamplePoints << endl;
    }
    SetUpCommInfo();
    if (m_updateOnInitialise)
    {
        v_Update(pFields, time);
    }
}

void FilterLagrangianPoints::v_ModifyVelocity(
    [[maybe_unused]] Array<OneD, NekDouble> gcoords,
    [[maybe_unused]] NekDouble time, Array<OneD, NekDouble> vel)
{
    if (m_frame.m_frameVelFunction.size() > 0)
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] -= m_frame.m_frameVel[i];
        }
    }
}

void FilterLagrangianPoints::ExtraPhysicsVars(
    std::vector<std::string> &extraVars)
{
    if (m_spacedim == 2)
    {
        extraVars.push_back("u_x");
        extraVars.push_back("u_y");
        extraVars.push_back("v_x");
        extraVars.push_back("v_y");
        extraVars.push_back("LW_z");
    }
    else if (m_spacedim == 3)
    {
        extraVars.push_back("u_x");
        extraVars.push_back("u_y");
        extraVars.push_back("u_z");
        extraVars.push_back("v_x");
        extraVars.push_back("v_y");
        extraVars.push_back("v_z");
        extraVars.push_back("w_x");
        extraVars.push_back("w_y");
        extraVars.push_back("w_z");
        extraVars.push_back("LW_x");
        extraVars.push_back("LW_y");
        extraVars.push_back("LW_z");
    }
}

void FilterLagrangianPoints::GetPhysicsData(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    std::vector<Array<OneD, NekDouble>> &PhysicsData)
{
    int npts = pFields[0]->GetTotPoints();
    for (int i = 0; i < m_spacedim; ++i)
    {
        PhysicsData.push_back(pFields[i]->UpdatePhys());
    }
    int offset = PhysicsData.size();
    for (int d = 0; d < m_spacedim; ++d)
    {
        Array<OneD, Array<OneD, NekDouble>> data(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            data[i] = Array<OneD, NekDouble>(npts, 0.);
        }
        if (m_spacedim == 2)
        {
            pFields[d]->PhysDeriv(pFields[d]->UpdatePhys(), data[0], data[1]);
        }
        else if (m_spacedim == 3)
        {
            pFields[d]->PhysDeriv(pFields[d]->UpdatePhys(), data[0], data[1],
                                  data[2]);
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            PhysicsData.push_back(data[i]);
        }
    }
    int vortexdim = m_spacedim == 3 ? 3 : 1;
    Array<OneD, Array<OneD, NekDouble>> vorticity(vortexdim);
    Array<OneD, Array<OneD, NekDouble>> Ldata(vortexdim);
    Array<OneD, Array<OneD, NekDouble>> temp(m_spacedim);
    Array<OneD, NekDouble> temps(npts, 0.);
    for (int d = 0; d < vortexdim; ++d)
    {
        Ldata[d]     = Array<OneD, NekDouble>(npts, 0.);
        vorticity[d] = Array<OneD, NekDouble>(npts, 0.);
    }
    for (int d = 0; d < m_spacedim; ++d)
    {
        temp[d] = Array<OneD, NekDouble>(npts, 0.);
    }
    if (m_spacedim == 2)
    {
        // 0 ux, 1 uy, 2 vx, 3 vy
        Vmath::Vsub(npts, PhysicsData[offset + 2], 1, PhysicsData[offset + 1],
                    1, vorticity[0], 1); // vx-uy
    }
    else if (m_spacedim == 3)
    {
        // 0 ux, 1 uy, 2 uz, 3 vx, 4 vy, 5 vz, 6 wx, 7 wy, 8 wz
        Vmath::Vsub(npts, PhysicsData[offset + 7], 1, PhysicsData[offset + 5],
                    1, vorticity[0], 1); // wy-vz
        Vmath::Vsub(npts, PhysicsData[offset + 2], 1, PhysicsData[offset + 6],
                    1, vorticity[1], 1); // uz-wx
        Vmath::Vsub(npts, PhysicsData[offset + 3], 1, PhysicsData[offset + 1],
                    1, vorticity[2], 1); // vx-uy
    }
    for (int d = 0; d < vortexdim; ++d)
    {
        if (m_spacedim == 2)
        {
            pFields[0]->PhysDeriv(vorticity[d], temp[0], temp[1]);
        }
        else if (m_spacedim == 3)
        {
            pFields[0]->PhysDeriv(vorticity[d], temp[0], temp[1], temp[2]);
        }
        for (int j = 0; j < m_spacedim; ++j)
        {
            pFields[0]->PhysDeriv(j, temp[j], temps);
            Vmath::Vadd(npts, temps, 1, Ldata[d], 1, Ldata[d], 1);
        }
        PhysicsData.push_back(Ldata[d]);
    }
}

void FilterLagrangianPoints::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    if (m_isMovablePoints || m_index == 0)
    {
        PartitionMobilePts(pFields);
    }
    std::vector<Array<OneD, NekDouble>> PhysicsData;
    GetPhysicsData(pFields, PhysicsData);
    m_frame.Update(time);
    EvaluateMobilePhys(pFields[0], PhysicsData, time);
    std::map<int, std::set<int>> callbackUpdateMobCoords;
    PassMobilePhysToStatic(callbackUpdateMobCoords);
    // output Lagrangian coordinates

    if (m_index % m_outputFrequency == 0 && m_staticPts != nullptr)
    {
        OutputStatPoints(time);
        if (m_outputL2Norm)
        {
            m_staticPts->OutputError();
        }
    }

    // output sample points
    if (m_index % m_outputSampleFrequency == 0 && !m_samplePointIDs.empty())
    {
        OutputSamplePoints(time);
    }
    // update Lagrangian coordinates
    if (m_isMovablePoints)
    {
        if (m_staticPts != nullptr)
        {
            m_staticPts->TimeAdvance(m_index + 1);
        }
        // return the updated global coordinates to mobile points
        PassStaticCoordsToMobile(callbackUpdateMobCoords);
    }
    ++m_index;
}

void FilterLagrangianPoints::OutputSamplePoints(NekDouble time)
{
    m_ofstreamSamplePoints << "ZONE T=\"" << time
                           << "\", I=" << m_samplePointIDs.size() << endl;
    for (auto &pid : m_samplePointIDs)
    {
        m_ofstreamSamplePoints << pid << " " << time << " ";
        Array<OneD, NekDouble> data(m_spacedim, 0.);
        m_staticPts->GetCoordsByPID(pid, data);
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_ofstreamSamplePoints << data[d] + m_frame.m_frameDisp[d] << " ";
        }
        m_staticPts->GetPhysicsByPID(pid, data);
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_ofstreamSamplePoints << data[d] + m_frame.m_frameVel[d] << " ";
        }
        m_ofstreamSamplePoints << endl;
    }
}

void FilterLagrangianPoints::OutputStatPoints(NekDouble time)
{
    boost::format filename(m_outputFile);
    filename % m_rank % time;
    std::map<std::string, NekDouble> params;
    if (m_frame.m_frameVelFunction.size())
    {
        std::vector<std::string> COORDS = {"x", "y", "z"};
        std::vector<std::string> VELOCI = {"u", "v", "w"};
        for (auto it : m_frame.m_frameVelFunction)
        {
            params[VELOCI[it.first]] = m_frame.m_frameVel[it.first];
        }
        for (auto it : m_frame.m_frameDispFunction)
        {
            params[COORDS[it.first]] = m_frame.m_frameDisp[it.first];
        }
    }
    m_staticPts->OutputData(filename.str(), m_outputL2Norm, params);
}

void FilterLagrangianPoints::v_Finalise(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    [[maybe_unused]] const NekDouble &time)
{
}

bool FilterLagrangianPoints::v_IsTimeDependent()
{
    return true;
}

} // namespace Nektar::SolverUtils
