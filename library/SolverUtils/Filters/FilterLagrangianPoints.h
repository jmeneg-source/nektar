///////////////////////////////////////////////////////////////////////////////
//
// File: FilterLagrangianPoints.h
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
// Description: base class for physics evaluation.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERLAGRANGIANPOINTS_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERLAGRANGIANPOINTS_H

#include "LibUtilities/BasicConst/NektarUnivTypeDefs.hpp"
#include "LibUtilities/BasicUtils/SharedArray.hpp"
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <SolverUtils/Filters/EvaluatePoints.h>
#include <SolverUtils/Filters/Filter.h>

namespace Nektar::SolverUtils
{

class FilterLagrangianPoints;

class StatLagrangianPoints : public StationaryPoints
{
public:
    friend class MemoryManager<StatLagrangianPoints>;

    /// Creates an instance of this class
    static StationaryPointsSharedPtr create(int rank, int dim, int intOrder,
                                            int idOffset, NekDouble dt,
                                            const std::vector<int> &Np,
                                            const std::vector<NekDouble> &Box,
                                            std::vector<std::string> extraVars)
    {
        StationaryPointsSharedPtr p =
            MemoryManager<StatLagrangianPoints>::AllocateSharedPtr(
                rank, dim, intOrder, idOffset, dt, Np, Box, extraVars);
        return p;
    }

    StatLagrangianPoints(int rank, int dim, int intOrder, int idOffset,
                         NekDouble dt, const std::vector<int> &Np,
                         const std::vector<NekDouble> &Box,
                         std::vector<std::string> extraVars);

    ~StatLagrangianPoints() override
    {
    }

protected:
    void v_OutputData(std::string filename, bool verbose,
                      std::map<std::string, NekDouble> &params) override;
    void v_OutputError() override;
    void v_TimeAdvance(int order) override;
    void v_GetCoords(int pid, Array<OneD, NekDouble> &gcoords) override;
    void v_SetCoords(int pid, const Array<OneD, NekDouble> &gcoords) override;
    void v_GetPhysics(int pid, Array<OneD, NekDouble> &data) override;
    void v_SetPhysics(int pid, const Array<OneD, NekDouble> &data) override;
    void v_ReSize(int Np) override;
    void v_AssignPoint(int id, int pid,
                       const Array<OneD, NekDouble> &gcoords) override;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_coords;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_velocity;
    Array<OneD, Array<OneD, NekDouble>> m_extraPhysics;
    std::vector<std::string> m_extraPhysVars;
    NekDouble m_dt;
    int m_idOffset;
    std::vector<int> m_N;
    int m_intOrder;
};

class FilterLagrangianPoints : public EvaluatePoints, public Filter
{
public:
    friend class MemoryManager<FilterLagrangianPoints>;

    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<EquationSystem> &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p =
            MemoryManager<FilterLagrangianPoints>::AllocateSharedPtr(
                pSession, pEquation, pParams);
        return p;
    }

    static std::string className;

    SOLVER_UTILS_EXPORT FilterLagrangianPoints(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<EquationSystem> &pEquation,
        const std::map<std::string, std::string> &pParams);
    SOLVER_UTILS_EXPORT ~FilterLagrangianPoints() override;

protected:
    void ExtraPhysicsVars(std::vector<std::string> &extraVars);
    SOLVER_UTILS_EXPORT void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;

    void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;

    void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    bool v_IsTimeDependent() override;
    void v_ModifyVelocity(Array<OneD, NekDouble> gcoords, NekDouble time,
                          Array<OneD, NekDouble> vel) override;
    void GetPhysicsData(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        std::vector<Array<OneD, NekDouble>> &PhysicsData);
    void OutputSamplePoints(NekDouble time);
    void OutputStatPoints(NekDouble time);
    std::vector<std::string> m_defaultValues;
    struct MovingFrame
    {
        std::map<int, LibUtilities::EquationSharedPtr> m_frameVelFunction;
        std::map<int, LibUtilities::EquationSharedPtr> m_frameDispFunction;
        std::map<int, NekDouble> m_frameVel;
        std::map<int, NekDouble> m_frameDisp;
        void Update(NekDouble time)
        {
            for (auto it : m_frameVelFunction)
            {
                m_frameVel[it.first] = it.second->Evaluate(0, 0, 0, time);
            }
            for (auto it : m_frameDispFunction)
            {
                m_frameDisp[it.first] = it.second->Evaluate(0, 0, 0, time);
            }
        }
    } m_frame;
    std::vector<NekDouble> m_box;
    std::map<std::string, std::string> v_params;
    unsigned int m_index = 0;
    unsigned int m_outputFrequency;
    std::string m_outputFile;
    bool m_outputL2Norm;
    bool m_isMovablePoints;
    unsigned int m_outputSampleFrequency;
    std::set<int> m_samplePointIDs;
    std::ofstream m_ofstreamSamplePoints;
};

} // namespace Nektar::SolverUtils

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERLAGRANGIANPOINTS_H */
