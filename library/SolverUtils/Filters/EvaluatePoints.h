///////////////////////////////////////////////////////////////////////////////
//
// File: EvaluatePoints.h
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

#ifndef NEKTAR_SOLVERUTILS_FILTERS_EVALUATEPOINTS_H
#define NEKTAR_SOLVERUTILS_FILTERS_EVALUATEPOINTS_H

#include "LibUtilities/BasicConst/NektarUnivTypeDefs.hpp"
#include "LibUtilities/BasicUtils/SharedArray.hpp"
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/EquationSystem.h>

#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar::SolverUtils
{

class MobilePoint;
class EvaluatePoints;
class StationaryPoints;

typedef std::shared_ptr<MobilePoint> MobilePointSharedPtr;
typedef std::shared_ptr<StationaryPoints> StationaryPointsSharedPtr;

class StationaryPoints
{
public:
    virtual ~StationaryPoints() = default;
    SOLVER_UTILS_EXPORT void OutputData(
        std::string filename, bool verbose,
        std::map<std::string, NekDouble> &params)
    {
        v_OutputData(filename, verbose, params);
    }
    SOLVER_UTILS_EXPORT void OutputError()
    {
        v_OutputError();
    }
    inline void ReSize(int Np)
    {
        v_ReSize(Np);
    }
    inline void AssignPoint(int id, int pid,
                            const Array<OneD, NekDouble> &gcoords)
    {
        v_AssignPoint(id, pid, gcoords);
    }
    inline void GetCoordsByPID(int pid, Array<OneD, NekDouble> &gcoords)
    {
        v_GetCoords(GlobalToLocal(pid), gcoords);
    }
    inline void SetCoordsByPID(int pid, const Array<OneD, NekDouble> &gcoords)
    {
        v_SetCoords(GlobalToLocal(pid), gcoords);
    }
    inline void GetPhysicsByPID(int pid, Array<OneD, NekDouble> &data)
    {
        v_GetPhysics(GlobalToLocal(pid), data);
    }
    inline void SetPhysicsByPID(int pid, const Array<OneD, NekDouble> &data)
    {
        v_SetPhysics(GlobalToLocal(pid), data);
    }
    inline void GetCoords(int pid, Array<OneD, NekDouble> &gcoords)
    {
        v_GetCoords(pid, gcoords);
    }
    inline void SetCoords(int pid, const Array<OneD, NekDouble> &gcoords)
    {
        v_SetCoords(pid, gcoords);
    }
    inline void GetPhysics(int pid, Array<OneD, NekDouble> &data)
    {
        v_GetPhysics(pid, data);
    }
    inline void SetPhysics(int pid, const Array<OneD, NekDouble> &data)
    {
        v_SetPhysics(pid, data);
    }
    inline void TimeAdvance(int order)
    {
        v_TimeAdvance(order);
    }
    inline int GetTotPoints()
    {
        return m_totPts;
    }
    inline int GetDim()
    {
        return m_dim;
    }
    inline int LocalToGlobal(int id)
    {
        return m_localIDToGlobal[id];
    }
    inline int GlobalToLocal(int id)
    {
        return m_globalIDToLocal[id];
    }

protected:
    virtual void v_OutputData(std::string filename, bool verbose,
                              std::map<std::string, NekDouble> &params)    = 0;
    virtual void v_OutputError()                                           = 0;
    virtual void v_TimeAdvance(int order)                                  = 0;
    virtual void v_GetCoords(int pid, Array<OneD, NekDouble> &gcoords)     = 0;
    virtual void v_SetCoords(int pid,
                             const Array<OneD, NekDouble> &gcoords)        = 0;
    virtual void v_GetPhysics(int pid, Array<OneD, NekDouble> &data)       = 0;
    virtual void v_SetPhysics(int pid, const Array<OneD, NekDouble> &data) = 0;
    virtual void v_ReSize(int Np)                                          = 0;
    virtual void v_AssignPoint(int id, int pid,
                               const Array<OneD, NekDouble> &gcoords);
    int m_dim;
    int m_totPts;
    std::map<int, int> m_localIDToGlobal; // size should be m_totPts
    std::map<int, int> m_globalIDToLocal; // size should be m_totPts
};

class MobilePoint
{
public:
    friend class MemoryManager<MobilePoint>;
    int m_sRank;
    Array<OneD, NekDouble> m_global;
    int m_eId;
    Array<OneD, NekDouble> m_local;
    Array<OneD, NekDouble> m_data;

    /// Creates an instance of this class
    static MobilePointSharedPtr create(int dim, int sRank,
                                       const Array<OneD, NekDouble> global)
    {
        MobilePointSharedPtr p =
            MemoryManager<MobilePoint>::AllocateSharedPtr(dim, sRank, global);
        return p;
    }

    MobilePoint(int dim, int sRank, const Array<OneD, NekDouble> global)
        : m_sRank(sRank)
    {
        m_global = Array<OneD, NekDouble>(dim);
        Vmath::Vcopy(dim, global, 1, m_global, 1);
        m_local = Array<OneD, NekDouble>(dim);
    }

    void SetLocalCoords(int eId, const Array<OneD, NekDouble> &local)
    {
        m_eId = eId;
        Vmath::Vcopy(m_local.size(), local, 1, m_local, 1);
    }

    void SetGlobalCoords(const Array<OneD, NekDouble> &global)
    {
        Vmath::Vcopy(m_global.size(), global, 1, m_global, 1);
    }

    void SetData(int nPhys, const Array<OneD, NekDouble> data)
    {
        if (m_data.size() < nPhys)
        {
            m_data = Array<OneD, NekDouble>(nPhys);
        }
        Vmath::Vcopy(nPhys, data, 1, m_data, 1);
    }

    inline Array<OneD, NekDouble> GetData()
    {
        return m_data;
    }

    inline Array<OneD, NekDouble> GetGlobalCoords()
    {
        return m_global;
    }

    inline Array<OneD, NekDouble> GetLocalCoords()
    {
        return m_local;
    }
};

class EvaluatePoints
{
public:
    SOLVER_UTILS_EXPORT EvaluatePoints();
    virtual ~EvaluatePoints();
    SOLVER_UTILS_EXPORT void PartitionMobilePts(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
    SOLVER_UTILS_EXPORT void EvaluateMobilePhys(
        const MultiRegions::ExpListSharedPtr &pField,
        std::vector<Array<OneD, NekDouble>> &PhysicsData, NekDouble time);
    SOLVER_UTILS_EXPORT void PassMobilePhysToStatic(
        std::map<int, std::set<int>> &callbackUpdateMobCoords);
    SOLVER_UTILS_EXPORT void PassStaticCoordsToMobile(
        std::map<int, std::set<int>> &callbackUpdateMobCoords);
    SOLVER_UTILS_EXPORT void CopyStaticPtsToMobile();
    SOLVER_UTILS_EXPORT void CopyMobilePtsToStatic();
    SOLVER_UTILS_EXPORT void Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time, std::vector<std::string> &defaultValues);
    SOLVER_UTILS_EXPORT void SetUpCommInfo();

protected:
    void UpdateMobileCoords(
        std::map<int, Array<OneD, NekDouble>> &globalCoords);
    void GatherMobilePhysics(std::map<int, Array<OneD, NekDouble>> &revData);
    void PartitionLocalPoints(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        std::set<int> &notfound);
    void PartitionExchangeNonlocalPoints(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        std::set<int> &notfound);
    void SyncColumnComm();
    NekDouble GetDefaultValue(const int i, const Array<OneD, const NekDouble> x,
                              const NekDouble time);
    void Pack2Int(const int &a, const int &b, double &d);
    void unPack2Int(int &a, int &b, const double &d);
    void Pack3Int(const int &a, const int &b, const int &c, double &d);
    void unPack3Int(int &a, int &b, int &c, const double &d);
    void Pack2Short(const int &a, const int &b, int &c);
    void unPack2Short(int &a, int &b, const int &c);
    virtual void v_ModifyVelocity(Array<OneD, NekDouble> gcoords,
                                  NekDouble time, Array<OneD, NekDouble> vel);
    int m_spacedim;
    int m_rank;
    int m_rowRank;
    int m_columnRank;
    bool m_isHomogeneous1D;
    int m_nPhysics;
    LibUtilities::CommSharedPtr m_comm;
    StationaryPointsSharedPtr m_staticPts;
    std::map<int, MobilePointSharedPtr> m_mobilePts;
    std::map<int, LibUtilities::EquationSharedPtr> m_defaultValueFunction;
    // m_recvMobInfo[t] = n, in the current stationary pts, there are n mobile
    // points from thread t; all threads; global
    // m_recvStatInfo[t] = n, in the
    // current mobile pts, there are n stationary points from thread t;
    // m_columnRank == 0 only; local
    std::map<int, int> m_recvMobInfo;
    std::map<int, int> m_recvStatInfo;
};

} // namespace Nektar::SolverUtils

#endif /* NEKTAR_SOLVERUTILS_FILTERS_EVALUATEPOINTS_H */
