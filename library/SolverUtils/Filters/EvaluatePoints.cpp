///////////////////////////////////////////////////////////////////////////////
//
// File: EvaluatePoints.cpp
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

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <SolverUtils/Filters/EvaluatePoints.h>
#include <boost/format.hpp>

namespace Nektar::SolverUtils
{

void StationaryPoints::v_AssignPoint(
    int id, int pid, [[maybe_unused]] const Array<OneD, NekDouble> &gcoords)
{
    m_localIDToGlobal[id]  = pid;
    m_globalIDToLocal[pid] = id;
}

void EvaluatePoints::Pack2Int(const int &a, const int &b, double &d)
{
    int *p = (int *)&d;
    p[0]   = a;
    p[1]   = b;
}

void EvaluatePoints::unPack2Int(int &a, int &b, const double &d)
{
    int *p = (int *)&d;
    a      = p[0];
    b      = p[1];
}

void EvaluatePoints::Pack3Int(const int &a, const int &b, const int &c,
                              double &d)
{
    int *p                = (int *)&d;
    unsigned short int *q = (unsigned short int *)&d;
    p[0]                  = a;
    q[2]                  = (unsigned short int)b;
    q[3]                  = (unsigned short int)c;
}

void EvaluatePoints::unPack3Int(int &a, int &b, int &c, const double &d)
{
    int *p                = (int *)&d;
    unsigned short int *q = (unsigned short int *)&d;
    a                     = p[0];
    b                     = q[2];
    c                     = q[3];
}

void EvaluatePoints::Pack2Short(const int &a, const int &b, int &c)
{
    short int *p = (short int *)&c;
    p[0]         = (short int)a;
    p[1]         = (short int)b;
}

void EvaluatePoints::unPack2Short(int &a, int &b, const int &c)
{
    short int *p = (short int *)&c;
    a            = (int)p[0];
    b            = (int)p[1];
}
/**
 *
 */
EvaluatePoints::EvaluatePoints()
{
}

/**
 *
 */
EvaluatePoints::~EvaluatePoints()
{
}

/**
 *
 */
void EvaluatePoints::Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const NekDouble &time,
    std::vector<std::string> &defaultValues)
{
    m_comm       = pFields[0]->GetComm();
    m_rank       = m_comm->GetRank();
    m_columnRank = m_comm->GetColumnComm()->GetRank();
    m_rowRank    = m_comm->GetRowComm()->GetRank();
    LibUtilities::SessionReaderSharedPtr session = pFields[0]->GetSession();
    session->MatchSolverInfo("Homogeneous", "1D", m_isHomogeneous1D, false);
    int coordim    = pFields[0]->GetGraph()->GetSpaceDimension();
    int numHomoDir = m_isHomogeneous1D ? 1 : 0;
    m_spacedim     = coordim + numHomoDir;
    ASSERTL0(m_rank <= 65535,
             "EvaluatePoints does not support more than 65536 threads.")

    for (size_t i = 0; i < defaultValues.size(); ++i)
    {
        try
        {
            m_defaultValueFunction[i] =
                MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
                    session->GetInterpreter(), defaultValues[i]);
        }
        catch (const std::runtime_error &)
        {
            NEKERROR(ErrorUtil::efatal,
                     "Error parsering expression" + defaultValues[i]);
        }
    }
}

/**
 *
 */
void EvaluatePoints::EvaluateMobilePhys(
    const MultiRegions::ExpListSharedPtr &pField,
    std::vector<Array<OneD, NekDouble>> &PhysicsData, NekDouble time)
{
    m_nPhysics  = PhysicsData.size();
    int coordim = pField->GetGraph()->GetSpaceDimension();
    Array<OneD, NekDouble> data(m_mobilePts.size() * m_nPhysics, 0.0);
    Array<OneD, NekDouble> physvals;
    Array<OneD, NekDouble> locCoord;
    if (m_isHomogeneous1D)
    {
        Array<OneD, const unsigned int> planes = pField->GetZIDs();
        int nPlanes    = pField->GetHomogeneousBasis()->GetZ().size();
        NekDouble lHom = pField->GetHomoLen();
        int nelmts     = pField->GetPlane(0)->GetExpSize();
        int indx       = -1;
        for (auto &x : m_mobilePts)
        {
            ++indx;
            locCoord  = x.second->GetLocalCoords();
            int expId = x.second->m_eId;
            if (expId < 0)
            {
                continue;
            }
            NekDouble BetaT = 2. * M_PI * fmod(locCoord[coordim], lHom) / lHom;
            std::vector<NekDouble> basis(planes.size(), 0.);
            for (size_t n = 0; n < planes.size(); ++n)
            {
                if (planes[n] == 0)
                {
                    basis[n] = 1.;
                }
                else if (planes[n] == 1)
                {
                    basis[n] = cos(0.5 * nPlanes * BetaT);
                }
                else if (planes[n] % 2 == 0)
                {
                    NekDouble phase = (planes[n] >> 1) * BetaT;
                    basis[n]        = cos(phase);
                }
                else
                {
                    NekDouble phase = (planes[n] >> 1) * BetaT;
                    basis[n]        = -sin(phase);
                }
            }
            LocalRegions::ExpansionSharedPtr elmt =
                pField->GetPlane(0)->GetExp(expId);
            for (int j = 0; j < m_nPhysics; ++j)
            {
                NekDouble value = 0.0;
                for (size_t n = 0; n < planes.size(); ++n)
                {
                    physvals = PhysicsData[j] +
                               pField->GetPhys_Offset(expId + n * nelmts);
                    NekDouble coeff = elmt->StdPhysEvaluate(locCoord, physvals);
                    value += basis[n] * coeff;
                }
                data[indx * m_nPhysics + j] = value;
            }
        }
        m_comm->GetColumnComm()->AllReduce(data, LibUtilities::ReduceSum);
    }
    else
    {
        for (int j = 0; j < m_nPhysics; ++j)
        {
            int indx = -1;
            for (auto &x : m_mobilePts)
            {
                ++indx;
                locCoord  = x.second->GetLocalCoords();
                int expId = x.second->m_eId;
                if (expId < 0)
                {
                    continue;
                }
                physvals = PhysicsData[j] + pField->GetPhys_Offset(expId);
                data[indx * m_nPhysics + j] =
                    pField->GetExp(expId)->StdPhysEvaluate(locCoord, physvals);
            }
        }
    }
    // copy data
    int count = 0;
    for (auto &x : m_mobilePts)
    {
        Array<OneD, NekDouble> datum;
        if (x.second->m_eId >= 0)
        {
            datum = data + count * m_nPhysics;
        }
        else
        {
            datum = Array<OneD, NekDouble>(m_nPhysics, 0.);
            if (m_columnRank == 0)
            {
                for (int j = 0; j < m_nPhysics; ++j)
                {
                    datum[j] =
                        GetDefaultValue(j, x.second->GetGlobalCoords(), time);
                }
            }
        }
        v_ModifyVelocity(x.second->GetGlobalCoords(), time, datum);
        x.second->SetData(m_nPhysics, datum);
        ++count;
    }
}

void EvaluatePoints::v_ModifyVelocity(
    [[maybe_unused]] Array<OneD, NekDouble> gcoords,
    [[maybe_unused]] NekDouble time,
    [[maybe_unused]] Array<OneD, NekDouble> vel)
{
}

// (pid, mobile thread, stationary thread, data)
void EvaluatePoints::GatherMobilePhysics(
    std::map<int, Array<OneD, NekDouble>> &revData)
{
    revData.clear();
    // send/receive data
    int nPackageSize = m_nPhysics + 1;
    if (m_columnRank == 0)
    {
        std::map<int, std::set<int>> threadToSend;
        for (const auto &p : m_mobilePts)
        {
            threadToSend[p.second->m_sRank].insert(p.first);
        }
        for (const auto &thread : threadToSend)
        {
            Array<OneD, NekDouble> dataToSend(
                thread.second.size() * nPackageSize, 0.);
            Array<OneD, NekDouble> tmp;
            int offset = 0;
            for (const int p : thread.second)
            {
                Pack3Int(p, m_rank, thread.first, dataToSend[offset]);
                Vmath::Vcopy(m_nPhysics, m_mobilePts[p]->GetData(), 1,
                             tmp = dataToSend + offset + 1, 1);
                offset += nPackageSize;
            }
            if (thread.first != m_rank)
            {
                m_comm->Send(thread.first, dataToSend);
            }
            else
            {
                revData[thread.first] = dataToSend;
            }
        }
    }
    for (const auto &p : m_recvMobInfo)
    {
        if (p.first != m_rank)
        {
            Array<OneD, NekDouble> dataToRecv(p.second * nPackageSize, 0.);
            m_comm->Recv(p.first, dataToRecv);
            revData[p.first] = dataToRecv;
        }
    }
}

void EvaluatePoints::SetUpCommInfo()
{
    // setup communication pair;
    int nRows   = m_comm->GetRowComm()->GetSize();
    int nColumn = m_comm->GetColumnComm()->GetSize();
    //(from short/int, to short/int, size/int)
    m_recvStatInfo.clear();
    int totPairs = 0;
    Array<OneD, int> data;
    if (m_columnRank == 0)
    {
        for (const auto &p : m_mobilePts)
        {
            m_recvStatInfo[p.second->m_sRank] += 1;
        }
        Array<OneD, int> npairs(nRows, 0);
        npairs[m_rowRank] = m_recvStatInfo.size();
        m_comm->GetRowComm()->AllReduce(npairs, LibUtilities::ReduceMax);
        int offset = 0;
        for (int t = 0; t < m_rowRank; ++t)
        {
            offset += npairs[t];
        }
        totPairs = offset;
        for (int t = m_rowRank; t < nRows; ++t)
        {
            totPairs += npairs[t];
        }
        data = Array<OneD, int>(totPairs * 2, 0);
        offset *= 2;
        for (const auto &p : m_recvStatInfo)
        {
            Pack2Short((short int)m_rowRank, (short int)p.first,
                       data[offset++]);
            data[offset++] = p.second;
        }
        m_comm->GetRowComm()->AllReduce(data, LibUtilities::ReduceSum);
    }
    // set up to-receive mobile points
    if (nColumn > 1)
    {
        m_comm->GetColumnComm()->AllReduce(totPairs, LibUtilities::ReduceMax);
        if (m_columnRank > 0)
        {
            data = Array<OneD, int>(totPairs * 2, 0);
        }
        m_comm->GetColumnComm()->AllReduce(data, LibUtilities::ReduceSum);
    }
    m_recvMobInfo.clear();
    for (int t = 0; t < totPairs; ++t)
    {
        int from, to;
        unPack2Short(from, to, data[2 * t]);
        if (to == m_rank)
        {
            m_recvMobInfo[(int)from] = data[2 * t + 1];
        }
    }
}

void EvaluatePoints::PartitionMobilePts(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    std::set<int> notfound;
    PartitionLocalPoints(pFields, notfound);
    PartitionExchangeNonlocalPoints(pFields, notfound);
    SyncColumnComm();
}

void EvaluatePoints::PassMobilePhysToStatic(
    std::map<int, std::set<int>> &callbackUpdateMobCoords)
{
    std::map<int, Array<OneD, NekDouble>> recvData;
    GatherMobilePhysics(recvData);
    //   unpack data
    for (auto const &data : recvData)
    {
        int offset = 0;
        while (offset < data.second.size())
        {
            int pid, mobt, statt;
            unPack3Int(pid, mobt, statt, data.second[offset++]);
            m_staticPts->SetPhysicsByPID(pid, data.second + offset);
            callbackUpdateMobCoords[data.first].insert(pid);
            offset += m_nPhysics;
        }
    }
}

void EvaluatePoints::PassStaticCoordsToMobile(
    std::map<int, std::set<int>> &callbackUpdateMobCoords)
{
    std::map<int, Array<OneD, NekDouble>> globalCoords;
    for (auto const &thread : callbackUpdateMobCoords)
    {
        int nDataSize = m_spacedim + 1, offset = 0;
        Array<OneD, NekDouble> data(thread.second.size() * nDataSize), tmp;
        for (const auto &p : thread.second)
        {
            Pack3Int(p, thread.first, m_rank, data[offset++]);
            m_staticPts->GetCoordsByPID(p, tmp = data + offset);
            offset += m_spacedim;
        }
        globalCoords[thread.first] = data;
    }
    UpdateMobileCoords(globalCoords);
}

void EvaluatePoints::PartitionLocalPoints(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    std::set<int> &notfound)
{
    if (m_isHomogeneous1D && m_columnRank != 0)
    {
        return;
    }
    // not repeat in local points
    for (auto &p : m_mobilePts)
    {
        Array<OneD, NekDouble> locCoords = p.second->GetLocalCoords();
        Array<OneD, NekDouble> gloCoords = p.second->GetGlobalCoords();
        // Determine the expansion and local coordinates
        if (m_isHomogeneous1D)
        {
            p.second->m_eId = pFields[0]->GetPlane(0)->GetExpIndex(
                gloCoords, locCoords, 0., false, p.second->m_eId);
        }
        else
        {
            p.second->m_eId = pFields[0]->GetExpIndex(gloCoords, locCoords, 0.,
                                                      false, p.second->m_eId);
        }
        if (p.second->m_eId < 0)
        {
            notfound.insert(p.first);
        }
    }
}

/***
 * Do global partition if mobile points are not found locally
 * Update the communication pair
 ***/
void EvaluatePoints::PartitionExchangeNonlocalPoints(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    std::set<int> &notfound)
{
    // Determine the unique process responsible for each history point
    // For points on a partition boundary, must select a single process
    int nExchange = 0, nTot = 0;
    Array<OneD, NekDouble> exchangeInfo;
    std::vector<int> npts(m_comm->GetRowComm()->GetSize(), 0);
    if (m_columnRank == 0 || !m_isHomogeneous1D)
    {
        npts[m_rowRank] = notfound.size();
        m_comm->GetRowComm()->AllReduce(npts, LibUtilities::ReduceMax);
        for (auto v : npts)
        {
            nTot += v;
        }
    }
    m_comm->AllReduce(nTot, LibUtilities::ReduceMax);
    if (nTot == 0)
    {
        return;
    }
    if (m_columnRank == 0 || !m_isHomogeneous1D)
    {
        // spread unfound points to all Row threads
        int nOffset = 0;
        for (int i = 0; i < m_rowRank; ++i)
        {
            nOffset += npts[i];
        }
        int nDataPackage = 1 + m_spacedim;
        int count        = nOffset * nDataPackage;
        Array<OneD, NekDouble> data(nTot * nDataPackage, 0.);
        for (auto &p : notfound)
        {
            Pack2Int(p, m_mobilePts[p]->m_sRank, data[count++]);
            for (int i = 0; i < m_spacedim; ++i)
            {
                data[count++] = m_mobilePts[p]->GetGlobalCoords()[i];
            }
        }
        m_comm->GetRowComm()->AllReduce(data, LibUtilities::ReduceSum);
        // global search
        std::vector<int> eIds(nTot, -1);
        std::vector<int> procIds(nTot, -1);
        Array<OneD, NekDouble> localCoords(nTot * m_spacedim, 0.);
        for (int p = 0; p < nTot; ++p)
        {
            if (nOffset <= p && p < nOffset + npts[m_rowRank])
            {
                continue;
            }
            Array<OneD, NekDouble> gloCoords = data + p * nDataPackage + 1;
            Array<OneD, NekDouble> locCoords = localCoords + p * m_spacedim;
            // Determine the expansion and local coordinates
            if (m_isHomogeneous1D && m_columnRank == 0)
            {
                eIds[p] =
                    pFields[0]->GetPlane(0)->GetExpIndex(gloCoords, locCoords);
            }
            else
            {
                eIds[p] = pFields[0]->GetExpIndex(gloCoords, locCoords);
            }
            if (eIds[p] >= 0)
            {
                procIds[p] = m_rowRank;
            }
        }
        m_comm->GetRowComm()->AllReduce(procIds, LibUtilities::ReduceMax);
        for (int p = 0; p < nTot; ++p)
        {
            if (procIds[p] >= 0)
            {
                ++nExchange;
                if (procIds[p] == m_rowRank)
                {
                    int pId, sRank;
                    unPack2Int(pId, sRank, data[p * nDataPackage]);
                    m_mobilePts[pId] =
                        MemoryManager<MobilePoint>::AllocateSharedPtr(
                            m_spacedim, sRank, data + p * nDataPackage + 1);
                    m_mobilePts[pId]->SetLocalCoords(
                        eIds[p], localCoords + p * m_spacedim);
                    m_recvStatInfo[sRank] += 1;
                }
            }
        }
        for (int p = nOffset; p < nOffset + npts[m_rowRank]; ++p)
        {
            if (procIds[p] >= 0)
            {
                int pId, sRank;
                unPack2Int(pId, sRank, data[p * nDataPackage]);
                m_mobilePts.erase(pId);
                if (m_recvStatInfo[sRank] == 1)
                {
                    m_recvStatInfo.erase(sRank);
                }
                else
                {
                    m_recvStatInfo[sRank] -= 1;
                }
            }
        }
        // build exchangeInfo [pid, from thread rank, to thread rank]
        int accumulate = 0;
        exchangeInfo   = Array<OneD, NekDouble>(nExchange, 0.);
        for (int t = 0, p = 0; t < npts.size(); ++t)
        {
            for (int i = 0; i < npts[t]; ++i, ++p)
            {
                if (procIds[p] >= 0)
                {
                    int pId, sRank;
                    unPack2Int(pId, sRank, data[p * nDataPackage]);
                    Pack3Int(sRank, t, procIds[p], exchangeInfo[accumulate++]);
                }
            }
        }
    }
    m_comm->AllReduce(nExchange, LibUtilities::ReduceMax);
    if (nExchange == 0)
    {
        return;
    }
    // sync column communication
    if (m_comm->GetColumnComm()->GetSize() > 1)
    {
        if (m_columnRank > 0)
        {
            exchangeInfo = Array<OneD, NekDouble>(nExchange, 0.);
        }
        m_comm->GetColumnComm()->AllReduce(exchangeInfo,
                                           LibUtilities::ReduceSum);
    }
    // update m_recvMobInfo
    for (int p = 0; p < nExchange; ++p)
    {
        int sRank, fromRank, toRank;
        unPack3Int(sRank, fromRank, toRank, exchangeInfo[p]);
        if (sRank == m_rank)
        {
            if (m_recvMobInfo[fromRank] == 1)
            {
                m_recvMobInfo.erase(fromRank);
            }
            else
            {
                m_recvMobInfo[fromRank] -= 1;
            }
            m_recvMobInfo[toRank] += 1;
        }
    }
}

void EvaluatePoints::SyncColumnComm()
{
    if (!m_isHomogeneous1D)
    {
        return;
    }
    if (m_columnRank == 0)
    {
        for (const auto &p : m_mobilePts)
        {
            p.second->m_local[2] = p.second->m_global[2];
        }
    }
    else
    {
        m_mobilePts.clear();
    }
    if (m_comm->GetColumnComm()->GetSize() == 1)
    {
        return;
    }
    int nPts = m_mobilePts.size();
    m_comm->GetColumnComm()->AllReduce(nPts, LibUtilities::ReduceMax);
    int nPackageSize = 1 + m_spacedim;
    Array<OneD, NekDouble> data(nPts * nPackageSize, 0.), tmp;
    if (m_columnRank == 0)
    {
        int count = 0;
        for (const auto &p : m_mobilePts)
        {
            Pack2Int(p.first, p.second->m_eId, data[count++]);
            Vmath::Vcopy(3, p.second->m_local, 1, tmp = data + count, 1);
            count += 3;
        }
    }
    m_comm->GetColumnComm()->AllReduce(data, LibUtilities::ReduceSum);
    if (m_columnRank != 0)
    {
        int offset = 0;
        for (int p = 0; p < nPts; ++p)
        {
            int pId, eId;
            Array<OneD, NekDouble> globalCoords(m_spacedim, 0.);
            unPack2Int(pId, eId, data[offset++]);
            m_mobilePts[pId] = MemoryManager<MobilePoint>::AllocateSharedPtr(
                m_spacedim, 0, globalCoords);
            m_mobilePts[pId]->SetLocalCoords(eId, data + offset);
            offset += m_spacedim;
        }
    }
}

NekDouble EvaluatePoints::GetDefaultValue(const int i,
                                          const Array<OneD, const NekDouble> x,
                                          const NekDouble time)
{
    NekDouble z = m_spacedim == 2 ? 0 : x[2];
    if (m_defaultValueFunction.find(i) == m_defaultValueFunction.end())
    {
        return 0.;
    }
    else
    {
        return m_defaultValueFunction[i]->Evaluate(x[0], x[1], z, time);
    }
}

void EvaluatePoints::UpdateMobileCoords(
    std::map<int, Array<OneD, NekDouble>> &globalCoords)
{
    if (m_columnRank == 0 || !m_isHomogeneous1D)
    {
        for (auto &iter : globalCoords)
        {
            if (iter.first != m_rank)
            {
                m_comm->Send(iter.first, iter.second);
            }
        }
        for (const auto &iter : m_recvStatInfo)
        {
            Array<OneD, NekDouble> data;
            if (iter.first != m_rank)
            {
                data = Array<OneD, NekDouble>(iter.second * (m_spacedim + 1));
                m_comm->Recv(iter.first, data);
            }
            else
            {
                data = globalCoords[iter.first];
            }
            for (int offset = 0; offset < data.size(); offset += m_spacedim)
            {
                int pid, mobt, stat;
                unPack3Int(pid, mobt, stat, data[offset++]);
                m_mobilePts[pid]->SetGlobalCoords(data + offset);
            }
        }
    }
}

void EvaluatePoints::CopyStaticPtsToMobile()
{
    ASSERTL0(m_columnRank == 0 || !m_isHomogeneous1D,
             "CopyStaticPtsToMobile shoule only be called by ColumnRank zero");
    m_mobilePts.clear();
    for (int p = 0; p < m_staticPts->GetTotPoints(); ++p)
    {
        Array<OneD, NekDouble> data(m_spacedim, 0.);
        m_staticPts->GetCoords(p, data);
        m_mobilePts[m_staticPts->LocalToGlobal(p)] =
            MemoryManager<MobilePoint>::AllocateSharedPtr(m_spacedim, m_rank,
                                                          data);
    }
}

void EvaluatePoints::CopyMobilePtsToStatic()
{
    ASSERTL0(m_columnRank == 0 || !m_isHomogeneous1D,
             "CopyStaticPtsToMobile shoule only be called by ColumnRank zero");
    m_staticPts->ReSize(m_mobilePts.size());
    int id = 0;
    for (const auto &p : m_mobilePts)
    {
        m_staticPts->AssignPoint(id, p.first, p.second->GetGlobalCoords());
    }
}

} // namespace Nektar::SolverUtils
