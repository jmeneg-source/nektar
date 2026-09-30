////////////////////////////////////////////////////////////////////////////////
//
//  File: SourcePoint.hpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: class and methods of curvature sampling point
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESH_OCTREE_SOURCEPOINT_H
#define NEKMESH_OCTREE_SOURCEPOINT_H

#include <array>

#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar::NekMesh
{

enum SPType
{
    eCBoundary,  // on a curved boundary
    eRCBoundary, // on a curved boundary with refinement
    ePBoundary,  // on a planar boundary (R=inf)
    eSrcPoint    // source point
};

/**
 * @brief base class of sizing point for octree construction
 *        these carry information the octree needs and have various types
 */
class SPBase
{
public:
    friend class MemoryManager<SPBase>;

    SPBase(std::array<NekDouble, 3> l)
    {
        m_loc = l;
    }

    virtual ~SPBase()
    {
    }

    SPType GetType()
    {
        return m_type;
    }

    std::array<NekDouble, 3> GetLoc()
    {
        return m_loc;
    }

    virtual NekDouble GetDelta()
    {
        return 0.0;
    }

    virtual NekDouble GetRDelta()
    {
        return 0.0;
    }

    virtual void SetDelta([[maybe_unused]] NekDouble i)
    {
    }

    virtual void GetCAD([[maybe_unused]] int &surf,
                        [[maybe_unused]] std::array<NekDouble, 2> &uv)
    {
    }

    bool HasDelta()
    {
        bool ret;
        if (m_type == eCBoundary || m_type == eRCBoundary)
        {
            ret = true;
        }
        else
        {
            ret = false;
        }
        return ret;
    }

    bool HasRDelta()
    {
        bool ret;
        if (m_type == eRCBoundary || m_type == eSrcPoint)
        {
            ret = true;
        }
        else
        {
            ret = false;
        }
        return ret;
    }

    bool Isboundary()
    {
        bool ret;
        if (m_type == eCBoundary || m_type == ePBoundary ||
            m_type == eRCBoundary)
        {
            ret = true;
        }
        else
        {
            ret = false;
        }
        return ret;
    }

protected:
    /// type
    SPType m_type;
    /// x,y,z location
    std::array<NekDouble, 3> m_loc;
};

typedef std::shared_ptr<SPBase> SPBaseSharedPtr;

/**
 * @brief class for a curvature based samlping Point
 */
class CPoint : public SPBase
{
public:
    friend class MemoryManager<CPoint>;

    /**
     * @brief constructor for a valid point (has radius of curvature)
     */
    CPoint(int i, std::array<NekDouble, 2> uv, std::array<NekDouble, 3> l,
           NekDouble d)
        : SPBase(l), sid(i), m_uv(uv), m_delta(d)
    {
        m_type = eCBoundary;
    }

    // Overloaded constructor with refined delta value
    CPoint(int i, std::array<NekDouble, 2> uv, std::array<NekDouble, 3> l,
           NekDouble d, NekDouble rd)
        : SPBase(l), sid(i), m_uv(uv), m_delta(d), m_rdelta(rd)
    {
        m_type = eRCBoundary;
    }

    ~CPoint() override {};

    /**
     * @brief get mesh spacing paramter
     */
    NekDouble GetDelta() override
    {
        return m_delta;
    }

    /**
     * @brief get mesh refinement spacing paramter
     */
    NekDouble GetRDelta() override
    {
        return m_rdelta;
    }

    /**
     * @brief gets the corresponding cad information for the point
     */
    void GetCAD(int &surf, std::array<NekDouble, 2> &uv) override
    {
        surf = sid;
        uv   = m_uv;
    }

    void SetRDelta(NekDouble i)
    {
        m_rdelta = i;
    }

    void SetDelta(NekDouble i) override
    {
        m_delta = i;
    }

private:
    /// surf id
    int sid;
    /// uv coord on surf
    std::array<NekDouble, 2> m_uv;
    NekDouble m_ti;
    /// delta parameter
    NekDouble m_delta;
    /// refinement delta parameter
    NekDouble m_rdelta;
};
typedef std::shared_ptr<CPoint> CPointSharedPtr;

/**
 * @brief class for a planar boundary based samlping Point
 */
class BPoint : public SPBase
{
public:
    friend class MemoryManager<BPoint>;

    /**
     * @brief constructor for a boundary point without delta
     */
    BPoint(int i, std::array<NekDouble, 2> uv, std::array<NekDouble, 3> l)
        : SPBase(l), sid(i), m_uv(uv)
    {
        m_type = ePBoundary;
    }

    ~BPoint() override {};

    /**
     * @brief gets the corresponding cad information for the point
     */
    void GetCAD(int &surf, std::array<NekDouble, 2> &uv) override
    {
        surf = sid;
        uv   = m_uv;
    }

    NekDouble GetDelta() override
    {
        NEKERROR(ErrorUtil::efatal, "Cannot retrieve delta from this type");
        return 0.0;
    }

    void SetDelta([[maybe_unused]] NekDouble i) override
    {
        NEKERROR(ErrorUtil::efatal,
                 "Cannot assign refinement delta to this type");
    }

    NekDouble GetRDelta() override
    {
        NEKERROR(ErrorUtil::efatal,
                 "Cannot retrieve refinment delta from this type");
        return 0.0;
    }

    void SetRDelta([[maybe_unused]] NekDouble i)
    {
        NEKERROR(ErrorUtil::efatal, "Cannot assign delta to this type");
    }

private:
    /// surf id
    int sid;
    /// uv coord on surf
    std::array<NekDouble, 2> m_uv;
    NekDouble m_ti;
};
typedef std::shared_ptr<BPoint> BPointSharedPtr;

/**
 * @brief class for a general source point
 */
class SrcPoint : public SPBase
{
public:
    friend class MemoryManager<SrcPoint>;

    /**
     * @brief constructor for a boundary point without delta
     */
    SrcPoint(std::array<NekDouble, 3> l, NekDouble d) : SPBase(l), m_rdelta(d)
    {
        m_type = eSrcPoint;
    }

    ~SrcPoint() override {};

    /**
     * @brief get mesh spacing paramter
     */
    NekDouble GetRDelta() override
    {
        return m_rdelta;
    }

    void SetRDelta(NekDouble i)
    {
        m_rdelta = i;
    }

    NekDouble GetDelta() override
    {
        NEKERROR(ErrorUtil::efatal, "Cannot retrieve delta from this type");
        return 0.0;
    }

    void SetDelta([[maybe_unused]] NekDouble i) override
    {
        NEKERROR(ErrorUtil::efatal, "Cannot assign delta to this type");
    }

    void GetCAD([[maybe_unused]] int &surf,
                [[maybe_unused]] std::array<NekDouble, 2> &uv) override
    {
        NEKERROR(ErrorUtil::efatal, "Cannot retrieve CAD from this type")
    }

private:
    /// refinement delta parameter
    NekDouble m_rdelta;
};
typedef std::shared_ptr<SrcPoint> SrcPointSharedPtr;

} // namespace Nektar::NekMesh

#endif
