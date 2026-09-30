///////////////////////////////////////////////////////////////////////////////
//
// File: ObjectPool.hpp
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
// Description: Fast geometry object allocator using boost::fast_pool_allocator.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_OBJECTPOOL_HPP
#define NEKTAR_LIBUTILITIES_OBJECTPOOL_HPP

#include <iostream>
#include <memory>
#include <vector>

#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>

namespace Nektar
{

template <typename DataType>
using PoolAllocator =
    boost::fast_pool_allocator<DataType,
                               boost::default_user_allocator_new_delete,
                               boost::details::pool::null_mutex>;

/**
 * @brief Generic object pool allocator/deallocator.
 *
 * This class provides an allocator based on the boost::fast_pool_allocator for
 * creating object pools.
 */
template <typename DataType> class ObjPoolManager
{
public:
    template <typename... Args> static DataType *Allocate(const Args &...args)
    {
        DataType *ptr = m_alloc.allocate();
        new (ptr) DataType(args...);
        return ptr;
    }

    static void Deallocate(DataType *ptr)
    {
        ptr->~DataType();
        return m_alloc.deallocate(ptr);
    }

    struct UniquePtrDeleter
    {
        void operator()(DataType *ptr) const
        {
            // ignore dealloc
            ObjPoolManager<DataType>::Deallocate(ptr);
        }
    };

    template <typename... Args>
    static std::unique_ptr<DataType, UniquePtrDeleter> AllocateUniquePtr(
        const Args &...args)
    {
        return std::unique_ptr<DataType, UniquePtrDeleter>(Allocate(args...));
    }

    static PoolAllocator<DataType> m_alloc;
};

template <typename T>
using unique_ptr_objpool =
    std::unique_ptr<T, typename ObjPoolManager<T>::UniquePtrDeleter>;

} // namespace Nektar

#endif
