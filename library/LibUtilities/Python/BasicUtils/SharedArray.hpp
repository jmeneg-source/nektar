///////////////////////////////////////////////////////////////////////////////
//
// File: SharedArray.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Python wrapper for ShareArray.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_PYTHON_BASICUTILS_SHAREDARRAY_HPP
#define NEKTAR_LIBUTILITIES_PYTHON_BASICUTILS_SHAREDARRAY_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/Python/NekPyConfig.hpp>

#include <type_traits>

namespace pybind11::detail
{

/// Template utility to determine whether @tparam T is a Nektar::Array<OneD, .>.
template <typename T> struct is_nekarray_oned : std::false_type
{
};

template <typename T>
struct is_nekarray_oned<Nektar::Array<Nektar::OneD, T>> : std::true_type
{
};

template <typename T>
struct is_nekarray_oned<const Nektar::Array<Nektar::OneD, T>> : std::true_type
{
};

/**
 * @brief Convert for Array<OneD, T> to Python list of objects for numeric
 * types, self-array types and shared_ptr types.
 *
 * @tparam Type  The overall type Array<OneD, T>
 * @tparam T     The data type T
 */
template <typename Type, typename T> struct nekarray_caster
{
public:
    PYBIND11_TYPE_CASTER(Type, const_name("Array<OneD, ") +
                                   handle_type_name<T>::name + const_name(">"));

    bool load(handle src, bool enable)
    {
        return load_impl(src, enable);
    }

    static handle cast(Nektar::Array<Nektar::OneD, T> const &arr,
                       return_value_policy policy, handle parent)
    {
        return cast_impl(arr, policy, parent);
    }

    template <typename U                                      = T,
              std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
    bool load_impl(handle src, bool)
    {
        if (!array_t<U>::check_(src))
        {
            return false;
        }

        auto arr = array_t<U, array::c_style | array::forcecast>::ensure(src);
        if (!arr)
        {
            return false;
        }

        if (arr.ndim() != 1)
        {
            return false;
        }

        using nonconst_t = typename std::remove_const<U>::type;
        value            = Nektar::Array<Nektar::OneD, T>(arr.shape(0),
                                                          (nonconst_t *)arr.data(),
                                                          (void *)src.ptr(), &decrement);

        // We increase the refcount on src so that even if on the Python side we
        // go out of scope, the data will still exist in at least one reference
        // until it is no longer used on the C++ side.
        src.inc_ref();

        return true;
    }

    static void decrement(void *objPtr)
    {
        if (!Py_IsInitialized())
        {
            // In deinitialisation phase, reference counters appear to not be
            // terribly robust; decremementing counters here can lead to
            // segfaults during program exit in some cases.
            return;
        }

        // Otherwise decrement reference counter.
        handle obj((PyObject *)objPtr);
        obj.dec_ref();
    }

    template <typename U                                      = T,
              std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
    static handle cast_impl(Nektar::Array<Nektar::OneD, T> const &arr,
                            return_value_policy, handle)
    {
        // Create a Python capsule to hold a pointer that contains a lightweight
        // copy of arr. Uhat way we guarantee Python will still have access to
        // the memory allocated inside arr even if arr is deallocated in C++.
        capsule c(new Nektar::Array<Nektar::OneD, U>(arr), [](void *ptr) {
            Nektar::Array<Nektar::OneD, U> *tmp =
                (Nektar::Array<Nektar::OneD, U> *)ptr;
            delete tmp;
        });

        // Create the NumPy array, passing the capsule. When we go out of scope,
        // c's reference count will have been reduced by 1, but array increases
        // the reference count when it assigns the base to the array.
        array_t<U> array({arr.size()}, {}, arr.data(), c);

        // This causes the array to be released without decreasing its reference
        // count, otherwise the array would be deallocated immediately when this
        // function returns.
        return array.release();
    }

    template <typename U = T,
              std::enable_if_t<
                  is_shared_ptr<typename std::remove_const<U>::type>::value ||
                      is_nekarray_oned<U>::value,
                  bool> = true>
    bool load_impl(handle src, bool)
    {
        if (!py::isinstance<py::list>(src))
        {
            return false;
        }

        py::list l               = py::cast<py::list>(src);
        const std::size_t nItems = l.size();

        using nonconst_t = typename std::remove_const<U>::type;
        auto tmparr      = Nektar::Array<Nektar::OneD, nonconst_t>(nItems);

        // We'll need to construct a temporary vector to hold each item in the
        // list.
        try
        {
            for (std::size_t i = 0; i < nItems; ++i)
            {
                tmparr[i] = py::cast<nonconst_t>(l[i]);
            }
        }
        catch (...)
        {
            return false;
        }

        value = tmparr;

        return true;
    }

    template <typename U = T,
              std::enable_if_t<
                  is_shared_ptr<typename std::remove_const<U>::type>::value ||
                      is_nekarray_oned<U>::value,
                  bool> = true>
    static handle cast_impl(Nektar::Array<Nektar::OneD, U> const &arr,
                            return_value_policy, handle)
    {
        py::list tmp;

        for (std::size_t i = 0; i < arr.size(); ++i)
        {
            tmp.append(arr[i]);
        }

        return tmp.release();
    }
};

template <typename Type>
struct type_caster<Nektar::Array<Nektar::OneD, Type>>
    : nekarray_caster<Nektar::Array<Nektar::OneD, Type>, Type>
{
};

} // namespace pybind11::detail

#endif
