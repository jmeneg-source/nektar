///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrix.hpp
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
// Description: Python wrapper for NekMatrix.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_PYTHON_LINEARALGERBA_NEKMATRIX_HPP
#define NEKTAR_LIBUTILITIES_PYTHON_LINEARALGERBA_NEKMATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <LibUtilities/Python/NekPyConfig.hpp>

namespace pybind11::detail
{

template <typename Type, typename T> struct standard_nekmatrix_caster
{
public:
    PYBIND11_TYPE_CASTER(Type, const_name("NekMatrix<T>"));

    bool load(handle src, bool)
    {
        // Perform some checks: the source should be an ndarray.
        if (!array_t<T>::check_(src))
        {
            return false;
        }

        // The array should be a c-style, contiguous array of type T.
        auto buf = array_t<T, array::c_style | array::forcecast>::ensure(src);
        if (!buf)
        {
            return false;
        }

        // There should be two dimensions.
        auto dims = buf.ndim();
        if (dims != 2)
        {
            return false;
        }

        // Copy data across from the Python array to C++.
        std::vector<size_t> shape = {buf.shape()[0], buf.shape()[1]};
        value = Nektar::NekMatrix<T, Nektar::StandardMatrixTag>(
            shape[0], shape[1], buf.data(), Nektar::eFULL, buf.data());

        return true;
    }

    // Conversion part 2 (C++ -> Python)
    static handle cast(
        const Nektar::NekMatrix<T, Nektar::StandardMatrixTag> &src,
        return_value_policy, handle)
    {
        using NMat = Nektar::NekMatrix<T, Nektar::StandardMatrixTag>;

        // Construct a new wrapper matrix to hold onto the data. Assign a
        // destructor so that the wrapper is cleaned up when the Python object
        // is deallocated.
        capsule c(new NMat(src.GetRows(), src.GetColumns(), src.GetPtr(),
                           Nektar::eWrapper),
                  [](void *ptr) {
                      NMat *mat = (NMat *)ptr;
                      delete mat;
                  });

        // Create the NumPy array, passing the capsule. When we go out of scope,
        // c's reference count will have been reduced by 1, but array increases
        // the reference count when it assigns the base to the array.
        array a({src.GetRows(), src.GetColumns()},
                {sizeof(T), src.GetRows() * sizeof(T)}, src.GetRawPtr(), c);

        // This causes the array to be released without decreasing its reference
        // count, which we do since if we just returned a, then the array would
        // be deallocated immediately when this function returns.
        return a.release();
    }
};

template <typename Type>
struct type_caster<Nektar::NekMatrix<Type, Nektar::StandardMatrixTag>>
    : standard_nekmatrix_caster<
          Nektar::NekMatrix<Type, Nektar::StandardMatrixTag>, Type>
{
};

} // namespace pybind11::detail

#endif
