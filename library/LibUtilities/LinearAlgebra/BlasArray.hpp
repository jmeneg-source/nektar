///////////////////////////////////////////////////////////////////////////////
//
// File: BlasArray.hpp
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
// Description: wrapper of functions around standard BLAS routines
// using Array's as calling arguments
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLASARRAY_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLASARRAY_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/TransF77.hpp>

// Translations for using Fortran version of blas
namespace Blas
{

extern "C"
{
    // -- BLAS Level 1:
    void F77NAME(daxpy)(const int &n, const double &alpha, const double *x,
                        const int &incx, const double *y, const int &incy);
}

/// \brief  BLAS level 1: y = alpha \a x plus \a y
static inline void Daxpy(const int &n, const double &alpha,
                         const Nektar::Array<Nektar::OneD, const double> &x,
                         const int &incx,
                         Nektar::Array<Nektar::OneD, double> &y,
                         const int &incy)
{
    ASSERTL1(static_cast<unsigned int>(n * incx) <= x.size() + x.GetOffset(),
             "Array out of bounds");
    ASSERTL1(static_cast<unsigned int>(n * incy) <= y.size() + y.GetOffset(),
             "Array out of bounds");

    F77NAME(daxpy)(n, alpha, &x[0], incx, &y[0], incy);
}

} // namespace Blas
#endif // NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLASARRAY_HPP
