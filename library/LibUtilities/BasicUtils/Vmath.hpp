///////////////////////////////////////////////////////////////////////////////
//
// File: Vmath.hpp
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
// Description: Collection of templated functions for vector mathematics
//
// Note: For those unfamiliar with the vector routines notation, it is
//       reverse polish notation (RPN).  For example:
//
//       In the function "Vvtvp()", it is "Vvt" means vector vector times,
//       which in infix notation is "v * v".  So "Vvtvp" is:
//
//       RPN:    vector vector times vector plus
//       Infix:  ( vector * vector ) + vector
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATH_HPP
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATH_HPP

#include <LibUtilities/LinearAlgebra/Blas.hpp>

namespace Vmath
{
/***************** Math routines  ***************/

/// \brief Fill a vector with a constant value
template <class T> inline void Fill(int n, const T alpha, T *x, const int incx)
{
    while (n--)
    {
        *x = alpha;
        x += incx;
    }
}

/// \brief Generates a number from ~Normal(0,1)
template <class T> T ran2(long *idum);

/// \brief Fills a vector with white noise.
template <class T>
void FillWhiteNoise(int n, const T eps, T *x, const int incx, int seed = 9999);

/// \brief Multiply vector z = x*y
template <class T>
inline void Vmul(int n, const T *x, const int incx, const T *y, const int incy,
                 T *z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1 && incz == 1)
    {
        while (--n)
        {
            *z = (*x) * (*y);
            ++x;
            ++y;
            ++z;
        }
    }
    else
    {
        while (--n)
        {
            *z = (*x) * (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
}

/// \brief Scalar multiply  y = alpha*x
template <class T>
inline void Smul(int n, const T alpha, const T *x, const int incx, T *y,
                 const int incy)
{
    ++n;
    if (incx == 1 && incy == 1)
    {
        while (--n)
        {
            *y = alpha * (*x);
            ++x;
            ++y;
        }
    }
    else
    {
        while (--n)
        {
            *y = alpha * (*x);
            x += incx;
            y += incy;
        }
    }
}

/// \brief Multiply vector z = x/y
template <class T>
inline void Vdiv(int n, const T *x, const int incx, const T *y, const int incy,
                 T *z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1)
    {
        while (--n)
        {
            *z = (*x) / (*y);
            ++x;
            ++y;
            ++z;
        }
    }
    else
    {
        while (--n)
        {
            *z = (*x) / (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
}

/// \brief Scalar multiply  y = alpha/x
template <class T>
inline void Sdiv(int n, const T alpha, const T *x, const int incx, T *y,
                 const int incy)
{
    ++n;
    if (incx == 1 && incy == 1)
    {
        while (--n)
        {
            *y = alpha / (*x);
            ++x;
            ++y;
        }
    }
    else
    {
        while (--n)
        {
            *y = alpha / (*x);
            x += incx;
            y += incy;
        }
    }
}

/// \brief Add vector z = x+y
template <class T>
inline void Vadd(int n, const T *x, const int incx, const T *y, const int incy,
                 T *z, const int incz)
{
    while (n--)
    {
        *z = (*x) + (*y);
        x += incx;
        y += incy;
        z += incz;
    }
}

/// \brief Add vector y = alpha + x
template <class T>
inline void Sadd(int n, const T alpha, const T *x, const int incx, T *y,
                 const int incy)
{
    ++n;
    if (incx == 1 && incy == 1)
    {
        while (--n)
        {
            *y = alpha + (*x);
            ++x;
            ++y;
        }
    }
    else
    {
        while (--n)
        {
            *y = alpha + (*x);
            x += incx;
            y += incy;
        }
    }
}

/// \brief Subtract vector z = x-y
template <class T>
inline void Vsub(int n, const T *x, const int incx, const T *y, const int incy,
                 T *z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1 && incz == 1)
    {
        while (--n)
        {
            *z = (*x) - (*y);
            ++x;
            ++y;
            ++z;
        }
    }
    else
    {
        while (--n)
        {
            *z = (*x) - (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
}

/// \brief Substract vector y = alpha - x
template <class T>
inline void Ssub(int n, const T alpha, const T *x, const int incx, T *y,
                 const int incy)
{
    ++n;
    if (incx == 1 && incy == 1)
    {
        while (--n)
        {
            *y = alpha - (*x);
            ++x;
            ++y;
        }
    }
    else
    {
        while (--n)
        {
            *y = alpha - (*x);
            x += incx;
            y += incy;
        }
    }
}

/// \brief Zero vector
template <class T> inline void Zero(int n, T *x, const int incx)
{
    if (incx == 1)
    {
        std::memset(x, '\0', n * sizeof(T));
    }
    else
    {
        T zero = 0;
        ++n;
        while (--n)
        {
            *x = zero;
            x += incx;
        }
    }
}

/// \brief Negate x = -x
template <class T> inline void Neg(int n, T *x, const int incx)
{
    while (n--)
    {
        *x = -(*x);
        x += incx;
    }
}

/// \brief log y = log(x)
template <class T>
inline void Vlog(int n, const T *x, const int incx, T *y, const int incy)
{
    while (n--)
    {
        *y = log(*x);
        x += incx;
        y += incy;
    }
}

/// \brief exp y = exp(x)
template <class T>
inline void Vexp(int n, const T *x, const int incx, T *y, const int incy)
{
    while (n--)
    {
        *y = exp(*x);
        x += incx;
        y += incy;
    }
}

/// \brief pow y = pow(x, f)
template <class T>
inline void Vpow(int n, const T *x, const int incx, const T f, T *y,
                 const int incy)
{
    while (n--)
    {
        *y = pow(*x, f);
        x += incx;
        y += incy;
    }
}

/// \brief sqrt y = sqrt(x)
template <class T>
inline void Vsqrt(int n, const T *x, const int incx, T *y, const int incy)
{
    while (n--)
    {
        *y = sqrt(*x);
        x += incx;
        y += incy;
    }
}

/// \brief vabs: y = |x|
template <class T>
inline void Vabs(int n, const T *x, const int incx, T *y, const int incy)
{
    while (n--)
    {
        *y = (*x > 0) ? *x : -(*x);
        x += incx;
        y += incy;
    }
}

/********** Triad  routines  ***********************/

/// \brief  vvtvp (vector times vector plus vector): z = w*x + y
template <class T>
inline void Vvtvp(int n, const T *w, const int incw, const T *x, const int incx,
                  const T *y, const int incy, T *z, const int incz)
{
    while (n--)
    {
        *z = (*w) * (*x) + (*y);
        w += incw;
        x += incx;
        y += incy;
        z += incz;
    }
}

/// \brief vvtvm (vector times vector minus vector): z = w*x - y
template <class T>
inline void Vvtvm(int n, const T *w, const int incw, const T *x, const int incx,
                  const T *y, const int incy, T *z, const int incz)
{
    while (n--)
    {
        *z = (*w) * (*x) - (*y);
        w += incw;
        x += incx;
        y += incy;
        z += incz;
    }
}

/// \brief  Svtvp (scalar times vector plus vector): z = alpha*x + y
template <class T>
inline void Svtvp(int n, const T alpha, const T *x, const int incx, const T *y,
                  const int incy, T *z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1 && incz == 1)
    {
        while (--n)
        {
            *z = alpha * (*x) + (*y);
            ++x;
            ++y;
            ++z;
        }
    }
    else
    {
        while (--n)
        {
            *z = alpha * (*x) + (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
}

/// \brief  Svtvm (scalar times vector minus vector): z = alpha*x - y
template <class T>
inline void Svtvm(int n, const T alpha, const T *x, const int incx, const T *y,
                  const int incy, T *z, const int incz)
{
    while (n--)
    {
        *z = alpha * (*x) - (*y);
        x += incx;
        y += incy;
        z += incz;
    }
}

/// \brief  vvtvvtp (vector times vector plus vector times vector):
// z = v*w + x*y
template <class T>
inline void Vvtvvtp(int n, const T *v, int incv, const T *w, int incw,
                    const T *x, int incx, const T *y, int incy, T *z, int incz)
{
    while (n--)
    {
        *z = (*v) * (*w) + (*x) * (*y);
        v += incv;
        w += incw;
        x += incx;
        y += incy;
        z += incz;
    }
}

/// \brief  vvtvvtm (vector times vector minus vector times vector):
// z = v*w - x*y
template <class T>
inline void Vvtvvtm(int n, const T *v, int incv, const T *w, int incw,
                    const T *x, int incx, const T *y, int incy, T *z, int incz)
{
    while (n--)
    {
        *z = (*v) * (*w) - (*x) * (*y);
        v += incv;
        w += incw;
        x += incx;
        y += incy;
        z += incz;
    }
}

/// \brief  Svtsvtp (scalar times vector plus scalar times vector):
// z = alpha*x + beta*y
template <class T>
inline void Svtsvtp(int n, const T alpha, const T *x, int incx, const T beta,
                    const T *y, int incy, T *z, int incz)
{
    while (n--)
    {
        *z = alpha * (*x) + beta * (*y);
        x += incx;
        y += incy;
        z += incz;
    }
}

/// \brief  Vstvpp (scalar times vector plus vector plus vector):
// z = alpha*w + x + y
template <class T>
inline void Vstvpp(int n, const T alpha, const T *v, int incv, const T *w,
                   int incw, const T *x, int incx, T *z, int incz)
{
    while (n--)
    {
        *z = alpha * (*v) + (*w) + (*x);
        v += incv;
        w += incw;
        x += incx;
        z += incz;
    }
}

/************ Misc routine from Veclib (and extras)  ************/

/// \brief Gather vector z[i] = x[y[i]]
template <class T, class I,
          typename = typename std::enable_if<std::is_floating_point_v<T> &&
                                             std::is_integral_v<I>>::type>
inline void Gathr(I n, const T *x, const I *y, T *z)
{
    while (n--)
    {
        *z++ = *(x + *y++);
    }
    return;
}

/// \brief Gather vector z[i] = sign*x[y[i]]
template <class T>
inline void Gathr(int n, const T sign, const T *x, const int *y, T *z)
{
    while (n--)
    {
        *z++ = sign * (*(x + *y++));
    }
    return;
}

/// \brief Gather vector z[i] = sign[i]*x[y[i]]
template <class T>
inline void Gathr(int n, const T *sign, const T *x, const int *y, T *z)
{
    while (n--)
    {
        *z++ = *(sign++) * (*(x + *y++));
    }
    return;
}

/// \brief Scatter vector z[y[i]] = x[i]
template <class T> inline void Scatr(int n, const T *x, const int *y, T *z)
{
    while (n--)
    {
        *(z + *(y++)) = *(x++);
    }
}

/// \brief Scatter vector z[y[i]] = sign*x[i]
template <class T>
inline void Scatr(int n, const T sign, const T *x, const int *y, T *z)
{
    while (n--)
    {
        *(z + *(y++)) = sign * (*(x++));
    }
}

/// \brief Scatter vector z[y[i]] = sign[i]*x[i]
template <class T>
inline void Scatr(int n, const T *sign, const T *x, const int *y, T *z)
{
    while (n--)
    {
        if (*sign)
        {
            *(z + *(y++)) = *(sign++) * (*(x++));
        }
        else
        {
            x++;
            y++;
            sign++;
        }
    }
}

/// \brief Assemble z[y[i]] += x[i]; z should be zero'd first
template <class T> inline void Assmb(int n, const T *x, const int *y, T *z)
{
    while (n--)
    {
        *(z + *(y++)) += *(x++);
    }
}

/// \brief Assemble z[y[i]] += sign*x[i]; z should be zero'd first
template <class T>
inline void Assmb(int n, const T sign, const T *x, const int *y, T *z)
{
    while (n--)
    {
        *(z + *(y++)) += sign * (*(x++));
    }
}

/// \brief Assemble z[y[i]] += sign[i]*x[i]; z should be zero'd first
template <class T>
inline void Assmb(int n, const T *sign, const T *x, const int *y, T *z)
{
    while (n--)
    {
        *(z + *(y++)) += *(sign++) * (*(x++));
    }
}

/************* Reduction routines  *****************/

/// \brief Subtract return sum(x)
template <class T> inline T Vsum(int n, const T *x, const int incx)
{

    T sum = 0;

    while (n--)
    {
        sum += (*x);
        x += incx;
    }

    return sum;
}

/// \brief Return the index of the maximum element in x
template <class T> inline int Imax(int n, const T *x, const int incx)
{

    int i, indx = (n > 0) ? 0 : -1;
    T xmax = *x;

    for (i = 0; i < n; i++)
    {
        if (*x > xmax)
        {
            xmax = *x;
            indx = i;
        }
        x += incx;
    }

    return indx;
}

/// \brief Return the maximum element in x -- called vmax to avoid
/// conflict with max
template <class T> inline T Vmax(int n, const T *x, const int incx)
{

    T xmax = *x;

    while (n--)
    {
        if (*x > xmax)
        {
            xmax = *x;
        }
        x += incx;
    }

    return xmax;
}

/// \brief Return the index of the maximum absolute element in x
template <class T> inline int Iamax(int n, const T *x, const int incx)
{

    int i, indx = (n > 0) ? 0 : -1;
    T xmax = *x;
    T xm;

    for (i = 0; i < n; i++)
    {
        xm = (*x > 0) ? *x : -*x;
        if (xm > xmax)
        {
            xmax = xm;
            indx = i;
        }
        x += incx;
    }

    return indx;
}

/// \brief Return the maximum absolute element in x
/// called vamax to avoid conflict with max
template <class T> inline T Vamax(int n, const T *x, const int incx)
{

    T xmax = *x;
    T xm;

    while (n--)
    {
        xm = (*x > 0) ? *x : -*x;
        if (xm > xmax)
        {
            xmax = xm;
        }
        x += incx;
    }
    return xmax;
}

/// \brief Return the index of the minimum element in x
template <class T> inline int Imin(int n, const T *x, const int incx)
{

    int i, indx = (n > 0) ? 0 : -1;
    T xmin = *x;

    for (i = 0; i < n; i++)
    {
        if (*x < xmin)
        {
            xmin = *x;
            indx = i;
        }
        x += incx;
    }

    return indx;
}

/// \brief Return the minimum element in x - called vmin to avoid
/// conflict with min
template <class T> inline T Vmin(int n, const T *x, const int incx)
{

    T xmin = *x;

    while (n--)
    {
        if (*x < xmin)
        {
            xmin = *x;
        }
        x += incx;
    }

    return xmin;
}

/// \brief Return number of NaN elements of x
template <class T> inline int Nnan(int n, const T *x, const int incx)
{

    int nNan = 0;

    while (n--)
    {
        if (*x != *x)
        {
            nNan++;
        }
        x += incx;
    }

    return nNan;
}

/// \brief dot product
template <class T> inline T Dot(int n, const T *w, const T *x)
{
    T sum = 0;

    while (n--)
    {
        sum += (*w) * (*x);
        ++w;
        ++x;
    }
    return sum;
}

/// \brief dot product
template <class T>
inline T Dot(int n, const T *w, const int incw, const T *x, const int incx)
{
    T sum = 0;

    while (n--)
    {
        sum += (*w) * (*x);
        w += incw;
        x += incx;
    }
    return sum;
}

/// \brief dot product
template <class T> inline T Dot2(int n, const T *w, const T *x, const int *y)
{
    T sum = 0;

    while (n--)
    {
        sum += (*y == 1 ? (*w) * (*x) : 0);
        ++w;
        ++x;
        ++y;
    }
    return sum;
}

/// \brief dot product
template <class T>
inline T Dot2(int n, const T *w, const int incw, const T *x, const int incx,
              const int *y, const int incy)
{
    T sum = 0;

    while (n--)
    {
        sum += (*y == 1 ? (*w) * (*x) : 0.0);
        w += incw;
        x += incx;
        y += incy;
    }
    return sum;
}

/********** Memory routines  ***********************/

// \brief copy one vector to another
template <class T>
inline void Vcopy(int n, const T *x, const int incx, T *y, const int incy)
{
    if (incx == 1 && incy == 1)
    {
        memcpy(y, x, n * sizeof(T));
    }
    else
    {
        while (n--)
        {
            *y = *x;
            x += incx;
            y += incy;
        }
    }
}

// \brief reverse the ordering of vector to another
template <class T>
inline void Reverse(int n, const T *x, const int incx, T *y, const int incy)
{
    int i;
    T store;

    // Perform element by element swaps in case x and y reference the same
    // array.
    int nloop = n / 2;

    // copy value in case of n is odd number
    y[nloop] = x[nloop];

    const T *x_end = x + (n - 1) * incx;
    T *y_end       = y + (n - 1) * incy;
    for (i = 0; i < nloop; ++i)
    {
        store  = *x_end;
        *y_end = *x;
        *y     = store;
        x += incx;
        y += incy;
        x_end -= incx;
        y_end -= incy;
    }
}

} // namespace Vmath
#endif // VECTORMATH_HPP
