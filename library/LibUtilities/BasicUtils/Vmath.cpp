///////////////////////////////////////////////////////////////////////////////
//
// File: Vmath.cpp
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
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Vmath.hpp>

namespace Vmath
{

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

/// \brief Generates a number from ~Normal(0,1)
template <class T> T ran2(long *idum)
/* ------------------------------------------------------------------------- *
 * Ran2 from NR 2e.  Returns a uniform random deviate between 0.0 &
 * 1.0 (exclusive of endpoints).  Call with idum a negative integer to
 * initialize; thereafter, do not alter idum between successive
 * deviates in a sequence.  RNMX should approximate the largest
 * floating value that is less than 1.
 * ------------------------------------------------------------------------- */
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy    = 0;
    static long iv[NTAB];
    T temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1)
        {
            *idum = 1;
        }
        else
        {
            *idum = -(*idum);
        }
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--)
        {
            k     = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0)
            {
                *idum += IM1;
            }
            if (j < NTAB)
            {
                iv[j] = *idum;
            }
        }
        iy = iv[0];
    }

    k     = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0)
    {
        *idum += IM1;
    }

    k     = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0)
    {
        idum2 += IM2;
    }

    j     = iy / NDIV;
    iy    = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1)
    {
        iy += IMM1;
    }

    if ((temp = AM * iy) > RNMX)
    {
        return RNMX;
    }
    else
    {
        return temp;
    }
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#ifdef NEKTAR_USE_THREAD_SAFETY
static std::mutex mutex;
#endif
template LIB_UTILITIES_EXPORT Nektar::NekDouble ran2(long *idum);
template LIB_UTILITIES_EXPORT Nektar::NekSingle ran2(long *idum);

/// \brief Fills a vector with white noise.
template <class T>
void FillWhiteNoise(int n, const T eps, T *x, const int incx, int outseed)
{
#ifdef NEKTAR_USE_THREAD_SAFETY
    // Protect the static vars here and in ran2
    std::scoped_lock l(mutex);
#endif

    // Define static variables for generating random numbers
    static int iset = 0;
    static T gset;
    static long seed = 0;

    // Bypass seed if outseed was specified
    if (outseed != 9999)
    {
        seed = long(outseed);
    }

    while (n--)
    {
        T fac, rsq, v1, v2;

        if (iset == 0)
        {
            do
            {
                v1  = 2.0 * ran2<T>(&seed) - 1.0;
                v2  = 2.0 * ran2<T>(&seed) - 1.0;
                rsq = v1 * v1 + v2 * v2;
            } while (rsq >= 1.0 || rsq == 0.0);
            fac  = sqrt(-2.0 * log(rsq) / rsq);
            gset = v1 * fac;
            iset = 1;
            *x   = eps * v2 * fac;
        }
        else
        {
            iset = 0;
            *x   = eps * gset;
        }
        x += incx;
    }
}
template LIB_UTILITIES_EXPORT void FillWhiteNoise(int n,
                                                  const Nektar::NekDouble eps,
                                                  Nektar::NekDouble *x,
                                                  const int incx, int outseed);
template LIB_UTILITIES_EXPORT void FillWhiteNoise(int n,
                                                  const Nektar::NekSingle eps,
                                                  Nektar::NekSingle *x,
                                                  const int incx, int outseed);

} // namespace Vmath
