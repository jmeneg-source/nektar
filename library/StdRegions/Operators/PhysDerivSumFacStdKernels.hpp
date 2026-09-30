///////////////////////////////////////////////////////////////////////////////
//
// File: PhysDerivSumFacStdKernels.hpp
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
// Description: Inner kernel functions leveraging the  AVX impelmentaiton of
// PHysDeriv. They are then used in a scalar manner within  StdRegions
//
///////////////////////////////////////////////////////////////////////////////

template <typename simd_type>
NEK_FORCE_INLINE static void PhysDerivTensor1DKernel(const unsigned int nq0,
                                                     const simd_type *in,
                                                     const simd_type *D0,
                                                     simd_type *out_d0)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // D0 * in
    for (unsigned int i = 0; i < nq0; ++i)
    { // Row index of D0 matrix

        simd_type prod_sum = 0.0;
        for (unsigned int k = 0; k < nq0; ++k)
        {                                    // Col index of D0, row index of IN
            simd_type v1 = D0[k * nq0 + i];  // Load 1x
            simd_type v2 = simd_type(in[k]); // Load 1x

            prod_sum.fma(v1, v2);
        }

        out_d0[i] = prod_sum; // Store 1x
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void PhysDerivTensor2DKernel(
    const unsigned int nq0, const unsigned int nq1, const simd_type *in,
    const simd_type *D0, const simd_type *D1, simd_type *out_d0,
    simd_type *out_d1, bool Deriv0 = true, bool Deriv1 = true)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // D0 * in
    if (Deriv0) // backwards compatibility with StdRegsion PhyTensorDeriv
    {
        for (unsigned int i = 0; i < nq0; ++i)
        { // Row index of D0 matrix
            for (unsigned int j = 0; j < nq1; ++j)
            { // Col index of IN matrix

                simd_type prod_sum = 0.0;
                for (unsigned int k = 0; k < nq0; ++k)
                { // Col index of D0, row index of IN
                    simd_type v1 = D0[k * nq0 + i]; // Load 1x
                    simd_type v2 = in[j * nq0 + k]; // Load 1x

                    prod_sum.fma(v1, v2);
                }

                out_d0[j * nq0 + i] = prod_sum; // Store 1x
            }
        }
    }

    // D1 * in
    if (Deriv1)
    {
        // in * D1^T
        for (unsigned int i = 0; i < nq0; ++i)
        { // row index for grid
            for (unsigned int j = 0; j < nq1; ++j)
            { // Column index for D1^T (row idx for D1)

                simd_type prod_sum = 0.0;
                for (unsigned int k = 0; k < nq1; ++k)
                {
                    simd_type v1 = in[k * nq0 + i]; // Load 1x
                    simd_type v2 = D1[k * nq1 + j]; // Load 1x

                    prod_sum.fma(v1, v2);
                }

                out_d1[j * nq0 + i] = prod_sum; // Store 1x
            }
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void SumDerivTensor2DKernel(
    const unsigned int nq0, const unsigned int nq1, const simd_type *in0,
    const simd_type *in1, const simd_type *D0, const simd_type *D1,
    simd_type *out, const typename simd_type::scalarType scale = 0.0,
    bool Deriv0 = true, bool Deriv1 = true)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // D0 * in
    if (Deriv0) // backwards compatibility with StdRegsion PhyTensorDeriv
    {
        for (unsigned int i = 0; i < nq0; ++i)
        { // Row index of D0 matrix
            for (unsigned int j = 0; j < nq1; ++j)
            { // Col index of IN matrix

                simd_type prod_sum = 0.0;
                for (unsigned int k = 0; k < nq0; ++k)
                { // Col index of D0, row index of IN
                    simd_type v1 = D0[i * nq0 + k];  // Load 1x
                    simd_type v2 = in0[j * nq0 + k]; // Load 1x

                    prod_sum.fma(v1, v2);
                }

                out[j * nq0 + i] *= simd_type(scale);
                out[j * nq0 + i] += prod_sum; // Store 1x
            }
        }
    }

    // D1 * in
    if (Deriv1)
    {
        // in * D1^T
        for (unsigned int j = 0; j < nq1; ++j)
        { // Column index for D1^T (row idx for D1)
            for (unsigned int i = 0; i < nq0; ++i)
            { // row index for grid

                simd_type prod_sum = 0.0;
                for (unsigned int k = 0; k < nq1; ++k)
                {
                    simd_type v1 = in1[k * nq0 + i]; // Load 1x
                    simd_type v2 = D1[j * nq1 + k];  // Load 1x

                    prod_sum.fma(v1, v2);
                }

                out[j * nq0 + i] += prod_sum; // Store 1x
            }
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void PhysDerivTensor3DKernel(
    const unsigned int nq0, const unsigned int nq1, const unsigned int nq2,
    const simd_type *in, const simd_type *D0, const simd_type *D1,
    const simd_type *D2, simd_type *out_d0, simd_type *out_d1,
    simd_type *out_d2, bool Deriv0 = true, bool Deriv1 = true,
    bool Deriv2 = true)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // Direction 0
    if (Deriv0)
    {
        for (unsigned int i = 0; i < nq0; ++i)
        {
            for (unsigned int j = 0; j < nq1 * nq2; ++j)
            {
                simd_type prod_sum = 0.0;
                for (unsigned int k = 0; k < nq0; ++k)
                {
                    simd_type v1 = D0[k * nq0 + i]; // Load 1x
                    simd_type v2 = in[j * nq0 + k]; // Load 1x

                    prod_sum.fma(v1, v2);
                }

                out_d0[j * nq0 + i] = prod_sum; // Store 1x
            }
        }
    }

    // Direction 1
    if (Deriv1)
    {
        for (unsigned int block = 0; block < nq2; ++block)
        {
            unsigned int start = block * nq0 * nq1;

            for (unsigned int i = 0; i < nq0; ++i)
            {
                for (unsigned int j = 0; j < nq1; ++j)
                {
                    simd_type prod_sum = 0.0;
                    for (unsigned int k = 0; k < nq1; ++k)
                    {
                        simd_type v1 = in[start + k * nq0 + i]; // Load 1x
                        simd_type v2 = D1[k * nq1 + j];         // Load 1x

                        prod_sum.fma(v1, v2);
                    }

                    out_d1[start + j * nq0 + i] = prod_sum; // Store 1x
                }
            }
        }
    }

    // Direction 2
    if (Deriv2)
    {
        for (unsigned int i = 0; i < nq0 * nq1; ++i)
        {
            for (unsigned int j = 0; j < nq2; ++j)
            {
                simd_type prod_sum = 0.0;
                for (unsigned int k = 0; k < nq2; ++k)
                {
                    simd_type v1 = simd_type(in[k * nq0 * nq1 + i]); // Load 1x
                    simd_type v2 = D2[k * nq2 + j];                  // Load 1x

                    prod_sum.fma(v1, v2);
                }

                out_d2[j * nq0 * nq1 + i] = prod_sum; // Store 1x
            }
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void SumDerivTensor3DKernel(
    const unsigned int nq0, const unsigned int nq1, const unsigned int nq2,
    const simd_type *in0, const simd_type *in1, const simd_type *in2,
    const simd_type *D0, const simd_type *D1, const simd_type *D2,
    simd_type *out, const typename simd_type::scalarType scale = 0.0,
    bool Deriv0 = true, bool Deriv1 = true, bool Deriv2 = true)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // Compared with PhysDerivTensor3DKernel, here we must multiply by the
    // transpose of D matrix

    // Direction 0
    if (Deriv0)
    {
        for (unsigned int p = 0; p < nq0; ++p)
        {
            unsigned int cnt_kji = 0, cnt_kj = 0;
            for (unsigned int k = 0; k < nq2; ++k)
            {
                for (unsigned int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    simd_type prod_sum = 0.0;
                    for (unsigned int i = 0; i < nq0; ++i, ++cnt_kji)
                    {
                        simd_type v1 = D0[p * nq0 + i]; // Load 1x
                        simd_type v2 = in0[cnt_kji];    // Load 1x

                        prod_sum.fma(v1, v2);
                    }
                    out[cnt_kj * nq0 + p] *= simd_type(scale);
                    out[cnt_kj * nq0 + p] += prod_sum; // Store 1x
                }
            }
        }
    }

    // Direction 1
    if (Deriv1)
    {
        for (unsigned int block = 0; block < nq2; ++block)
        {
            unsigned int start = block * nq0 * nq1;

            for (unsigned int i = 0; i < nq0; ++i)
            {
                for (unsigned int j = 0; j < nq1; ++j)
                {
                    simd_type prod_sum = 0.0;
                    for (unsigned int k = 0; k < nq1; ++k)
                    {
                        simd_type v1 = in1[start + k * nq0 + i]; // Load 1x
                        simd_type v2 = D1[j * nq1 + k];          // Load 1x

                        prod_sum.fma(v1, v2);
                    }

                    out[start + j * nq0 + i] += prod_sum; // Store 1x
                }
            }
        }
    }

    // Direction 2
    if (Deriv2)
    {
        unsigned int cnt_hi = 0;
        for (unsigned int h = 0; h < nq1; ++h)
        {
            for (unsigned int i = 0; i < nq0; ++i, ++cnt_hi)
            {
                for (unsigned int j = 0; j < nq2; ++j)
                {
                    simd_type prod_sum = 0.0;
                    for (unsigned int k = 0; k < nq2; ++k)
                    {
                        simd_type v1 = in2[k * nq0 * nq1 + cnt_hi]; // Load 1x
                        simd_type v2 = D2[j * nq2 + k];             // Load 1x

                        prod_sum.fma(v1, v2);
                    }

                    out[j * nq0 * nq1 + cnt_hi] += prod_sum; // Store 1x
                }
            }
        }
    }
}
