///////////////////////////////////////////////////////////////////////////////
//
// File: LinearAdvectionDiffusionReactionKernels.hpp
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
// Description: Matrix-free operator kernels for 2D and 3D shape types
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_MF_LINEARADVECTIONDIFFUSIONREACTION_KERNELS_H
#define NEKTAR_LIBRARY_MF_LINEARADVECTIONDIFFUSIONREACTION_KERNELS_H

namespace Nektar::MatrixFree
{

using namespace tinysimd;
using vec_t = simd<NekDouble>;

// The dimension and shape kernels. NOTE: They are NOT duplicate
// templated version based on the array size like the
// operators. HOWEVER, they are forced to be INLINED. The inlining is
// critical so that when used in the templated version of the operator
// that loop unrolling occurs.

// The seven shape kernels where the work gets done.
#if defined(SHAPE_TYPE_TRI)

template <bool DEFORMED>
NEK_FORCE_INLINE static void AdvectionTriKernel(
    const size_t nq0, const size_t nq1, const vec_t *advVel_ptr,
    const vec_t *df_ptr, const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    const std::vector<vec_t, allocator<vec_t>> &deriv0,
    const std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
    constexpr auto ndf  = 4;
    constexpr auto nvel = 2;

    vec_t vx, vy;
    vec_t df0, df1, df2, df3;

    // Precompute Laplacian metricsp
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
    }

    // Step 4a: Construct Laplacian metrics
    for (size_t j = 0, cnt = 0; j < nq1; ++j)
    {
        // Collapsed coordinate mapping
        vec_t h1j = m_h1[j];
        for (size_t i = 0; i < nq0; ++i, ++cnt)
        {
            // Collapsed coordinate mapping
            vec_t h0i = m_h0[i];

            // Set deformed derivative factors
            if (DEFORMED)
            {
                df0 = df_ptr[cnt * ndf];
                df1 = df_ptr[cnt * ndf + 1];
                df2 = df_ptr[cnt * ndf + 2];
                df3 = df_ptr[cnt * ndf + 3];
            }

            // Get advection velocity
            vx = advVel_ptr[cnt * nvel];
            vy = advVel_ptr[cnt * nvel + 1];

            // Get derivatives
            vec_t d0 = deriv0[cnt];
            vec_t d1 = deriv1[cnt];

            // Use metrics to reduce operations
            vec_t metric00 = df0 + df1 * h0i;
            metric00       = metric00 * h1j;
            vec_t metric01 = df1;

            vec_t metric10 = df2 + df3 * h0i;
            metric10       = metric10 * h1j;
            vec_t metric11 = df3;
            // Use (faster) metrics to reduce operations
            // vec_t metric0 = h1j * d0;
            // vec_t metric1 = metric0 * h0i + d1;

            // Do x-advection
            vec_t tmp = metric00 * d0;
            tmp.fma(metric01, d1);
            tmp = tmp * vx;

            // Do y-advection and sum-up
            vec_t tmp2 = metric10 * d0;
            tmp2.fma(metric11, d1);
            tmp.fma(tmp2, vy);
            bwd[cnt] = bwd[cnt] + tmp; // Appends to bwd
        }
    }
}

#elif defined(SHAPE_TYPE_QUAD)

template <bool DEFORMED>
NEK_FORCE_INLINE static void AdvectionQuadKernel(
    const size_t nq0, const size_t nq1, const vec_t *advVel_ptr,
    const vec_t *df_ptr, const std::vector<vec_t, allocator<vec_t>> &deriv0,
    const std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
    constexpr auto ndf  = 4;
    constexpr auto nvel = 2;

    vec_t vx, vy;
    vec_t df0, df1, df2, df3;

    // Precompute Laplacian metricsp
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
    }

    // Apply Laplacian metrics on all quad points
    for (size_t j = 0, cnt = 0; j < nq1; ++j)
    {
        for (size_t i = 0; i < nq0; ++i, ++cnt)
        {
            // Set deformed derivative factors
            if (DEFORMED)
            {
                df0 = df_ptr[cnt * ndf];
                df1 = df_ptr[cnt * ndf + 1];
                df2 = df_ptr[cnt * ndf + 2];
                df3 = df_ptr[cnt * ndf + 3];
            }

            // Get advection velocity
            vx = advVel_ptr[cnt * nvel];
            vy = advVel_ptr[cnt * nvel + 1];

            // Get derivatives
            vec_t d0 = deriv0[cnt];
            vec_t d1 = deriv1[cnt];

            // Do x-advection
            vec_t tmp = df0 * d0;
            tmp.fma(df1, d1);
            tmp = tmp * vx;

            // Do y-advection and sum-up
            vec_t tmp2 = df2 * d0;
            tmp2.fma(df3, d1);
            tmp.fma(tmp2, vy);

            bwd[cnt] = bwd[cnt] + tmp; // Appends to bwd
        }
    }
}

#elif defined(SHAPE_TYPE_HEX)

template <bool DEFORMED>
NEK_FORCE_INLINE static void AdvectionHexKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const vec_t *advVel_ptr, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &deriv0,
    const std::vector<vec_t, allocator<vec_t>> &deriv1,
    const std::vector<vec_t, allocator<vec_t>> &deriv2,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
    constexpr auto ndf  = 9;
    constexpr auto nvel = 3;

    vec_t vx, vy, vz;
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    // Precompute deriv factors
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    // Compute advection term
    for (size_t k = 0, cnt = 0; k < nq2; k++)
    {
        for (size_t j = 0; j < nq1; j++)
        {
            for (size_t i = 0; i < nq0; i++, ++cnt)
            {
                // Set deformed derivative factors
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];
                    df4 = df_ptr[cnt * ndf + 4];
                    df5 = df_ptr[cnt * ndf + 5];
                    df6 = df_ptr[cnt * ndf + 6];
                    df7 = df_ptr[cnt * ndf + 7];
                    df8 = df_ptr[cnt * ndf + 8];
                }

                // Get advection velocity
                vx = advVel_ptr[cnt * nvel];
                vy = advVel_ptr[cnt * nvel + 1];
                vz = advVel_ptr[cnt * nvel + 2];

                // Get derivatives
                vec_t d0 = deriv0[cnt];
                vec_t d1 = deriv1[cnt];
                vec_t d2 = deriv2[cnt];

                // Do x-advection
                vec_t tmp = df0 * d0;
                tmp.fma(df1, d1);
                tmp.fma(df2, d2);
                tmp = tmp * vx;

                // Do y-advection
                vec_t tmp2 = df3 * d0;
                tmp2.fma(df4, d1);
                tmp2.fma(df5, d2);
                tmp.fma(tmp2, vy);

                // Do z-advection and sum-up
                tmp2 = df6 * d0;
                tmp2.fma(df7, d1);
                tmp2.fma(df8, d2);
                tmp.fma(tmp2, vz);

                bwd[cnt] = bwd[cnt] + tmp; // Appends to bwd
            }
        }
    }
}

#elif defined(SHAPE_TYPE_TET)

template <bool DEFORMED>
NEK_FORCE_INLINE static void AdvectionTetKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const vec_t *advVel_ptr, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    const std::vector<vec_t, allocator<vec_t>> &m_h2,
    const std::vector<vec_t, allocator<vec_t>> &m_h3,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
    constexpr auto ndf  = 9;
    constexpr auto nvel = 3;

    vec_t vx, vy, vz;
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    // Precompute deriv factors
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    // Compute advection term
    for (size_t k = 0, cnt = 0; k < nq2; ++k)
    {
        vec_t h3 = m_h3[k];
        for (size_t j = 0; j < nq1; ++j)
        {
            vec_t h1   = m_h1[j];
            vec_t h2   = m_h2[j];
            vec_t h2h3 = h2 * h3;
            vec_t h1h3 = h1 * h3;

            for (int i = 0; i < nq0; ++i, ++cnt)
            {
                vec_t h0 = m_h0[i];
                // Set deformed derivative factors
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];
                    df4 = df_ptr[cnt * ndf + 4];
                    df5 = df_ptr[cnt * ndf + 5];
                    df6 = df_ptr[cnt * ndf + 6];
                    df7 = df_ptr[cnt * ndf + 7];
                    df8 = df_ptr[cnt * ndf + 8];
                }

                // Get advection velocity
                vx = advVel_ptr[cnt * nvel];
                vy = advVel_ptr[cnt * nvel + 1];
                vz = advVel_ptr[cnt * nvel + 2];

                // Get derivatives
                vec_t d0 = deriv0[cnt];
                vec_t d1 = deriv1[cnt];
                vec_t d2 = deriv2[cnt];

                // Use metrics to reduce operations
                vec_t metric0 = h2h3 * d0;
                vec_t metric1 = metric0 * h0;
                vec_t metric2 = metric1 + d2;
                metric2.fma(h1h3, d1);
                metric1.fma(h3, d1);

                // Do x-advection
                vec_t tmp = df0 * metric0;
                tmp.fma(df1, metric1);
                tmp.fma(df2, metric2);
                tmp = tmp * vx;

                // Do y-advection
                vec_t tmp2 = df3 * metric0;
                tmp2.fma(df4, metric1);
                tmp2.fma(df5, metric2);
                tmp.fma(tmp2, vy);

                // Do z-advection and sum-up
                tmp2 = df6 * metric0;
                tmp2.fma(df7, metric1);
                tmp2.fma(df8, metric2);
                tmp.fma(tmp2, vz);

                bwd[cnt] = bwd[cnt] + tmp; // Appends to bwd
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PRISM)

template <bool DEFORMED>
NEK_FORCE_INLINE static void AdvectionPrismKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const vec_t *advVel_ptr, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
    constexpr auto ndf  = 9;
    constexpr auto nvel = 3;

    vec_t vx, vy, vz;
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    // Precompute deriv factors
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    // Compute advection term
    for (size_t k = 0, cnt = 0; k < nq2; ++k)
    {
        vec_t h1 = m_h1[k];
        for (size_t j = 0; j < nq1; ++j)
        {
            for (size_t i = 0; i < nq0; ++i, ++cnt)
            {
                vec_t h0 = m_h0[i];
                // Set deformed derivative factors
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];
                    df4 = df_ptr[cnt * ndf + 4];
                    df5 = df_ptr[cnt * ndf + 5];
                    df6 = df_ptr[cnt * ndf + 6];
                    df7 = df_ptr[cnt * ndf + 7];
                    df8 = df_ptr[cnt * ndf + 8];
                }

                // Get advection velocity
                vx = advVel_ptr[cnt * nvel];
                vy = advVel_ptr[cnt * nvel + 1];
                vz = advVel_ptr[cnt * nvel + 2];

                // Get derivatives
                vec_t d0 = deriv0[cnt];
                vec_t d1 = deriv1[cnt];
                vec_t d2 = deriv2[cnt];

                // Use metrics to reduce operations
                vec_t metric0 = h1 * d0;
                vec_t metric2 = h0 * metric0 + d2;

                // Do x-advection
                vec_t tmp = df0 * metric0;
                tmp.fma(df1, d1);
                tmp.fma(df2, metric2);
                tmp = tmp * vx;

                // Do y-advection
                vec_t tmp2 = df3 * metric0;
                tmp2.fma(df4, d1);
                tmp2.fma(df5, metric2);
                tmp.fma(tmp2, vy);

                // Do z-advection and sum-up
                tmp2 = df6 * metric0;
                tmp2.fma(df7, d1);
                tmp2.fma(df8, metric2);
                tmp.fma(tmp2, vz);

                bwd[cnt] = bwd[cnt] + tmp; // Appends to bwd
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PYR)

template <bool DEFORMED>
NEK_FORCE_INLINE static void AdvectionPyrKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const vec_t *advVel_ptr, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    const std::vector<vec_t, allocator<vec_t>> &m_h2,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
    constexpr auto ndf  = 9;
    constexpr auto nvel = 3;

    vec_t vx, vy, vz;
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    // Precompute deriv factors
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    // Compute advection term
    for (size_t k = 0, cnt = 0; k < nq2; ++k)
    {
        vec_t h2 = m_h2[k];
        for (size_t j = 0; j < nq1; ++j)
        {
            vec_t h1   = m_h1[j];
            vec_t h1h2 = h1 * h2;
            for (size_t i = 0; i < nq0; ++i, cnt++)
            {
                vec_t h0   = m_h0[i];
                vec_t h0h2 = h0 * h2;

                // Set deformed derivative factors
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];
                    df4 = df_ptr[cnt * ndf + 4];
                    df5 = df_ptr[cnt * ndf + 5];
                    df6 = df_ptr[cnt * ndf + 6];
                    df7 = df_ptr[cnt * ndf + 7];
                    df8 = df_ptr[cnt * ndf + 8];
                }

                // Get advection velocity
                vx = advVel_ptr[cnt * nvel];
                vy = advVel_ptr[cnt * nvel + 1];
                vz = advVel_ptr[cnt * nvel + 2];

                // Get derivatives
                vec_t d0 = deriv0[cnt];
                vec_t d1 = deriv1[cnt];
                vec_t d2 = deriv2[cnt];

                // Use metrics to reduce operations
                vec_t metric0 = h2 * d0;
                vec_t metric1 = h2 * d1;
                vec_t metric2 = d2;
                metric2.fma(h0h2, d0);
                metric2.fma(h1h2, d1);

                // Do x-advection
                vec_t tmp = df0 * metric0;
                tmp.fma(df1, metric1);
                tmp.fma(df2, metric2);
                tmp = tmp * vx;

                // Do y-advection
                vec_t tmp2 = df3 * metric0;
                tmp2.fma(df4, metric1);
                tmp2.fma(df5, metric2);
                tmp.fma(tmp2, vy);

                // Do z-advection and sum-up
                tmp2 = df6 * metric0;
                tmp2.fma(df7, metric1);
                tmp2.fma(df8, metric2);
                tmp.fma(tmp2, vz);

                bwd[cnt] = bwd[cnt] + tmp; // Appends to bwd
            }
        }
    }
}

#endif // SHAPE_TYPE

#if defined(SHAPE_TYPE_TRI)
template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void GetLinearAdvectionDiffusionReaction2DHalfSpace(
    const size_t nq0, const size_t nq1, const Array<OneD, const NekDouble> &z0,
    const Array<OneD, const NekDouble> &z1,
    std::vector<vec_t, allocator<vec_t>> &h0,
    std::vector<vec_t, allocator<vec_t>> &h1)
{
    h0.resize(nq0);
    h1.resize(nq1);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int j = 0; j < nq1; ++j)
    {
        h1[j] = 2.0 / (1 - z1[j]);
    }
}

#elif defined(SHAPE_TYPE_TET) || defined(SHAPE_TYPE_PRISM) ||                  \
    defined(SHAPE_TYPE_PYR)
template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void GetLinearAdvectionDiffusionReaction3DHalfSpace(
    const size_t nq0, [[maybe_unused]] const size_t nq1, const size_t nq2,
    const Array<OneD, const NekDouble> &z0,
    [[maybe_unused]] const Array<OneD, const NekDouble> &z1,
    const Array<OneD, const NekDouble> &z2,
    std::vector<vec_t, allocator<vec_t>> &h0,
    std::vector<vec_t, allocator<vec_t>> &h1,
    [[maybe_unused]] std::vector<vec_t, allocator<vec_t>> &h2,
    [[maybe_unused]] std::vector<vec_t, allocator<vec_t>> &h3)
{
#if defined(SHAPE_TYPE_TET)
    h0.resize(nq0);
    h1.resize(nq1);
    h2.resize(nq1);
    h3.resize(nq2);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int j = 0; j < nq1; ++j)
    {
        h1[j] = 0.5 * (1 + z1[j]);
        h2[j] = 2.0 / (1 - z1[j]);
    }

    for (int k = 0; k < nq2; ++k)
    {
        h3[k] = 2.0 / (1 - z2[k]);
    }

#elif defined(SHAPE_TYPE_PRISM)

    h0.resize(nq0);
    h1.resize(nq2);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int k = 0; k < nq2; ++k)
    {
        h1[k] = 2.0 / (1 - z2[k]);
    }

#elif defined(SHAPE_TYPE_PYR)

    h0.resize(nq0);
    h1.resize(nq1);
    h2.resize(nq2);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int j = 0; j < nq1; ++j)
    {
        h1[j] = 0.5 * (1 + z1[j]);
    }

    for (int k = 0; k < nq2; ++k)
    {
        h2[k] = 2.0 / (1 - z2[k]);
    }
#endif
}
#endif

// The dimension kernels which select the shape kernel.
#if defined(SHAPE_DIMENSION_2D)

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED>
NEK_FORCE_INLINE static void Advection2DKernel(
    const size_t nq0, const size_t nq1, const vec_t *advVel_ptr,
    const vec_t *df_ptr,
    [[maybe_unused]] const std::vector<vec_t, allocator<vec_t>> &h0,
    [[maybe_unused]] const std::vector<vec_t, allocator<vec_t>> &h1,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
#if defined(SHAPE_TYPE_TRI)
    AdvectionTriKernel<DEFORMED>(nq0, nq1, advVel_ptr, df_ptr, h0, h1, deriv0,
                                 deriv1, bwd);
#elif defined(SHAPE_TYPE_QUAD)
    AdvectionQuadKernel<DEFORMED>(nq0, nq1, advVel_ptr, df_ptr, deriv0, deriv1,
                                  bwd);
#endif
}

#elif defined(SHAPE_DIMENSION_3D)

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED>
NEK_FORCE_INLINE static void Advection3DKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const vec_t *advVel_ptr, const vec_t *df_ptr,
    [[maybe_unused]] const std::vector<vec_t, allocator<vec_t>> &h0,
    [[maybe_unused]] const std::vector<vec_t, allocator<vec_t>> &h1,
    [[maybe_unused]] const std::vector<vec_t, allocator<vec_t>> &h2,
    [[maybe_unused]] const std::vector<vec_t, allocator<vec_t>> &h3,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2,
    std::vector<vec_t, allocator<vec_t>> &bwd)
{
#if defined(SHAPE_TYPE_HEX)
    AdvectionHexKernel<DEFORMED>(nq0, nq1, nq2, advVel_ptr, df_ptr, deriv0,
                                 deriv1, deriv2, bwd);
#elif defined(SHAPE_TYPE_TET)
    AdvectionTetKernel<DEFORMED>(nq0, nq1, nq2, advVel_ptr, df_ptr, h0, h1, h2,
                                 h3, deriv0, deriv1, deriv2, bwd);
#elif defined(SHAPE_TYPE_PRISM)
    AdvectionPrismKernel<DEFORMED>(nq0, nq1, nq2, advVel_ptr, df_ptr, h0, h1,
                                   deriv0, deriv1, deriv2, bwd);
#elif defined(SHAPE_TYPE_PYR)
    AdvectionPyrKernel<DEFORMED>(nq0, nq1, nq2, advVel_ptr, df_ptr, h0, h1, h2,
                                 deriv0, deriv1, deriv2, bwd);
#endif
}

#endif // SHAPE_DIMENSION

} // namespace Nektar::MatrixFree

#endif
