///////////////////////////////////////////////////////////////////////////////
//
// File: PhysInterp1DScaledKernels.hpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_MF_PHYSINTERP1DSCALED_KERNELS_HPP
#define NEKTAR_LIBRARY_MF_PHYSINTERP1DSCALED_KERNELS_HPP

#include <LibUtilities/BasicUtils/NekInline.hpp>

namespace Nektar::MatrixFree
{

using namespace tinysimd;
using vec_t = simd<NekDouble>;

// Workspace - used to dynamically get the workspace size needed for
// temporary memory.
#if defined(SHAPE_DIMENSION_1D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysInterp1DScaled1DWorkspace(
    [[maybe_unused]] const size_t nq_in0, [[maybe_unused]] const size_t nq_out0)

{
}

#elif defined(SHAPE_DIMENSION_2D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysInterp1DScaled2DWorkspace(
    [[maybe_unused]] const size_t nq_in0, [[maybe_unused]] const size_t nq_in1,
    [[maybe_unused]] const size_t nq_out0,
    [[maybe_unused]] const size_t nq_out1, size_t &wsp0Size)
{
    wsp0Size = std::max(wsp0Size, nq_in0 * nq_out0);
}

#elif defined(SHAPE_DIMENSION_3D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysInterp1DScaled3DWorkspace(
    [[maybe_unused]] const size_t nq_in0, [[maybe_unused]] const size_t nq_in1,
    [[maybe_unused]] const size_t nq_in2, [[maybe_unused]] const size_t nq_out0,
    [[maybe_unused]] const size_t nq_out1,
    [[maybe_unused]] const size_t nq_out2, size_t &wsp0Size, size_t &wsp1Size)
{

    wsp0Size =
        std::max(wsp0Size, nq_out0 * nq_in1 * nq_in2); // nq_in1 == nq_in2
    wsp1Size =
        std::max(wsp1Size, nq_out0 * nq_out1 * nq_in2); // nq_out0 == nq_out1
}

#endif // SHAPE_DIMENSION

// The dimension kernels which select the shape kernel.
#if defined(SHAPE_DIMENSION_1D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysInterp1DScaled1DKernel(
    const size_t nq_in0, const size_t nq_out0,
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &basis0,
    std::vector<vec_t, allocator<vec_t>> &out)
{
    for (int i = 0; i < nq_out0; ++i)
    {
        vec_t tmp = in[0] * basis0[i]; // Load 2x

        for (int p = 1; p < nq_in0; ++p)
        {
            tmp.fma(in[p], basis0[p * nq_out0 + i]); // Load 2x
        }

        out[i] = tmp; // Store 1x
    }
}

#elif defined(SHAPE_DIMENSION_2D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysInterp1DScaled2DKernel(
    const size_t nq_in0, const size_t nq_in1, const size_t nq_out0,
    const size_t nq_out1, const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &basis0,
    const std::vector<vec_t, allocator<vec_t>> &basis1,
    std::vector<vec_t, allocator<vec_t>> &wsp,
    std::vector<vec_t, allocator<vec_t>> &out)
{
    for (int i = 0, cnt_iq = 0; i < nq_out0; ++i)
    {
        for (int q = 0, cnt_pq = 0; q < nq_in1; ++q, ++cnt_iq)
        {
            vec_t tmp = in[cnt_pq] * basis0[i]; // Load 2x
            ++cnt_pq;
            for (int p = 1; p < nq_in0; ++p, ++cnt_pq)
            {
                tmp.fma(in[cnt_pq], basis0[p * nq_out0 + i]); // Load 2x
            }
            wsp[cnt_iq] = tmp; // Store 1x
        }
    }

    for (int j = 0, cnt_ij = 0; j < nq_out1; ++j)
    {
        for (int i = 0, cnt_iq = 0; i < nq_out0; ++i, ++cnt_ij)
        {
            vec_t tmp = wsp[cnt_iq] * basis1[j]; // Load 2x
            ++cnt_iq;
            for (int q = 1; q < nq_in1; ++q, ++cnt_iq)
            {
                tmp.fma(wsp[cnt_iq], basis1[q * nq_out1 + j]); // Load 2x
            }
            out[cnt_ij] = tmp; // Store 1x
        }
    }
}

#elif defined(SHAPE_DIMENSION_3D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysInterp1DScaled3DKernel(
    const size_t nq_in0, const size_t nq_in1, const size_t nq_in2,
    const size_t nq_out0, const size_t nq_out1, const size_t nq_out2,
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &basis0,
    const std::vector<vec_t, allocator<vec_t>> &basis1,
    const std::vector<vec_t, allocator<vec_t>> &basis2,
    std::vector<vec_t, allocator<vec_t>> &sum_irq, // nq_out0 * nq_in2 * nq_in1
    std::vector<vec_t, allocator<vec_t>> &sum_jir, // nq_out1 * nq_out0 * nq_in2
    std::vector<vec_t, allocator<vec_t>> &out)
{
    for (int i = 0, cnt_irq = 0; i < nq_out0; ++i)
    {
        for (int r = 0, cnt_rqp = 0; r < nq_in2; ++r)
        {
            for (int q = 0; q < nq_in1; ++q, ++cnt_irq)
            {
                vec_t tmp = in[cnt_rqp] * basis0[i];
                ++cnt_rqp;

                for (int p = 1; p < nq_in0; ++p, ++cnt_rqp)
                {
                    tmp.fma(in[cnt_rqp], basis0[p * nq_out0 + i]);
                }

                sum_irq[cnt_irq] = tmp;
            }
        }
    }

    for (int j = 0, cnt_jir = 0; j < nq_out1; ++j)
    {
        for (int i = 0, cnt_irq = 0; i < nq_out0; ++i)
        {
            for (int r = 0; r < nq_in2; ++r, ++cnt_jir)
            {
                vec_t tmp = sum_irq[cnt_irq] * basis1[j];
                ++cnt_irq;

                for (int q = 1; q < nq_in1; ++q)
                {
                    tmp.fma(sum_irq[cnt_irq++], basis1[q * nq_out1 + j]);
                }

                sum_jir[cnt_jir] = tmp;
            }
        }
    }

    for (int k = 0, cnt_kji = 0; k < nq_out2; ++k)
    {
        for (int j = 0, cnt_jir = 0; j < nq_out1; ++j)
        {
            for (int i = 0; i < nq_out0; ++i, ++cnt_kji)
            {
                vec_t tmp = sum_jir[cnt_jir] * basis2[k];
                ++cnt_jir;

                for (int r = 1; r < nq_in2; ++r)
                {
                    tmp.fma(sum_jir[cnt_jir++], basis2[r * nq_out2 + k]);
                }

                out[cnt_kji] = tmp;
            }
        }
    }
}

#endif // SHAPE_DIMENSION

} // namespace Nektar::MatrixFree

#endif
