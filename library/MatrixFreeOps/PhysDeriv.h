///////////////////////////////////////////////////////////////////////////////
//
// File: PhysDeriv.h
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

#ifndef MF_PHYSDERIV_H
#define MF_PHYSDERIV_H

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"
#include "PhysDerivKernels.hpp"

namespace Nektar::MatrixFree
{

// As each opertor has seven shapes over three dimension to get to the
// "work" each operator uses a series of preprocessor directives based
// on the dimension and shape type so to limit the code
// inclusion. This keeps the library size as minimal as possible while
// removing runtime conditionals.

// The preprocessor directives, SHAPE_DIMENSION_?D and SHAPE_TYPE_*
// are constructed by CMake in the CMakeLists.txt file which uses an
// implementation file.  See the CMakeLists.txt files for more
// details.
template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED = false>
struct PhysDerivTemplate : public PhysDeriv, public Helper<SHAPE_TYPE, DEFORMED>
{
    PhysDerivTemplate(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
        : PhysDeriv(basis, nElmt), Helper<SHAPE_TYPE, DEFORMED>(basis, nElmt)
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<PhysDerivTemplate<SHAPE_TYPE, DEFORMED>>(basis,
                                                                         nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output) final
    {
        const auto nq0 = this->m_nq[0];
#if defined(SHAPE_DIMENSION_2D)
        const auto nq1 = this->m_nq[1];
#elif defined(SHAPE_DIMENSION_3D)
        const auto nq1 = this->m_nq[1];
        const auto nq2 = this->m_nq[2];
#endif
#include "SwitchPoints.h"
    }

    // There must be separate 1D, 2D, and 3D operators because the
    // helper base class is based on the shape type dim. Thus indices
    // for the basis must be respected.

    // Further there are duplicate 1D, 2D, and 3D operators so to
    // allow for compile time array sizes to be part of the template
    // and thus gain loop unrolling. The only difference is the
    // location of the size value, in the template or the function:
    // foo<int bar>() {} vs foo() { int bar = ... }

#if defined(SHAPE_DIMENSION_1D)

    // Non-size based operator.
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 0, "PhysDerivTemplate::Operator1D: Cannot "
                                    "call 1D routine with no output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator1D: Operator "
                                     "not set up for other dimensions.")

        const auto nq0 = this->m_nq[0];

        const auto nqTot    = nq0;
        const auto nqBlocks = nqTot * vec_t::width;

        const auto ndf        = output.size();
        constexpr int max_ndf = 3;

        auto *inptr = &input[0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Call 1D kernel
        const vec_t *df_ptr   = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf] = {}; // max_ndf is a constexpr

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_ndf] = {}; // max_ndf is a constexpr

        for (int d = 0; d < ndf; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            // Get the basic derivative
            PhysDerivTensor1DKernel(nq0, tmpIn, this->m_D[0], tmpOut[0]);

            // Calculate physical derivative
            PhysDeriv1DKernel<DEFORMED>(nq0, ndf, df_ptr, df_tmp, tmpOut);

            // De-interleave and store data
            for (int d = 0; d < ndf; ++d)
            {
                // de-interleave and store data
                // deinterleave_store(tmpOut[d], nqTot, locField);
                // std::copy(locField, locField + nqBlocks, out_ptr[d]);
                deinterleave_unalign_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }

            inptr += nqBlocks;
            df_ptr += dfSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            // Get the basic derivative
            PhysDerivTensor1DKernel(nq0, tmpIn, this->m_D[0], tmpOut[0]);

            // Calculate physical derivative
            PhysDeriv1DKernel<DEFORMED>(nq0, ndf, df_ptr, df_tmp, tmpOut);

            // De-interleave and store data
            for (int d = 0; d < ndf; ++d)
            {
                // de-interleave and store data
                deinterleave_store(tmpOut[d], nqTot, locField);
                std::copy(locField, locField + acturalSize, out_ptr[d]);
            }
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nq0>
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 0, "PhysDerivTemplate::Operator1D: Cannot "
                                    "call 1D routine with no output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator1D: Operator "
                                     "not set up for other dimensions.")

        constexpr auto nqTot    = nq0;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        const auto ndf        = output.size();
        constexpr int max_ndf = 3;

        auto *inptr = &input[0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        PhysDeriv1DWorkspace<SHAPE_TYPE>(nq0);

        // Call 1D kernel
        const vec_t *df_ptr   = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf] = {}; // max_ndf is a constexpr

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_ndf] = {}; // max_ndf is a constexpr

        for (int d = 0; d < ndf; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqBlocks];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            // Get the basic derivative
            PhysDerivTensor1DKernel(nq0, tmpIn, this->m_D[0], tmpOut[0]);

            // Calculate physical derivative
            PhysDeriv1DKernel<DEFORMED>(nq0, ndf, df_ptr, df_tmp, tmpOut);

            // De-interleave and store data
            for (int d = 0; d < ndf; ++d)
            {
                // de-interleave and store data
                // deinterleave_store(tmpOut[d], nqTot, locField);
                // std::copy(locField, locField + nqBlocks, out_ptr[d]);
                deinterleave_unalign_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }

            inptr += nqBlocks;
            df_ptr += dfSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            // Get the basic derivative
            PhysDerivTensor1DKernel(nq0, tmpIn, this->m_D[0], tmpOut[0]);

            // Calculate physical derivative
            PhysDeriv1DKernel<DEFORMED>(nq0, ndf, df_ptr, df_tmp, tmpOut);

            // De-interleave and store data
            for (int d = 0; d < ndf; ++d)
            {
                // de-interleave and store data
                deinterleave_store(tmpOut[d], nqTot, locField);
                std::copy(locField, locField + acturalSize, out_ptr[d]);
            }
        }
    }

#elif defined(SHAPE_DIMENSION_2D)

    // Non-size based operator.
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 1, "PhysDerivTemplate::Operator2D: Cannot "
                                    "call 2D routine with one output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator2D: Operator "
                                     "not set up for 3D coordinates.");

        const auto nq0 = this->m_nq[0];
        const auto nq1 = this->m_nq[1];

        const auto nqTot    = nq0 * nq1;
        const auto nqBlocks = nqTot * vec_t::width;

        const auto outdim        = output.size();
        constexpr int max_outdim = 3;

        const auto ndf        = 2 * outdim;
        constexpr int max_ndf = 2 * max_outdim;

        auto *inptr = &input[0];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_outdim] = {}; // max_outdim is a constexpr

        for (int d = 0; d < outdim; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        const vec_t *df_ptr   = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf] = {}; // max_ndf is a constexpr

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            // Results written to out_d0, out_d1
            PhysDerivTensor2DKernel(nq0, nq1, tmpIn, this->m_D[0], this->m_D[1],
                                    tmpOut[0], tmpOut[1]);
            // Calculate physical derivative
            PhysDeriv2DKernel<DEFORMED>(nq0, nq1, outdim, this->m_Z[0],
                                        this->m_Z[1], df_ptr, df_tmp, tmpOut);

            inptr += nqBlocks;
            df_ptr += dfSize;
            // de-interleave and store data
            for (int d = 0; d < outdim; ++d)
            {
                // de-interleave and store data
                // deinterleave_store(tmpOut[d], nqTot, locField);
                // std::copy(locField, locField + nqBlocks, out_ptr[d]);
                deinterleave_unalign_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            // Results written to out_d0, out_d1
            PhysDerivTensor2DKernel(nq0, nq1, tmpIn, this->m_D[0], this->m_D[1],
                                    tmpOut[0], tmpOut[1]);
            // Calculate physical derivative
            PhysDeriv2DKernel<DEFORMED>(nq0, nq1, outdim, this->m_Z[0],
                                        this->m_Z[1], df_ptr, df_tmp, tmpOut);

            // De-interleave and store data
            for (int d = 0; d < outdim; ++d)
            {
                // de-interleave and store data
                deinterleave_store(tmpOut[d], nqTot, locField);
                std::copy(locField, locField + acturalSize, out_ptr[d]);
            }
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nq0, int nq1>
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 1, "PhysDerivTemplate::Operator2D: Cannot "
                                    "call 2D routine with one output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator2D: Operator "
                                     "not set up for 3D coordinates.");

        constexpr auto nqTot    = nq0 * nq1;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        const auto outdim        = output.size();
        constexpr int max_outdim = 3;

        const auto ndf        = 2 * outdim;
        constexpr int max_ndf = 2 * max_outdim;

        auto *inptr = &input[0];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_outdim] = {}; // max_outdim is a constexpr

        for (int d = 0; d < outdim; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        const vec_t *df_ptr   = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf] = {}; // max_ndf is a constexpr

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqBlocks];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            // Results written to out_d0, out_d1
            PhysDerivTensor2DKernel(nq0, nq1, tmpIn, this->m_D[0], this->m_D[1],
                                    tmpOut[0], tmpOut[1]);
            // Calculate physical derivative
            PhysDeriv2DKernel<DEFORMED>(nq0, nq1, outdim, this->m_Z[0],
                                        this->m_Z[1], df_ptr, df_tmp, tmpOut);

            inptr += nqBlocks;
            df_ptr += dfSize;
            // de-interleave and store data
            for (int d = 0; d < outdim; ++d)
            {
                // de-interleave and store data
                // deinterleave_store(tmpOut[d], nqTot, locField);
                // std::copy(locField, locField + nqBlocks, out_ptr[d]);
                deinterleave_unalign_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            // Results written to out_d0, out_d1
            PhysDerivTensor2DKernel(nq0, nq1, tmpIn, this->m_D[0], this->m_D[1],
                                    tmpOut[0], tmpOut[1]);
            // Calculate physical derivative
            PhysDeriv2DKernel<DEFORMED>(nq0, nq1, outdim, this->m_Z[0],
                                        this->m_Z[1], df_ptr, df_tmp, tmpOut);

            // De-interleave and store data
            for (int d = 0; d < outdim; ++d)
            {
                // de-interleave and store data
                deinterleave_store(tmpOut[d], nqTot, locField);
                std::copy(locField, locField + acturalSize, out_ptr[d]);
            }
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    // Non-size based operator.
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() == 3, "PhysDerivTemplate::Operator3D: Cannot "
                                     "call 3D routine with 1 or 2 outputs.");

        const auto nq0 = this->m_nq[0];
        const auto nq1 = this->m_nq[1];
        const auto nq2 = this->m_nq[2];

        const auto nqTot    = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;

        constexpr auto ndf = 9;

        const auto *inptr = input.data();
        auto *outptr_d0   = &output[0][0];
        auto *outptr_d1   = &output[1][0];
        auto *outptr_d2   = &output[2][0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        PhysDeriv3DWorkspace<SHAPE_TYPE>(nq0, nq1, nq2, wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            tmpIn(nqTot), tmpd0(nqTot), tmpd1(nqTot), tmpd2(nqTot);

        const vec_t *df_ptr = &((*this->m_df)[0]);
        vec_t df_tmp[ndf]   = {}; // ndf is a constexpr

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            // Results written to out_d0, out_d1, out_d2
            PhysDerivTensor3DKernel(nq0, nq1, nq2, tmpIn, this->m_D[0],
                                    this->m_D[1], this->m_D[2], tmpd0, tmpd1,
                                    tmpd2);
            // Calculate physical derivative
            PhysDeriv3DKernel<DEFORMED>(
                nq0, nq1, nq2, this->m_Z[0], this->m_Z[1], this->m_Z[2], df_ptr,
                df_tmp, wsp0, wsp1, tmpd0, tmpd1, tmpd2);

            // de-interleave and store data
            deinterleave_unalign_store(tmpd0, nqTot, outptr_d0);
            deinterleave_unalign_store(tmpd1, nqTot, outptr_d1);
            deinterleave_unalign_store(tmpd2, nqTot, outptr_d2);

            inptr += nqBlocks;
            df_ptr += dfSize;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            // Results written to out_d0, out_d1, out_d2
            PhysDerivTensor3DKernel(nq0, nq1, nq2, tmpIn, this->m_D[0],
                                    this->m_D[1], this->m_D[2], tmpd0, tmpd1,
                                    tmpd2);
            // Calculate physical derivative
            PhysDeriv3DKernel<DEFORMED>(
                nq0, nq1, nq2, this->m_Z[0], this->m_Z[1], this->m_Z[2], df_ptr,
                df_tmp, wsp0, wsp1, tmpd0, tmpd1, tmpd2);

            // De-interleave and store data
            deinterleave_store(tmpd0, nqTot, locField);
            std::copy(locField, locField + acturalSize, outptr_d0);
            deinterleave_store(tmpd1, nqTot, locField);
            std::copy(locField, locField + acturalSize, outptr_d1);
            deinterleave_store(tmpd2, nqTot, locField);
            std::copy(locField, locField + acturalSize, outptr_d2);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nq0, int nq1, int nq2>
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() == 3, "PhysDerivTemplate::Operator3D: Cannot "
                                     "call 3D routine with 1 or 2 outputs.");

        constexpr auto nqTot    = nq0 * nq1 * nq2;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        constexpr auto ndf = 9;

        const auto *inptr = input.data();
        auto *outptr_d0   = &output[0][0];
        auto *outptr_d1   = &output[1][0];
        auto *outptr_d2   = &output[2][0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        PhysDeriv3DWorkspace<SHAPE_TYPE>(nq0, nq1, nq2, wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpd0(nqTot),
            tmpd1(nqTot), tmpd2(nqTot), wsp0(wsp0Size), wsp1(wsp1Size);

        const vec_t *df_ptr = &((*this->m_df)[0]);
        vec_t df_tmp[ndf]   = {}; // ndf is a constexpr

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqBlocks];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            // Results written to out_d0, out_d1, out_d2
            PhysDerivTensor3DKernel(nq0, nq1, nq2, tmpIn, this->m_D[0],
                                    this->m_D[1], this->m_D[2], tmpd0, tmpd1,
                                    tmpd2);
            // Calculate physical derivative
            PhysDeriv3DKernel<DEFORMED>(
                nq0, nq1, nq2, this->m_Z[0], this->m_Z[1], this->m_Z[2], df_ptr,
                df_tmp, wsp0, wsp1, tmpd0, tmpd1, tmpd2);

            // de-interleave and store data
            deinterleave_unalign_store(tmpd0, nqTot, outptr_d0);
            deinterleave_unalign_store(tmpd1, nqTot, outptr_d1);
            deinterleave_unalign_store(tmpd2, nqTot, outptr_d2);

            inptr += nqBlocks;
            df_ptr += dfSize;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            // Results written to out_d0, out_d1, out_d2
            PhysDerivTensor3DKernel(nq0, nq1, nq2, tmpIn, this->m_D[0],
                                    this->m_D[1], this->m_D[2], tmpd0, tmpd1,
                                    tmpd2);
            // Calculate physical derivative
            PhysDeriv3DKernel<DEFORMED>(
                nq0, nq1, nq2, this->m_Z[0], this->m_Z[1], this->m_Z[2], df_ptr,
                df_tmp, wsp0, wsp1, tmpd0, tmpd1, tmpd2);

            // De-interleave and store data
            deinterleave_store(tmpd0, nqTot, locField);
            std::copy(locField, locField + acturalSize, outptr_d0);
            deinterleave_store(tmpd1, nqTot, locField);
            std::copy(locField, locField + acturalSize, outptr_d1);
            deinterleave_store(tmpd2, nqTot, locField);
            std::copy(locField, locField + acturalSize, outptr_d2);
        }
    }

#endif

private:
};

} // namespace Nektar::MatrixFree

#endif
