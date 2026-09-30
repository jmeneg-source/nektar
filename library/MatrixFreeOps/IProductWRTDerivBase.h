///////////////////////////////////////////////////////////////////////////////
//
// File: IProductWRTDerivBase.h
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

#ifndef NEKTAR_LIBRARY_MF_IPRODUCTWRTDERIVBASE_H
#define NEKTAR_LIBRARY_MF_IPRODUCTWRTDERIVBASE_H

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"

#include "IProductKernels.hpp"
#include "IProductWRTDerivBaseKernels.hpp"

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
struct IProductWRTDerivBaseTemplate : public IProductWRTDerivBase,
                                      public Helper<SHAPE_TYPE, DEFORMED>
{
    IProductWRTDerivBaseTemplate(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : IProductWRTDerivBase(basis, nElmt),
          Helper<SHAPE_TYPE, DEFORMED>(basis, nElmt)
    {
        constexpr auto DIM = LibUtilities::ShapeTypeDimMap[SHAPE_TYPE];

        if (DIM == 1)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(SHAPE_TYPE,
                                                            this->m_nm[0]);
        }
        else if (DIM == 2)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1]);
        }
        else if (DIM == 3)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1], this->m_nm[2]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<
            IProductWRTDerivBaseTemplate<SHAPE_TYPE, DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, Array<OneD, NekDouble>> &input,
                    Array<OneD, NekDouble> &output) final
    {
        const int nm0 = this->m_nm[0];
        const int nq0 = this->m_nq[0];
#if defined(SHAPE_DIMENSION_2D)
        const int nm1 = this->m_nm[1];
        const int nq1 = this->m_nq[1];
#elif defined(SHAPE_DIMENSION_3D)
        const int nm1 = this->m_nm[1];
        const int nm2 = this->m_nm[2];
        const int nq1 = this->m_nq[1];
        const int nq2 = this->m_nq[2];
#endif
#include "SwitchNodesPoints.h"
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
    void operator1D(const Array<OneD, Array<OneD, NekDouble>> &input,
                    Array<OneD, NekDouble> &output)
    {
        ASSERTL1(input.size() > 0, "IProductWRTDerivBaseTemplate::Operator1D: "
                                   "Cannot call 1D routine with no input.");
        ASSERTL1(input.size() <= 3, "IProductWRTDerivBaseTemplate::Operator1D: "
                                    "Operator not set up for other dimensions.")

        const auto nm0 = this->m_nm[0];
        const auto nq0 = this->m_nq[0];

        const auto nqTot    = nq0;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        const auto ndf        = input.size();
        constexpr int max_ndf = 3;

        auto *outptr = output.data();

        // Workspace for kernels - also checks preconditions
        IProduct1DWorkspace<SHAPE_TYPE>(nm0, nq0);

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;

        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf * nqTot;
        }

        std::vector<vec_t, allocator<vec_t>>
            tmpIn[max_ndf]; // max_ndf is a constexpr
        std::vector<vec_t, allocator<vec_t>> tmp0(nqTot), tmpOut(m_nmTot);

        const vec_t *jac_ptr = &((*this->m_jac)[0]);
        const vec_t *df_ptr  = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf];     // max_ndf is a constexpr
        NekDouble *inptr[max_ndf]; // max_ndf is a constexpr

        for (int d = 0; d < ndf; ++d)
        {
            tmpIn[d].resize(nqTot);
            inptr[d] = const_cast<NekDouble *>(input[d].data());
        }

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load and transpose data
            for (int d = 0; d < ndf; ++d)
            {
                // Load data to aligned storage and interleave it
                // std::copy(inptr[d], inptr[d] + nqBlocks, locField);
                // load_interleave(locField, nqTot, tmpIn[d]);
                load_unalign_interleave(inptr[d], nqTot, tmpIn[d]);
            }

            IProductWRTDerivBase1DKernel<SHAPE_TYPE, DEFORMED>(
                nqTot, ndf, df_ptr, df_tmp, tmpIn, tmp0);

            // IP DB0
            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmp0, this->m_dbdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            for (int d = 0; d < ndf; ++d)
            {
                inptr[d] += nqBlocks;
            }
            outptr += nmBlocks;
            df_ptr += dfSize;
            jac_ptr += dJSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            for (int d = 0; d < ndf; ++d)
            {
                std::copy(inptr[d], inptr[d] + acturalSize, locField);
                load_interleave(locField, nqTot, tmpIn[d]);
            }

            IProductWRTDerivBase1DKernel<SHAPE_TYPE, DEFORMED>(
                nqTot, ndf, df_ptr, df_tmp, tmpIn, tmp0);

            // IP DB0
            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmp0, this->m_dbdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nm0, int nq0>
    void operator1D(const Array<OneD, Array<OneD, NekDouble>> &input,
                    Array<OneD, NekDouble> &output)
    {
        ASSERTL1(input.size() > 0, "IProductWRTDerivBaseTemplate::Operator1D: "
                                   "Cannot call 1D routine with no input.");
        ASSERTL1(input.size() <= 3, "IProductWRTDerivBaseTemplate::Operator1D: "
                                    "Operator not set up for other dimensions.")

        constexpr auto nqTot    = nq0;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        const auto ndf        = input.size();
        constexpr int max_ndf = 3;

        auto *outptr = output.data();

        // Workspace for kernels - also checks preconditions
        IProduct1DWorkspace<SHAPE_TYPE>(nm0, nq0);

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;

        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf * nqTot;
        }

        std::vector<vec_t, allocator<vec_t>>
            tmpIn[max_ndf]; // max_ndf is a constexpr
        std::vector<vec_t, allocator<vec_t>> tmp0(nqTot), tmpOut(m_nmTot);

        const vec_t *jac_ptr = &((*this->m_jac)[0]);
        const vec_t *df_ptr  = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf];     // max_ndf is a constexpr
        NekDouble *inptr[max_ndf]; // max_ndf is a constexpr

        for (int d = 0; d < ndf; ++d)
        {
            tmpIn[d].resize(nqTot);
            inptr[d] = const_cast<NekDouble *>(input[d].data());
        }

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqTot * vec_t::width];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load and transpose data
            for (int d = 0; d < ndf; ++d)
            {
                // Load data to aligned storage and interleave it
                // std::copy(inptr[d], inptr[d] + nqBlocks, locField);
                // load_interleave(locField, nqTot, tmpIn[d]);
                load_unalign_interleave(inptr[d], nqTot, tmpIn[d]);
            }

            IProductWRTDerivBase1DKernel<SHAPE_TYPE, DEFORMED>(
                nqTot, ndf, df_ptr, df_tmp, tmpIn, tmp0);

            // IP DB0
            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmp0, this->m_dbdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            for (int d = 0; d < ndf; ++d)
            {
                inptr[d] += nqBlocks;
            }
            outptr += nmBlocks;
            df_ptr += dfSize;
            jac_ptr += dJSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            for (int d = 0; d < ndf; ++d)
            {
                std::copy(inptr[d], inptr[d] + acturalSize, locField);
                load_interleave(locField, nqTot, tmpIn[d]);
            }

            IProductWRTDerivBase1DKernel<SHAPE_TYPE, DEFORMED>(
                nqTot, ndf, df_ptr, df_tmp, tmpIn, tmp0);

            // IP DB0
            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmp0, this->m_dbdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

#elif defined(SHAPE_DIMENSION_2D)

    // Non-size based operator.
    void operator2D(const Array<OneD, Array<OneD, NekDouble>> &input,
                    Array<OneD, NekDouble> &output)
    {
        ASSERTL1(input.size() > 1, "IProductWRTDerivBaseTemplate::Operator2D: "
                                   "Cannot call 2D routine with one output.");
        ASSERTL1(input.size() <= 3, "IProductWRTDerivBaseTemplate::Operator2D: "
                                    "Operator not set up for 3D coordinates.");

        const auto nm0 = this->m_nm[0];
        const auto nm1 = this->m_nm[1];

        const auto nq0 = this->m_nq[0];
        const auto nq1 = this->m_nq[1];

        const auto nqTot    = nq0 * nq1;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        const auto indim        = input.size();
        constexpr int max_indim = 3;

        const auto ndf        = 2 * indim;
        constexpr int max_ndf = 2 * max_indim;

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        auto *outptr = output.data();

        std::vector<vec_t, allocator<vec_t>>
            tmpIn[max_indim], // max_indim is a constexpr
            tmp0(nqTot), tmp1(nqTot), tmpOut(m_nmTot);

        const NekDouble *inptr[max_indim]; // max_indim is a constexpr

        for (int d = 0; d < indim; ++d)
        {
            tmpIn[d].resize(nqTot);
            inptr[d] = &input[d][0];
        }

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;

        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf * nqTot;
        }

        // Workspace for kernels
        size_t wsp0Size = 0;
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size);

        const vec_t *jac_ptr = &((*this->m_jac)[0]);
        const vec_t *df_ptr  = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf]; // max_ndf is a constexpr

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load and transpose data
            for (int d = 0; d < indim; ++d)
            {
                // Load data to aligned storage and interleave it
                // std::copy(inptr[d], inptr[d] + nqBlocks, locField);
                // load_interleave(locField, nqTot, tmpIn[d]);
                load_unalign_interleave(inptr[d], nqTot, tmpIn[d]);
            }

            fusedKernel2D(nm0, nm1, nq0, nq1, indim, correct, jac_ptr, df_ptr,
                          tmpIn, df_tmp, tmp0, tmp1, wsp0, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            for (int d = 0; d < indim; ++d)
            {
                inptr[d] += nqBlocks;
            }
            outptr += nmBlocks;
            df_ptr += dfSize;
            jac_ptr += dJSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            for (int d = 0; d < indim; ++d)
            {
                std::copy(inptr[d], inptr[d] + acturalSize, locField);
                load_interleave(locField, nqTot, tmpIn[d]);
            }

            fusedKernel2D(nm0, nm1, nq0, nq1, indim, correct, jac_ptr, df_ptr,
                          tmpIn, df_tmp, tmp0, tmp1, wsp0, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nm0, int nm1, int nq0, int nq1>
    void operator2D(const Array<OneD, Array<OneD, NekDouble>> &input,
                    Array<OneD, NekDouble> &output)
    {
        ASSERTL1(input.size() > 1, "IProductWRTDerivBaseTemplate::Operator2D: "
                                   "Cannot call 2D routine with one output.");
        ASSERTL1(input.size() <= 3, "IProductWRTDerivBaseTemplate::Operator2D: "
                                    "Operator not set up for 3D coordinates.");

        constexpr auto nqTot    = nq0 * nq1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        const auto indim        = input.size();
        constexpr int max_indim = 3;

        const auto ndf        = 2 * indim;
        constexpr int max_ndf = 2 * max_indim;

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        auto *outptr = output.data();

        std::vector<vec_t, allocator<vec_t>>
            tmpIn[max_indim], // max_indim is a constexpr
            tmp0(nqTot), tmp1(nqTot), tmpOut(m_nmTot);

        const NekDouble *inptr[max_indim]; // max_ndf is a constexpr

        for (int d = 0; d < indim; ++d)
        {
            tmpIn[d].resize(nqTot);
            inptr[d] = &input[d][0];
        }

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;

        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf * nqTot;
        }

        // Workspace for kernels
        size_t wsp0Size = 0;
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size);

        const vec_t *jac_ptr = &((*this->m_jac)[0]);
        const vec_t *df_ptr  = &((*this->m_df)[0]);
        vec_t df_tmp[max_ndf]; // max_ndf is a constexpr

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqTot * vec_t::width];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load and transpose data
            for (int d = 0; d < indim; ++d)
            {
                // Load data to aligned storage and interleave it
                // std::copy(inptr[d], inptr[d] + nqBlocks, locField);
                // load_interleave(locField, nqTot, tmpIn[d]);
                load_unalign_interleave(inptr[d], nqTot, tmpIn[d]);
            }

            fusedKernel2D(nm0, nm1, nq0, nq1, indim, correct, jac_ptr, df_ptr,
                          tmpIn, df_tmp, tmp0, tmp1, wsp0, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            for (int d = 0; d < indim; ++d)
            {
                inptr[d] += nqBlocks;
            }
            outptr += nmBlocks;
            df_ptr += dfSize;
            jac_ptr += dJSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            for (int d = 0; d < indim; ++d)
            {
                std::copy(inptr[d], inptr[d] + acturalSize, locField);
                load_interleave(locField, nqTot, tmpIn[d]);
            }

            fusedKernel2D(nm0, nm1, nq0, nq1, indim, correct, jac_ptr, df_ptr,
                          tmpIn, df_tmp, tmp0, tmp1, wsp0, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

    NEK_FORCE_INLINE void fusedKernel2D(
        const size_t nm0, const size_t nm1, const size_t nq0, const size_t nq1,
        const size_t indim, const bool correct, const vec_t *jac_ptr,
        const vec_t *df_ptr, const std::vector<vec_t, allocator<vec_t>> *tmpIn,
        vec_t *df_tmp, std::vector<vec_t, allocator<vec_t>> &tmp0,
        std::vector<vec_t, allocator<vec_t>> &tmp1,
        std::vector<vec_t, allocator<vec_t>> &wsp0,
        std::vector<vec_t, allocator<vec_t>> &tmpOut)
    {
        IProductWRTDerivBase2DKernel<SHAPE_TYPE, DEFORMED>(
            nq0, nq1, indim, df_ptr, df_tmp, this->m_Z[0], this->m_Z[1], tmpIn,
            tmp0, tmp1);

        // IP DB0 B1
        IProduct2DKernel<SHAPE_TYPE, false, false, DEFORMED>(
            nm0, nm1, nq0, nq1, correct, tmp0, this->m_dbdata[0],
            this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
            tmpOut);

        // IP DB1 B0
        IProduct2DKernel<SHAPE_TYPE, false, true, DEFORMED>(
            nm0, nm1, nq0, nq1, correct, tmp1, this->m_bdata[0],
            this->m_dbdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
            tmpOut);
    }

#elif defined(SHAPE_DIMENSION_3D)

    // Non-size based operator.
    void operator3D(const Array<OneD, Array<OneD, NekDouble>> &input,
                    Array<OneD, NekDouble> &output)
    {
        ASSERTL1(input.size() == 3,
                 "IProductWRTDerivBaseTemplate::Operator3D: Cannot call 3D "
                 "routine with 1 or 2 outputs.");

        const auto nm0 = this->m_nm[0];
        const auto nm1 = this->m_nm[1];
        const auto nm2 = this->m_nm[2];

        const auto nq0 = this->m_nq[0];
        const auto nq1 = this->m_nq[1];
        const auto nq2 = this->m_nq[2];

        constexpr auto ndf  = 9u;
        const auto nqTot    = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        auto *inptr0 = input[0].data();
        auto *inptr1 = input[1].data();
        auto *inptr2 = input[2].data();
        auto *outptr = output.data();

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;

        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf * nqTot;
        }

        // Workspace for kernels
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size), tmpIn0(nqTot), tmpIn1(nqTot), tmpIn2(nqTot),
            tmp0(nqTot), tmp1(nqTot), tmp2(nqTot), tmpOut(m_nmTot);

        const vec_t *jac_ptr = &((*this->m_jac)[0]);
        const vec_t *df_ptr  = &((*this->m_df)[0]);

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // load_interleave(inptr0, nqTot, tmpIn0);
            // load_interleave(inptr1, nqTot, tmpIn1);
            // load_interleave(inptr2, nqTot, tmpIn2);
            load_unalign_interleave(inptr0, nqTot, tmpIn0);
            load_unalign_interleave(inptr1, nqTot, tmpIn1);
            load_unalign_interleave(inptr2, nqTot, tmpIn2);

            fusedKernel3D(nm0, nm1, nm2, nq0, nq1, nq2, correct, jac_ptr,
                          df_ptr, tmpIn0, tmpIn1, tmpIn2, tmp0, tmp1, tmp2,
                          wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr0 += nqBlocks;
            inptr1 += nqBlocks;
            inptr2 += nqBlocks;
            outptr += nmBlocks;
            df_ptr += dfSize;
            jac_ptr += dJSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr0, inptr0 + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn0);
            std::copy(inptr1, inptr1 + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn1);
            std::copy(inptr2, inptr2 + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn2);

            fusedKernel3D(nm0, nm1, nm2, nq0, nq1, nq2, correct, jac_ptr,
                          df_ptr, tmpIn0, tmpIn1, tmpIn2, tmp0, tmp1, tmp2,
                          wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nm0, int nm1, int nm2, int nq0, int nq1, int nq2>
    void operator3D(const Array<OneD, Array<OneD, NekDouble>> &input,
                    Array<OneD, NekDouble> &output)
    {
        ASSERTL1(input.size() == 3,
                 "IProductWRTDerivBaseTemplate::Operator3D: Cannot call 3D "
                 "routine with 1 or 2 outputs.");

        constexpr auto ndf      = 9u;
        constexpr auto nqTot    = nq0 * nq1 * nq2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        auto *inptr0 = input[0].data();
        auto *inptr1 = input[1].data();
        auto *inptr2 = input[2].data();
        auto *outptr = output.data();

        // Get size of jacobian factor block
        auto dJSize = 1u;
        auto dfSize = ndf;

        if (DEFORMED)
        {
            dJSize = nqTot;
            dfSize = ndf * nqTot;
        }

        // Workspace for kernels
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size), tmpIn0(nqTot), tmpIn1(nqTot), tmpIn2(nqTot),
            tmp0(nqTot), tmp1(nqTot), tmp2(nqTot), tmpOut(m_nmTot);

        const vec_t *jac_ptr = &((*this->m_jac)[0]);
        const vec_t *df_ptr  = &((*this->m_df)[0]);

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqTot * vec_t::width];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load and transpose data
            // load_interleave(inptr0, nqTot, tmpIn0);
            // load_interleave(inptr1, nqTot, tmpIn1);
            // load_interleave(inptr2, nqTot, tmpIn2);
            load_unalign_interleave(inptr0, nqTot, tmpIn0);
            load_unalign_interleave(inptr1, nqTot, tmpIn1);
            load_unalign_interleave(inptr2, nqTot, tmpIn2);

            fusedKernel3D(nm0, nm1, nm2, nq0, nq1, nq2, correct, jac_ptr,
                          df_ptr, tmpIn0, tmpIn1, tmpIn2, tmp0, tmp1, tmp2,
                          wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr0 += nqBlocks;
            inptr1 += nqBlocks;
            inptr2 += nqBlocks;
            outptr += nmBlocks;
            df_ptr += dfSize;
            jac_ptr += dJSize;
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr0, inptr0 + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn0);
            std::copy(inptr1, inptr1 + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn1);
            std::copy(inptr2, inptr2 + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn2);

            fusedKernel3D(nm0, nm1, nm2, nq0, nq1, nq2, correct, jac_ptr,
                          df_ptr, tmpIn0, tmpIn1, tmpIn2, tmp0, tmp1, tmp2,
                          wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

    NEK_FORCE_INLINE void fusedKernel3D(
        const size_t nm0, const size_t nm1, const size_t nm2, const size_t nq0,
        const size_t nq1, const size_t nq2, const bool correct,
        const vec_t *jac_ptr, const vec_t *df_ptr,
        const std::vector<vec_t, allocator<vec_t>> &tmpIn0,
        const std::vector<vec_t, allocator<vec_t>> &tmpIn1,
        const std::vector<vec_t, allocator<vec_t>> &tmpIn2,
        std::vector<vec_t, allocator<vec_t>> &tmp0,
        std::vector<vec_t, allocator<vec_t>> &tmp1,
        std::vector<vec_t, allocator<vec_t>> &tmp2,
        std::vector<vec_t, allocator<vec_t>> &wsp0,
        std::vector<vec_t, allocator<vec_t>> &wsp1,
        std::vector<vec_t, allocator<vec_t>> &wsp2,
        std::vector<vec_t, allocator<vec_t>> &tmpOut)
    {
        IProductWRTDerivBase3DKernel<SHAPE_TYPE, DEFORMED>(
            nq0, nq1, nq2, df_ptr, this->m_Z[0], this->m_Z[1], this->m_Z[2],
            tmpIn0, tmpIn1, tmpIn2, tmp0, tmp1, tmp2);

        // IP DB0 B1 B2
        IProduct3DKernel<SHAPE_TYPE, false, false, DEFORMED>(
            nm0, nm1, nm2, nq0, nq1, nq2, correct, tmp0, this->m_dbdata[0],
            this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
            this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

        // IP B0 DB1 B2
        IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
            nm0, nm1, nm2, nq0, nq1, nq2, correct, tmp1, this->m_bdata[0],
            this->m_dbdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
            this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

        // IP B0 B1 DB2
        IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
            nm0, nm1, nm2, nq0, nq1, nq2, correct, tmp2, this->m_bdata[0],
            this->m_bdata[1], this->m_dbdata[2], this->m_w[0], this->m_w[1],
            this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);
    }

#endif // SHAPE_DIMENSION

private:
    int m_nmTot;
};

} // namespace Nektar::MatrixFree

#endif
