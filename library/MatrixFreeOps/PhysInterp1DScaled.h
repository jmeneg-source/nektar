///////////////////////////////////////////////////////////////////////////////
//
// File: PhysInterp1DScaled.h
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

#ifndef NEKTAR_LIBRARY_MF_PHYSINTERP1DSCALED_H
#define NEKTAR_LIBRARY_MF_PHYSINTERP1DSCALED_H

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"
#include "PhysInterp1DScaledKernels.hpp"

// As each opertor has seven shapes over three dimension to get to the
// "work" each operator uses a series of preprocessor directives based
// on the dimension and shape type so to limit the code
// inclusion. This keeps the library size as minimal as possible while
// removing runtime conditionals.

// The preprocessor directives, SHAPE_DIMENSION_?D and SHAPE_TYPE_*
// are constructed by CMake in the CMakeLists.txt file which uses an
// implementation file.  See the CMakeLists.txt files for more
// details.
namespace Nektar::MatrixFree
{

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED = false>
struct PhysInterp1DScaledTemplate : public PhysInterp1DScaled,
                                    public Helper<SHAPE_TYPE>
{
    PhysInterp1DScaledTemplate(std::vector<LibUtilities::BasisSharedPtr> basis,
                               int nElmt)
        : PhysInterp1DScaled(basis, nElmt), Helper<SHAPE_TYPE>(basis, nElmt)
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<PhysInterp1DScaledTemplate<SHAPE_TYPE>>(basis,
                                                                        nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output) final
    {
        // Following the same variable naming as in SwitchNodesPoints.h
        // Both input and output number of points of the kernels are in the
        // physical space
        // nm0 is replaced by nq_in0
        // nq0 is replaced by nq_out0
        const auto nm0 = this->m_nq[0];
        const auto nq0 = this->m_enhancednq[0];
#if defined(SHAPE_DIMENSION_2D)
        // nm1 is replaced by nq_in1
        // nq1 is replaced by nq_out1
        // Following the same variable naming as in SwitchNodesPoints.h
        const auto nm1 = this->m_nq[1];
        const auto nq1 = this->m_enhancednq[1];
#elif defined(SHAPE_DIMENSION_3D)
        // nm2 is replaced by nq_in2
        // nq2 is replaced by nq_out2
        // Following the same variable naming as in SwitchNodesPoints.h
        const auto nm1 = this->m_nq[1];
        const auto nm2 = this->m_nq[2];
        const auto nq1 = this->m_enhancednq[1];
        const auto nq2 = this->m_enhancednq[2];
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
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nq_in0  = this->m_nq[0];
        const auto nq_out0 = this->m_enhancednq[0];

        const auto nb_out = nq_out0 * vec_t::width;
        const auto nb_in  = nq_in0 * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        PhysInterp1DScaled1DWorkspace<SHAPE_TYPE>(nq_in0, nq_out0);

        std::vector<vec_t, allocator<vec_t>> tmpIn(nq_in0), tmpOut(nq_out0);

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nb_out * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nb_in, locField);
            // load_interleave(locField, nq_in0, tmpIn);
            load_unalign_interleave(inptr, nq_in0, tmpIn);

            PhysInterp1DScaled1DKernel<SHAPE_TYPE>(nq_in0, nq_out0, tmpIn,
                                                   this->m_I[0], tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, nq_out0, locField);
            // std::copy(locField, locField + nb_out, outptr);
            deinterleave_unalign_store(tmpOut, nq_out0, outptr);

            inptr += nb_in;
            outptr += nb_out;
        }
        // last block
        {
            int acturalSize = nb_in - this->m_nPads * nq_in0;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nq_in0, tmpIn);

            PhysInterp1DScaled1DKernel<SHAPE_TYPE>(nq_in0, nq_out0, tmpIn,
                                                   this->m_I[0], tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nq_out0, locField);
            acturalSize = nb_out - this->m_nPads * nq_out0;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nq_in0, int nq_out0>
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nOutTot = nq_out0;
        constexpr auto nInTot  = nq_in0;
        constexpr auto nb_out  = nOutTot * vec_t::width;
        const auto nb_in       = nInTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        PhysInterp1DScaled1DWorkspace<SHAPE_TYPE>(nq_in0, nq_out0);

        std::vector<vec_t, allocator<vec_t>> tmpIn(nInTot), tmpOut(nOutTot);

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nb_out];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nb_in, locField);
            // load_interleave(locField, m_nInTot, tmpIn);
            load_unalign_interleave(inptr, nInTot, tmpIn);

            PhysInterp1DScaled1DKernel<SHAPE_TYPE>(nq_in0, nq_out0, tmpIn,
                                                   this->m_I[0], tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, nOutTot, locField);
            // std::copy(locField, locField + nb_out, outptr);
            deinterleave_unalign_store(tmpOut, nOutTot, outptr);

            inptr += nb_in;
            outptr += nb_out;
        }
        // last block
        {
            int acturalSize = nb_in - this->m_nPads * nInTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nInTot, tmpIn);

            PhysInterp1DScaled1DKernel<SHAPE_TYPE>(nq_in0, nq_out0, tmpIn,
                                                   this->m_I[0], tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nOutTot, locField);
            acturalSize = nb_out - this->m_nPads * nOutTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

#elif defined(SHAPE_DIMENSION_2D)

    // Non-size based operator.
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nq_in0 = this->m_nq[0];
        const auto nq_in1 = this->m_nq[1];

        const auto nq_out0 = this->m_enhancednq[0];
        const auto nq_out1 = this->m_enhancednq[1];

        const auto nInTot  = nq_in0 * nq_in1;
        const auto nOutTot = nq_out0 * nq_out1;
        const auto nb_out  = nOutTot * vec_t::width;
        const auto nb_in   = nInTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        PhysInterp1DScaled2DWorkspace<SHAPE_TYPE>(nq_in0, nq_in1, nq_out0,
                                                  nq_out1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(nInTot),
            tmpOut(nOutTot);

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nb_out * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nb_in, locField);
            // load_interleave(locField, m_nInTot, tmpIn);
            load_unalign_interleave(inptr, nInTot, tmpIn);

            PhysInterp1DScaled2DKernel<SHAPE_TYPE>(nq_in0, nq_in1, nq_out0,
                                                   nq_out1, tmpIn, this->m_I[0],
                                                   this->m_I[1], wsp0, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, nOutTot, locField);
            // std::copy(locField, locField + nb_out, outptr);
            deinterleave_unalign_store(tmpOut, nOutTot, outptr);

            inptr += nb_in;
            outptr += nb_out;
        }
        // last block
        {
            int acturalSize = nb_in - this->m_nPads * nInTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nInTot, tmpIn);

            PhysInterp1DScaled2DKernel<SHAPE_TYPE>(nq_in0, nq_in1, nq_out0,
                                                   nq_out1, tmpIn, this->m_I[0],
                                                   this->m_I[1], wsp0, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nOutTot, locField);
            acturalSize = nb_out - this->m_nPads * nOutTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nq_in0, int nq_in1, int nq_out0, int nq_out1>
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nInTot  = nq_in0 * nq_in1;
        constexpr auto nOutTot = nq_out0 * nq_out1;
        constexpr auto nb_out  = nOutTot * vec_t::width;
        const auto nb_in       = nInTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        PhysInterp1DScaled2DWorkspace<SHAPE_TYPE>(nq_in0, nq_in1, nq_out0,
                                                  nq_out1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(nInTot),
            tmpOut(nOutTot);

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nb_out];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nb_in, locField);
            // load_interleave(locField, m_nInTot, tmpIn);
            load_unalign_interleave(inptr, nInTot, tmpIn);

            PhysInterp1DScaled2DKernel<SHAPE_TYPE>(nq_in0, nq_in1, nq_out0,
                                                   nq_out1, tmpIn, this->m_I[0],
                                                   this->m_I[1], wsp0, tmpOut);
            // de-interleave and store data
            // deinterleave_store(tmpOut, nOutTot, locField);
            // std::copy(locField, locField + nb_out, outptr);
            deinterleave_unalign_store(tmpOut, nOutTot, outptr);

            inptr += nb_in;
            outptr += nb_out;
        }
        // last block
        {
            int acturalSize = nb_in - this->m_nPads * nInTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nInTot, tmpIn);

            PhysInterp1DScaled2DKernel<SHAPE_TYPE>(nq_in0, nq_in1, nq_out0,
                                                   nq_out1, tmpIn, this->m_I[0],
                                                   this->m_I[1], wsp0, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nOutTot, locField);
            acturalSize = nb_out - this->m_nPads * nOutTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    // Non-size based operator.
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nq_in0 = this->m_nq[0];
        const auto nq_in1 = this->m_nq[1];
        const auto nq_in2 = this->m_nq[2];

        const auto nq_out0 = this->m_enhancednq[0];
        const auto nq_out1 = this->m_enhancednq[1];
        const auto nq_out2 = this->m_enhancednq[2];

        const auto nInTot  = nq_in0 * nq_in1 * nq_in2;
        const auto nOutTot = nq_out0 * nq_out1 * nq_out2;
        const auto nb_out  = nOutTot * vec_t::width;
        const auto nb_in   = nInTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        PhysInterp1DScaled3DWorkspace<SHAPE_TYPE>(nq_in0, nq_in1, nq_in2,
                                                  nq_out0, nq_out1, nq_out2,
                                                  wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            tmpIn(nInTot), tmpOut(nOutTot);

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nb_out * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nb_in, locField);
            // load_interleave(locField, m_nInTot, tmpIn);
            load_unalign_interleave(inptr, nInTot, tmpIn);

            PhysInterp1DScaled3DKernel<SHAPE_TYPE>(
                nq_in0, nq_in1, nq_in2, nq_out0, nq_out1, nq_out2, tmpIn,
                this->m_I[0], this->m_I[1], this->m_I[2], wsp0, wsp1, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, nOutTot, locField);
            // std::copy(locField, locField + nb_out, outptr);
            deinterleave_unalign_store(tmpOut, nOutTot, outptr);

            inptr += nb_in;
            outptr += nb_out;
        }
        // last block
        {
            int acturalSize = nb_in - this->m_nPads * nInTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nInTot, tmpIn);

            PhysInterp1DScaled3DKernel<SHAPE_TYPE>(
                nq_in0, nq_in1, nq_in2, nq_out0, nq_out1, nq_out2, tmpIn,
                this->m_I[0], this->m_I[1], this->m_I[2], wsp0, wsp1, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nOutTot, locField);
            acturalSize = nb_out - this->m_nPads * nOutTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nq_in0, int nq_in1, int nq_in2, int nq_out0, int nq_out1,
              int nq_out2>
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nInTot  = nq_in0 * nq_in1 * nq_in2;
        constexpr auto nOutTot = nq_out0 * nq_out1 * nq_out2;
        constexpr auto nb_out  = nOutTot * vec_t::width;
        const auto nb_in       = nInTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        PhysInterp1DScaled3DWorkspace<SHAPE_TYPE>(nq_in0, nq_in1, nq_in2,
                                                  nq_out0, nq_out1, nq_out2,
                                                  wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            tmpIn(nInTot), tmpOut(nOutTot);

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nb_out];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nb_in, locField);
            // load_interleave(locField, m_nInTot, tmpIn);
            load_unalign_interleave(inptr, nInTot, tmpIn);

            PhysInterp1DScaled3DKernel<SHAPE_TYPE>(
                nq_in0, nq_in1, nq_in2, nq_out0, nq_out1, nq_out2, tmpIn,
                this->m_I[0], this->m_I[1], this->m_I[2], wsp0, wsp1, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, nOutTot, locField);
            // std::copy(locField, locField + nb_out, outptr);
            deinterleave_unalign_store(tmpOut, nOutTot, outptr);

            inptr += nb_in;
            outptr += nb_out;
        }
        // last block
        {
            int acturalSize = nb_in - this->m_nPads * nInTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nInTot, tmpIn);

            PhysInterp1DScaled3DKernel<SHAPE_TYPE>(
                nq_in0, nq_in1, nq_in2, nq_out0, nq_out1, nq_out2, tmpIn,
                this->m_I[0], this->m_I[1], this->m_I[2], wsp0, wsp1, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nOutTot, locField);
            acturalSize = nb_out - this->m_nPads * nOutTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

#endif // SHAPE_DIMENSION

private:
};

} // namespace Nektar::MatrixFree

#endif
