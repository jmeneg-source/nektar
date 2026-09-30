///////////////////////////////////////////////////////////////////////////////
//
// File: BwdTransSumFacStdKernels.hpp
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

template <typename simd_type>
NEK_FORCE_INLINE static void BwdTransSegKernel(const unsigned int nm0,
                                               const unsigned int nq0,
                                               const simd_type *basis0,
                                               const simd_type *in,
                                               simd_type *out)
{
    for (unsigned int i = 0; i < nq0; ++i)
    {
        simd_type tmp = in[0] * basis0[i]; // Load 2x

        for (unsigned int p = 1; p < nm0; ++p)
        {
            tmp.fma(in[p], basis0[p * nq0 + i]); // Load 2x
        }

        out[i] = tmp; // Store 1x
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void BwdTransTriKernel(
    const unsigned int nm0, const unsigned int nm1, const unsigned int nq0,
    const unsigned int nq1, const bool isModified, const simd_type *basis0,
    const simd_type *basis1, simd_type *p_sums, const simd_type *in,
    simd_type *out) // wsp : nm0
{
    for (unsigned int eta1 = 0, eta_idx = 0; eta1 < nq1; ++eta1)
    {
        for (unsigned int p = 0, mode = 0; p < nm0; ++p)
        {
            simd_type p_sum = 0.0;

            for (unsigned int q = 0; q < (nm1 - p); ++q, ++mode)
            {
                p_sum.fma(basis1[mode * nq1 + eta1], in[mode]);
            }

            p_sums[p] = p_sum; // Store 1x
        }

        // Already have q_sum at each quadrature point in eta 1 for
        // each mode, p.  From this assemble the tensor produce of
        // each quadrature point, eta1
        for (unsigned int eta0 = 0; eta0 < nq0; ++eta0, ++eta_idx)
        {
            simd_type p_sum = 0.0;
            for (unsigned int p = 0; p < nm0; ++p)
            {
                p_sum.fma(p_sums[p], basis0[p * nq0 + eta0]); // Load 2x
            }

            if (isModified)
            {
                // p_sum += coef * basis0 * basis1
                p_sum.fma(in[1] * basis0[nq0 + eta0], basis1[nq1 + eta1]);
            }

            out[eta_idx] = p_sum; // Store 1x
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void BwdTransQuadKernel(
    const unsigned int nm0, const unsigned int nm1, const unsigned int nq0,
    const unsigned int nq1, const simd_type *basis0, const simd_type *basis1,
    simd_type *wsp, const simd_type *in, simd_type *out) // wsp : nq0 * nm1
{
    for (unsigned int i = 0, cnt_iq = 0; i < nq0; ++i)
    {
        for (unsigned int q = 0, cnt_pq = 0; q < nm1; ++q, ++cnt_iq)
        {
            simd_type tmp = in[cnt_pq] * basis0[i]; // Load 2x
            ++cnt_pq;
            for (unsigned int p = 1; p < nm0; ++p, ++cnt_pq)
            {
                tmp.fma(in[cnt_pq], basis0[p * nq0 + i]); // Load 2x
            }
            wsp[cnt_iq] = tmp; // Store 1x
        }
    }

    for (unsigned int j = 0, cnt_ij = 0; j < nq1; ++j)
    {
        for (unsigned int i = 0, cnt_iq = 0; i < nq0; ++i, ++cnt_ij)
        {
            simd_type tmp = wsp[cnt_iq] * basis1[j]; // Load 2x
            ++cnt_iq;
            for (unsigned int q = 1; q < nm1; ++q, ++cnt_iq)
            {
                tmp.fma(wsp[cnt_iq], basis1[q * nq1 + j]); // Load 2x
            }
            out[cnt_ij] = tmp; // Store 1x
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void BwdTransHexKernel(
    const unsigned int nm0, const unsigned int nm1, const unsigned int nm2,
    const unsigned int nq0, const unsigned int nq1, const unsigned int nq2,
    const simd_type *basis0, const simd_type *basis1, const simd_type *basis2,
    simd_type *sum_irq, // nq0 * nm2 * nm1
    simd_type *sum_jir, // nq1 * nq0 * nm2
    const simd_type *in, simd_type *out)
{
    for (unsigned int i = 0, cnt_irq = 0; i < nq0; ++i)
    {
        for (unsigned int r = 0, cnt_rqp = 0; r < nm2; ++r)
        {
            for (unsigned int q = 0; q < nm1; ++q, ++cnt_irq)
            {
                simd_type tmp = in[cnt_rqp] * basis0[i];
                ++cnt_rqp;

                for (unsigned int p = 1; p < nm0; ++p, ++cnt_rqp)
                {
                    tmp.fma(in[cnt_rqp], basis0[p * nq0 + i]);
                }

                sum_irq[cnt_irq] = tmp;
            }
        }
    }

    for (unsigned int j = 0, cnt_jir = 0; j < nq1; ++j)
    {
        for (unsigned int i = 0, cnt_irq = 0; i < nq0; ++i)
        {
            for (unsigned int r = 0; r < nm2; ++r, ++cnt_jir)
            {
                simd_type tmp = sum_irq[cnt_irq] * basis1[j];
                ++cnt_irq;

                for (unsigned int q = 1; q < nm1; ++q)
                {
                    tmp.fma(sum_irq[cnt_irq++], basis1[q * nq1 + j]);
                }

                sum_jir[cnt_jir] = tmp;
            }
        }
    }

    for (unsigned int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        for (unsigned int j = 0, cnt_jir = 0; j < nq1; ++j)
        {
            for (unsigned int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                simd_type tmp = sum_jir[cnt_jir] * basis2[k];
                ++cnt_jir;

                for (unsigned int r = 1; r < nm2; ++r)
                {
                    tmp.fma(sum_jir[cnt_jir++], basis2[r * nq2 + k]);
                }

                out[cnt_kji] = tmp;
            }
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void BwdTransTetKernel(
    const unsigned int nm0, const unsigned int nm1, const unsigned int nm2,
    const unsigned int nq0, const unsigned int nq1, const unsigned int nq2,
    const bool isModified, const simd_type *basis0, const simd_type *basis1,
    const simd_type *basis2,
    simd_type *fpq, // nm0 * nm1
    simd_type *fp,  // nm0
    const simd_type *in, simd_type *out)
{
    for (unsigned int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        unsigned int cnt_pq = 0, cnt_pqr = 0, mode = 0;
        for (unsigned int p = 0; p < nm0; ++p)
        {
            for (unsigned int q = 0; q < nm1 - p; ++q, ++cnt_pq)
            {
                simd_type prod =
                    in[cnt_pqr] * basis2[k + nq2 * mode]; // Load 2x
                ++mode;
                ++cnt_pqr;

                for (unsigned int r = 1; r < nm2 - p - q;
                     ++r, ++mode, ++cnt_pqr)
                {
                    prod.fma(in[cnt_pqr], basis2[k + nq2 * mode]); // Load 2x
                }

                fpq[cnt_pq] = prod; // Store 1x
            }

            // increment mode in case order1!=order2
            for (unsigned int q = nm1 - p; q < nm2 - p; ++q)
            {
                mode += nm2 - p - q;
            }
        }

        for (unsigned int j = 0; j < nq1; ++j)
        {
            mode = cnt_pq = 0;
            for (unsigned int p = 0; p < nm0; ++p)
            {
                simd_type prod =
                    fpq[cnt_pq] * basis1[mode * nq1 + j]; // Load 2x
                ++cnt_pq;

                for (unsigned int q = 1; q < nm1 - p; ++q, ++cnt_pq)
                {
                    prod.fma(fpq[cnt_pq],
                             basis1[(mode + q) * nq1 + j]); // Load 2x
                }

                fp[p] = prod; // Store 1x
                mode += nm1 - p;
            }

            for (unsigned int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                simd_type tmp = basis0[i] * fp[0]; // Load 2x

                for (unsigned int p = 1; p < nm0; ++p)
                {
                    tmp.fma(basis0[p * nq0 + i], fp[p]); // Load 2x
                }

                if (isModified)
                {
                    // top vertex
                    //
                    // sum += inarray[1] * base2[nquad2 + k] * (
                    //     base0[i] * base1[nquad1+j] +
                    //     base0[nquad0+i] * base1[j] +
                    //     base0[nquad0+i] * base1[nquad1+j]);

                    simd_type tmp1 = basis0[i] * basis1[nq1 + j]; // Load 2x
                    tmp1.fma(basis0[nq0 + i], basis1[j]);         // Load 2x
                    tmp1.fma(basis0[nq0 + i], basis1[nq1 + j]);   // Load 2x
                    tmp1 = tmp1 * basis2[nq2 + k];                // Load 1x

                    simd_type inarray1 = simd_type(in[1]); // Load 1x
                    tmp.fma(tmp1, inarray1);

                    // bottom vertex
                    //
                    // sum += inarray[order2] * base2[k] * (
                    //     base0[nquad0+i] * base1[nquad1+j]);
                    tmp1     = basis0[nq0 + i] * basis1[nq1 + j]; // Load 2x
                    tmp1     = tmp1 * basis2[k];                  // Load 1x
                    inarray1 = simd_type(in[nm2]);                // Load 1x
                    tmp.fma(inarray1, tmp1);

                    // singular edge
                    for (unsigned int r = 1; r < nm2 - 1; ++r)
                    {
                        // sum += inarray[order2+r] * base2[(r+1)*nquad2+k] *
                        //     base1[nquad1+j] * base0[nquad0+i];
                        tmp1     = basis1[nq1 + j] * basis0[nq0 + i]; // Load 2x
                        tmp1     = tmp1 * basis2[(r + 1) * nq2 + k];  // Load 1x
                        inarray1 = simd_type(in[nm2 + r]);            // Load 1x
                        tmp.fma(inarray1, tmp1);
                        // multiply by (1-a)/2
                    }
                }

                out[cnt_kji] = tmp; // Store 1x
            }
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void BwdTransPrismKernel(
    const unsigned int nm0, const unsigned int nm1, const unsigned int nm2,
    const unsigned int nq0, const unsigned int nq1, const unsigned int nq2,
    const bool isModified, const simd_type *basis0, const simd_type *basis1,
    const simd_type *basis2,
    simd_type *fpq, // nm0 * nm1
    simd_type *fp,  // nm0
    const simd_type *in, simd_type *out)
{
    for (unsigned int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        unsigned int mode_pqr = 0, mode_pq = 0, mode_pr = 0;
        for (unsigned int p = 0; p < nm0; ++p)
        {
            for (unsigned int q = 0; q < nm1; ++q, ++mode_pq)
            {
                simd_type prod = 0.0;
                for (unsigned int r = 0; r < nm2 - p; ++r, ++mode_pqr)
                {
                    prod.fma(in[mode_pqr],
                             basis2[(mode_pr + r) * nq2 + k]); // Load 2x
                }

                fpq[mode_pq] = prod; // Store 1x
            }

            mode_pr += nm2 - p;
        }

        for (unsigned int j = 0; j < nq1; ++j)
        {
            mode_pq = 0;
            for (unsigned int p = 0; p < nm0; ++p)
            {
                simd_type prod = 0.0;
                for (unsigned int q = 0; q < nm1; ++q, ++mode_pq)
                {
                    prod.fma(fpq[mode_pq], basis1[q * nq1 + j]); // Load 2x
                }
                fp[p] = prod; // Store 1x
            }

            for (unsigned int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                simd_type val_kji = 0.0;
                for (unsigned int p = 0; p < nm0; ++p)
                {
                    val_kji.fma(fp[p], basis0[p * nq0 + i]); // Load 2x
                }

                if (isModified)
                {
                    simd_type basis_2 = basis2[nq2 + k]; // Load 1x
                    simd_type basis_0 = basis0[nq0 + i]; // Load 1x

                    for (unsigned int q = 0; q < nm1; ++q)
                    {
                        simd_type coef_0q1 =
                            simd_type(in[q * nm2 + 1]);          // Load 1x
                        simd_type basis_1 = basis1[q * nq1 + j]; // Load 1x
                        val_kji.fma(basis_2 * basis_1, basis_0 * coef_0q1);
                    }
                }
                out[cnt_kji] = val_kji; // Store 1x
            }
        }
    }
}

template <typename simd_type>
NEK_FORCE_INLINE static void BwdTransPyrKernel(
    const unsigned int nm0, const unsigned int nm1, const unsigned int nm2,
    const unsigned int nq0, const unsigned int nq1, const unsigned int nq2,
    const bool isModified, const simd_type *basis0, const simd_type *basis1,
    const simd_type *basis2,
    simd_type *fpq, // nm0 * nm1
    simd_type *fp,  // nm0
    const simd_type *in, simd_type *out)
{
    for (unsigned int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        unsigned int cnt_pqr = 0, mode_pqr = 0, mode_pq = 0;
        for (unsigned int p = 0; p < nm0; ++p)
        {
            for (unsigned int q = 0; q < p; ++q, ++mode_pq)
            {
                simd_type prod = 0.0;
                for (unsigned int r = 0; r < nm2 - p;
                     ++r, ++mode_pqr, ++cnt_pqr)
                {
                    // Load 2x
                    prod.fma(in[cnt_pqr], basis2[mode_pqr * nq2 + k]);
                }
                fpq[mode_pq] = prod; // Store 1x
            }

            for (unsigned int q = p; q < nm1; ++q, ++mode_pq)
            {
                simd_type prod = 0.0;
                for (unsigned int r = 0; r < nm2 - q;
                     ++r, ++mode_pqr, ++cnt_pqr)
                {
                    // Load 2x
                    prod.fma(in[cnt_pqr], basis2[mode_pqr * nq2 + k]);
                }

                fpq[mode_pq] = prod; // Store 1x
            }

            // increment mode in case nm2>nm1
            for (unsigned int q = nm1; q < nm2; ++q)
            {
                mode_pqr += nm2 - q;
            }
        }

        for (unsigned int j = 0; j < nq1; ++j)
        {
            mode_pq = 0;
            for (unsigned int p = 0; p < nm0; ++p)
            {
                simd_type prod = 0.0;
                for (unsigned int q = 0; q < nm1; ++q, ++mode_pq)
                {
                    prod.fma(fpq[mode_pq], basis1[q * nq1 + j]); // Load 2x
                }
                fp[p] = prod; // Store 1x
            }

            for (unsigned int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                simd_type val_kji = 0.0;
                for (unsigned int p = 0; p < nm0; ++p)
                {
                    val_kji.fma(fp[p], basis0[p * nq0 + i]); // Load 2x
                }

                if (isModified)
                {
                    // top vertex
                    //
                    // sum += inarray[1] * base2[nquad2 + k] * (
                    //     base0[i] * base1[nquad1+j] +
                    //     base0[nquad0+i] * base1[j] +
                    //     base0[nquad0+i] * base1[nquad1+j]);
                    simd_type tmp1 = basis0[i] * basis1[nq1 + j]; // Load 2x
                    tmp1.fma(basis0[nq0 + i], basis1[j]);         // Load 2x
                    tmp1.fma(basis0[nq0 + i], basis1[nq1 + j]);   // Load 2x
                    tmp1 = tmp1 * basis2[nq2 + k];                // Load 1x

                    simd_type inarray1 = simd_type(in[1]); // Load 1x
                    val_kji.fma(tmp1, inarray1);
                }

                out[cnt_kji] = val_kji; // Store 1x
            }
        }
    }
}
