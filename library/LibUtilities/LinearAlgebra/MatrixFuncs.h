///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixFuncs.h
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
// Description: Matrix functions that depend on storage policy.  Putting
// methods in these separate classes makes it easier to use them from
// normal, scaled, and block matrices.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_FUNCS_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_FUNCS_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>

namespace Nektar
{
struct LIB_UTILITIES_EXPORT BandedMatrixFuncs
{
    /// \brief Calculates and returns the storage size required.
    ///
    /// This method assumes that the matrix will be used with LU factorizationa
    /// and allocates additional storage as appropriate.
    static unsigned int GetRequiredStorageSize(unsigned int totalRows,
                                               unsigned int totalColumns,
                                               unsigned int subDiags,
                                               unsigned int superDiags);

    static unsigned int CalculateNumberOfDiags(unsigned int totalRows,
                                               unsigned int diags);

    static unsigned int CalculateNumberOfRows(unsigned int totalRows,
                                              unsigned int subDiags,
                                              unsigned int superDiags);

    static unsigned int CalculateIndex(unsigned int totalRows,
                                       unsigned int totalColumns,
                                       unsigned int row, unsigned int column,
                                       unsigned int sub, unsigned int super);

    static std::tuple<unsigned int, unsigned int> Advance(
        const unsigned int totalRows, const unsigned int totalColumns,
        const unsigned int curRow, const unsigned int curColumn);
};

struct LIB_UTILITIES_EXPORT FullMatrixFuncs
{

    static unsigned int GetRequiredStorageSize(unsigned int rows,
                                               unsigned int columns);
    static unsigned int CalculateIndex(unsigned int totalRows,
                                       unsigned int totalColumns,
                                       unsigned int curRow,
                                       unsigned int curColumn);

    static std::tuple<unsigned int, unsigned int> Advance(
        const unsigned int totalRows, const unsigned int totalColumns,
        const unsigned int curRow, const unsigned int curColumn);

    template <typename DataType>
    static void Invert(unsigned int rows, unsigned int columns,
                       Array<OneD, DataType> &data, const char transpose)
    {
        ASSERTL0(rows == columns, "Only square matrices can be inverted.");
        ASSERTL0(transpose == 'N',
                 "Only untransposed matrices may be inverted.");

        int m    = rows;
        int n    = columns;
        int info = 0;
        Array<OneD, int> ipivot(n);
        Array<OneD, DataType> work(n);

        if (std::is_floating_point_v<DataType>)
        {
            switch (sizeof(DataType))
            {
                case sizeof(NekDouble):
                    break;
                case sizeof(NekSingle):
                    break;
                default:
                    ASSERTL0(
                        false,
                        "Invert DataType is neither NekDouble nor NekSingle");
                    break;
            }
        }
        else
        {
            ASSERTL0(false,
                     "FullMatrixFuncs::Invert DataType is not floating point");
        }

        Lapack::DoSgetrf(m, n, data.data(), m, ipivot.data(), info);

        if (info < 0)
        {
            std::string message =
                "ERROR: The " + std::to_string(-info) +
                "th parameter had an illegal parameter for dgetrf";
            ASSERTL0(false, message.c_str());
        }
        else if (info > 0)
        {
            std::string message = "ERROR: Element u_" + std::to_string(info) +
                                  std::to_string(info) + " is 0 from dgetrf";
            ASSERTL0(false, message.c_str());
        }

        Lapack::DoSgetri(n, data.data(), n, ipivot.data(), work.data(), n,
                         info);

        if (info < 0)
        {
            std::string message =
                "ERROR: The " + std::to_string(-info) +
                "th parameter had an illegal parameter for dgetri";
            ASSERTL0(false, message.c_str());
        }
        else if (info > 0)
        {
            std::string message = "ERROR: Element u_" + std::to_string(info) +
                                  std::to_string(info) + " is 0 from dgetri";
            ASSERTL0(false, message.c_str());
        }
    }

    // add a method to do the pseudo inverse of a full matrix, using SVD routine
    // from Lapack.
    template <typename DataType>
    static void PseudoInverse(unsigned int rows, unsigned int columns,
                              Array<OneD, DataType> &data)
    {
        int m      = rows;
        int n      = columns;
        int lda    = m;
        int ldu    = m;
        int ldvt   = n;
        int info   = 0;
        int min_mn = std::min(m, n);
        // Copy input data since SVD destroys it
        Array<OneD, DataType> a(data);
        // SVD outputs
        Array<OneD, DataType> s(min_mn, 0.0);
        Array<OneD, DataType> u(ldu * m, 0.0);
        Array<OneD, DataType> vt(ldvt * n, 0.0);
        // Workspace query
        DataType wkopt = 0.0;
        int lwork      = -1;
        // SVD job: all singular vectors
        char jobu  = 'A';
        char jobvt = 'A';
        //
        // DGESVD computes the singular value decomposition (SVD) of a real
        // M-by-N matrix A, optionally computing the left and/or right singular
        // vectors. The SVD is written
        //
        //      A = U * SIGMA * transpose(V)
        //
        // where SIGMA is an M-by-N matrix which is zero except for its
        // min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
        // V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
        // are the singular values of A; they are real and non-negative, and
        // are returned in descending order.  The first min(m,n) columns of
        // U and V are the left and right singular vectors of A.
        //
        // Note that the routine returns V**T, not V.
        //
        if (sizeof(DataType) == sizeof(NekDouble))
        {
            Lapack::Dgesvd(jobu, jobvt, m, n, (double *)a.data(), lda,
                           (double *)s.data(), (double *)u.data(), ldu,
                           (double *)vt.data(), ldvt, (double *)&wkopt, lwork,
                           info);
        }
        else if (sizeof(DataType) == sizeof(NekSingle))
        {
            Lapack::Sgesvd(jobu, jobvt, m, n, (float *)a.data(), lda,
                           (float *)s.data(), (float *)u.data(), ldu,
                           (float *)vt.data(), ldvt, (float *)&wkopt, lwork,
                           info);
        }
        else
        {
            ASSERTL0(
                false,
                "PseudoInverse DataType is neither NekDouble nor NekSingle");
        }
        lwork = static_cast<int>(wkopt);
        Array<OneD, DataType> work(lwork);
        if (sizeof(DataType) == sizeof(NekDouble))
        {
            Lapack::Dgesvd(jobu, jobvt, m, n, (double *)a.data(), lda,
                           (double *)s.data(), (double *)u.data(), ldu,
                           (double *)vt.data(), ldvt, (double *)work.data(),
                           lwork, info);
        }
        else if (sizeof(DataType) == sizeof(NekSingle))
        {
            Lapack::Sgesvd(jobu, jobvt, m, n, (float *)a.data(), lda,
                           (float *)s.data(), (float *)u.data(), ldu,
                           (float *)vt.data(), ldvt, (float *)work.data(),
                           lwork, info);
        }
        if (info < 0)
        {
            std::string message = "ERROR: The " + std::to_string(-info) +
                                  "th parameter had an illegal value for gesvd";
            ASSERTL0(false, message.c_str());
        }
        else if (info > 0)
        {
            std::string message = "ERROR: SVD did not converge in gesvd";
            ASSERTL0(false, message.c_str());
        }
        // Compute pseudo-inverse: A^+ = V * S^+ * U^T
        // S^+ is diagonal with 1/s_i for s_i > tol, 0 otherwise
        DataType tol = std::numeric_limits<DataType>::epsilon() *
                       std::max(m, n) * std::abs(s[0]);
        Array<OneD, DataType> sinv(min_mn, 0.0);
        for (int i = 0; i < min_mn; ++i)
        {
            if (std::abs(s[i]) > tol)
            {
                sinv[i] = 1.0 / s[i];
            }
            else
            {
                sinv[i] = 0.0;
            }
        }
        // temp = S^+ * U^T (size min_mn x m)
        Array<OneD, DataType> temp(min_mn * m, 0.0);
        // scale each row of U^T by sinv
        for (int i = 0; i < min_mn; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                temp[i + j * min_mn] = sinv[i] * u[j + i * ldu];
            }
        }
        // Use BLAS GEMM: data = V (n x min_mn) * temp (min_mn x m)
        data           = Array<OneD, DataType>(n * m, 0.0);
        char transa    = 'T';
        char transb    = 'N';
        DataType alpha = 1.0;
        DataType beta  = 0.0;
        if (sizeof(DataType) == sizeof(NekDouble))
        {
            Blas::Gemm(transa, transb, n, m, min_mn, alpha, (double *)vt.data(),
                       ldvt, (double *)temp.data(), min_mn, beta,
                       (double *)data.data(), n);
        }
        else if (sizeof(DataType) == sizeof(NekSingle))
        {
            Blas::Gemm(transa, transb, n, m, min_mn, (float)alpha,
                       (float *)vt.data(), ldvt, (float *)temp.data(), min_mn,
                       (float)beta, (float *)data.data(), n);
        }
    }

    template <typename DataType>
    static void EigenSolve(unsigned int n, const Array<OneD, const DataType> &A,
                           Array<OneD, DataType> &EigValReal,
                           Array<OneD, DataType> &EigValImag,
                           Array<OneD, DataType> &EigVecs)
    {
        int lda = n, info = 0;
        DataType dum;
        char uplo = 'N';

        if (EigVecs.size() != 0) // calculate Right Eigen Vectors
        {
            int lwork = 4 * lda;
            Array<OneD, DataType> work(4 * lda);
            char lrev = 'V';
            Lapack::DoSgeev(uplo, lrev, lda, A.data(), lda, EigValReal.data(),
                            EigValImag.data(), &dum, 1, EigVecs.data(), lda,
                            &work[0], lwork, info);
        }
        else
        {
            int lwork = 3 * lda;
            Array<OneD, DataType> work(3 * lda);
            char lrev = 'N';
            Lapack::DoSgeev(uplo, lrev, lda, A.data(), lda, EigValReal.data(),
                            EigValImag.data(), &dum, 1, &dum, 1, &work[0],
                            lwork, info);
        }
        ASSERTL0(info == 0, "Info is not zero");
    }
};

struct LIB_UTILITIES_EXPORT TriangularMatrixFuncs
{
    static unsigned int GetRequiredStorageSize(unsigned int rows,
                                               unsigned int columns);
};

struct LIB_UTILITIES_EXPORT UpperTriangularMatrixFuncs
    : public TriangularMatrixFuncs
{
    static unsigned int CalculateIndex(unsigned int curRow,
                                       unsigned int curColumn);

    static std::tuple<unsigned int, unsigned int> Advance(
        const unsigned int totalRows, const unsigned int totalColumns,
        const unsigned int curRow, const unsigned int curColumn);
};

struct LIB_UTILITIES_EXPORT LowerTriangularMatrixFuncs
    : public TriangularMatrixFuncs
{
    static unsigned int CalculateIndex(unsigned int totalColumns,
                                       unsigned int curRow,
                                       unsigned int curColumn);

    static std::tuple<unsigned int, unsigned int> Advance(
        const unsigned int totalRows, const unsigned int totalColumns,
        const unsigned int curRow, const unsigned int curColumn,
        char transpose = 'N');
};

/// \internal
/// Symmetric matrices use upper triangular packed storage.
struct LIB_UTILITIES_EXPORT SymmetricMatrixFuncs : private TriangularMatrixFuncs
{
    using TriangularMatrixFuncs::GetRequiredStorageSize;

    static unsigned int CalculateIndex(unsigned int curRow,
                                       unsigned int curColumn);

    template <typename DataType>
    static void Invert(unsigned int rows, unsigned int columns,
                       Array<OneD, DataType> &data)
    {
        ASSERTL0(rows == columns, "Only square matrices can be inverted.");

        int n    = columns;
        int info = 0;
        Array<OneD, int> ipivot(n);
        Array<OneD, DataType> work(n);

        Lapack::DoSsptrf('U', n, data.data(), ipivot.data(), info);

        if (info < 0)
        {
            std::string message =
                "ERROR: The " + std::to_string(-info) +
                "th parameter had an illegal parameter for dsptrf";
            ASSERTL0(false, message.c_str());
        }
        else if (info > 0)
        {
            std::string message = "ERROR: Element u_" + std::to_string(info) +
                                  std::to_string(info) + " is 0 from dsptrf";
            ASSERTL0(false, message.c_str());
        }

        Lapack::DoSsptri('U', n, data.data(), ipivot.data(), work.data(), info);

        if (info < 0)
        {
            std::string message =
                "ERROR: The " + std::to_string(-info) +
                "th parameter had an illegal parameter for dsptri";
            ASSERTL0(false, message.c_str());
        }
        else if (info > 0)
        {
            std::string message = "ERROR: Element u_" + std::to_string(info) +
                                  std::to_string(info) + " is 0 from dsptri";
            ASSERTL0(false, message.c_str());
        }
    }

    static std::tuple<unsigned int, unsigned int> Advance(
        const unsigned int totalRows, const unsigned int totalColumns,
        const unsigned int curRow, const unsigned int curColumn);
};

struct LIB_UTILITIES_EXPORT DiagonalMatrixFuncs
{
    static std::tuple<unsigned int, unsigned int> Advance(
        const unsigned int totalRows, const unsigned int totalColumns,
        const unsigned int curRow, const unsigned int curColumn);

    template <typename DataType>
    static void Invert(unsigned int rows, unsigned int columns,
                       Array<OneD, DataType> &data)
    {
        ASSERTL0(rows == columns, "Only square matrices can be inverted.");
        for (unsigned int i = 0; i < rows; ++i)
        {
            data[i] = 1.0 / data[i];
        }
    }

    static unsigned int GetRequiredStorageSize(unsigned int rows,
                                               unsigned int columns);

    static unsigned int CalculateIndex(unsigned int row, unsigned int col);
};

struct LIB_UTILITIES_EXPORT TriangularBandedMatrixFuncs
{
    static unsigned int GetRequiredStorageSize(unsigned int rows,
                                               unsigned int columns,
                                               unsigned int nSubSuperDiags);
};

struct LIB_UTILITIES_EXPORT UpperTriangularBandedMatrixFuncs
    : public TriangularBandedMatrixFuncs
{
};

struct LIB_UTILITIES_EXPORT LowerTriangularBandedMatrixFuncs
    : public TriangularBandedMatrixFuncs
{
};

/// \internal
/// Symmetric banded matrices use upper triangular banded packed storage.
struct LIB_UTILITIES_EXPORT SymmetricBandedMatrixFuncs
    : private TriangularBandedMatrixFuncs
{
    using TriangularBandedMatrixFuncs::GetRequiredStorageSize;

    static unsigned int CalculateIndex(unsigned int curRow,
                                       unsigned int curColumn,
                                       unsigned int nSuperDiags);
};

} // namespace Nektar

#endif // NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_FUNCS_H
