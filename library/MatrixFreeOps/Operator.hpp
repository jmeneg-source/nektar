///////////////////////////////////////////////////////////////////////////////
//
// File: Operator.hpp
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

#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "MatrixFreeDeclspec.h"

#include <iostream>

#include <Collections/Operator.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/NekInline.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/SimdLib/tinysimd.hpp>

namespace Nektar::MatrixFree
{

class Operator;

typedef std::shared_ptr<Operator> OperatorSharedPtr;

using vec_t = tinysimd::simd<NekDouble>;

using OperatorFactory =
    LibUtilities::NekFactory<std::string, Operator,
                             std::vector<LibUtilities::BasisSharedPtr>, int>;

MATRIXFREE_EXPORT OperatorFactory &GetOperatorFactory();

/// Helper function, get operator string
MATRIXFREE_EXPORT std::string GetOpstring(LibUtilities::ShapeType shape,
                                          bool deformed = false);

/// Operator base class
class Operator
{
public:
    virtual ~Operator() = default;

    /// This operator requires derivative factors.
    MATRIXFREE_EXPORT virtual bool NeedsDF()
    {
        return false;
    }

    /// This operator requires Jacobian.
    MATRIXFREE_EXPORT virtual bool NeedsJac()
    {
        return false;
    }

    MATRIXFREE_EXPORT virtual void SetUpBdata(
        [[maybe_unused]] std::vector<LibUtilities::BasisSharedPtr> &basis) = 0;
    MATRIXFREE_EXPORT virtual void SetUpDBdata(
        [[maybe_unused]] std::vector<LibUtilities::BasisSharedPtr> &basis) = 0;
    MATRIXFREE_EXPORT virtual void SetUpInterp1D(
        [[maybe_unused]] std::vector<LibUtilities::BasisSharedPtr> &basis,
        [[maybe_unused]] NekDouble factor) = 0;
    MATRIXFREE_EXPORT virtual void SetUpZW(
        [[maybe_unused]] std::vector<LibUtilities::BasisSharedPtr> &basis) = 0;
    MATRIXFREE_EXPORT virtual void SetUpD(
        [[maybe_unused]] std::vector<LibUtilities::BasisSharedPtr> &basis) = 0;

    MATRIXFREE_EXPORT virtual void SetDF(
        const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>
            &df) = 0;

    MATRIXFREE_EXPORT virtual void SetJac(
        const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>
            &jac) = 0;
};

// Base class for backwards transform operator.
class BwdTrans : virtual public Operator
{
public:
    BwdTrans(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : m_basis(basis), m_nElmt(nElmt)
    {
    }

    ~BwdTrans() override = default;

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0; // Abstract Method

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

// Base class for PhysInterp1DScaled operator.
class PhysInterp1DScaled : virtual public Operator
{
public:
    PhysInterp1DScaled(std::vector<LibUtilities::BasisSharedPtr> basis,
                       int nElmt)
        : m_basis(basis), m_nElmt(nElmt)
    {
        // Since the class is constructed before the
        // PhysInterp1DScaled_MatrixFree class inside collections, we are
        // initializing the m_scalingFactor member with the default value of 1.5
        // and upon construction of the aforemmentioned class, it will be
        // updated to the correct value.
        m_scalingFactor = 1.5;
    }

    ~PhysInterp1DScaled() override = default;

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0; // Abstract Method

    // Function to initialize the scaling factor that is used to increase the
    // number of quadrature points in each direction of the flow
    NEK_FORCE_INLINE void SetScalingFactor(NekDouble scalingFactor)
    {
        m_scalingFactor = scalingFactor;
    }

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;

    // Scaling factor to enhance the polynomial space
    NekDouble m_scalingFactor;
};

// Base class for product operator.
class IProduct : virtual public Operator
{
public:
    IProduct(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : m_basis(basis), m_nElmt(nElmt)
    {
    }

    ~IProduct() override = default;

    bool NeedsJac() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0;

protected:
    /// Vector of tensor product basis directions
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

// Base class for product WRT derivative base operator.
class IProductWRTDerivBase : virtual public Operator
{
public:
    IProductWRTDerivBase(std::vector<LibUtilities::BasisSharedPtr> basis,
                         int nElmt)
        : m_basis(basis), m_nElmt(nElmt)
    {
    }

    ~IProductWRTDerivBase() override = default;

    bool NeedsDF() final
    {
        return true;
    }

    bool NeedsJac() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, Array<OneD, NekDouble>> &input,
        Array<OneD, NekDouble> &output) = 0;

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

// Base class for physical derivatives operator.
class PhysDeriv : virtual public Operator
{
public:
    PhysDeriv(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : m_basis(basis), m_nElmt(nElmt)
    {
    }

    ~PhysDeriv() override = default;

    bool NeedsDF() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, Array<OneD, NekDouble>> &output) = 0;

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

// Base class for the Helmholtz base operator.
class Helmholtz : virtual public Operator
{
public:
    Helmholtz(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : m_basis(basis), m_nElmt(nElmt), m_lambda(1.0),
          m_isConstVarDiff(false), m_isVarDiff(false)
    {
        int n          = m_basis.size();
        m_constVarDiff = Array<OneD, NekDouble>(n * (n + 1) / 2);
        int tp         = 1;
        for (int bn = 0; bn < n; ++bn)
        {
            tp *= m_basis[bn]->GetNumPoints();
        }
        switch (n)
        {
            case 2:
                m_varD00 = Array<OneD, NekDouble>(tp);
                m_varD01 = Array<OneD, NekDouble>(tp);
                m_varD11 = Array<OneD, NekDouble>(tp);
                break;
            case 3:
                m_varD00 = Array<OneD, NekDouble>(tp);
                m_varD01 = Array<OneD, NekDouble>(tp);
                m_varD11 = Array<OneD, NekDouble>(tp);
                m_varD02 = Array<OneD, NekDouble>(tp);
                m_varD12 = Array<OneD, NekDouble>(tp);
                m_varD22 = Array<OneD, NekDouble>(tp);
                break;
            default:
                break;
        }
    }

    ~Helmholtz() override = default;

    bool NeedsDF() final
    {
        return true;
    }

    bool NeedsJac() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0;

    NEK_FORCE_INLINE void SetLambda(NekDouble lambda)
    {
        m_lambda = lambda;
    }

    NEK_FORCE_INLINE void SetConstVarDiffusion(Array<OneD, NekDouble> diff)
    {
        m_isConstVarDiff = true;

        int n = m_basis.size();

        for (int i = 0; i < n * (n + 1) / 2; ++i)
        {
            m_constVarDiff[i] = diff[i];
        }
    }

    NEK_FORCE_INLINE void SetVarDiffusion(
        [[maybe_unused]] Array<OneD, NekDouble> diff)
    {
        m_isVarDiff      = true;
        m_isConstVarDiff = false;

        int n  = m_basis.size();
        int tp = 1;

        for (int bn = 0; bn < n; ++bn)
        {
            tp *= m_basis[bn]->GetNumPoints();
        }

        // fixed values for testing!
        for (int i = 0; i < tp; ++i)
        {
            switch (n)
            {
                case 2:
                    m_varD00[i] = diff[0];
                    m_varD01[i] = diff[1];
                    m_varD11[i] = diff[2];
                    break;
                case 3:
                    m_varD00[i] = diff[0];
                    m_varD01[i] = diff[1];
                    m_varD11[i] = diff[2];
                    m_varD02[i] = diff[3];
                    m_varD12[i] = diff[4];
                    m_varD22[i] = diff[5];
                    break;
                default:
                    break;
            }
        }
    }

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
    NekDouble m_lambda;
    bool m_isConstVarDiff;
    Array<OneD, NekDouble> m_constVarDiff;
    bool m_isVarDiff;
    Array<OneD, NekDouble> m_varD00;
    Array<OneD, NekDouble> m_varD01;
    Array<OneD, NekDouble> m_varD11;
    Array<OneD, NekDouble> m_varD02;
    Array<OneD, NekDouble> m_varD12;
    Array<OneD, NekDouble> m_varD22;
};

// Base class for the LinearAdvectionDiffusionReaction base operator.
class LinearAdvectionDiffusionReaction : virtual public Operator
{
public:
    LinearAdvectionDiffusionReaction(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : m_basis(basis), m_nElmt(nElmt), m_lambda(1.0),
          m_isConstVarDiff(false), m_isVarDiff(false)
    {
        int n          = m_basis.size();
        m_constVarDiff = Array<OneD, NekDouble>(n * (n + 1) / 2);
        int tp         = 1;
        for (int bn = 0; bn < n; ++bn)
        {
            tp *= m_basis[bn]->GetNumPoints();
        }

        switch (n)
        {
            case 2:
                m_varD00 = Array<OneD, NekDouble>(tp);
                m_varD01 = Array<OneD, NekDouble>(tp);
                m_varD11 = Array<OneD, NekDouble>(tp);
                break;
            case 3:
                m_varD00 = Array<OneD, NekDouble>(tp);
                m_varD01 = Array<OneD, NekDouble>(tp);
                m_varD11 = Array<OneD, NekDouble>(tp);
                m_varD02 = Array<OneD, NekDouble>(tp);
                m_varD12 = Array<OneD, NekDouble>(tp);
                m_varD22 = Array<OneD, NekDouble>(tp);
                break;
            default:
                break;
        }
    }

    ~LinearAdvectionDiffusionReaction() override = default;

    bool NeedsDF() final
    {
        return true;
    }

    bool NeedsJac() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0;

    NEK_FORCE_INLINE void SetLambda(NekDouble lambda)
    {
        m_lambda = lambda;
    }

    NEK_FORCE_INLINE void SetConstVarDiffusion(Array<OneD, NekDouble> diff)
    {
        m_isConstVarDiff = true;

        int n = m_basis.size();

        for (int i = 0; i < n * (n + 1) / 2; ++i)
        {
            m_constVarDiff[i] = diff[i];
        }
    }

    NEK_FORCE_INLINE void SetVarDiffusion(
        [[maybe_unused]] Array<OneD, NekDouble> diff)
    {
        m_isVarDiff      = true;
        m_isConstVarDiff = false;

        int n  = m_basis.size();
        int tp = 1;

        for (int bn = 0; bn < n; ++bn)
        {
            tp *= m_basis[bn]->GetNumPoints();
        }

        // fixed values for testing!
        for (int i = 0; i < tp; ++i)
        {
            switch (n)
            {
                case 2:
                    m_varD00[i] = diff[0];
                    m_varD01[i] = diff[1];
                    m_varD11[i] = diff[2];
                    break;
                case 3:
                    m_varD00[i] = diff[0];
                    m_varD01[i] = diff[1];
                    m_varD11[i] = diff[2];
                    m_varD02[i] = diff[3];
                    m_varD12[i] = diff[4];
                    m_varD22[i] = diff[5];
                    break;
                default:
                    break;
            }
        }
    }

    NEK_FORCE_INLINE void SetAdvectionVelocities(
        const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>
            &advVel)
    {
        m_advVel = advVel;
    }

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
    NekDouble m_lambda;
    bool m_isConstVarDiff;
    Array<OneD, NekDouble> m_constVarDiff;
    bool m_isVarDiff;
    Array<OneD, NekDouble> m_varD00;
    Array<OneD, NekDouble> m_varD01;
    Array<OneD, NekDouble> m_varD11;
    Array<OneD, NekDouble> m_varD02;
    Array<OneD, NekDouble> m_varD12;
    Array<OneD, NekDouble> m_varD22;
    std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>> m_advVel;
};

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED = false>
class Helper : virtual public Operator
{
    static const size_t DIM{LibUtilities::ShapeTypeDimMap[SHAPE_TYPE]};

protected:
    Helper(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : Operator()
    {
        if (nElmt % vec_t::width == 0) // No padding or already padded
        {
            // Calculate number of 'blocks', i.e. meta-elements
            m_nBlocks = nElmt / vec_t::width;
            m_nPads   = 0;
        }
        else // Need padding internally
        {
            // Calculate number of 'blocks', i.e. meta-elements
            m_nBlocks = nElmt / vec_t::width + 1;
            m_nPads   = vec_t::width - (nElmt % vec_t::width);
        }
        for (int i = 0; i < DIM; ++i)
        {
            m_nm[i] = basis[i]->GetNumModes();
            m_nq[i] = basis[i]->GetNumPoints();
        }
    }

    // Depending on element dimension, set up basis information,
    // inside vectorised environment.
    void SetUpBdata(std::vector<LibUtilities::BasisSharedPtr> &basis) final
    {
        for (int i = 0; i < DIM; ++i)
        {
            const Array<OneD, const NekDouble> bdata = basis[i]->GetBdata();

            m_bdata[i].resize(bdata.size());
            for (auto j = 0; j < bdata.size(); ++j)
            {
                m_bdata[i][j] = bdata[j];
            }
        }
    }

    // Depending on element dimension, set up derivative of basis
    // information, inside vectorised environment.
    void SetUpDBdata(std::vector<LibUtilities::BasisSharedPtr> &basis) final
    {
        for (int i = 0; i < DIM; ++i)
        {
            const Array<OneD, const NekDouble> dbdata = basis[i]->GetDbdata();

            m_dbdata[i].resize(dbdata.size());
            for (auto j = 0; j < dbdata.size(); ++j)
            {
                m_dbdata[i][j] = dbdata[j];
            }
        }
    }

    // Depending on element dimension, set up 1D interpolation matrix in
    // basis data inside vectorised environment.
    void SetUpInterp1D(std::vector<LibUtilities::BasisSharedPtr> &basis,
                       NekDouble factor) final
    {
        // Is this updated at each time-step?
        // Depending on element dimension, set up interpolation matrices,
        // inside vectorised environment.

        int npt0 = basis[0]->GetNumPoints();

        for (int i = 0; i < DIM; ++i)
        {
            int npt         = basis[i]->GetNumPoints();
            m_enhancednq[i] = (npt0 - npt == 1) ? (int)(npt0 * factor - 1)
                                                : (int)(npt * factor);

            LibUtilities::PointsKey PointsKeyIn(m_nq[i],
                                                basis[i]->GetPointsType());
            LibUtilities::PointsKey PointsKeyOut(m_enhancednq[i],
                                                 basis[i]->GetPointsType());
            auto I = LibUtilities::PointsManager()[PointsKeyIn]
                         ->GetI(PointsKeyOut)
                         ->GetPtr();

            m_I[i].resize(I.size());
            for (int j = 0; j < I.size(); ++j)
            {
                m_I[i][j] = I[j];
            }
        }
    }

    void SetUpZW(std::vector<LibUtilities::BasisSharedPtr> &basis) final
    {

        // Depending on element dimension, set up quadrature,
        // inside vectorised environment.
        for (int i = 0; i < DIM; ++i)
        {
            const Array<OneD, const NekDouble> w = basis[i]->GetW();
            NekDouble fac                        = 1.0;
            if (basis[i]->GetPointsType() ==
                LibUtilities::eGaussRadauMAlpha1Beta0)
            {
                fac = 0.5;
            }
            else if (basis[i]->GetPointsType() ==
                     LibUtilities::eGaussRadauMAlpha2Beta0)
            {
                fac = 0.25;
            }

            m_w[i].resize(w.size());
            for (auto j = 0; j < w.size(); ++j)
            {
                m_w[i][j] = fac * w[j];
            }

            auto Z = basis[i]->GetZ();
            m_Z[i].resize(Z.size());
            for (int j = 0; j < Z.size(); ++j)
            {
                m_Z[i][j] = Z[j];
            }
        }
    }

    void SetUpD(std::vector<LibUtilities::BasisSharedPtr> &basis) final
    {

        // Depending on element dimension, set up quadrature,
        // inside vectorised environment.
        for (int i = 0; i < DIM; ++i)
        {
            auto D = basis[i]->GetD()->GetPtr();
            m_D[i].resize(D.size());
            for (int j = 0; j < D.size(); ++j)
            {
                m_D[i][j] = D[j];
            }
        }
    }

    void SetDF(
        const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>
            &df) final
    {
        m_df = df;
    }

    void SetJac(
        const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>
            &jac) final
    {
        m_jac = jac;
    }

    int m_nBlocks;
    int m_nPads;
    std::array<int, DIM> m_nm, m_nq, m_enhancednq;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_bdata;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_dbdata;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM>
        m_D; // Derivatives
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM>
        m_Z; // Zeroes
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM>
        m_w; // Weights
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM>
        m_I; // Interpolation Matrices
    std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>
        m_df; // Chain rule function deriviatives for each element (00, 10,
              // 20, 30...)
    std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>> m_jac;
};
} // namespace Nektar::MatrixFree

#endif
