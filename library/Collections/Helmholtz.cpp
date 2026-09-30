///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz.cpp
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
// Description: Helmholtz operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Collection.h>
#include <Collections/IProduct.h>
#include <Collections/MatrixFreeBase.h>
#include <Collections/Operator.h>
#include <MatrixFreeOps/Operator.hpp>

using namespace std;

namespace Nektar::Collections
{

using LibUtilities::eHexahedron;
using LibUtilities::eNodalPrism;
using LibUtilities::eNodalTet;
using LibUtilities::eNodalTri;
using LibUtilities::ePrism;
using LibUtilities::ePyramid;
using LibUtilities::eQuadrilateral;
using LibUtilities::eSegment;
using LibUtilities::eTetrahedron;
using LibUtilities::eTriangle;

/**
 * @brief Helmholtz help class to calculate the size of the collection
 * that is given as an input and as an output to the Helmholtz Operator. The
 * size evaluation takes into account that the evaluation of the Helmholtz
 * operator takes input from the coeff space and gives the output in the coeff
 * space too.
 */
class Helmholtz_Helper : virtual public Operator
{
protected:
    Helmholtz_Helper()
    {
        // expect input to be number of elements by the number of coefficients
        m_inputSize = m_numElmt * m_stdExp->GetNcoeffs();

        // expect output to be number of elements by the number of coefficients
        // computation is from coeff space to coeff space
        m_outputSize = m_inputSize;
    }
};

/**
 * @brief Helmholtz operator using LocalRegions implementation.
 */
class Helmholtz_NoCollection final : virtual public Operator, Helmholtz_Helper
{
public:
    OPERATOR_CREATE(Helmholtz_NoCollection)

    ~Helmholtz_NoCollection() final = default;

    void operator()(const Array<OneD, const NekDouble> &entry0,
                    Array<OneD, NekDouble> &entry1,
                    [[maybe_unused]] Array<OneD, NekDouble> &entry2,
                    [[maybe_unused]] Array<OneD, NekDouble> &entry3,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        unsigned int nmodes = m_expList[0]->GetNcoeffs();
        unsigned int nphys  = m_expList[0]->GetTotPoints();
        Array<OneD, NekDouble> tmp;

        for (int n = 0; n < m_numElmt; ++n)
        {
            // Restrict varcoeffs to size of element
            StdRegions::VarCoeffMap varcoeffs = StdRegions::NullVarCoeffMap;
            if (m_varcoeffs.size())
            {
                varcoeffs =
                    StdRegions::RestrictCoeffMap(m_varcoeffs, n * nphys, nphys);
            }

            StdRegions::StdMatrixKey mkey(
                StdRegions::eHelmholtz, (m_expList)[n]->DetShapeType(),
                *(m_expList)[n], m_factors, varcoeffs);
            m_expList[n]->GeneralMatrixOp(entry0 + n * nmodes,
                                          tmp = entry1 + n * nmodes, mkey);
        }
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

    void UpdateFactors(StdRegions::FactorMap factors) override
    {
        m_factors = factors;
    }

    void UpdateVarcoeffs(StdRegions::VarCoeffMap &varcoeffs) override
    {
        m_varcoeffs = varcoeffs;
    }

protected:
    int m_dim;
    int m_coordim;
    vector<LocalRegions::ExpansionSharedPtr> m_expList;
    StdRegions::FactorMap m_factors;
    StdRegions::VarCoeffMap m_varcoeffs;

private:
    Helmholtz_NoCollection(vector<LocalRegions::ExpansionSharedPtr> pCollExp,
                           CoalescedGeomDataSharedPtr pGeomData,
                           StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), Helmholtz_Helper()
    {
        m_expList = pCollExp;
        m_dim     = pCollExp[0]->GetNumBases();
        m_coordim = pCollExp[0]->GetCoordim();

        m_factors   = factors;
        m_varcoeffs = StdRegions::NullVarCoeffMap;
    }
};

/// Factory initialisation for the Helmholtz_NoCollection operators
OperatorKey Helmholtz_NoCollection::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment, eHelmholtz, eNoCollection, false),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eHelmholtz, eNoCollection, false),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eNodalTri, eHelmholtz, eNoCollection, true),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eNoCollection, false),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eHelmholtz, eNoCollection, false),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eNodalTet, eHelmholtz, eNoCollection, true),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eHelmholtz, eNoCollection, false),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eHelmholtz, eNoCollection, false),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eNodalPrism, eHelmholtz, eNoCollection, true),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eHelmholtz, eNoCollection, false),
        Helmholtz_NoCollection::create, "Helmholtz_NoCollection_Hex")};

/**
 * @brief Helmholtz operator using LocalRegions implementation.
 */
class Helmholtz_IterPerExp final : virtual public Operator, Helmholtz_Helper
{
public:
    OPERATOR_CREATE(Helmholtz_IterPerExp)

    ~Helmholtz_IterPerExp() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        const int nCoeffs = m_stdExp->GetNcoeffs();
        const int nPhys   = m_stdExp->GetTotPoints();

        ASSERTL1(input.size() >= m_numElmt * nCoeffs,
                 "input array size is insufficient");
        ASSERTL1(output.size() >= m_numElmt * nCoeffs,
                 "output array size is insufficient");

        Array<OneD, NekDouble> tmpphys, t1;
        Array<OneD, Array<OneD, NekDouble>> dtmp(3);
        Array<OneD, Array<OneD, NekDouble>> tmp(3);

        tmpphys = wsp;
        for (int i = 1; i < m_coordim + 1; ++i)
        {
            dtmp[i - 1] = wsp + i * nPhys;
            tmp[i - 1]  = wsp + (i + m_coordim) * nPhys;
        }

        for (int i = 0; i < m_numElmt; ++i)
        {
            m_stdExp->BwdTrans(input + i * nCoeffs, tmpphys);

            // local derivative
            m_stdExp->PhysDeriv(tmpphys, dtmp[0], dtmp[1], dtmp[2]);

            // determine mass matrix term
            if (m_isDeformed)
            {
                Vmath::Vmul(nPhys, m_jac + i * nPhys, 1, tmpphys, 1, tmpphys,
                            1);
            }
            else
            {
                Vmath::Smul(nPhys, m_jac[i], tmpphys, 1, tmpphys, 1);
            }

            m_stdExp->IProductWRTBase(tmpphys, t1 = output + i * nCoeffs);
            Vmath::Smul(nCoeffs, m_lambda, output + i * nCoeffs, 1,
                        t1 = output + i * nCoeffs, 1);

            if (m_isDeformed)
            {
                // calculate full derivative
                for (int j = 0; j < m_coordim; ++j)
                {
                    Vmath::Vmul(nPhys,
                                m_derivFac[j * m_dim].origin() + i * nPhys, 1,
                                &dtmp[0][0], 1, &tmp[j][0], 1);

                    for (int k = 1; k < m_dim; ++k)
                    {
                        Vmath::Vvtvp(
                            nPhys,
                            m_derivFac[j * m_dim + k].origin() + i * nPhys, 1,
                            &dtmp[k][0], 1, &tmp[j][0], 1, &tmp[j][0], 1);
                    }
                }

                if (m_HasVarCoeffDiff)
                {
                    // calculate dtmp[i] = dx/dxi sum_j diff[0][j] tmp[j]
                    //                   + dy/dxi sum_j diff[1][j] tmp[j]
                    //                   + dz/dxi sum_j diff[2][j] tmp[j]

                    // First term
                    Vmath::Smul(nPhys, m_diff[0][0], &tmp[0][0], 1, &tmpphys[0],
                                1);
                    for (int l = 1; l < m_coordim; ++l)
                    {
                        Vmath::Svtvp(nPhys, m_diff[0][l], &tmp[l][0], 1,
                                     &tmpphys[0], 1, &tmpphys[0], 1);
                    }

                    for (int j = 0; j < m_dim; ++j)
                    {
                        Vmath::Vmul(nPhys, m_derivFac[j].origin() + i * nPhys,
                                    1, &tmpphys[0], 1, &dtmp[j][0], 1);
                    }

                    // Second and third terms
                    for (int k = 1; k < m_coordim; ++k)
                    {
                        Vmath::Smul(nPhys, m_diff[k][0], &tmp[0][0], 1,
                                    &tmpphys[0], 1);
                        for (int l = 1; l < m_coordim; ++l)
                        {
                            Vmath::Svtvp(nPhys, m_diff[k][l], &tmp[l][0], 1,
                                         &tmpphys[0], 1, &tmpphys[0], 1);
                        }

                        for (int j = 0; j < m_dim; ++j)
                        {
                            Vmath::Vvtvp(nPhys,
                                         m_derivFac[j + k * m_dim].origin() +
                                             i * nPhys,
                                         1, &tmpphys[0], 1, &dtmp[j][0], 1,
                                         &dtmp[j][0], 1);
                        }
                    }
                }
                else
                {
                    // calculate dx/dxi tmp[0] + dy/dxi tmp[1]
                    //                         + dz/dxi tmp[2]
                    for (int j = 0; j < m_dim; ++j)
                    {
                        Vmath::Vmul(nPhys, m_derivFac[j].origin() + i * nPhys,
                                    1, &tmp[0][0], 1, &dtmp[j][0], 1);

                        for (int k = 1; k < m_coordim; ++k)
                        {
                            Vmath::Vvtvp(nPhys,
                                         m_derivFac[j + k * m_dim].origin() +
                                             i * nPhys,
                                         1, &tmp[k][0], 1, &dtmp[j][0], 1,
                                         &dtmp[j][0], 1);
                        }
                    }
                }

                // calculate Iproduct WRT Std Deriv
                for (int j = 0; j < m_dim; ++j)
                {
                    // multiply by Jacobian
                    Vmath::Vmul(nPhys, m_jac + i * nPhys, 1, dtmp[j], 1,
                                dtmp[j], 1);

                    m_stdExp->IProductWRTDerivBase(j, dtmp[j], tmp[0]);
                    Vmath::Vadd(nCoeffs, tmp[0], 1, output + i * nCoeffs, 1,
                                t1 = output + i * nCoeffs, 1);
                }
            }
            else
            {
                // calculate full derivative
                for (int j = 0; j < m_coordim; ++j)
                {
                    Vmath::Smul(nPhys, m_derivFac[j * m_dim][i], &dtmp[0][0], 1,
                                &tmp[j][0], 1);

                    for (int k = 1; k < m_dim; ++k)
                    {
                        Vmath::Svtvp(nPhys, m_derivFac[j * m_dim + k][i],
                                     &dtmp[k][0], 1, &tmp[j][0], 1, &tmp[j][0],
                                     1);
                    }
                }

                if (m_HasVarCoeffDiff)
                {
                    // calculate dtmp[i] = dx/dxi sum_j diff[0][j] tmp[j]
                    //                   + dy/dxi sum_j diff[1][j] tmp[j]
                    //                   + dz/dxi sum_j diff[2][j] tmp[j]

                    // First term
                    Vmath::Smul(nPhys, m_diff[0][0], &tmp[0][0], 1, &tmpphys[0],
                                1);
                    for (int l = 1; l < m_coordim; ++l)
                    {
                        Vmath::Svtvp(nPhys, m_diff[0][l], &tmp[l][0], 1,
                                     &tmpphys[0], 1, &tmpphys[0], 1);
                    }

                    for (int j = 0; j < m_dim; ++j)
                    {
                        Vmath::Smul(nPhys, m_derivFac[j][i], &tmpphys[0], 1,
                                    &dtmp[j][0], 1);
                    }

                    // Second and third terms
                    for (int k = 1; k < m_coordim; ++k)
                    {
                        Vmath::Smul(nPhys, m_diff[k][0], &tmp[0][0], 1,
                                    &tmpphys[0], 1);
                        for (int l = 1; l < m_coordim; ++l)
                        {
                            Vmath::Svtvp(nPhys, m_diff[k][l], &tmp[l][0], 1,
                                         &tmpphys[0], 1, &tmpphys[0], 1);
                        }

                        for (int j = 0; j < m_dim; ++j)
                        {
                            Vmath::Svtvp(nPhys, m_derivFac[j + k * m_dim][i],
                                         &tmpphys[0], 1, &dtmp[j][0], 1,
                                         &dtmp[j][0], 1);
                        }
                    }
                }
                else
                {
                    // calculate dx/dxi tmp[0] + dy/dxi tmp[2]
                    //                         + dz/dxi tmp[3]
                    for (int j = 0; j < m_dim; ++j)
                    {
                        Vmath::Smul(nPhys, m_derivFac[j][i], &tmp[0][0], 1,
                                    &dtmp[j][0], 1);

                        for (int k = 1; k < m_coordim; ++k)
                        {
                            Vmath::Svtvp(nPhys, m_derivFac[j + k * m_dim][i],
                                         &tmp[k][0], 1, &dtmp[j][0], 1,
                                         &dtmp[j][0], 1);
                        }
                    }
                }

                // calculate Iproduct WRT Std Deriv
                for (int j = 0; j < m_dim; ++j)
                {
                    // multiply by Jacobian
                    Vmath::Smul(nPhys, m_jac[i], dtmp[j], 1, dtmp[j], 1);

                    m_stdExp->IProductWRTDerivBase(j, dtmp[j], tmp[0]);
                    Vmath::Vadd(nCoeffs, tmp[0], 1, output + i * nCoeffs, 1,
                                t1 = output + i * nCoeffs, 1);
                }
            }
        }
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

    /**
     * @brief Check the validity of supplied constant factors.
     *
     * @param factors Map of factors
     */
    void UpdateFactors(StdRegions::FactorMap factors) override
    {
        // If match previous factors, nothing to do.
        if (m_factors == factors)
        {
            return;
        }

        m_factors = factors;

        // Check Lambda constant of Helmholtz operator
        auto x = factors.find(StdRegions::eFactorLambda);
        ASSERTL1(
            x != factors.end(),
            "Constant factor not defined: " +
                std::string(
                    StdRegions::ConstFactorTypeMap[StdRegions::eFactorLambda]));
        m_lambda = x->second;

        // If varcoeffs not supplied, nothing else to do.
        m_HasVarCoeffDiff = false;
        auto d            = factors.find(StdRegions::eFactorCoeffD00);
        if (d == factors.end())
        {
            return;
        }

        m_diff = Array<OneD, Array<OneD, NekDouble>>(m_coordim);
        for (int i = 0; i < m_coordim; ++i)
        {
            m_diff[i] = Array<OneD, NekDouble>(m_coordim, 0.0);
        }

        for (int i = 0; i < m_coordim; ++i)
        {
            d = factors.find(m_factorCoeffDef[i][i]);
            if (d != factors.end())
            {
                m_diff[i][i] = d->second;
            }
            else
            {
                m_diff[i][i] = 1.0;
            }

            for (int j = i + 1; j < m_coordim; ++j)
            {
                d = factors.find(m_factorCoeffDef[i][j]);
                if (d != factors.end())
                {
                    m_diff[i][j] = m_diff[j][i] = d->second;
                }
            }
        }
        m_HasVarCoeffDiff = true;
    }

    /**
     * @brief Check the validity of supplied variable coefficients.
     * Note that the m_varcoeffs are not implemented for the IterPerExp
     * operator.
     * There exists a specialised implementation for variable diffusion in the
     * above function UpdateFactors.
     *
     * @param varcoeffs Map of varcoeffs
     */
    void UpdateVarcoeffs(StdRegions::VarCoeffMap &varcoeffs) override
    {
        m_varcoeffs = varcoeffs;
    }

protected:
    Array<TwoD, const NekDouble> m_derivFac;
    Array<OneD, const NekDouble> m_jac;
    int m_dim;
    int m_coordim;
    StdRegions::FactorMap m_factors;
    StdRegions::VarCoeffMap m_varcoeffs;
    NekDouble m_lambda;
    bool m_HasVarCoeffDiff;
    Array<OneD, Array<OneD, NekDouble>> m_diff;
    const StdRegions::ConstFactorType m_factorCoeffDef[3][3] = {
        {StdRegions::eFactorCoeffD00, StdRegions::eFactorCoeffD01,
         StdRegions::eFactorCoeffD02},
        {StdRegions::eFactorCoeffD01, StdRegions::eFactorCoeffD11,
         StdRegions::eFactorCoeffD12},
        {StdRegions::eFactorCoeffD02, StdRegions::eFactorCoeffD12,
         StdRegions::eFactorCoeffD22}};

private:
    Helmholtz_IterPerExp(vector<LocalRegions::ExpansionSharedPtr> pCollExp,
                         CoalescedGeomDataSharedPtr pGeomData,
                         StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), Helmholtz_Helper()
    {
        m_dim     = pCollExp[0]->GetShapeDimension();
        m_coordim = pCollExp[0]->GetCoordim();
        int nqtot = m_stdExp->GetTotPoints();

        m_derivFac = pGeomData->GetDerivFactors(pCollExp);
        m_jac      = pGeomData->GetJac(pCollExp);
        m_wspSize  = (2 * m_coordim + 1) * nqtot;

        m_lambda          = 1.0;
        m_HasVarCoeffDiff = false;
        m_factors         = StdRegions::NullFactorMap;
        m_varcoeffs       = StdRegions::NullVarCoeffMap;
        this->UpdateFactors(factors);
    }
};

/// Factory initialisation for the Helmholtz_IterPerExp operators
OperatorKey Helmholtz_IterPerExp::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment, eHelmholtz, eIterPerExp, false),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eHelmholtz, eIterPerExp, false),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eNodalTri, eHelmholtz, eIterPerExp, true),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eIterPerExp, false),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eHelmholtz, eIterPerExp, false),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eNodalTet, eHelmholtz, eIterPerExp, true),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eHelmholtz, eIterPerExp, false),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eHelmholtz, eIterPerExp, false),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eNodalPrism, eHelmholtz, eIterPerExp, true),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eHelmholtz, eIterPerExp, false),
        Helmholtz_IterPerExp::create, "Helmholtz_IterPerExp_Hex")};

/**
 * @brief Helmholtz operator using matrix free operators.
 */
class Helmholtz_MatrixFree final : virtual public Operator,
                                   MatrixFreeBase,
                                   Helmholtz_Helper
{
public:
    OPERATOR_CREATE(Helmholtz_MatrixFree)

    ~Helmholtz_MatrixFree() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output0,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        (*m_oper)(input, output0);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

    /**
     *
     */
    void UpdateFactors(StdRegions::FactorMap factors) override
    {
        if (factors == m_factors)
        {
            return;
        }

        m_factors = factors;

        // Set lambda for this call
        auto x = factors.find(StdRegions::eFactorLambda);
        ASSERTL1(
            x != factors.end(),
            "Constant factor not defined: " +
                std::string(
                    StdRegions::ConstFactorTypeMap[StdRegions::eFactorLambda]));
        m_oper->SetLambda(x->second);

        // set constant diffusion coefficients
        bool isConstVarDiff         = false;
        Array<OneD, NekDouble> diff = Array<OneD, NekDouble>(6, 0.0);
        diff[0] = diff[2] = diff[5] = 1.0;

        auto xd00 = factors.find(StdRegions::eFactorCoeffD00);
        if (xd00 != factors.end() && xd00->second != 1.0)
        {
            isConstVarDiff = true;
            diff[0]        = xd00->second;
        }

        auto xd01 = factors.find(StdRegions::eFactorCoeffD01);
        if (xd01 != factors.end() && xd01->second != 0.0)
        {
            isConstVarDiff = true;
            diff[1]        = xd01->second;
        }

        auto xd11 = factors.find(StdRegions::eFactorCoeffD11);
        if (xd11 != factors.end() && xd11->second != 1.0)
        {
            isConstVarDiff = true;
            diff[2]        = xd11->second;
        }

        auto xd02 = factors.find(StdRegions::eFactorCoeffD02);
        if (xd02 != factors.end() && xd02->second != 0.0)
        {
            isConstVarDiff = true;
            diff[3]        = xd02->second;
        }

        auto xd12 = factors.find(StdRegions::eFactorCoeffD12);
        if (xd12 != factors.end() && xd12->second != 0.0)
        {
            isConstVarDiff = true;
            diff[4]        = xd12->second;
        }

        auto xd22 = factors.find(StdRegions::eFactorCoeffD22);
        if (xd22 != factors.end() && xd22->second != 1.0)
        {
            isConstVarDiff = true;
            diff[5]        = xd22->second;
        }

        if (isConstVarDiff)
        {
            m_oper->SetConstVarDiffusion(diff);
        }

        // set random here for fn variable diffusion
        auto k = factors.find(StdRegions::eFactorTau);
        if (k != factors.end() && k->second != 0.0)
        {
            m_oper->SetVarDiffusion(diff);
        }
    }

    /**
     * @brief Check the validity of supplied variable coefficients.
     * Note that the m_varcoeffs are not implemented for the MatrixFree
     * operator.
     * There exists a specialised implementation for variable diffusion in the
     * above function UpdateFactors.
     *
     * @param varcoeffs Map of varcoeffs
     */
    void UpdateVarcoeffs(StdRegions::VarCoeffMap &varcoeffs) override
    {
        m_varcoeffs = varcoeffs;
    }

private:
    std::shared_ptr<MatrixFree::Helmholtz> m_oper;
    unsigned int m_nmtot;
    StdRegions::FactorMap m_factors;
    StdRegions::VarCoeffMap m_varcoeffs;

    Helmholtz_MatrixFree(vector<LocalRegions::ExpansionSharedPtr> pCollExp,
                         CoalescedGeomDataSharedPtr pGeomData,
                         StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors),
          MatrixFreeBase(pCollExp[0]->GetNcoeffs(), pCollExp[0]->GetNcoeffs(),
                         pCollExp.size()),
          Helmholtz_Helper()
    {

        m_nmtot = m_numElmt * pCollExp[0]->GetNcoeffs();

        const auto dim = pCollExp[0]->GetShapeDimension();

        // Basis vector.
        std::vector<LibUtilities::BasisSharedPtr> basis(dim);
        for (auto i = 0; i < dim; ++i)
        {
            basis[i] = pCollExp[0]->GetBasis(i);
        }

        // Get shape type
        auto shapeType = pCollExp[0]->DetShapeType();

        // Generate operator string and create operator.
        std::string op_string = "Helmholtz";
        op_string += MatrixFree::GetOpstring(shapeType, m_isDeformed);
        auto oper = MatrixFree::GetOperatorFactory().CreateInstance(
            op_string, basis, pCollExp.size());

        // Set Jacobian
        oper->SetJac(pGeomData->GetJacInterLeave(pCollExp, m_nElmtPad));

        // Store derivative factor
        oper->SetDF(pGeomData->GetDerivFactorsInterLeave(pCollExp, m_nElmtPad));

        oper->SetUpBdata(basis);
        oper->SetUpDBdata(basis);
        oper->SetUpZW(basis);
        oper->SetUpD(basis);

        m_oper = std::dynamic_pointer_cast<MatrixFree::Helmholtz>(oper);
        ASSERTL0(m_oper, "Failed to cast pointer.");

        // Set factors
        m_factors   = StdRegions::NullFactorMap;
        m_varcoeffs = StdRegions::NullVarCoeffMap;
        this->UpdateFactors(factors);
    }
};

/// Factory initialisation for the Helmholtz_MatrixFree operators
OperatorKey Helmholtz_MatrixFree::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Hex"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Tet"),
};

} // namespace Nektar::Collections
