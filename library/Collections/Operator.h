///////////////////////////////////////////////////////////////////////////////
//
// File: Operator.h
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
// Description: Operator top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_COLLECTIONS_OPERATOR_H
#define NEKTAR_LIBRARY_COLLECTIONS_OPERATOR_H

#include <Collections/CollectionsDeclspec.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LocalRegions/Expansion.h>
#include <SpatialDomains/Geometry.h>

#define OPERATOR_CREATE(cname)                                                 \
    static OperatorKey m_type;                                                 \
    static OperatorKey m_typeArr[];                                            \
    friend class MemoryManager<cname>;                                         \
    static OperatorSharedPtr create(                                           \
        std::vector<LocalRegions::ExpansionSharedPtr> pCollExp,                \
        std::shared_ptr<CoalescedGeomData> GeomData,                           \
        StdRegions::FactorMap factors)                                         \
    {                                                                          \
        return MemoryManager<cname>::AllocateSharedPtr(pCollExp, GeomData,     \
                                                       factors);               \
    }

namespace Nektar::Collections
{

class CoalescedGeomData;
typedef std::shared_ptr<CoalescedGeomData> CoalescedGeomDataSharedPtr;

enum OperatorType
{
    eBwdTrans,
    eHelmholtz,
    eLinearAdvectionDiffusionReaction,
    eIProductWRTBase,
    eIProductWRTDerivBase,
    ePhysDeriv,
    ePhysInterp1DScaled,
    SIZE_OperatorType
};

const char *const OperatorTypeMap[] = {"BwdTrans",
                                       "Helmholtz",
                                       "LinearAdvectionDiffusionReaction",
                                       "IProductWRTBase",
                                       "IProductWRTDerivBase",
                                       "PhysDeriv",
                                       "PhysInterp1DScaled"};

const char *const OperatorTypeMap1[] = {
    "BwdTrans",   "Helmholtz", "LinearADR",         "IPWrtBase",
    "IPWrtDBase", "PhysDeriv", "PhysInterp1DScaled"};

enum ImplementationType
{
    eNoImpType,
    eNoCollection,
    eIterPerExp,
    eStdMat,
    eSumFac,
    eMatrixFree,
    SIZE_ImplementationType
};

const char *const ImplementationTypeMap[] = {"NoImplementationType",
                                             "NoCollection",
                                             "IterPerExp",
                                             "StdMat",
                                             "SumFac",
                                             "MatrixFree"};

const char *const ImplementationTypeMap1[] = {
    "NoImplementationType",
    "IterLocExp", // formerly "NoCollection"
    "IterStdExp", // formerly "IterPerExp"
    "StdMat    ",
    "SumFac    ",
    "MatFree   " // formerly "MatrixFree"
};

typedef bool ExpansionIsNodal;

class Operator;

/// Key for describing an Operator
typedef std::tuple<LibUtilities::ShapeType, OperatorType, ImplementationType,
                   ExpansionIsNodal>
    OperatorKey;

/// Operator factory definition
typedef Nektar::LibUtilities::NekFactory<
    OperatorKey, Operator, std::vector<LocalRegions::ExpansionSharedPtr>,
    CoalescedGeomDataSharedPtr, StdRegions::FactorMap>
    OperatorFactory;

/// Returns the singleton Operator factory object
OperatorFactory &GetOperatorFactory();

typedef std::map<OperatorType, ImplementationType> OperatorImpMap;

/// simple Operator Implementation Map generator
OperatorImpMap SetFixedImpType(ImplementationType defaultType);

/// Base class for operators on a collection of elements
class Operator
{
public:
    /// Constructor
    Operator(std::vector<LocalRegions::ExpansionSharedPtr> pCollExp,
             std::shared_ptr<CoalescedGeomData> GeomData,
             StdRegions::FactorMap factors);

    virtual ~Operator() = default;

    /// Perform operation
    COLLECTIONS_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output0, Array<OneD, NekDouble> &output1,
        Array<OneD, NekDouble> &output2,
        Array<OneD, NekDouble> &wsp = NullNekDouble1DArray) = 0;

    COLLECTIONS_EXPORT virtual void operator()(
        int dir, const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output,
        Array<OneD, NekDouble> &wsp = NullNekDouble1DArray) = 0;

    /// Update the supplied factor map
    COLLECTIONS_EXPORT virtual void UpdateFactors(
        [[maybe_unused]] StdRegions::FactorMap factors)
    {
        ASSERTL0(false, "This method needs to be re-implemented in derived "
                        "operator class.");
    }

    /// Update the supplied variable coefficients
    COLLECTIONS_EXPORT virtual void UpdateVarcoeffs(
        [[maybe_unused]] StdRegions::VarCoeffMap &varcoeffs)
    {
        ASSERTL0(false,
                 "This method needs to be re-implemented in derived "
                 "operator class. Make sure it is implemented for the operator"
                 " in the .opt file");
    }

    /// Get the size of the required workspace
    unsigned int GetWspSize()
    {
        return m_wspSize;
    }

    /// Get number of elements
    unsigned int GetNumElmt()
    {
        return m_numElmt;
    }

    /// Get expansion pointer
    StdRegions::StdExpansionSharedPtr GetExpSharedPtr()
    {
        return m_stdExp;
    }

    /*
     * Return the input size for this collection.
     * Optionally return the size for the opposite (Phys or Coeff) space.
     */
    inline unsigned int GetInputSize(void)
    {
        return m_inputSize;
    }

    /*
     * Return the output size for this collection.
     */
    inline unsigned int GetOutputSize(void)
    {
        return m_outputSize;
    }

    /*
     * Return the phys size for this collection.
     */
    inline unsigned int GetPhysSize(void)
    {
        return m_numElmt * m_stdExp->GetTotPoints();
    }

    /*
     * Return the coeff size for this collection.
     */
    inline unsigned int GetCoeffSize(void)
    {
        return m_numElmt * m_stdExp->GetNcoeffs();
    }

protected:
    bool m_isDeformed;
    StdRegions::StdExpansionSharedPtr m_stdExp;
    /// number of elements that the operator is applied on
    unsigned int m_numElmt;
    unsigned int m_nqe;
    unsigned int m_wspSize;
    /// number of modes or quadrature points that are passed as input to an
    /// operator
    unsigned int m_inputSize;
    /// number of modes or quadrature points  that are taken as output from an
    /// operator
    unsigned int m_outputSize;
};

/// Shared pointer to an Operator object
typedef std::shared_ptr<Operator> OperatorSharedPtr;

/// Less-than comparison operator for OperatorKey objects
bool operator<(OperatorKey const &p1, OperatorKey const &p2);

/// Stream output operator for OperatorKey objects
std::ostream &operator<<(std::ostream &os, OperatorKey const &p);

} // namespace Nektar::Collections

#endif
