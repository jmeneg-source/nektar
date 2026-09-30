///////////////////////////////////////////////////////////////////////////////
//
// File: Collection.h
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
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_COLLECTIONS_COLLECTION_H
#define NEKTAR_LIBRARY_COLLECTIONS_COLLECTION_H

#include <vector>

#include <Collections/CoalescedGeomData.h>
#include <Collections/CollectionsDeclspec.h>
#include <Collections/Operator.h>
#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <LocalRegions/Expansion.h>
#include <SpatialDomains/Geometry.h>

namespace Nektar::Collections
{

/**
 * @brief Collection
 */
class Collection
{
public:
    COLLECTIONS_EXPORT Collection(
        std::vector<LocalRegions::ExpansionSharedPtr> pCollExp,
        OperatorImpMap &impTypes);

    /// Update the factor map in operator @param opType by input @param factors
    COLLECTIONS_EXPORT void UpdateFactors(const OperatorType opType,
                                          StdRegions::FactorMap factors);

    COLLECTIONS_EXPORT void UpdateVarcoeffs(
        const OperatorType opType,
        StdRegions::VarCoeffMap &varcoeffs = StdRegions::NullVarCoeffMap);

    COLLECTIONS_EXPORT void Initialise(
        const OperatorType opType,
        StdRegions::FactorMap factors = StdRegions::NullFactorMap);

    inline void ApplyOperator(const OperatorType &op,
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &output);

    inline void ApplyOperator(const OperatorType &op,
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &output0,
                              Array<OneD, NekDouble> &output1);

    inline void ApplyOperator(const OperatorType &op,
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &output0,
                              Array<OneD, NekDouble> &output1,
                              Array<OneD, NekDouble> &output2);

    inline void ApplyOperator(const OperatorType &op, int dir,
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &output);

    inline bool HasOperator(const OperatorType &op);

    inline OperatorSharedPtr GetOpSharedPtr(const OperatorType &op)
    {
        return m_ops[op];
    }

    inline CoalescedGeomDataSharedPtr GetGeomSharedPtr()
    {
        return m_geomData;
    }
    inline std::vector<LocalRegions::ExpansionSharedPtr> &GetExpVector()
    {
        return m_collExp;
    }

    inline int GetInputSize(const OperatorType &op)
    {
        return m_ops[op]->GetInputSize();
    }

    inline int GetOutputSize(const OperatorType &op)
    {
        return m_ops[op]->GetOutputSize();
    }

    inline int GetPhysSize(const OperatorType &op)
    {
        return m_ops[op]->GetPhysSize();
    }

    inline int GetCoeffSize(const OperatorType &op)
    {
        return m_ops[op]->GetCoeffSize();
    }

    /// Return the number of elements in collection attached to @param op
    inline int GetNumElmt(const OperatorType &op)
    {
        return m_ops[op]->GetNumElmt();
    }

protected:
    std::unordered_map<OperatorType, OperatorSharedPtr> m_ops;
    CoalescedGeomDataSharedPtr m_geomData;
    // store details for initialisation on call rather than default
    // initialisation
    std::vector<LocalRegions::ExpansionSharedPtr> m_collExp;
    OperatorImpMap m_impTypes;
};

typedef std::vector<Collection> CollectionVector;
typedef std::shared_ptr<CollectionVector> CollectionVectorSharedPtr;

/**
 *
 */
inline void Collection::ApplyOperator(
    const OperatorType &op, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &output)
{
    Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
    (*m_ops[op])(inarray, output, NullNekDouble1DArray, NullNekDouble1DArray,
                 wsp);
}

/**
 *
 */
inline void Collection::ApplyOperator(
    const OperatorType &op, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &output0, Array<OneD, NekDouble> &output1)
{
    Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
    (*m_ops[op])(inarray, output0, output1, NullNekDouble1DArray, wsp);
}

/**
 *
 */
inline void Collection::ApplyOperator(
    const OperatorType &op, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &output0, Array<OneD, NekDouble> &output1,
    Array<OneD, NekDouble> &output2)
{
    Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
    (*m_ops[op])(inarray, output0, output1, output2, wsp);
}

/**
 *
 */
inline void Collection::ApplyOperator(
    const OperatorType &op, int dir,
    const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &output)
{
    Array<OneD, NekDouble> wsp(m_ops[op]->GetWspSize());
    (*m_ops[op])(dir, inarray, output, wsp);
}

inline bool Collection::HasOperator(const OperatorType &op)
{
    return (m_ops.find(op) != m_ops.end());
}

} // namespace Nektar::Collections

#endif
