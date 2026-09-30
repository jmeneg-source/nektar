///////////////////////////////////////////////////////////////////////////////
//
// File: GJPStabilisation.h
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
// Description: GJP data
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_GJPSTABILISATION_H
#define NEKTAR_LIBS_MULTIREGIONS_GJPSTABILISATION_H

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/DisContField.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/MultiRegionsDeclspec.h>

namespace Nektar::MultiRegions
{

enum GJPFormulation
{
    eGJPNoFormulation,
    eGJPExplicit,
    eGJPImplicit,
    eGJPSemiImplicit
};

class GJPStabilisation
{
public:
    MULTI_REGIONS_EXPORT GJPStabilisation(ExpListSharedPtr field);

    MULTI_REGIONS_EXPORT ~GJPStabilisation(){};

    MULTI_REGIONS_EXPORT void Apply(
        const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray,
        const Array<OneD, NekDouble> &pUnorm = NullNekDouble1DArray,
        const NekDouble scale                = 1.0) const;

    Array<OneD, Array<OneD, NekDouble>> &GetTraceNormals(void)
    {
        return m_traceNormals;
    }

    int GetNumTracePts(void) const
    {
        return m_dgfield->GetTrace()->GetTotPoints();
    }

    bool IsSemiImplicit() const
    {
        return (m_formulation == eGJPSemiImplicit);
    }

    bool IsExplicit() const
    {
        return (m_formulation == eGJPExplicit);
    }

    bool IsImplicit() const
    {
        return (m_formulation == eGJPImplicit);
    }

    MULTI_REGIONS_EXPORT Array<OneD, NekDouble> GetTraceWeightVarFactors(void);

private:
    unsigned m_coordDim;
    unsigned m_traceDim;
    unsigned m_nLocTracePts;
    GJPFormulation m_formulation = eGJPNoFormulation;
    static std::string GJPStabilisationLookupIds[];

    // Trace normals
    Array<OneD, Array<OneD, NekDouble>> m_traceNormals;

    /// DG expansion for projection evalaution along trace
    MultiRegions::ExpListSharedPtr m_dgfield;

    /// Scale factor for phys values along trace involving the local
    /// normals and geometric factors
    Array<OneD, Array<OneD, NekDouble>> m_scalTrace;

    Array<OneD, NekDouble> m_locTraceWeights;

    /// phys offset in trace expannsion of each trace as we loop over elmts
    std::vector<unsigned> m_traceOffset;
    /// npoints in local trace expannsion in dir 0
    std::vector<unsigned> m_locTracePts0;
    /// npoints in local trace expannsion in dir 1
    std::vector<unsigned> m_locTracePts1;
    /// list of the number of traces over an element;
    std::vector<unsigned> m_ntrace;
    /// local trace and multiregion dg trace  if different (i.e variable p and
    /// BC trace)
    std::map<int, std::pair<LocalRegions::ExpansionSharedPtr,
                            LocalRegions::ExpansionSharedPtr>>
        m_interpTrace;

    std::vector<bool> m_traceFwd;

    std::vector<std::pair<int, Array<OneD, DNekMatSharedPtr>>>
        m_StdDBaseOnTraceMat;

    std::vector<NekDouble> m_locEdgeScale;

    void ConstructLocalTraceJump(const int dir,
                                 const Array<OneD, const NekDouble> &in,
                                 Array<OneD, NekDouble> &store) const;

    void ConstructLocalTraceJumpSI(const int dir,
                                   const Array<OneD, const NekDouble> &Fwd,
                                   const Array<OneD, const NekDouble> &Bwd,
                                   Array<OneD, NekDouble> &store) const;

    void IProductwrtStdDerivBaseOnTraceMat(int i, Array<OneD, NekDouble> &in,
                                           Array<OneD, NekDouble> &out) const;

    void StdDerivOnTraceFromModes(int i, Array<OneD, NekDouble> &in,
                                  Array<OneD, NekDouble> &out) const;

    void TraceJumpFromLocTraceNormDeriv(Array<OneD, NekDouble> &normderiv,
                                        Array<OneD, NekDouble> &Fwd,
                                        Array<OneD, NekDouble> &Bwd) const;
};

typedef std::shared_ptr<GJPStabilisation> GJPStabilisationSharedPtr;

} // namespace Nektar::MultiRegions
#endif // GJP
