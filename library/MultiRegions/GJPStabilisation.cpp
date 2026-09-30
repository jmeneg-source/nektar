///////////////////////////////////////////////////////////////////////////////
//
// File: GJPStabilisation.cpp
//
// For mre information, please see: http://www.nektar.info
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

#include <MultiRegions/GJPStabilisation.h>

namespace Nektar::MultiRegions
{
std::string GJPStabilisation::GJPStabilisationLookupIds[3] = {
    LibUtilities::SessionReader::RegisterEnumValue(
        "GJPStabilisation", "Explicit", eExplicitGJPStabilisation),
    LibUtilities::SessionReader::RegisterEnumValue(
        "GJPStabilisation", "SemiImplicit", eSemiImplicitGJPStabilisation),
    LibUtilities::SessionReader::RegisterEnumValue(
        "GJPStabilisation", "Implicit", eFullImplicitGJPStabilisation),
};
GJPStabilisation::GJPStabilisation(ExpListSharedPtr pField)
{
    LibUtilities::SessionReaderSharedPtr session = pField->GetSession();

    bool test;

    session->MatchSolverInfo("GJPStabilisation", "Explicit", test, false);
    if (test)
    {
        m_formulation = eGJPExplicit;
    }

    session->MatchSolverInfo("GJPStabilisation", "Implicit", test, false);
    if (test)
    {
        m_formulation = eGJPImplicit;

        ASSERTL0(session->MatchSolverInfo("GlobalSysSoln", "IterativeFull"),
                 "To use GJP Fully implicit stabilisation you must use a "
                 "Iterative Full solver");
    }

    session->MatchSolverInfo("GJPStabilisation", "SemiImplicit", test, false);
    if (test)
    {
        m_formulation = eGJPSemiImplicit;
    }
    ASSERTL0(m_formulation != eGJPNoFormulation,
             "Need a valid formualtion type for GradientJumpStabilisation: "
             "Explicit, Implicit, SemiImplicit");

    // Call GetTrace on the initialising field will set up
    // DG. Store a copy so that if we make a soft copy of
    // this class we can re-used this field for operators.
    pField->GetTrace();
    m_dgfield = pField;

    m_coordDim = m_dgfield->GetCoordim(0);
    m_traceDim = m_dgfield->GetShapeDimension() - 1;

    // set up trace normals but would be better if could use
    // equation system definition
    m_traceNormals = Array<OneD, Array<OneD, NekDouble>>(m_coordDim);
    for (int i = 0; i < m_coordDim; ++i)
    {
        m_traceNormals[i] =
            Array<OneD, NekDouble>(m_dgfield->GetTrace()->GetNpoints());
    }
    m_dgfield->GetTrace()->GetNormals(m_traceNormals);

    m_scalTrace = Array<OneD, Array<OneD, NekDouble>>(m_traceDim + 1);

    MultiRegions::ExpListSharedPtr dgtrace = m_dgfield->GetTrace();

    const std::shared_ptr<LocalRegions::ExpansionVector> exp =
        m_dgfield->GetExp();

    Array<OneD, Array<OneD, NekDouble>> dfactors[3];
    Array<OneD, NekDouble> e_tmp;

    m_nLocTracePts = 0;
    for (unsigned e = 0; e < (*exp).size(); ++e)
    {
        for (unsigned t = 0; t < (*exp)[e]->GetNtraces(); ++t)
        {
            m_nLocTracePts += (*exp)[e]->GetLocTraceExp(t)->GetTotPoints();
        }
    }

    m_scalTrace = Array<OneD, Array<OneD, NekDouble>>(m_traceDim + 1);
    for (int i = 0; i < m_traceDim + 1; ++i)
    {
        m_scalTrace[i] = Array<OneD, NekDouble>(m_nLocTracePts);
    }

    int cnt         = 0;
    int offset_phys = 0;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> dbasis;
    Array<OneD, Array<OneD, Array<OneD, unsigned int>>> traceToCoeffMap;

    Array<OneD, int> sign, sign1;
    NekDouble h, p;
    std::map<unsigned, std::pair<NekDouble, unsigned>> hpscale;

    MultiRegions::DisContFieldSharedPtr dgfield =
        std::dynamic_pointer_cast<MultiRegions::DisContField>(m_dgfield);

    for (int e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        LocalRegions::ExpansionSharedPtr elmt = (*exp)[e];

        elmt->NormalTraceDerivFactors(dfactors[0], dfactors[1], dfactors[2]);

        for (int n = 0; n < elmt->GetNtraces(); ++n, ++cnt)
        {
            // collect offset for traces in dgtrace
            unsigned eid = dgfield->GetTraceElmtId(e, n);
            m_traceOffset.push_back(dgtrace->GetPhys_Offset(eid));
            LocalRegions::ExpansionSharedPtr LocTraceExp =
                elmt->GetLocTraceExp(n);
            unsigned LocTracepts = LocTraceExp->GetTotPoints();
            if (LocTracepts != dgtrace->GetExp(eid)->GetTotPoints())
            {
                m_interpTrace[cnt] =
                    std::make_pair(LocTraceExp, dgtrace->GetExp(eid));
            }

            m_locTracePts0.push_back(LocTraceExp->GetNumPoints(0));
            if (m_traceDim == 1)
            {
                m_locTracePts1.push_back(1);
            }
            else
            {
                m_locTracePts1.push_back(LocTraceExp->GetNumPoints(1));
            }
            m_traceFwd.push_back(dgfield->IsLeftAdjacentTrace(e, n));

            // gather h and p average informatoin
            elmt->TraceNormLen(n, h, p);
            unsigned edgeid = elmt->GetTraceExp(n)->GetGeom()->GetGlobalID();
            if (hpscale.count(edgeid))
            {
                auto hp         = hpscale[edgeid];
                hpscale[edgeid] = std::pair<NekDouble, unsigned>(
                    0.5 * (hp.first + h), (int)(0.5 * (hp.second + p)));
            }
            else
            {
                hpscale[edgeid] = std::pair<NekDouble, unsigned>(h, p);
            }

            int nptrace = dfactors[0][n].size();
            for (int i = 0; i < m_traceDim + 1; ++i)
            {
                Vmath::Smul(nptrace, 1.0, dfactors[i][n], 1,
                            e_tmp = m_scalTrace[i] + offset_phys, 1);
            }
            offset_phys += nptrace;
        }
        m_ntrace.push_back(elmt->GetNtraces());
    }
    m_nLocTracePts = offset_phys;

    m_locTraceWeights = Array<OneD, NekDouble>(m_nLocTracePts);
    // set up scale factors
    for (int e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        LocalRegions::ExpansionSharedPtr elmt = (*exp)[e];
        for (int n = 0; n < elmt->GetNtraces(); ++n, ++cnt)
        {
            unsigned edgeid = elmt->GetTraceExp(n)->GetGeom()->GetGlobalID();

            ASSERTL1(hpscale.count(edgeid), "Scale has not been defined");
            auto hp     = hpscale[edgeid];
            NekDouble h = hp.first;
            unsigned p  = hp.second;
            NekDouble jumpScal =
                (p == 1) ? 0.02 * h * h : 0.8 * pow(p + 1, -4.0) * h * h;

            m_locEdgeScale.push_back(jumpScal);
        }
    }

    m_locTraceWeights = Array<OneD, NekDouble>(m_nLocTracePts, 1.0);

    //  Generate array of quadrature and Jacobian for loc Trace
    unsigned offset = 0;
    cnt             = 0;
    for (int e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        for (unsigned j = 0; j < m_ntrace[e]; ++j)
        {
            unsigned eid =
                std::dynamic_pointer_cast<MultiRegions::DisContField>(m_dgfield)
                    ->GetTraceElmtId(e, j);

            if (m_interpTrace.count(cnt)) // manage interpolated case
            {
                m_interpTrace[cnt].first->MultiplyByQuadratureMetric(
                    m_locTraceWeights + offset,
                    e_tmp = m_locTraceWeights + offset);
            }
            else
            {
                dgtrace->GetExp(eid)->MultiplyByQuadratureMetric(
                    m_locTraceWeights + offset,
                    e_tmp = m_locTraceWeights + offset);

                // reverse data if necessary
                m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                    m_dgfield->GetExp(e)->GetTraceOrient(j), e_tmp, e_tmp,
                    m_locTracePts0[cnt], m_locTracePts1[cnt], false);
            }
            offset += m_locTracePts0[cnt] * m_locTracePts1[cnt];
            cnt++;
        }
    }
    // Assemble list of Matrix Product
    Array<OneD, DNekMatSharedPtr> TraceMat;

    m_dgfield->GetExp(0)->StdDerivBaseOnTraceMat(TraceMat);

    int nelmt = 1;
    Array<OneD, const LibUtilities::BasisSharedPtr> base_sav =
        m_dgfield->GetExp(0)->GetBase();

    for (int n = 1; n < m_dgfield->GetExpSize(); ++n)
    {
        const Array<OneD, const LibUtilities::BasisSharedPtr> &base =
            m_dgfield->GetExp(n)->GetBase();

        // check to see if same element expansion as previous matrix and if
        // so can reused
        int i;
        for (i = 0; i < base.size(); ++i)
        {
            if (base[i] != base_sav[i])
            {
                break;
            }
        }

        if (i == base.size())
        {
            nelmt++;
        }
        else
        {
            // save previous block of data.
            m_StdDBaseOnTraceMat.push_back(
                std::pair<int, Array<OneD, DNekMatSharedPtr>>(nelmt, TraceMat));

            // start new block
            m_dgfield->GetExp(n)->StdDerivBaseOnTraceMat(TraceMat);
            nelmt    = 1;
            base_sav = m_dgfield->GetExp(n)->GetBase();
        }
    }

    // save latest block of data.
    m_StdDBaseOnTraceMat.push_back(
        std::pair<int, Array<OneD, DNekMatSharedPtr>>(nelmt, TraceMat));
}

void GJPStabilisation::Apply(const Array<OneD, NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray,
                             const Array<OneD, NekDouble> &pUnorm,
                             NekDouble scale) const
{
    int ncoeffs   = m_dgfield->GetNcoeffs();
    int nTracePts = m_dgfield->GetTrace()->GetTotPoints();
    LibUtilities::Timer timer, timer1;
    timer.Start();
    timer1.Start();

    Array<OneD, NekDouble> FilterCoeffs(ncoeffs), CoeffsTmp(ncoeffs);
    Array<OneD, NekDouble> Fwd(nTracePts, 0.0), Bwd(nTracePts, 0.0);
    Array<OneD, NekDouble> Store(m_nLocTracePts);
    Array<OneD, NekDouble> unorm;

    if (pUnorm == NullNekDouble1DArray)
    {
        unorm = Array<OneD, NekDouble>(nTracePts, 1.0);
    }
    else
    {
        unorm = pUnorm;
    }

    Array<OneD, NekDouble> dudn(m_nLocTracePts, 0.0);

    timer1.Stop();
    timer1.AccumulateRegion("GJP:Init", 10);
    timer1.Start();

    // Fwd Trans
    m_dgfield->FwdTransLocalElmt(inarray, FilterCoeffs);

    // Construct local derivaitves on trace
    for (int n = 0; n < m_coordDim; ++n)
    {
        StdDerivOnTraceFromModes(n, FilterCoeffs, Store);
        Vmath::Vvtvp(m_nLocTracePts, Store, 1, m_scalTrace[n], 1, dudn, 1, dudn,
                     1);
    }

    timer1.Stop();
    timer1.AccumulateRegion("GJP:FwdTrans + StdDerivFromModes", 10);
    timer1.Start();

    if (m_formulation == eGJPSemiImplicit)
    {
        // want to put Fwd vals on bwd trace and vice versa
        TraceJumpFromLocTraceNormDeriv(dudn, Bwd, Fwd);
    }
    else
    {
        TraceJumpFromLocTraceNormDeriv(dudn, Fwd, Bwd);
        Vmath::Vadd(nTracePts, Fwd, 1, Bwd, 1, Fwd, 1);
    }
    timer1.Stop();
    timer1.AccumulateRegion("GJP:TransJumpFromLocDeriv", 10);
    timer1.Start();

    if (m_formulation == eGJPSemiImplicit)
    {
        Vmath::Vmul(nTracePts, unorm, 1, Fwd, 1, Fwd, 1);
        Vmath::Vmul(nTracePts, unorm, 1, Bwd, 1, Bwd, 1);
        // Evaluate trace inner product with respect to derivative of the
        // basis
        ConstructLocalTraceJumpSI(0, Fwd, Bwd, Store);
        timer1.Stop();
        timer1.AccumulateRegion("GJP:Construct Trace ", 10);
        timer1.Start();
        IProductwrtStdDerivBaseOnTraceMat(0, Store, FilterCoeffs);
        timer1.Stop();
        timer1.AccumulateRegion("GJP:Deriv on Trace", 10);

        timer1.Start();

        for (int i = 0; i < m_traceDim; ++i)
        {
            ConstructLocalTraceJumpSI(i + 1, Fwd, Bwd, Store);
            timer1.Stop();
            timer1.AccumulateRegion("GJP:Construct Trace ", 10);

            timer1.Start();
            IProductwrtStdDerivBaseOnTraceMat(i + 1, Store, CoeffsTmp);
            timer1.Stop();
            timer1.AccumulateRegion("GJP:Deriv on Trace", 10);

            timer1.Start();
            Vmath::Vadd(ncoeffs, CoeffsTmp, 1, FilterCoeffs, 1, FilterCoeffs,
                        1);
        }
    }
    else
    {
        Vmath::Vmul(nTracePts, unorm, 1, Fwd, 1, Fwd, 1);

        // Evaluate trace inner product with respect to derivitive of the
        // basis
        ConstructLocalTraceJump(0, Fwd, Store);
        IProductwrtStdDerivBaseOnTraceMat(0, Store, FilterCoeffs);

        for (int i = 0; i < m_traceDim; ++i)
        {
            ConstructLocalTraceJump(i + 1, Fwd, Store);
            IProductwrtStdDerivBaseOnTraceMat(i + 1, Store, CoeffsTmp);
            Vmath::Vadd(ncoeffs, CoeffsTmp, 1, FilterCoeffs, 1, FilterCoeffs,
                        1);
        }
    }

    Vmath::Svtvp(ncoeffs, scale, FilterCoeffs, 1, outarray, 1, outarray, 1);
    timer.Stop();
    // Elapsed time
    timer.AccumulateRegion("GJP:Total", 10);
}

/**
 * construct the gradient jump on local trace using input
 * 'in' which  is in trace orientation
 */
void GJPStabilisation::ConstructLocalTraceJump(
    const int dir, const Array<OneD, const NekDouble> &in,
    Array<OneD, NekDouble> &store) const
{
    LibUtilities::Timer timer;
    unsigned offset = 0;
    unsigned cnt    = 0;
    Array<OneD, NekDouble> tmp;

    for (unsigned e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        for (unsigned j = 0; j < m_ntrace[e]; ++j)
        {
            unsigned locTracePts = m_locTracePts0[cnt] * m_locTracePts1[cnt];

            NekDouble jumpScal = m_locEdgeScale[cnt];

            if (m_interpTrace.count(cnt)) // variable p and BCs
            {
                auto it = m_interpTrace.find(cnt);

                StdRegions::Orientation orient =
                    m_dgfield->GetExp(e)->GetTraceOrient(j);

                // interpolate to new space
                it->second.first->PhysInterp(
                    it->second.second, in + m_traceOffset[cnt],
                    tmp = store + offset,
                    orient >= StdRegions::eDir1FwdDir2_Dir2FwdDir1);

                // reorientate
                m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                    orient, tmp, tmp, m_locTracePts0[cnt], m_locTracePts1[cnt],
                    false);
            }
            else
            {
                tmp = store + offset;
                timer.Start();
                // Reverse data if necessary - using mapping from trace to local
                // trace
                m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                    m_dgfield->GetExp(e)->GetTraceOrient(j),
                    in + m_traceOffset[cnt], tmp, m_locTracePts0[cnt],
                    m_locTracePts1[cnt], false);
                timer.Stop();
                timer.AccumulateRegion("GJPStab::Construct Trace:: Reorient",
                                       10);
            }

            // multiply by local factors in local trace space
            for (unsigned i = 0; i < locTracePts; ++i)
            {
                tmp[i] *= jumpScal * m_scalTrace[dir][offset + i] *
                          m_locTraceWeights[offset + i];
            }

            offset += locTracePts;
            cnt++;
        }
    }
}

void GJPStabilisation::ConstructLocalTraceJumpSI(
    const int dir, const Array<OneD, const NekDouble> &Fwd,
    const Array<OneD, const NekDouble> &Bwd,
    Array<OneD, NekDouble> &store) const
{
    unsigned offset = 0;
    unsigned cnt    = 0;
    Array<OneD, NekDouble> tmp;
    LibUtilities::Timer timer;

    for (unsigned e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        for (unsigned j = 0; j < m_ntrace[e]; ++j)
        {
            unsigned locTracePts = m_locTracePts0[cnt] * m_locTracePts1[cnt];

            NekDouble jumpScal = m_locEdgeScale[cnt];

            if (m_traceFwd[cnt])
            {
                if (m_interpTrace.count(cnt)) // variable p and BCs
                {
                    auto it = m_interpTrace.find(cnt);

                    StdRegions::Orientation orient =
                        m_dgfield->GetExp(e)->GetTraceOrient(j);

                    // interpolate to new space
                    it->second.first->PhysInterp(
                        it->second.second, Fwd + m_traceOffset[cnt],
                        tmp = store + offset,
                        orient >= StdRegions::eDir1FwdDir2_Dir2FwdDir1);

                    // reorientate
                    m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                        orient, tmp, tmp, m_locTracePts0[cnt],
                        m_locTracePts1[cnt], false);
                }
                else
                {
                    // Reverse data if necessary - using mapping from trace to
                    // local trace
                    timer.Start();
                    m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                        m_dgfield->GetExp(e)->GetTraceOrient(j),
                        Fwd + m_traceOffset[cnt], tmp = store + offset,
                        m_locTracePts0[cnt], m_locTracePts1[cnt], false);
                    timer.Stop();
                    timer.AccumulateRegion(
                        "GJPStab::Construct Trace:: Reorient", 10);
                }
            }
            else
            {
                if (m_interpTrace.count(cnt)) // variable p and BCs
                {
                    auto it = m_interpTrace.find(cnt);

                    StdRegions::Orientation orient =
                        m_dgfield->GetExp(e)->GetTraceOrient(j);

                    // interpolate to new space
                    it->second.first->PhysInterp(
                        it->second.second, Fwd + m_traceOffset[cnt],
                        tmp = store + offset,
                        orient >= StdRegions::eDir1FwdDir2_Dir2FwdDir1);

                    // reorientate
                    m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                        orient, tmp, tmp, m_locTracePts0[cnt],
                        m_locTracePts1[cnt], false);
                }
                else
                {
                    // Reverse data if necessary - using mapping from trace to
                    // local trace
                    timer.Start();
                    m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                        m_dgfield->GetExp(e)->GetTraceOrient(j),
                        Bwd + m_traceOffset[cnt], tmp = store + offset,
                        m_locTracePts0[cnt], m_locTracePts1[cnt], false);
                    timer.Stop();
                    timer.AccumulateRegion(
                        "GJPStab::Construct Trace:: Reorient", 10);
                }
            }
            // multiply by local factors in local trace space
            for (unsigned i = 0; i < locTracePts; ++i)
            {
                tmp[i] *= jumpScal * m_scalTrace[dir][offset + i] *
                          m_locTraceWeights[offset + i];
            }

            offset += locTracePts;
            cnt++;
        }
    }
    ASSERTL1(offset <= 2 * m_dgfield->GetTrace()->GetTotPoints(),
             "Stroage is not large enough");
}

/* Given the local trace point of a function whcih include the integration
 * weights return the Inner product with respect to the Std Derivative in the
 * dir direction */
void GJPStabilisation::IProductwrtStdDerivBaseOnTraceMat(
    int dir, Array<OneD, NekDouble> &in, Array<OneD, NekDouble> &out) const
{
    unsigned cnt  = 0;
    unsigned cnt1 = 0;

    for (auto &it : m_StdDBaseOnTraceMat)
    {
        unsigned modes    = it.second[dir]->GetRows();
        unsigned tracepts = it.second[dir]->GetColumns();

        Blas::Dgemm('N', 'N', modes, it.first, tracepts, 1.0,
                    &(it.second[dir]->GetPtr())[0], modes, &in[0] + cnt,
                    tracepts, 0.0, &out[0] + cnt1, modes);

        cnt += tracepts * it.first;
        cnt1 += modes * it.first;
    }
}

void GJPStabilisation::TraceJumpFromLocTraceNormDeriv(
    Array<OneD, NekDouble> &normderiv, Array<OneD, NekDouble> &Fwd,
    Array<OneD, NekDouble> &Bwd) const
{
    unsigned offset = 0;
    unsigned cnt    = 0;
    Array<OneD, NekDouble> tmp, tmp1;

    for (unsigned e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        for (unsigned j = 0; j < m_ntrace[e]; ++j)
        {
            unsigned locTracePts = m_locTracePts0[cnt] * m_locTracePts1[cnt];

            tmp = normderiv + offset;
            ASSERTL0(offset < normderiv.size(), "Issue");
            // Reverse data if necessary
            m_dgfield->GetExp(e)->ReOrientTracePhysVals(
                m_dgfield->GetExp(e)->GetTraceOrient(j), tmp, tmp,
                m_locTracePts0[cnt], m_locTracePts1[cnt]);

            if (m_traceFwd[cnt])
            {
                tmp1 = Fwd + m_traceOffset[cnt];
            }
            else
            {
                tmp1 = Bwd + m_traceOffset[cnt];
            }

            if (m_interpTrace.count(cnt)) // variable p and BCs
            {
                auto it = m_interpTrace.find(cnt);
                it->second.second->PhysInterp(it->second.first, tmp, tmp1);
            }
            else
            {
                Vmath::Vcopy(locTracePts, tmp, 1, tmp1, 1);
            }

            offset += locTracePts;
            cnt++;
        }
    }

    // globally assemble Fwd and Bwd traces;
    m_dgfield->PeriodicBwdCopy(Fwd, Bwd);
    m_dgfield->FillBwdWithBoundCond(Fwd, Bwd, true);
    m_dgfield->GetTraceMap()->GetAssemblyCommDG()->PerformExchange(Fwd, Bwd);
}

/* Given the modes of an element on input  return the std derivative in
 * direction dir on the trace at local tracepoints
 */
void GJPStabilisation::StdDerivOnTraceFromModes(
    int dir, Array<OneD, NekDouble> &in, Array<OneD, NekDouble> &out) const
{
    unsigned cnt  = 0;
    unsigned cnt1 = 0;

    for (auto &it : m_StdDBaseOnTraceMat)
    {
        unsigned modes    = it.second[dir]->GetRows();
        unsigned tracepts = it.second[dir]->GetColumns();

        Blas::Dgemm('T', 'N', tracepts, it.first, modes, 1.0,
                    &(it.second[dir]->GetPtr())[0], modes, &in[0] + cnt, modes,
                    0.0, &out[0] + cnt1, tracepts);

        cnt += modes * it.first;
        cnt1 += tracepts * it.first;
    }
}

Array<OneD, NekDouble> GJPStabilisation::GetTraceWeightVarFactors(void)
{

    Array<OneD, NekDouble> returnval(m_dgfield->GetNumElmts() * 6, 0.0);

    unsigned cnt = 0;
    for (unsigned e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        for (unsigned n = 0; n < m_dgfield->GetExp(e)->GetNtraces(); ++n)
        {
            returnval[6 * e + n] = m_locEdgeScale[cnt++];
        }
    }
    return returnval;
}

} // namespace Nektar::MultiRegions
