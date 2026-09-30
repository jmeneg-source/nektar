///////////////////////////////////////////////////////////////////////////////
//
// File: FilterReynoldsStresses.cpp
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
// Description: Append Reynolds stresses to the average fields
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Filters/FilterReynoldsStresses.h>

namespace Nektar::SolverUtils
{

std::string FilterReynoldsStresses::className =
    GetFilterFactory().RegisterCreatorFunction("ReynoldsStresses",
                                               FilterReynoldsStresses::create);

/**
 * @class FilterReynoldsStresses
 *
 * @brief Append Reynolds stresses to the average fields
 *
 * This class appends the average fields with the Reynolds stresses of the form
 * \f$ \overline{u' v'} \f$.
 *
 * For the default case, this is achieved by calculating
 * \f$ C_{n} = \Sigma_{i=1}^{n} (u_i - \bar{u}_n)(v_i - \bar{v}_n)\f$
 * using the recursive relation:
 *
 * \f[ C_{n} = C_{n-1} + \frac{n}{n-1} (u_n - \bar{u}_n)(v_n - \bar{v}_n) \f]
 *
 * The FilterSampler base class then divides the result by n, leading
 * to the Reynolds stress.
 *
 * It is also possible to perform the averages using an exponential moving
 *  average, in which case either the moving average parameter \f$ \alpha \f$
 * or the time constant \f$ \tau \f$ must be prescribed.
 */
FilterReynoldsStresses::FilterReynoldsStresses(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::shared_ptr<SolverUtils::EquationSystem> &pEquation,
    const std::map<std::string, std::string> &pParams)
    : FilterFieldConvert(pSession, pEquation, pParams)
{
    // Load sampling frequency
    auto it = pParams.find("SampleFrequency");
    if (it == pParams.end())
    {
        m_sampleFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_sampleFrequency = round(equ.Evaluate());
    }

    // Check if we should use moving average
    it = pParams.find("MovingAverage");
    if (it == pParams.end())
    {
        m_movAvg = false;
    }
    else
    {
        std::string sOption = it->second.c_str();
        m_movAvg            = (boost::iequals(sOption, "true")) ||
                   (boost::iequals(sOption, "yes"));
    }

    // Check if highOrder is switched on
    it = pParams.find("ScaleNumModes");
    if (it == pParams.end())
    {
        m_Scale = false;
    }
    else
    {
        m_Scale = true;
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_ScaleNumModes = equ.Evaluate();
    }

    if (!m_movAvg)
    {
        return;
    }

    // Load alpha parameter for moving average
    it = pParams.find("alpha");
    if (it == pParams.end())
    {
        it = pParams.find("tau");
        if (it == pParams.end())
        {
            ASSERTL0(false, "MovingAverage needs either alpha or tau.");
        }
        else
        {
            // Load time constant
            LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
            NekDouble tau = equ.Evaluate();
            // Load delta T between samples
            NekDouble dT;
            m_session->LoadParameter("TimeStep", dT);
            dT = dT * m_sampleFrequency;
            // Calculate alpha
            m_alpha = dT / (tau + dT);
        }
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_alpha = equ.Evaluate();
        // Check if tau was also defined
        it = pParams.find("tau");
        if (it != pParams.end())
        {
            NEKERROR(ErrorUtil::efatal,
                     "Cannot define both alpha and tau in MovingAverage");
        }
    }
    // Check bounds of m_alpha
    ASSERTL0(m_alpha > 0 && m_alpha < 1, "Alpha out of bounds.");
}

void FilterReynoldsStresses::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    size_t dim          = pFields.size() - 1;
    size_t nExtraFields = (dim + 1) * dim / 2;
    size_t origFields   = pFields.size();
    size_t nqtot        = pFields[0]->GetTotPoints();
    bool waveSpace      = pFields[0]->GetWaveSpace();

    if (m_Scale)
    {
        // Getting mesh
        SpatialDomains::MeshGraphSharedPtr graph = pFields[0]->GetGraph();

        // Initialising high order m_pFieldsInterp
        int nvariables = m_session->GetVariables().size();
        m_pFieldsScaled =
            Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);

        const SpatialDomains::ExpansionInfoMap expInfo =
            graph->GetExpansionInfo(m_session->GetVariable(0));

        SpatialDomains::ExpansionInfoMapShPtr expInfoScaled = MemoryManager<
            SpatialDomains::ExpansionInfoMap>::AllocateSharedPtr();

        for (auto expIt = expInfo.begin(); expIt != expInfo.end(); ++expIt)
        {
            int expSpDim = expIt->second->m_basisKeyVector.size();

            std::vector<LibUtilities::BasisKey> BKeyVector;
            std::vector<int> oldPts(expSpDim);

            for (int i = 0; i < expSpDim; ++i)
            {
                LibUtilities::BasisKey bkey =
                    expIt->second->m_basisKeyVector[i];
                oldPts[i] = bkey.GetNumPoints();
            }

            for (int i = 0; i < expSpDim; ++i)
            {
                LibUtilities::BasisKey bkeyold =
                    expIt->second->m_basisKeyVector[i];
                int newNumModes =
                    static_cast<int>(m_ScaleNumModes * bkeyold.GetNumModes());

                int npts;
                // 1D
                if (i == 0)
                {
                    npts = static_cast<int>(oldPts[0] * m_ScaleNumModes);
                }
                else // 2D and 3D
                {
                    npts =
                        (oldPts[0] - oldPts[i] == 1)
                            ? static_cast<int>(oldPts[0] * m_ScaleNumModes - 1)
                            : static_cast<int>(oldPts[i] * m_ScaleNumModes);
                }

                const LibUtilities::PointsKey pkey(npts,
                                                   bkeyold.GetPointsType());
                LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(),
                                               newNumModes, pkey);
                BKeyVector.push_back(bkeynew);
            }

            (*expInfoScaled)[expIt->first] =
                MemoryManager<SpatialDomains::ExpansionInfo>::AllocateSharedPtr(
                    expIt->second->m_geomPtr, BKeyVector);
        }

        graph->SetExpansionInfo("Highorder", expInfoScaled);

        for (int i = 0; i < pFields.size(); ++i)
        {
            if (waveSpace)
            {
                ASSERTL0(pFields[0]->GetExpType() == MultiRegions::e3DH1D,
                         "Nummodes scaling for Reynolds stresses is "
                         "implemented for 3DH1D expansions only.");

                int npointsZ;
                NekDouble LhomZ;
                bool useFFT, homogen_dealiasing;

                m_session->LoadParameter("HomModesZ", npointsZ);
                m_session->LoadParameter("LZ", LhomZ);
                m_session->MatchSolverInfo("USEFFT", "FFTW", useFFT, false);
                m_session->MatchSolverInfo("DEALIASING", "True",
                                           homogen_dealiasing, false);

                const LibUtilities::PointsKey PkeyZ(
                    npointsZ, LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey BkeyZ(LibUtilities::eFourier,
                                                   npointsZ, PkeyZ);

                m_pFieldsScaled[i] =
                    MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                        AllocateSharedPtr(m_session, BkeyZ, LhomZ, useFFT,
                                          homogen_dealiasing, graph,
                                          "Highorder", Collections::eNoImpType);
            }
            else
            {
                m_pFieldsScaled[i] =
                    MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(
                        m_session, graph, false, "Highorder",
                        Collections::eNoImpType);
            }
        }

        nqtot = m_pFieldsScaled[0]->GetTotPoints();
        WARNINGL0(nqtot != pFields[0]->GetTotPoints(),
                  "The scaled number of modes did not increase in the Reynolds "
                  "Stress filter, "
                  "please increase scaled factor further.")
    }

    // Allocate storage
    m_fields.resize(origFields + nExtraFields);
    m_delta.resize(dim);
    // Initialising average fields
    for (size_t n = 0; n < m_fields.size(); ++n)
    {
        m_fields[n] = Array<OneD, NekDouble>(nqtot, 0.0);
    }
    // Initialising fluctuating fields
    for (size_t n = 0; n < m_delta.size(); ++n)
    {
        m_delta[n] = Array<OneD, NekDouble>(nqtot, 0.0);
    }

    // Initialise output arrays
    if (m_Scale)
    {
        FilterFieldConvert::v_Initialise(m_pFieldsScaled, time);
    }
    else
    {
        FilterFieldConvert::v_Initialise(pFields, time);
    }

    // Update m_fields if using restart file
    if (m_numSamples)
    {
        for (size_t j = 0; j < m_fields.size(); ++j)
        {
            pFields[0]->BwdTrans(m_outFields[j], m_fields[j]);
            if (pFields[0]->GetWaveSpace())
            {
                pFields[0]->HomogeneousBwdTrans(nqtot, m_fields[j],
                                                m_fields[j]);
            }
        }
    }
}

void FilterReynoldsStresses::v_FillVariablesName(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    size_t dim        = pFields.size() - 1;
    size_t origFields = pFields.size();

    // Fill name of variables
    for (size_t n = 0; n < origFields; ++n)
    {
        m_variables.push_back(pFields[n]->GetSession()->GetVariable(n));
    }
    for (int i = 0; i < dim; ++i)
    {
        for (int j = i; j < dim; ++j)
        {
            std::string var = pFields[i]->GetSession()->GetVariable(i) +
                              pFields[j]->GetSession()->GetVariable(j);
            m_variables.push_back(var);
        }
    }
}

void FilterReynoldsStresses::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    [[maybe_unused]] const NekDouble &time)
{
    size_t i, j, n;
    size_t dim         = pFields.size() - 1;
    bool waveSpace     = pFields[0]->GetWaveSpace();
    NekDouble nSamples = (NekDouble)m_numSamples;
    size_t nq          = pFields[0]->GetTotPoints();
    size_t ncoeffs     = pFields[0]->GetNcoeffs();

    if (m_Scale)
    {
        nq      = m_pFieldsScaled[0]->GetTotPoints();
        ncoeffs = m_pFieldsScaled[0]->GetNcoeffs();
    }

    // For moving average, take first sample as initial vector
    NekDouble alpha = m_alpha;
    if (m_numSamples == 1)
    {
        alpha = 1.0;
    }

    // Define auxiliary constants for averages
    NekDouble facOld, facAvg, facStress, facDelta;
    if (m_movAvg)
    {
        facOld    = 1.0 - alpha;
        facAvg    = alpha;
        facStress = alpha;
        facDelta  = 1.0;
    }
    else
    {
        facOld    = 1.0;
        facAvg    = 1.0;
        facStress = nSamples / (nSamples - 1);
        facDelta  = 1.0 / nSamples;
    }

    Array<OneD, NekDouble> vel(nq);
    Array<OneD, NekDouble> tmp(nq);

    // Update original velocities in phys space and calculate (\bar{u} - u_n)
    for (n = 0; n < dim; ++n)
    {
        if (waveSpace)
        {
            if (m_Scale)
            {
                // BwdTrans into phys space before interpolating onto high-order
                // grid
                Array<OneD, NekDouble> phys(pFields[n]->GetTotPoints());
                pFields[n]->HomogeneousBwdTrans(pFields[n]->GetTotPoints(),
                                                pFields[n]->GetPhys(), phys);
                pFields[0]->PhysInterp1DScaled(m_ScaleNumModes, phys, vel);
            }
            else
            {
                pFields[n]->HomogeneousBwdTrans(nq, pFields[n]->GetPhys(), vel);
            }
        }
        else
        {
            if (m_Scale)
            {
                // Interpolate phys-field by ScaleNumModes
                pFields[0]->PhysInterp1DScaled(m_ScaleNumModes,
                                               pFields[n]->GetPhys(), vel);
            }
            else
            {
                vel = pFields[n]->GetPhys();
            }
        }

        Vmath::Svtsvtp(nq, facAvg, vel, 1, facOld, m_fields[n], 1, m_fields[n],
                       1);
        Vmath::Svtvm(nq, facDelta, m_fields[n], 1, vel, 1, m_delta[n], 1);
    }

    // Update pressure (directly to outFields)
    if (m_Scale)
    {
        Array<OneD, NekDouble> wsp1(nq);
        Array<OneD, NekDouble> wsp2(ncoeffs);

        if (waveSpace)
        {
            Array<OneD, NekDouble> phys(pFields[dim]->GetTotPoints());
            pFields[dim]->HomogeneousBwdTrans(pFields[dim]->GetTotPoints(),
                                              pFields[dim]->GetPhys(), phys);

            // Interpolate phys-field by ScaleNumModes
            pFields[dim]->PhysInterp1DScaled(m_ScaleNumModes, phys, wsp1);

            m_pFieldsScaled[dim]->FwdTransLocalElmt(wsp1, wsp2);
            Vmath::Svtsvtp(m_outFields[dim].size(), facAvg, wsp2, 1, facOld,
                           m_outFields[dim], 1, m_outFields[dim], 1);
        }
        else
        {
            // Interpolate phys-field by ScaleNumModes
            pFields[dim]->PhysInterp1DScaled(m_ScaleNumModes,
                                             pFields[dim]->GetPhys(), wsp1);
            m_pFieldsScaled[dim]->FwdTransLocalElmt(wsp1, wsp2);
            Vmath::Svtsvtp(m_outFields[dim].size(), facAvg, wsp2, 1, facOld,
                           m_outFields[dim], 1, m_outFields[dim], 1);
        }
    }
    else
    {
        Vmath::Svtsvtp(m_outFields[dim].size(), facAvg,
                       pFields[dim]->GetCoeffs(), 1, facOld, m_outFields[dim],
                       1, m_outFields[dim], 1);
    }

    // Ignore Reynolds stress for first sample (its contribution is zero)
    if (m_numSamples == 1)
    {
        return;
    }

    // Calculate C_{n} = facOld * C_{n-1} + facStress * deltaI * deltaJ
    for (i = 0, n = dim + 1; i < dim; ++i)
    {
        for (j = i; j < dim; ++j, ++n)
        {
            Vmath::Vmul(nq, m_delta[i], 1, m_delta[j], 1, tmp, 1);
            Vmath::Svtsvtp(nq, facStress, tmp, 1, facOld, m_fields[n], 1,
                           m_fields[n], 1);
        }
    }
}

void FilterReynoldsStresses::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const NekDouble &time)
{
    size_t dim = pFields.size() - 1;

    m_fieldMetaData["NumberOfFieldDumps"] = std::to_string(m_numSamples);

    // Set wavespace to false, as calculations were performed in physical space
    bool waveSpace = pFields[0]->GetWaveSpace();
    pFields[0]->SetWaveSpace(false);

    // Forward transform and put into m_outFields (except pressure)
    for (size_t i = 0; i < m_fields.size(); ++i)
    {
        if (i != dim)
        {
            if (m_Scale)
            {
                m_pFieldsScaled[0]->FwdTransLocalElmt(m_fields[i],
                                                      m_outFields[i]);
            }
            else
            {
                pFields[0]->FwdTransLocalElmt(m_fields[i], m_outFields[i]);
            }
        }
    }

    // Restore waveSpace
    pFields[0]->SetWaveSpace(waveSpace);
}

NekDouble FilterReynoldsStresses::v_GetScale()
{
    if (m_movAvg)
    {
        return 1.0;
    }
    else
    {
        return 1.0 / m_numSamples;
    }
}

void FilterReynoldsStresses::v_OutputField(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, int dump)
{
    NekDouble scale = v_GetScale();
    for (int n = 0; n < m_outFields.size(); ++n)
    {
        Vmath::Smul(m_outFields[n].size(), scale, m_outFields[n], 1,
                    m_outFields[n], 1);
    }

    // Generating high order field
    if (m_Scale)
    {
        CreateFields(m_pFieldsScaled);
    }
    else
    {
        CreateFields(pFields);
    }

    // Determine new file name
    std::stringstream tmpOutname;
    std::string outname;
    int dot            = m_outputFile.find_last_of('.');
    std::string name   = m_outputFile.substr(0, dot);
    std::string ext    = m_outputFile.substr(dot, m_outputFile.length() - dot);
    std::string suffix = v_GetFileSuffix();

    if (dump == -1) // final dump
    {
        tmpOutname << name << suffix << ext;
    }
    else
    {
        tmpOutname << name << "_" << dump << suffix << ext;
    }
    outname = Filter::SetupOutput(ext, tmpOutname.str());
    m_modules[m_modules.size() - 1]->RegisterConfig("outfile", outname);

    // Run field process.
    for (int n = 0; n < SIZE_ModulePriority; ++n)
    {
        ModulePriority priority = static_cast<ModulePriority>(n);
        for (int i = 0; i < m_modules.size(); ++i)
        {
            if (m_modules[i]->GetModulePriority() == priority)
            {
                m_modules[i]->Process(m_vm);
            }
        }
    }

    // Empty m_f to save memory
    m_f->ClearField();

    if (dump != -1) // not final dump so rescale
    {
        for (int n = 0; n < m_outFields.size(); ++n)
        {
            Vmath::Smul(m_outFields[n].size(), 1.0 / scale, m_outFields[n], 1,
                        m_outFields[n], 1);
        }
    }
}

} // namespace Nektar::SolverUtils
