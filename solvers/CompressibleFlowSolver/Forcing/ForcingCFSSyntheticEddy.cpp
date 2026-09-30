///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingCFSSyntheticEddy.cpp
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
// Description: Derived base class - Synthetic turbulence forcing for the
//              Compressible solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/Forcing/ForcingCFSSyntheticEddy.h>

namespace Nektar::SolverUtils
{

std::string ForcingCFSSyntheticEddy::className =
    GetForcingFactory().RegisterCreatorFunction(
        "CFSSyntheticTurbulence", ForcingCFSSyntheticEddy::create,
        "Compressible Synthetic Eddy Turbulence Forcing (Generation)");

ForcingCFSSyntheticEddy::ForcingCFSSyntheticEddy(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation), ForcingSyntheticEddy(pSession, pEquation)
{
}

/**
 * @brief Apply forcing term if an eddy left the box of eddies and
 *        update the eddies positions.
 *
 * @param pFields   Pointer to fields.
 * @param inarray   Given fields. The fields are in in physical space.
 * @param outarray  Calculated solution after forcing term being applied
 *                  in physical space.
 * @param time      time.
 */
void ForcingCFSSyntheticEddy::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble &time)
{
    // Number of Variables
    int nVars = pFields.size();
    // Total number of coefficients
    unsigned int nqTot = pFields[0]->GetTotPoints();

    // Only calculates the forcing in the first time step and when an eddy
    // leaves the synthetic eddy region (box).
    if (m_calcForcing)
    {
        CalculateForcing(pFields);
        m_calcForcing = false;
    }

    for (size_t i = 0; i < nVars; ++i)
    {
        Vmath::Vadd(nqTot, m_Forcing[i], 1, outarray[i], 1, outarray[i], 1);
    }

    // Update eddies position inside the box.
    UpdateEddiesPositions();
}

/**
 * @brief Apply forcing term if an eddy left the box of eddies and
 *        update the eddies positions.
 *
 * @param pFields   Pointer to fields.
 * @param inarray   Given fields. The fields are in in physical space.
 * @param outarray  Calculated solution after forcing term being applied
 *                  in coeficient space.
 * @param time      time.
 */
void ForcingCFSSyntheticEddy::v_ApplyCoeff(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    // Number of Variables
    int nVars = fields.size();
    // Total number of coefficients
    unsigned int nCoeff = fields[0]->GetNcoeffs();

    // Only calculates the forcing in the first time step and when an eddy
    // leaves the box
    if (m_calcForcing)
    {
        CalculateForcing(fields);
        m_calcForcing = false;
    }

    Array<OneD, NekDouble> forcingCoeff(nCoeff);
    for (size_t i = 0; i < nVars; ++i)
    {
        fields[i]->FwdTrans(m_Forcing[i], forcingCoeff);
        Vmath::Vadd(nCoeff, forcingCoeff, 1, outarray[i], 1, outarray[i], 1);
    }

    // Check for update: implicit solver
    if (m_currTime != time)
    {
        UpdateEddiesPositions();
        m_currTime = time;
    }
}

/**
 * @brief Calculate forcing term.
 *
 * @param pFfields  Pointer to fields.
 */
void ForcingCFSSyntheticEddy::CalculateForcing(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of Variables
    int nVars = pFields.size();
    // Define parameters
    Array<OneD, NekDouble> rho(nqTot), temperature(nqTot), pressure(nqTot);
    // Velocity fluctuation module
    NekDouble velFlucMod;
    // Variable converter
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(m_session,
                                                                    m_spacedim);

    // Compute Stochastic Signal
    Array<OneD, Array<OneD, NekDouble>> stochasticSignal;
    stochasticSignal = ComputeStochasticSignal(pFields);

    // Compute rho and mach mean inside the box of eddies
    std::pair<NekDouble, NekDouble> rhoMachMean;
    rhoMachMean = ComputeRhoMachMean(pFields);

    // Compute velocity fluctuation
    Array<OneD, Array<OneD, NekDouble>> velFluc;
    velFluc = ComputeVelocityFluctuation(pFields, stochasticSignal);

    // Compute density fluctuation
    Array<OneD, Array<OneD, NekDouble>> rhoFluc;
    rhoFluc = ComputeDensityFluctuation(pFields, velFluc, rhoMachMean);

    // Compute characteristic convective turbulent time
    Array<OneD, Array<OneD, NekDouble>> convTurbTime;
    convTurbTime = ComputeCharConvTurbTime(pFields);

    // Compute smoothing factor
    Array<OneD, Array<OneD, NekDouble>> smoothFac;
    smoothFac = ComputeSmoothingFactor(pFields, convTurbTime);

    // Get physical data
    Array<OneD, Array<OneD, NekDouble>> physFields(nVars);
    for (size_t i = 0; i < nVars; ++i)
    {
        physFields[i] = pFields[i]->GetPhys();
    }

    // Calculate parameters
    m_varConv->GetTemperature(physFields, temperature);
    m_varConv->GetPressure(physFields, pressure);
    m_varConv->GetRhoFromPT(pressure, temperature, rho);

    // Check if eddies left the box (reinjected). Note that the member
    // m_eddiesIDForcing is populate with the eddies that left the box.
    // If any left it is going to be empty.
    if (!m_eddiesIDForcing.empty())
    {
        // Forcing must stop applying for eddies that left the box
        RemoveEddiesFromForcing(pFields);

        // Update Forcing term which are going to be applied until
        // the eddy leave the leave the domain
        // Select the eddies to apply the forcing. Superposition.
        for (auto &n : m_eddiesIDForcing)
        {
            for (size_t i = 0; i < nqTot; ++i)
            {
                if (m_mask[i])
                {
                    // density term
                    m_ForcingEddy[n][0][i] =
                        (rhoFluc[n][i] * smoothFac[0][i]) / convTurbTime[0][i];
                    // Update forcing
                    m_Forcing[0][i] += m_ForcingEddy[n][0][i];
                    //  velocity term
                    for (size_t j = 0; j < m_spacedim; ++j)
                    {
                        m_ForcingEddy[n][j + 1][i] =
                            (rho[i] * velFluc[n][j * nqTot + i] *
                             smoothFac[j][i]) /
                            convTurbTime[j][i];
                        // Update forcing
                        m_Forcing[j + 1][i] += m_ForcingEddy[n][j + 1][i];
                    }
                    // energy term
                    velFlucMod =
                        velFluc[n][i] * velFluc[n][i] +
                        velFluc[n][nqTot + i] * velFluc[n][nqTot + i] +
                        velFluc[n][2 * nqTot + i] * velFluc[n][2 * nqTot + i];

                    m_ForcingEddy[n][nVars - 1][i] =
                        (0.5 * rho[i] * velFlucMod * smoothFac[0][i]) /
                        convTurbTime[0][i];
                    // Update forcing
                    m_Forcing[nVars - 1][i] += m_ForcingEddy[n][nVars - 1][i];
                }
            }
        }
        // delete eddies
        m_eddiesIDForcing.erase(m_eddiesIDForcing.begin(),
                                m_eddiesIDForcing.end());
    }
    else
    {
        NEKERROR(ErrorUtil::efatal, "Failed: Eddies ID vector is empty.");
    }
}

/**
 * @brief Compute density fluctuation for the source term
 *
 * @param pFfields  Pointer to fields.
 */
Array<OneD, Array<OneD, NekDouble>> ForcingCFSSyntheticEddy::
    ComputeDensityFluctuation(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &velFluc,
        std::pair<NekDouble, NekDouble> rhoMachMean)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Density fluctuation
    Array<OneD, Array<OneD, NekDouble>> rhoFluc(m_N);

    for (auto &n : m_eddiesIDForcing)
    {
        rhoFluc[n] = Array<OneD, NekDouble>(nqTot, 0.0);

        for (size_t i = 0; i < nqTot; ++i)
        {
            if (m_mask[i])
            {
                // only need velocity fluctation in the x-direction
                rhoFluc[n][i] = rhoMachMean.first * (m_gamma - 1) *
                                rhoMachMean.second * rhoMachMean.second *
                                (velFluc[n][i] / m_Ub);
            }
        }
    }

    return rhoFluc;
}

/**
 * @brief Calculate the density and mach number mean
 *        inside the box of eddies.
 *
 * @param  pFfields  Pointer to fields.
 * @return           Pair with density and mach nember means.
 */
std::pair<NekDouble, NekDouble> ForcingCFSSyntheticEddy::ComputeRhoMachMean(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of Variables
    int nVars = pFields.size();
    // <rho> and <mach>^{2}
    NekDouble rhoMean = 0.0, machMean = 0.0;
    // Counter for the mean
    int count = 0;

    // Get physical field
    Array<OneD, Array<OneD, NekDouble>> physFields(nVars);
    for (size_t i = 0; i < nVars; ++i)
    {
        physFields[i] = pFields[i]->GetPhys();
    }

    // Define parameters
    Array<OneD, NekDouble> soundSpeed(nqTot), mach(nqTot), rho(nqTot),
        temperature(nqTot), pressure(nqTot);
    // Calculate parameters

    m_varConv->GetTemperature(physFields, temperature);
    m_varConv->GetPressure(physFields, pressure);
    m_varConv->GetRhoFromPT(pressure, temperature, rho);
    m_varConv->GetSoundSpeed(physFields, soundSpeed);
    m_varConv->GetMach(physFields, soundSpeed, mach);

    // Sum
    for (size_t i = 0; i < nqTot; ++i)
    {
        if (m_mask[i])
        {
            rhoMean += rho[i];
            machMean += mach[i];
            count += 1;
        }
    }

    // Density mean
    rhoMean = rhoMean / count;
    // Mach number mean
    machMean = machMean / count;

    return std::make_pair(rhoMean, machMean);
}

} // namespace Nektar::SolverUtils
