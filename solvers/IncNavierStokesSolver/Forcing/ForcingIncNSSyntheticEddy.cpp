///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingIncNSSyntheticEddy.cpp
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
//              Incompressible solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingIncNSSyntheticEddy.h>

namespace Nektar::SolverUtils
{

std::string ForcingIncNSSyntheticEddy::className =
    GetForcingFactory().RegisterCreatorFunction(
        "IncNSSyntheticTurbulence", ForcingIncNSSyntheticEddy::create,
        "Inc NS Synthetic Eddy Turbulence Forcing (Generation)");

ForcingIncNSSyntheticEddy::ForcingIncNSSyntheticEddy(
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
void ForcingIncNSSyntheticEddy::v_Apply(
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

    // Incompressible version
    // Only velocities u,v,w: nVars - 1
    for (size_t i = 0; i < (nVars - 1); ++i)
    {
        Vmath::Vadd(nqTot, m_Forcing[i], 1, outarray[i], 1, outarray[i], 1);
    }

    // Update eddies position inside the box.
    UpdateEddiesPositions();
}

/**
 * @brief Calculate forcing term.
 *
 * @param pFfields  Pointer to fields.
 */
void ForcingIncNSSyntheticEddy::CalculateForcing(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();

    // Compute Stochastic Signal
    Array<OneD, Array<OneD, NekDouble>> stochasticSignal;
    stochasticSignal = ComputeStochasticSignal(pFields);

    // Compute velocity flsuctuation
    Array<OneD, Array<OneD, NekDouble>> velFluc;
    velFluc = ComputeVelocityFluctuation(pFields, stochasticSignal);

    // Compute characteristic convective turbulent time
    Array<OneD, Array<OneD, NekDouble>> convTurbTime;
    convTurbTime = ComputeCharConvTurbTime(pFields);

    // Compute smoothing factor
    Array<OneD, Array<OneD, NekDouble>> smoothFac;
    smoothFac = ComputeSmoothingFactor(pFields, convTurbTime);

    // Check if eddies left the box. Note that the member m_eddiesIDForcing
    // is populate with the eddies that left the box. If any left it is going
    // to be empty.
    if (!m_eddiesIDForcing.empty())
    {
        // Forcing must stop applying for eddies that left the box.
        RemoveEddiesFromForcing(pFields);

        // Update Forcing term which are going to be applied until
        // the eddy leave the domain.
        // Select the eddies to apply the forcing. Superposition.
        for (auto &n : m_eddiesIDForcing)
        {
            //  velocity term
            for (size_t j = 0; j < m_spacedim; ++j)
            {
                for (size_t i = 0; i < nqTot; ++i)
                {
                    if (m_mask[i])
                    {
                        m_ForcingEddy[n][j][i] +=
                            ((velFluc[n][j * nqTot + i] * smoothFac[j][i]) /
                             convTurbTime[j][i]);
                        // Update forcing
                        m_Forcing[j][i] += m_ForcingEddy[n][j][i];
                    }
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

} // namespace Nektar::SolverUtils
