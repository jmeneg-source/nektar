///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingSyntheticEddy.cpp
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
// Description: Derived base class - Synthetic turbulence generation.
//              This code implements the Synthetic Eddy Method (SEM).
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingSyntheticEddy.h>
#include <ctime>
#include <fstream>
#include <iomanip>

namespace Nektar::SolverUtils
{

std::string ForcingSyntheticEddy::className =
    GetForcingFactory().RegisterCreatorFunction(
        "SyntheticTurbulence", ForcingSyntheticEddy::create,
        "Synthetic Eddy Turbulence Method");

ForcingSyntheticEddy::ForcingSyntheticEddy(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

/**
 * @brief Read input from xml file and initialise the class members.
 *        The main parameters are the characteristic lengths, Reynolds
 *        stresses and the synthetic eddy region (box of eddies).
 *
 * @param pFields           Pointer to fields.
 * @param pNumForcingField  Number of forcing fields.
 * @param pForce            Xml element describing the mapping.
 */
void ForcingSyntheticEddy::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const unsigned int &pNumForcingFields,
    [[maybe_unused]] const TiXmlElement *pForce)
{
    // Space dimension
    m_spacedim = pFields[0]->GetGraph()->GetMeshDimension();

    if (m_spacedim != 3)
    {
        NEKERROR(Nektar::ErrorUtil::efatal,
                 "Sythetic eddy method "
                 "only supports fully three-dimensional simulations");
    }

    // Get gamma parameter
    m_session->LoadParameter("Gamma", m_gamma, 1.4);

    const TiXmlElement *elmtInfTurb;

    // Reynolds Stresses
    elmtInfTurb = pForce->FirstChildElement("ReynoldsStresses");
    ASSERTL0(elmtInfTurb, "Requires ReynoldsStresses tag specifying function "
                          "name which prescribes the reynodls stresses.");

    std::string reyStressesFuncName = elmtInfTurb->GetText();
    ASSERTL0(m_session->DefinesFunction(reyStressesFuncName),
             "Function '" + reyStressesFuncName +
                 "' is not defined in the session.");

    // Reynolds stresses variables. Do not change the order of the variables.
    std::vector<std::string> reyStressesVars = {"r00", "r10", "r20",
                                                "r11", "r21", "r22"};

    for (size_t i = 0; i < reyStressesVars.size(); ++i)
    {
        std::string varStr = reyStressesVars[i];
        if (m_session->DefinesFunction(reyStressesFuncName, varStr))
        {
            m_R[i] = m_session->GetFunction(reyStressesFuncName, varStr);
        }
        else
        {
            NEKERROR(Nektar::ErrorUtil::efatal,
                     "Missing parameter '" + varStr +
                         "' in the FUNCTION NAME = " + reyStressesFuncName +
                         ".");
        }
    }

    // Characteristic length scales
    elmtInfTurb = pForce->FirstChildElement("CharLengthScales");
    ASSERTL0(elmtInfTurb, "Requires CharLengthScales tag specifying function "
                          "name which prescribes the characteristic "
                          "length scales.");

    std::string charLenScalesFuncName = elmtInfTurb->GetText();
    ASSERTL0(m_session->DefinesFunction(charLenScalesFuncName),
             "Function '" + charLenScalesFuncName +
                 "' is not defined in the session.");

    // Characteristic length scale variables
    // Do not change the order of the variables
    std::vector<std::string> lenScalesVars = {"l00", "l10", "l20", "l01", "l11",
                                              "l21", "l02", "l12", "l22"};
    LibUtilities::EquationSharedPtr clsAux;
    m_l = Array<OneD, NekDouble>(m_spacedim * m_spacedim, 0.0);

    for (size_t i = 0; i < lenScalesVars.size(); ++i)
    {
        std::string varStr = lenScalesVars[i];
        if (m_session->DefinesFunction(charLenScalesFuncName, varStr))
        {
            clsAux = m_session->GetFunction(charLenScalesFuncName, varStr);
            m_l[i] = (NekDouble)clsAux->Evaluate();
        }
        else
        {
            NEKERROR(Nektar::ErrorUtil::efatal,
                     "Missing parameter '" + varStr +
                         "' in the FUNCTION NAME = " + charLenScalesFuncName +
                         ".");
        }
    }

    // Read box of eddies parameters
    m_rc = Array<OneD, NekDouble>(m_spacedim, 0.0);
    // Array<OneD, NekDouble> m_lyz(m_spacedim - 1, 0.0);
    m_lyz       = Array<OneD, NekDouble>(m_spacedim - 1, 0.0);
    elmtInfTurb = pForce->FirstChildElement("BoxOfEddies");
    ASSERTL0(elmtInfTurb,
             "Unable to find BoxOfEddies tag. in SyntheticTurbulence forcing");

    if (elmtInfTurb)
    {
        std::stringstream boxStream;
        std::string boxStr = elmtInfTurb->GetText();
        boxStream.str(boxStr);
        size_t countVar = 0;
        for (size_t i = 0; i < (2 * m_spacedim - 1); ++i)
        {
            boxStream >> boxStr;
            if (i < m_spacedim)
            {
                m_rc[i] = std::stod(boxStr);
            }
            else
            {
                m_lyz[i - m_spacedim] = std::stod(boxStr);
            }
            countVar += 1;
        }
        if (countVar != (2 * m_spacedim - 1))
        {
            NEKERROR(Nektar::ErrorUtil::efatal,
                     "Missing parameter in the BoxOfEddies tag");
        }
    }

    // Read sigma
    elmtInfTurb = pForce->FirstChildElement("Sigma");
    ASSERTL0(elmtInfTurb,
             "Unable to find Sigma tag. in SyntheticTurbulence forcing");
    std::string sigmaStr = elmtInfTurb->GetText();
    m_sigma              = std::stod(sigmaStr);

    // Read bulk velocity
    elmtInfTurb = pForce->FirstChildElement("BulkVelocity");
    ASSERTL0(elmtInfTurb,
             "Unable to find BulkVelocity tag. in SyntheticTurbulence forcing");
    std::string bVelStr = elmtInfTurb->GetText();
    m_Ub                = std::stod(bVelStr);

    // Read flag to check if the run is a test case
    elmtInfTurb          = pForce->FirstChildElement("TestCase");
    const char *tcaseStr = (elmtInfTurb) ? elmtInfTurb->GetText() : "NoName";
    m_tCase              = (strcmp(tcaseStr, "ChanFlow3D") == 0) ? true : false;

    // Set Cholesky decomposition of the Reynolds Stresses in the domain
    SetCholeskyReyStresses(pFields);
    // Compute reference lengths
    ComputeRefLenghts();
    // Compute Xi maximum
    ComputeXiMax();
    // Set Number of Eddies
    SetNumberOfEddies();
    // Set mask
    SetBoxOfEddiesMask(pFields);
    // Set Forcing for each eddy
    InitialiseForcingEddy(pFields);
    // Check for test case
    if (!m_tCase)
    {
        // Compute initial location of the eddies in the box
        ComputeInitialRandomLocationOfEddies();
    }
    else
    {
        // Compute initial location of the eddies for the test case
        ComputeInitialLocationTestCase();
    }

    // Seed to generate random positions for the eddies
    srand(time(nullptr));

    // Initialise member from the base class
    m_Forcing = Array<OneD, Array<OneD, NekDouble>>(pFields.size());
    for (int i = 0; i < pFields.size(); ++i)
    {
        m_Forcing[i] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }

    // Initialise eddies ID vector
    for (size_t n = 0; n < m_N; ++n)
    {
        m_eddiesIDForcing.push_back(n);
    }
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
void ForcingSyntheticEddy::v_Apply(
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble &time)
{
}

/**
 * @brief Apply forcing term if an eddy left the box of eddies and
 *        update the eddies positions.
 *
 * @param pFields   Pointer to fields.
 * @param inarray   Given fields. The fields are in in physical space.
 * @param outarray  Calculated solution after forcing term being applied
 *                  in coefficient space.
 * @param time      time.
 */
void ForcingSyntheticEddy::v_ApplyCoeff(
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble &time)
{
}

/**
 * @brief Compute characteristic convective turbulent time.
 *
 * @param pFields  Pointer to fields.
 */
Array<OneD, Array<OneD, NekDouble>> ForcingSyntheticEddy::
    ComputeCharConvTurbTime(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Characteristic convective turbulent time
    Array<OneD, Array<OneD, NekDouble>> convTurbTime(m_spacedim);

    for (size_t k = 0; k < m_spacedim; ++k)
    {
        convTurbTime[k] = Array<OneD, NekDouble>(nqTot, 0.0);
        for (size_t i = 0; i < nqTot; ++i)
        {
            NekDouble convTurbLength = m_xiMaxMin * m_lref[0];
            // 3*k because of the structure of the m_l parameter
            // to obtain lxk.
            if ((m_l[3 * k] > m_xiMaxMin * m_lref[0]) && (m_mask[i]))
            {
                convTurbLength = m_l[3 * k];
            }
            convTurbTime[k][i] = convTurbLength / m_Ub;
        }
    }

    return convTurbTime;
}

/**
 * @brief Compute smoothing factor to avoid strong variations
 *        of the source term across the domain.
 *
 * @param pFields       Pointer to fields.
 * @param convTurbTime  Convective turbulent time.
 */
Array<OneD, Array<OneD, NekDouble>> ForcingSyntheticEddy::
    ComputeSmoothingFactor(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, Array<OneD, NekDouble>> convTurbTime)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Total number of quadrature points of each element
    int nqe;
    // Characteristic convective turbulent time
    Array<OneD, Array<OneD, NekDouble>> smoothFac(m_spacedim);
    // Counter
    int count = 0;
    // module
    NekDouble mod;
    // Create Array with zeros
    for (size_t i = 0; i < m_spacedim; ++i)
    {
        smoothFac[i] = Array<OneD, NekDouble>(nqTot, 0.0);
    }

    for (size_t e = 0; e < nElmts; ++e)
    {
        nqe = pFields[0]->GetExp(e)->GetTotPoints();

        Array<OneD, NekDouble> coords0(nqe, 0.0);
        Array<OneD, NekDouble> coords1(nqe, 0.0);
        Array<OneD, NekDouble> coords2(nqe, 0.0);

        pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

        for (size_t i = 0; i < nqe; ++i)
        {
            if (m_mask[count + i])
            {
                mod = (coords0[i] - m_rc[0] + m_lref[0]) *
                      (coords0[i] - m_rc[0] + m_lref[0]);

                smoothFac[0][count + i] =
                    exp((-0.5 * M_PI * mod) /
                        (convTurbTime[0][count + i] *
                         convTurbTime[0][count + i] * m_Ub * m_Ub));
                smoothFac[1][count + i] =
                    exp((-0.5 * M_PI * mod) /
                        (convTurbTime[1][count + i] *
                         convTurbTime[1][count + i] * m_Ub * m_Ub));
                smoothFac[2][count + i] =
                    exp((-0.5 * M_PI * mod) /
                        (convTurbTime[2][count + i] *
                         convTurbTime[2][count + i] * m_Ub * m_Ub));
            }
        }
        count += nqe;
    }

    return smoothFac;
}

/**
 * @brief Calculate velocity fluctuation for the source term
 *
 * @param pFields           Pointer to fields.
 * @param stochasticSignal  Stochastic signal.
 * @return                  Velocity fluctuation.
 */
Array<OneD, Array<OneD, NekDouble>> ForcingSyntheticEddy::
    ComputeVelocityFluctuation(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, Array<OneD, NekDouble>> stochasticSignal)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Velocity fluctuation vector
    Array<OneD, Array<OneD, NekDouble>> velFluc(m_N);
    // Control loop for the m_Cholesky (Cholesky decomposition matrix)
    int l;

    for (auto &n : m_eddiesIDForcing)
    {
        velFluc[n] = Array<OneD, NekDouble>(nqTot * m_spacedim, 0.0);
        for (size_t k = 0; k < m_spacedim; ++k)
        {
            for (size_t j = 0; j < k + 1; ++j)
            {
                for (size_t i = 0; i < nqTot; ++i)
                {
                    if (m_mask[i])
                    {
                        l = k + j * (2 * m_spacedim - j - 1) * 0.5;
                        velFluc[n][k * nqTot + i] +=
                            m_Cholesky[i][l] *
                            stochasticSignal[n][j * nqTot + i];
                    }
                }
            }
        }
    }

    return velFluc;
}

/**
 * @brief Compute stochastic signal.
 *
 * @param pFields  Pointer to fields.
 * @return         Stochastic signal.
 */
Array<OneD, Array<OneD, NekDouble>> ForcingSyntheticEddy::
    ComputeStochasticSignal(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Total number of quadrature points of each element
    int nqe;
    // Stochastic Signal vector
    Array<OneD, Array<OneD, NekDouble>> stochasticSignal(m_N);
    // Random numbers: -1 and 1
    Array<OneD, Array<OneD, int>> epsilonSign;

    // Generate only for the new eddies after the first time step
    epsilonSign = GenerateRandomOneOrMinusOne();

    // Calculate the stochastic signal for the eddies
    for (auto &n : m_eddiesIDForcing)
    {
        stochasticSignal[n] = Array<OneD, NekDouble>(nqTot * m_spacedim, 0.0);

        // Evaluate the function at interpolation points for each component
        for (size_t j = 0; j < m_spacedim; ++j)
        {
            // Count the number of quadrature points
            int nqeCount = 0;

            for (size_t e = 0; e < nElmts; ++e)
            {
                nqe = pFields[0]->GetExp(e)->GetTotPoints();

                Array<OneD, NekDouble> coords0(nqe, 0.0);
                Array<OneD, NekDouble> coords1(nqe, 0.0);
                Array<OneD, NekDouble> coords2(nqe, 0.0);

                pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

                // i: degrees of freedom, j: direction, n: eddies
                for (size_t i = 0; i < nqe; ++i)
                {
                    if (m_mask[nqeCount + i])
                    {
                        stochasticSignal[n][j * nqTot + nqeCount + i] =
                            epsilonSign[j][n] *
                            ComputeGaussian((coords0[i] - m_eddyPos[n][0]) /
                                                m_lref[0],
                                            m_xiMax[j * m_spacedim + 0],
                                            ComputeConstantC(0, j)) *
                            ComputeGaussian((coords1[i] - m_eddyPos[n][1]) /
                                                m_lref[1],
                                            m_xiMax[j * m_spacedim + 1],
                                            ComputeConstantC(1, j)) *
                            ComputeGaussian((coords2[i] - m_eddyPos[n][2]) /
                                                m_lref[2],
                                            m_xiMax[j * m_spacedim + 2],
                                            ComputeConstantC(2, j));
                    }
                }
                nqeCount += nqe;
            }
        }
    }

    return stochasticSignal;
}

/**
 * @brief Update the position of the eddies for every time step.
 */
void ForcingSyntheticEddy::UpdateEddiesPositions()
{
    NekDouble dt = m_session->GetParameter("TimeStep");

    for (size_t n = 0; n < m_N; ++n)
    {
        // Check if any eddy went through the outlet plane (box)
        if ((m_eddyPos[n][0] + m_Ub * dt) < (m_rc[0] + m_lref[0]))
        {
            m_eddyPos[n][0] = m_eddyPos[n][0] + m_Ub * dt;
        }
        else
        {
            // Generate a new one in the inlet plane
            m_eddyPos[n][0] = m_rc[0] - m_lref[0];
            // Check if test case
            if (!m_tCase)
            {
                m_eddyPos[n][1] =
                    (m_rc[1] - 0.5 * m_lyz[0]) +
                    (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * (m_lyz[0]);
                m_eddyPos[n][2] =
                    (m_rc[2] - 0.5 * m_lyz[1] + 0.5 * m_lref[2]) +
                    (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * (m_lyz[1]);
            }
            else
            {
                // same place (center) for the test case
                m_eddyPos[n][1] = m_rc[1];
                m_eddyPos[n][2] = m_rc[2];
            }
            m_eddiesIDForcing.push_back(n);
            m_calcForcing = true;
        }
    }
}

/**
 * @brief Remove eddy from forcing term
 *
 * @param pFields  Pointer to fields.
 */
void ForcingSyntheticEddy::RemoveEddiesFromForcing(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of Variables
    int nVars = pFields.size();

    for (auto &n : m_eddiesIDForcing)
    {
        for (size_t j = 0; j < nVars; ++j)
        {
            for (size_t i = 0; i < nqTot; ++i)
            {
                m_Forcing[j][i] -= m_ForcingEddy[n][j][i];
            }
        }
    }
}

/**
 * @brief Calculate distribution of eddies in the box.
 */
void ForcingSyntheticEddy::ComputeInitialRandomLocationOfEddies()
{
    m_eddyPos = Array<OneD, Array<OneD, NekDouble>>(m_N);

    for (size_t n = 0; n < m_N; ++n)
    {
        m_eddyPos[n] = Array<OneD, NekDouble>(m_spacedim);
        // Generate randomly eddies inside the box
        m_eddyPos[n][0] =
            (m_rc[0] - m_lref[0]) +
            (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * 2 * m_lref[0];
        m_eddyPos[n][1] =
            (m_rc[1] - 0.5 * m_lyz[0]) +
            (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * m_lyz[0];
        m_eddyPos[n][2] =
            (m_rc[2] - 0.5 * m_lyz[1]) +
            (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * m_lyz[1];
    }
}

/**
 * @brief Compute standard Gaussian with zero mean.
 *
 * @param coord  Coordianate.
 * @return       Gaussian value for the coordinate.
 */
NekDouble ForcingSyntheticEddy::ComputeGaussian(NekDouble coord,
                                                NekDouble xiMaxVal,
                                                NekDouble constC)
{
    NekDouble epsilon = 1e-6;
    if (abs(coord) <= (xiMaxVal + epsilon))
    {
        return ((1.0 / (m_sigma * sqrt(2.0 * M_PI * constC))) *
                exp(-0.5 * (coord / m_sigma) * (coord / m_sigma)));
    }
    else
    {
        return 0.0;
    }
}

/**
 * @brief Compute constant C for the gaussian funcion.
 *
 * @param row  index for the rows of the matrix.
 * @param col  index for the columns of the matrix.
 * @return     Value of C.
 */
NekDouble ForcingSyntheticEddy::ComputeConstantC(int row, int col)
{
    NekDouble sizeLenScale = m_xiMax[col * m_spacedim + row];

    // Integration
    NekDouble sum  = 0.0;
    NekDouble step = 0.025;
    NekDouble xi0  = -1;
    NekDouble xif  = 1;
    while (xi0 < xif)
    {
        sum += (ComputeGaussian(xi0 + step, sizeLenScale) *
                    ComputeGaussian(xi0 + step, sizeLenScale) +
                ComputeGaussian(xi0, sizeLenScale) *
                    ComputeGaussian(xi0, sizeLenScale));
        xi0 += step;
    }

    return (0.5 * 0.5 * step * sum);
}

/**
 * @brief Generate random 1 or -1 values to be use to compute
 *        the stochastic signal.
 * @return ramdom 1 or -1 values.
 */
Array<OneD, Array<OneD, int>> ForcingSyntheticEddy::
    GenerateRandomOneOrMinusOne()
{
    Array<OneD, Array<OneD, int>> epsilonSign(m_spacedim);

    // j: component of the stochastic signal
    // n: eddy
    for (size_t j = 0; j < m_spacedim; ++j)
    {
        epsilonSign[j] = Array<OneD, int>(m_N, 0.0);
        for (auto &n : m_eddiesIDForcing)
        {
            // Convert to -1 or to 1
            // Check if test case
            if (!m_tCase)
            {
                epsilonSign[j][n] =
                    ((NekSingle(std::rand()) / NekSingle(RAND_MAX)) <= 0.5) ? -1
                                                                            : 1;
            }
            else
            {
                // Positive values only for the test case
                epsilonSign[j][n] = 1;
            }
        }
    }

    return epsilonSign;
}

/**
 * @brief Set box of eddies mask to be used to seprate the
 *        degrees of freedom (coordinates) inside and outside
 *        the box of eddies.
 *
 * @param pFields  Pointer to fields.
 */
void ForcingSyntheticEddy::SetBoxOfEddiesMask(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Total number of quadrature points of each element
    int nqe;
    // Mask
    m_mask = Array<OneD, int>(nqTot, 0); // 0 for ouside, 1 for inside
    // Counter
    int count = 0;

    for (size_t e = 0; e < nElmts; ++e)
    {
        nqe = pFields[0]->GetExp(e)->GetTotPoints();

        Array<OneD, NekDouble> coords0(nqe, 0.0);
        Array<OneD, NekDouble> coords1(nqe, 0.0);
        Array<OneD, NekDouble> coords2(nqe, 0.0);

        pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

        for (size_t i = 0; i < nqe; ++i)
        {
            if (InsideBoxOfEddies(coords0[i], coords1[i], coords2[i]))
            {
                // 0 for outside, 1 for inside
                m_mask[count + i] = 1;
            }
        }
        count += nqe;
    }
}

/**
 * @brief Initialise Forcing term for each eddy.
 */
void ForcingSyntheticEddy::InitialiseForcingEddy(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of Variables
    int nVars     = pFields.size();
    m_ForcingEddy = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_N);

    for (size_t i = 0; i < m_N; ++i)
    {
        m_ForcingEddy[i] = Array<OneD, Array<OneD, NekDouble>>(nVars);
        for (size_t j = 0; j < nVars; ++j)
        {
            m_ForcingEddy[i][j] = Array<OneD, NekDouble>(nqTot);
            for (size_t k = 0; k < nqTot; ++k)
            {
                m_ForcingEddy[i][j][k] = 0.0;
            }
        }
    }
}

/**
 * @brief Check it point is inside the box of eddies.
 *
 * @param coord0  coordinate in the x-direction.
 * @param coord1  coordinate in the y-direction.
 * @param coord2  coordinate in the z-direction.
 * @return        true or false
 */
bool ForcingSyntheticEddy::InsideBoxOfEddies(NekDouble coord0, NekDouble coord1,
                                             NekDouble coord2)
{
    NekDouble tol = 1e-6;
    if ((coord0 >= (m_rc[0] - m_lref[0] - m_lref[0])) &&
        (coord0 <= (m_rc[0] + m_lref[0] + tol)) &&
        (coord1 >= (m_rc[1] - 0.5 * m_lyz[0] - tol)) &&
        (coord1 <= (m_rc[1] + 0.5 * m_lyz[0] + tol)) &&
        (coord2 >= (m_rc[2] - 0.5 * m_lyz[1] - tol)) &&
        (coord2 <= (m_rc[2] + 0.5 * m_lyz[1] + tol)))
    {
        return true;
    }

    return false;
}

/**
 * @brief Calculates the reference lenghts
 */
void ForcingSyntheticEddy::ComputeRefLenghts()
{
    m_lref    = Array<OneD, NekDouble>(m_spacedim, 0.0);
    m_lref[0] = m_l[0];
    m_lref[1] = m_l[1];
    m_lref[2] = m_l[2];

    // The l_{x}^{ref}, l_{y}^{ref} and l_{z}^{ref}
    // are the maximum among the velocity components
    // over the domain.
    for (size_t i = 0; i < m_spacedim; ++i)
    {
        for (size_t j = 1; j < m_spacedim; ++j)
        {
            if (m_l[m_spacedim * j + i] > m_lref[i])
            {
                m_lref[i] = m_l[m_spacedim * j + i];
            }
        }
    }
}

/**
 * @brief Calculates the \f$\xi_{max}\f$.
 */
void ForcingSyntheticEddy::ComputeXiMax()
{
    NekDouble value;
    m_xiMax = Array<OneD, NekDouble>(m_spacedim * m_spacedim, 0.0);

    for (size_t i = 0; i < m_spacedim; i++)
    {
        for (size_t j = 0; j < m_spacedim; j++)
        {
            value = (m_l[m_spacedim * j + i] / m_lref[i]);
            if (value > m_xiMaxMin)
            {
                m_xiMax[m_spacedim * j + i] = value;
            }
            else
            {
                m_xiMax[m_spacedim * j + i] = m_xiMaxMin;
            }
        }
    }
}

/**
 * @brief Calculates the Cholesky decomposition of the Reynolds Stresses
 *        in each degree of freedom of the mesh.
 *
 * @param pFields  Pointer to fields.
 */
void ForcingSyntheticEddy::SetCholeskyReyStresses(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Total number of quadrature points of each element
    int nqe;
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Count the degrees of freedom
    int nqeCount = 0;
    // Block diagonal size
    int diagSize = m_spacedim * (m_spacedim + 1) * 0.5;
    // Cholesky decomposition matrix for the synthetic eddy region (box)
    m_Cholesky = Array<OneD, Array<OneD, NekDouble>>(nqTot);

    for (size_t e = 0; e < nElmts; ++e)
    {
        nqe = pFields[0]->GetExp(e)->GetTotPoints();

        Array<OneD, NekDouble> coords0(nqe, 0.0);
        Array<OneD, NekDouble> coords1(nqe, 0.0);
        Array<OneD, NekDouble> coords2(nqe, 0.0);

        // Coordinates (for each degree of freedom) for the element k.
        pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

        // Evaluate Cholesky decomposition
        for (size_t i = 0; i < nqe; ++i)
        {
            int l = 0;
            // Size of Cholesky decomposition matrix - aux vector
            Array<OneD, NekDouble> A(diagSize, 0.0);

            while (l < diagSize)
            {
                // Variable to compute the Cholesky decomposition for each
                // degree of freedom
                A[l] = m_R[l]->Evaluate(coords0[i], coords1[i], coords2[i]);
                l++;
            }
            int info = 0;
            Lapack::Dpptrf('L', m_spacedim, A.data(), info);
            if (info < 0)
            {
                std::string message =
                    "ERROR: The " + std::to_string(-info) +
                    "th parameter had an illegal parameter for dpptrf";
                NEKERROR(ErrorUtil::efatal, message.c_str());
            }
            m_Cholesky[nqeCount + i] = Array<OneD, NekDouble>(diagSize);
            for (size_t l = 0; l < diagSize; ++l)
            {
                m_Cholesky[nqeCount + i][l] = A[l];
            }
        }
        nqeCount += nqe;
    }
}

/**
 * @brief Calculate the number of eddies that are going to be
 *        injected in the synthetic eddy region (box).
 */
void ForcingSyntheticEddy::SetNumberOfEddies()
{
    m_N = int((m_lyz[0] * m_lyz[1]) /
              (4 * m_lref[m_spacedim - 2] * m_lref[m_spacedim - 1])) +
          1;
}

/**
 * @brief Place eddies in specific locations in the test case
 *        geometry for consistency and comparison.
 *
 *        This function was design for a three-dimensional
 *        channel flow test case (ChanFlow3d_infTurb).
 *        It is only called for testing purposes (unit test).
 */
void ForcingSyntheticEddy::ComputeInitialLocationTestCase()
{
    m_N       = 3; // Redefine number of eddies
    m_eddyPos = Array<OneD, Array<OneD, NekDouble>>(m_N);

    // First eddy (center)
    m_eddyPos[0]    = Array<OneD, NekDouble>(m_spacedim);
    m_eddyPos[0][0] = (m_rc[0] - m_lref[0]) + 0.6 * 2 * m_lref[0];
    m_eddyPos[0][1] = m_rc[1];
    m_eddyPos[0][2] = m_rc[2];

    // Second eddy (top)
    m_eddyPos[1]    = Array<OneD, NekDouble>(m_spacedim);
    m_eddyPos[1][0] = (m_rc[0] - m_lref[0]) + 0.7 * 2 * m_lref[0];
    m_eddyPos[1][1] = (m_rc[1] - 0.5 * m_lyz[0]) + 0.9 * m_lyz[0];
    m_eddyPos[1][2] = m_rc[2];

    // Third eddy (bottom)
    m_eddyPos[2]    = Array<OneD, NekDouble>(m_spacedim);
    m_eddyPos[2][0] = (m_rc[0] - m_lref[0]) + 0.8 * 2 * m_lref[0];
    m_eddyPos[2][1] = (m_rc[1] - 0.5 * m_lyz[0]) + 0.1 * m_lyz[0];
    m_eddyPos[2][2] = m_rc[2];
}

} // namespace Nektar::SolverUtils
