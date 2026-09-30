///////////////////////////////////////////////////////////////////////////////
//
// File NoSchemeTimeIntegrationScheme.cpp
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
// Description: implementation of time integration scheme NoScheme class
// This integration scheme does nothing
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/TimeIntegration/NoSchemeTimeIntegrationScheme.h>

namespace Nektar::LibUtilities
{

// Access Methods
std::string NoSchemeTimeIntegrationScheme::v_GetName() const
{
    return std::string("NoIntegration");
}

std::string NoSchemeTimeIntegrationScheme::v_GetVariant() const
{
    return std::string("");
}

size_t NoSchemeTimeIntegrationScheme::v_GetOrder() const
{
    return 0;
}

std::vector<NekDouble> NoSchemeTimeIntegrationScheme::v_GetFreeParams() const
{
    return std::vector<NekDouble>(0);
}

NekDouble NoSchemeTimeIntegrationScheme::v_GetTimeStability() const
{
    return 0.;
}

TimeIntegrationSchemeType NoSchemeTimeIntegrationScheme::
    v_GetIntegrationSchemeType() const
{
    return eNoTimeIntegrationSchemeType;
}

size_t NoSchemeTimeIntegrationScheme::v_GetNumIntegrationPhases() const
{
    return 0;
}

/**
 * @brief Worker method to initialize the integration scheme.
 */
void NoSchemeTimeIntegrationScheme::v_InitializeScheme(
    [[maybe_unused]] const NekDouble deltaT, ConstDoubleArray &y_0,
    [[maybe_unused]] const NekDouble time,
    [[maybe_unused]] const TimeIntegrationSchemeOperators &op)
{
    m_doubleArray = y_0;
}

/**
 * @brief Worker method that actually does the time integration.
 */
ConstDoubleArray &NoSchemeTimeIntegrationScheme::v_TimeIntegrate(
    [[maybe_unused]] const size_t timestep,
    [[maybe_unused]] const NekDouble delta_t)
{
    return m_doubleArray;
}

TripleArray &NoSchemeTimeIntegrationScheme::v_UpdateSolutionVector()
{
    return NullNekDoubleTensorOfArray3D;
}

/**
 * @brief Worker method to print details on the integration scheme
 */
void NoSchemeTimeIntegrationScheme::v_print(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl;
}

void NoSchemeTimeIntegrationScheme::v_printFull(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl;
}

} // namespace Nektar::LibUtilities
