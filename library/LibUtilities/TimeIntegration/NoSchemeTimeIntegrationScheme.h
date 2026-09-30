///////////////////////////////////////////////////////////////////////////////
//
// File: NoSchemeTimeIntegrationScheme.h
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
// Description: Header file of time integration scheme NoScheme class
// This integration scheme does nothing
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NOSCHEME_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME
#define NEKTAR_LIB_UTILITIES_NOSCHEME_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar::LibUtilities
{

/**
 * @brief Base class for No time integration schemes.
 */
class NoSchemeTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    LUE NoSchemeTimeIntegrationScheme(std::string variant, unsigned int order,
                                      std::vector<NekDouble> freeParams)
        : TimeIntegrationScheme(variant, order, freeParams)
    {
    }

    ~NoSchemeTimeIntegrationScheme() override = default;

    static NoTimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order,
        std::vector<NekDouble> freeParams)
    {
        NoTimeIntegrationSchemeSharedPtr p =
            MemoryManager<NoSchemeTimeIntegrationScheme>::AllocateSharedPtr(
                variant, order, freeParams);
        return p;
    }

    // Access methods
    LUE std::string v_GetName() const override;
    // Values stored by each integration phase.
    LUE std::string v_GetVariant() const override;
    LUE size_t v_GetOrder() const override;
    LUE std::vector<NekDouble> v_GetFreeParams() const override;

    LUE NekDouble v_GetTimeStability() const override;

    LUE TimeIntegrationSchemeType v_GetIntegrationSchemeType() const override;

    LUE size_t v_GetNumIntegrationPhases() const override;

    // Gets the solution Vector
    inline const TripleArray &v_GetSolutionVector() const override
    {
        return NullNekDoubleTensorOfArray3D;
    }

    // Sets the solution Vector
    inline void v_SetSolutionVector(
        [[maybe_unused]] const size_t Offset,
        [[maybe_unused]] const DoubleArray &y) override
    {
    }

    // The worker methods
    LUE void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE ConstDoubleArray &v_TimeIntegrate(const size_t timestep,
                                          const NekDouble delta_t) override;

    LUE void v_print(std::ostream &os) const override;
    LUE void v_printFull(std::ostream &os) const override;
    LUE TripleArray &v_UpdateSolutionVector() override;

    static std::string className;

protected:
    DoubleArray m_doubleArray;

}; // end class TimeIntegrationScheme

} // namespace Nektar::LibUtilities

#endif
