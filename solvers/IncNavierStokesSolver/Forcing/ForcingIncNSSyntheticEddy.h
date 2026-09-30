///////////////////////////////////////////////////////////////////////////////
//
// File:  ForcingIncNSSyntheticEddy.h
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
//              Incompressible solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGINCNSSYNTHETICEDDY
#define NEKTAR_SOLVERUTILS_FORCINGINCNSSYNTHETICEDDY

#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/Forcing/ForcingSyntheticEddy.h>
#include <string>

namespace Nektar::SolverUtils
{
class ForcingIncNSSyntheticEddy

    : virtual public SolverUtils::Forcing,
      virtual public SolverUtils::ForcingSyntheticEddy
{
public:
    friend class MemoryManager<ForcingIncNSSyntheticEddy>;

    /// Creates an instance of this class
    static ForcingSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
    {
        ForcingSharedPtr p =
            MemoryManager<ForcingIncNSSyntheticEddy>::AllocateSharedPtr(
                pSession, pEquation);
        p->InitObject(pFields, pNumForcingFields, pForce);
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    // Apply forcing term
    void v_Apply(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                 const Array<OneD, Array<OneD, NekDouble>> &inarray,
                 Array<OneD, Array<OneD, NekDouble>> &outarray,
                 const NekDouble &time) override;

    /// Calculate Forcing
    void CalculateForcing(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

private:
    ForcingIncNSSyntheticEddy(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation);
    ~ForcingIncNSSyntheticEddy(void) override {};
};

} // namespace Nektar::SolverUtils

#endif
