///////////////////////////////////////////////////////////////////////////////
//
// File: UndefinedInOutflow.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: UndefinedInOutflow class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/UndefinedInOutflow.h>

namespace Nektar
{

std::string UndefinedInOutflow::className =
    GetBoundaryFactory().RegisterCreatorFunction(
        "NoUserDefined", UndefinedInOutflow::create, "No boundary condition");

UndefinedInOutflow::UndefinedInOutflow(
    Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
    const LibUtilities::SessionReaderSharedPtr pSession,
    PulseWavePressureAreaSharedPtr pressureArea)
    : PulseWaveBoundary(pVessel, pSession, pressureArea)
{
}

void UndefinedInOutflow::v_DoBoundary(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &A_0,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &beta,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &alpha,
    [[maybe_unused]] const NekDouble time, [[maybe_unused]] int omega,
    [[maybe_unused]] int offset, [[maybe_unused]] int n)
{
}

} // namespace Nektar
