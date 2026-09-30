///////////////////////////////////////////////////////////////////////////////
//
// File: Diffusion.cpp
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
// Description: Abstract base class for diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar::SolverUtils
{
DiffusionFactory &GetDiffusionFactory()
{
    static DiffusionFactory instance;
    return instance;
}

void Diffusion::v_DiffuseCoeffs(
    [[maybe_unused]] const std::size_t nConvectiveFields,
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    NEKERROR(ErrorUtil::efatal, "v_DiffuseCoeffs not defined");
}

const Array<OneD, const Array<OneD, NekDouble>> &Diffusion::v_GetTraceNormal()
{
    NEKERROR(ErrorUtil::efatal, "v_GetTraceNormal not defined");
    return NullNekDoubleArrayOfArray;
}

void Diffusion::v_DiffuseCalcDerivative(
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] TensorOfArray3D<NekDouble> &qfields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    NEKERROR(ErrorUtil::efatal, "Not defined for function DiffuseVolumeFLux.");
}

void Diffusion::v_DiffuseVolumeFlux(
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] TensorOfArray3D<NekDouble> &qfields,
    [[maybe_unused]] TensorOfArray3D<NekDouble> &VolumeFlux,
    [[maybe_unused]] Array<OneD, int> &nonZeroIndex)
{
    NEKERROR(ErrorUtil::efatal, "Not defined for function DiffuseVolumeFLux.");
}

void Diffusion::v_DiffuseTraceFlux(
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] TensorOfArray3D<NekDouble> &qfields,
    [[maybe_unused]] TensorOfArray3D<NekDouble> &VolumeFlux,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &TraceFlux,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pBwd,
    [[maybe_unused]] Array<OneD, int> &nonZeroIndex)
{
    NEKERROR(ErrorUtil::efatal, "Not defined function DiffuseTraceFLux.");
}

} // namespace Nektar::SolverUtils
