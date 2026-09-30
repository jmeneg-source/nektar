///////////////////////////////////////////////////////////////////////////////
//
// File: StdPointExp.h
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
// Description: Definition of a Point expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGIONS_STDPOINTEXP_H
#define NEKTAR_LIBS_STDREGIONS_STDPOINTEXP_H

#include <StdRegions/StdExpansion0D.h>

namespace Nektar::StdRegions
{
class StdPointExp : virtual public StdExpansion0D
{
public:
    STD_REGIONS_EXPORT StdPointExp(const LibUtilities::BasisKey &Ba);
    STD_REGIONS_EXPORT StdPointExp()                     = default;
    STD_REGIONS_EXPORT StdPointExp(const StdPointExp &T) = default;
    STD_REGIONS_EXPORT ~StdPointExp() override           = default;

protected:
    //----------------------------
    // Evaluations Methods
    //---------------------------
    STD_REGIONS_EXPORT void v_GetCoords(
        Array<OneD, NekDouble> &coords_0, Array<OneD, NekDouble> &coords_1,
        Array<OneD, NekDouble> &coords_2) override;

    //----------------------------
    // Helper functions
    //---------------------------
    STD_REGIONS_EXPORT LibUtilities::ShapeType v_DetShapeType() const override;

    STD_REGIONS_EXPORT void PhysTensorDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);

    //-----------------------------
    // Transforms
    //-----------------------------
    STD_REGIONS_EXPORT void v_BwdTrans(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //----------------------------
    // Inner product functions
    //----------------------------
    STD_REGIONS_EXPORT void v_IProductWRTBase(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    STD_REGIONS_EXPORT void v_IProductWRTDerivBase(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------
    // Evaluations Methods
    //---------------------------
    STD_REGIONS_EXPORT DNekMatSharedPtr
    v_GenMatrix(const StdMatrixKey &mkey) override;
    STD_REGIONS_EXPORT DNekMatSharedPtr
    v_CreateStdMatrix(const StdMatrixKey &mkey) override;

private:
    int v_GetNverts() const final
    {
        return 1;
    }

    int v_NumBndryCoeffs() const final
    {
        return 0;
    }

    int v_NumDGBndryCoeffs() const final
    {
        return 0;
    }

    int v_GetTraceNcoeffs([[maybe_unused]] const int i) const final
    {
        return 0;
    }

    int v_GetTraceIntNcoeffs([[maybe_unused]] const int i) const final
    {
        return 0;
    }

    int v_GetTraceNumPoints([[maybe_unused]] const int i) const final
    {
        return 0;
    }

    int v_GetVertexMap([[maybe_unused]] int localVertexId,
                       [[maybe_unused]] bool useCoeffPacking = false) override
    {
        ASSERTL2(localVertexId == 0, "Only single point in StdPointExp!");
        return 0;
    }
};

// type defines for use of PointExp in a boost vector
typedef std::shared_ptr<StdPointExp> StdPointExpSharedPtr;
} // namespace Nektar::StdRegions

#endif // STDPOINTEXP_H
