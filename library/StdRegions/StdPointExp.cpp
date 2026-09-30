///////////////////////////////////////////////////////////////////////////////
//
// File: StdPointExp.cpp
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

#include <StdRegions/StdPointExp.h>

using namespace std;

namespace Nektar::StdRegions
{

StdPointExp::StdPointExp(const LibUtilities::BasisKey &Ba)
    : StdExpansion(Ba.GetNumModes(), 1, Ba),
      StdExpansion0D(Ba.GetNumModes(), Ba)
{
}

LibUtilities::ShapeType StdPointExp::v_DetShapeType() const
{
    return LibUtilities::ePoint;
}

void StdPointExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                              [[maybe_unused]] Array<OneD, NekDouble> &coords_1,
                              [[maybe_unused]] Array<OneD, NekDouble> &coords_2)
{
    Vmath::Vcopy(GetNumPoints(0), (m_base[0]->GetZ()).data(), 1, &coords_0[0],
                 1);
}

void StdPointExp::PhysTensorDeriv(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &outarray)
{
    int nquad          = GetTotPoints();
    DNekMatSharedPtr D = m_base[0]->GetD();

    if (inarray.data() == outarray.data())
    {
        Array<OneD, NekDouble> wsp(nquad);
        CopyArray(inarray, wsp);
        Blas::Dgemv('N', nquad, nquad, 1.0, &(D->GetPtr())[0], nquad, &wsp[0],
                    1, 0.0, &outarray[0], 1);
    }
    else
    {
        Blas::Dgemv('N', nquad, nquad, 1.0, &(D->GetPtr())[0], nquad,
                    &inarray[0], 1, 0.0, &outarray[0], 1);
    }
}

void StdPointExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray)
{
    int nquad = m_base[0]->GetNumPoints();

    if (m_base[0]->Collocation())
    {
        Vmath::Vcopy(nquad, inarray, 1, outarray, 1);
    }
    else
    {
        Blas::Dgemv('N', nquad, m_base[0]->GetNumModes(), 1.0,
                    (m_base[0]->GetBdata()).data(), nquad, &inarray[0], 1, 0.0,
                    &outarray[0], 1);
    }
}

// Inner product
void StdPointExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    Array<OneD, const NekDouble> w = m_base[0]->GetW();
    outarray[0]                    = inarray[0] * w[0];
}

void StdPointExp::v_IProductWRTDerivBase(
    [[maybe_unused]] const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    ASSERTL1(dir >= 0 && dir < 1, "input dir is out of range");
    Array<OneD, const NekDouble> w = m_base[0]->GetW();
    outarray[0]                    = inarray[0] * w[0];
}

DNekMatSharedPtr StdPointExp::v_GenMatrix(const StdMatrixKey &mkey)
{
    DNekMatSharedPtr Mat;
    MatrixType mattype;

    switch (mattype = mkey.GetMatrixType())
    {
        case eFwdTrans:
        {
            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
            StdMatrixKey iprodkey(eIProductWRTBase, v_DetShapeType(), *this);
            DNekMat &Iprod = *GetStdMatrix(iprodkey);
            StdMatrixKey imasskey(eInvMass, v_DetShapeType(), *this);
            DNekMat &Imass = *GetStdMatrix(imasskey);

            (*Mat) = Imass * Iprod;
        }
        break;
        default:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
        }
        break;
    }

    return Mat;
}

DNekMatSharedPtr StdPointExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

} // namespace Nektar::StdRegions
