///////////////////////////////////////////////////////////////////////////////
//
// File: IncBoundaryConditions.h
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
// Description: Base class for boundary conditions.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_INCBOUNDARYCONDITIONS_H
#define NEKTAR_SOLVERS_INCBOUNDARYCONDITIONS_H
#include <IncNavierStokesSolver/BoundaryConditions/IncBaseCondition.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
class IncBoundaryConditions;
typedef std::shared_ptr<IncBoundaryConditions> IncBoundaryConditionsSharedPtr;

class IncBoundaryConditions
{
public:
    IncBoundaryConditions();

    void Initialize(const LibUtilities::SessionReaderSharedPtr pSession,
                    Array<OneD, MultiRegions::ExpListSharedPtr> pFields);

    void Update(const Array<OneD, const Array<OneD, NekDouble>> &fields,
                const Array<OneD, const Array<OneD, NekDouble>> &Adv,
                std::map<std::string, NekDouble> &params);

protected:
    std::map<int, IncBaseConditionSharedPtr> m_bounds;
    static std::set<std::string> m_BndType;
    int m_spacedim;
    int m_bnd_dim;
};

} // namespace Nektar

#endif
