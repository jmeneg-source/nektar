///////////////////////////////////////////////////////////////////////////////
//
// File: UnsteadyReactionDiffusion.cpp
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
// Description: Unsteady reaction-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/UnsteadyReactionDiffusion.h>
#include <iomanip>
#include <iostream>

using namespace std;

namespace Nektar
{
string UnsteadyReactionDiffusion::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "UnsteadyReactionDiffusion", UnsteadyReactionDiffusion::create);

UnsteadyReactionDiffusion::UnsteadyReactionDiffusion(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadyDiffusion(pSession, pGraph)
{
}

/**
 * @brief Initialisation object for the unsteady reaction-diffusion problem.
 */
void UnsteadyReactionDiffusion::v_InitObject(bool DeclareFields)
{
    UnsteadyDiffusion::v_InitObject(DeclareFields);

    // Forcing terms
    m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                           m_fields, m_fields.size());

    // Reset OdeRhs functor
    m_ode.DefineOdeRhs(&UnsteadyReactionDiffusion::DoOdeRhs, this);
}

/**
 * @brief Compute the right-hand side for the unsteady reaction diffusion
 * problem.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyReactionDiffusion::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    if (m_explicitDiffusion)
    {
        UnsteadyDiffusion::DoOdeRhs(inarray, outarray, time);
    }
    else
    {
        // RHS should be set to zero.
        for (int i = 0; i < outarray.size(); ++i)
        {
            Vmath::Zero(outarray[i].size(), &outarray[i][0], 1);
        }
    }

    // Add forcing terms for reaction.
    for (auto &x : m_forcing)
    {
        // set up non-linear terms
        x->Apply(m_fields, inarray, outarray, time);
    }
}

} // namespace Nektar
