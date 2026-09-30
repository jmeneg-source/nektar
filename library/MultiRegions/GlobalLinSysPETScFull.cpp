///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysPETScFull.cpp
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
// Description: GlobalLinSysPETScFull definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList.h>
#include <MultiRegions/GlobalLinSysPETScFull.h>

#include "petscao.h"
#include "petscis.h"

using namespace std;

namespace Nektar::MultiRegions
{
/**
 * @class GlobalLinSysPETScFull
 */

/**
 * Registers the class with the Factory.
 */
string GlobalLinSysPETScFull::className =
    GetGlobalLinSysFactory().RegisterCreatorFunction(
        "PETScFull", GlobalLinSysPETScFull::create, "PETSc Full Matrix.");

/// Constructor for full direct matrix solve.
GlobalLinSysPETScFull::GlobalLinSysPETScFull(
    const GlobalLinSysKey &pLinSysKey, const std::weak_ptr<ExpList> &pExp,
    const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    : GlobalLinSys(pLinSysKey, pExp, pLocToGloMap),
      GlobalLinSysPETSc(pLinSysKey, pExp, pLocToGloMap)
{
    const int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

    int i, j, n, cnt, gid1, gid2, loc_lda;
    NekDouble sign1, sign2, value;
    DNekScalMatSharedPtr loc_mat;

    // CALCULATE REORDERING MAPPING
    CalculateReordering(pLocToGloMap->GetGlobalToUniversalMap(),
                        pLocToGloMap->GetGlobalToUniversalMapUnique(),
                        pLocToGloMap);

    // SET UP VECTORS AND MATRIX
    SetUpMatVec(pLocToGloMap->GetNumGlobalCoeffs(), nDirDofs);

    // SET UP SCATTER OBJECTS
    SetUpScatter();

    // CONSTRUCT KSP OBJECT
    LibUtilities::SessionReaderSharedPtr pSession =
        m_expList.lock()->GetSession();
    string variable                    = pLocToGloMap->GetVariable();
    NekDouble IterativeSolverTolerance = NekConstants::kNekIterativeTol;
    if (pSession->DefinesGlobalSysSolnInfo(variable,
                                           "IterativeSolverTolerance"))
    {
        IterativeSolverTolerance = boost::lexical_cast<double>(
            pSession->GetGlobalSysSolnInfo(variable, "IterativeSolverTolerance")
                .c_str());
    }
    else if (pSession->DefinesParameter("IterativeSolverTolerance"))
    {
        IterativeSolverTolerance =
            pSession->GetParameter("IterativeSolverTolerance");
    }
    SetUpSolver(IterativeSolverTolerance);

    // POPULATE MATRIX
    for (n = cnt = 0; n < m_expList.lock()->GetNumElmts(); ++n)
    {
        loc_mat = GetBlock(n);
        loc_lda = loc_mat->GetRows();

        for (i = 0; i < loc_lda; ++i)
        {
            gid1  = pLocToGloMap->GetLocalToGlobalMap(cnt + i) - nDirDofs;
            sign1 = pLocToGloMap->GetLocalToGlobalSign(cnt + i);
            if (gid1 >= 0)
            {
                int gid1ro = m_reorderedMap[gid1];
                for (j = 0; j < loc_lda; ++j)
                {
                    gid2 =
                        pLocToGloMap->GetLocalToGlobalMap(cnt + j) - nDirDofs;
                    sign2 = pLocToGloMap->GetLocalToGlobalSign(cnt + j);
                    if (gid2 >= 0)
                    {
                        int gid2ro = m_reorderedMap[gid2];
                        value      = sign1 * sign2 * (*loc_mat)(i, j);
                        MatSetValue(m_matrix, gid1ro, gid2ro, value,
                                    ADD_VALUES);
                    }
                }
            }
        }
        cnt += loc_lda;
    }

    // ASSEMBLE MATRIX
    MatAssemblyBegin(m_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_matrix, MAT_FINAL_ASSEMBLY);
}

/**
 * Solve the linear system using a full global matrix system.
 */
void GlobalLinSysPETScFull::v_Solve(
    const Array<OneD, const NekDouble> &pLocInput,
    Array<OneD, NekDouble> &pLocOutput,
    const AssemblyMapSharedPtr &pLocToGloMap,
    const Array<OneD, const NekDouble> &pDirForcing)
{
    m_locToGloMap = pLocToGloMap;

    bool dirForcCalculated = (bool)pDirForcing.size();
    int nDirDofs           = pLocToGloMap->GetNumGlobalDirBndCoeffs();
    int nGlobDofs          = pLocToGloMap->GetNumGlobalCoeffs();
    int nLocDofs           = pLocToGloMap->GetNumLocalCoeffs();

    int nDirTotal                                  = nDirDofs;
    std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
    expList->GetComm()->GetRowComm()->AllReduce(nDirTotal,
                                                LibUtilities::ReduceSum);

    if (nDirTotal)
    {
        Array<OneD, NekDouble> rhs(nLocDofs);

        // Calculate the Dirichlet forcing
        if (dirForcCalculated)
        {
            // Assume pDirForcing is in local space
            ASSERTL0(
                pDirForcing.size() >= nLocDofs,
                "DirForcing is not of sufficient size. Is it in local space?");
            Vmath::Vsub(nLocDofs, pLocInput, 1, pDirForcing, 1, rhs, 1);
        }
        else
        {
            // Calculate initial condition and Dirichlet forcing and subtract it
            // from the rhs
            expList->GeneralMatrixOp(m_linSysKey, pLocOutput, rhs);

            // Iterate over all the elements computing Robin BCs where
            // necessary
            for (auto &r : m_robinBCInfo) // add robin mass matrix
            {
                RobinBCInfoSharedPtr rBC;
                Array<OneD, NekDouble> rhsloc;

                int n      = r.first;
                int offset = expList->GetCoeff_Offset(n);

                LocalRegions::ExpansionSharedPtr vExp = expList->GetExp(n);
                // Add local matrix contribution
                for (rBC = r.second; rBC; rBC = rBC->next)
                {
                    vExp->AddRobinTraceContribution(
                        rBC->m_robinID, rBC->m_robinPrimitiveCoeffs,
                        pLocOutput + offset, rhsloc = rhs + offset);
                }
            }
            Vmath::Vsub(nLocDofs, pLocInput, 1, rhs, 1, rhs, 1);
        }

        Array<OneD, NekDouble> diff(nLocDofs);

        // Solve for perturbation from initial guess in pOutput
        SolveLinearSystem(nGlobDofs, rhs, diff, pLocToGloMap, nDirDofs);

        // Add back initial and boundary condition
        Vmath::Vadd(nLocDofs, diff, 1, pLocOutput, 1, pLocOutput, 1);
    }
    else
    {
        SolveLinearSystem(nGlobDofs, pLocInput, pLocOutput, pLocToGloMap,
                          nDirDofs);
    }
}

/**
 * @brief Apply matrix-vector multiplication using local approach and
 * the assembly map.
 *
 * @param input   Vector input.
 * @param output  Result of multiplication.
 */
void GlobalLinSysPETScFull::v_DoMatrixMultiply(
    const Array<OneD, const NekDouble> &input, Array<OneD, NekDouble> &output)
{
    std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();

    int nLocDofs = m_locToGloMap->GetNumLocalCoeffs();

    Array<OneD, NekDouble> tmp(nLocDofs);
    Array<OneD, NekDouble> tmp1(nLocDofs);

    m_locToGloMap->GlobalToLocal(input, tmp);

    // Perform matrix-vector operation A*d_i
    expList->GeneralMatrixOp(m_linSysKey, tmp, tmp1);

    m_locToGloMap->Assemble(tmp1, output);
}

} // namespace Nektar::MultiRegions
