///////////////////////////////////////////////////////////////////////////////
//
// File: PreconditionerLOR.cpp
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
// Description: PreconditionerLOR definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysDirectFull.h>
#include <MultiRegions/GlobalLinSysIterativeFull.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/PreconditionerLOR.h>
#include <SpatialDomains/LinearMeshGraph.hpp>
#include <SpatialDomains/MeshGraphIO.h>

#ifdef NEKTAR_USE_MPI
#include <MultiRegions/GlobalLinSysXxtFull.h>
#endif

#ifdef NEKTAR_USING_PETSC
#include <MultiRegions/GlobalLinSysPETScFull.h>
#endif

#include <cmath>

namespace Nektar::MultiRegions
{
/**
 * Registers the class with the Factory.
 */
std::string PreconditionerLOR::className =
    GetPreconFactory().RegisterCreatorFunction("LOR", PreconditionerLOR::create,
                                               "LOR Preconditioning");

/**
 * @class PreconditionerLOR
 *
 * This class implements LOR preconditioning for the conjugate gradient matrix
 * solver.
 */
PreconditionerLOR::PreconditionerLOR(
    const std::shared_ptr<GlobalLinSys> &plinsys,
    const AssemblyMapSharedPtr &pLocToGloMap)
    : Preconditioner(plinsys, pLocToGloMap), m_verboseIter(false)
{
}

void PreconditionerLOR::v_InitObject()
{
}

void PreconditionerLOR::v_BuildPreconditioner()
{
    GlobalSysSolnType solvertype = m_locToGloMap.lock()->GetGlobalSysSolnType();

    ASSERTL0(solvertype == eIterativeFull, "Only works for full iterative");

    // Grab explist
    auto expList = ((m_linsys.lock())->GetLocMat()).lock();

    std::shared_ptr<LibUtilities::SessionReader> session =
        expList->GetSession();

    m_nsplit = expList->GetExp(0)->GetBasisNumModes(0) - 1;

    std::string var = m_locToGloMap.lock()->GetVariable();

    m_equiSpaced = false;
    if (session->DefinesGlobalSysSolnInfo(var, "LORPointDistribution"))
    {
        if (boost::iequals("GLL", session->GetGlobalSysSolnInfo(
                                      var, "LORPointDistribution")))
        {
            m_equiSpaced = false;
        }
        else if (boost::iequals("Equispaced", session->GetGlobalSysSolnInfo(
                                                  var, "LORPointDistribution")))
        {
            m_equiSpaced = true;
        }
    }

    m_useSimplex = false;
    if (session->DefinesGlobalSysSolnInfo(var, "LORSimplex"))
    {
        if (boost::iequals("True",
                           session->GetGlobalSysSolnInfo(var, "LORSimplex")))
        {
            m_useSimplex = true;
        }
    }

    // Setup session !m_equispaced since true means GLL distribution
    std::map<int, std::pair<int, std::vector<int>>> coeffmap;
    SpatialDomains::LinearMeshGraph linear(expList->GetGraph());
    linear.CreateLinearGraph(m_nsplit, coeffmap, !m_equiSpaced, m_useSimplex);
    m_lor_graph = linear.GetLinearGraph();

    std::string preconType;

    if (session->DefinesGlobalSysSolnInfo(var, "LORSolverType"))
    {
        m_slvType = session->GetGlobalSysSolnInfo(var, "LORSolverType");
    }
    else
    {
        NEKERROR(
            ErrorUtil::efatal,
            "Need to define LORSolverType "
            "(DirectFull, IterativeFull, Preconditioner, PETScFull, XxtFull)"
            "with LOR Preconditioner");
    }

    if (session->GetComm()->GetRank() == 0)
    {
        std::cout << "LOR Solve Type: " << m_slvType << std::endl;
    }

    std::string LinSysIterSolver;
    if (session->DefinesGlobalSysSolnInfo(var, "LORLinSysIterSolver"))
    {
        LinSysIterSolver =
            session->GetGlobalSysSolnInfo(var, "LORLinSysIterSolver");
    }
    else
    {
        LinSysIterSolver = "ConjugateGradient";
    }

    if (session->DefinesGlobalSysSolnInfo(var, "LORPreconditioner"))
    {
        preconType = session->GetGlobalSysSolnInfo(var, "LORPreconditioner");
    }
    else
    {
        preconType = "Diagonal";
    }

    if (boost::iequals(m_slvType, "Preconditioner"))
    {
        session->SetGlobalSysSolnInfo(var, "GlobalSysSoln", "IterativeFull");
    }
    else
    {
        session->SetGlobalSysSolnInfo(var, "GlobalSysSoln", m_slvType.c_str());

        if (boost::iequals(m_slvType, "IterativeFull"))
        {
            m_verboseIter =
                expList->GetSession()->DefinesCmdLineArgument("verbose");
        }
    }

    session->SetGlobalSysSolnInfo(var, "Preconditioner", preconType.c_str());

    // setup Fixed iteration if difined in GlobalSysSoln section
    std::string prevFixedIter;
    if (session->DefinesGlobalSysSolnInfo(var, "LORFixedIterations"))
    {

        if (session->DefinesGlobalSysSolnInfo(var, "NekLinSysFixedIterations"))
        {
            prevFixedIter =
                session->GetGlobalSysSolnInfo(var, "NekLinSysFixedIterations");
        }

        session->SetGlobalSysSolnInfo(
            var, "NekLinSysFixedIterations",
            session->GetGlobalSysSolnInfo(var, "LORFixedIterations"));
    }

    std::string prevIterTol;
    if (session->DefinesGlobalSysSolnInfo(var, "LORIterativeTolerance"))
    {

        if (session->DefinesGlobalSysSolnInfo(var, "IterativeSolverTolerance"))
        {
            prevIterTol =
                session->GetGlobalSysSolnInfo(var, "IterativeSolverTolerance");
        }

        session->SetGlobalSysSolnInfo(
            var, "IterativeSolverTolerance",
            session->GetGlobalSysSolnInfo(var, "LORIterativeTolerance"));
    }

    session->SetGlobalSysSolnInfo(var, "LinSysIterSolver",
                                  LinSysIterSolver.c_str());

    // Create contfield
    m_lor_field = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        session, m_lor_graph, var);

    // Set up high-order to LOR map
    m_ho2lor = Array<OneD, int>(m_lor_field->GetNcoeffs());
    std::map<int, int> offset;

    // check if the coeffmap has been setup correctly in the LOR space
    for (int e = 0; e < expList->GetNumElmts(); ++e)
    {
        offset[expList->GetExp(e)->GetGeom()->GetGlobalID()] =
            expList->GetCoeff_Offset(e);
    }

    int cnt = 0;
    for (int e = 0; e < m_lor_field->GetNumElmts(); ++e)
    {
        int gid = m_lor_field->GetExp(e)->GetGeom()->GetGlobalID();
        int id  = coeffmap[gid].first;

        for (int i = 0; i < coeffmap[gid].second.size(); ++i)
        {
            m_ho2lor[cnt++] = offset[id] + coeffmap[gid].second[i];
        }
    }
    ASSERTL1(cnt == m_lor_field->GetNcoeffs(),
             "Error in setting up high order to lower order map");

    // for preconditioner have a homogeneous system so need to set boundary
    // conditions to homogeneous too
    m_lor_field->SetBCsToHomogeneous();

    GlobalLinSysKey preconKey(m_linsys.lock()->GetKey().GetMatrixType(),
                              m_lor_field->GetLocalToGlobalMap(),
                              m_linsys.lock()->GetKey().GetConstFactors());

    if (boost::iequals(m_slvType, "DirectFull"))
    {
        m_LORLinSys = MemoryManager<GlobalLinSysDirectFull>::AllocateSharedPtr(
            preconKey, m_lor_field, m_lor_field->GetLocalToGlobalMap());
    }
    else if (boost::iequals(m_slvType, "XxtFull"))
    {
#ifdef NEKTAR_USE_MPI
        m_LORLinSys = MemoryManager<GlobalLinSysXxtFull>::AllocateSharedPtr(
            preconKey, m_lor_field, m_lor_field->GetLocalToGlobalMap());
#else
        ASSERTL0(false, "Nektar++ has not been compiled with MPI support.");
#endif
    }
    else if (boost::iequals(m_slvType, "IterativeFull"))
    {
        m_LORLinSys =
            MemoryManager<GlobalLinSysIterativeFull>::AllocateSharedPtr(
                preconKey, m_lor_field, m_lor_field->GetLocalToGlobalMap());
    }
    else if (boost::iequals(m_slvType, "PETScFull"))
    {
#ifdef NEKTAR_USING_PETSC
        m_LORLinSys = MemoryManager<GlobalLinSysPETScFull>::AllocateSharedPtr(
            preconKey, m_lor_field, m_lor_field->GetLocalToGlobalMap());
#else
        ASSERTL0(false, "Nektar++ has not been compiled with PETSc support.");
#endif
    }
    else if (boost::iequals(m_slvType, "Preconditioner"))
    {
        AssemblyMapSharedPtr l2gmap = m_lor_field->GetLocalToGlobalMap();

        m_LORLinSys =
            MemoryManager<GlobalLinSysIterativeFull>::AllocateSharedPtr(
                preconKey, m_lor_field, l2gmap);

        std::dynamic_pointer_cast<GlobalLinSysIterativeFull>(m_LORLinSys)
            ->Initialise(l2gmap->GetNumGlobalCoeffs(), l2gmap,
                         l2gmap->GetNumGlobalDirBndCoeffs());
    }
    else
    {
        NEKERROR(ErrorUtil::efatal,
                 "LORSlvType must be specified as one of "
                 "DirectFull, IterativeFull, Preconditioner, PETScFull");
    }

    // reset GlobalSysSolnType
    session->SetGlobalSysSolnInfo(var, "GlobalSysSoln",
                                  GlobalSysSolnTypeMap[solvertype]);

    // set fixed iterations if required;
    if (prevFixedIter.size())
    {
        session->SetGlobalSysSolnInfo(var, "NekLinSysFixedIterations",
                                      prevFixedIter);
    }

    // reset iterative tolerance
    if (prevIterTol.size())
    {
        session->SetGlobalSysSolnInfo(var, "IterativeSolverTolerance",
                                      prevIterTol);
    }

    CreateInvMultiplicity();

    // store ndir for easy access in DoPreconditioner
    m_nDir = m_locToGloMap.lock()->GetNumGlobalDirBndCoeffs();
}

/**
 *
 */
void PreconditionerLOR::v_DoPreconditioner(const Array<OneD, NekDouble> &pInput,
                                           Array<OneD, NekDouble> &pOutput,
                                           const bool &isLocal)
{
    auto expList = ((m_linsys.lock())->GetLocMat()).lock();

    int nLocal     = expList->GetNcoeffs();
    int lor_nLocal = m_lor_field->GetNcoeffs();

    Array<OneD, NekDouble> tmp, solve(std::max(nLocal, lor_nLocal), 0.0);
    Array<OneD, NekDouble> Local(nLocal);

    if (isLocal)
    {
        // Input array is in local (non-assembled format) and to reproduce the
        // globally assembled problem the LOR preconditioner works best if we
        // take the averaged value so assemble and then divide by the
        // multiplicity.
        //
        // Note: Not entirely clear to me why we cannot take the input value
        m_locToGloMap.lock()->Assemble(pInput, solve);
        Vmath::Vmul(m_invMultiplicity.size() - m_nDir,
                    m_invMultiplicity + m_nDir, 1, solve + m_nDir, 1,
                    tmp = solve + m_nDir, 1);
        Vmath::Zero(m_nDir, solve, 1);
        expList->GlobalToLocal(solve, Local);
    }
    else // recover local forcing coefficients
    {
        Vmath::Vmul(m_invMultiplicity.size() - m_nDir,
                    m_invMultiplicity + m_nDir, 1, pInput, 1,
                    tmp = solve + m_nDir, 1);
        Vmath::Zero(m_nDir, solve, 1);

        expList->GlobalToLocal(solve, Local);
    }

    StdRegions::MatrixType mtype = (m_equiSpaced)
                                       ? StdRegions::eEquiSpacedToCoeffs
                                       : StdRegions::eGLLToCoeffs;
    GlobalMatrixKey key(mtype, m_locToGloMap.lock(),
                        m_linsys.lock()->GetKey().GetConstFactors());

    // Transform inner product to equispaced polynomial  with transpose multiply
    expList->MultiplyByBlockMatrix(key, Local, solve, true);

    // Get force for LOR system
    Array<OneD, NekDouble> fce(lor_nLocal);
    // makes inner system symmetric by dividing by linear mesh elemental
    // multiplicity
    Vmath::Vmul(nLocal, m_invLinMeshMultiplicity, 1, solve, 1, solve, 1);
    Vmath::Gathr(lor_nLocal, solve, m_ho2lor, fce);

    // solve linear sys solve for low order system
    AssemblyMapSharedPtr l2gmap = m_lor_field->GetLocalToGlobalMap();
    Vmath::Zero(lor_nLocal, solve, 1);

    if (boost::iequals(m_slvType, "DirectFull") ||
        boost::iequals(m_slvType, "IterativeFull") ||
        boost::iequals(m_slvType, "XxtFull") ||
        boost::iequals(m_slvType, "PETScFull"))
    {
        m_LORLinSys->SolveLinearSystem(l2gmap->GetNumGlobalCoeffs(), fce, solve,
                                       l2gmap,
                                       l2gmap->GetNumGlobalDirBndCoeffs());
    }
    else
    {
        std::dynamic_pointer_cast<GlobalLinSysIterativeFull>(m_LORLinSys)
            ->DoPreconditionerFlag(fce, solve, true);
    }

    // Project solution back to p-expansion and in Lagrange space
    Vmath::Scatr(lor_nLocal, solve, m_ho2lor, Local);

    if (isLocal)
    {
        // Transform residual to  polynomial

        expList->MultiplyByBlockMatrix(key, Local, pOutput);

        // enforce zero Bcs - should be done with a mask on local variables.
        expList->LocalToGlobal(pOutput, solve);
        Vmath::Zero(m_nDir, solve, 1);
        expList->GlobalToLocal(solve, pOutput);
    }
    else
    {
        // Transform residual to polynomial
        expList->MultiplyByBlockMatrix(key, Local, solve);
        expList->LocalToGlobal(solve, Local);
        Vmath::Vcopy(pOutput.size(), Local + m_nDir, 1, pOutput, 1);
    }
}

/**
 * Create inv multiplicity
 */
void PreconditionerLOR::CreateInvMultiplicity(void)
{
    auto asmMap      = m_locToGloMap.lock();
    unsigned nGlobal = asmMap->GetNumGlobalCoeffs();
    unsigned nLocal  = asmMap->GetNumLocalCoeffs();

    // Construct a mask array
    m_invMultiplicity = Array<OneD, NekDouble>(nGlobal, 0.0);

    const Array<OneD, const int> &map = asmMap->GetLocalToGlobalMap();

    for (unsigned i = 0; i < nLocal; ++i)
    {
        m_invMultiplicity[map[i]] += 1.0;
    }

    asmMap->UniversalAssemble(m_invMultiplicity);
    Vmath::Sdiv(nGlobal, 1.0, m_invMultiplicity, 1, m_invMultiplicity, 1);

    auto lor_nLocal          = m_ho2lor.size();
    m_invLinMeshMultiplicity = Array<OneD, NekDouble>(nLocal, 0.0);

    for (unsigned i = 0; i < lor_nLocal; ++i)
    {
        m_invLinMeshMultiplicity[m_ho2lor[i]] += 1.0;
    }
    for (unsigned i = 0; i < nLocal; ++i)
    {
        m_invLinMeshMultiplicity[i] = 1.0 / m_invLinMeshMultiplicity[i];
    }
}

} // namespace Nektar::MultiRegions
