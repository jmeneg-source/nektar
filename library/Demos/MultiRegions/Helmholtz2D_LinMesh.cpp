///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz2D_LinMesh.cpp
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
// Description: Solve a Helmholtz 2D problem by generating a linear
// mesh and interpolating the forcing to the linear mesh then
// interpolating the linear mesh back to the high order mesh
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField.h>
#include <SpatialDomains/LinearMeshGraph.hpp>
#include <SpatialDomains/MeshGraphIO.h>

using namespace Nektar;

// #define TIMING
#ifdef TIMING
#include <time.h>
#define Timing(s)                                                              \
    fprintf(stdout, "%s Took %g seconds\n", s,                                 \
            (clock() - st) / (double)CLOCKS_PER_SEC);                          \
    st = clock();
#else
#define Timing(s) /* Nothing */
#endif

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContFieldSharedPtr Exp, LinExp;
    int i, nq, coordim;
    Array<OneD, NekDouble> fce;
    Array<OneD, NekDouble> xc0, xc1, xc2;
    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;

    if (argc < 2)
    {
        std::cerr << "Usage: Helmholtz2D_LinMesh meshfile" << std::endl;
        return 1;
    }

    try
    {
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateDefault(vSession);

        //----------------------------------------------
        // Read in mesh from input file and generate Linear graph
        SpatialDomains::MeshGraphSharedPtr graph2D =
            SpatialDomains::MeshGraphIO::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
            vSession, graph2D, vSession->GetVariable(0));

        int nsplit = Exp->GetExp(0)->GetBasisNumModes(0) - 1;
        std::map<int, std::pair<int, std::vector<int>>> CoeffMap;

        SpatialDomains::LinearMeshGraph linGraph2D(graph2D);
        linGraph2D.CreateLinearGraph(nsplit, CoeffMap, false, false);
        auto linGraph2DShPtr = linGraph2D.GetLinearGraph();

        LinExp = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
            vSession, linGraph2DShPtr, vSession->GetVariable(0));
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionInfoMap &expansions =
            graph2D->GetExpansionInfo();
        LibUtilities::BasisKey bkey0 =
            expansions.begin()->second->m_basisKeyVector[0];

        if (vSession->GetComm()->GetRank() == 0)
        {
            std::cout << "Solving 2D Helmholtz: " << std::endl;
            std::cout << "         Communication: "
                      << vSession->GetComm()->GetType() << std::endl;
            std::cout << "         Solver type  : "
                      << vSession->GetSolverInfo("GlobalSysSoln") << std::endl;
            std::cout << "         Lambda       : "
                      << factors[StdRegions::eFactorLambda] << std::endl;
            std::cout << "         No. modes    : " << bkey0.GetNumModes()
                      << std::endl;
            std::cout << "         nsplit       : " << nsplit << std::endl;
            std::cout << std::endl;
        }
        //----------------------------------------------

        Timing("Read files and define exp ..");

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        coordim = Exp->GetCoordim(0);
        nq      = Exp->GetTotPoints();

        xc0 = Array<OneD, NekDouble>(nq, 0.0);
        xc1 = Array<OneD, NekDouble>(nq, 0.0);
        xc2 = Array<OneD, NekDouble>(nq, 0.0);

        switch (coordim)
        {
            case 2:
                Exp->GetCoords(xc0, xc1);
                break;
            case 3:
                Exp->GetCoords(xc0, xc1, xc2);
                break;
            default:
                ASSERTL0(false, "Coordim not valid");
                break;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function and interpolate to linear mesh
        fce = Array<OneD, NekDouble>(nq);
        LibUtilities::EquationSharedPtr ffunc =
            vSession->GetFunction("Forcing", 0);
        ffunc->Evaluate(xc0, xc1, xc2, fce);

        // Fwd & bwd Trans
        Array<OneD, NekDouble> fcecoeffs(nq);
        Exp->FwdTransLocalElmt(fce, fcecoeffs);
        Exp->BwdTrans(fcecoeffs, fce);

        // interpolate phys values to equispaced simplex points
        Array<OneD, NekDouble> equi(LinExp->GetNcoeffs()), tmp;
        int offset0 = 0;
        int offset1 = 0;
        for (auto &elmt : *(Exp->GetExp()))
        {
            ASSERTL0(offset1 < equi.size(),
                     "offset1 larger than vector size of equi");
            elmt->PhysInterpToSimplexEquiSpaced(
                fce + offset0, tmp = equi + offset1, nsplit + 1);
            offset0 += elmt->GetTotPoints();
            offset1 += (elmt->DetShapeType() == LibUtilities::eQuadrilateral)
                           ? (nsplit + 1) * (nsplit + 1)
                           : (nsplit + 1) * (nsplit + 2) / 2;
        }

        // extract high order to low order mapping
        std::map<int, int> offset;
        for (int e = 0; e < Exp->GetNumElmts(); ++e)
        {
            offset[Exp->GetExp(e)->GetGeom()->GetGlobalID()] =
                Exp->GetCoeff_Offset(e);
        }

        Array<OneD, int> ho2lor = Array<OneD, int>(LinExp->GetNcoeffs());
        int cnt                 = 0;
        for (int e = 0; e < LinExp->GetNumElmts(); ++e)
        {
            int gid = LinExp->GetExp(e)->GetGeom()->GetGlobalID();
            int id  = CoeffMap[gid].first;

            for (int i = 0; i < CoeffMap[gid].second.size(); ++i)
            {
                ho2lor[cnt++] = offset[id] + CoeffMap[gid].second[i];
            }
        }

        // fill cooef space
        Vmath::Gathr(LinExp->GetNcoeffs(), equi, ho2lor,
                     LinExp->UpdateCoeffs());
        // evaluate at phys points of linear expansion
        LinExp->BwdTrans(LinExp->GetCoeffs(), LinExp->UpdatePhys());
        //----------------------------------------------
        Timing("Define forcing ..");

        //----------------------------------------------
        // Helmholtz solution taking physical forcing after setting
        // initial condition to zero
        Vmath::Zero(LinExp->GetNcoeffs(), LinExp->UpdateCoeffs(), 1);
        LinExp->HelmSolve(LinExp->GetPhys(), LinExp->UpdateCoeffs(), factors,
                          varcoeffs);
        //----------------------------------------------
        Timing("Helmholtz Solve ..");

        //----------------------------------------------
        // Project solution back to p-expansion and put in phys space
        Vmath::Scatr(LinExp->GetNcoeffs(), LinExp->GetCoeffs(), ho2lor, equi);
        offset0 = offset1 = 0;
        for (auto &elmt : *(Exp->GetExp()))
        {
            elmt->EquiSpacedToPhys(nsplit + 1, equi + offset0,
                                   tmp = Exp->UpdatePhys() + offset1);

            offset0 += (elmt->DetShapeType() == LibUtilities::eQuadrilateral)
                           ? (nsplit + 1) * (nsplit + 1)
                           : (nsplit + 1) * (nsplit + 2) / 2;
            offset1 += elmt->GetTotPoints();
        }
        Exp->FwdTransLocalElmt(Exp->GetPhys(), Exp->UpdateCoeffs());
        //----------------------------------------------
        Timing("Project points back ..");

        //-----------------------------------------------
        // Write solution and graph to file
        std::string out = vSession->GetSessionName() + "_LinMesh.fld";
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
            LinExp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());

        for (i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            LinExp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        fld->Write(out, FieldDef, FieldData);

        std::string outname = vSession->GetSessionName() + "_linMesh.xml";
        auto graphIO =
            SpatialDomains::GetMeshGraphIOFactory().CreateInstance("Xml");
        graphIO->SetMeshGraph(linGraph2DShPtr);
        graphIO->WriteGeometry(outname);
        //-----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol =
            vSession->GetFunction("ExactSolution", 0);

        if (ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution
            ex_sol->Evaluate(xc0, xc1, xc2, fce);
            //--------------------------------------------

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Exp->GetPhys(), fce);
            NekDouble vL2Error   = Exp->L2(Exp->GetPhys(), fce);
            NekDouble vH1Error   = Exp->H1(Exp->GetPhys(), fce);
            if (vSession->GetComm()->GetRank() == 0)
            {
                std::cout << "L infinity error: " << vLinfError << std::endl;
                std::cout << "L 2 error:        " << vL2Error << std::endl;
                std::cout << "H 1 error:        " << vH1Error << std::endl;
            }
            //--------------------------------------------
        }
        //----------------------------------------------
    }
    catch (const std::runtime_error &)
    {
        std::cout << "Caught an error" << std::endl;
        return 1;
    }

    vSession->Finalise();

    return 0;
}
