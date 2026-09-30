///////////////////////////////////////////////////////////////////////////////
//
// File: NonlinearPeregrine.h
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
// Description: Nonlinear Boussinesq equations of Peregrine in conservative
//              variables
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_NONLINEARPEREGRINE_H
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_NONLINEARPEREGRINE_H

#include <ShallowWaterSolver/EquationSystems/NonlinearSWE.h>

namespace Nektar
{

enum ProblemType
{
    eGeneral,        ///< No problem defined - Default Inital data
    eSolitaryWave,   ///< First order Laitone solitary wave
    SIZE_ProblemType ///< Length of enum list
};

const char *const ProblemTypeMap[] = {"General", "SolitaryWave"};

class NonlinearPeregrine : public NonlinearSWE
{
public:
    friend class MemoryManager<NonlinearPeregrine>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<NonlinearPeregrine>::AllocateSharedPtr(pSession,
                                                                 pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

protected:
    NonlinearPeregrine(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const SpatialDomains::MeshGraphSharedPtr &pGraph);

    ~NonlinearPeregrine() override = default;

    void v_InitObject(bool DeclareFields = true) override;

    void v_SetInitialConditions(NekDouble initialtime      = 0.0,
                                bool dumpInitialConditions = true,
                                const int domain           = 0) override;

    void v_DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                    Array<OneD, Array<OneD, NekDouble>> &outarray,
                    const NekDouble time) override;

    void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    void LaitoneSolitaryWave(NekDouble amp, NekDouble d, NekDouble time,
                             NekDouble x_offset);

    // Dispersive parts
    void WCESolve(Array<OneD, NekDouble> &fce, NekDouble lambda);

    void SetBoundaryConditionsForcing(
        Array<OneD, Array<OneD, NekDouble>> &inarray, NekDouble time);

    void WallBoundaryForcing(int bcRegion, int cnt,
                             Array<OneD, Array<OneD, NekDouble>> &inarray);

    void NumericalFluxForcing(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, NekDouble> &numfluxX, Array<OneD, NekDouble> &numfluxY);

    void SetBoundaryConditionsContVariables(
        const Array<OneD, const NekDouble> &inarray, NekDouble time);

    void WallBoundaryContVariables(int bcRegion, int cnt,
                                   const Array<OneD, const NekDouble> &inarray);

    void NumericalFluxConsVariables(
        const Array<OneD, const NekDouble> &physfield,
        Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY);

private:
    StdRegions::ConstFactorMap m_factors;
    ProblemType m_problemType;
    NekDouble m_const_depth;
};

} // namespace Nektar

#endif
