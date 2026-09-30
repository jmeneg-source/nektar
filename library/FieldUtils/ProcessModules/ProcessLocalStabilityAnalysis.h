////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLocalStabilityAnalysis.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Local linear stability analysis of compressible flow.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSLOCALSTABILITYANALYSIS
#define FIELDUTILS_PROCESSLOCALSTABILITYANALYSIS

#include "ProcessBoundaryExtract.h"

#include <string>

namespace Nektar::FieldUtils
{

/**
 * @brief This processing module calculates the wall shear stress and adds it
 * as an extra-field to the output file, and writes it to a surface output file.
 */

// Define the

extern "C"
{
    void F77NAME(helloworld)();
    void F77NAME(copse3d)(const long int &option, const NekDouble *finalValue1,
                          const long int &numStep1,
                          const NekDouble *finalValue2,
                          const long int &numStep2);
}

class ProcessLocalStabilityAnalysis : public ProcessBoundaryExtract
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessLocalStabilityAnalysis>::AllocateSharedPtr(
            f);
    }
    static ModuleKey className;

    ProcessLocalStabilityAnalysis(FieldSharedPtr f);
    ~ProcessLocalStabilityAnalysis() override;

    /// Write mesh to output file.
    void v_Process(po::variables_map &vm) override;

    std::string v_GetModuleName() override
    {
        return "ProcessLocalStabilityAnalysis";
    }

    std::string v_GetModuleDescription() override
    {
        return "Performing local (linear) stability analysis (LST)";
    }

    // test
    inline void call_hello()
    {
        F77NAME(helloworld)();
    }

    inline void call_lst(const long int &option, const NekDouble *finalValue1,
                         const long int &numStep1, const NekDouble *finalValue2,
                         const long int &numStep2)
    {
        F77NAME(copse3d)(option, finalValue1, numStep1, finalValue2, numStep2);
    }

protected:
private:
    void GenerateInputFiles(); // Generate the PSEREAD.NAM for the LST solver

    // PSEREAD parameters
    int m_NSTART        = 1;
    int m_ITSE          = 0;
    int m_INEUTRL       = 0;    // flag for spanwise wavenumber/frequency sweep
    NekDouble m_BETALIN = 0;    // spanwise wavenumber, unit [1/m]
    NekDouble m_FREQLIN = 2000; // frequency of disturbances, unit [1/s]
    int m_NGLOBAL       = 25;   // no. of points to be used in global analysis.
    int m_NLOCAL        = 121;  // no. of points to be used in local analysis.
    int m_ITOLRAY       = 5;    // tolerance in Rayleigh iteration
    std::string m_LABEL = "\'localStabilityAnalysis\'";
    int m_ITSPRES       = 0;   // history output, 0 ,1, 2
    int m_IEIGL         = 1;   // output the eigen function
    NekDouble m_RGAS  = 287.0; // in physical unit rather than the non-dim value
    NekDouble m_GAMMA = 1.4;
    NekDouble m_PRANDTL = 0.72;
};
} // namespace Nektar::FieldUtils

#endif
