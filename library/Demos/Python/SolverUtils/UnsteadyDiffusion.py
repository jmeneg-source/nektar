###############################################################################
##
## File: UnsteadyDiffusion.py
##
## For more information, please see: http://www.nektar.info
##
## The MIT License
##
## Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
## Department of Aeronautics, Imperial College London (UK), and Scientific
## Computing and Imaging Institute, University of Utah (USA).
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included
## in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
## OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Description: Example solver using the UnsteadySystem Python bindings.
##
###############################################################################

import sys
import numpy as np
from NekPy.LibUtilities import SessionReader
from NekPy.StdRegions import ConstFactorMap, ConstFactorType
from NekPy.SpatialDomains import MeshGraphIO
from NekPy.SolverUtils import EquationSystem, UnsteadySystem

class Diffusion(UnsteadySystem):
    def __init__(self, session, graph):
        UnsteadySystem.__init__(self, session, graph)

    def InitObject(self, declare_field):
        UnsteadySystem.InitObject(self, declare_field)

        # Load diffusion parameter from session file
        if session.DefinesParameter("epsilon"):
            self.epsilon = session.GetParameter("epsilon")
        else:
            self.epsilon = 1.0

        # Create a factors map
        self.factors = ConstFactorMap()

        # Define implicit solution function
        self.ode.DefineImplicitSolve(self.DoImplicitSolve)

    def DoImplicitSolve(self, fields, time, l):
        outfields = []

        lamb = 1.0 / self.epsilon / l
        self.factors[ConstFactorType.FactorLambda] = lamb
        self.factors[ConstFactorType.FactorTau] = 1.0

        for i, f in enumerate(fields):
            self.fields[i].coeffs = self.fields[i].HelmSolve(-lamb * f, self.factors)
            outfields.append(self.fields[i].BwdTrans(self.fields[i].coeffs))

        return outfields

EquationSystem.Register("Diffusion", Diffusion)
            
if __name__ == '__main__':
    session = SessionReader.CreateInstance(sys.argv)
    graph = MeshGraphIO.Read(session)

    # Create equation system from the factory.
    eqsys = EquationSystem.Create("Diffusion", session, graph)

    # Initialise, print summary, and solve the system
    eqsys.InitObject(True)
    eqsys.DoInitialise(True)
    eqsys.DoSolve()
    eqsys.WriteFld("test.fld")

    for i, var in enumerate(session.GetVariables()):
        print("L 2 error (variable {:s}) : {:}".format(var, eqsys.L2Error(i)))
        print("L inf error (variable {:s}) : {:}".format(var, eqsys.LinfError(i)))
