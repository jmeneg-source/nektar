###############################################################################
##
## File: Helmholtz.py
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
## Description: Example solver using the EquationSystem Python bindings.
##
###############################################################################

import sys
import numpy as np
from NekPy.LibUtilities import SessionReader, NekError
from NekPy.StdRegions import ConstFactorMap, ConstFactorType
from NekPy.SpatialDomains import MeshGraphIO
from NekPy.SolverUtils import EquationSystem, Filter

class Helmholtz(EquationSystem):
    def __init__(self, session, graph):
        super(Helmholtz, self).__init__(session, graph)

    def InitObject(self, declare_field):
        super(Helmholtz, self).InitObject(declare_field)

        self.lamb = session.GetParameter("Lambda")
        self.factors = ConstFactorMap()
        self.factors[ConstFactorType.FactorLambda] = self.lamb

        # Get the forcing function
        func = self.GetFunction("Forcing")

        # Print this out
        print("Forcing function:")
        for v in session.GetVariables():
            print(" - variable {:s}: {:s}".format(v, func.Describe(v)))

        func.Evaluate(session.GetVariables(), self.fields)

    def DoSolve(self):
        for field in self.GetFields():
            field.coeffs = field.HelmSolve(field.phys, self.factors)
            field.SetPhysState(False)

EquationSystem.Register("Helmholtz", Helmholtz)
            
if __name__ == '__main__':
    session = SessionReader.CreateInstance(sys.argv)
    graph = MeshGraphIO.Read(session)

    # Manually set projection
    session.SetSolverInfo("Projection", "Continuous")

    # Create equation system from the factory.
    eqsys = EquationSystem.Create("Helmholtz", session, graph)
    # Alternative: eqsys = Helmholtz(session, graph)

    # Initialise, print summary, and solve the system
    eqsys.InitObject(True)
    eqsys.DoSolve()

    for i, var in enumerate(session.GetVariables()):
        print("L 2 error (variable {:s}) : {:}".format(var, eqsys.L2Error(i)))
        print("L inf error (variable {:s}) : {:}".format(var, eqsys.LinfError(i)))
