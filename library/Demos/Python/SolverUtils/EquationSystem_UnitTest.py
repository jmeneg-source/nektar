###############################################################################
##
## File: EquationSystem_UnitTest.py
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
## Description: Unit tests for the EquationSystem class.
##
###############################################################################

import math, sys, io, os, unittest, argparse
import numpy as np

from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraphIO
from NekPy.SolverUtils import EquationSystem
from UnitTestUtils import SuppressStream

# Create an EquationSystem class to test registration of classes from Python.
class TestingEquationSystem(EquationSystem):
    def __init__(self, session, graph):
        # Call super's constructor.
        super(TestingEquationSystem, self).__init__(session, graph)
        self.init_called = True
        self.InitObject(True)

    def InitObject(self, declare_field):
        super(TestingEquationSystem, self).InitObject(declare_field)
        self.initobject_called = True

    def DoInitialise(self, dump_initial_conditions):
        super(TestingEquationSystem, self).DoInitialise(dump_initial_conditions)
        self.initialise_called = True

    def DoSolve(self):
        self.solve_called = True

    def EvaluateExactSolution(self, field, time):
        tmp = super(TestingEquationSystem, self).EvaluateExactSolution(field, time)
        self.exact_called = True
        return tmp

# Register the EquationSystem class with the factory.
EquationSystem.Register("Testing", TestingEquationSystem)

class TestEquationSystem(unittest.TestCase):
    def setUp(self):
        # Create session and meshgraph
        self.session = SessionReader.CreateInstance(['EquationSystem_UnitTest', nektar_filename])
        self.graph = MeshGraphIO.Read(self.session)

        assert self.graph.GetSpaceDimension() == 2, "Input file should be 2D"

    def testInit(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        self.assertEqual(eqsys.init_called, True)
        self.assertEqual(eqsys.initobject_called, True)

    def testDoInitialise(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        eqsys.DoInitialise(True)
        self.assertEqual(eqsys.initialise_called, True)
        
    def testDoSolve(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        eqsys.DoSolve()
        self.assertEqual(eqsys.solve_called, True)

    def testEvaluateExactSolution(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)

        for v in range(eqsys.GetNvariables()):
            exact_sol = eqsys.EvaluateExactSolution(v, 0)
            self.assertEqual(eqsys.exact_called, True)

            # Try to test against function from session
            func = self.session.GetFunction("ExactSolution", v)
            x, y = eqsys.GetFields()[v].GetCoords()
            exact_eval = func.Evaluate(x, y, np.zeros(len(x)), 0.0)

            # These should be equal!
            self.assertEqual(math.isclose(np.max(np.abs(exact_sol - exact_eval)), 0.0), True)

            # Try to test against function from SessionFunction
            exact_sessionfunc = eqsys.GetFunction("ExactSolution").Evaluate(
                self.session.GetVariables()[v:v+1], eqsys.fields[v:v+1])
            self.assertEqual(math.isclose(np.max(np.abs(exact_sol - exact_sessionfunc[0].phys)), 0.0), True)

    def testGetFields(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        variables = self.session.GetVariables()
        fields = eqsys.GetFields()
        self.assertEqual(len(variables), len(fields))

        # Make sure that fields was initialised properly.
        self.assertGreater(len(variables), 0)
        self.assertGreater(fields[0].GetNpoints(), 0)

    def testGetNvariables(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        self.assertEqual(len(self.session.GetVariables()), eqsys.GetNvariables())

    def testGetNpoints(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        self.assertEqual(eqsys.GetNpoints(), eqsys.fields[0].GetNpoints())

    def testTime(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        eqsys.SetTime(1.0)
        self.assertEqual(eqsys.GetTime(), 1.0)

        eqsys.time = 2.0
        self.assertEqual(eqsys.time, 2.0)

    def testTimeStep(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)
        eqsys.SetTimeStep(1.0)
        self.assertEqual(eqsys.GetTimeStep(), 1.0)

        eqsys.timestep = 2.0
        self.assertEqual(eqsys.timestep, 2.0)

    def testSteps(self):
        eqsys = EquationSystem.Create("Testing", self.session, self.graph)

        eqsys.SetSteps(100)
        self.assertEqual(eqsys.GetSteps(), 100)

        eqsys.steps = 200
        self.assertEqual(eqsys.steps, 200)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('unittest_args', nargs='*')

    args, unknown = parser.parse_known_args()
    nektar_filename = args.filename
    sys.argv[1:] = args.unittest_args + unknown
    unittest.main()
