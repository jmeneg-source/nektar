###############################################################################
##
## File: Filter_UnitTest.py
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
## Description: Unit tests for the Filter class.
##
###############################################################################

import sys, io, os, unittest, argparse
from NekPy.LibUtilities import SessionReader, NekError
from NekPy.SpatialDomains import MeshGraph, MeshGraphIO
from NekPy.SolverUtils import EquationSystem, Filter
from UnitTestUtils import SuppressStream

# Create a filter class to test registration of classes from Python.
class TestingFilter(Filter):
    def __init__(self, session, eqsys, params):
        # Call super's constructor.
        super(TestingFilter, self).__init__(session, eqsys)

        # Check we've been passed a parameter.
        assert params['param'] == 'testing'
        self.param = params['param']

    def Initialise(self, expLists, time):
        # Counter to check number of times Update function is called.
        self.num = 0

    def Update(self, expLists, time):
        # Update the counter
        self.num += 1

    def Finalise(self, expLists, time):
        # Assert the number is correct
        assert self.num == 10

    def IsTimeDependent(self):
        return True

class TestFilter(unittest.TestCase):
    def setUp(self):
        # Create session and meshgraph
        self.session = SessionReader.CreateInstance(['filtertest', nektar_filename])
        self.graph = MeshGraphIO.Read(self.session)

        # Create our 'dummy' equationsystem
        self.eqsys = EquationSystem.Create("Dummy", self.session, self.graph)

    def testCreateFilter(self):
        # Create a simple checkpoint filter. Parameters are passed as Python
        # keyword arguments.
        f = Filter.Create("Checkpoint", self.session, self.eqsys,
                          OutputFrequency=1, OutputFile="testing")
        self.assertEqual(f.IsTimeDependent(), True)

    def testUnsuppliedConfig(self):
        try:
            f = Filter.Create("Checkpoint", self.session, self.eqsys)
        except NekError:
            pass

    def testGetFields(self):
        # Get list of fields from EquationSystem.
        f = Filter.Create("Checkpoint", self.session, self.eqsys,
                          OutputFrequency=1, OutputFile="testing")
        fields = self.eqsys.GetFields()
        self.assertEqual(len(fields) > 0, True)

    def testRunFilter(self):
        f = Filter.Create("Checkpoint", self.session, self.eqsys,
                          OutputFrequency=1, OutputFile="testing")

        with SuppressStream(sys.stdout):
            # Initialise the filter
            fields = self.eqsys.GetFields()
            f.Initialise(fields, 0.0)

            time = 0.0
            for i in range(10):
                time += 1.0
                f.Update(fields, time)

            # Finalise
            f.Finalise(fields, time)

    def testInheritance(self):
        # Register the filter and create
        Filter.Register("test", TestingFilter)
        f2 = Filter.Create("test", self.session, self.eqsys, param='testing')

        fields = self.eqsys.GetFields()
        f2.Initialise(fields, 0.0)

        for i in range(10):
            f2.Update(fields, 0.0)

        f2.Finalise(fields, 0.0)
        self.assertEqual(f2.num, 10)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('unittest_args', nargs='*')

    args, unknown = parser.parse_known_args()
    nektar_filename = args.filename
    sys.argv[1:] = args.unittest_args + unknown
    unittest.main()
