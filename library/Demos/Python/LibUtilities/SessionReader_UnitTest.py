###############################################################################
##
## File: SessionReader_UnitTest.py
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
## Description: Unit tests for the SessionReader class.
##
###############################################################################

from NekPy.LibUtilities import SessionReader
import unittest

class TestInterpreter(unittest.TestCase):
    def setUp(self):
        session_name = ["NekPy_SessionReader_UnitTest.py", "newsquare_2x2.xml"]
        self.session = SessionReader.CreateInstance(session_name)
        self.session.InitSession()

    # Check that session parameters can be retrieved in a dict
    def testReadParamsMap(self):
        params = self.session.GetParameters()
        # Check that a dict is returned
        self.assertTrue(
            isinstance(params, dict),
            "GetParameters returned object of type {}".format(type(params)),
        )
        # Check that dict isn't empty
        self.assertTrue(params, "Expected non-empty parameter dict")

    # Check that session variables can be retrieved in a list
    def testReadVarsList(self):
        vars = self.session.GetVariables()
        # Check that a list is returned
        self.assertTrue(
            isinstance(vars, list),
            "GetVariables returned object of type {}".format(type(vars)),
        )
        # Check that list isn't empty
        self.assertTrue(vars, "Expected non-empty list of variables")


if __name__ == "__main__":
    unittest.main()
