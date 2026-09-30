import vtk, sys

from NekPy.SolverUtils import Filter
from NekPy.FieldUtils import Field, OutputModule

class TestVTKFilter(Filter):
    def __init__(self, session, eqsys, params):
        # Call super's constructor.
        super(TestVTKFilter, self).__init__(session, eqsys)

        self.field = Field([], force_output=True)
        self.num = 0
        self.module = OutputModule.Create(
            "vtu", self.field, outfile="output_0.vtu")

    def Initialise(self, fields, time):
        self.field.SetupFromExpList(fields)

    def Update(self, fields, time):
        self.num += 1
        if self.num % 10 != 0:
            return

        self.field.SetupFromExpList(fields)
        self.module.RegisterConfig("outfile", "output_%03d.vtu" % (self.num / 10, ))
        self.module.Run()

    def Finalise(self, fields, time):
        # Print for testing purposes
        sys.stdout.write("TestVTKFilter: %d\n" % self.num)
        sys.stdout.flush()

    def IsTimeDependent(self):
        return True

# Register filter with factory.
Filter.Register("TestVTKFilter", TestVTKFilter)
