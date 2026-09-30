import vtk, sys
from NekPy.FieldUtils import Field, OutputModule

field = Field([], force_output=True)
output_num = 0
output_vtk = None

def filter_initialise(fields, time):
    # Generally this is bad form but we're doing this for proof of concept!!
    global field, output_vtk

    field.SetupFromExpList(fields)
    output_vtk = OutputModule.Create("vtu", field, outfile="output_0.vtu")

def filter_update(fields, time):
    global field, output_num, output_vtk

    output_num += 1
    if output_num % 10 != 0:
        return

    field.SetupFromExpList(fields)

    output_vtk.RegisterConfig("outfile", "output_%03d.vtu" % (output_num / 10, ))
    output_vtk.Run()

    # now here you could e.g. extract a slice using VTK and then dump it out
    # to a file instead of dumping whole file, or do whatever postprocessing
    # you would like to

def filter_finalise(fields, time):
    # print out for testing!
    sys.stdout.write("FilterPython_Function: %d\n" % output_num)
    sys.stdout.flush()
