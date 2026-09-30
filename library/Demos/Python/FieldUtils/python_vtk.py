import sys, vtk
from NekPy.FieldUtils import Field, InputModule, OutputModule

field = Field(sys.argv)

# Create input and output modules. Output module deliberately sets an empty
# string, which prevents output from actually being written.
InputModule.Create("xml", field, infile={"xml": sys.argv[1]}).Run()
output_module = OutputModule.Create("vtu", field, outfile="")
output_module.Run()

# From the output module we can now grab a VTK handle.
grid = output_module.GetVtkGrid()

# Grab the bounding box and the centre point
center = grid.GetCenter()
bounds = grid.GetBounds()
ncells = grid.GetNumberOfCells()

# Print out some broad information about the mesh for testing
print("Number of cells:", ncells)
print("xmin:", bounds[0])
print("xmax:", bounds[1])
print("ymin:", bounds[2])
print("ymax:", bounds[3])
print("zmin:", bounds[4])
print("zmax:", bounds[5])
print("cx:", center[0])
print("cy:", center[1])
print("cz:", center[2])
