# -f -e --no-equispaced -m pointdatatofld outflow.pts outflow.xml outflow.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, force_output=True, error=True, no_equispaced=True)

InputModule.Create("pts",  field, "outflow.pts").Run()
InputModule.Create("xml",  field, "outflow.xml").Run()
ProcessModule.Create("pointdatatofld", field).Run()
OutputModule.Create("fld", field, "outflow.fld").Run()
