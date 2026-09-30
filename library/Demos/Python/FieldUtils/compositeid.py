# -f -e -m addcompositeid compositeid.xml compositeid.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, force_output=True, error=True)

InputModule.Create("xml",  field, "compositeid.xml").Run()
ProcessModule.Create("addcompositeid", field).Run()
OutputModule.Create("fld", field, "compositeid.fld").Run()
