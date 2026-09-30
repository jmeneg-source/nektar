# -f -r 0,0.5,0,0.5,0,0.5 cube_prismhex.xml out.xml
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, force_output=True)

InputModule.Create("xml", field,  "cube_prismhex.xml",
                   range="0,0.5,0,0.5,0,0.5").Run()
OutputModule.Create("xml", field, "out.xml").Run()
