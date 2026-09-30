**Overview**

This is an example of a 2-DoF self-propelled flapping airfoil.

**Execution**

To execute this test case, type: 

`IncNavierStokesSolver SelfPropelledAirfoil.xml -v`

The -v option is an optional argument which provides a verbose output. The first time you run this example it will generate an optimisation .opt (`SelfPropelledAirfoil.opt`) file which it will reuse on subsequent runs.

**Output**

The run will produce an output file `SelfPropelledAirfoil.fld` and two checkpoint files `SelfPropelledAirfoil_0.chk` and `SelfPropelledAirfoil_1.chk` which can be postprocessed using FieldConvert, i.e.

`FieldConvert SelfPropelledAirfoil.xml SelfPropelledAirfoil.fld SelfPropelledAirfoil.plt`

_Note if you want to edit the session .xml file you should copy it to another file name otherwise ctest might failed for this test._
