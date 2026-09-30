**Overview**
This example is to solve the 2D Kovasznay flow at Reynolds number Re=40. The case is run for `NUMMODES=3` or at P=2. This case uses steady boundary conditions for velocity and pressure.

**Execution**

To execute this test case, type: 

`IncNavierStokesSolver KovasznayFlow2D.xml -v`

The `-v` option is an optional argument which provides a verbose output. The first time you run this example it will generate an optimisation .opt (`KovasznayFlow2D.opt`) file which it will reuse on subsequent runs.
 
_Note if you want to edit the session .xml file you should copy it to another file name otherwise ctest might failed for this test._  

**Output**

The run will produce an output file `KovasznayFlow2D.fld` and two checkpoint files `KovasznayFlow2D_0.chk` and `KovasznayFlow2D_1.chk` which can be postprocessed using FieldConvert, i.e.

`FieldConvert KovasznayFlow2D.xml KovasznayFlow2D.fld KovasznayFlow2D.vtu`
