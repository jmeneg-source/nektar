Overview
--------
In this example, it will be demonstrated how the Advection equation can be solved on a
one-dimensional domain.


The input for this example is given in the example file Advection1D_WeakDG_GLL_LAGRANGE.xml

Execution
---------
To execute this test case run: \
ADRSolver Advection1D_WeakDG_GLL_LAGRANGE.xml

Output
------
To visualise the output, we can convert it into either Tecplot or VTK formats: \
FieldConvert Advection1D_WeakDG_GLL_LAGRANGE.xml Advection1D_WeakDG_GLL_LAGRANGE.fld Advection1D_WeakDG_GLL_LAGRANGE.dat \
FieldConvert Advection1D_WeakDG_GLL_LAGRANGE.xml Advection1D_WeakDG_GLL_LAGRANGE.fld Advection1D_WeakDG_GLL_LAGRANGE.vtu 

Note if you want to edit the Advection1D_WeakDG_GLL_LAGRANGE.xml file you should copy it to another file name otherwise ctest might failed for this test.
