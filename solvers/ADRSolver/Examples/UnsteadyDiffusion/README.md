Overview
--------
In this example, it will be demonstrated how the Unsteady Reaction-Diffusion equation can be solved on a
two-dimensional domain.


The input for this example is given in the example file ImDiffusion_VarCoeff.xml

Execution
---------
To execute this test case run: \
ADRSolver ImDiffusion_VarCoeff.xml

Output
------
To visualise the output, we can convert it into either Tecplot or VTK formats: \
FieldConvert ImDiffusion_VarCoeff.xml ImDiffusion_VarCoeff.fld ImDiffusion_VarCoeff.dat \
FieldConvert ImDiffusion_VarCoeff.xml ImDiffusion_VarCoeff.fld ImDiffusion_VarCoeff.vtu 

Note if you want to edit the ImDiffusion_VarCoeff.xml file you should copy it to another file name otherwise ctest might failed for this test.
