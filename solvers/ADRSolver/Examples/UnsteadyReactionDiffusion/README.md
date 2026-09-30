Overview
--------
In this example, it will be demonstrated how the Unsteady Reaction-Diffusion equation can be solved on a
two-dimensional domain.


The input for this example is given in the example file ReactionDiffusion2D.xml

Execution
---------
To execute this test case run: \
ADRSolver ReactionDiffusion2D.xml

Output
------
To visualise the output, we can convert it into either Tecplot or VTK formats: \
FieldConvert ReactionDiffusion2D.xml ReactionDiffusion2D.fld ReactionDiffusion2D.dat \
FieldConvert ReactionDiffusion2D.xml ReactionDiffusion2D.fld ReactionDiffusion2D.vtu 

Note if you want to edit the ReactionDiffusion2D.xml file you should copy it to another file name otherwise ctest might failed for this test.
