**Overview**

This is an example of performing a transient growth stability analysis for a flow past a Backward-Facing Step at Re = 500. The mesh used triangular and quadrilateral elements. The case is run at P=7.

The parameter `EvolutionOperator` must be `TransientGrowth`. The `Driver` is set up as `Arpack` to solve the eigenproblem. This also means that you need the `NEKTAR_USE_ARPACK` turned on while compiling.

**Execution**

A base flow is required, which is specified under `Function` named `BaseFlow` in the session file. The base flow file in the folder is `BackwardFacingStep_TG.bse`. The initial guess is specified through the `InitialConditions` function, and the file is `BackwardFacingStep_TG.rst`. 

To execute this test case, type: 

`IncNavierStokesSolver BackwardFacingStep_TG.xml -v`

The `-v` option is an optional argument which provides a verbose output. The first time you run this example it will generate an optimisation .opt (`BackwardFacingStep_TG.opt`) file which it will reuse on subsequent runs.

**Output**

The run will produce an output file `BackwardFacingStep_TG.fld` and two checkpoint files `BackwardFacingStep_TG_0.chk` and `BackwardFacingStep_TG_1.chk` which can be postprocessed using FieldConvert, i.e.

`FieldConvert BackwardFacingStep_TG.xml BackwardFacingStep_TG.fld BackwardFacingStep_TG.vtu`

Alongwith the `.chk` files, the eigenvectors are also dumped in `BackwardFacingStep_TG_eig_0.fld` and `BackwardFacingStep_TG_eig_1.fld` which can be post-processed in a similar manner. The run also produces the eigenvalues and eigenvectors computed at every iteration in the file called `BackwardFacingStep_TG.evl`. 

_Note if you want to edit the session .xml file you should copy it to another file name otherwise ctest might failed for this test._