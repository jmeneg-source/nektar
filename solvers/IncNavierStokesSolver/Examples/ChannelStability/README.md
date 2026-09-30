**Overview**

This is an example of performing a direct stability analysis on a Poiseuille flow (flow is parallel to the walls between two infinite plates in 2D) case at Re = 7500. We are interested in computing the leading eigenvalue of the system using the Arnoldi method. The mesh used triangular and quadrilateral elements. The case is run at P=10.

The parameter `EvolutionOperator` must be `Direct` to consider the linearised Navier-Stokes equations. The `Driver` is set up as `ModifiedArnoldi` to solve the eigenproblem.

**Execution**

A base flow is required, which is specified under `Function` named `BaseFlow` in the session file. The base flow file in the folder is `ChannelStability.bse`. The initial guess is specified through the `InitialConditions` function, and the file is `ChannelStability.rst`.

To execute this test case, type: 

`IncNavierStokesSolver ChannelStability.xml -v`

The -v option is an optional argument which provides a verbose output. The first time you run this example it will generate an optimisation .opt (`ChannelStability.opt`) file which it will reuse on subsequent runs.

**Output**

The run will produce an output file `ChannelStability.fld` and two checkpoint files `ChannelStability_0.chk` and `ChannelStability_1.chk` which can be postprocessed using FieldConvert, i.e.

`FieldConvert ChannelStability.xml ChannelStability.fld ChannelStability.vtu`

Alongwith the `.chk` files, the eigenvectors are also dumped in `ChannelStability_eig_0.fld` and `ChannelStability_eig_1.fld` which can be post-processed in a similar manner. The run also produces the eigenvalues and eigenvectors computed at every iteration in the file called `ChannelStability.evl`. 

_Note if you want to edit the session .xml file you should copy it to another file name otherwise ctest might failed for this test._
