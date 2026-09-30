Examples
========
This folder has brief examples for using the IncNavierStokesSolver for 2D and 3D flow cases. The examples are soft linked from the `solvers/IncNavierStokesSolver/Tests` folder, and should not be modified here in order to keep the testing pipeline working.

Laminar Channel Flow 2D
----------
Flow through a 2D channel at Reynolds number 1 with fixed boundary conditions. Folder: `ChannelFlow2D`

Laminar Channel Flow 3D
---------
Flow through a 3D channel at Reynolds number 1 with fixed boundary conditions. This example uses a 3D mesh in an HDF5 format. The mesh contains hexahedrons, prisms, tetrahedrons and pyramid elements. Folder: `ChannelFlow3D`.

Laminar Channel Flow Quasi-3D
---------
Flow through a 3D channel at Reynolds number 1 with fixed boundary conditions. This example makes use of a pure spectral discretisation in one direction as the case is geometrically homogeneous in one direction. Folder: `ChannelFlowQuasi3D`.

2D Direct Stability Analysis of the Channel Flow
--------------
Perform a direct stability analysis on a Poiseuille flow (flow is parallel to the walls between two infinite plates in 2D) case at Re = 7500. Folder: `ChannelStability`.

2D Adjoint stability analysis of the channel flow
--------------
Perform a adjoint stability analysis on a Poiseuille flow (flow is parallel to the walls between two infinite plates in 2D) case at Re = 7500. Folder: `ChannelStabilityAdjoint`.

Kovasznay Flow 2D
--------------
Solving the 2D Kovasznay flow at Reynolds number Re=40. Folder: `KovasznayFlow2D`.

2D Transient Growth Analysis of a flow past a Backward-Facing Step
--------------
Perform a transient growth stability analysis. This will help in understanding the effects of separation caused by abrupt changes in geometry for a flow past a Backward-Facing Step at Re = 500. Folder: `BackwardFacingStep_TG`.

Self-propelled flapping airfoil
--------------
Two-degree-of-freedom (2-DoF) self-propelled flapping airfoil in an incompressible viscous flow case at Ref=310. Folder:
`SelfPropelledAirfoil`.
