**Overview**

This is an example of a laminar channel flow in 3D using quadrilateral elements. This example makes use of a pure spectral discretisation in one direction as the case is geometrically homogeneous in one direction. 

This case reuses the 2D files in `solvers/IncNavierStokesSolver/Examples/ChannelFlow2D`, with a Fourier expansion in the third direction. The parameter `HomModesZ` specifies the number of Fourier modes to use in the homogeneous direction. The `LZ` parameter specifies the physical length of the domain in that direction.

**Execution**

To execute this test case, type: 

`IncNavierStokesSolver ChannelFlowQuasi3D.xml -v`

The -v option is an optional argument which provides a verbose output. The first time you run this example it will generate an optimisation .opt (`ChannelFlowQuasi3D.opt`) file which it will reuse on subsequent runs.

**Output**

The run will produce an output file `ChannelFlowQuasi3D.fld` and two checkpoint files `ChannelFlowQuasi3D_0.chk` and `ChannelFlowQuasi3D_1.chk` which can be postprocessed using FieldConvert, i.e.

`FieldConvert ChannelFlowQuasi3D.xml ChannelFlowQuasi3D.fld ChannelFlowQuasi3D.vtu`

_Note if you want to edit the session .xml file you should copy it to another file name otherwise ctest might failed for this test._
