**Overview**

This is an example of a laminar channel flow in 2D using quadrilateral elements. The case is run for `NUMMODES=3` or at P=2. 

**Execution**

To execute this test case, type: 

`IncNavierStokesSolver ChannelFlow2D.xml -v`

The -v option is an optional argument which provides a verbose output. The first time you run this example it will generate an optimisation .opt (`ChannelFlow2D.opt`) file which it will reuse on subsequent runs.

**Output**

The run will produce an output file `ChannelFlow2D.fld` and two checkpoint files `ChannelFlow2D_0.chk` and `ChannelFlow2D_1.chk` which can be postprocessed using FieldConvert, i.e.

`FieldConvert ChannelFlow2D.xml ChanneFlow2D.fld ChannelFlow2D.vtu`

_Note if you want to edit the session .xml file you should copy it to another file name otherwise ctest might failed for this test._