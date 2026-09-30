**Overview**
This is an example of a laminar channel flow using all types of 3D elemets (Hex, Prism, Pyramid and Tetrahedron). 

The mesh file is in Hdf5 format and so you need to have compiled Nektar with the `NEKTAR_USE_HDF5` and `NEKTAR_USE_MPI` options. To view the file you can run

`h5dump ChannelFlow3D.nekg` 

Or you can convert the input into XML ASCII format using 

`NekMesh ChannelFlow3D.xml  ChannelFlow3D_xml.xml:xml:uncompress`

**Execution**

To execute this test case type: 

`IncNavierStokesSolver ChannelFlow3D.xml -v`

or to run in parallel 

`mpirun -n 2 IncNavierStokesSolver ChannelFlow3D.xml -v`  

(Note this will automatically select an iterative solver)

The `-v` option is an optional argument which provides a verbose output. The first time you run this example it will generate an optimisation `.opt` file which it will reuse on subseuquent runs.  

_Note if you want to edit the session .xml file you should copy it to another file name otherwise ctest might failed for this test._  

**Output**

The run will produce an output file `ChannelFlow3D.fld` and two checkpoint files `ChannelFlow3D_0.chk` and `ChannelFlow3D_1.chk` which can be postprocessed using FieldConvert, i.e.

`FieldConvert ChannelFlow3D.xml ChanneFlow3D.fld ChannelFlow3D.vtu`
