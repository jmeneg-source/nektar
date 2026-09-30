Overview
--------
In this example, it will be demonstrated how the Helmholtz equation can be solved on a
two-dimensional domain

Execute
-------
To execute this run case go to the path stated within Helmholtz2D_modal.xml and run: \
ADRSolver Helmholtz2D_modal.xml

Output
------
To visualise the output, we can convert it into either Tecplot or VTK formats: \
FieldConvert Helmholtz2D_modal.xml Helmholtz2D_modal.fld Helmholtz2D_modal.dat \
FieldConvert Helmholtz2D_modal.xml Helmholtz2D_modal.fld Helmholtz2D_modal.vtu

Note if you want to edit the Helmholtz2D_modal.xml file you should copy it to another file name otherwise ctest might failed for this test.
