<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=4-9, 20 Fourier modes in parallel with HDF5 output - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM_hdf5.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM_hdf5.xml</file>
        <file description="Restart File">KovaFlow_3DH1D_P5_20modes_MVM_hdf5.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">1.55744e-05</value>
            <value variable="v" tolerance="1e-11">9.39817e-06</value>
            <value variable="w" tolerance="1e-11">8.12806e-06</value>
            <value variable="p" tolerance="1e-10">9.12185e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">3.3134e-05</value>
            <value variable="v" tolerance="1e-10">2.51101e-05</value>
            <value variable="w" tolerance="1e-10">1.70084e-05</value>
            <value variable="p" tolerance="1e-09">0.000138221</value>
        </metric>
    </metrics>
</test>
