<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.09988e-06</value>
            <value variable="v" tolerance="1e-11">1.41301e-06</value>
            <value variable="w" tolerance="1e-11">1.23744e-06</value>
	    <value variable="p" tolerance="1e-10">2.96068e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">3.02806e-06</value>
            <value variable="v" tolerance="1e-11">2.22791e-06</value>
            <value variable="w" tolerance="1e-11">1.98813e-06</value>
	    <value variable="p" tolerance="1e-10">6.84108e-05</value>
        </metric>
    </metrics>
</test>
