<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=8, 16 Fourier modes, using explicit mapping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P8_16modes_Mapping-explicit.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P8_16modes_Mapping-explicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">2.07269e-05</value>
            <value variable="v" tolerance="2e-13">1.77651e-08</value>
            <value variable="w" tolerance="2e-12">1.53263e-07</value>
	    <value variable="p" tolerance="1e-10">1.5383e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">1.64559e-05</value>
            <value variable="v" tolerance="2e-12">2.01647e-07</value>
            <value variable="w" tolerance="2e-12">3.49113e-07</value>
	    <value variable="p" tolerance="1e-10">8.01198e-06</value>
        </metric>
    </metrics>
</test>
