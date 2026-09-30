<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes (MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM_Deal.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM_Deal.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.17793e-06</value>
            <value variable="v" tolerance="1e-11">1.34849e-06</value>
            <value variable="w" tolerance="1e-11">1.18944e-06</value>
	    <value variable="p" tolerance="1e-10">2.89563e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">3.08536e-06</value>
            <value variable="v" tolerance="1e-11">2.22781e-06</value>
            <value variable="w" tolerance="1e-11">1.84095e-06</value>
	    <value variable="p" tolerance="1e-11">6.83072e-05</value>
        </metric>
    </metrics>
</test>
