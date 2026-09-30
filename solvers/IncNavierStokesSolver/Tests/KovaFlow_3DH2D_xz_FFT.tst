<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 2D, flow in xz, FFT</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH2D_xz_FFT.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH2D_xz_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.02685e-07</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">1.15255e-08</value>
	    <value variable="p" tolerance="1e-9">3.85049e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.67617e-07</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">1.91761e-08</value>
	    <value variable="p" tolerance="1e-9">4.77718e-06</value>
        </metric>
    </metrics>
</test>
