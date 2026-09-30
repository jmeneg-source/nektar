<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Kovasznay flow on domain with all elements and some non-regular elements with Semi-Implicit GJP stabilisation</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Kovasznay3D_AllElmts_GJP.xml</parameters>
    <files>
        <file description="Session File">Kovasznay3D_AllElmts_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">0.0189978</value>
            <value variable="v" tolerance="1e-7">0.0043984</value>
            <value variable="w" tolerance="1e-7">0.00213548</value>
            <value variable="p" tolerance="1e-7">0.0266785</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0.0465685</value>
            <value variable="v" tolerance="1e-7">0.0154147</value>
            <value variable="w" tolerance="1e-7">0.0120201</value>
            <value variable="p" tolerance="1e-7">0.101036</value>
        </metric>
    </metrics>
</test>
