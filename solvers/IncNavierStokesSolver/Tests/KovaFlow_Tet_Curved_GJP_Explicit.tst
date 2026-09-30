<?xml version="1.0" encoding="utf-8"?>
<test>
3   <description>3D Kovasznay flow on curved Tetrahedrons with explicit GJP stabilisation</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_Tet_Curved.xml -I GJPStabilisation=Explicit -I GJPNormalVelocity=True</parameters>
    <files>
        <file description="Session File">KovaFlow_Tet_Curved.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.000539915</value>
            <value variable="v" tolerance="1e-8">0.000245073</value>
            <value variable="w" tolerance="1e-10">0.000163579</value>
            <value variable="p" tolerance="1e-8">0.00400963</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">0.00209877</value>
            <value variable="v" tolerance="1e-8">0.000873068</value>
            <value variable="w" tolerance="1e-8">0.000513662</value>
            <value variable="p" tolerance="1e-8">0.0121881</value>
        </metric>
    </metrics>
</test>
