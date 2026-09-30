<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, tets at P=2 with Gradient Stabilisation </description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_Tet_GJP.xml</parameters>
    <files>
        <file description="Session File">Advection3D_Tet_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">0.0983789</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0.532112</value>
        </metric>
    </metrics>
</test>

