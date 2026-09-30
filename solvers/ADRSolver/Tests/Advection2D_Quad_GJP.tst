<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady DG advection, quad elmts P=4 with Gradient Stabilisation (two forced curved elmt) </description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_Quad_GJP.xml</parameters>
    <files>
        <file description="Session File">Advection2D_Quad_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">0.00349066</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0.0383412</value>
        </metric>
    </metrics>
</test>
