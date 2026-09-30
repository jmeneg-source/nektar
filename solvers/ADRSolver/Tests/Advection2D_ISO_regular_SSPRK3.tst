<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady advection GLL_LAGRANGE, P=10, homogeneous Dirichlet bcs, regular elements, SSPRK3</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_ISO_regular_SSPRK3.xml</parameters>
    <files>
        <file description="Session File">Advection2D_ISO_regular_SSPRK3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.00114112 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.00522123 </value>
        </metric>
    </metrics>
</test>
