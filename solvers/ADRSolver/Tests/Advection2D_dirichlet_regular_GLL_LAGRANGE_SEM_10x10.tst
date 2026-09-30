<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady advection GLL_LAGRANGE_SEM, P=3, Dirichlet bcs, regular elements</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_dirichlet_regular_GLL_LAGRANGE_SEM_10x10.xml</parameters>
    <files>
        <file description="Session File">Advection2D_dirichlet_regular_GLL_LAGRANGE_SEM_10x10.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.000222548 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.000814097 </value>
        </metric>
    </metrics>
</test>
