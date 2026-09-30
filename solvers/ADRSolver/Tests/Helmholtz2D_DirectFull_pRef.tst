<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D Helmholtz/Steady Diffusion with local p-refinement </description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz2D_DirectFull_pRef.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_DirectFull_pRef.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8"> 0.000126986 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">  0.000262014 </value>
        </metric>
    </metrics>
</test>
