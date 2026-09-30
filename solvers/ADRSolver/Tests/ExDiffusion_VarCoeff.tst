<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D unsteady DG explicit diffusion, variable coeffs.</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_VarCoeff.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_VarCoeff.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">7.51811e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.0606e-07</value>
        </metric>
    </metrics>
</test>
