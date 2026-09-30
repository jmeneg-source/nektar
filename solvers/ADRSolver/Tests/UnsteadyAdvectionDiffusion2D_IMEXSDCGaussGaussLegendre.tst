<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG explicit-implicit advection-diffusion</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion2D_IMEXSDCGaussGaussLegendre.xml</parameters>
    <processes>1</processes>
    <files>
        <file description="Session File"> UnsteadyAdvectionDiffusion2D_IMEXSDCGaussGaussLegendre.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-12">2.71598e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-12">1.5227e-08</value>
        </metric>
    </metrics>
</test>
