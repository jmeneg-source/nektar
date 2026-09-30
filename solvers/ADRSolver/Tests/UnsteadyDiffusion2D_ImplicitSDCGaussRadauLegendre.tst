<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyDiffusion2D_ImplicitSDCGaussRadauLegendre.xml</parameters>
    <processes>1</processes>
    <files>
        <file description="Session File"> UnsteadyDiffusion2D_ImplicitSDCGaussRadauLegendre.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-12">1.78681e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-12">3.01586e-10</value>
        </metric>
    </metrics>
</test>
