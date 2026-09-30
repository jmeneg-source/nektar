<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG explicit advection</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvection2D_ExplicitSDCGaussRadauLegendre.xml</parameters>
    <processes>1</processes>
    <files>
        <file description="Session File"> UnsteadyAdvection2D_ExplicitSDCGaussRadauLegendre.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-14">6.65957e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-14">7.60936e-14</value>
        </metric>
    </metrics>
</test>
