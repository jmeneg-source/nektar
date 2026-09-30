<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG explicit advection</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvection2D_ExplicitSDCGaussLobattoLegendre.xml</parameters>
    <processes>1</processes>
    <files>
        <file description="Session File"> UnsteadyAdvection2D_ExplicitSDCGaussLobattoLegendre.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-13">7.3708e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-13">6.2686e-14</value>
        </metric>
    </metrics>
</test>
