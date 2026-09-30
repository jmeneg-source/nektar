<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG explicit advection</description>
    <executable>ADRSolver</executable>
    <parameters>--npt 8 PararealDriverUnsteadyAdvection2D.xml</parameters>
    <processes>8</processes>
    <files>
        <file description="Session File"> PararealDriverUnsteadyAdvection2D.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-12">1.84946e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-12">1.40477e-08</value>
        </metric>
    </metrics>
</test>
