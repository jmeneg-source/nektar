<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG explicit advection</description>
    <executable>ADRSolver</executable>
    <parameters>--npt 8 PFASSTDriverUnsteadyAdvection2D.xml</parameters>
    <processes>8</processes>
    <files>
        <file description="Session File"> PFASSTDriverUnsteadyAdvection2D.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-12">1.04397e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-12">2.75569e-09</value>
        </metric>
    </metrics>
</test>
