<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion </description>
    <executable python="true">UnsteadyDiffusion.py</executable>
    <parameters> ImDiffusion_m6.xml</parameters>
    <files>
        <file description="Session File">../../../../../solvers/ADRSolver/Tests/ImDiffusion_m6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08"> 0.003957190532100898 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08"> 0.03203139672729776 </value>
        </metric>
    </metrics>
</test>
