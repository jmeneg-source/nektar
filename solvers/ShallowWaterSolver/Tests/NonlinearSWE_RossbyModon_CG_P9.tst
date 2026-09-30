<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Rossby modon, CG, P=9</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>NonlinearSWE_RossbyModon_CG_P9.xml</parameters>
    <files>
        <file description="Session File">NonlinearSWE_RossbyModon_CG_P9.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="h" tolerance="1e-8">27.8415</value>
            <value variable="hu" tolerance="1e-8">0.600778</value>
            <value variable="hv" tolerance="1e-8">0.152893</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h" tolerance="1e-8">1.15912</value>
            <value variable="hu" tolerance="1e-8">0.295309</value>
            <value variable="hv" tolerance="1e-8">0.0482965</value>
        </metric>
    </metrics>
</test>


