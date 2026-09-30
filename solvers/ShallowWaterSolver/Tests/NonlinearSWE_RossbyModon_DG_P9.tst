<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Rossby modon, DG, P=9</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>NonlinearSWE_RossbyModon_DG_P9.xml</parameters>
    <files>
        <file description="Session File">NonlinearSWE_RossbyModon_DG_P9.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="h" tolerance="1e-12">27.8415</value>
            <value variable="hu" tolerance="1e-12">0.60092</value>
            <value variable="hv" tolerance="1e-12">0.153466</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h" tolerance="1e-12">1.1588</value>
            <value variable="hu" tolerance="1e-12">0.3005</value>
            <value variable="hv" tolerance="1e-12">0.0484094</value>
        </metric>
    </metrics>
</test>


