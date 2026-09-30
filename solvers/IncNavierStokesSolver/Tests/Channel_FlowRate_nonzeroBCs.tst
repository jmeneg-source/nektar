<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=4 with flow rate and non-zero BCs</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Channel_FlowRate_nonzeroBCs.xml</parameters>
    <files>
            <file description="Session File">Channel_FlowRate_nonzeroBCs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.146542</value>
            <value variable="v" tolerance="1e-6">0.0651023</value>
            <value variable="p" tolerance="1e-6">0.299874</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.5</value>
            <value variable="v" tolerance="1e-6">0.0811353</value>
            <value variable="p" tolerance="1e-6">0.451181</value>
        </metric>
    </metrics>
</test>


