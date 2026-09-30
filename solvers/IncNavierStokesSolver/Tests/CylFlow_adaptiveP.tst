<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D cylinder flow simulation using adaptive polynomial order</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_adaptiveP.xml</parameters>
    <files>
        <file description="Session File">CylFlow_adaptiveP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-4">21.4162</value>
            <value variable="v" tolerance="1e-6">0.911294</value>
            <value variable="p" tolerance="1e-4">0.8743</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">1.76893</value>
            <value variable="v" tolerance="1e-6">0.837797</value>
            <value variable="p" tolerance="1e-5">1.29058</value>
        </metric>
    </metrics>
</test>

