<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=3 with explicit GJP and Normal velocity using Semi-Implicit VCS</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m3_GJP.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m3_GJP.xml</file>
        <file description="Session File">KovaFlow_m3.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.342624</value>
            <value variable="v" tolerance="1e-6">0.0260572</value>
            <value variable="p" tolerance="1e-6">0.321985</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.322039</value>
            <value variable="v" tolerance="1e-6">0.0257921</value>
            <value variable="p" tolerance="1e-6">0.500148</value>
        </metric>
    </metrics>
</test>
