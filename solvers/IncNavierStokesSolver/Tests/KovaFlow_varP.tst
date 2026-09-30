<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow variable P</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_varP.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_varP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.38342e-05</value>
            <value variable="v" tolerance="1e-12">6.97199e-06</value>
	    <value variable="p" tolerance="1e-12">0.000950975</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000144107</value>
            <value variable="v" tolerance="1e-12">1.01063e-05</value>
	    <value variable="p" tolerance="1e-12">0.00283351</value>
        </metric>
    </metrics>
</test>
