<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8 Weak Pressure VSS</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_ConOBC_VCSWeakPress.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_ConOBC_VCSWeakPress.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.30721e-08</value>
            <value variable="v" tolerance="1e-11">5.73276e-09</value>
	    <value variable="p" tolerance="1e-11">1.12259e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">7.08359e-08</value>
            <value variable="v" tolerance="1e-11">4.61745e-08</value>
	    <value variable="p" tolerance="1e-10">2.9248e-07</value>
        </metric>
    </metrics>
</test>
