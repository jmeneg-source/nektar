<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8 Weak Pressure VSS</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_HOBC_VCSWeakPress.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_HOBC_VCSWeakPress.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.30589e-08</value>
            <value variable="v" tolerance="1e-11">5.67503e-09</value>
	    <value variable="p" tolerance="1e-11">1.07923e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">7.07099e-08</value>
            <value variable="v" tolerance="1e-11">4.59269e-08</value>
	    <value variable="p" tolerance="1e-10">2.92718e-007</value>
        </metric>
    </metrics>
</test>
