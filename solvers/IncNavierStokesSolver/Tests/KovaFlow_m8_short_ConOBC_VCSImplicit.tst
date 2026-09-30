<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8 with convective link outflow BC and Implicit scheme</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_ConOBC_VCSWeakPress.xml -I SolverType=VCSImplicit</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_ConOBC_VCSWeakPress.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">2.30703e-08</value>
            <value variable="v" tolerance="1e-10">5.73407e-09</value>
            <value variable="p" tolerance="1e-9">1.15756e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">7.08275e-08</value>
            <value variable="v" tolerance="1e-9">4.61323e-08</value>
            <value variable="p" tolerance="1e-8">2.93157e-07</value>
        </metric>
    </metrics>
</test>
