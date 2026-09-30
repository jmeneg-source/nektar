<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Linear stability (Mod. Arnoldi): Channel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChannelStability.xml</parameters>
    <files>
        <file description="Session File">ChannelStability.xml</file>
        <file description="Session File">ChannelStability.bse</file>
        <file description="Session File">ChannelStability.rst</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.00031,0.0349782</value>
            <value tolerance="0.001">1.00031,-0.0349782</value>
        </metric>
    </metrics>
</test>
