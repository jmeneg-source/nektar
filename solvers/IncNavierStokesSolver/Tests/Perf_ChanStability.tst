<?xml version="1.0" encoding="utf-8" ?>
<test runs="20">
    <description>Linear stability (Mod. Arnoldi): Channel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Perf_ChanStability.xml</parameters>
    <files>
        <file description="Session File">Perf_ChanStability.xml</file>
        <file description="Session File">Perf_ChanStability.bse</file>
        <file description="Session File">Perf_ChanStability.rst</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.00013,0.0149907</value>
            <value tolerance="0.001">1.00013,-0.0149907</value>
        </metric>
        <metric type="ExecutionTime" id="1">
            <regex>^.*Execute\s*(\d+\.?\d*).*</regex>
            <value tolerance="8e-1" hostname="42.debian-bullseye-performance-build-and-test">9.3</value>
        </metric>
    </metrics>
</test>
