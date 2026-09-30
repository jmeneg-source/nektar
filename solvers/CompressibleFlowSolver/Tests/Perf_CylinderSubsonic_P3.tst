<?xml version="1.0" encoding="utf-8"?>
<test runs="10">
    <description>Euler, Subsonic Cylinder P=3</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Perf_CylinderSubsonic_P3.xml</parameters>
    <files>
        <file description="Session File">Perf_CylinderSubsonic_P3.xml</file>
        <file description="Restart File">Perf_CylinderSubsonic_P3.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">36.4222</value>
            <value variable="rhou" tolerance="1e-12">3.6391</value>
            <value variable="rhov" tolerance="1e-12">1.29803</value>
            <value variable="E" tolerance="1e-12">4.29634</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.33193</value>
            <value variable="rhou" tolerance="1e-12">0.195761</value>
            <value variable="rhov" tolerance="1e-12">0.114451</value>
            <value variable="E" tolerance="1e-12">0.163198</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="1e0" hostname="42.debian-bullseye-performance-build-and-test">16.5</value>
        </metric>
    </metrics>
</test>
