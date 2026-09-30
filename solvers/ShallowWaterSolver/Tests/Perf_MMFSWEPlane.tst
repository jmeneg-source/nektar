<?xml version="1.0" encoding="utf-8"?>
<test runs="20">
    <description>MMF SWE solver, DG, P=4</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>Perf_MMFSWEPlane.xml</parameters>
    <files>
        <file description="Session File">Perf_MMFSWEPlane.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-5">0.0303929</value>
            <value variable="u" tolerance="1e-5">0.0343446</value>
            <value variable="v" tolerance="1e-5">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-5">0.032063</value>
            <value variable="u" tolerance="1e-5">0.041687</value>
            <value variable="v" tolerance="1e-5">0.0</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="5e-1" hostname="42.debian-bullseye-performance-build-and-test">7.3</value>
        </metric>
    </metrics>
</test>


