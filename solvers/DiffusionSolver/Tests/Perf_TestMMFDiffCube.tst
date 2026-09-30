<?xml version="1.0" encoding="utf-8"?>
<test runs="10">
    <description> 2D MMF implicit diffusion </description>
    <executable>MMFDiffusion</executable>
    <parameters>Perf_TestMMFDiffCube.xml</parameters>
    <files>
        <file description="Session File">Perf_TestMMFDiffCube.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">0.000339924</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">0.000529421</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="2e0" hostname="42.debian-bullseye-performance-build-and-test">33.6845</value>
        </metric>
    </metrics>
</test>
