<?xml version="1.0" encoding="utf-8"?>
<test runs="30">
    <description>Double Bifurcation, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>Perf_TwoBifurcations.xml</parameters>
    <files>
        <file description="Session File">Perf_TwoBifurcations.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">95.5631</value>
            <value variable="u" tolerance="1e-12">18.8865</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">6.40844</value>
            <value variable="u" tolerance="1e-12">6.05588</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="6e-1" hostname="42.debian-bullseye-performance-build-and-test">6.05743</value>
        </metric>
    </metrics>
</test>
