<?xml version="1.0" encoding="utf-8"?>
<test runs="50">
    <description> 3D Helmholtz/Steady Diffusion Reaction P=5 All Element Types </description>
    <executable>ADRSolver</executable>
    <parameters>CubeAllElements.xml</parameters>
    <files>
        <file description="Session File">CubeAllElements.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.000455961</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00297549</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="2e-1" hostname="42.debian-bullseye-performance-build-and-test">3.61641</value>
        </metric>
    </metrics>
</test>
