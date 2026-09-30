<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=4, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16Implicit_P4.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16Implicit_P4.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-04">1.67364e-03</value>
            <value variable="rhou" tolerance="1e-04">3.21887e-03</value>
            <value variable="rhov" tolerance="1e-04">3.35898e-03</value>
            <value variable="E" tolerance="1e-04">7.08994e-03</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-04">3.63568e-03</value>
            <value variable="rhou" tolerance="1e-04">8.48122e-03</value>
            <value variable="rhov" tolerance="1e-04">8.07722e-03</value>
            <value variable="E" tolerance="1e-03">1.52236e-02</value>
        </metric>
    </metrics>
</test>

