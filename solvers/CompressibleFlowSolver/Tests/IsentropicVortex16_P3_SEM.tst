<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=3, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P3_SEM.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P3_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0111668</value>
            <value variable="rhou" tolerance="1e-12">0.0222945</value>
            <value variable="rhov" tolerance="1e-12">0.0202924</value>
            <value variable="E" tolerance="1e-12">0.065215</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0126154</value>
            <value variable="rhou" tolerance="1e-12">0.0278622</value>
            <value variable="rhov" tolerance="1e-12">0.0228415</value>
            <value variable="E" tolerance="1e-12">0.0913806</value>
        </metric>
    </metrics>
</test>


