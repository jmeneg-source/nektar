<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=4, WeakDG</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_WeakDG_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_WeakDG_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-10">0.000119397</value>
            <value variable="rhou" tolerance="1e-10">0.000401715</value>
            <value variable="rhov" tolerance="1e-10">0.000387776</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-09">0.00113732</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-10">0.000212639</value>
            <value variable="rhou" tolerance="1e-10">0.000904271</value>
            <value variable="rhov" tolerance="1e-10">0.000819096</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-10">0.00240946</value>
        </metric>
    </metrics>
</test>


