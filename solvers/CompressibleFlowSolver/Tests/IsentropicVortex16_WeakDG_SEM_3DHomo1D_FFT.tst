<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=4, WeakDG, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_WeakDG_SEM_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_WeakDG_SEM_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-11">4.24655e-06</value>
            <value variable="rhou" tolerance="1e-11">9.19532e-06</value>
            <value variable="rhov" tolerance="1e-11">8.62834e-06</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-10">3.11402e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-11">8.1427e-06</value>
            <value variable="rhou" tolerance="1e-10">2.15656e-05</value>
            <value variable="rhov" tolerance="1e-10">1.86306e-05</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-10">7.35728e-05</value>
        </metric>
    </metrics>
</test>


