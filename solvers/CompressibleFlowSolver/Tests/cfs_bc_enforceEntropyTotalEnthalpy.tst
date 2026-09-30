<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test different boundary condition enforcements: entropy-total enthalpy </description>
    <executable>CompressibleFlowSolver </executable>
    <parameters>square_mix.xml session_enforceEntropyTotalEnthalpy.xml </parameters>
    <files>
        <file description="Session File">square_mix.xml</file>
        <file description="Session File">session_enforceEntropyTotalEnthalpy.xml</file>
    </files>
    <metrics>
        <metric type="Linf" id="1">
            <value variable="rho"  tolerance="1e-12">1.0</value>
            <value variable="rhou" tolerance="1e-12">1.0</value>
            <value variable="rhov" tolerance="4e-12">0.0</value>
            <value variable="E"    tolerance="1e-12">3.29018</value>
        </metric>
        <metric type="L2" id="2">
            <value variable="rho"  tolerance="1e-12">2.0</value>
            <value variable="rhou" tolerance="1e-12">2.0</value>
            <value variable="rhov" tolerance="3e-12">0.0</value>
            <value variable="E"    tolerance="1e-12">6.58036</value>
        </metric>
 
    </metrics>
</test>
