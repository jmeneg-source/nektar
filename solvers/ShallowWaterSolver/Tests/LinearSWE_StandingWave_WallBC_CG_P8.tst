<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Standing Wave, CG, P=8</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>LinearSWE_StandingWave_WallBC_CG_P8.xml</parameters>
    <files>
        <file description="Session File">LinearSWE_StandingWave_WallBC_CG_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-12">9.67385e-09</value>
            <value variable="u" tolerance="1e-12">5.10534e-07</value>
            <value variable="v" tolerance="1e-12">5.10534e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-12">1.50797e-09</value>
            <value variable="u" tolerance="1e-12">3.83906e-08</value>
            <value variable="v" tolerance="1e-12">3.83906e-08</value>
        </metric>
    </metrics>
</test>


