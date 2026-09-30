<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Standing Wave, DG, P=4</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>LinearSWE_StandingWave_WallBC_DG_P4.xml</parameters>
    <files>
        <file description="Session File">LinearSWE_StandingWave_WallBC_DG_P4.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-12">4.44046e-05</value>
            <value variable="u" tolerance="1e-12">0.000104914</value>
            <value variable="v" tolerance="1e-12">0.000104914</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-12">8.19871e-06</value>
            <value variable="u" tolerance="1e-12">1.53915e-05</value>
            <value variable="v" tolerance="1e-12">1.53915e-05</value>
        </metric>
    </metrics>
</test>


