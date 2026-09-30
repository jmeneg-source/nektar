<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Solitary Wave, DG, P=6</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>NonlinearSWE_Peregrine_SolitaryWave_DG_P6.xml</parameters>
    <files>
        <file description="Session File">NonlinearSWE_Peregrine_SolitaryWave_DG_P6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="h"  tolerance="1e-8">54.0898</value>
            <value variable="hu" tolerance="1e-8">4.62569</value>
            <value variable="hv" tolerance="1e-8">0.0052069</value>
            <value variable="z"  tolerance="1e-8">0.0823785</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h"  tolerance="1e-8">1.10008</value>
            <value variable="hu" tolerance="1e-8">0.344453</value>
            <value variable="hv" tolerance="1e-8">0.00270664</value>
            <value variable="z"  tolerance="1e-8">0.00837449</value>
        </metric>
    </metrics>
</test>


