<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Solitary Wave, DG, P=6</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>NonlinearSWE_Peregrine_SolitaryWave_DG_P6_implicit.xml</parameters>
    <files>
        <file description="Session File">NonlinearSWE_Peregrine_SolitaryWave_DG_P6_implicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="h"  tolerance="1e-8">54.0874</value>
            <value variable="hu" tolerance="1e-8">4.62012</value>
            <value variable="hv" tolerance="1e-8">0.000315187</value>
            <value variable="z"  tolerance="1e-8">0.334792</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="h"  tolerance="1e-8">1.1002</value>
            <value variable="hu" tolerance="1e-8">0.341875</value>
            <value variable="hv" tolerance="1e-8">0.000356266</value>
            <value variable="z"  tolerance="1e-8">0.0324715</value>
        </metric>
    </metrics>
</test>


