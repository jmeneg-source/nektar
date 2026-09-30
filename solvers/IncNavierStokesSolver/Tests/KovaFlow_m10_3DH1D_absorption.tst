<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=9</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m10_3DH1D_absorption.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m10_3DH1D_absorption.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-07">0.0102828</value>
            <value variable="v" tolerance="1e-07">0.0106737</value>
            <value variable="w" tolerance="1e-08">0.00110054</value>
            <value variable="p" tolerance="1e-08">0.0895775</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-07">0.0231446</value>
            <value variable="v" tolerance="1e-07">0.0447882</value>
            <value variable="w" tolerance="1e-08">0.00442845</value>
            <value variable="p" tolerance="1e-06">0.224952</value>
        </metric>
    </metrics>
</test>
