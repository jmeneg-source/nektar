<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Hexahedral elements, P=3, Successive RHS(5), par(2)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch Hex_channel_m3_srhs.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Hex_channel_m3_srhs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.00209513</value>
            <value variable="v" tolerance="1e-8">0.00137204</value>
            <value variable="w" tolerance="1e-8">0.01153110</value>
            <value variable="p" tolerance="1e-8">0.38012600</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">0.00495176</value>
            <value variable="v" tolerance="1e-8">0.00372357</value>
            <value variable="w" tolerance="1e-8">0.07484400</value>
            <value variable="p" tolerance="1e-8">0.90268500</value>
        </metric>
    </metrics>
</test>
