<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Hex Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_Kovasnay_SubStep.xml</parameters>
    <files>
        <file description="Session File">Hex_Kovasnay_SubStep.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.0525213 </value>
            <value variable="v" tolerance="1e-6">0.00347786</value>
            <value variable="w" tolerance="1e-6">0.00339647</value>
            <value variable="p" tolerance="1e-6">0.0136371 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0637019 </value>
            <value variable="v" tolerance="1e-6">0.00775727</value>
            <value variable="w" tolerance="1e-6">0.00629852</value>
            <value variable="p" tolerance="1e-6">0.0471015 </value>
        </metric>
    </metrics>
</test>
