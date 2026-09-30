<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_SubStep.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_SubStep.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.00627866</value>
            <value variable="v" tolerance="1e-6">0.00199549</value>
            <value variable="w" tolerance="1e-6">0.00192544</value>
            <value variable="p" tolerance="1e-6">0.00639417</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0223768</value>
            <value variable="v" tolerance="1e-6">0.00706921</value>
            <value variable="w" tolerance="1e-6">0.00628175</value>
            <value variable="p" tolerance="2e-6">0.038014</value>
        </metric>
    </metrics>
</test>
