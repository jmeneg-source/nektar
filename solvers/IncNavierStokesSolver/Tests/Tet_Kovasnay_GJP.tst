<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using GJP Stabilisation and dealiasing</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_GJP.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.000873329</value>
            <value variable="v" tolerance="1e-8">0.000449329</value>
            <value variable="w" tolerance="1e-10">0.000159493</value>
            <value variable="p" tolerance="1e-8">0.00205201</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">0.00370413</value>
            <value variable="v" tolerance="1e-8">0.00336431</value>
            <value variable="w" tolerance="1e-8">0.000684325</value>
            <value variable="p" tolerance="1e-8">0.020534</value>
        </metric>
    </metrics>
</test>
