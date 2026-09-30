<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using LORPreconditioner only for all Vars</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_PreconditionerLOR_all.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_PreconditionerLOR_all.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.000620708</value>
            <value variable="v" tolerance="1e-6">0.00013943</value>
            <value variable="w" tolerance="1e-6">0.000127092</value>
            <value variable="p" tolerance="1e-6">0.000415308</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00210742</value>
            <value variable="v" tolerance="1e-6">0.000611576</value>
            <value variable="w" tolerance="1e-6">0.000443681</value>
            <value variable="p" tolerance="1e-6">0.0024698</value>
        </metric>
    </metrics>
</test>
