<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using (Semi-Implicit) GJP Stabilisation and dealiasing with the Implicit VCS</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_GJP.xml -I SolverType=VCSImplicit</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-4">0.000871731</value>
            <value variable="v" tolerance="1e-5">0.000453683</value>
            <value variable="w" tolerance="1e-5">0.000159363</value>
            <value variable="p" tolerance="1e-4">0.00205074</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-4">0.0036831</value>
            <value variable="v" tolerance="1e-4">0.00342033</value>
            <value variable="w" tolerance="1e-4">0.000675693</value>
            <value variable="p" tolerance="1e-4">0.0204599</value>
        </metric>
    </metrics>
</test>
