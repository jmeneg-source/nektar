<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Test Lagrangian points tracking in a laminar Channel Flow 3D</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>PointsTracking3D.xml </parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">PointsTracking3D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="Lx" tolerance="1e-6">0.346482</value>
            <value variable="Ly" tolerance="1e-6">0.590326</value>
            <value variable="Lz" tolerance="1e-6">0.5</value>
            <value variable="Lu" tolerance="1e-6">0</value>
            <value variable="Lv" tolerance="1e-6">0</value>
            <value variable="Lw" tolerance="1e-6">0.174796</value>
            <value variable="Lx" tolerance="1e-6">0.346482</value>
            <value variable="Ly" tolerance="1e-6">0.590326</value>
            <value variable="Lz" tolerance="1e-6">0.657319</value>
            <value variable="Lu" tolerance="1e-6">0</value>
            <value variable="Lv" tolerance="1e-6">0</value>
            <value variable="Lw" tolerance="1e-6">0.174796</value>
            <value variable="u" tolerance="1e-6">0</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">0</value>
        </metric>
    </metrics>
</test>


