<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of NS equation with plunging, NACA0012 airfoil</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>NACA0012_Re400.xml PlungingAirfoil_3DH1D_parBCs.xml --npz 2</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">NACA0012_Re400.xml</file>
        <file description="Session File">PlungingAirfoil_3DH1D_parBCs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.197324</value>
            <value variable="v" tolerance="1e-6">0.182629</value>
            <value variable="w" tolerance="1e-6">0.0</value>
            <value variable="p" tolerance="1e-4">13.3313</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">1.14017</value>
            <value variable="v" tolerance="1e-5">1.28725</value>
            <value variable="w" tolerance="1e-6">0.0</value>
            <value variable="p" tolerance="1e-5">7.20143</value>
        </metric>
    </metrics>
</test>
