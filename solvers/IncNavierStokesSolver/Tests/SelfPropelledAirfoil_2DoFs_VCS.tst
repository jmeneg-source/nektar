<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Formulation of a velocity correction scheme for a 2-DoFs self-propelled NACA0012 airfoil</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>SelfPropelledAirfoil_2DoFs_VCS.xml</parameters>
    <files>
        <file description="Session File">SelfPropelledAirfoil_2DoFs_VCS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">0.462121</value>
            <value variable="v" tolerance="5e-4">0.459917</value>
            <value variable="p" tolerance="5e+0">43325.3</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-4">5.20019</value>
            <value variable="v" tolerance="5e-4">5.54633</value>
            <value variable="p" tolerance="5e+0">15037</value>
        </metric>
    </metrics>
</test>
