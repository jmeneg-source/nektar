<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Formulation of a  pressure-decomposed velocity correction scheme for a 1-DoF self-propelled NACA0012 airfoil</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>SelfPropelledAirfoil_1DoF_PD-VCS.xml</parameters>
    <files>
        <file description="Session File">SelfPropelledAirfoil_1DoF_PD-VCS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">1.0496052e-01</value>
            <value variable="v" tolerance="5e-4">1.0544453e-01</value>
            <value variable="p" tolerance="5e-2">7.5141928e+00</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-4">1.9880973e+00</value>
            <value variable="v" tolerance="5e-4">2.5251975e+00</value>
            <value variable="p" tolerance="5e-2">7.8161777e+00</value>
        </metric>
    </metrics>
</test>
