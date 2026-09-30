<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Formulation of a pressure-decomposed velocity correction scheme for a 2-DoFs self-propelled NACA0012 airfoil</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>SelfPropelledAirfoil_2DoFs_PD-VCS.xml</parameters>
    <files>
        <file description="Session File">SelfPropelledAirfoil_2DoFs_PD-VCS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">6.1462253e-02</value>
            <value variable="v" tolerance="5e-4">6.2272365e-02</value>
            <value variable="p" tolerance="5e-2">5.2780029e+01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-4">1.1165357e+00</value>
            <value variable="v" tolerance="5e-4">1.3776587e+00</value>
            <value variable="p" tolerance="5e-2">1.8841437e+01</value>
        </metric>
    </metrics>
</test>
