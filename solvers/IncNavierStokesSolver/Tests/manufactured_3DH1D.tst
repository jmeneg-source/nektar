<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Manufactured solution varying in space and time; solved via 2.5D simulation</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>manufactured_3DH1D.xml</parameters>
    <files>
        <file description="Session File">manufactured_3DH1D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">7.39590e-06</value>
            <value variable="v" tolerance="1e-11">7.39590e-06</value>
            <value variable="w" tolerance="1e-11">7.89264e-06</value>
            <value variable="p" tolerance="1e-09">0.000124706</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">5.70178e-06</value>
            <value variable="v" tolerance="1e-11">5.70178e-06</value>
            <value variable="w" tolerance="1e-10">1.00807e-05</value>
            <value variable="p" tolerance="1e-10">9.87479e-05</value>
        </metric>
    </metrics>
</test>
