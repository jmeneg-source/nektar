<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D flexible cylinder flow simulation using "MovingBody" module</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_MovBody.xml</parameters>
    <files>
        <file description="Session File">CylFlow_MovBody.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-04">78.4858</value>
            <value variable="v" tolerance="1e-05">5.28416</value>
            <value variable="w" tolerance="1e-08">0.00921949</value>
            <value variable="p" tolerance="1e-03">237.459</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-05">1.44442</value>
            <value variable="v" tolerance="1e-06">0.670923</value>
            <value variable="w" tolerance="1e-08">0.004443</value>
            <value variable="p" tolerance="1e-05">5.69645</value>
        </metric>
    </metrics>
</test>
