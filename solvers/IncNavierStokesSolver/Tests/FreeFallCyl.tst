<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of a free falling circular cylinder</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>FreeFallCyl.xml FreeFallCylc.xml --set-start-time 0 --set-start-chknumber 0</parameters>
    <files>
        <file description="Mesh File">FreeFallCyl.xml</file>
        <file description="Session File">FreeFallCylc.xml</file>
        <file description="initial condition File">FreeFallCyl.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">2.9125363e-01</value>
            <value variable="v" tolerance="5e-4">3.0579505e-01</value>
            <value variable="p" tolerance="5e-2">1.9414285e-02</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-3">1.5042508e-01</value>
            <value variable="v" tolerance="5e-4">7.9868420e-02</value>
            <value variable="p" tolerance="5e-2">7.6541462e-03</value>
        </metric>
    </metrics>
</test>
