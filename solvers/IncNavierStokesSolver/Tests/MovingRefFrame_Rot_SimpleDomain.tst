<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of NS equation with rotation, simple domain without body</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>MovingRefFrame_Rot_SimpleDomain.xml</parameters>
    <files>
        <file description="Session File">MovingRefFrame_Rot_SimpleDomain.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-7">0.000433915</value>
            <value variable="v" tolerance="5e-7">0.000332251</value>
            <value variable="p" tolerance="5e-5">0.0329402</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-7">0.000849685</value>
            <value variable="v" tolerance="5e-7">0.000433858</value>
            <value variable="p" tolerance="5e-5">0.0100733</value>
        </metric>
    </metrics>
</test>
