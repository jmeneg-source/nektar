<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of NS equation with rotation, simple domain without body</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>MovingRefFrame_Rot_naca0012.xml</parameters>
    <files>
        <file description="Session File">MovingRefFrame_Rot_naca0012.xml</file>
        <file description="Restart File">MovingRefFrame_Rot_naca0012.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">29.9834</value>
            <value variable="v" tolerance="5e-5">0.548735</value>
            <value variable="p" tolerance="5e-3">0.857797</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-4">1.33746</value>
            <value variable="v" tolerance="5e-4">0.947605</value>
            <value variable="p" tolerance="5e-3">1.98366</value>
        </metric>
    </metrics>
</test>
