<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Channel flow file solution</description>
    <executable>ReviewSolution</executable>
    <parameters>ReviewSolution_2D.xml</parameters>
    <files>
        <file description="Session File">ReviewSolution_2D.xml</file>
        <file description="Session File">solution_0.chk</file>
        <file description="Session File">solution_1.chk</file>
        <file description="Session File">solution_2.chk</file>
        <file description="Session File">solution_3.chk</file>
        <file description="Session File">solution_4.chk</file>
        <file description="Session File">solution_5.chk</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">0</value>
            <value variable="v" tolerance="1e-08">0</value>
            <value variable="p" tolerance="1e-08">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">0</value>
            <value variable="v" tolerance="1e-08">0</value>
            <value variable="p" tolerance="1e-08">0</value>
        </metric>
    </metrics>
</test>
