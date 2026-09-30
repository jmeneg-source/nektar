<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NonlinSysIterDemo for a small nonlinear problem of dimension 4</description>
    <executable>NonlinSysIterDemo</executable>
    <parameters>NonlinSysIterDemo_LinearSys_NewtonScale.xml</parameters>
    <files>
        <file description="Session File">NonlinSysIterDemo_LinearSys_NewtonScale.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-06">1.000000e+00</value>
            <value variable="v" tolerance="1e-9">1.0e-9</value>
        </metric>
    </metrics>
</test>
