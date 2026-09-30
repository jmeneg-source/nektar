<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interp field to a box of points (also calculate cp and cp0)</description>
    <executable>FieldConvert</executable>
    <parameters>-f -e -m interpfield:fromxml=intake0.xml:fromfld=intake0.fld intake1.xml intake.vtu</parameters>
    <files>
        <file description="Session File">intake0.xml</file>
        <file description="Session File">intake1.xml</file>
        <file description="Fld File">intake0.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">26.7415</value>
            <value variable="y" tolerance="1e-4">22.3235</value>
            <value variable="rho" tolerance="1e-5">8.46494</value>
            <value variable="rhou" tolerance="1e-5">8.14898</value>
            <value variable="rhov" tolerance="1e-6">0.341581</value>
            <value variable="E" tolerance="1e-5">5.82780</value>
            <value variable="u" tolerance="1e-5">7.51314</value>
            <value variable="v" tolerance="1e-6">0.169977</value>
            <value variable="p" tolerance="1e-6">0.776144</value>
            <value variable="T" tolerance="1e-5">7.74504</value>
            <value variable="s" tolerance="5e-7">0.00449409</value>
            <value variable="a" tolerance="1e-5">2.54447</value>
            <value variable="Mach" tolerance="1e-4">22.4282</value>
            <value variable="Sensor" tolerance="5e-3">95.8186</value>
        </metric>

<metric type="Linf" id="2">
            <value variable="x" tolerance="1e-5">7.00000</value>
            <value variable="y" tolerance="1e-5">5.00000</value>
            <value variable="rho" tolerance="1e-5">6.09411</value>
            <value variable="rhou" tolerance="1e-5">4.13321</value>
            <value variable="rhov" tolerance="1e-6"> 0.508248</value>
            <value variable="E" tolerance="1e-5">3.78123</value>
            <value variable="u" tolerance="1e-5">1.18404</value>
            <value variable="v" tolerance="1e-6">0.152246</value>
            <value variable="p" tolerance="1e-6">0.96063</value>
            <value variable="T" tolerance="1e-5">2.15795</value>
            <value variable="s" tolerance="5e-7">0.00793072</value>
            <value variable="a" tolerance="1e-6">0.493238</value>
            <value variable="Mach" tolerance="1e-5">3.40483</value>
            <value variable="Sensor" tolerance="1e-5">0.982668</value>
        </metric>
    </metrics>
</test>

