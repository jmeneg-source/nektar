<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Generate vortex-induced velocity </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m vortexinducedvelocity:vortex=vortexfilament.dat  Helmholtz.xml dummy.fld</parameters>
    <files>
        <file description="Session File">Helmholtz.xml</file>
        <file description="vortex parameters file">vortexfilament.dat</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">1.1547</value>
            <value variable="y" tolerance="1e-4">1.1547</value>
            <value variable="u" tolerance="2e-3">0.443226</value>
            <value variable="v" tolerance="1e-2">0.52802</value>
            <value variable="w" tolerance="1e-6">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">1.0</value>
            <value variable="y" tolerance="1e-4">1.0</value>
            <value variable="u" tolerance="4.1e-2">0.693835</value>
            <value variable="v" tolerance="1.1e-2">0.740137</value>
            <value variable="w" tolerance="1e-6">0</value>
        </metric>
    </metrics>
</test>
