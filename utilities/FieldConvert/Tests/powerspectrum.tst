<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3DH1D powerspectrum </description>
    <executable>FieldConvert</executable>
    <parameters> -e -f -m powerspectrum:vars=w:box=0,0.5,0,0.5 chan3DH1D.xml chan3DH1D_0.chk chan3DH1D.plt</parameters>
    <files>
        <file description="Session File">chan3DH1D.xml</file>
        <file description="Session File">chan3DH1D_0.chk</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="E" tolerance="1e-8">0.00221566</value>
            <value variable="x" tolerance="1e-6">0.577350</value>
            <value variable="y" tolerance="1e-6">0.577350</value>
            <value variable="z" tolerance="1e-6">0.522913</value>
            <value variable="u" tolerance="1e-6">0.183207</value>
            <value variable="v" tolerance="1e-6">0.816497</value>
            <value variable="w" tolerance="1e-6">0.842564</value>
            <value variable="p" tolerance="1e-6">0</value>
            <value variable="regions" tolerance="1e-6">0.5</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-5">1</value>
            <value variable="y" tolerance="1e-5">1</value>
            <value variable="z" tolerance="1e-6">0.875</value>
            <value variable="u" tolerance="1e-6">0.25</value>
            <value variable="v" tolerance="1e-5">1.70711</value>
            <value variable="w" tolerance="1e-5">1.62067</value>
            <value variable="p" tolerance="1e-5">0</value>
            <value variable="regions" tolerance="1e-5">1</value>
        </metric>
    </metrics>
</test>
