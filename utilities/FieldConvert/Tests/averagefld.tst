<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3DH1D powerspectrum </description>
    <executable>FieldConvert</executable>
    <parameters> -e -f -m averagefld:inputfld=chan3DH1D_%d.chk:range=0,1,4 chan3DH1D.xml phaseavg.fld </parameters>
    <files>
        <file description="Session File">chan3DH1D.xml</file>
        <file description="Session File">chan3DH1D_0.chk</file>
        <file description="Session File">chan3DH1D_1.chk</file>
        <file description="Session File">chan3DH1D_2.chk</file>
        <file description="Session File">chan3DH1D_3.chk</file>
        <file description="Session File">chan3DH1D_4.chk</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.577350</value>
            <value variable="y" tolerance="1e-6">0.577350</value>
            <value variable="z" tolerance="1e-6">0.522913</value>
            <value variable="u" tolerance="1e-6">0.182574</value>
            <value variable="v" tolerance="1e-5">0.163299</value>
            <value variable="w" tolerance="1e-6">0.143376</value>
            <value variable="p" tolerance="1e-6">0.923760</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-5">1</value>
            <value variable="y" tolerance="1e-5">1</value>
            <value variable="z" tolerance="1e-6">0.875</value>
            <value variable="u" tolerance="1e-6">0.25</value>
            <value variable="v" tolerance="1e-6">0.341421</value>
            <value variable="w" tolerance="1e-6">0.35345</value>
            <value variable="p" tolerance="1e-5">1.6</value>
        </metric>
    </metrics>
</test>
