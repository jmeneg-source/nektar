<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fully 3D Channel Flow. Synthetic turbulence generation is introduced in the flow field using forcing tag. </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow3D_infTurb.xml</parameters>
    <files>
        <file description="Session File">ChanFlow3D_infTurb.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">3.64515</value>
            <value variable="v" tolerance="1e-7">0.0182432</value>
            <value variable="w" tolerance="1e-7">0.033767</value>
            <value variable="p" tolerance="1e-7">0.0201723</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">1.32574</value>
            <value variable="v" tolerance="1e-7">0.0871106</value>
            <value variable="w" tolerance="1e-5">0.179135</value>
            <value variable="p" tolerance="1e-7">0.0363918</value>
        </metric>
    </metrics>
</test>


