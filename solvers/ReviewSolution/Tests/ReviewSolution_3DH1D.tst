<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D, equation solution</description>
    <executable>ReviewSolution</executable>
    <parameters>ReviewSolution_3DH1D.xml</parameters>
    <files>
        <file description="Session File">ReviewSolution_3DH1D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">2.14538e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">3.90343e-08</value>
        </metric>
    </metrics>
</test>
