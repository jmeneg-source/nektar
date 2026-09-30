<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>HDG Helmholtz3D Homogeneous 1D</description>
    <executable>HDGHelmholtz3DHomo1D</executable>
    <parameters>Helmholtz3D_Homo1D.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Homo1D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">0.00550986</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">0.02504420</value>
        </metric>
    </metrics>
</test>


