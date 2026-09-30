<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Testing 3D homogeneous 1D derivatives</description>
    <executable>Deriv3DHomo1D</executable>
    <parameters>Deriv3D_Homo1D.xml</parameters>
    <files>
        <file description="Session File">Deriv3D_Homo1D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="dudx" tolerance="1e-12">0</value>
            <value variable="dvdy" tolerance="1e-08">0.00575299</value>
            <value variable="dwdz" tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="dudx" tolerance="1e-12">0</value>
            <value variable="dvdy" tolerance="1e-08">0.00464427</value>
            <value variable="dwdz" tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>


