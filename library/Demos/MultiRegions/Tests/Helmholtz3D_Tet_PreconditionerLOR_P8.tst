<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D for Tets, LOR Preconditioner, P=8</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Tet_PreconditionerLOR_P8.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Tet_PreconditionerLOR_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">0.000000217955</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">0.00000190495</value>
        </metric>
    </metrics>
</test>
