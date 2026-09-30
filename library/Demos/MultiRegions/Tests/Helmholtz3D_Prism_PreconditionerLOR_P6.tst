<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 3D for Prisms, LOR Preconditioner, P=6</description>
    <executable>Helmholtz3D</executable>
    <parameters>Helmholtz3D_Prism_PreconditionerLOR_P6.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Prism_PreconditionerLOR_P6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">0.000000249129</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">0.00000420692</value>
        </metric>
    </metrics>
</test>
