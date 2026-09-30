<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 2D for Quads, LOR Preconditioner, P=6</description>
    <executable>Helmholtz2D</executable>
    <parameters>Helmholtz2D_Quad16_PreconditionerLOR_P6.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_Quad16_PreconditionerLOR_P6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-4">0.000000249129</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-4">0.00000420692</value>
        </metric>
    </metrics>
</test>
