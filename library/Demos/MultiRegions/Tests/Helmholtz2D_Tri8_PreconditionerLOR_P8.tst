<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 2D for Triangles, LOR Preconditioner, P=8</description>
    <executable>Helmholtz2D</executable>
    <parameters>Helmholtz2D_Tri8_PreconditionerLOR_P8.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_Tri8_PreconditionerLOR_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-6">0.000562378</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.000857352</value>
        </metric>
    </metrics>
</test>
