<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG linear expansion P=7</description>
    <executable>Helmholtz2D_LinMesh</executable>
    <parameters>Helmholtz2D_P7_LinSol.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_P7_LinSol.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.78788e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.33227e-15</value>
        </metric>
    </metrics>
</test>


