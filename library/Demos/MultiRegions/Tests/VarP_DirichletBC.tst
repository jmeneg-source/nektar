<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D for Tets with Variable P</description>
    <executable>Helmholtz3D</executable>
    <parameters>VarP_DirichletBC.xml</parameters>
    <files>
        <file description="Session File">VarP_DirichletBC.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.0</value>
        </metric>
    </metrics>
</test>


