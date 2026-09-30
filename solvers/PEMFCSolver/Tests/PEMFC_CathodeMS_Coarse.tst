<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Steady cathode two-species PEMFC system (PEMFCCathodeMS), coarse 24x48 parallel flow field, P=3, Picard iteration with Butler-Volmer kinetics on a frozen penalized flow</description>
    <executable>PEMFCSolver</executable>
    <parameters>PEMFC_CathodeMS_Coarse_mesh.xml PEMFC_CathodeMS_Coarse.xml</parameters>
    <files>
        <file description="Mesh File">PEMFC_CathodeMS_Coarse_mesh.xml</file>
        <file description="Session File">PEMFC_CathodeMS_Coarse.xml</file>
        <file description="Advection velocity field">velocity.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="wO2" tolerance="1e-6">0.00282351</value>
            <value variable="wH2O" tolerance="1e-6">0.00230188</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="wO2" tolerance="1e-6">0.228</value>
            <value variable="wH2O" tolerance="1e-6">0.229603</value>
        </metric>
    </metrics>
</test>
