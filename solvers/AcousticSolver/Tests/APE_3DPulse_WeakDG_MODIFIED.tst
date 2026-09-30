<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>AcousticSolver</executable>
    <parameters>APE_3DPulse_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_3DPulse_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">7.33275</value>
            <value variable="u" tolerance="1e-4">0.0101853</value>
            <value variable="v" tolerance="1e-4">0.0101853</value>
            <value variable="w" tolerance="1e-4">0.0101853</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-3">21.8578</value>
            <value variable="u" tolerance="1e-5">0.037688</value>
            <value variable="v" tolerance="1e-5">0.037688</value>
            <value variable="w" tolerance="1e-5">0.037688</value>
        </metric>
    </metrics>
</test>

