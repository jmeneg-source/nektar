<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Adjoint stability (Mod. Arnoldi): Channel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChannelStabilityAdjoint.xml</parameters>
    <files>
        <file description="Session File">ChannelStabilityAdjoint.xml</file>
	<file description="Session File">ChannelStabilityAdjoint.bse</file>
	<file description="Session File">ChannelStabilityAdjoint.rst</file>
    </files>
    <metrics>
        <!-- <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">2.58928</value>
            <value variable="v" tolerance="1e-6">0.00336105</value>
            <value variable="p" tolerance="1e-6">0.00596128</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">1</value>
            <value variable="v" tolerance="1e-6">0.00284298</value>
            <value variable="p" tolerance="1e-6">0.00593204</value>
        </metric> -->
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.00031,0.0349782</value>
            <value tolerance="0.001">1.00031,-0.0349782</value>
        </metric>
    </metrics>
</test>


