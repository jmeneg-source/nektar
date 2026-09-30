<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fully 3D Compressible Channel Flow. Synthetic turbulence generation is introduced in the flow field using forcing tag. This test case uses the explicit solver. </description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ChanFlow3D_infTurbExpl.xml</parameters>
    <files>
        <file description="Session File">ChanFlow3D_infTurbExpl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="1e-7">0.0496257</value>
            <value variable="rhou" tolerance="1e-7">0.0510006</value>
            <value variable="rhov" tolerance="1e-8">0.000046682</value>
            <value variable="rhow" tolerance="1e-8">0.0000572015</value>
            <value variable="E"    tolerance="1e-5">8.88776</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="1e-7">0.0140208</value>
            <value variable="rhou" tolerance="1e-7">0.0172187</value>
            <value variable="rhov" tolerance="1e-8">0.000614856</value>
            <value variable="rhow" tolerance="1e-8">0.000399892</value>
            <value variable="E"    tolerance="1e-5">2.51459</value>
        </metric>
    </metrics>
</test>


