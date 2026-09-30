<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fully 3D Compressible Channel Flow. Synthetic turbulence generation is introduced in the flow field using forcing tag. This test cases uses the implicit solver. </description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ChanFlow3D_infTurbImpl.xml</parameters>
    <files>
        <file description="Session File">ChanFlow3D_infTurbImpl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="1e-7">0.0496257</value>
            <value variable="rhou" tolerance="1e-7">0.0510022</value>
            <value variable="rhov" tolerance="1e-9">0.0000771330</value>
            <value variable="rhow" tolerance="1e-8">0.000135185</value>
            <value variable="E"    tolerance="1e-5">8.88777</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="1e-7">0.0140176</value>
            <value variable="rhou" tolerance="1e-7">0.0173112</value>
            <value variable="rhov" tolerance="1e-8">0.000770763</value>
            <value variable="rhow" tolerance="1e-8">0.000988695</value>
            <value variable="E"    tolerance="1e-5">2.51265</value>
        </metric>
    </metrics>
</test>


