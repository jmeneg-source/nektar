<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_LDG_MODIFIED_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_LDG_MODIFIED_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-05">5.50775</value>
            <value variable="rhou" tolerance="1e-02">2016.35</value>
            <value variable="rhov" tolerance="1e-04">24.5452</value>
            <value variable="rhow" tolerance="1e-04">24.4789</value>
            <value variable="E" tolerance="1e+01">6.27023e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-05">0.297685</value>
            <value variable="rhou" tolerance="1e-04">86.9862</value>
            <value variable="rhov" tolerance="1e-04">15.3531</value>
            <value variable="rhow" tolerance="1e-05">1.00004</value>
            <value variable="E" tolerance="1e+00">277498</value>
        </metric>
    </metrics>
</test>


