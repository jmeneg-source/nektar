<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-09">0.000683625</value>
            <value variable="rhou" tolerance="1e-04">68.0612</value>
            <value variable="rhov" tolerance="1e-06">0.250790</value>
             <value variable="rhow" tolerance="1e-11">9.74017e-06</value>
            <value variable="E" tolerance="1e-01">24775.4</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-09">0.00207327</value>
            <value variable="rhou" tolerance="1e-04">83.3215</value>
            <value variable="rhov" tolerance="1e-06">0.752264</value>
             <value variable="rhow" tolerance="1e-10">3.39885e-05</value>
            <value variable="E" tolerance="1e-01">18732.3</value>
        </metric>
    </metrics>
</test>


