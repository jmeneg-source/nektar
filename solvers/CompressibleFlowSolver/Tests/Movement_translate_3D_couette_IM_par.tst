<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NS, 3D Couette flow with non-conformal translate, exact solution, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Movement_translate_3D_couette_IM.xml</parameters>
    <processes>6</processes>
    <files>
        <file description="Session File">Movement_translate_3D_couette_IM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">1.41721e-06</value>
            <value variable="rhou" tolerance="1e-6">1.74764e-07</value>
            <value variable="rhov" tolerance="1e-6">2.04018e-06</value>
            <value variable="rhow" tolerance="1e-6">3.85054e-13</value>
            <value variable="E" tolerance="1e-6">2.52549e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">3.37801e-06</value>
            <value variable="rhou" tolerance="1e-6">6.63742e-07</value>
            <value variable="rhov" tolerance="1e-6">5.16763e-06</value>
            <value variable="rhow" tolerance="1e-6">2.23365e-12</value>
            <value variable="E" tolerance="2e-6">5.68531e-06</value>
        </metric>
    </metrics>
</test>
