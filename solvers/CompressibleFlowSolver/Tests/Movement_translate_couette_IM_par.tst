<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NS, Couette flow with non-conformal translate, exact solution, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Movement_translate_couette_IM.xml</parameters>
    <processes>6</processes>
    <files>
        <file description="Session File">Movement_translate_couette_IM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">2.92134e-07</value>
            <value variable="rhou" tolerance="1e-6">9.12755e-08</value>
            <value variable="rhov" tolerance="1e-6">4.45173e-07</value>
            <value variable="E" tolerance="1e-6">5.24223e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">5.96995e-07</value>
            <value variable="rhou" tolerance="1e-6">2.91595e-07</value>
            <value variable="rhov" tolerance="1e-6">8.17435e-07</value>
            <value variable="E" tolerance="2e-6">9.08052e-07</value>
        </metric>
    </metrics>
</test>
