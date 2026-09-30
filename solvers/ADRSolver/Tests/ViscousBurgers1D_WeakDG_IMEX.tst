<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady WeakDG viscous Burgers MODIFIED IMEX, P=10</description>
    <executable>ADRSolver</executable>
    <parameters>ViscousBurgers1D_WeakDG_IMEX.xml</parameters>
    <files>
        <file description="Session File">ViscousBurgers1D_WeakDG_IMEX.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.33284</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.0</value>
        </metric>
    </metrics>
</test>
