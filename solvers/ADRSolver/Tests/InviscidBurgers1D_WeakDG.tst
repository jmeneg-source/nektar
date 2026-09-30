<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady WeakDG inviscidBurgers MODIFIED, P=10</description>
    <executable>ADRSolver</executable>
    <parameters>InviscidBurgers1D_WeakDG.xml</parameters>
    <files>
        <file description="Session File">InviscidBurgers1D_WeakDG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.36679</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.28739</value>
        </metric>
    </metrics>
</test>
