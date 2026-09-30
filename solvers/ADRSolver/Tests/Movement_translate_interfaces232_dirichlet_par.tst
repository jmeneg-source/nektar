<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D advection with 3 zones and 2 non-conformal interfaces with periodic BCs on 4 processes with central zone translating</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_translate_interfaces232_dirichlet.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Movement_translate_interfaces232_dirichlet.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">1.55223e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">3.20845e-07</value>
        </metric>
    </metrics>
</test>


