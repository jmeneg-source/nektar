<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Channel Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_Tri_Lagrange.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_Tri_Lagrange.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.44199e-12</value>
            <value variable="v" tolerance="1e-12">4.18302e-13</value>
            <value variable="p" tolerance="1e-8">3.31598e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.24133e-11</value>
            <value variable="v" tolerance="1e-12">6.22878e-12</value>
            <value variable="p" tolerance="1e-8">5.98182e-10</value>
        </metric>
    </metrics>
</test>
