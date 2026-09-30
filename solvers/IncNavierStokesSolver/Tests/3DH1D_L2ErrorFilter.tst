<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Computes L2norm for 3DH1D expansions</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>3DH1D_L2ErrorFilter.xml</parameters>
    <files>
        <file description="Session File">3DH1D_L2ErrorFilter.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">3.14159e+00</value>
            <value variable="v" tolerance="1e-11">3.14159e+00</value>
            <value variable="w" tolerance="1e-11">3.14159e+00</value>
            <value variable="theta" tolerance="1e-11">3.14159e+00</value>
            <value variable="p" tolerance="1e-09">0.000000e+00</value>
        </metric>
    </metrics>
</test>
