<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3DH1D cylinder flow, P=3, HomModes=8, cylinder defined via IB SPM function</description>    
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow3DH1D_SPM.xml</parameters>
    <files>
        <file description="Session File">CylFlow3DH1D_SPM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-04">63.3723</value>
            <value variable="v" tolerance="1e-05">2.00052</value>
            <value variable="w" tolerance="1e-11">1.05046e-06</value>
            <value variable="p" tolerance="1e-04">19.3327</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-04">1.7799</value>
            <value variable="v" tolerance="1e-05">1.07304</value>
            <value variable="w" tolerance="1e-12">3.13475e-07</value>
            <value variable="p" tolerance="1e-04">2.98637</value>
        </metric>
    </metrics>
</test>
