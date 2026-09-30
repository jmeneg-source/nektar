<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Laminar Channel Flow 3DH1D, P=4, 12 Fourier modes (MVM) with flowrate control and forcing direction explicitly defined</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-13">3.65972e-08</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="2e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-13">4.57503e-08</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="2e-12">0</value>
        </metric>
    </metrics>
</test>
