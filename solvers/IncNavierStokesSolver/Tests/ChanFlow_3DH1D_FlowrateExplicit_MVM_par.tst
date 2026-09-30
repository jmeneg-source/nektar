<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Laminar Channel Flow 3DH1D, P=4, 12 Fourier modes (MVM) with flowrate control and forcing direction explicitly defined</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ChanFlow_3DH1D_FlowrateExplicit_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-13">3.66285e-08</value>
            <value variable="v" tolerance="1e-12">9.60875e-11</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-10">2.76791e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="3e-13">4.67694e-08</value>
            <value variable="v" tolerance="1e-12">2.21341e-10</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-10">3.62202e-07</value>
        </metric>
    </metrics>
</test>
