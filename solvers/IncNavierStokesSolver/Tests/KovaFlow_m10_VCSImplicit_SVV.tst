<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay flow test case for using SVV with VCSImplicit</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m10_VCSImplicit.xml -I SpectralVanishingViscosity=DGKernel</parameters>
    <files>
        <file description="Session File">KovaFlow_m10_VCSImplicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">7.12999e-07</value>
            <value variable="v" tolerance="1e-08">4.33648e-07</value>
            <value variable="p" tolerance="1e-08">2.46092e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-07">1.55041e-06</value>
            <value variable="v" tolerance="1e-08">8.44357e-07</value>
            <value variable="p" tolerance="1e-07">4.57509e-06</value>
        </metric>
    </metrics>
</test>
