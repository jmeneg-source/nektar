<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=9 with large dt and taylor-hood expansion</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m10_taylorHood_VCSImplicitLoc.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m10_taylorHood_VCSImplicitLoc.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">1.54223e-08</value>
            <value variable="v" tolerance="1e-08">1.33973e-08</value>
            <value variable="p" tolerance="1e-08">2.01951e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">4.73025e-08</value>
            <value variable="v" tolerance="1e-08">3.11874e-08</value>
            <value variable="p" tolerance="1e-07">2.41298e-07</value>
        </metric>
    </metrics>
</test>
