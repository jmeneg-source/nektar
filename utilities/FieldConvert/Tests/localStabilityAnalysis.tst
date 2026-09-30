<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test localStabilityAnalysis module</description>
    <executable>FieldConvert</executable>
    <parameters>-f -e -m localStabilityAnalysis:opt=6:finVal1=0.05:numStep1=200 localStabilityAnalysis_mesh.xml localStabilityAnalysis_session.xml localStabilityAnalysis.fld</parameters>
    <files>
        <file description="Session File">localStabilityAnalysis_mesh.xml</file>
        <file description="Session File">localStabilityAnalysis_session.xml</file>
        <file description="Session File">localStabilityAnalysis.fld</file>
        <file description="Baseflow File">PSEMEAN.PRO</file>
    </files>

    <metrics>
        <metric type="Regex" id="1">
            <regex>
                \s*(-?\d\.\d*)\s*(-?\d\.\d*)\s*(-?\d\.\d*E[+-]?\d*)
            </regex>
            <matches>
                <match>
                    <!--Root 1-->
                    <field>0.069805</field>
                    <field>0.014125</field>
                    <field>0.20000000E+05</field>
                    <!--Root 2-->
                    <field>0.029448</field>
                    <field>-0.016068</field>
                    <field>0.20000000E+05</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
