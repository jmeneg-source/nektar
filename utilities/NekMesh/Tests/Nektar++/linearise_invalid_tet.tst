<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Remove curvature from invalid tetrahedrons</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list -m linearise:invalid -m jac:list linearise_invalid_tet.xml linearise_invalid_tet-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">linearise_invalid_tet.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">3</field>
                </match>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
