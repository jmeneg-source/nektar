<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Scaling a 3D mesh</description>
    <executable>NekMesh</executable>
    <parameters>-m scale:scaleX=2.0:scaleY=2.0:scaleZ=2.0 -v scaleHO.xml
        scaleHO_scaled.xml:xml:test:stats</parameters>
    <files>
        <file description="Input File">scaleHO.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>^.*Lower mesh extent\s*:\s*(-?\d+\.?\d*)\s*(-?\d+\.?\d*)\s*
                (-?\d+\.?\d*)</regex>
            <matches>
                <match>
                    <field id="0">-1</field>
                    <field id="1">-1</field>
                    <field id="2">0</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="1">
            <regex>^.*Upper mesh extent\s*:\s*(-?\d+\.?\d*)\s*(-?\d+\.?\d*)\s*
                (-?\d+\.?\d*)</regex>
            <matches>
                <match>
                    <field id="0">3</field>
                    <field id="1">1</field>
                    <field id="2">0.6</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
