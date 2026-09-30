<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Create output PNG via Python VTK interface</description>
    <executable python="true">python_vtk.py</executable>
    <parameters>Tet_channel_m3.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_m3.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^(?:xmin|xmax|ymin|ymax|zmin|zmax|cx|cy|cz): ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf)</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-12">0.0</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">1.0000000000000002</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">0.0</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">1.0000000000000002</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">0.0</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">1.0000000000000002</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">0.5000000000000001</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">0.5000000000000001</field>
                </match>
                <match>
                    <field id="0" tolerance="1e-12">0.5000000000000001</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
