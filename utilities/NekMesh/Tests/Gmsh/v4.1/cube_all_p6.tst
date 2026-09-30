<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.1) linear mesh of cube, all elements, order 6</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list cube_all_p6.msh cube_all-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">cube_all_p6.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
