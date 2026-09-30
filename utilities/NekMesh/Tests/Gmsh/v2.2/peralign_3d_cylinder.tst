<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3d mesh with prism layers with periodic boundary using peralign</description>
    <executable>NekMesh</executable>
    <parameters>-m peralign:surf1=4:surf2=5:dir=z:orient -m jac:list peralign_3d_cylinder.msh cyl-out.xml:xml:test</parameters>
    <files>
        <file description="Input File"> peralign_3d_cylinder.msh </file>
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
        <metric type="nowarning" id="2"/>
    </metrics>
</test>
