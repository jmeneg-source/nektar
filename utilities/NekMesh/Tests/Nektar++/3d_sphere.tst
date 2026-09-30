<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test quality of half sphere</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list:quality 3d_sphere.xml 3d_sphere-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">3d_sphere.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>.*Integration of Jacobian: ([0-9]*\.[0-9]+)%?</regex>
            <matches>
                <match>
                    <field id="0" tolerance="2e-2">99.89</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
