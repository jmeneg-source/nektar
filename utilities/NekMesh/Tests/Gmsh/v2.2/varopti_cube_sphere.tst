<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Variational optimiser test on all tet cube/sphere</description>
    <executable>NekMesh</executable>
    <parameters>varopti_cube_sphere.msh varopti_cube_sphere-out.xml:xml:test -v -m varopti:hyperelastic:maxiter=5:nq=4</parameters>
    <files>
        <file description="Input File">varopti_cube_sphere.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>^.*Worst at end\s*: (-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?)</regex>
            <matches>
                <match>
                    <field id="0" tolerance="6e-2">5.652371e-01</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
