<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> NekMesh ProjectCAD with Shaft HO - checks if it can reconstruct HO-infromation </description>
    <executable>NekMesh</executable>
    <parameters> -m projectcad:file=00_projectcad_shaft.stp:order=3:surfopti=1:cLength=1:ho -m jac:list 01_projectcad_shaftHO.xml 01_projectcad_shaftHO-out.xml:xml:test -v </parameters>
    <files>
        <file description="Input File">01_projectcad_shaftHO.xml</file>
        <file description="Input File">00_projectcad_shaft.stp</file>
    </files>
    <metrics>
        # Mesh Quality
        <metric type="regex" id="1">
            <regex>.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>

        # Main CAD Reconstruction
        <metric type="regex" id="2">
            <regex>.*\s+Vertices No CAD \(Includes PlanarSurf\)\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>

        <metric type="regex" id="3">
            <regex>.*\s+Edges without CADObject\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>

        <metric type="regex" id="4">
            <regex>.*\s+Faces without CADObject\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>

        # Detailed CAD Reconstruction
        <metric type="regex" id="7">
            <regex>.*\s+Vertices CAD3 or more Surf\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0">8</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="8">
            <regex>.*\s+Vertices- CADVertex \(also have CADSu\)\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0">8</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="9">
            <regex>.*\s+Vertices- CADCurve\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0">44</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="10">
            <regex>.*\s+Edges CADCurve\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                <field id="0" comparison="greater">48</field>
                </match>
            </matches>
        </metric>
        # HO Surf Module 
        <metric type="regex" id="10">
            <regex>.*Curved edges inserted\s+N\s*=\s*(\d+)</regex>
            <matches>
                <match>
                    <field id="0">342</field>
                </match>
            </matches>
        </metric>
    </metrics> 
</test>
