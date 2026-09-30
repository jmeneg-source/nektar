<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> NekMesh ProjectCAD with Shaft HO - picky case for CADCurve association that can result in invalid J </description>
    <executable>NekMesh</executable>
    <parameters> -m projectcad:file=03_projectcad_NACA0012.stp:order=4 -m jac:list 03_projectcad_NACA0012.xml 03_projectcad_NACA0012-out.xml:xml:test -v </parameters>
    <files>
        <file description="Input File">03_projectcad_NACA0012.xml</file>
        <file description="Input File">03_projectcad_NACA0012.stp</file>
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
                    <field id="0">14</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="8">
            <regex>.*\s+Vertices- CADVertex \(also have CADSu\)\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0">14</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="9">
            <regex>.*\s+Vertices- CADCurve\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                    <field id="0" comparison="greater">349</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="10">
            <regex>.*\s+Edges CADCurve\s+N=\s+(\d+)</regex>
            <matches>
                <match>
                <field id="0" comparison="greater">356</field>
                </match>
            </matches>
        </metric>
        # HO Surf Module 
        <metric type="regex" id="10">
            <regex>.*Curved edges inserted\s+N\s*=\s*(\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics> 
</test>
