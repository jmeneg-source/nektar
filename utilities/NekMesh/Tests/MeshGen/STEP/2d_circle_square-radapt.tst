<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D internal circle inside a square with r-adaption along the circumference</description>
    <executable>NekMesh</executable>
    <parameters> -m varopti:linearelastic:nq=2:radaptcurves=1-4:radaptscale=0.5:radaptrad=0.0015:subiter=20:maxiter=500:restol=1e-5:numthreads=2
    -m varopti:linearelastic:nq=2:radaptcurves=1-4:radaptscale=0.4:subiter=10:maxiter=500:restol=1e-5:numthreads=2
    -m varopti:hyperelastic:nq=5:numthreads=2 -m jac:list
    2d_circle_square-radapt.mcf 2d_circle_square-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">2d_circle_square-radapt.mcf</file>
        <file description="Input File 2">2d_circle_square.stp</file>
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
