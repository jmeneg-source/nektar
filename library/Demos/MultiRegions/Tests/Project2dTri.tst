<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Projection h-type convergence with 6 modes in 2D (triangles, continuous)</description>
    <executable>Project2D</executable>
    <parameters>Project_blank.xml -t -n 6</parameters>
    <files>
        <file description="Blank Session File">Project_blank.xml</file>
    </files>

    <metrics>
        <metric type="REGEX" id="1">
            <regex>
                Gradient: ([+-]?\d+(?:\.\d*)?)
            </regex>
            <matches>
                <match>
                    <field tolerance="1e-3" id="1">-6.0027621894</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>


