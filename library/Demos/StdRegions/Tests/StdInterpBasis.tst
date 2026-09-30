<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdInterpBasis Hexahedron Modified basis P=7 Q=8</description>
        <executable>StdInterpBasis</executable>
        <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 7 7 7 -p 8 8 8</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpBasis Prism Orthonormal basis P=6 Q=7</description>
        <executable>StdInterpBasis</executable>
        <parameters>-s prism -b Ortho_A Ortho_A Ortho_B -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.89694e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.77636e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpBasis Pyramid Modified basis P=7 Q=8</description>
        <executable>StdInterpBasis</executable>
        <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 7 7 7 -p 8 8 8</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpBasis Quadrilateral Lagrange basis P=6 Q=7</description>
        <executable>StdInterpBasis</executable>
        <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpBasis Tet Modified basis P=4 Q=5</description>
        <executable>StdInterpBasis</executable>
        <parameters> -s tetrahedron -b Modified_A Modified_B Modified_C -o 4 4 4 -p 5 5 5</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.89694e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.77636e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpBasis Triangle Modified basis P=7 Q=8</description>
        <executable>StdInterpBasis</executable>
        <parameters>-s triangle -b Modified_A Modified_B -o 7 7 -p 8 8</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpBasis Triangle Ortho basis P=7 Q=8</description>
        <executable>StdInterpBasis</executable>
        <parameters>-s triangle -b Ortho_A Ortho_B -o 7 7 -p 8 8</parameters>
        <metrics>
            <metric type="Linf" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="L2" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
</tests>
