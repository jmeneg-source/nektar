<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdInterpDeriv Hex Lagrange basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s hexahedron -b GLL_Lagrange GLL_Lagrange GLL_Lagrange -o 7 7 7 -p 8 8 8 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Hex Mod basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s hexahedron -b Modified_A Modified_A Modified_A -o 7 7 7 -p 8 8 8 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Hex Ortho basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s hexahedron -b Ortho_A Ortho_A Ortho_A -o 7 7 7 -p 8 8 8 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Prism Mod basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s prism -b Modified_A Modified_A Modified_B  -o 7 7 7 -p 8 8 8 -P GaussLobattoLegendre GaussLobattoLegendre GaussRadauMAlpha1Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Prism Ortho basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s prism -b Ortho_A Ortho_A Ortho_B -o 7 7 7 -p 8 8 8 -P GaussLobattoLegendre GaussLobattoLegendre GaussRadauMAlpha1Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Pyr Mod basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 7 7 7 -p 8 8 8 -P GaussLobattoLegendre GaussLobattoLegendre GaussRadauMAlpha2Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Pyramid Ortho basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s pyramid -b Ortho_A Ortho_A OrthoPyr_C -o 7 7 7 -p 8 8 8 -P GaussLobattoLegendre GaussLobattoLegendre GaussRadauMAlpha2Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Quadrilateral Lagrange basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 7 7 -p 8 8 -P GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Quadrilateral Modified basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters>-s quadrilateral -b Modified_A Modified_A -o 7 7 -p 8 8 -P GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Quadrilateral Orth basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters>-s quadrilateral -b Ortho_A Ortho_A -o 7 7 -p 8 8 -P GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Seg Lagrange basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s segment -b GLL_Lagrange -o 7 -p 8 -P GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Seg Mod basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s segment -b Modified_A -o 7 -p 8 -P GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Seg Orth basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s segment -b Ortho_A -o 7 -p 8 -P GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Tet Mod basis P=8 Q=9</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s tetrahedron -b Modified_A Modified_B Modified_C -o 8 8 8 -p 9 9 9 -P GaussLobattoLegendre GaussRadauMAlpha1Beta0 GaussRadauMAlpha2Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Tet Orth basis P=8 Q=9</description>
        <executable>StdInterpDeriv</executable>
        <parameters> -s tetrahedron -b Ortho_A Ortho_B Ortho_C -o 8 8 8 -p 9 9 9 -P GaussLobattoLegendre GaussRadauMAlpha1Beta0 GaussRadauMAlpha2Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Tri Modified basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters>-s triangle -b Modified_A Modified_B -o 7 7 -p 8 8 -P GaussLobattoLegendre GaussRadauMAlpha1Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpDeriv Tri Ortho basis P=7 Q=8</description>
        <executable>StdInterpDeriv</executable>
        <parameters>-s triangle -b Ortho_A Ortho_B -o 7 7 -p 8 8 -P GaussLobattoLegendre GaussRadauMAlpha1Beta0</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-7">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-7">0</value>
            </metric>
        </metrics>
    </test>
</tests>
