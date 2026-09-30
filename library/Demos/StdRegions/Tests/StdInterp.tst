<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdInterp Hex Modified basis P=6 Q=7</description>
        <executable>StdInterp</executable>
        <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 6 6 6 -p 7 7 7 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
        <metrics>
            <metric type="Linf" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="L2" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterp Prism Modified basis P=6 Q=7</description>
        <executable>StdInterp</executable>
        <parameters>-s prism -b Modified_A Modified_A Modified_B -o 6 6 6 -p 7 7 7 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
        <metrics>
            <metric type="Linf" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="L2" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterpBasis Prism Orthonormal basis P=6 Q=7</description>
        <executable>StdInterp</executable>
        <parameters>-s prism -b Ortho_A Ortho_A Ortho_B -o 6 6 6 -p 7 7 7 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
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
        <description>StdInterp Pyramid Modified basis P=6 Q=7</description>
        <executable>StdInterp</executable>
        <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 7 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
        <metrics>
            <metric type="Linf" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="L2" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterp Quadrilateral Lagrange basis P=6 Q=7</description>
        <executable>StdInterp</executable>
        <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7 -P GaussGaussLegendre GaussGaussLegendre</parameters>
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
        <description>StdInterp Seg Modified_A basis P=7 Q=8</description>
        <executable>StdInterp</executable>
        <parameters>-s segment -b Modified_A -o 7 -p 8 -P GaussGaussLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.81216e-16</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.33227e-15</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterp Tetrahedron Modified basis P=6 Q=7</description>
        <executable>StdInterp</executable>
        <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.89776e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="2e-12">4.32543e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterp Triangle Modified basis P=6 Q=7</description>
        <executable>StdInterp</executable>
        <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7 -P GaussGaussLegendre GaussGaussLegendre </parameters>
        <metrics>
            <metric type="Linf" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="L2" id="2">
                <value tolerance="1e-12">0</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdInterp Triangle Orthonormal basis P=7n Q=8</description>
        <executable>StdInterp</executable>
        <parameters>-s triangle -b Ortho_A Ortho_B -o 7 7 -p 8 8 -P GaussGaussLegendre GaussGaussLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.81216e-16</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.33227e-15</value>
            </metric>
        </metrics>
    </test>
</tests>
