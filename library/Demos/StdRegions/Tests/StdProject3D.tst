<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdProject3D Hexahedron Chebyshev basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b chebyshev chebyshev chebyshev -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">4.11256e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">4.1743e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Hexahedron Lagrange basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b GLL_Lagrange GLL_Lagrange GLL_Lagrange -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">3.06382e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">3.41061e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Hexahedron Legendre basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b legendre legendre legendre -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.20528e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.84217e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Hexahedron Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-8">4.76767e-11</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-8">8.0945e-11</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Hexahedron Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b Ortho_A Ortho_A Ortho_A -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.20528e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.84217e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Prism Modified basis P=6 Q=7 with alternative weights</description>
        <executable>StdProject</executable>
        <parameters>-s prism -b Modified_A Modified_A Modified_B -o 6 6 6 -p 7 7 7 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-11">2.08237e-12</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-11">4.65228e-12</value>
            </metric>
            <metric type="Regex" id="3">
                <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
                <matches>
                    <match>
                        <field id="0" tolerance="1e-11">7.38076e-13</field>
                    </match>
                </matches>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Prism Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s prism -b Modified_A Modified_A Modified_B -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-11">2.08237e-12</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-11">4.65228e-12</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Prism Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
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
        <description>StdProject3D Pyramid Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 7 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-11">8.8186e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-11">2.41681e-12</value>
            </metric>
            <metric type="Regex" id="3">
                <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
                <matches>
                    <match>
                        <field id="0" tolerance="1e-11">2.49289e-12</field>
                    </match>
                </matches>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Pyramid Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 6</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-9">3.2756e-12</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-9">9.72786e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Pyramid Orthogonal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s pyramid -b Ortho_A Ortho_A OrthoPyr_C -o 6 6 6 -p 7 7 6</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-9">5.32907e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-9">5.32907e-15</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Tetrahedron Modified basis P=6 Q=7, alternative quadrature weights</description>
        <executable>StdProject</executable>
        <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.89257e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="4e-12">7.77156e-13</value>
            </metric>
            <metric type="Regex" id="3">
                <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
                <matches>
                    <match>
                        <field id="0" tolerance="1e-11">2.15439e-13</field>
                    </match>
                </matches>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Tetrahedron Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.89776e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="2e-12">2.55795e-12</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject3D Tetrahedron Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s tetrahedron -b Ortho_A Ortho_B Ortho_C -o 6 6 6 -p 7 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">6.072e-16</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.66454e-15</value>
            </metric>
        </metrics>
    </test>
</tests>
