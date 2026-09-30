<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdProject2D Quadrilateral Fourier basis P=6 Q=8</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b Fourier Fourier -o 6 6 -p 8 8</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.06892e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.77636e-15</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Quadrilateral Fourier Single Mode basis P=2 Q=2</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b FourierSingleMode FourierSingleMode -o 2 2 -p 2 2</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">3.59678e-16</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">8.88178e-16</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Quadrilateral Lagrange basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">4.62859e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.13163e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Quadrilateral Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b Modified_A Modified_A -o 6 6 -p 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">5.08794e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="2e-12">1.1795e-12</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Quadrilateral Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b Ortho_A Ortho_A -o 6 6 -p 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.85686e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.06581e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Triangle Modified basis P=6 Q=7, alternative quadrature weights</description>
        <executable>StdProject</executable>
        <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7 -P GaussLobattoLegendre GaussLobattoLegendre</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">3.24325e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">3.86358e-14</value>
            </metric>
            <metric type="Regex" id="3">
                <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
                <matches>
                    <match>
                        <field id="0" tolerance="1e-11">3.99614e-12</field>
                    </match>
                </matches>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Triangle Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.18252e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">5.15143e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Triangle Nodal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-n NodalTriElec -o 6 6 -p 7 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.78107e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.5099e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject2D Triangle Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s triangle -b Ortho_A Ortho_B -o 6 6 -p 7 7</parameters>
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
