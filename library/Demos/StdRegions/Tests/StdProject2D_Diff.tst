<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdProject_Diff2D Quadrilateral Fourier basis P=6 Q=8</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b Fourier Fourier -o 6 6 -p 8 8 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.54556e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.84217e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff2D Quadrilateral Lagrange basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">4.14376e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.13687e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff2D Quadrilateral Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b Modified_A Modified_A -o 6 6 -p 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-11">1.2852e-12</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-11">6.05382e-12</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff2D Quadrilateral Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s quadrilateral -b Ortho_A Ortho_A -o 6 6 -p 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">3.27917e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.13687e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff2D Triangle Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.0808e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.74891e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff2D Triangle Nodal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-n NodalTriElec -o 6 6 -p 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.06405e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">9.59233e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff2D Triangle Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s triangle -b Ortho_A Ortho_B -o 6 6 -p 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">6.45129e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">4.9738e-14</value>
            </metric>
        </metrics>
    </test>
</tests>
