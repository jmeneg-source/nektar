<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdProject_Diff3D Hexahedron Chebyshev basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b chebyshev chebyshev chebyshev -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-8">5.62024e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-8">1.20508e-11</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Hexahedron Lagrange basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b GLL_Lagrange GLL_Lagrange GLL_Lagrange -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">5.11314e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-11">1.11413e-11</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Hexahedron Legendre basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b legendre legendre legendre -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">4.88883e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="4e-12">1.13687e-11</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Hexahedron Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-8">3.12754e-10</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-8">5.69116e-10</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Hexahedron Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s hexahedron -b Ortho_A Ortho_A Ortho_A -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">4.88883e-13</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="4e-12">1.13687e-11</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Prism Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s prism -b Modified_A Modified_A Modified_B -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-8">1.52731e-11</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-8">4.03908e-11</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Prism Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s prism -b Ortho_A Ortho_A Ortho_B -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">2.56024e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.98428e-13</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Pyramid Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 6 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-9">4.99628e-12</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-9">2.5608e-11</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Tetrahedron Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-11">2.1782e-12</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-11">5.08749e-12</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff3D Tetrahedron Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s tetrahedron -b Ortho_A Ortho_B Ortho_C -o 6 6 6 -p 7 7 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">8.49152e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">9.23706e-14</value>
            </metric>
        </metrics>
    </test>
</tests>
