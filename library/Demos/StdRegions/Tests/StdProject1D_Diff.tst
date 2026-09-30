<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdProject_Diff1D Segment Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s segment -b Modified_A -o 6 -p 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">1.8355e-14</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">4.9738e-14</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject_Diff1D Segment Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s Segment -b Ortho_A -o 6 -p 7 -d</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">3.65229e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">4.88498e-15</value>
            </metric>
        </metrics>
    </test>
</tests>
