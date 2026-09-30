<?xml version="1.0" encoding="utf-8"?>
<tests>
    <test>
        <description>StdProject1D Segment Modified basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s segment -b Modified_A -o 6 -p 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">5.07435e-15</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">5.10703e-15</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject1D Segment Orthonormal basis P=6 Q=7</description>
        <executable>StdProject</executable>
        <parameters>-s Segment -b Ortho_A -o 6 -p 7</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">5.37715e-16</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">1.16573e-15</value>
            </metric>
        </metrics>
    </test>
    <test>
        <description>StdProject1D Segment single mode Fourier basis P=2 Q=2</description>
        <executable>StdProject</executable>
        <parameters>-s segment -b FourierSingleMode -o 2 -p 2</parameters>
        <metrics>
            <metric type="L2" id="1">
                <value tolerance="1e-12">0</value>
            </metric>
            <metric type="Linf" id="2">
                <value tolerance="1e-12">2.77556e-17</value>
            </metric>
        </metrics>
    </test>
</tests>
