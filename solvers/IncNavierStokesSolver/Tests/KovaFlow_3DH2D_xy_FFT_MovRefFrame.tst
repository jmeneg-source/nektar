<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 2D, flow in xy, MovingReferenceFrame, FFT</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH2D_xy_FFT_MovRefFrame.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH2D_xy_FFT_MovRefFrame.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="3.5e-5">3.41544e-05</value>
            <value variable="v" tolerance="4.46e-5">4.45434e-05</value>
            <value variable="w" tolerance="1e-9">2.11462e-17</value>
	    <value variable="p" tolerance="5e-8">3.95241e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="7.1e-5">5.2044e-05</value>
            <value variable="v" tolerance="1.1e-5">0.000107335</value>
            <value variable="w" tolerance="1e-9">9.28156e-17</value>
            <value variable="p" tolerance="5e-8">7.20077e-05</value>
        </metric>
    </metrics>
</test>
