<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, adaptive P, 16 Fourier modes, using explicit mapping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_adaptive_16modes_FFTW_Mapping.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_adaptive_16modes_FFTW_Mapping.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">5.96501e-05</value>
            <value variable="v" tolerance="1e-10">1.65234e-05</value>
            <value variable="w" tolerance="1e-10">1.78794e-05</value>
	    <value variable="p" tolerance="1e-09">0.000368607</value>
        </metric>
        <metric type="Linf" id="2">
	    <value variable="u" tolerance="1e-10">9.04205e-05</value>
            <value variable="v" tolerance="1e-10">2.48215e-05</value>
            <value variable="w" tolerance="1e-10">2.89005e-05</value>
	    <value variable="p" tolerance="1e-09">0.000805191</value>
        </metric>
    </metrics>
</test>
