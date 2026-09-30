<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, Check MVM/FFTW consistency</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_MVM_FFTW_Consistency.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_MVM_FFTW_Consistency.xml</file>
        <file description="Restart file">KovaFlow_3DH1D_MVM_FFTW_Consistency.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.19175e-06</value>
            <value variable="v" tolerance="1e-11">1.34946e-06</value>
            <value variable="w" tolerance="1e-11">1.27995e-06</value>
	    <value variable="p" tolerance="1e-10">2.89141e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">3.14595e-06</value>
            <value variable="v" tolerance="1e-11">2.25562e-06</value>
            <value variable="w" tolerance="1e-11">1.81458e-06</value>
	    <value variable="p" tolerance="2e-10">6.87311e-05</value>
        </metric>
    </metrics>
</test>
