<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D  output with a range restriction and  hdf5 fld  files with variable p</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e KovaFlow_varP.xml:xml:range="0,1,0.1,0.9"  KovaFlow_varP_hdf5.fld KovaFlow_varP.plt</parameters>
    <files>
	<file description="Mesh File">KovaFlow_varP.xml</file>
	<file description="Field File">KovaFlow_varP_hdf5.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.612372</value>
            <value variable="y" tolerance="1e-6">0.707107</value>
            <value variable="u" tolerance="1e-6">1.46348</value>
            <value variable="v" tolerance="1e-8">0.122973</value>
            <value variable="p" tolerance="1e-7">0.420724</value>
        </metric>
    </metrics>
</test>

