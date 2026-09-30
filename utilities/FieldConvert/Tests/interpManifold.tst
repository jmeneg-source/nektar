<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> interpolation on a manifold </description>
    <executable>FieldConvert</executable>
    <parameters> -e -f -v -m interppoints:fromxml=naca0012_3D_bnd.xml:fromfld=naca0012_3D_bnd.fld:topts=naca0012_3D_interp.csv:distTolerance=1E-2 naca0012_3D_interp.dat </parameters>
    <files>
        <file description="Session File">naca0012_3D_bnd.xml</file>
        <file description="Field File">naca0012_3D_bnd.fld</file>
        <file description="Field File">naca0012_3D_interp.csv</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">0.54064</value>
            <value variable="y" tolerance="1e-4">0.150715</value>
            <value variable="z" tolerance="1e-4">0.195195</value>
            <value variable="Shear_x" tolerance="1e-4">0.130942</value>
            <value variable="Shear_y" tolerance="1e-4">0.17828</value>
            <value variable="Shear_z" tolerance="1e-4">0.0128847</value>
            <value variable="Shear_mag" tolerance="1e-4">0.221541</value>
            <value variable="Norm_x" tolerance="1e-4">0.351291</value>
            <value variable="Norm_y" tolerance="1e-4">0.936267</value>
            <value variable="Norm_z" tolerance="1e-4">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">0.966252</value>
            <value variable="y" tolerance="1e-4">0.260036</value>
            <value variable="z" tolerance="1e-4">0.195195</value>
            <value variable="Shear_x" tolerance="1e-4">0.714302</value>
            <value variable="Shear_y" tolerance="1e-4">1.14682</value>
            <value variable="Shear_z" tolerance="1e-4">0.0390763</value>
            <value variable="Shear_mag" tolerance="1e-4">1.22686</value>
            <value variable="Norm_x" tolerance="1e-4">0.999588</value>
            <value variable="Norm_y" tolerance="1e-4">0.999984</value>
            <value variable="Norm_z" tolerance="1e-4">0</value>
        </metric>
    </metrics>
</test>
