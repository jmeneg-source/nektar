<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process fld directory without P0000000.fld </description>
    <executable>FieldConvert</executable>
    <parameters> -e -f naca0012_3D_bnd.xml naca0012_3D_bnd.fld naca0012_3D_bnd.plt </parameters>
    <files>
        <file description="Session File">naca0012_3D_bnd.xml</file>
        <file description="Field File">naca0012_3D_bnd.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-4">1.77227</value>
            <value variable="y" tolerance="1e-4">0.493971</value>
            <value variable="z" tolerance="1e-4">9.2236</value>
            <value variable="Shear_x" tolerance="1e-4">0.312949</value>
            <value variable="Shear_y" tolerance="1e-4">0.32099</value>
            <value variable="Shear_z" tolerance="1e-4">0.0578149</value>
            <value variable="Shear_mag" tolerance="1e-4">0.452696</value>
            <value variable="Norm_x" tolerance="1e-4">0.987202</value>
            <value variable="Norm_y" tolerance="1e-4">3.03874</value>
            <value variable="Norm_z" tolerance="1e-4">4.96156e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-4">0.966252</value>
            <value variable="y" tolerance="1e-4">0.021514</value>
            <value variable="z" tolerance="1e-4">5</value>
            <value variable="Shear_x" tolerance="1e-4">0.725261</value>
            <value variable="Shear_y" tolerance="1e-4">1.11074</value>
            <value variable="Shear_z" tolerance="1e-4">0.808733</value>
            <value variable="Shear_mag" tolerance="1e-4">1.57498</value>
            <value variable="Norm_x" tolerance="1e-4">0.999636</value>
            <value variable="Norm_y" tolerance="1e-4">0.992511</value>
            <value variable="Norm_z" tolerance="1e-4">3.57464e-14</value>
        </metric>
    </metrics>
</test>
