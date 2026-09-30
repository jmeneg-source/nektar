<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test the zero-homo-plane FieldConvert module</description>
    <executable>FieldConvert</executable>
    <parameters>-m zero-homo-plane:planeid=0 square_mix.xml session_zero-homo-plane.xml zero-homo-plane_pre.fld zero-homo-plane_post.fld -f -e</parameters>
    <files>
        <file description="Session File">square_mix.xml</file>
        <file description="Session File">session_zero-homo-plane.xml</file>
        <file description="Session File">zero-homo-plane_pre.fld</file>
    </files>
    <metrics>
        
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1.63299</value>
            <value variable="y" tolerance="1e-6">1.63299</value>
            <value variable="z" tolerance="1e-6">0.661438</value>
            <value variable="u" tolerance="1e-6">0.0</value>
            <value variable="v" tolerance="1e-6">0.0</value>
            <value variable="w" tolerance="1e-6">0.0</value>
            <value variable="p" tolerance="1e-6">0.0</value> 
        </metric>

        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">2.0</value>
            <value variable="y" tolerance="1e-6">2.0</value>
            <value variable="z" tolerance="1e-6">0.75</value>
            <value variable="u" tolerance="1e-6">0.0</value>
            <value variable="v" tolerance="1e-6">0.0</value>
            <value variable="w" tolerance="1e-6">0.0</value>
            <value variable="p" tolerance="1e-6">0.0</value>
        </metric>
        
    </metrics>
</test>
