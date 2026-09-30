<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D nodal tri and quad vtu output </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e Nodal2DTriQuad.xml Nodal2DTriQuad.fld Nodal2DTriQuad.vtu</parameters>
    <files>
        <file description="Session File">Nodal2DTriQuad.xml</file>
	<file description="Session File">Nodal2DTriQuad.fld</file>
    </files>
     <metrics>
       <metric type="L2" id="1">
         <value variable="x" tolerance="1e-4">13.5641</value>
         <value variable="y" tolerance="1e-4">12.7083</value>
         <value variable="u" tolerance="1e-4">2.63222</value>
       </metric>
     </metrics>
</test>

