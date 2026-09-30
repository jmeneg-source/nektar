<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Test of creating and writing out meshes using NekPy </description>
    <executable python="true"> CreateMesh.py </executable>
    <parameters></parameters>
    <files></files>
    <metrics>
	<metric type="regex" id="1">
	    <regex>^.*Test (.*)</regex>
            <matches>
                <match>
                    <field id="0">successful!</field>
                </match>
            </matches>			
	</metric>
    </metrics>
</test>
