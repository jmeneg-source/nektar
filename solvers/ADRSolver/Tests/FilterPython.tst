<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test of the Python filter</description>
    <executable>ADRSolver</executable>
    <parameters>-v FilterPython.xml</parameters>
    <files>
        <file description="Session File">FilterPython.xml</file>
        <file description="Session File">FilterPython_Function.py</file>
        <file description="Session File">FilterPython_Class.py</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^TestVTKFilter: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">200</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="2">
            <regex>^FilterPython_Function: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">200</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
