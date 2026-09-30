<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit tests for the SessionReader class</description>
    <executable python="true">SessionReader_UnitTest.py</executable>
    <parameters>-v</parameters>
    <files>
        <file description="Session File">../../MultiRegions/Tests/newsquare_2x2.xml</file>
    </files>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testReadParamsMap</function>
            <function>testReadVarsList</function>
        </metric>
    </metrics>
</test>