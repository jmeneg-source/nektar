<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the SolverUtils::Filter class.</description>
    <executable python="true">Filter_UnitTest.py</executable>
    <parameters>-v Tet_channel_m3.xml</parameters>
    <files>
        <file description="Session File">../../../../../solvers/IncNavierStokesSolver/Tests/Tet_channel_m3.xml</file>
    </files>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testCreateFilter</function>
            <function>testGetFields</function>
            <function>testInheritance</function>
            <function>testRunFilter</function>
            <function>testUnsuppliedConfig</function>
        </metric>
    </metrics>
</test>
