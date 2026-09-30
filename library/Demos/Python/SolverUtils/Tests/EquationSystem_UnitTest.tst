<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the SolverUtils::EquationSystem class.</description>
    <executable python="true">EquationSystem_UnitTest.py</executable>
    <parameters>-v TriQuadChannel.xml</parameters>
    <files>
        <file description="Session File">../../../../../solvers/IncNavierStokesSolver/Tests/TriQuadChannel.xml</file>
    </files>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testDoInitialise</function>
            <function>testDoSolve</function>
            <function>testEvaluateExactSolution</function>
            <function>testGetFields</function>
            <function>testGetNpoints</function>
            <function>testGetNvariables</function>
            <function>testInit</function>
            <function>testSteps</function>
            <function>testTime</function>
            <function>testTimeStep</function>
        </metric>
    </metrics>
</test>
