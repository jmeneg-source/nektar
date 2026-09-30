<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process halfmodetofourier with w being transformed to imaginary mode </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m halfmodetofourier:realmodetoimag=2  TriQuad.xml TriQuad.fld  TriQuad.plt</parameters>
    <files>
        <file description="Session File">TriQuad.xml</file>
	<file description="input fld File">TriQuad.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x"   tolerance="1e-05">1.15470</value>
            <value variable="y"   tolerance="1e-05">1.15470</value>
            <value variable="z"   tolerance="1e-06">0.935414</value>
            <value variable="u"   tolerance="1e-05">3.00785</value>
            <value variable="v"   tolerance="1e-06">0.643477</value>
            <value variable="w"   tolerance="1e-06">0.816497</value>
            <value variable="p"   tolerance="1e-05">1.91485</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x"   tolerance="1e-05">1.</value>
            <value variable="y"   tolerance="1e-05">1.</value>
            <value variable="z"   tolerance="1e-06">0.75</value>
            <value variable="u"   tolerance="1e-05">6.28319</value>
            <value variable="v"   tolerance="1e-05">1.</value>
            <value variable="w"   tolerance="1e-05">1.</value>
            <value variable="p"   tolerance="1e-05">4.</value>
        </metric>
    </metrics>
</test>

