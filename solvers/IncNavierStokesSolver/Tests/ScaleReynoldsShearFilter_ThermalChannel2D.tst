<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Scaling the number of modes in the Reynolds stress filter</description>
        <segment type="sequential">
            <executable>IncNavierStokesSolver</executable>
            <parameters>ScaleReynoldsShearFilter_ThermalChannel2D.xml</parameters>
            <processes> 1 </processes>
        </segment>
        <segment type="sequential">
            <executable>../../utilities/FieldConvert/FieldConvert</executable>
            <parameters>-m printfldnorms ScaleReynoldsShearFilter_ThermalChannel2D.xml ScaleReynoldsShearFilter_ThermalChannel2D_stress.fld stdout</parameters>
            <processes> 1 </processes>
        </segment>
    <files>
        <file description="Session File">ScaleReynoldsShearFilter_ThermalChannel2D.xml</file>
        <file description="Session File">ScaleReynoldsShearFilter_ThermalChannel2D.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.83126</value>
            <value variable="v" tolerance="1e-12">0.0581544</value>
            <value variable="theta" tolerance="1e-12">1.43323</value>
            <value variable="p" tolerance="1e-12">283.566</value>
            <value variable="uu" tolerance="1e-12">2.27537</value>
            <value variable="uv" tolerance="1e-12">0.0127426</value>
            <value variable="utheta" tolerance="1e-12">0.00254648</value>
            <value variable="vv" tolerance="1e-9">3.3398</value>
            <value variable="vtheta" tolerance="1e-9">0.226775</value>
            <value variable="thetatheta" tolerance="1e-9">0.0155594</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.01319</value>
            <value variable="v" tolerance="1e-12">0.0585532</value>
            <value variable="theta" tolerance="1e-12">1</value>
            <value variable="p" tolerance="1e-12">1.6809</value>
            <value variable="uu" tolerance="1e-12">1.46165</value>
            <value variable="uv" tolerance="1e-12">0.0129115</value>
            <value variable="utheta" tolerance="1e-12">0.00229161</value>
            <value variable="vv" tolerance="1e-9">2.43021</value>
            <value variable="vtheta" tolerance="1e-9">0.158068</value>
            <value variable="thetatheta" tolerance="1e-9">0.0103278</value>
        </metric>
    </metrics>
</test>
