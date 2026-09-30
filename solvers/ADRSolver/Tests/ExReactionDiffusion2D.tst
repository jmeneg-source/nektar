<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D Explicit Reaction Diffusion with bodyforce</description>
    <executable>ADRSolver</executable>
    <parameters>ExReactionDiffusion2D.xml</parameters>
    <files>
        <file description="Session File">ExReactionDiffusion2D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.04893e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.39248e-08</value>
        </metric>
    </metrics>
</test>
