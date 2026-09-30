<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction with local p-refinement using mixed TYPE methods </description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_Cube_pRef_MixedType.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Cube_pRef_MixedType.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9"> 2.65743e-05 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9"> 0.000118796 </value>
        </metric>
    </metrics>
</test>
