<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady CG implicit diffusion, variable coeffs. </description>
    <executable>ADRSolver</executable>
    <parameters> ImDiffusion_VarCoeff.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_VarCoeff.xml </file>
    </files>
    <metrics>
   <metric type="L2" id="1">
	   <value variable="u" tolerance="1e-10">7.57494e-08 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">1.21038e-07 </value>
        </metric>
    </metrics>
</test>
