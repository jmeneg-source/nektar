<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Courtemanche Cell model</description>
    <executable>PrePacing</executable>
    <parameters>CellMLTest.xml</parameters>
    <files>
        <file description="Session File">CellMLTest.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>
                ^#\s([\w]*)\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)
            </regex>
            <matches>
                <match>
                    <field>V</field>
                    <field tolerance="1e-06">-8.09679</field>
                </match>
                <match>
                    <field>m</field>
                    <field tolerance="1e-06">0.987074</field>
                </match>
                <match>
                    <field>h</field>
                    <field tolerance="1e-06">2.03226e-178</field>
                </match>
                <match>
                    <field>j</field>
                    <field tolerance="1e-06">3.67039e-12</field>
                </match>
                <match>
                    <field>oa</field>
                    <field tolerance="1e-06">0.670337</field>
                </match>
                <match>
                    <field>oi</field>
                    <field tolerance="1e-06">0.00186459</field>
                </match>
                <match>
                    <field>ua</field>
                    <field tolerance="1e-06">0.910572</field>
                </match>
                <match>
                    <field>ui</field>
                    <field tolerance="1e-06">0.993321</field>
                </match>
                <match>
                    <field>xr</field>
                    <field tolerance="1e-06">0.211011</field>
                </match>
                <match>
                    <field>xs</field>
                    <field tolerance="1e-06">0.0843768</field>
                </match>
                <match>
                    <field>d</field>
                    <field tolerance="1e-06">0.562857</field>
                </match>
                <match>
                    <field>f</field>
                    <field tolerance="1e-06">0.678725</field>
                </match>
                <match>
                    <field>f_Ca</field>
                    <field tolerance="1e-06">0.344585</field>
                </match>
                <match>
                    <field>U</field>
                    <field tolerance="1e-06">0.000418484</field>
                </match>
                <match>
                    <field>V</field>
                    <field tolerance="1e-06">2.5819e-11</field>
                </match>
                <match>
                    <field>W</field>
                    <field tolerance="1e-06">0.944047</field>
                </match>
                <match>
                    <field>Na_i</field>
                    <field tolerance="1e-06">11.1714</field>
                </match>
                <match>
                    <field>Ca_i</field>
                    <field tolerance="1e-06">0.000661552</field>
                </match>
                <match>
                    <field>K_i</field>
                    <field tolerance="1e-06">138.991</field>
                </match>
                <match>
                    <field>Ca_rel</field>
                    <field tolerance="1e-06">0.204121</field>
                </match>
                <match>
                    <field>Ca_up</field>
                    <field tolerance="1e-06">1.58659</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>




