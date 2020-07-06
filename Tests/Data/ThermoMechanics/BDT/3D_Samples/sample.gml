<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>sample</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="p_1"/>
        <point id="1" x="0.02" y="0" z="0" name="p_2"/>
        <point id="2" x="0.02" y="0.01" z="0" name="p_3"/>
        <point id="3" x="0" y="0.01" z="0"  name="p_4"/>

        <point id="4" x="0" y="0" z="0.001" name="p_5"/>
        <point id="5" x="0.02" y="0" z="0.001" name="p_6"/>
        <point id="6" x="0.02" y="0.01" z="0.001" name="p_7"/>
        <point id="7" x="0" y="0.01" z="0.001"  name="p_8"/>
        
        <point id="8" x="0.02" y="0.005" z="0"  name="p_9"/>
        <point id="9" x="0.02" y="0.005" z="0.001"  name="p_10"/>

        <point id="10" x="0.02" y="0.0" z="0.0005"  name="p_11"/>
        <point id="11" x="0.02" y="0.01" z="0.0005"  name="p_12"/>


    </points>

    <surfaces>
        <surface id="0" name="left">
            <element p1="0" p2="4" p3="3"/>
            <element p1="4" p2="3" p3="7"/>
        </surface>
        <surface id="1" name="right">
            <element p1="1" p2="5" p3="2"/>
            <element p1="5" p2="2" p3="6"/>
        </surface>
        <surface id="2" name="top">
            <element p1="3" p2="2" p3="7"/>
            <element p1="2" p2="7" p3="6"/>
        </surface>
        <surface id="3" name="bottom">
            <element p1="0" p2="1" p3="4"/>
            <element p1="1" p2="4" p3="5"/>
        </surface>
        <surface id="4" name="front">
            <element p1="0" p2="1" p3="3"/>
            <element p1="1" p2="3" p3="2"/>
        </surface>
        <surface id="5" name="back">
            <element p1="4" p2="5" p3="7"/>
            <element p1="5" p2="7" p3="6"/>
        </surface>
    </surfaces>

    <polylines>
    <polyline id="0" name="fix_y">
        <pnt>8</pnt>
        <pnt>9</pnt>
    </polyline>
    <polyline id="1" name="fix_z">
        <pnt>10</pnt>
        <pnt>11</pnt>
    </polyline>

    </polylines>


    <!--polylines>
    <polyline id="0" name="top">
        <pnt>2</pnt>
        <pnt>3</pnt>
    </polyline>
    <polyline id="1" name="bottom">
        <pnt>0</pnt>
        <pnt>1</pnt>
    </polyline>
    <polyline id="2" name="left">
        <pnt>3</pnt>
        <pnt>0</pnt>
    </polyline>
    <polyline id="3" name="right">
        <pnt>1</pnt>
        <pnt>2</pnt>
    </polyline>
    </polylines-->

</OpenGeoSysGLI>
