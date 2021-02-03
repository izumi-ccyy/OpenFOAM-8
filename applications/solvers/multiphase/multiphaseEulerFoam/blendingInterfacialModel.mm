<map version="freeplane 1.8.0">
<!--To view this file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net -->
<node TEXT="blendingInterfacialModel" FOLDED="false" ID="ID_1031683301" CREATED="1612366302616" MODIFIED="1612366317207" STYLE="oval">
<font SIZE="18"/>
<hook NAME="MapStyle">
    <properties edgeColorConfiguration="#808080ff,#ff0000ff,#0000ffff,#00ff00ff,#ff00ffff,#00ffffff,#7c0000ff,#00007cff,#007c00ff,#7c007cff,#007c7cff,#7c7c00ff" fit_to_viewport="false"/>

<map_styles>
<stylenode LOCALIZED_TEXT="styles.root_node" STYLE="oval" UNIFORM_SHAPE="true" VGAP_QUANTITY="24.0 pt">
<font SIZE="24"/>
<stylenode LOCALIZED_TEXT="styles.predefined" POSITION="right" STYLE="bubble">
<stylenode LOCALIZED_TEXT="default" ICON_SIZE="12.0 pt" COLOR="#000000" STYLE="fork">
<font NAME="SansSerif" SIZE="10" BOLD="false" ITALIC="false"/>
</stylenode>
<stylenode LOCALIZED_TEXT="defaultstyle.details"/>
<stylenode LOCALIZED_TEXT="defaultstyle.attributes">
<font SIZE="9"/>
</stylenode>
<stylenode LOCALIZED_TEXT="defaultstyle.note" COLOR="#000000" BACKGROUND_COLOR="#ffffff" TEXT_ALIGN="LEFT"/>
<stylenode LOCALIZED_TEXT="defaultstyle.floating">
<edge STYLE="hide_edge"/>
<cloud COLOR="#f0f0f0" SHAPE="ROUND_RECT"/>
</stylenode>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.user-defined" POSITION="right" STYLE="bubble">
<stylenode LOCALIZED_TEXT="styles.topic" COLOR="#18898b" STYLE="fork">
<font NAME="Liberation Sans" SIZE="10" BOLD="true"/>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.subtopic" COLOR="#cc3300" STYLE="fork">
<font NAME="Liberation Sans" SIZE="10" BOLD="true"/>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.subsubtopic" COLOR="#669900">
<font NAME="Liberation Sans" SIZE="10" BOLD="true"/>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.important">
<icon BUILTIN="yes"/>
</stylenode>
</stylenode>
<stylenode LOCALIZED_TEXT="styles.AutomaticLayout" POSITION="right" STYLE="bubble">
<stylenode LOCALIZED_TEXT="AutomaticLayout.level.root" COLOR="#000000" STYLE="oval" SHAPE_HORIZONTAL_MARGIN="10.0 pt" SHAPE_VERTICAL_MARGIN="10.0 pt">
<font SIZE="18"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,1" COLOR="#0033ff">
<font SIZE="16"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,2" COLOR="#00b439">
<font SIZE="14"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,3" COLOR="#990000">
<font SIZE="12"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,4" COLOR="#111111">
<font SIZE="10"/>
</stylenode>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,5"/>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,6"/>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,7"/>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,8"/>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,9"/>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,10"/>
<stylenode LOCALIZED_TEXT="AutomaticLayout.level,11"/>
</stylenode>
</stylenode>
</map_styles>
</hook>
<hook NAME="AutomaticEdgeColor" COUNTER="2" RULE="ON_BRANCH_CREATION"/>
<node TEXT="blendingInterfacialModel.H" POSITION="right" ID="ID_859280618" CREATED="1612366331509" MODIFIED="1612366345106">
<edge COLOR="#ff0000"/>
<node TEXT="private data" ID="ID_605005455" CREATED="1612366365393" MODIFIED="1612366368619">
<node TEXT="phase1_" ID="ID_451639107" CREATED="1612366370114" MODIFIED="1612366394019">
<node TEXT="Reference to phase 1" ID="ID_1832646816" CREATED="1612366452189" MODIFIED="1612366456970"/>
</node>
<node TEXT="phase2_" ID="ID_1652929615" CREATED="1612366395882" MODIFIED="1612366401324">
<node TEXT="Reference to phase 2" ID="ID_1847156470" CREATED="1612366523446" MODIFIED="1612366523446"/>
</node>
<node TEXT="blending_" ID="ID_665910182" CREATED="1612366402518" MODIFIED="1612366408967">
<node TEXT="Blending model" ID="ID_1238418097" CREATED="1612366526612" MODIFIED="1612366526612"/>
</node>
<node TEXT="model_" ID="ID_949995368" CREATED="1612366409289" MODIFIED="1612366411852">
<node TEXT="Model for region with no obvious dispersed phase" ID="ID_1855208482" CREATED="1612366529883" MODIFIED="1612366529883"/>
</node>
<node TEXT="model1In2_" ID="ID_654968193" CREATED="1612366412606" MODIFIED="1612366420052">
<node TEXT="Model for dispersed phase 1 in continuous phase 2" ID="ID_1487787915" CREATED="1612366534132" MODIFIED="1612366534132"/>
</node>
<node TEXT="model2In1_" ID="ID_1216579370" CREATED="1612366420444" MODIFIED="1612366424315">
<node TEXT="Model for dispersed phase 2 in continuous phase 1" ID="ID_282465810" CREATED="1612366538415" MODIFIED="1612366538415"/>
</node>
<node TEXT="correctFixedFluxBCs_" ID="ID_1453667251" CREATED="1612366426285" MODIFIED="1612366436867">
<node TEXT="If true set coefficients and forces to 0 at fixed-flux BCs" ID="ID_1334281608" CREATED="1612366542815" MODIFIED="1612366542815"/>
</node>
</node>
<node TEXT="private member functions" FOLDED="true" ID="ID_1977576819" CREATED="1612366379805" MODIFIED="1612366386219">
<node TEXT="calculateBlendingCoeffs()" ID="ID_1631595684" CREATED="1612366546189" MODIFIED="1612366556839">
<node TEXT="Calculate the blending coefficients" ID="ID_535011069" CREATED="1612366680219" MODIFIED="1612366680219"/>
</node>
<node TEXT="correctFixedFluxBCs()" ID="ID_613674019" CREATED="1612366557169" MODIFIED="1612366579726">
<node TEXT="Correct coeff/value on fixed flux boundary conditions" ID="ID_989408058" CREATED="1612366675851" MODIFIED="1612366675851"/>
</node>
<node TEXT="evaluate()" ID="ID_1944605532" CREATED="1612366580896" MODIFIED="1612366653571"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions, one for tmp, one for hashPtrTable
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Return the blended coeff/value" ID="ID_83585701" CREATED="1612366665799" MODIFIED="1612366665799"/>
<node TEXT="Return the blended coeff/value" ID="ID_1040002250" CREATED="1612366667481" MODIFIED="1612366667481"/>
</node>
</node>
<node TEXT="public" ID="ID_1386505012" CREATED="1612366879350" MODIFIED="1612366881587">
<node TEXT="Typename" ID="ID_240660860" CREATED="1612366883761" MODIFIED="1612366885995"/>
<node TEXT="constructors" FOLDED="true" ID="ID_1609896155" CREATED="1612366886630" MODIFIED="1612366932615"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Construct from two phases, blending method and three models" ID="ID_1693379648" CREATED="1612366921380" MODIFIED="1612366921380"/>
<node TEXT="Construct from the model table, dictionary and pairs" ID="ID_1401120077" CREATED="1612366925948" MODIFIED="1612366925948"/>
</node>
<node TEXT="destructor" FOLDED="true" ID="ID_1539915016" CREATED="1612366891625" MODIFIED="1612366898158">
<node TEXT="~BlendedInterfacialModel()" ID="ID_226096439" CREATED="1612366904304" MODIFIED="1612366904304"/>
</node>
<node TEXT="member functions" ID="ID_1916436222" CREATED="1612366939862" MODIFIED="1612366946187">
<node TEXT="K()" ID="ID_1420115795" CREATED="1612366948746" MODIFIED="1612366963695">
<node TEXT="Return the blended force coefficient" ID="ID_761377839" CREATED="1612367047932" MODIFIED="1612367047932"/>
</node>
<node TEXT="K(residulaAlpha)" ID="ID_1630784442" CREATED="1612366964154" MODIFIED="1612366975112">
<node TEXT="Return the blended force coefficient with a specified residual alpha" ID="ID_1257676744" CREATED="1612367056028" MODIFIED="1612367056028"/>
</node>
<node TEXT="Kf()" ID="ID_1783900670" CREATED="1612366975733" MODIFIED="1612366988888">
<node TEXT="Return the face blended force coefficient" ID="ID_1688762778" CREATED="1612367060600" MODIFIED="1612367060600"/>
</node>
<node TEXT="F()" ID="ID_610803990" CREATED="1612366989198" MODIFIED="1612366994820">
<node TEXT="Return the blended force" ID="ID_1660312374" CREATED="1612367064592" MODIFIED="1612367064592"/>
</node>
<node TEXT="Ff()" ID="ID_1411472740" CREATED="1612366995126" MODIFIED="1612366998239">
<node TEXT="Return the face blended force" ID="ID_580137829" CREATED="1612367069594" MODIFIED="1612367069594"/>
</node>
<node TEXT="D()" ID="ID_1320837781" CREATED="1612366998506" MODIFIED="1612367001479">
<node TEXT="Return the blended diffusivity" ID="ID_1428266291" CREATED="1612367073764" MODIFIED="1612367073764"/>
</node>
<node TEXT="mixture()" ID="ID_1886359815" CREATED="1612367001794" MODIFIED="1612367004595">
<node TEXT="Return the list of individual species that are transferred" ID="ID_47441601" CREATED="1612367079027" MODIFIED="1612367079027"/>
</node>
<node TEXT="dmdtf()" ID="ID_1651588002" CREATED="1612367004925" MODIFIED="1612367008227">
<node TEXT="Return the blended mass transfer rate" ID="ID_1894664779" CREATED="1612367083261" MODIFIED="1612367083261"/>
</node>
<node TEXT="species()" ID="ID_756221333" CREATED="1612367008449" MODIFIED="1612367011535">
<node TEXT="Return the list of individual species that are transferred" ID="ID_165838053" CREATED="1612367087359" MODIFIED="1612367087359"/>
</node>
<node TEXT="dmidtf()" ID="ID_1408638653" CREATED="1612367011721" MODIFIED="1612367034395">
<node TEXT="Return the blended mass transfer rates for individual species" ID="ID_125204099" CREATED="1612367092103" MODIFIED="1612367092103"/>
</node>
<node TEXT="writeData()" ID="ID_1695326481" CREATED="1612367034641" MODIFIED="1612367042971">
<node TEXT="Dummy write for regIOobject" ID="ID_1308926824" CREATED="1612367095984" MODIFIED="1612367095984"/>
</node>
</node>
<node TEXT="member operator" ID="ID_970507252" CREATED="1612367107872" MODIFIED="1612367111351">
<node TEXT="operator=()" ID="ID_934683240" CREATED="1612367124925" MODIFIED="1612367135926">
<node TEXT="Disallow default bitwise assignment" ID="ID_1638869144" CREATED="1612367137472" MODIFIED="1612367137472"/>
</node>
</node>
</node>
</node>
<node TEXT="blendingInterfacialModel.C" POSITION="right" ID="ID_304549322" CREATED="1612366353282" MODIFIED="1612366361471">
<edge COLOR="#0000ff"/>
</node>
</node>
</map>
