<map version="freeplane 1.8.0">
<!--To view this file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net -->
<node TEXT="phaseSystem.H" FOLDED="false" ID="ID_1450468438" CREATED="1611134120853" MODIFIED="1611134143945" STYLE="oval">
<font SIZE="18"/>
<hook NAME="MapStyle">
    <properties edgeColorConfiguration="#808080ff,#ff0000ff,#0000ffff,#00ff00ff,#ff00ffff,#00ffffff,#7c0000ff,#00007cff,#007c00ff,#7c007cff,#007c7cff,#7c7c00ff" fit_to_viewport="false" show_note_icons="true"/>

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
<hook NAME="AutomaticEdgeColor" COUNTER="7" RULE="ON_BRANCH_CREATION"/>
<node TEXT="public" POSITION="right" ID="ID_1861990465" CREATED="1611134157309" MODIFIED="1611134512461">
<edge COLOR="#ff0000"/>
<node TEXT="public typedefs" ID="ID_1911934725" CREATED="1611134233208" MODIFIED="1611134524776">
<node TEXT="momentumTransferTable" ID="ID_1815302653" CREATED="1611134246985" MODIFIED="1611134260050"/>
<node TEXT="heatTransferTable" ID="ID_1169473602" CREATED="1611134261105" MODIFIED="1611134266311"/>
<node TEXT="specieTransferTable" ID="ID_413537823" CREATED="1611134267276" MODIFIED="1611134279047"/>
<node TEXT="phaseModelList" ID="ID_1008170813" CREATED="1611134279820" MODIFIED="1611134289041"/>
<node TEXT="phaseModelPartialList" ID="ID_1414880174" CREATED="1611134299204" MODIFIED="1611134314970"/>
<node TEXT="phasePairTable" ID="ID_844615315" CREATED="1611134315381" MODIFIED="1611134329655"/>
<node TEXT="dmdtTable" ID="ID_1992912802" CREATED="1611134332302" MODIFIED="1611134336837"/>
<node TEXT="dmidtfTable" ID="ID_1444756432" CREATED="1611134337639" MODIFIED="1611134357470"/>
</node>
</node>
<node TEXT="protected" POSITION="right" ID="ID_101882061" CREATED="1611134182789" MODIFIED="1611134517941">
<edge COLOR="#0000ff"/>
<node TEXT="protected typedefs" ID="ID_1617093736" CREATED="1611134371147" MODIFIED="1611134527711">
<node TEXT="dictTable" ID="ID_1958986782" CREATED="1611134396400" MODIFIED="1611134408601"/>
<node TEXT="blendingMethodTable" ID="ID_1355043544" CREATED="1611134381120" MODIFIED="1611134394621"/>
<node TEXT="surfaceTensionModelTable" ID="ID_244522032" CREATED="1611134412811" MODIFIED="1611134423738"/>
<node TEXT="aspectRatioModelTable" ID="ID_282704732" CREATED="1611134430200" MODIFIED="1611134443689"/>
<node TEXT="cAlphaTable" ID="ID_1343719262" CREATED="1611134444505" MODIFIED="1611134450905"/>
</node>
<node TEXT="protected data" ID="ID_45133563" CREATED="1611134493002" MODIFIED="1611134530161">
<node TEXT="mesh_" ID="ID_1167191567" CREATED="1611134540588" MODIFIED="1611134544939"/>
<node TEXT="referencePhaseName_" ID="ID_1246019262" CREATED="1611134545509" MODIFIED="1611134568557"/>
<node TEXT="phaseModels_" ID="ID_578022762" CREATED="1611134568800" MODIFIED="1611134573205"/>
<node TEXT="movingPhaseModels_" ID="ID_1324280705" CREATED="1611134573384" MODIFIED="1611134581767"/>
<node TEXT="stationaryPhaseModels_" ID="ID_529755993" CREATED="1611134582261" MODIFIED="1611134591181"/>
<node TEXT="anisothermalPhaseModels_" ID="ID_1494237370" CREATED="1611134591472" MODIFIED="1611134610667"/>
<node TEXT="multiComponentPhaseModels_" ID="ID_1782781722" CREATED="1611134611076" MODIFIED="1611134623481"/>
<node TEXT="phasePairs_" ID="ID_1320197260" CREATED="1611134624021" MODIFIED="1611134632490"/>
<node TEXT="phi_" ID="ID_710148763" CREATED="1611134632756" MODIFIED="1611134636623">
<node TEXT="Total volumetric flux" ID="ID_946382642" CREATED="1611134658988" MODIFIED="1611134662109"/>
</node>
<node TEXT="dpdt_" ID="ID_1485099577" CREATED="1611134640904" MODIFIED="1611134643894"/>
<node TEXT="MRF_" ID="ID_1412774964" CREATED="1611134645384" MODIFIED="1611134647950"/>
<node TEXT="blendingMethods_" ID="ID_312556119" CREATED="1611134648293" MODIFIED="1611134675941"/>
<node TEXT="cAlphas_" FOLDED="true" ID="ID_773775563" CREATED="1611134677220" MODIFIED="1611134681229">
<node TEXT="Interface compression coefficients" ID="ID_973304035" CREATED="1611134704212" MODIFIED="1611134711621"/>
</node>
<node TEXT="deltaN_" FOLDED="true" ID="ID_484256301" CREATED="1611134718841" MODIFIED="1611134733457">
<node TEXT="Stabilisation for normalisation of the interface normal" ID="ID_1472530693" CREATED="1611134735057" MODIFIED="1611134741053"/>
</node>
<node TEXT="surfaceTensionModels_" ID="ID_696949979" CREATED="1611134792421" MODIFIED="1611134802235"/>
<node TEXT="aspectRatioModels_" ID="ID_1902449267" CREATED="1611134803424" MODIFIED="1611134808898"/>
<node TEXT="fillFields_" ID="ID_1404050848" CREATED="1611134817589" MODIFIED="1611134830810"/>
</node>
<node TEXT="protected member functions" ID="ID_1271932990" CREATED="1611134929925" MODIFIED="1611134940321">
<node TEXT="calcPhi()" FOLDED="true" ID="ID_36785264" CREATED="1611134943209" MODIFIED="1611140149388" BACKGROUND_COLOR="#ffff99">
<node TEXT="Calculate and return the mixture flux" ID="ID_936739669" CREATED="1611134976553" MODIFIED="1611134980227"/>
</node>
<node TEXT="generatePairs()" FOLDED="true" ID="ID_428338018" CREATED="1611134981556" MODIFIED="1611140157377" BACKGROUND_COLOR="#ffff99">
<node TEXT="Generate pairs" ID="ID_1902646126" CREATED="1611135154252" MODIFIED="1611135164354"/>
</node>
<node TEXT="createSubModels()" ID="ID_1670810133" CREATED="1611134994573" MODIFIED="1611139745204" BACKGROUND_COLOR="#00ffff">
<node TEXT="Generate pairs and sub-model tables" ID="ID_511473553" CREATED="1611135122537" MODIFIED="1611135152272"/>
</node>
<node TEXT="generatePairsAndSubModels" ID="ID_1074577567" CREATED="1611135001941" MODIFIED="1611139745204" BACKGROUND_COLOR="#00ffff">
<node TEXT="Generate pairs and sub-model tables" ID="ID_503700419" CREATED="1611135106405" MODIFIED="1611135107774"/>
<node TEXT="Generate pairs and blended sub-model tables" ID="ID_1431659720" CREATED="1611135108120" MODIFIED="1611135121527"/>
<node TEXT="Generate pairs and two-sided sub-model tables" ID="ID_108921543" CREATED="1611135188213" MODIFIED="1611135192145"/>
</node>
<node TEXT="sumAlphaMoving()" FOLDED="true" ID="ID_60774567" CREATED="1611135200801" MODIFIED="1611140166052" BACKGROUND_COLOR="#ffff99">
<node TEXT="Return the sum of the phase fractions of the moving phases" ID="ID_1934877155" CREATED="1611135239329" MODIFIED="1611135241393"/>
</node>
<node TEXT="setMixtureU()" FOLDED="true" ID="ID_868652333" CREATED="1611135211360" MODIFIED="1611140172924" BACKGROUND_COLOR="#ffff99">
<node TEXT="Re-normalise the velocity of the phases around the specified mixture mean" ID="ID_1653762123" CREATED="1611135243488" MODIFIED="1611135269300"/>
</node>
<node TEXT="setMixturePhi()" ID="ID_3948302" CREATED="1611135222981" MODIFIED="1611140178790" BACKGROUND_COLOR="#ffff99">
<node TEXT="Re-normalise the flux of the phases around the specified mixture mean" ID="ID_807894114" CREATED="1611135273080" MODIFIED="1611135282170"/>
</node>
<node TEXT="functions required for interface compression" ID="ID_153268146" CREATED="1611135313785" MODIFIED="1611135331681">
<node TEXT="nHatfv()" ID="ID_505414119" CREATED="1611135334301" MODIFIED="1611140205613" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="nHatf()" ID="ID_376795922" CREATED="1611135344961" MODIFIED="1611140205613" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="correctContactAngle()" ID="ID_449202687" CREATED="1611135353296" MODIFIED="1611140205613" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="K()" ID="ID_1407656221" CREATED="1611135380801" MODIFIED="1611140205613" BACKGROUND_COLOR="#ffff99"/>
</node>
<node TEXT="functions required by twoPhaseSystem" ID="ID_722098846" CREATED="1611135416708" MODIFIED="1611140517234"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      defined as virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Kd()" ID="ID_1882234470" CREATED="1611135434425" MODIFIED="1611135446229">
<node TEXT="Return the drag coefficient for phase pair" ID="ID_436775323" CREATED="1611135447653" MODIFIED="1611135450693"/>
</node>
<node TEXT="Vm()" ID="ID_1569523092" CREATED="1611135463417" MODIFIED="1611135465419">
<node TEXT="Return the virtual mass coefficient for phase pair" ID="ID_36885446" CREATED="1611135477337" MODIFIED="1611135479997"/>
</node>
<node TEXT="Kdf()" ID="ID_1786521354" CREATED="1611135466905" MODIFIED="1611135470098">
<node TEXT="Return the face drag coefficient for phase pair" ID="ID_1632460493" CREATED="1611135485885" MODIFIED="1611135488361"/>
</node>
</node>
</node>
</node>
<node TEXT="Public" POSITION="right" ID="ID_1529041872" CREATED="1611135502032" MODIFIED="1611135505962">
<edge COLOR="#ff00ff"/>
<node TEXT="declaration" ID="ID_529736966" CREATED="1611135507629" MODIFIED="1611135906950">
<node TEXT="phaseSystem" ID="ID_936367680" CREATED="1611135910229" MODIFIED="1611135924947">
<node TEXT="Runtime type information" ID="ID_1514099818" CREATED="1611136038936" MODIFIED="1611136044565"/>
</node>
<node TEXT="propertiesName" ID="ID_1523158328" CREATED="1611135925536" MODIFIED="1611140126862" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="declareRunTimeSelectionTable" ID="ID_505294684" CREATED="1611135930885" MODIFIED="1611135945270"/>
<node TEXT="phaseSystem()" ID="ID_712420015" CREATED="1611135945685" MODIFIED="1611140227701" BACKGROUND_COLOR="#ffff99">
<node TEXT="constructors" ID="ID_160788601" CREATED="1611135989637" MODIFIED="1611136008084"/>
</node>
<node TEXT="New" FOLDED="true" ID="ID_1248053942" CREATED="1611135970839" MODIFIED="1611140447252" BACKGROUND_COLOR="#ff9900">
<node TEXT="selectors" ID="ID_1234331931" CREATED="1611135995705" MODIFIED="1611136004095"/>
</node>
<node TEXT="~phaseSystem()" ID="ID_1012281728" CREATED="1611135974857" MODIFIED="1611140227701" BACKGROUND_COLOR="#ffff99">
<node TEXT="destructor" ID="ID_101687693" CREATED="1611136018220" MODIFIED="1611136037708"/>
</node>
</node>
<node TEXT="member functions" ID="ID_1426835561" CREATED="1611136052741" MODIFIED="1611136058889">
<node TEXT="access" ID="ID_1442275744" CREATED="1611136060968" MODIFIED="1611136161501"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      defined as inline functions
    </p>
  </body>
</html>

</richcontent>
<node TEXT="" ID="ID_403249965" CREATED="1611136567005" MODIFIED="1611136567005">
<hook NAME="FirstGroupNode"/>
</node>
<node TEXT="mesh()" ID="ID_1141870226" CREATED="1611136064800" MODIFIED="1611139994976" BACKGROUND_COLOR="#99ff99"/>
<node TEXT="phases()" ID="ID_1029130090" CREATED="1611136168540" MODIFIED="1611139994976" BACKGROUND_COLOR="#99ff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions: const and non-const
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="movingPhases()" ID="ID_1324486329" CREATED="1611136195104" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions: const and non-const
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="stationaryPhases()" ID="ID_810310745" CREATED="1611136200600" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions: const and non-const
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="anisothermalPhases()" ID="ID_424722533" CREATED="1611136212933" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions: const and non-const
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="multiComponentPhases()" ID="ID_1756163998" CREATED="1611136222728" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions: const and non-const
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="phasePairs()" ID="ID_1446552678" CREATED="1611136234704" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"/>
<node TEXT="ptherPhase()" ID="ID_1707994226" CREATED="1611136251376" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"/>
<node TEXT="phi()" ID="ID_1660554586" CREATED="1611136260796" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions: const and non-const
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="dpdt()" ID="ID_811239906" CREATED="1611136263477" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      two definitions: const and non-const
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="MRF()" ID="ID_1174922546" CREATED="1611136274428" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"/>
<node TEXT="fvOptions()" ID="ID_1904212085" CREATED="1611136276952" MODIFIED="1611139994977" BACKGROUND_COLOR="#99ff99"/>
<node TEXT="" ID="ID_273533533" CREATED="1611136567002" MODIFIED="1611136567005">
<hook NAME="SummaryNode"/>
<hook NAME="AlwaysUnfoldedNode"/>
<node TEXT="defined in phaseSystemI.H" ID="ID_1409562776" CREATED="1611136567006" MODIFIED="1611136657426">
<font BOLD="false"/>
</node>
</node>
</node>
<node TEXT="sub-model lookup" ID="ID_1910387793" CREATED="1611136686941" MODIFIED="1611136699970">
<node TEXT="foundSubModel()" ID="ID_406384363" CREATED="1611136701741" MODIFIED="1611139624149" BACKGROUND_COLOR="#00ffff"><richcontent TYPE="DETAILS">

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
<node TEXT="Check availability of a sub model for a given phase pair" ID="ID_1721560936" CREATED="1611136797861" MODIFIED="1611136812174"/>
<node TEXT="Check availability of a sub model between two phases" ID="ID_1226228767" CREATED="1611136812641" MODIFIED="1611136816554"/>
</node>
<node TEXT="lookupSubModel()" ID="ID_868515926" CREATED="1611136713910" MODIFIED="1611139624149" BACKGROUND_COLOR="#00ffff"><richcontent TYPE="DETAILS">

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
<node TEXT="Return a sub model between a phase pair" ID="ID_509198654" CREATED="1611136824233" MODIFIED="1611136827678"/>
<node TEXT="Return a sub model between two phases" ID="ID_1372294155" CREATED="1611136827825" MODIFIED="1611136835178"/>
</node>
<node TEXT="foundBlendedSubModel()" ID="ID_1039497673" CREATED="1611136764248" MODIFIED="1611139624149" BACKGROUND_COLOR="#00ffff">
<node TEXT="Check availability of a blended sub model for a given phase pair" ID="ID_360852293" CREATED="1611136846585" MODIFIED="1611136853321"/>
</node>
<node TEXT="lookupBlendedSubModel()" ID="ID_1559028251" CREATED="1611136742518" MODIFIED="1611139624149" BACKGROUND_COLOR="#00ffff">
<node TEXT="Return a blended sub model between a phase pair" ID="ID_633243305" CREATED="1611136854793" MODIFIED="1611136857775"/>
</node>
</node>
<node TEXT="field construction" ID="ID_1704663455" CREATED="1611136927497" MODIFIED="1611136935919">
<node TEXT="fillFields()" ID="ID_550052949" CREATED="1611136937465" MODIFIED="1611139675700" BACKGROUND_COLOR="#00ffff">
<node TEXT="Fill up gaps in a phase-indexed list of fields with zeros" ID="ID_1016521243" CREATED="1611137023905" MODIFIED="1611137026686"/>
<node TEXT="Fill up gaps in a phase-indexed table of fields with zeros" ID="ID_382717094" CREATED="1611137045101" MODIFIED="1611137063642"/>
</node>
</node>
<node TEXT="properties" ID="ID_860397392" CREATED="1611137067513" MODIFIED="1611137071303">
<node TEXT="rho()" FOLDED="true" ID="ID_1432112446" CREATED="1611137072170" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99">
<node TEXT="Return the mixture density" ID="ID_633095318" CREATED="1611137146729" MODIFIED="1611137150043"/>
</node>
<node TEXT="U()" FOLDED="true" ID="ID_1575047611" CREATED="1611137086702" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99">
<node TEXT="Return the mixture velocity" ID="ID_1036632064" CREATED="1611137154938" MODIFIED="1611137168384"/>
</node>
<node TEXT="E()" FOLDED="true" ID="ID_1129079783" CREATED="1611137090433" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99">
<node TEXT="Return the aspect-ratio for a pair" ID="ID_236534992" CREATED="1611137169481" MODIFIED="1611137172561"/>
</node>
<node TEXT="sigma()" FOLDED="true" ID="ID_819717637" CREATED="1611137092201" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99"><richcontent TYPE="DETAILS">

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
<node TEXT="Return the surface tension coefficient for a pair" ID="ID_1552367620" CREATED="1611137186801" MODIFIED="1611137197751"/>
<node TEXT="Return the surface tension coefficient for a pair on a patch" ID="ID_1177302624" CREATED="1611137198345" MODIFIED="1611137207900"/>
</node>
<node TEXT="nearInterface()" FOLDED="true" ID="ID_1447916056" CREATED="1611137105241" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99">
<node TEXT="Indicator of the proximity of the interface, field values are 1 near and 0 away for the interface" ID="ID_463060695" CREATED="1611137209494" MODIFIED="1611137229494"/>
</node>
<node TEXT="dmdtf()" FOLDED="true" ID="ID_522990057" CREATED="1611137116257" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      defined as virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Return the mass transfer rate for an interface" ID="ID_1877560237" CREATED="1611137245590" MODIFIED="1611137258465"/>
</node>
<node TEXT="dmdts()" FOLDED="true" ID="ID_1096658048" CREATED="1611137127634" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      defined as virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Return the mass transfer rates for each phase" ID="ID_439410465" CREATED="1611137259525" MODIFIED="1611137268905"/>
</node>
<node TEXT="incompressible()" FOLDED="true" ID="ID_124457509" CREATED="1611137132938" MODIFIED="1611140264766" BACKGROUND_COLOR="#ffff99">
<node TEXT="Return incompressibility" ID="ID_1704932324" CREATED="1611137270205" MODIFIED="1611137273190"/>
</node>
</node>
<node TEXT="Transfer" ID="ID_230377253" CREATED="1611137286830" MODIFIED="1611137819212"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      defined as virtual
    </p>
    <ul>
      <li>
        f means for face based algorithm
      </li>
      <li>
        s means it's a list
      </li>
    </ul>
  </body>
</html>

</richcontent>
<node TEXT="momentumTransfer()" FOLDED="true" ID="ID_1250333864" CREATED="1611137650254" MODIFIED="1611137668016">
<node TEXT="Return the momentum transfer matrices for the cell-based algorithm" ID="ID_247038376" CREATED="1611137669046" MODIFIED="1611137693287"/>
</node>
<node TEXT="momentumTransferf()" FOLDED="true" ID="ID_1927645106" CREATED="1611137710870" MODIFIED="1611137719423">
<node TEXT="Return the momentum transfer matrices for the face-based algorithm" ID="ID_1134143408" CREATED="1611138185662" MODIFIED="1611138193659"/>
</node>
<node TEXT="AFfs()" FOLDED="true" ID="ID_828176415" CREATED="1611137722618" MODIFIED="1611137726151">
<node TEXT="Return the implicit force coefficients for the face-based algorithm" ID="ID_387218864" CREATED="1611138201750" MODIFIED="1611138239745"/>
</node>
<node TEXT="phiFs()" FOLDED="true" ID="ID_775637807" CREATED="1611137727366" MODIFIED="1611137738695">
<node TEXT="Return the force fluxes for the cell-based algorithm" ID="ID_1115995952" CREATED="1611138240858" MODIFIED="1611138250887"/>
</node>
<node TEXT="phiFfs()" FOLDED="true" ID="ID_1300987994" CREATED="1611137739553" MODIFIED="1611137758524">
<node TEXT="Return the force fluxes for the face-based algorithm" FOLDED="true" ID="ID_1590713756" CREATED="1611138259144" MODIFIED="1611138259144">
<node TEXT="virtual PtrList&lt;surfaceScalarField&gt; phiFfs" ID="ID_20603473" CREATED="1611138259144" MODIFIED="1611138259144"/>
</node>
</node>
<node TEXT="phiKdPhis()" FOLDED="true" ID="ID_1326083237" CREATED="1611137846163" MODIFIED="1611137855315">
<node TEXT="Return the force fluxes for the cell-based algorithm" ID="ID_1292344073" CREATED="1611138272612" MODIFIED="1611138272612"/>
</node>
<node TEXT="phiKdPhifs()" FOLDED="true" ID="ID_1110907067" CREATED="1611137856081" MODIFIED="1611137867391">
<node TEXT="Return the force fluxes for the face-based algorithm" ID="ID_27674799" CREATED="1611138285196" MODIFIED="1611138285196"/>
</node>
<node TEXT="KdUByAs()" FOLDED="true" ID="ID_110024133" CREATED="1611137867894" MODIFIED="1611137884304">
<node TEXT="Return the explicit part of the drag force" ID="ID_1340107870" CREATED="1611138169130" MODIFIED="1611138184457"/>
</node>
<node TEXT="implicitPhasesPressure()" FOLDED="true" ID="ID_285567292" CREATED="1611137891282" MODIFIED="1611140289094" BACKGROUND_COLOR="#ffff99"><richcontent TYPE="DETAILS">

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
<node TEXT="Returns true if the phase pressure is treated implicitly in the phase fraction equation" ID="ID_524076076" CREATED="1611137926478" MODIFIED="1611137949183"/>
<node TEXT="Returns true if the phase pressure is treated implicitly in the phase fraction equation for any phase" ID="ID_1192277971" CREATED="1611137949653" MODIFIED="1611137956407"/>
</node>
<node TEXT="DByAfs()" FOLDED="true" ID="ID_528648945" CREATED="1611137966510" MODIFIED="1611137973171">
<node TEXT="Return the phase diffusivity divided by the momentum central coefficient" ID="ID_877159223" CREATED="1611138145722" MODIFIED="1611138168146"/>
</node>
<node TEXT="partialElimination()" FOLDED="true" ID="ID_1870801425" CREATED="1611137974242" MODIFIED="1611137980492">
<node TEXT="Solve the drag system for the new velocities and fluxes" ID="ID_687701213" CREATED="1611138132967" MODIFIED="1611138138039"/>
</node>
<node TEXT="partialEliminationf()" FOLDED="true" ID="ID_765640149" CREATED="1611137980858" MODIFIED="1611137991051">
<node TEXT="Solve the drag system for the new fluxes" ID="ID_1867648403" CREATED="1611138115406" MODIFIED="1611138131084"/>
</node>
<node TEXT="ddtCorrByAs()" FOLDED="true" ID="ID_1555871790" CREATED="1611137991610" MODIFIED="1611137998815">
<node TEXT="Return the flux corrections for the cell-based algorithm" ID="ID_1560177633" CREATED="1611138098601" MODIFIED="1611138105067"/>
</node>
<node TEXT="heatTransfer()" FOLDED="true" ID="ID_1200775483" CREATED="1611137999171" MODIFIED="1611138011796">
<node TEXT="Return the heat transfer matrices" ID="ID_72269356" CREATED="1611138063738" MODIFIED="1611138068411"/>
</node>
<node TEXT="specieTransfer()" FOLDED="true" ID="ID_47603591" CREATED="1611138012802" MODIFIED="1611138018519">
<node TEXT="Return the specie transfer matrices" ID="ID_1265975038" CREATED="1611138050638" MODIFIED="1611138061841"/>
</node>
<node TEXT="surfaceTension()" ID="ID_1336193551" CREATED="1611138019186" MODIFIED="1611140289093" BACKGROUND_COLOR="#ffff99"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      not virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Return the surface tension force" ID="ID_1762813242" CREATED="1611138035966" MODIFIED="1611138048917"/>
</node>
</node>
<node TEXT="Evolution" ID="ID_1715914849" CREATED="1611137312881" MODIFIED="1611137596863"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      defined as virtual; solve and correct
    </p>
  </body>
</html>

</richcontent>
<node TEXT="solve()" FOLDED="true" ID="ID_923328542" CREATED="1611137460285" MODIFIED="1611140051272" BACKGROUND_COLOR="#ff99ff">
<node TEXT="Solve for the phase fractions" ID="ID_281401809" CREATED="1611137571289" MODIFIED="1611137574158"/>
</node>
<node TEXT="correct()" FOLDED="true" ID="ID_394372308" CREATED="1611137465725" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99">
<node TEXT="Correct the fluid properties other than those listed below" ID="ID_1466200910" CREATED="1611137555398" MODIFIED="1611137558366"/>
</node>
<node TEXT="correctContinuityError()" ID="ID_441915662" CREATED="1611137470073" MODIFIED="1611140364277" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="correctKinematics()" ID="ID_1233396352" CREATED="1611137480990" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="correctThermo()" ID="ID_1106634332" CREATED="1611137490525" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="correctReactions()" ID="ID_220371915" CREATED="1611140374056" MODIFIED="1611140386017" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="correctSpecies()" ID="ID_1646669768" CREATED="1611137496377" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="correctTurbulence()" ID="ID_1775636394" CREATED="1611137522361" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="correctEnergyTransport()" ID="ID_279169362" CREATED="1611137531130" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99"/>
</node>
<node TEXT="IO" ID="ID_1228496289" CREATED="1611137329345" MODIFIED="1611137445307"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      defined as virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="read()" ID="ID_564768683" CREATED="1611137340593" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99">
<node TEXT="Read base phaseProperties dictionary" ID="ID_439216723" CREATED="1611137356905" MODIFIED="1611137360874"/>
</node>
</node>
</node>
</node>
<node TEXT="other definitions" POSITION="right" ID="ID_404951138" CREATED="1611138321935" MODIFIED="1611138332848">
<edge COLOR="#00ffff"/>
<node TEXT="byDt()" ID="ID_2513212" CREATED="1611138333963" MODIFIED="1611140364278" BACKGROUND_COLOR="#ffff99"><richcontent TYPE="DETAILS">

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
<node TEXT="for volScarlarField" ID="ID_1674095355" CREATED="1611138361135" MODIFIED="1611138368427"/>
<node TEXT="for surfaceScalarField" ID="ID_1130952232" CREATED="1611138368898" MODIFIED="1611138378971"/>
</node>
</node>
<node TEXT="#include" FOLDED="true" POSITION="right" ID="ID_964177636" CREATED="1611138395418" MODIFIED="1611138402351">
<edge COLOR="#7c0000"/>
<node TEXT="phaseSystemI.H" ID="ID_1920513622" CREATED="1611138402739" MODIFIED="1611138410288"/>
</node>
<node TEXT="remarks" POSITION="right" ID="ID_1790549476" CREATED="1611139770417" MODIFIED="1611139775360">
<edge COLOR="#00007c"/>
<node TEXT="phaseSystem.C" ID="ID_1771530383" CREATED="1611139776051" MODIFIED="1611139954194" BACKGROUND_COLOR="#ffff99"/>
<node TEXT="phaseSystem.H" ID="ID_519676628" CREATED="1611139831731" MODIFIED="1611139836474"/>
<node TEXT="phaseSystemI.H" ID="ID_63305327" CREATED="1611139791676" MODIFIED="1611139966176" BACKGROUND_COLOR="#99ff99"/>
<node TEXT="phaseSystemNew.C" ID="ID_1812007754" CREATED="1611139838513" MODIFIED="1611140433929" BACKGROUND_COLOR="#ff9900"/>
<node TEXT="phaseSystemSolve.C" ID="ID_1397336956" CREATED="1611139850360" MODIFIED="1611139947425" BACKGROUND_COLOR="#ff99ff"/>
<node TEXT="phaseSystemTemplates.C" ID="ID_444832457" CREATED="1611139800427" MODIFIED="1611139881348" BACKGROUND_COLOR="#00ffff"/>
</node>
</node>
</map>
