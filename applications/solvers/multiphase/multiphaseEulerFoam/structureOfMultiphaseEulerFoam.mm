<map version="freeplane 1.8.0">
<!--To view this file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net -->
<node TEXT="multiphaseEulerFoam" FOLDED="false" ID="ID_919720003" CREATED="1611342930810" MODIFIED="1611342944105" STYLE="oval">
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
<hook NAME="AutomaticEdgeColor" COUNTER="11" RULE="ON_BRANCH_CREATION"/>
<node TEXT="functionObjects" POSITION="right" ID="ID_26080841" CREATED="1611342946940" MODIFIED="1611342991016">
<edge COLOR="#ff0000"/>
<node TEXT="make" ID="ID_1824772266" CREATED="1611343179516" MODIFIED="1611343181200"/>
<node TEXT="phaseForces" ID="ID_1713292859" CREATED="1611343181675" MODIFIED="1611343340526"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      This functionObject calculates and outputs the blended interfacial forces acting on a given phase, i.e. drag, virtual mass, lift, wall-lubrication and turbulent dispersion. Note that it works only in run-time processing mode and in combination with the reactingEulerFoam solvers.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="phaseMap" ID="ID_244770064" CREATED="1611343190743" MODIFIED="1611343378017"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      This functionObject writes the phase-fraction map field alpha.map with incremental value ranges for each phase e.g., with values 0 for water, 1 for air, 2 for oil etc.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="sizeDistribution" ID="ID_1931314155" CREATED="1611343194887" MODIFIED="1611343448165"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      This function object calculates and outputs volume-averaged information about the size distribution of the dispersed phase, such as the number density function or its moments. It is designed to be used exclusively with the population balance modeling functionality of the reactingEulerFoam solvers. It can be applied to a specific cellZone or the entire domain. The function type determines whether the density function and its moments are based on the number of dispersed phase elements in a size group or their total volume.
    </p>
  </body>
</html>

</richcontent>
</node>
</node>
<node TEXT="include" POSITION="right" ID="ID_1427443934" CREATED="1611342991574" MODIFIED="1611342995700">
<edge COLOR="#0000ff"/>
<node TEXT="createRDeltaTf.H" ID="ID_814964704" CREATED="1611343474042" MODIFIED="1611343479813"/>
<node TEXT="setRDeltaTf.H" ID="ID_943664726" CREATED="1611343480227" MODIFIED="1611343492921"/>
</node>
<node TEXT="interfacialCompositionModels" POSITION="right" ID="ID_591441573" CREATED="1611342996038" MODIFIED="1611343005888">
<edge COLOR="#00ff00"/>
<node TEXT="diffusiveMassTransferModels" ID="ID_1351213013" CREATED="1611343509802" MODIFIED="1611343519992">
<node TEXT="diffusiveMassTransferModel" ID="ID_1391783198" CREATED="1611343612979" MODIFIED="1611343623748">
<node TEXT="diffusiveMassTransferModel.C" ID="ID_1148611250" CREATED="1611343777270" MODIFIED="1611343783334"/>
<node TEXT="diffusiveMassTransferModel.H" ID="ID_427253553" CREATED="1611343783883" MODIFIED="1611343787258"/>
<node TEXT="diffusiveMassTransferModelNew.C" ID="ID_1943825764" CREATED="1611343787892" MODIFIED="1611343795285"/>
</node>
<node TEXT="Frossling" ID="ID_701949611" CREATED="1611343624371" MODIFIED="1611343698789"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Frossling correlation for turbulent mass transfer from the surface of a sphere to the surrounding fluid.
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Frossling.C" ID="ID_1904995240" CREATED="1611343739086" MODIFIED="1611343746425"/>
<node TEXT="Frossling.H" ID="ID_343405915" CREATED="1611343747008" MODIFIED="1611343750496"/>
</node>
<node TEXT="sphericalDiffusiveMassTransfer" ID="ID_120577878" CREATED="1611343629503" MODIFIED="1611343718421"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Model which applies an analytical solution for mass transfer from the surface of a sphere to the fluid within the sphere.
    </p>
  </body>
</html>

</richcontent>
<node TEXT="sphericalDiffusiveMassTransfer.C" ID="ID_1325168856" CREATED="1611343805007" MODIFIED="1611343808088"/>
<node TEXT="sphericalDiffusiveMassTransfer.H" ID="ID_851067668" CREATED="1611343808376" MODIFIED="1611343815453"/>
</node>
</node>
<node TEXT="interfaceCompositionModels" ID="ID_1198251783" CREATED="1611343521715" MODIFIED="1611343531860">
<node TEXT="interfaceCompositionModel" ID="ID_213024467" CREATED="1611343840963" MODIFIED="1611343878998"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Generic base class for interface composition models. These models describe the composition in phase 1 of the supplied pair at the interface with phase 2.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="Henry" ID="ID_434898189" CREATED="1611343880828" MODIFIED="1611344059322"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Henry's law for gas solubility in liquid. The concentration of a dissolved species in the liquid is proportional to its partial pressure in the gas. A dimensionless solubility, $k$, is given for each species. This is the ratio of the concentration of the species in the liquid to the corresponding concentration in the gas; i.e., $k =&nbsp;&nbsp;c_{i,liq}/c_{i,gas}$. Mixing in the gas is assumed to be ideal.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="nonRandomTwoLiquid" ID="ID_1939509718" CREATED="1611343924895" MODIFIED="1611344093873"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Non ideal law for the mixing of two species. A separate composition model is given for each species. The composition of a species is equal to the value given by the model, scaled by the species fraction in the bulk of the other phase, and multiplied by the activity coefficient for that species. The gas behaviour is assumed ideal; i.e. the fugacity coefficient is taken as equal to 1.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="Raoult" ID="ID_1113899952" CREATED="1611343935384" MODIFIED="1611344145570"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Raoult's law of ideal mixing. A separate composition model is given for each species. The composition of a species is equal to the value given by the model scaled by the species fraction in the bulk of the other phase.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="saturated" ID="ID_1214325357" CREATED="1611343941962" MODIFIED="1611344167909"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Model which uses a saturation pressure model for a single species to calculate the interface composition.
    </p>
  </body>
</html>

</richcontent>
</node>
</node>
<node TEXT="Make" ID="ID_1537803715" CREATED="1611343532323" MODIFIED="1611343533773"/>
<node TEXT="saturationModels" ID="ID_23027508" CREATED="1611343534815" MODIFIED="1611343540161">
<node TEXT="saturationModel" ID="ID_571056571" CREATED="1611344223516" MODIFIED="1611344229265"/>
<node TEXT="Antoine" ID="ID_1893305579" CREATED="1611344230508" MODIFIED="1611344236513"/>
<node TEXT="AntoineExtended" ID="ID_285647210" CREATED="1611344237468" MODIFIED="1611344245669"/>
<node TEXT="ArdenBuck" ID="ID_1865832464" CREATED="1611344246004" MODIFIED="1611344251309"/>
<node TEXT="constantSaturationConditions" ID="ID_1332745945" CREATED="1611344251791" MODIFIED="1611344260539"/>
<node TEXT="function1" ID="ID_1588301834" CREATED="1611344264052" MODIFIED="1611344266168"/>
<node TEXT="polynomial" ID="ID_1673034674" CREATED="1611344266988" MODIFIED="1611344270113"/>
</node>
<node TEXT="surfaceTensionModels" ID="ID_1243222347" CREATED="1611343542166" MODIFIED="1611343566041">
<node TEXT="surfaceTensionModel" ID="ID_1082059199" CREATED="1611344316837" MODIFIED="1611344321797"/>
<node TEXT="constantSurfaceTensionCoefficient" ID="ID_1397588241" CREATED="1611344324859" MODIFIED="1611344364742"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Constant value surface tension model.
    </p>
  </body>
</html>

</richcontent>
<node TEXT="constantSurfaceTensionCoefficient.C" ID="ID_506410159" CREATED="1611344355674" MODIFIED="1611344373298"/>
<node TEXT="constantSurfaceTensionCoefficient.H" ID="ID_1822917556" CREATED="1611344373654" MODIFIED="1611344378962"/>
</node>
</node>
</node>
<node TEXT="interfacialModels" POSITION="right" ID="ID_1746247967" CREATED="1611343010475" MODIFIED="1611343014688">
<edge COLOR="#ff00ff"/>
<node TEXT="aspectRatioModels" ID="ID_1922384465" CREATED="1611351132217" MODIFIED="1611351139940">
<node TEXT="aspectRatioModel" ID="ID_888568580" CREATED="1611351800334" MODIFIED="1611351815200"/>
<node TEXT="constantAspectRatio" ID="ID_897466769" CREATED="1611351815770" MODIFIED="1611351824204"/>
<node TEXT="TomiyamaAspectRatio" ID="ID_1083118244" CREATED="1611351828285" MODIFIED="1611351837720"/>
<node TEXT="VakhrushevEfremov" ID="ID_730492015" CREATED="1611351838470" MODIFIED="1611351850189"/>
<node TEXT="Wellek" ID="ID_218589646" CREATED="1611351851217" MODIFIED="1611351859296"/>
</node>
<node TEXT="dragModels" ID="ID_555226416" CREATED="1611351142643" MODIFIED="1611351146612">
<node TEXT="dragModel" ID="ID_1360348715" CREATED="1611351861494" MODIFIED="1611351870301"/>
<node TEXT="aerosolDrag" ID="ID_1438990928" CREATED="1611351887610" MODIFIED="1611351910625"/>
<node TEXT="AttouFerschneider" ID="ID_970289238" CREATED="1611351912490" MODIFIED="1611351922371"/>
<node TEXT="Beetstra" ID="ID_508364578" CREATED="1611351923869" MODIFIED="1611351928531"/>
<node TEXT="Ergun" ID="ID_1227250150" CREATED="1611351929701" MODIFIED="1611351933917"/>
<node TEXT="Gibilaro" ID="ID_1167361364" CREATED="1611351934315" MODIFIED="1611351939446"/>
<node TEXT="GidaspowErgunWenYu" ID="ID_1536762313" CREATED="1611351941719" MODIFIED="1611351952180"/>
<node TEXT="GidaspowSchillerNaumann" ID="ID_912859938" CREATED="1611351952482" MODIFIED="1611351961913"/>
<node TEXT="IshiiZuber" ID="ID_871301906" CREATED="1611351969213" MODIFIED="1611351970344"/>
<node TEXT="Lain" ID="ID_1937167760" CREATED="1611351971434" MODIFIED="1611351975492"/>
<node TEXT="SchillerNaumann" ID="ID_1118414996" CREATED="1611351976575" MODIFIED="1611351987857"/>
<node TEXT="segregated" ID="ID_348935783" CREATED="1611351988014" MODIFIED="1611351996921"/>
<node TEXT="SyamlalOBrien" ID="ID_1708672278" CREATED="1611351997811" MODIFIED="1611352008700"/>
<node TEXT="Tenneti" ID="ID_363434786" CREATED="1611352008934" MODIFIED="1611352014689"/>
<node TEXT="timeScaleFilteredDrag" ID="ID_1775722854" CREATED="1611352014986" MODIFIED="1611352023885"/>
<node TEXT="TomiyamaAnalytic" ID="ID_1271929955" CREATED="1611352024058" MODIFIED="1611352032309"/>
<node TEXT="TomiyamaCorrelated" ID="ID_709008412" CREATED="1611352032474" MODIFIED="1611352040340"/>
<node TEXT="TomiyamaKataokaZunSakaguchi" ID="ID_384871524" CREATED="1611352040519" MODIFIED="1611352049992"/>
<node TEXT="WenYu" ID="ID_1432887946" CREATED="1611352052029" MODIFIED="1611352055970"/>
</node>
<node TEXT="heatTransferModels" ID="ID_1789391182" CREATED="1611351147786" MODIFIED="1611351152348">
<node TEXT="heatTransferModel" ID="ID_262585139" CREATED="1611351728691" MODIFIED="1611351734348"/>
<node TEXT="constantNu" ID="ID_1855410543" CREATED="1611351735634" MODIFIED="1611351738967"/>
<node TEXT="Gunn" ID="ID_981764016" CREATED="1611351739273" MODIFIED="1611351742203"/>
<node TEXT="RanzMarshall" ID="ID_1952311370" CREATED="1611351742418" MODIFIED="1611351749557"/>
<node TEXT="sphericalHeatTransfer" ID="ID_1140499536" CREATED="1611351751688" MODIFIED="1611351761119"/>
<node TEXT="timeScaleFilteredHeatTransfer" ID="ID_860792456" CREATED="1611351761797" MODIFIED="1611351773560"/>
</node>
<node TEXT="liftModels" ID="ID_255915873" CREATED="1611351153710" MODIFIED="1611351157489">
<node TEXT="liftModel" ID="ID_455842364" CREATED="1611351649446" MODIFIED="1611351654892"/>
<node TEXT="constantLiftCoefficient" ID="ID_977160227" CREATED="1611351655407" MODIFIED="1611351664473"/>
<node TEXT="LegendreMagnaudet" ID="ID_672575966" CREATED="1611351666062" MODIFIED="1611351683956"/>
<node TEXT="Moraga" ID="ID_1848593010" CREATED="1611351688578" MODIFIED="1611351692776"/>
<node TEXT="noLift" ID="ID_1428889432" CREATED="1611351694590" MODIFIED="1611351697680"/>
<node TEXT="TomiyamaLift" ID="ID_1525806060" CREATED="1611351697970" MODIFIED="1611351705820"/>
<node TEXT="wallDampedLift" ID="ID_178934439" CREATED="1611351706122" MODIFIED="1611351713604"/>
</node>
<node TEXT="Make" ID="ID_1528813338" CREATED="1611351158689" MODIFIED="1611351161651"/>
<node TEXT="phaseTransferModels" ID="ID_1511552280" CREATED="1611351161970" MODIFIED="1611351167749">
<node TEXT="phaseTransferModel" ID="ID_780022971" CREATED="1611351622790" MODIFIED="1611351629801"/>
<node TEXT="deposition" ID="ID_522119637" CREATED="1611351630214" MODIFIED="1611351632161"/>
<node TEXT="reactionDriven" ID="ID_149166398" CREATED="1611351632514" MODIFIED="1611351638269"/>
</node>
<node TEXT="swarmCorrections" ID="ID_612918309" CREATED="1611351171325" MODIFIED="1611351180177">
<node TEXT="swarmCorrection" ID="ID_1855380449" CREATED="1611351573062" MODIFIED="1611351584228"/>
<node TEXT="noSwarm" ID="ID_1022367719" CREATED="1611351584630" MODIFIED="1611351591172"/>
<node TEXT="TomiyamaSwarm" ID="ID_1851665494" CREATED="1611351591518" MODIFIED="1611351603572"/>
</node>
<node TEXT="turbulentDispersionModels" ID="ID_1063897565" CREATED="1611351181522" MODIFIED="1611351192148">
<node TEXT="turbulentDispersionModel" ID="ID_685555639" CREATED="1611351481930" MODIFIED="1611351496439"/>
<node TEXT="Burns" ID="ID_199081885" CREATED="1611351502038" MODIFIED="1611351505232"/>
<node TEXT="constantTurbulentDispersionCoefficient" ID="ID_72158253" CREATED="1611351506222" MODIFIED="1611351517961"/>
<node TEXT="Gosman" ID="ID_455380031" CREATED="1611351522574" MODIFIED="1611351526549"/>
<node TEXT="LopezDeBertodano" ID="ID_1462406937" CREATED="1611351526873" MODIFIED="1611351536224"/>
<node TEXT="noTurbulentDispersion" ID="ID_1762704041" CREATED="1611351537869" MODIFIED="1611351547361"/>
</node>
<node TEXT="virtualMassModels" ID="ID_588149717" CREATED="1611351193661" MODIFIED="1611351205319">
<node TEXT="virtualMassModel" ID="ID_1626528853" CREATED="1611351400425" MODIFIED="1611351411256"/>
<node TEXT="constantVirtualMassCoefficient" ID="ID_537180720" CREATED="1611351374174" MODIFIED="1611351385440"/>
<node TEXT="Lamb" ID="ID_899454194" CREATED="1611351389534" MODIFIED="1611351391741"/>
<node TEXT="noVirtualMass" ID="ID_453780699" CREATED="1611351392082" MODIFIED="1611351396633"/>
</node>
<node TEXT="wallDampingModels" ID="ID_1902385799" CREATED="1611351207530" MODIFIED="1611351212043">
<node TEXT="wallDampingModel" ID="ID_576150090" CREATED="1611351362071" MODIFIED="1611351366913"/>
<node TEXT="cosine" ID="ID_1103839527" CREATED="1611351354790" MODIFIED="1611351357055"/>
<node TEXT="linear" ID="ID_951629113" CREATED="1611351357350" MODIFIED="1611351359043"/>
<node TEXT="sine" ID="ID_1337940846" CREATED="1611351359354" MODIFIED="1611351361704"/>
</node>
<node TEXT="wallDependentModel" ID="ID_1563234934" CREATED="1611351213066" MODIFIED="1611351223941">
<node TEXT="wallDependentModel.C" ID="ID_874155992" CREATED="1611351329430" MODIFIED="1611351332104"/>
<node TEXT="wallDependentModel.H" ID="ID_1826852944" CREATED="1611351332669" MODIFIED="1611351335440"/>
</node>
<node TEXT="wallLubricationModels" ID="ID_657402582" CREATED="1611351226462" MODIFIED="1611351233548">
<node TEXT="wallLubricationModel" ID="ID_1837894988" CREATED="1611351255082" MODIFIED="1611351262963"/>
<node TEXT="Antal" ID="ID_1540340845" CREATED="1611351263822" MODIFIED="1611351267005"/>
<node TEXT="Frank" ID="ID_111144315" CREATED="1611351267350" MODIFIED="1611351268989"/>
<node TEXT="noWallLubrication" ID="ID_1476809992" CREATED="1611351270064" MODIFIED="1611351280540"/>
<node TEXT="TomiyamaWallLubrication" ID="ID_45787378" CREATED="1611351282082" MODIFIED="1611351298409"/>
</node>
</node>
<node TEXT="multiphaseCompressibleMomentumTransportModels" POSITION="right" ID="ID_191255006" CREATED="1611343021714" MODIFIED="1611343052935">
<edge COLOR="#00ffff"/>
<node TEXT="derivedPatchFields" ID="ID_593399612" CREATED="1611350265873" MODIFIED="1611350285972">
<node TEXT="alphatFixedDmdtfWallBoilingWallFunction" ID="ID_453698046" CREATED="1611350611310" MODIFIED="1611350613411">
<node TEXT="alphatFixedDmdtfWallBoilingWallFunction.C" ID="ID_1692397347" CREATED="1611350614461" MODIFIED="1611350624505"/>
<node TEXT="alphatFixedDmdtfWallBoilingWallFunction.H" ID="ID_486909930" CREATED="1611350624817" MODIFIED="1611350627887"/>
</node>
<node TEXT="alphatPhaseChangeWallFunction" ID="ID_1550133602" CREATED="1611350640107" MODIFIED="1611350642083">
<node TEXT="alphatPhaseChangeWallFunction.C" ID="ID_491194350" CREATED="1611350643124" MODIFIED="1611350648323"/>
<node TEXT="alphatPhaseChangeWallFunction.H" ID="ID_1122548561" CREATED="1611350648949" MODIFIED="1611350652129"/>
</node>
<node TEXT="alphatPhaseJayatillekeWallFunction" ID="ID_1902608932" CREATED="1611350662565" MODIFIED="1611350664135">
<node TEXT="alphatPhaseJayatillekeWallFunction.C" ID="ID_1398764439" CREATED="1611350664877" MODIFIED="1611350677745"/>
<node TEXT="alphatPhaseJayatillekeWallFunction.H" ID="ID_152826328" CREATED="1611350671592" MODIFIED="1611350675022"/>
</node>
<node TEXT="alphatWallBoilingWallFunction" ID="ID_742875864" CREATED="1611350687465" MODIFIED="1611350689132">
<node TEXT="alphatWallBoilingWallFunction.C" ID="ID_1432955426" CREATED="1611350689809" MODIFIED="1611350692828"/>
<node TEXT="alphatWallBoilingWallFunction.H" ID="ID_966159729" CREATED="1611350693197" MODIFIED="1611350699791"/>
</node>
<node TEXT="copiedFixedValue" ID="ID_1626041439" CREATED="1611350709149" MODIFIED="1611350710143">
<node TEXT="copiedFixedValue.C" ID="ID_1818960211" CREATED="1611350711121" MODIFIED="1611350713839"/>
<node TEXT="copiedFixedValue.H" ID="ID_989113720" CREATED="1611350714150" MODIFIED="1611350716919"/>
</node>
<node TEXT="fixedMultiPhaseHeatFlux" ID="ID_830821379" CREATED="1611350726130" MODIFIED="1611350727630">
<node TEXT="fixedMultiPhaseHeatFlux.C" ID="ID_394897924" CREATED="1611350728141" MODIFIED="1611350730483"/>
<node TEXT="fixedMultiPhaseHeatFlux.H" ID="ID_1090951555" CREATED="1611350730789" MODIFIED="1611350732779"/>
</node>
<node TEXT="JohnsonJacksonParticleSlip" ID="ID_702513059" CREATED="1611350734414" MODIFIED="1611350744839">
<node TEXT="JohnsonJacksonParticleSlip.C" ID="ID_627764732" CREATED="1611350745757" MODIFIED="1611350749117"/>
<node TEXT="JohnsonJacksonParticleSlip.H" ID="ID_782494397" CREATED="1611350749429" MODIFIED="1611350752431"/>
</node>
<node TEXT="JohnsonJacksonParticleTheta" ID="ID_1658142594" CREATED="1611350763881" MODIFIED="1611350765059">
<node TEXT="JohnsonJacksonParticleTheta.C" ID="ID_64418402" CREATED="1611350765445" MODIFIED="1611350767491"/>
<node TEXT="JohnsonJacksonParticleTheta.H" ID="ID_1562890781" CREATED="1611350767793" MODIFIED="1611350770228"/>
</node>
<node TEXT="wallBoilingSubModels" ID="ID_1153599724" CREATED="1611350783398" MODIFIED="1611350784967">
<node TEXT="departureDiameterModels" ID="ID_797759711" CREATED="1611350786053" MODIFIED="1611350808811">
<node TEXT="departureDiameterModel" ID="ID_1538840168" CREATED="1611350955553" MODIFIED="1611350957688"/>
<node TEXT="KocamustafaogullariIshiiDepartureDiameter" ID="ID_1567630573" CREATED="1611350957934" MODIFIED="1611350971446"/>
<node TEXT="TolubinskiKostanchuk" ID="ID_1373897094" CREATED="1611350971753" MODIFIED="1611350981391"/>
</node>
<node TEXT="departureFrequencyModels" ID="ID_567458290" CREATED="1611350816442" MODIFIED="1611350825675">
<node TEXT="Cole" ID="ID_1459322505" CREATED="1611350987350" MODIFIED="1611350988963"/>
<node TEXT="departureFrequencyModel" ID="ID_1455364484" CREATED="1611350989250" MODIFIED="1611351003004"/>
<node TEXT="KocamustafaogullariIshiiDepartureFrequency" ID="ID_1678553532" CREATED="1611351003482" MODIFIED="1611351005272"/>
</node>
<node TEXT="nucleationSiteModels" ID="ID_475448131" CREATED="1611350826917" MODIFIED="1611350847569">
<node TEXT="KocamustafaogullariIshiiNucleationSite" ID="ID_1708771719" CREATED="1611351018706" MODIFIED="1611351020640"/>
<node TEXT="LemmertChawla" ID="ID_304518785" CREATED="1611351020886" MODIFIED="1611351028412"/>
<node TEXT="nucleationSiteModel" ID="ID_1138609599" CREATED="1611351028562" MODIFIED="1611351038466"/>
</node>
<node TEXT="partitioningModels" ID="ID_949644684" CREATED="1611350848006" MODIFIED="1611350850503">
<node TEXT="cosine" ID="ID_1669776564" CREATED="1611351054885" MODIFIED="1611351060772"/>
<node TEXT="Lavieville" ID="ID_1309804984" CREATED="1611351061165" MODIFIED="1611351070856"/>
<node TEXT="linear" ID="ID_1828782740" CREATED="1611351071906" MODIFIED="1611351074720"/>
<node TEXT="partitioningModel" ID="ID_573827546" CREATED="1611351075906" MODIFIED="1611351084939"/>
<node TEXT="phaseFraction" ID="ID_1045357418" CREATED="1611351087782" MODIFIED="1611351094924"/>
</node>
</node>
</node>
<node TEXT="kineticTheoryModels" ID="ID_1270681089" CREATED="1611350286453" MODIFIED="1611350302443">
<node TEXT="conductivityModel" ID="ID_1503018319" CREATED="1611350433601" MODIFIED="1611350447087"/>
<node TEXT="frictionalStressModel" ID="ID_6615139" CREATED="1611350454668" MODIFIED="1611350470871"/>
<node TEXT="granularPressureModel" ID="ID_1120638865" CREATED="1611350471217" MODIFIED="1611350483572"/>
<node TEXT="kineticTheoryModel" ID="ID_458880170" CREATED="1611350488489" MODIFIED="1611350496900"/>
<node TEXT="radialModel" ID="ID_1339680671" CREATED="1611350497277" MODIFIED="1611350503635"/>
<node TEXT="viscosityModel" ID="ID_623669797" CREATED="1611350505035" MODIFIED="1611350512018"/>
</node>
<node TEXT="Make" ID="ID_1667009981" CREATED="1611350303386" MODIFIED="1611350309272"/>
<node TEXT="phasePressureModel" FOLDED="true" ID="ID_834643236" CREATED="1611350310261" MODIFIED="1611350422044"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Particle-particle phase-pressure RAS model
    </p>
  </body>
</html>

</richcontent>
<node TEXT="phasePressureModel.C" ID="ID_574978509" CREATED="1611350310261" MODIFIED="1611350398294"/>
<node TEXT="phasePressureModel.H" ID="ID_1654073359" CREATED="1611350399504" MODIFIED="1611350403795"/>
</node>
<node TEXT="multiphaseCompressibleMomentumTransportModels.C" ID="ID_590734385" CREATED="1611350321385" MODIFIED="1611350337107"/>
</node>
<node TEXT="multiphaseEulerFoam" POSITION="right" ID="ID_355277845" CREATED="1611343053368" MODIFIED="1611343074913">
<edge COLOR="#7c0000"/>
<node TEXT="Make" ID="ID_1533155572" CREATED="1611346129369" MODIFIED="1611346132627"/>
<node TEXT="multiphaseSystems" ID="ID_1469635566" CREATED="1611346133865" MODIFIED="1611346139631">
<node TEXT="Make" ID="ID_1287747960" CREATED="1611346265279" MODIFIED="1611346268439"/>
<node TEXT="multiphaseSystems.C" ID="ID_346293187" CREATED="1611346268758" MODIFIED="1611346294087"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      typdef phase systems
    </p>
  </body>
</html>

</richcontent>
</node>
</node>
<node TEXT="pU" ID="ID_1679000118" CREATED="1611346141697" MODIFIED="1611346144475">
<node TEXT="pEqn.H" ID="ID_442259848" CREATED="1611346230169" MODIFIED="1611346234215"/>
<node TEXT="UEqns.H" ID="ID_534941714" CREATED="1611346235380" MODIFIED="1611346240990"/>
</node>
<node TEXT="pUf" ID="ID_1333268372" CREATED="1611346144961" MODIFIED="1611346147891">
<node TEXT="pEqn.H" ID="ID_502465653" CREATED="1611346243104" MODIFIED="1611346245931"/>
<node TEXT="UEqns.H" ID="ID_1476041552" CREATED="1611346246329" MODIFIED="1611346251122"/>
</node>
<node TEXT="Allwclean" ID="ID_919128733" CREATED="1611346148965" MODIFIED="1611346153667"/>
<node TEXT="Allwmake" ID="ID_1958754096" CREATED="1611346154013" MODIFIED="1611346156459"/>
<node TEXT="CourantNo.H" ID="ID_1264740177" CREATED="1611346157496" MODIFIED="1611346162603"/>
<node TEXT="createFieldRefs.H" ID="ID_751373251" CREATED="1611346163042" MODIFIED="1611346172815"/>
<node TEXT="createFields.H" ID="ID_1711094555" CREATED="1611346173661" MODIFIED="1611346184611"/>
<node TEXT="EEqns.H" ID="ID_310509915" CREATED="1611346187061" MODIFIED="1611346191359"/>
<node TEXT="multiphaseEulerFoam.C" ID="ID_500526852" CREATED="1611346192718" MODIFIED="1611346201912"/>
<node TEXT="setRDeltaT.H" ID="ID_1613852027" CREATED="1611346207784" MODIFIED="1611346216191"/>
<node TEXT="YEqns.H" ID="ID_1935403341" CREATED="1611346217661" MODIFIED="1611346221455"/>
</node>
<node TEXT="multiphaseThermophysicalTransportModels" POSITION="right" ID="ID_143273305" CREATED="1611343075178" MODIFIED="1611343096964">
<edge COLOR="#00007c"/>
<node TEXT="Make" ID="ID_1716308322" CREATED="1611346109685" MODIFIED="1611346112595"/>
<node TEXT="multiphaseThermophysicalTransportModels.C" ID="ID_1877598682" CREATED="1611346080196" MODIFIED="1611346084798"/>
<node TEXT="rhoReactionMultiphaseThermophysicalTransportModels.C" ID="ID_348212764" CREATED="1611346100301" MODIFIED="1611346104663"/>
</node>
<node TEXT="phaseSystems" POSITION="right" ID="ID_1397540131" CREATED="1611343097274" MODIFIED="1611343108364">
<edge COLOR="#007c00"/>
<node TEXT="alphaContactAngle" ID="ID_692330220" CREATED="1611344431200" MODIFIED="1611344520862"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Contact-angle boundary condition for multi-phase interface-capturing simulations.&nbsp;&nbsp;Used in conjunction with phaseSystem.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="BlendedInterfacialModel" ID="ID_466453374" CREATED="1611344439780" MODIFIED="1611344447066">
<node TEXT="blendingMethods" ID="ID_61948590" CREATED="1611344539140" MODIFIED="1611344546222">
<node TEXT="blendingMethod" ID="ID_142441308" CREATED="1611344587020" MODIFIED="1611344596120"/>
<node TEXT="hyperbolic" ID="ID_869143724" CREATED="1611344597780" MODIFIED="1611344602078"/>
<node TEXT="linear" ID="ID_1195285502" CREATED="1611344602663" MODIFIED="1611344604798"/>
<node TEXT="noBlending" ID="ID_703370052" CREATED="1611344605748" MODIFIED="1611344612470"/>
</node>
<node TEXT="BlendedInterfacialModel.C" ID="ID_1286550280" CREATED="1611344549073" MODIFIED="1611344563237"/>
<node TEXT="BlendedInterfacialModel.H" ID="ID_859739567" CREATED="1611344564532" MODIFIED="1611344571143"/>
</node>
<node TEXT="diameterModels" ID="ID_1694851591" CREATED="1611344448385" MODIFIED="1611344451925">
<node TEXT="diameterModel" ID="ID_1100115800" CREATED="1611344662120" MODIFIED="1611344740914"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Abstract base-class for dispersed-phase particle diameter models.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="constantDiameter" ID="ID_1928728533" CREATED="1611344671520" MODIFIED="1611344770169"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Constant dispersed-phase particle diameter model.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="IATE" ID="ID_1895048466" CREATED="1611344678840" MODIFIED="1611344777510"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      IATE (Interfacial Area Transport Equation) bubble diameter model.
    </p>
  </body>
</html>

</richcontent>
<node TEXT="IATEsources" ID="ID_174202464" CREATED="1611344780187" MODIFIED="1611344790114"/>
<node TEXT="IATE.C" ID="ID_1891162915" CREATED="1611344790528" MODIFIED="1611344794471"/>
<node TEXT="IATE.H" ID="ID_694077320" CREATED="1611344794852" MODIFIED="1611344798493"/>
</node>
<node TEXT="isothermalDiameter" ID="ID_485978694" CREATED="1611344681408" MODIFIED="1611344846788"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Isothermal dispersed-phase particle diameter model.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="linearTsubDiameter" ID="ID_271064027" CREATED="1611344695004" MODIFIED="1611344859054"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Vapour bubble diameter model for modelling of condensation of vapour bubbles. Calculates bubble diameter as a function of liquid phase subcooling.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="sphericalDiameter" ID="ID_1976392181" CREATED="1611344710216" MODIFIED="1611344897941"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Base class for models which represent spherical diameter models, providing a common implementation of surface area per unit volume
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="velocityGroup" ID="ID_1237919516" CREATED="1611344719084" MODIFIED="1611344955070"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      This diameterModel is intended for use with a populationBalanceModel in order to simulate polydispersed bubbly or particulate flows. It can hold any number of sizeGroups from which the Sauter mean diameter is calculated. It can also be used as a diameterModel without a populationBalance and would then behave like a constantDiameter model. In this case, some arbitrary name must be entered for the populationBalance keyword.
    </p>
  </body>
</html>

</richcontent>
<node TEXT="sizeGroup" ID="ID_680100964" CREATED="1611344901184" MODIFIED="1611344904742"/>
<node TEXT="velocityGroup.C" ID="ID_520423924" CREATED="1611344907068" MODIFIED="1611344912118"/>
<node TEXT="velocityGroup.H" ID="ID_841119549" CREATED="1611344980496" MODIFIED="1611344984761"/>
<node TEXT="velocityGroupI.H" ID="ID_26979537" CREATED="1611344985281" MODIFIED="1611344995638"/>
</node>
</node>
<node TEXT="Make" ID="ID_36795805" CREATED="1611344456514" MODIFIED="1611344457598"/>
<node TEXT="phaseModel" ID="ID_365125857" CREATED="1611344458506" MODIFIED="1611344461755">
<node TEXT="phaseModel" ID="ID_1926269986" CREATED="1611345180638" MODIFIED="1611345185847"/>
<node TEXT="AnisothermalPhaseModel" ID="ID_1759946081" CREATED="1611345129377" MODIFIED="1611345141535"/>
<node TEXT="InertPhaseModel" ID="ID_961490484" CREATED="1611345142540" MODIFIED="1611345148067"/>
<node TEXT="IsothermalPhaseModel" ID="ID_1107529327" CREATED="1611345152172" MODIFIED="1611345160198"/>
<node TEXT="MovingPhaseModel" ID="ID_1407775655" CREATED="1611345160581" MODIFIED="1611345167595"/>
<node TEXT="MultiComponentPhaseModel" ID="ID_983874362" CREATED="1611345168356" MODIFIED="1611345177003"/>
<node TEXT="PurePhaseModel" ID="ID_1360505690" CREATED="1611345186965" MODIFIED="1611345194635"/>
<node TEXT="ReactingPhaseModel" ID="ID_951222427" CREATED="1611345196065" MODIFIED="1611345201926"/>
<node TEXT="StationaryPhaseModel" ID="ID_1508807484" CREATED="1611345202885" MODIFIED="1611345212134"/>
<node TEXT="ThermoPhaseModel" ID="ID_540685803" CREATED="1611345218069" MODIFIED="1611345224502"/>
</node>
<node TEXT="phasePair" ID="ID_1043750722" CREATED="1611344462867" MODIFIED="1611344465512">
<node TEXT="orderedPhasePair" ID="ID_1418594899" CREATED="1611345242549" MODIFIED="1611345248959"/>
<node TEXT="phasePair" ID="ID_1902993220" CREATED="1611345249233" MODIFIED="1611345256880"/>
<node TEXT="phasePairKey" ID="ID_1917889152" CREATED="1611345301153" MODIFIED="1611345325718"/>
</node>
<node TEXT="phaseSystem" ID="ID_1848557646" CREATED="1611344465903" MODIFIED="1611344472003">
<node TEXT="phaseSystem.C" ID="ID_275573378" CREATED="1611345643438" MODIFIED="1611345647780"/>
<node TEXT="phaseSystem.H" ID="ID_139937214" CREATED="1611345648297" MODIFIED="1611345650778"/>
<node TEXT="phaseSystemI.H" ID="ID_1443894876" CREATED="1611345651273" MODIFIED="1611345656027"/>
<node TEXT="phaseSystemNew.C" ID="ID_1574137696" CREATED="1611345656445" MODIFIED="1611345661215"/>
<node TEXT="phaseSystemSolve.C" ID="ID_1836755468" CREATED="1611345663065" MODIFIED="1611345667523"/>
<node TEXT="phaseSystemTemplates.C" ID="ID_1969741391" CREATED="1611345667941" MODIFIED="1611345677059"/>
</node>
<node TEXT="phaseSystems" ID="ID_136434564" CREATED="1611344472448" MODIFIED="1611344478562">
<node TEXT="InterfaceCompositionPhaseChangePhaseSystem" ID="ID_1478220673" CREATED="1611345697025" MODIFIED="1611345925015"/>
<node TEXT="MomentumTransferPhaseSystem" ID="ID_954694069" CREATED="1611345925562" MODIFIED="1611345937718"/>
<node TEXT="OneResistanceHeatTransferPhaseSystem" ID="ID_48346229" CREATED="1611345939049" MODIFIED="1611345955995"/>
<node TEXT="PhaseTransferPhaseSystem" ID="ID_1589608933" CREATED="1611345956260" MODIFIED="1611345967694"/>
<node TEXT="PopulationBalancePhaseSystem" ID="ID_1932344099" CREATED="1611345968468" MODIFIED="1611345977946"/>
<node TEXT="ThermalPhaseChangePhaseSystem" ID="ID_367529229" CREATED="1611345978184" MODIFIED="1611345996694"/>
<node TEXT="TwoResistanceHeatTransferPhaseSystem" ID="ID_231964067" CREATED="1611345997163" MODIFIED="1611345998950"/>
</node>
<node TEXT="populationBalanceModel" ID="ID_999455388" CREATED="1611344478892" MODIFIED="1611344487466">
<node TEXT="populationBalanceModel" ID="ID_1028689397" CREATED="1611345375925" MODIFIED="1611345456581"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      Class that solves the univariate population balance equation by means of a class method (also called sectional or discrete method). The internal coordinate is set to the particle volume, so the equation is based on a transport equation of the volume-based number density function. The discretization is done using the fixed pivot technique of Kumar and Ramkrishna (1996). The source terms are written in a way that particle number and mass are preserved. Coalescence (aggregation), breakup, drift (growth and surface loss) as well as nucleation are supported. For the discrete breakup term two recipes are available, depending on the model choice. For models which state a total breakup rate and a separate daughter size distribution function, the formulation of Kumar and Ramkrishna (1996) is applied which is applicable for binary and multiple breakup events. The second formulation is given by Liao et al. (2018). It is useful for binary breakup models which give the breakup rate between a sizeGroup pair directly, without an explicit expression for the daughter size distribution. The drift term is implemented using a finite difference upwind scheme. Although it is diffusive, it ensures a stable and number-conservative solution.
    </p>
  </body>
</html>

</richcontent>
</node>
<node TEXT="binaryBreakupModels" ID="ID_1496833591" CREATED="1611345523554" MODIFIED="1611345544280"/>
<node TEXT="breakupModels" ID="ID_2329277" CREATED="1611345465329" MODIFIED="1611345558494"/>
<node TEXT="coalescenceModels" ID="ID_794344495" CREATED="1611345558904" MODIFIED="1611345564794"/>
<node TEXT="daughterSizeDistributionModels" ID="ID_1895461762" CREATED="1611345566008" MODIFIED="1611345580534"/>
<node TEXT="driftModels" ID="ID_1864540441" CREATED="1611345585852" MODIFIED="1611345592723"/>
<node TEXT="nucleationModels" ID="ID_1765273842" CREATED="1611345593881" MODIFIED="1611345600146"/>
</node>
</node>
<node TEXT="Allwclean" POSITION="right" ID="ID_453049896" CREATED="1611343110586" MODIFIED="1611343116528">
<edge COLOR="#7c007c"/>
</node>
<node TEXT="Allwmake" POSITION="right" ID="ID_1492595504" CREATED="1611343117339" MODIFIED="1611343120175">
<edge COLOR="#007c7c"/>
</node>
<node TEXT="debug.txt" POSITION="right" ID="ID_28909695" CREATED="1611343121386" MODIFIED="1611343177879">
<edge COLOR="#7c7c00"/>
</node>
</node>
</map>
