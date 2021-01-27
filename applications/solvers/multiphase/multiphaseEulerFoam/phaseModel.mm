<map version="freeplane 1.8.0">
<!--To view this file, download free mind mapping software Freeplane from http://freeplane.sourceforge.net -->
<node TEXT="phaseModel" FOLDED="false" ID="ID_1078109624" CREATED="1611781851876" MODIFIED="1611781882087" STYLE="oval">
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
<node TEXT="phaseModel.H" POSITION="right" ID="ID_1080589005" CREATED="1611781883259" MODIFIED="1611781887244">
<edge COLOR="#ff0000"/>
<node TEXT="private data" ID="ID_1103304700" CREATED="1611781893604" MODIFIED="1611781978776">
<node TEXT="fluid_" ID="ID_1423232192" CREATED="1611782066187" MODIFIED="1611782071308">
<node TEXT="Reference to the phaseSystem to which this phase belongs" ID="ID_722403937" CREATED="1611782167569" MODIFIED="1611782167569"/>
</node>
<node TEXT="name_" ID="ID_406455817" CREATED="1611782071592" MODIFIED="1611782074988">
<node TEXT="Name of phase" ID="ID_1876335287" CREATED="1611782162853" MODIFIED="1611782162853"/>
</node>
<node TEXT="index_" ID="ID_1717731267" CREATED="1611782075344" MODIFIED="1611782078340">
<node TEXT="Index of phase" ID="ID_1153729630" CREATED="1611782159594" MODIFIED="1611782159594"/>
</node>
<node TEXT="residualAlpha_" ID="ID_653437841" CREATED="1611782078711" MODIFIED="1611782087461">
<node TEXT="Return the residual phase-fraction for given phase, used to stabilize the phase momentum as the phase-fraction -&gt; 0" ID="ID_434545633" CREATED="1611782127720" MODIFIED="1611782147214"/>
</node>
<node TEXT="alphaMax_" ID="ID_1157007010" CREATED="1611782087860" MODIFIED="1611782091081">
<node TEXT="Optional maximum phase-fraction (e.g. packing limit)" ID="ID_504823735" CREATED="1611782120256" MODIFIED="1611782120256"/>
</node>
<node TEXT="diameterModel_" ID="ID_693647628" CREATED="1611782092499" MODIFIED="1611782101829">
<node TEXT="Diameter model" ID="ID_1853702952" CREATED="1611782110923" MODIFIED="1611782113329"/>
</node>
</node>
<node TEXT="public" ID="ID_1950744957" CREATED="1611781979168" MODIFIED="1611781983249">
<node TEXT="constructors" ID="ID_1896190843" CREATED="1611782180303" MODIFIED="1611784657665" BACKGROUND_COLOR="#ccffff">
<node TEXT="phaseModel()" ID="ID_861472065" CREATED="1611782200840" MODIFIED="1611782215894"/>
<node TEXT="clone()" ID="ID_1333792793" CREATED="1611782204920" MODIFIED="1611782209050"/>
</node>
<node TEXT="selectors" ID="ID_255034811" CREATED="1611782185856" MODIFIED="1611784657665" BACKGROUND_COLOR="#ccffff">
<node TEXT="New()" ID="ID_281045526" CREATED="1611782225556" MODIFIED="1611782229644"/>
<node TEXT="iNew" ID="ID_871598255" CREATED="1611782232284" MODIFIED="1611782237193">
<node TEXT="Return a pointer to a new phase created on freestore from Istream" ID="ID_438652831" CREATED="1611782244868" MODIFIED="1611782252263"/>
</node>
</node>
<node TEXT="destructors" ID="ID_1883449909" CREATED="1611782189535" MODIFIED="1611784657664" BACKGROUND_COLOR="#ccffff">
<node TEXT="~phaseModel()" ID="ID_708591632" CREATED="1611782269369" MODIFIED="1611782269369"/>
</node>
<node TEXT="member functions" ID="ID_604920680" CREATED="1611782274472" MODIFIED="1611782278745">
<node TEXT="name()" ID="ID_129384413" CREATED="1611782284511" MODIFIED="1611784646807" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return the name of this phase" ID="ID_1271540625" CREATED="1611782346237" MODIFIED="1611782346237"/>
</node>
<node TEXT="keyword()" ID="ID_1683691862" CREATED="1611782287175" MODIFIED="1611784646807" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return the name of the phase for use as the keyword in PtrDictionary" ID="ID_1827654005" CREATED="1611782350513" MODIFIED="1611782350513"/>
</node>
<node TEXT="index()" ID="ID_812322530" CREATED="1611782290244" MODIFIED="1611784646808" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return the index of the phase" ID="ID_1122509749" CREATED="1611782354345" MODIFIED="1611782354345"/>
</node>
<node TEXT="fluid()" ID="ID_670446757" CREATED="1611782294903" MODIFIED="1611784646808" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return the system to which this phase belongs" ID="ID_634578027" CREATED="1611782358518" MODIFIED="1611782358518"/>
</node>
<node TEXT="residualAlpha()" ID="ID_560475430" CREATED="1611782301283" MODIFIED="1611784646808" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return the residual phase-fraction for given phase, used to stabilize the phase momentum as the phase-fraction -&gt; 0" ID="ID_1005064568" CREATED="1611782362878" MODIFIED="1611782376022"/>
</node>
<node TEXT="alphaMax()" ID="ID_1568973747" CREATED="1611782312519" MODIFIED="1611784646808" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return the maximum phase-fraction (e.g. packing limit)" ID="ID_508870069" CREATED="1611782380869" MODIFIED="1611782380869"/>
</node>
<node TEXT="d()" ID="ID_305749214" CREATED="1611782315784" MODIFIED="1611784646808" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return the Sauter-mean diameter" ID="ID_406001339" CREATED="1611782385905" MODIFIED="1611782385905"/>
</node>
<node TEXT="dPtr()" FOLDED="true" ID="ID_390993866" CREATED="1611782320192" MODIFIED="1611784646807" BACKGROUND_COLOR="#ccffff">
<node TEXT="Return const-reference to diameterModel of the phase" ID="ID_749945942" CREATED="1611782392658" MODIFIED="1611782392658"/>
</node>
<node TEXT="correct()" FOLDED="true" ID="ID_1707312341" CREATED="1611782324416" MODIFIED="1611784635495" BACKGROUND_COLOR="#ccffff"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the phase properties" ID="ID_285119998" CREATED="1611782399830" MODIFIED="1611782399830"/>
</node>
<node TEXT="correctContinuityError()" ID="ID_1623134047" CREATED="1611782412395" MODIFIED="1611782583098"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the continuity error" ID="ID_196502339" CREATED="1611782510157" MODIFIED="1611782510157"/>
</node>
<node TEXT="correctKinematics()" ID="ID_1584522839" CREATED="1611782427535" MODIFIED="1611782606790"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the kinematics" ID="ID_1410869801" CREATED="1611782514522" MODIFIED="1611782514522"/>
</node>
<node TEXT="correctThermo()" ID="ID_1465163809" CREATED="1611782437364" MODIFIED="1611782609351"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the thermodynamics" ID="ID_1664893441" CREATED="1611782522870" MODIFIED="1611782522870"/>
</node>
<node TEXT="correctReactions()" ID="ID_1541906251" CREATED="1611782443504" MODIFIED="1611782611995"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the reactions" ID="ID_776043159" CREATED="1611782526390" MODIFIED="1611782526390"/>
</node>
<node TEXT="correctSpecies()" ID="ID_1955891849" CREATED="1611782450140" MODIFIED="1611782614191"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the species concentrations" ID="ID_1862433854" CREATED="1611782529937" MODIFIED="1611782529937"/>
</node>
<node TEXT="correctTurbulence()" ID="ID_1124310278" CREATED="1611782461044" MODIFIED="1611782617132"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the turbulence" ID="ID_1583663963" CREATED="1611782533386" MODIFIED="1611782533386"/>
</node>
<node TEXT="correctEnergyTransport()" ID="ID_1520554080" CREATED="1611782468715" MODIFIED="1611782619371"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Correct the energy transport" ID="ID_1969649077" CREATED="1611782536966" MODIFIED="1611782536966"/>
</node>
<node TEXT="correctInflowOutflow()" FOLDED="true" ID="ID_508783473" CREATED="1611782488352" MODIFIED="1611784611621" BACKGROUND_COLOR="#ccffff">
<node TEXT="Ensure that the flux at inflow/outflow BCs is preserved" ID="ID_409589251" CREATED="1611782542586" MODIFIED="1611782542586"/>
</node>
<node TEXT="read()" FOLDED="true" ID="ID_1138910354" CREATED="1611782496832" MODIFIED="1611784625310" BACKGROUND_COLOR="#ccffff"><richcontent TYPE="DETAILS">

<html>
  <head>
    
  </head>
  <body>
    <p>
      virtual
    </p>
  </body>
</html>

</richcontent>
<node TEXT="Read phase properties dictionary" ID="ID_400142216" CREATED="1611782549026" MODIFIED="1611782549026"/>
</node>
<node TEXT="Density variation and compressibility" ID="ID_277354207" CREATED="1611782634820" MODIFIED="1611782637641">
<node TEXT="incompressible()" ID="ID_261228446" CREATED="1611782718488" MODIFIED="1611782725654">
<node TEXT="Return true if the phase is incompressible otherwise false" ID="ID_987702005" CREATED="1611782791366" MODIFIED="1611782791366"/>
</node>
<node TEXT="isochoric()" ID="ID_430205115" CREATED="1611782726056" MODIFIED="1611782730817">
<node TEXT="Return true if the phase is constant density otherwise false" ID="ID_458239789" CREATED="1611782796838" MODIFIED="1611782796838"/>
</node>
<node TEXT="divU()" ID="ID_839480967" CREATED="1611782731159" MODIFIED="1611782735573">
<node TEXT="Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))" ID="ID_965040130" CREATED="1611782801465" MODIFIED="1611782801465"/>
</node>
<node TEXT="divU(divU)" ID="ID_1823776794" CREATED="1611782736004" MODIFIED="1611782781523">
<node TEXT="Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))" ID="ID_241278719" CREATED="1611782806469" MODIFIED="1611782806469"/>
</node>
</node>
<node TEXT="Thermo" ID="ID_1788257293" CREATED="1611782647376" MODIFIED="1611782652362">
<node TEXT="thermo()" ID="ID_1748791085" CREATED="1611782815004" MODIFIED="1611782817037">
<node TEXT="Return the thermophysical model" ID="ID_282860280" CREATED="1611782840194" MODIFIED="1611782840194"/>
</node>
<node TEXT="thermoRef()" ID="ID_470734436" CREATED="1611782818068" MODIFIED="1611782821394">
<node TEXT="Access the thermophysical model" ID="ID_1244610972" CREATED="1611782845805" MODIFIED="1611782845805"/>
</node>
<node TEXT="rho()" ID="ID_1587404578" CREATED="1611782822320" MODIFIED="1611782824130">
<node TEXT="Return the density field" ID="ID_1783865242" CREATED="1611782849673" MODIFIED="1611782849673"/>
</node>
<node TEXT="isothermal()" ID="ID_1560034303" CREATED="1611782824608" MODIFIED="1611782828869">
<node TEXT="Return whether the phase is isothermal" ID="ID_1331061704" CREATED="1611782853613" MODIFIED="1611782853613"/>
</node>
<node TEXT="heEqn()" ID="ID_987818243" CREATED="1611782829132" MODIFIED="1611782834550">
<node TEXT="Return the enthalpy equation" ID="ID_1108105411" CREATED="1611782858469" MODIFIED="1611782858469"/>
</node>
</node>
<node TEXT="Species" ID="ID_1790212272" CREATED="1611782653893" MODIFIED="1611782659028">
<node TEXT="pure()" ID="ID_1148212513" CREATED="1611782910621" MODIFIED="1611782912601">
<node TEXT="Return whether the phase is pure (i.e., not multi-component)" ID="ID_1153098725" CREATED="1611782954134" MODIFIED="1611782954134"/>
</node>
<node TEXT="YiEqn()" ID="ID_978235920" CREATED="1611782912912" MODIFIED="1611782918265">
<node TEXT="Return the species fraction equation" ID="ID_1492090940" CREATED="1611782957774" MODIFIED="1611782957774"/>
</node>
<node TEXT="Y()" FOLDED="true" ID="ID_1456555948" CREATED="1611782920533" MODIFIED="1611782922133">
<node TEXT="Return the species mass fractions" ID="ID_1380031101" CREATED="1611782966007" MODIFIED="1611782966007"/>
</node>
<node TEXT="Y(name)" ID="ID_35353390" CREATED="1611783001836" MODIFIED="1611783006662">
<node TEXT="Return a species mass fraction by name" ID="ID_273551901" CREATED="1611783008185" MODIFIED="1611783008185"/>
</node>
<node TEXT="YRef()" ID="ID_112327175" CREATED="1611782922360" MODIFIED="1611782926425">
<node TEXT="Access the species mass fractions" ID="ID_162472245" CREATED="1611783014610" MODIFIED="1611783014610"/>
</node>
<node TEXT="YActive()" ID="ID_675642777" CREATED="1611782926684" MODIFIED="1611782932177">
<node TEXT="Return the active species mass fractions" ID="ID_403173479" CREATED="1611783018822" MODIFIED="1611783018822"/>
</node>
<node TEXT="YActiveRef()" ID="ID_1036672133" CREATED="1611782932480" MODIFIED="1611782943017">
<node TEXT="Access the active species mass fractions" ID="ID_1165061634" CREATED="1611783022969" MODIFIED="1611783022969"/>
</node>
<node TEXT="R()" ID="ID_213865188" CREATED="1611782943388" MODIFIED="1611782947942">
<node TEXT="Return the fuel consumption rate matrix" ID="ID_1898441970" CREATED="1611783027914" MODIFIED="1611783027914"/>
</node>
</node>
<node TEXT="Momentum" ID="ID_261806333" CREATED="1611782660476" MODIFIED="1611782661857">
<node TEXT="stationary()" ID="ID_1670824145" CREATED="1611783041656" MODIFIED="1611783045041">
<node TEXT="Return whether the phase is stationary" ID="ID_273036191" CREATED="1611783128694" MODIFIED="1611783128694"/>
</node>
<node TEXT="UEqn()" ID="ID_490961232" CREATED="1611783045425" MODIFIED="1611783051574">
<node TEXT="Return the momentum equation" ID="ID_1705803080" CREATED="1611783132754" MODIFIED="1611783132754"/>
</node>
<node TEXT="UfEqn()" ID="ID_1404805463" CREATED="1611783051849" MODIFIED="1611783056285">
<node TEXT="Return the momentum equation for the face-based algorithm" ID="ID_745055876" CREATED="1611783136610" MODIFIED="1611783136610"/>
</node>
<node TEXT="U()" ID="ID_309517097" CREATED="1611783057500" MODIFIED="1611783059062">
<node TEXT="Return the velocity" ID="ID_1244718411" CREATED="1611783140650" MODIFIED="1611783140650"/>
</node>
<node TEXT="URef()" ID="ID_1304842852" CREATED="1611783059440" MODIFIED="1611783061901">
<node TEXT="Access the velocity" ID="ID_1546413156" CREATED="1611783144770" MODIFIED="1611783144770"/>
</node>
<node TEXT="phi()" ID="ID_473226680" CREATED="1611783063180" MODIFIED="1611783070413">
<node TEXT="Return the volumetric flux" ID="ID_563077129" CREATED="1611783148674" MODIFIED="1611783148674"/>
</node>
<node TEXT="phiRef()" ID="ID_1784661774" CREATED="1611783071489" MODIFIED="1611783075489">
<node TEXT="Access the volumetric flux" ID="ID_1050283816" CREATED="1611783154806" MODIFIED="1611783154806"/>
</node>
<node TEXT="alphaPhi()" ID="ID_1420800559" CREATED="1611783075761" MODIFIED="1611783079571">
<node TEXT="Return the volumetric flux of the phase" ID="ID_552901861" CREATED="1611783160255" MODIFIED="1611783160255"/>
</node>
<node TEXT="alphaPhiRef()" ID="ID_1031382713" CREATED="1611783080512" MODIFIED="1611783087110">
<node TEXT="Access the volumetric flux of the phase" ID="ID_602149818" CREATED="1611783164910" MODIFIED="1611783164910"/>
</node>
<node TEXT="alphaRhoPhi()" ID="ID_1554949506" CREATED="1611783088921" MODIFIED="1611783093401">
<node TEXT="Return the mass flux of the phase" ID="ID_684766953" CREATED="1611783168914" MODIFIED="1611783168914"/>
</node>
<node TEXT="alphaRhoPhiRef()" ID="ID_565462271" CREATED="1611783094045" MODIFIED="1611783103969">
<node TEXT="Access the mass flux of the phase" ID="ID_1918557045" CREATED="1611783173855" MODIFIED="1611783173855"/>
</node>
<node TEXT="DUDt()" ID="ID_1402586714" CREATED="1611783104304" MODIFIED="1611783108746">
<node TEXT="Return the substantive acceleration" ID="ID_1247496730" CREATED="1611783178839" MODIFIED="1611783178839"/>
</node>
<node TEXT="DUDtf()" ID="ID_1469490227" CREATED="1611783108964" MODIFIED="1611783112905">
<node TEXT="Return the substantive acceleration on the faces" ID="ID_132142389" CREATED="1611783184126" MODIFIED="1611783184126"/>
</node>
<node TEXT="continuityError()" ID="ID_134379268" CREATED="1611783114433" MODIFIED="1611783121628">
<node TEXT="Return the continuity error" ID="ID_1942017412" CREATED="1611783192266" MODIFIED="1611783192266"/>
</node>
<node TEXT="K()" ID="ID_1635727041" CREATED="1611783198680" MODIFIED="1611783206439">
<node TEXT="Return the phase kinetic energy" ID="ID_350925349" CREATED="1611783208269" MODIFIED="1611783214335"/>
</node>
</node>
<node TEXT="Transport" ID="ID_771369655" CREATED="1611782669388" MODIFIED="1611782670101">
<node TEXT="kappaEff()" ID="ID_1340162375" CREATED="1611783370861" MODIFIED="1611783376289">
<node TEXT="Effective thermal turbulent diffusivity for temperature of mixture for patch [W/m/" ID="ID_1266754182" CREATED="1611783395966" MODIFIED="1611783405019"/>
</node>
<node TEXT="k()" ID="ID_1201223380" CREATED="1611783377469" MODIFIED="1611783378974">
<node TEXT="Return the turbulent kinetic energy" ID="ID_1897513337" CREATED="1611783409274" MODIFIED="1611783409274"/>
</node>
<node TEXT="pPrime()" ID="ID_69098156" CREATED="1611783379293" MODIFIED="1611783384010">
<node TEXT="Return the phase-pressure&apos;, (derivative of phase-pressure w.r.t. phase-fraction)" ID="ID_1141273632" CREATED="1611783420810" MODIFIED="1611783436319"/>
</node>
</node>
</node>
</node>
</node>
<node TEXT="phaseModel.C" POSITION="right" ID="ID_535523031" CREATED="1611781887543" MODIFIED="1611784684343" BACKGROUND_COLOR="#ccffff">
<edge COLOR="#0000ff"/>
<richcontent TYPE="DETAILS">

<html>
  

  <head>

  </head>
  <body>
  </body>
</html>

</richcontent>
</node>
</node>
</map>
