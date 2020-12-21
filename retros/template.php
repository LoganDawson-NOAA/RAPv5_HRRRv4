<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
<title>EMC RAPv5/HRRRv4 Retrospectives</title>
<link rel="stylesheet" type="text/css" href="../../settings/rap_hrrr.css">

<meta http-equiv="content-type" content="text/html; charset=utf-8" />

<!--<meta http-equiv="Content-Style-Type" content="text/css" />-->

<!--- IMPORTANT PART THAT WILL BE MODIFIED ONCE THE JAVASCRIPT REFRESH PART GETS ADDED IN -->
	
<!--<meta http-equiv="refresh" content="600" />-->

<!--<link rel="stylesheet" href="http://www.emc.ncep.noaa.gov/mmb/mmbpll/css/href.css" type="text/css" media="all" />-->

<script src="../../settings/looper.js" type="text/javascript"></script>
<script src="../../settings/slidingBox.js" type="text/javascript"></script>
<script src="../../settings/keycommands.js" type="text/javascript"></script>
<script src="https://d3js.org/d3.v4.min.js"></script>

</head>

<?php
$randomtoken = base64_encode( openssl_random_pseudo_bytes(32));
$_SESSION['csrfToken']=$randomtoken;
?>

<body>

<!-- Head element -->
<!-- <div class="page-top">
        <span><a style="color:#ffffff">RAPv5/HRRRv4 Graphics Comparison</a></span>
</div> -->

<!-- div for each dropdown menu-->
<div id="divpps">
  <select id="pps"></select>
  <div id="afterpps"></div>
</div>
<div id="divcsa">
  <select id="csa"></select>
  <div id="aftercsa"></div>
</div>
<div id="divsf">
  <select id="sf"></select>
  <div id="aftersf"></div>
</div>
<div id="divcld">
  <select id="cld"></select>
  <div id="aftercld"></div>
</div>
<div id="divupa">
  <select id="upa"></select>
  <div id="afterupa"></div>
</div>
<div id="divdom">
  <select id="dom"></select>
  <div id="afterdom"></div>
</div>
<div id="divcyc">
  <select id="cyc"></select>
  <div id="aftercyc"></div>
</div>
<div id="divmod">
  <select id="mod"></select>
  <div id="aftermod"></div>
</div>

<table id="prodMenu">

<tr><td>

<div id="mainDiv">

<div id="pageContents">

 <div id="keyCodeMenu" style="display:none;opacity:.9;filter:alpha(opacity=90);background-color:#FFFFFF;position:absolute;left:700px;top:480px;border:1px solid #DDDDDD;">
 <div style="width:200px;background-color:#EEEEEE;border-bottom:1px solid #997033;"><a href="javascript:open_key()"><b>Close</b> (or press 1 again)</a></div>
 <div style="display:none" id="loadingBox">Images loading, please wait . . .</div>

Type in these keys to control the <br />

         animation from your keyboard. <br />
         You do not need to press the shift key. <br />
         <table> 

            <tr>
               <td>
                  <b>p</b> Play<br />
                  <b>o</b> Stop<br />
                  <b>+</b> Increase Speed<br />
                  <b>-</b> Decrease Speed<br />
                  <b>]</b> Forward<br />
                  <b>[</b> Reverse<br />
               </td>

               <td>
                  <b>r</b> Rock (on/off)<br />
                  <b>></b> Forward Step<br />
                  <b>&#60;</b> Back Step<br /> 
                  <b>d</b> Dwell (on/off)<br />
                  <b>w</b> increase Dwell<br /> 
                  <b>q</b> decrease Dwell<br />
               </td>
            </tr>
         </table>
         <b>m</b> Turn on/off rollover times
      </div>

   <!--PAGE CONTENTS -->

<form name="form" method="post" action="">
<input type="hidden" id="fh" name="csrfToken" value='<?php echo($_SESSION["csrfToken"]) ?>'/>

     <table border="0" cellspacing="0" cellpadding="0" width="125%" bgcolor="#D5D5D5" style="border-bottom:1px solid #DDDDDD">
        <tr>
	  <td style="padding:0 2px 0 2px;text-align:center;">
            <input value="Start" onclick="start_play();" name="button3" type="button" class="mapbutton" style="float:left" />
            <input value="Stop" onclick="stop_play();" name="button" type="button" class="mapbutton" style="float:right" />
          </td>

          <td bgcolor="" style="border-left:1px solid #BBBBBB;text-align:center;"><font style="font-size:10px;"><b>FRM</b></font></td>
	  <td>
            <input value="" name="frame" size="6" type="text" class="maptext" />
	  </td>

          <td bgcolor="" style="border-left:1px solid #BBBBBB;text-align:center;"><font style="font-size:10px;"><b>SPD</b></font></td>
	  <td> 
            <input value="&lt;&lt;" onclick="delay=delay*inc; show_delay();" class="mapbutton2" name="button2" type="button" />
            <input value="&gt;&gt;" onclick="delay=delay/inc; show_delay();" class="mapbutton2" name="button2" type="button" />
	    <input value="" name="dly" size="2" type="text" class="maptext" />img/sec
	  </td>

          <td bgcolor="" style="border-left:1px solid #BBBBBB;text-align:center;"><font style="font-size:10px;"><b>STEP (Keyboard Arrows)</b></font></td>
          <td style="text-align:center">
            <input value=" &lt; " onclick="backstep();" name="button2" type="button" class="mapbutton2" />
            <input value=" &gt; " onclick="forwardstep();" name="button2" type="button" class="mapbutton2" />
          </td>

	</tr>
    </table>

</table>



<table id="prodMenu2" border="0" cellspacing="0" cellpadding="0" width="125%" bgcolor="#D5D5D5" style="border-bottom:1px solid #DDDDDD; display: none">
       <tr>        
          <td bgcolor="" style="border-left:1px solid #BBBBBB;text-align:center;"><font style="font-size:10px;"><b>DWELL</b></font></td>
          <td>
           <label><input type="checkbox" name="dwell" onclick="show_dwell();" checked="checked" />Dwell</label> 
            <input type="button" value=" - " onclick="decDwell(); show_dwell();" class="mapbutton2" />
            <input type="button" value=" + " onclick="incDwell(); show_dwell();" class="mapbutton2" />
            <input type="text" size="2" name="dwl" class="maptext" />sec
          </td>

          <td bgcolor="" style="border-left:1px solid #BBBBBB;text-align:center;"><font style="font-size:10px;"><b>STEP (Keyboard Arrows)</b></font></td>
	  <td style="text-align:center">
            <input value=" &lt; " onclick="backstep();" name="button2" type="button" class="mapbutton2" />
            <input value=" &gt; " onclick="forwardstep();" name="button2" type="button" class="mapbutton2" />
	  </td>

	  <td><label><input type="checkbox" checked="checked" name="nostep"/>On/Off</label></td>

	</tr>
</table>

<div class="mapbox">
<center>
<img src="https://www.emc.ncep.noaa.gov/users/meg/rapv5_hrrrv4/retros/%RETRO_PERIOD%/%CASE%/hrrr/%DAY1_00Z_CYCLE%/comparerefc_conus_f00.png" name="animation" alt="Loading..." ondblclick="document.location.href = this.src" />
</center>
<br />
<div id="text1" style="position: absolute; left:225px;">
<span>
</span>
</div>
 </div>




<!--  Delay (ms): <INPUT TYPE=text VALUE="" NAME="delay" SIZE=6> //-->
</form>

<script type='text/javascript'>

//allow the use of arrow keys for advancing images
window.addEventListener("keydown", function(e) {
    if([37,39].indexOf(e.keyCode) > -1) {
        e.preventDefault();
        }
});
document.onkeydown = checkKey;
function checkKey(e) {
  e = e || window.event;
  if (e.keyCode == '37') {
     backstep();
    }
  if (e.keyCode == '39') {
      forwardstep();
    }
  }

//variables and functions for dropping menu if not already dropped and 
//retracting menu if it is already dropped when clicked
var dropped=0;
var droppedpps=0;
var droppedcsa=0;
var droppedsf=0;
var droppedcld=0;
var droppedupa=0;
var droppeddom=0;
var droppedcyc=0;
var droppedmod=0;

function change() {
  if (dropped==0) {
    afterdate.attr("class","select-styled active");
    uldate.style("display","block")
    dropped=1;
   }
   else {
    dropped=0;
    uldate.style("display","none");
    afterdate.attr("class","select-styled");
    showdate.attr("class","select");};
  subul.style("display","none");
};
function changepps() {
  if (droppedpps==0) {
    afterpps.attr("class","select-styled active");
    ulpps.style("display","block")
    droppedpps=1;
   }
   else {
    droppedpps=0;
    ulpps.style("display","none");
    afterpps.attr("class","select-styled");
    showpps.attr("class","selectpps");};
};
function changecsa() {
  if (droppedcsa==0) {
    aftercsa.attr("class","select-styled active");
    ulcsa.style("display","block")
    droppedcsa=1;
   }
   else {
    droppedcsa=0;
    ulcsa.style("display","none");
    aftercsa.attr("class","select-styled");
    showcsa.attr("class","selectpps");};
};
function changesf() {
  if (droppedsf==0) {
    aftersf.attr("class","select-styled active");
    ulsf.style("display","block")
    droppedsf=1;
   }
   else {
    droppedsf=0;
    ulsf.style("display","none");
    aftersf.attr("class","select-styled");
    showsf.attr("class","selectpps");};
};
function changecld() {
  if (droppedcld==0) {
    aftercld.attr("class","select-styled active");
    ulcld.style("display","block")
    droppedcld=1;
   }
   else {
    droppedcld=0;
    ulcld.style("display","none");
    aftercld.attr("class","select-styled");
    showcld.attr("class","selectpps");};
};
function changeupa() {
  if (droppedupa==0) {
    afterupa.attr("class","select-styled active");
    ulupa.style("display","block")
    droppedupa=1;
   }
   else {
    droppedupa=0;
    ulupa.style("display","none");
    afterupa.attr("class","select-styled");
    showupa.attr("class","selectpps");};
};
function changedom() {
  if (droppeddom==0) {
    afterdom.attr("class","select-styled active");
    uldom.style("display","block")
    droppeddom=1;
   }
   else {
    droppeddom=0;
    uldom.style("display","none");
    afterdom.attr("class","select-styled");
    showdom.attr("class","selectdom");};
};

function changecyc() {
  if (droppedcyc==0) {
    aftercyc.attr("class","select-styled active");
    ulcyc.style("display","block")
    droppedcyc=1;
   }
   else {
    droppedcyc=0;
    ulcyc.style("display","none");
    aftercyc.attr("class","select-styled");
    showcyc.attr("class","selectcycs");};
};

function changemod() {
  if (droppedmod==0) {
    aftermod.attr("class","select-styled active");
    ulmod.style("display","block")
    droppedmod=1;
   }
   else {
    droppedmod=0;
    ulmod.style("display","none");
    aftermod.attr("class","select-styled");
    showmod.attr("class","selectmods");};
};


//a bunch of arrays and variables to keep track of which image is currently being looped
var path="https://www.emc.ncep.noaa.gov/users/meg/rapv5_hrrrv4/retros/%RETRO_PERIOD%/%CASE%";
var field="refc";
var dom="conus";
var mod="hrrr";
var cyc="%DAY1_00Z_CYCLE%";
var imax=37;
var datelist=[""];

var ppslist=["Total Precip","3-h Precip","Precip Type","Accumulated Snow Depth","6-hr Change in Snow Depth"];
var ppsdata=["qpf","qpf3","ptype","snow","snow6"];
var imaxespps=[37,37,37,7,7];

var csalist=["Composite Reflectivity","1-km AGL Reflectivty","1-h Max 2-5 km UH","Run-Total Max 2-5 km UH","Surface-Based CAPE/CIN","Mixed Layer CAPE/CIN","Most Unstable CAPE/CIN","0-6 km Vertical Wind Shear","0-3 km Storm Relative Helicity","0-1 km Storm Relative Helicity"];
var csadata=["refc","ref1km","uh25","uh25_accum","sbcape","mlcape","mucape","shr6km","srh3km","srh1km"];  
var imaxescsa=[37,37,37,37,37,37,37,37,37,37];

var sflist=["Sea Level Pressure","2-m Temperature","Skin Temperature","2-m Dewpoint Temperature","10-m Wind Speed","PBL Height"];
var sfdata=["slp","2mt","tsfc","2mdew","10mwind","hpbl"];
var imaxessf=[37,37,37,37,37,37];

var cldlist=["Precipitable Water","Ceiling","Visibility"];
var clddata=["pw","ceil","vis"];
var imaxescld=[37,37,37];

//var upalist=["850-hPa Theta-e","850-hPa Temperatures","850-hPa Heights/Winds","700-hPa Omega & RH","500-hPa Heights/Winds/Vorticity","250-hPa Winds"];
//var upadata=["850t","850temp","850wind","700","500","250wind"];
//var imaxesupa=[37,37,37,37,37,37];
var upalist=["N/A"];
var upadata=[""];
var imaxesupa=[0];

var domlist=["CONUS","Southern Plains","Central Plains","Northern Plains","Midwest","Northeast","Southeast"];
var domdata=["conus","splains","cplains","nplains","midwest","northeast","southeast"];

var cyclist=["12Z %DATE_M1%","18Z %DATE_M1%","00Z %DATE%","03Z %DATE%","06Z %DATE%","09Z %DATE%","12Z %DATE%","15Z %DATE%","18Z %DATE%","21Z %DATE%","00Z %DATE_P1%"];
var cycdata=["%ymd_m1%12","%ymd_m1%18","%ymd%00","%ymd%03","%ymd%06","%ymd%09","%ymd%12","%ymd%15","%ymd%18","%ymd%21","%ymd_p1%00"];
//var cyclist=["2019092300","2019080712","2019080700","2019080612","2019080500"]
//var cycdata=["..\/cyc\/rap","..\/cycm1\/rap","..\/cycm2\/rap","..\/cycm3\/rap","..\/cycm4\/rap"]
//var cycdata=["cyc\/nam","cycm1\/nam","cycm2\/nam","cycm3\/nam","cycm4\/nam","cycm5\/nam"]
///home/people/emc/www/htdocs/mmb/bblake/fv3

//var modlist=["HRRR","RAP"]
//var moddata=[".\/hrrr",".\/rap"]
var modlist=["HRRR"];
var moddata=[".\/hrrr"];


//function for loading images
function loadimages() {
  window.pauseOnStart = true;
  window.pauseWhere = 0;

  if (mod == ".\/hrrr" || mod == "hrrr"){ 
    var imaxespps=[37,37,37,7,7];
    var imaxescsa=[37,37,37,37,37,37,37,37,37];
    var imaxessf=[37,37,37,37,37,37];
    var imaxescld=[37,37,37];
//  var imaxesupa=[37,37,37,37,37,37];
    var imaxesupa=[0];
  }
  else if (mod == ".\/rap" || mod == "rap"){
    var imaxespps=[40,40,40,7,7];
    var imaxescsa=[40,40,40,40,40,40,40,40,40];
    var imaxessf=[40,40,40,40,40,40];
    var imaxescld=[40,40,40];
//  var imaxesupa=[40,40,40,40,40,40];
    var imaxesupa=[0];
  };


  if (imax==61) {
    var fhours=["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60"];
  }
  if (imax==60) {
    var fhours=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60"];
  }
  if (imax==37) {
    var fhours=["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"];
  }
  if (imax==40) {
    var fhours=["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39"];
  }
  if (imax==29) {
    var fhours=["00","01","02","03","04","05","06","07","08","09","10","11","12","15","18","21","24","27","30","33","36","39","42","45","48","51","54","57","60"];
  }
  if (imax==22) {
    var fhours=["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21"];
  }
  if (imax==21) {
    var fhours=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21"];
  }
  if (imax==28) {
    var fhours=["01","02","03","04","05","06","07","08","09","10","11","12","15","18","21","24","27","30","33","36","39","42","45","48","51","54","57","60"];
  }
  if (imax==20) {
    var fhours=["03","06","09","12","15","18","21","24","27","30","33","36","39","42","45","48","51","54","57","60"];
  }
  if (imax==19) {
    var fhours=["06","09","12","15","18","21","24","27","30","33","36","39","42","45","48","51","54","57","60"];
  }
  if (imax==16) {
    var fhours=["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15"];
  }
  if (imax==7) {
    var fhours=["00","06","12","18","24","30","36"];
  }
  if (imax==0) {
    var fhours=["00"];
  }
  window.temp_list=[];
  window.temp_list = new Array(imax);
  for (i=0;i<imax;i++) {
  window.temp_list[i]=path+'/'+mod+'/'+cyc+'/compare'+field+'_'+dom+'_f'+String(fhours[i])+'.png'
  }
  initialize_looper();
};

//initial call for loading images
loadimages();

//precipitation, p-type, snow menu creation
var hidepps=d3.select("#pps")
  .attr("class","select-hidden");
hidepps.selectAll("option")
  .data(ppslist)
  .enter().append("option")
    .text(function(d) { return d; });
hidepps.property("value","precipitation|p-type");
var showpps=d3.select("#divpps")
  .attr("class","selectpps")
  .style("left","5px");
var afterpps=d3.select("#afterpps")
  .attr("class","select-styled")
  .text("Precipitation|P-Type|Snow")
  .on("click",changepps);
afterpps.property("value","precipitation|p-type|snow");
var ulpps=d3.select("#divpps").append("ul")
  .attr("class","uldate");
ulpps.selectAll("li")
  .data(ppslist)
  .enter().append("li")
    .attr("class","mainlipps")
    .text(function(d,i) { return ppslist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
ulpps.style("display","none");
ulpps.selectAll(".mainlipps").on("click",function(d) {
  for (i=0;i<ppslist.length;i++) {
    if (ppslist[i]==d) {
      field=ppsdata[i];
      imax=imaxespps[i];
    }
   }
    droppedpps=0;
    ulpps.style("display","none");
    afterpps.attr("class","select-styled");
    showpps.attr("class","selectpps");
  loadimages();
})

//convection, severe, aviation menu creation
var hidecsa=d3.select("#csa")
  .attr("class","select-hidden");
hidecsa.selectAll("option")
  .data(csalist)
  .enter().append("option")
    .text(function(d) { return d; });
hidecsa.property("value","convection|severe");
var showcsa=d3.select("#divcsa")
  .attr("class","selectpps")
  .style("left","215px");
var aftercsa=d3.select("#aftercsa")
  .attr("class","select-styled")
  .text("Convection|Severe Wx")
  .on("click",changecsa);
aftercsa.property("value","convection|severe");
var ulcsa=d3.select("#divcsa").append("ul")
  .attr("class","uldate");
ulcsa.selectAll("li")
  .data(csalist)
  .enter().append("li")
    .attr("class","mainlipps")
    .text(function(d,i) { return csalist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
ulcsa.style("display","none");
ulcsa.selectAll(".mainlipps").on("click",function(d) {
  for (i=0;i<csalist.length;i++) {
    if (csalist[i]==d) {
      field=csadata[i];
      imax=imaxescsa[i];
    }
   }
    droppedcsa=0;
    ulcsa.style("display","none");
    aftercsa.attr("class","select-styled");
    showcsa.attr("class","selectpps");
  loadimages();
})

//surface fields menu creation

var hidesf=d3.select("#sf")
  .attr("class","select-hidden");
hidesf.selectAll("option")
  .data(sflist)
  .enter().append("option")
    .text(function(d) { return d; });
hidesf.property("value","surface");
var showsf=d3.select("#divsf")
  .attr("class","selectpps")
  .style("left","425px");
var aftersf=d3.select("#aftersf")
  .attr("class","select-styled")
  .text("Surface Fields")
  .on("click",changesf);
aftersf.property("value","surface fields");
var ulsf=d3.select("#divsf").append("ul")
  .attr("class","uldate");
ulsf.selectAll("li")
  .data(sflist)
  .enter().append("li")
    .attr("class","mainlipps")
    .text(function(d,i) { return sflist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
ulsf.style("display","none");
ulsf.selectAll(".mainlipps").on("click",function(d) {
  for (i=0;i<sflist.length;i++) {
    if (sflist[i]==d) {
      field=sfdata[i];
      imax=imaxessf[i];
    }
   }
    droppedsf=0;
    ulsf.style("display","none");
    aftersf.attr("class","select-styled");
    showsf.attr("class","selectpps");
  loadimages();
})


//cloud variables
var hidecld=d3.select("#cld")
  .attr("class","select-hidden");
hidecld.selectAll("option")
  .data(cldlist)
  .enter().append("option")
    .text(function(d) { return d; });
hidecld.property("value","Cloud");
var showcld=d3.select("#divcld")
  .attr("class","selectpps")
  .style("left","635px");
var aftercld=d3.select("#aftercld")
  .attr("class","select-styled")
  .text("Cloud|Sfc Fluxes")
  .on("click",changecld);
aftercld.property("value","Cloud");
var ulcld=d3.select("#divcld").append("ul")
  .attr("class","uldate");
ulcld.selectAll("li")
  .data(cldlist)
  .enter().append("li")
    .attr("class","mainlipps")
    .text(function(d,i) { return cldlist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
ulcld.style("display","none");
ulcld.selectAll(".mainlipps").on("click",function(d) {
  for (i=0;i<cldlist.length;i++) {
    if (cldlist[i]==d) {
      field=clddata[i];
      imax=imaxescld[i];
    }
   }
    droppedcld=0;
    ulcld.style("display","none");
    aftercld.attr("class","select-styled");
    showcld.attr("class","selectpps");
  loadimages();
})

//upper level variables

var hideupa=d3.select("#upa")
  .attr("class","select-hidden");
hideupa.selectAll("option")
  .data(upalist)
  .enter().append("option")
    .text(function(d) { return d; });
hideupa.property("value","Upper Air");
var showupa=d3.select("#divupa")
  .attr("class","selectpps")
  .style("left","845px");
var afterupa=d3.select("#afterupa")
  .attr("class","select-styled")
  .text("Upper Air Fields")
  .on("click",changeupa);
afterupa.property("value","Upper Air");
var ulupa=d3.select("#divupa").append("ul")
  .attr("class","uldate");
ulupa.selectAll("li")
  .data(upalist)
  .enter().append("li")
    .attr("class","mainlipps")
    .text(function(d,i) { return upalist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
ulupa.style("display","none");
ulupa.selectAll(".mainlipps").on("click",function(d) {
  for (i=0;i<upalist.length;i++) {
    if (upalist[i]==d) {
      field=upadata[i];
      imax=imaxesupa[i];
    }
   }
    droppedupa=0;
    ulupa.style("display","none");
    afterupa.attr("class","select-styled");
    showupa.attr("class","selectpps");
  loadimages();
})

//region menu creation
var hidedom=d3.select("#dom")
  .attr("class","select-hidden");
hidedom.selectAll("option")
  .data(domlist)
  .enter().append("option")
    .text(function(d) { return d; });
hidedom.property("value","region");
var showdom=d3.select("#divdom")
  .attr("class","selectdom")
  .style("left","655px");
var afterdom=d3.select("#afterdom")
  .attr("class","select-styled")
  .text("Region")
  .on("click",changedom);
afterdom.property("value","region");
var uldom=d3.select("#divdom").append("ul")
  .attr("class","uldate");
uldom.selectAll("li")
  .data(domlist)
  .enter().append("li")
    .attr("class","mainli")
    .text(function(d,i) { return domlist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
uldom.style("display","none");
uldom.selectAll(".mainli").on("click",function(d) {
  for (i=0;i<domlist.length;i++) {
    if (domlist[i]==d) {
      dom=domdata[i];
    }
   }
    droppeddom=0;
    uldom.style("display","none");
    afterdom.attr("class","select-styled");
    showdom.attr("class","selectdom");
  loadimages();
})

//cycle menu creation
var hidecyc=d3.select("#cyc")
  .attr("class","select-hidden");
hidecyc.selectAll("option")
  .data(cyclist)
  .enter().append("option")
    .text(function(d) { return d; });
hidecyc.property("value","cycle");
var showcyc=d3.select("#divcyc")
  .attr("class","selectcycs")
  .style("left","445px");
var aftercyc=d3.select("#aftercyc")
  .attr("class","select-styled")
  .text("Cycle")
  .on("click",changecyc);
aftercyc.property("value","cycle");
var ulcyc=d3.select("#divcyc").append("ul")
  .attr("class","uldate");
ulcyc.selectAll("li")
  .data(cyclist)
  .enter().append("li")
    .attr("class","mainli")
    .text(function(d,i) { return cyclist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
ulcyc.style("display","none");
ulcyc.selectAll(".mainli").on("click",function(d) {
  for (i=0;i<cyclist.length;i++) {
    if (cyclist[i]==d) {
      cyc=cycdata[i];
    }
   }
    droppedcyc=0;
    ulcyc.style("display","none");
    aftercyc.attr("class","select-styled");
    showcyc.attr("class","selectcycs");
  loadimages();
})

//model menu creation
var hidemod=d3.select("#mod")
  .attr("class","select-hidden");
hidemod.selectAll("option")
  .data(modlist)
  .enter().append("option")
    .text(function(d) { return d; });
hidemod.property("value","model");
var showmod=d3.select("#divmod")
  .attr("class","selectmods")
  .style("left","225px");
var aftermod=d3.select("#aftermod")
  .attr("class","select-styled")
  .text("Model")
  .on("click",changemod);
aftermod.property("value","model");
var ulmod=d3.select("#divmod").append("ul")
  .attr("class","uldate");
ulmod.selectAll("li")
  .data(modlist)
  .enter().append("li")
    .attr("class","mainli")
    .text(function(d,i) { return modlist[i]; })
    .attr("id",function(d) { return "MD"+String(d);});
ulmod.style("display","none");
ulmod.selectAll(".mainli").on("click",function(d) {
  for (i=0;i<modlist.length;i++) {
    if (modlist[i]==d) {
      mod=moddata[i];
    }
   }
    droppedmod=0;
    ulmod.style("display","none");
    aftermod.attr("class","select-styled");
    showmod.attr("class","selectmods");
  loadimages();
})

 
</script>

</div>
</div>

</body>
</html>
