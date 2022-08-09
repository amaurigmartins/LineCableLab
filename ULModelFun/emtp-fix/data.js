// Author: Jesus Morales
//	Last change: 01/12/2020

/* This script drives the data tab (data.htm) of the Line/Cable data device. 

   This script is divided in 9 sections:
   Section 1. Generic EMTP methods, onopen, onload, and onclose
   Section 2. General functions, control the options displayed depending on selections and inputs
   Section 3. Functions related to soil options
   Section 4. Functions related to overhead line with single wires
   Section 5. Functions related to overhead line with bundle of conductors
   Section 6. Additional options for overhead lines (database, and consider midspan and hollow)
   Section 7. Functions related to single-core cables
   Section 8. Functions related to pipe-type cables
   Section 9. Other functions (line length and additional table definition)

   In addition to this script, the function DrawCables.js is used to draw the lines/cables.
   All calls to the functions in DrawCables.js are done from the htm file and not from this script.
*/

var myglobal
var designfilepath
var mydatafile

var phases_single_wire, phases_bundle, phases_single_core, phases_pipe_cable
var previousphases_single_wire, previousphases_bundle, previousphases_single_core, previousphases_pipe_cable

var SoilLayers_handle                     //SoilLayers grid handle
var SoilLayers_handle_Set = false;

var Overhead_line_data_handle             //Overhead line data grid handle
var Overhead_line_data_handle_Set = false;

var Bundle_data_handle                    //Overhead bundle data grid handle
var Bundle_data_handle_Set = false;

var Single_core_main_handle               //Single-core cable main grid handle
var Single_core_main_handle_Set = false;
var Single_core_data_handle               //Single-core cable data grid handle
var Single_core_data_handle_Set = false;

var Pipe_cable_main_handle                //Pipe-type cable main grid handle
var Pipe_cable_main_handle_Set = false;
var Pipe_cable_data_handle                //Pipe-type cable data grid handle
var Pipe_cable_data_handle_Set = false;

var units
var linecable_length, soil_resistivity, soil_permeability, soil_permittivity
var old_linecable_length, old_soil_resistivity, old_soil_permeability, old_soil_permittivity

var old_selectorState
var Rin, Rout, cable, conductor

var n_soil_layers, previous_n_soil_layers
var previous_multilayerSoil_checked



//*************************************************************************************
//  Section 1. Generic EMTP methods, onopen, onload and, onclose
//************************************************************************************ */

// DESCRIPTION: Called first, when opening the html page
function onopen(global){

   //initialize resolution issues
   standard_onload_GlobalDataInit(global); //this function is in standard_onload
  
   //*Set the directory of the design, used for file dialog initdir
   designfilepath = global.designfilepath; 
   project_folder_name = global.project_folder_name; 
   mydatafile = global.SPFile();   //Give an SPFile object
   myglobal = global;

   return true;
   
}


//+ DESCRIPTION: Called second. Called only once when we move to this tab the first time
// It first checks for correct input values, then call the corresponding functions for setting the data tables. 
function LineCableData_onload(){

   phases_single_wire = parseInt(Nphases_single_wire.value);  // auxiliar variable used all over the code 
   if (isNaN(phases_single_wire)) {    // check if number of conductors saved in FormData is valid
      phases_single_wire = 0;          
      Nphases_single_wire.value = '0';   
   }
   previousphases_single_wire = phases_single_wire;

   phases_bundle = parseInt(Nphases_bundle.value);  // auxiliar variable used all over the code 
   if (isNaN(phases_bundle)) {   // check if number of conductors saved in FormData is valid
      phases_bundle = 0;         
      Nphases_bundle.value = '0';   
   }
   previousphases_bundle = phases_bundle;

   phases_single_core = parseInt(Nphases_single_core.value);  // auxiliar variable used all over the code 
   if (isNaN(phases_single_core)) {   // check if number of conductors saved in FormData is valid
      phases_single_core = 0;        
      Nphases_single_core.value = '0';   
   }
   previousphases_single_core = phases_single_core;

   //* Number of pipe type cables
   phases_pipe_cable = parseInt(Nphases_pipe_cable.value);  // auxiliar variable used all over the code 
   if (isNaN(phases_pipe_cable)) {   // check if number of conductors saved in FormData is valid
      phases_pipe_cable = 0;        
      Nphases_pipe_cable.value = '0';   
   }
   previousphases_pipe_cable = phases_pipe_cable;

   previous_multilayerSoil_checked = MultilayerSoil.checked;

   n_soil_layers = parseFloat(SoilLayers.value);      // auxiliar variable used all over the code
   if( isNaN(n_soil_layers) || n_soil_layers <= 2 ) {  // check if value saved in FormData is valid
      SoilLayers.value = '2';   // if invalid, reset to default 
      n_soil_layers = 2;
   }
   previous_n_soil_layers = SoilLayers.value;

   soil_resistivity = parseFloat(Resistivity.value);  // auxiliar variable 
   if( isNaN(soil_resistivity) || soil_resistivity <= 0 ) { 
      Resistivity.value = '100';   // if invalid, reset to default 
   }
   old_soil_resistivity = Resistivity.value;

   soil_permeability = parseFloat(soil_mu.value);  // auxiliar variable 
   if( isNaN(soil_permeability) || soil_permeability <= 0 ) {    // check if value saved in FormData is valid
      soil_mu.value = '1';     // if invalid, reset to default 
   }
   old_soil_permeability = soil_mu.value;

   soil_permittivity = parseFloat(soil_epsilon.value);  // auxiliar variable 
   if( isNaN(soil_permittivity) || soil_permittivity <= 0 ) { // check if value saved in FormData is valid
      soil_epsilon.value = '1';   // if invalid, reset to default 
   }
   old_soil_permittivity = soil_epsilon.value;

   linecable_length = parseFloat(Length.value);   // auxiliar variable 
   if( isNaN(linecable_length) || linecable_length <= 0 ) {   // check if value saved in FormData is valid
      Length.value = '100';     // if invalid, reset to default 
   }
   old_linecable_length = Length.value;
   
   old_selectorState = SelectedType.value;    // save the initial selection made (line/cable/both)

   set_LineLengthUnits();        // Set the label of units (metric or english) for the line_length input  
   setOptions_MainTable();       // Set the options in the main table 
   set_ListOfTables();           // Set the list of tables containing data (list of checkboxes)

   // Since this function is only called at opening, check all the existing checkboxes from the list of tables
   if(MultilayerSoil.checked){
      document.getElementById("SoilLayersbox").checked = true;
   }
   if(phases_single_wire>0){
      document.getElementById("Overheadsinglebox").checked = true;
   }
   if(phases_bundle>0){
      document.getElementById("Overheadbundlebox").checked = true;
   }
   if(phases_single_core>0){
      document.getElementById("Singlecorebox").checked = true;
   }
   if(phases_pipe_cable>0){
      document.getElementById("Pipetypebox").checked = true;
   }

   // If at least one table exists, show the list of tables (all its items)
   if(phases_single_wire+phases_bundle+phases_single_core+phases_pipe_cable > 0){
      document.getElementById("List_of_Tables_section").style.display = "block" 
      document.getElementById("List_of_Tables").style.display = "block"  // show list
   }else{
      document.getElementById("List_of_Tables_section").style.display = "none" 
   }

   set_DataTables();             // Set the dhtmlx grids
   ShowHide_dataTables();
   ShowHide_canvas();
   
}


//+ DESCRIPTION: This function is called when the HTML page is closed. It is the final data check.
function onclose(global){

   //*type of data
   if(SelectedType.value=="none"){
      global.FirstDataTab_Conductor_Set=false
      return true //We accept to change the tab
   }else{
      // Nphases_single_wire Nphases_bundle  Nphases_single_core Nphases_pipe_cable
      TotalNumberOfWires=  parseInt(Nphases_single_wire.value)+
                           parseInt(Nphases_bundle.value)+
                           parseInt(Nphases_single_core.value)+
                           parseInt(Nphases_pipe_cable.value);

      if(TotalNumberOfWires==0){
         if( SelectedType.value == "pul_parameters" & pul_parameters_file.value == '' ){
            global.FirstDataTab_Conductor_Set=false
            SelectedType.value='none'
            changeInMainSelector()
            return true
         }
      }                     
      global.FirstDataTab_Conductor_Set=true;
   }
   //*Between tabs only, does not go to dwj
   global.SelectedType = SelectedType.value;

   //*+Nlayer soil, number of layers value is tested here
   // will also call and reset all tables and data
   if(test_NlayersSoil()){ 
      //found an error, fixed with a valid number of levels in the soil
      //The grid will be tested below
      return false 
   }
   
   //* n_soil_layers is global in this file, defined on top
   global.SoilLayers = n_soil_layers;
   if (MultilayerSoil.checked) {

      //*Test grid for soils
      if( !testgrids(SoilLayers_handle,true) ){
         // true means that blank values are not allowed
         alert('Empty cells or errors in table "Soil layers"!')
         return false
      }
      save_Clean_Grid_Data_Table(SoilLayers_handle,SoilLayersDataGrid,SoilLayersGrid_nlines);
      global.SoilLayersDataGrid = SoilLayersDataGrid
      
      //* Verify that conductors are not placed between medium borders
      var Nrows = SoilLayers_handle.getRowsNum();
      var Boundaries = [];
      var Y_axis = 0;  // initialize at zero
      for(var row=0; row<Nrows-1; row++){
         var idrow = SoilLayers_handle.getRowId(row)  //currently tested row
         var Thickness = SoilLayers_handle.cells(idrow,1).getValue();
         Y_axis = Y_axis - Thickness;
         Boundaries[row] = Y_axis;  // These are the boundaries between soil layers
      }

      if (phases_single_core > 0) { //-

         // 1. Check that inner radius is smaller than outer radius for each conductor
         for(var row=0; row<Single_core_main_handle.getRowsNum(); row++){
            idrow = Single_core_main_handle.getRowId(row)  //currently tested row

            var cable_n = Single_core_main_handle.cells(idrow,0).getValue();
            var CableYaxisCenter = Single_core_main_handle.cells(idrow,3).getValue();
            var Rout = Single_core_main_handle.cells(idrow,4).getValue() * 0.01;

            // Only check for underground cables
            if(CableYaxisCenter<0){

               var CableYaxis_min = CableYaxisCenter - Rout;
               var diam = 2*Rout;
               var CableYaxis_max = CableYaxis_min + diam;

               for(var h=0; h<Boundaries.length; h++){
                  var n_boundary = Boundaries[h];
                  if(n_boundary<CableYaxis_max && n_boundary>CableYaxis_min){
                     var n_layer = h+1;
                     var next_layer = h+2;
                     alert('Single-core cable '+cable_n+' is placed at the boundary between soil layers '+n_layer+' and '+next_layer+'. Please correct.')
                     return
                  }
               }

            }
            
         }
   
      }

      if (phases_pipe_cable > 0) {

         // 1. Check that inner radius is smaller than outer radius for each conductor
         for(var row=0; row<Pipe_cable_main_handle.getRowsNum(); row++){
            idrow = Pipe_cable_main_handle.getRowId(row)  //currently tested row

            var cable_n = Pipe_cable_main_handle.cells(idrow,0).getValue();
            var CableYaxisCenter = Pipe_cable_main_handle.cells(idrow,3).getValue();
            var Rout = Pipe_cable_main_handle.cells(idrow,4).getValue() * 0.01;

            // Only check for underground cables
            if(CableYaxisCenter<0){

               var CableYaxis_min = CableYaxisCenter - Rout;
               var diam = 2*Rout;
               var CableYaxis_max = CableYaxis_min + diam;

               for(var h=0; h<Boundaries.length; h++){
                  var n_boundary = Boundaries[h];
                  if(n_boundary<CableYaxis_max && n_boundary>CableYaxis_min){
                     var n_layer = h+1;
                     var next_layer = h+2;
                     alert('Pipe-type cable '+cable_n+' is placed at the boundary between soil layers '+n_layer+' and '+next_layer+'. Please correct.')
                     return
                  }
               }

            }
            
         }
   
      }

   }  //-end of code for testing soil

   //+BUNDLES: number of bundles test
   if(testNphases_bundle()){ 
      //found an error, fixed with a valid number of phases
      // this will call the DHTMLX table creation functions
     return false 
   }
   global.Nphases_bundle = phases_bundle;
   //*BUNDLES: test of bundle grid
   if (global.Nphases_bundle > 0) {
      if( !testgrids(Bundle_data_handle,true) ){
         // true means that blank values are not allowed
         alert('Empty cells or errors in table "Overhead line bundled (B) conductors"!')
         return false
      }
      save_Clean_Grid_Data_Table(Bundle_data_handle,BundleDataGrid,BundleDataGrid_nlines)
      global.BundleDataGrid = BundleDataGrid
   }

   //+SNGLE WIRES: Single-wire (W) conductors
   if(testNphases_single_wire()){ 
      //found an error, fixed with a valid number of phases
      return false 
   }
   //*SINGLE WIRES: test of single wires grid
   global.Nphases_single_wire = phases_single_wire;
   if (global.Nphases_single_wire > 0) {
      if( !testgrids(Overhead_line_data_handle,true) ){
         // true means that blank values are not allowed
         alert('Empty cells or errors in table "Overhead line single-wire (W) conductors"!')
         return false
      }
      save_Clean_Grid_Data_Table(Overhead_line_data_handle,OverheadLineDataGrid,OverheadLineDataGrid_nlines);
      global.OverheadLineDataGrid = OverheadLineDataGrid
   }
   

   //+Single-core Cable number test
   if ( testNphases_single_core() ) { 
      //found an error, fixed with a valid number of phases
      return false
   }
   global.Nphases_single_core = phases_single_core;
   //*Single-core: test grid. There are 2 grids
   if (global.Nphases_single_core > 0) {
      if( !testgrids(Single_core_main_handle,true) ){
         // true means that blank values are not allowed
         alert('Empty cells or errors in table "Single-core (SC) cable main data"!')
         return false
      }
      if( !testgrids(Single_core_data_handle,true) ){
         // true means that blank values are not allowed
         alert('Empty cells or errors in table "Single-core (SC) cable conductors/insulators data"!')
         return false
      }
      save_Clean_Grid_Data_Table(Single_core_main_handle,SingleCoreMainGrid,SingleCoreMainGrid_nlines)
      save_Clean_Grid_Data_Table(Single_core_data_handle,SingleCoreDataGrid,SingleCoreDataGrid_nlines)
   }

   //+Check for inconsistent data in single core cables
   if (global.Nphases_single_core > 0) {

      //* 1. Check that inner radius is smaller than outer radius for each conductor
      for(var row=0; row<Single_core_data_handle.getRowsNum(); row++){
         idrow = Single_core_data_handle.getRowId(row)  //currently tested row
         Rin = Single_core_data_handle.cells(idrow,3).getValue();
         Rout = Single_core_data_handle.cells(idrow,4).getValue();
         if(Rout<Rin){
            cable = Single_core_data_handle.cells(idrow,0).getValue();
            conductor = Single_core_data_handle.cells(idrow,1).getValue();
            Single_core_data_handle.selectRowById(idrow);
            alert('Inconsistent data: '+'\n"Inner radius" is larger than "Outer radius" for: \nsingle-core cable: '+cable+'\nconductor: '+conductor)
            return false
         }
      }

      //* 2. Check that all outer radius of conductors are smaller than outermost radius of cable
      for(var row=0; row<Single_core_data_handle.getRowsNum(); row++){

         idrow = Single_core_data_handle.getRowId(row)  //currently tested row
         cable =      Single_core_data_handle.cells(idrow,0).getValue();
         conductor =  Single_core_data_handle.cells(idrow,1).getValue();
         Rout =       Single_core_data_handle.cells(idrow,4).getValue();

         // Grab the cable outermost radius
         idrow_main_grid = Single_core_main_handle.getRowId(cable-1)
         outerCableRadius = Single_core_main_handle.cells(idrow_main_grid,4).getValue();

         if(outerCableRadius<Rout){
            Single_core_data_handle.selectRowById(idrow);
            Single_core_main_handle.selectRowById(idrow_main_grid);
            alert('Inconsistent data: '+'\n"Outer radius" of single-core cable: '+cable+', conductor: '+conductor+' ('+Rout+') is larger than the outermost cable radius ('+outerCableRadius+')')
            return false
         }
      }

   }
   
   //+Pipe-type cables: test the number
   if ( testNphases_pipe_cable() ) { 
      //found an error, fixed with a valid number of phases
      return false
   }

   //*Pipe-type cables: test grids
   global.Nphases_pipe_cable = phases_pipe_cable;
   // Check for inconsistent data
   if (global.Nphases_pipe_cable > 0) {

      //* Remove the tag "pipe" from the column Conductor
      for (var i=0; i<Pipe_cable_data_handle.getRowsNum(); i++){
         id_row = Pipe_cable_data_handle.getRowId(i)
         value = Pipe_cable_data_handle.cells(id_row,1).getValue(); 
         Pipe_cable_data_handle.cells(id_row,1).setValue(parseInt(value));   
      }

      //* 1. Check that each inner radius is smaller than each outer radius
      for(var row=0; row<Pipe_cable_data_handle.getRowsNum(); row++){
         idrow = Pipe_cable_data_handle.getRowId(row)  //currently tested row
         Rin = Pipe_cable_data_handle.cells(idrow,5).getValue();
         Rout = Pipe_cable_data_handle.cells(idrow,6).getValue();
         if(Rout<Rin){
            cable = Pipe_cable_data_handle.cells(idrow,0).getValue();
            conductor = Pipe_cable_data_handle.cells(idrow,1).getValue();
            Pipe_cable_data_handle.selectRowById(idrow);
            alert('Inconsistent data: '+'\n"Inner radius" is larger than "Outer radius" for: \npipe-type cable: '+cable+'\nconductor: '+conductor)
            return false
         }
      }

      //* 2. Check that all outer radius of conductors are smaller than outermost radius of cable
      for(var row=0; row<Pipe_cable_data_handle.getRowsNum(); row++){
         idrow = Pipe_cable_data_handle.getRowId(row)  //currently tested row
         cable =      Pipe_cable_data_handle.cells(idrow,0).getValue();
         conductor =  Pipe_cable_data_handle.cells(idrow,1).getValue();
         Rout =       Pipe_cable_data_handle.cells(idrow,6).getValue()*1;

         // Grab the cable outermost radius
         idrow_main_grid = Pipe_cable_main_handle.getRowId(cable-1)
         outerCableRadius = Pipe_cable_main_handle.cells(idrow_main_grid,4).getValue()*1;

         if(outerCableRadius<Rout){
            Pipe_cable_data_handle.selectRowById(idrow);
            Pipe_cable_main_handle.selectRowById(idrow_main_grid);
            alert('Inconsistent data: '+'\n"Outer radius" of pipe-type cable: '+cable+', conductor: '+conductor+' ('+Rout+') is larger than the outermost cable radius ('+outerCableRadius+')')
            return false
         }
      }

      if( !testgrids(Pipe_cable_main_handle,true) ){
         // true means that blank values are not allowed
         alert('Empty cells or errors in table "Pipe-type (PT) cable main data "!')
         return false
      }
      if( !testgrids(Pipe_cable_data_handle,true) ){
         // true means that blank values are not allowed
         alert('Empty cells or errors in table "Pipe-type (PT) cable conductors/insulators data"!')
         return false
      }
      save_Clean_Grid_Data_Table(Pipe_cable_main_handle,PipeCableMainGrid,PipeCableMainGrid_nlines)
      save_Clean_Grid_Data_Table(Pipe_cable_data_handle,PipeCableDataGrid,PipeCableDataGrid_nlines)
   }

   if(SelectedType.value=="pul_parameters"){
      if(pul_parameters_file.value == ''){
         alert('ERROR, No file selected!')
         return false
      }else{
         global.pul_parameters_file = pul_parameters_file.value;
      }
   }

   return true;
}






//*****************************************************************************************************
// Section 2. General functions, control the options displayed depending on the selection and inputs
//*************************************************************************************************** */


//+ DESCRIPTION: Changes the available options in the main table according to the selection (overhead-line/underground-cable/both)
// For example, if cables selected, will hide options related to overhead lines.
// If lines, cables or both is selected but the numer of phases is 0, hide soil and length options.
// This function is called by onload(), onchange in main selector and after changes in number of phases.
function setOptions_MainTable(){

   // ************************************************************************************************************************
   // 1. Check if "pul_parameters" selected. If yes, hide the complete main table and only show the button for selecting file.
   // ************************************************************************************************************************

   if (SelectedType.value=="pul_parameters"){
      
      document.getElementById("enter_pul_parameters_option").style.display = "";
      document.getElementById("Maintable").style.display = "none";

      if(pul_parameters_file.value != ''){
         SelectedFileDisplay.innerHTML = 'File selected: '+pul_parameters_file.value;
      }

      return

   }else{
      document.getElementById("enter_pul_parameters_option").style.display = "none";
      document.getElementById("Maintable").style.display = "";
   }

   // ************************************************************************************************************************
   // 2. Show/Hide line/cable/both options depending on main selector  
   // ************************************************************************************************************************

   if (SelectedType.value=="lines"){
      var Style1 = "";
      var Style2 = "none";
      var Style3 = "none";
   }else if(SelectedType.value=="cable"){
      var Style1 = "none";
      var Style2 = "";
      var Style3 = "none";
   }else if(SelectedType.value=="both"){
      var Style1 = "";
      var Style2 = "";
      var Style3 = "";
   }else if(SelectedType.value=="none"){
      document.getElementById("Maintable").style.display = "none";  // hide complete table
      return
   }

   // Lines
   document.getElementById("OverheadTitleR").style.display = Style1;
   document.getElementById("OverheadDatabaseR").style.display = Style1;
   document.getElementById("singlewireConducsR").style.display = Style1;
   document.getElementById("BundleConducsR").style.display = Style1;
   document.getElementById("OverheadCarR").style.display = Style1;
   document.getElementById("MidspanR").style.display = Style1;
   document.getElementById("HollowR").style.display = Style1;

   document.getElementById("line_cableBorderR").style.display = Style3;

   // Cables
   document.getElementById("CabletitleR").style.display = Style2;
   document.getElementById("SinglecoreR").style.display = Style2;
   document.getElementById("PipetypeR").style.display = Style2;
   document.getElementById("StrandedR").style.display = Style2;


   // ************************************************************************************************************************
   // 3. Check if the total number of lines and cables is zero. If yes, hide soil and line-length options. 
   // ************************************************************************************************************************

   var TotalN_conductors = phases_bundle + phases_single_wire + phases_single_core + phases_pipe_cable;
   if (TotalN_conductors==0){
      var Style="none";
   }else{
      var Style="";
   }

   // Soil 
   document.getElementById("soiltitleR").style.display = Style;
   document.getElementById("ResistivityR").style.display = Style;
   document.getElementById("RelativepermeabilityR").style.display = Style;
   document.getElementById("RelativepermittivityR").style.display = Style;
   soil_lengthBorderR.style.display = Style;
   document.getElementById("soil_cableBorderR").style.display = Style;
   MultilayerSoil_section.style.display = Style;       
   document.getElementById("MultilayerSoil_InputSection").style.display = Style;  

   // Line length 
   document.getElementById("LengthtitleR").style.display = Style;
   document.getElementById("lengthunitsR").style.display = Style;
   document.getElementById("lengthvalueR").style.display = Style;

   if (TotalN_conductors==0){  return  }
   
   // ************************************************************************************************************************
   // 4. Check if multilayer soil box is checked and set options. 
   // ************************************************************************************************************************

   if(MultilayerSoil.checked){
      var Style1 = "none"
      var Style2 = ""
   }else{
      var Style1 = ""
      var Style2 = "none"
   }

   // Uniform soil
   ResistivityR.style.display = Style1;
   RelativepermeabilityR.style.display = Style1;
   RelativepermittivityR.style.display = Style1;

   // Multilayer soil
   MultilayerSoil_InputSection.style.display = Style2;

}


// +DESCRIPTION: This function changes the displayed units (km/miles) for the line length.
//              Called from unload function and from units selector (on change)
function set_LineLengthUnits(){

   var unit_selector = document.getElementById("UnitChoice");
   var unit_label    = document.getElementById("length_km/miles_label");
   var unit_choice   = unit_selector.options[unit_selector.selectedIndex].value;

   //change unit label text
   if(unit_choice == "English"){
      unit_label.innerHTML = "miles";
   }else if (unit_choice == "Metric"){
      unit_label.innerHTML = "km";
   }
   
}
//+ Callback when changing units
function onChangingUnits(){
   reset_soil_layers_grid();
   reset_single_wire_grid();

   reset_bundle_grid();

   reset_single_core_main_grid();

   reset_pipe_cable_main_grid();

   initializeDrawing();
   set_LineLengthUnits()
}

//+


// +DESCRIPTION: Sets the list of tables with checkboxes, called by onload(), testNphases_single_wire(), 
// testNphases_bundle(), testNphases_single_core() and testNphases_pipe_cable()
function set_ListOfTables(){

   // ******************************************************************************************************************
   // 1. Before refreshing the list, rememeber which checkboxes were checked to set them to same state in the new list
   // ******************************************************************************************************************

   var OldStatus_ListofTables_section = document.getElementById("List_of_Tables_section").style.display;
   if(OldStatus_ListofTables_section == 'none'){
      var ListofTablesExistBefore = false;  // variable used to know if the list is new or existed before
   }else{
      var ListofTablesExistBefore = true;
   }

   var oldListofTablesState = document.getElementById("List_of_Tables").style.display;

   var SoilLayersCheckbox = document.getElementById("SoilLayersbox")
   if(SoilLayersCheckbox != null){  // It is null if did not exist before in the list
      OldMultilaterSoilState = document.getElementById("SoilLayersbox").checked;
   }

   var SingleWiresCheckbox = document.getElementById("Overheadsinglebox")
   if(SingleWiresCheckbox != null){  // It is null if did not exist before in the list
      OldOverheadsinglebox = document.getElementById("Overheadsinglebox").checked;
   }

   var BundleCheckbox = document.getElementById("Overheadbundlebox")
   if(BundleCheckbox != null){
      OldBundlebox = document.getElementById("Overheadbundlebox").checked;
   }

   var SingleCoreCheckbox = document.getElementById("Singlecorebox")
   if(SingleCoreCheckbox != null){
      OldSingleCorebox = document.getElementById("Singlecorebox").checked;
   }

   var PipeCheckbox = document.getElementById("Pipetypebox")
   if(PipeCheckbox != null){
      OldPipebox = document.getElementById("Pipetypebox").checked;
   }

   // ********************************************************************************
   // 2. Create code for the new list of tables 
   // ********************************************************************************

   var HTML = ''
   HTML += '<table>'+'\n'

   if(MultilayerSoil.checked){
      HTML += '    <tr>'+'\n'
      HTML += '        <label id="LabelSoil" style="cursor: pointer" title="Display the data table for the soil layers">' + '\n'
      HTML += '            <input type="checkbox" id="SoilLayersbox">Soil layers' + '\n'
      HTML += '        </label><br>' + '\n'
      HTML += '    </tr>'+'\n'
   }
   
   if(phases_single_wire>0){
      HTML += '    <tr>'+'\n'
      HTML += '        <label id="LabelSingle" style="cursor: pointer" title="Display the data table for single-wire conductors">' + '\n'
      HTML += '            <input type="checkbox" id="Overheadsinglebox">Overhead single-wire conductors' + '\n'
      HTML += '        </label><br>' + '\n'
      HTML += '    </tr>'+'\n'
   }

   if(phases_bundle>0){
      HTML += '    <tr>'+'\n'
      HTML += '        <label id="LabelBundle" style="cursor: pointer" title="Display the data table for the bundle conductors">' + '\n'
      HTML += '            <input type="checkbox" id="Overheadbundlebox">Overhead bundle conductors' + '\n'
      HTML += '        </label><br>' + '\n'
      HTML += '    </tr>'+'\n'
   }

   if(phases_single_core>0){
      HTML += '    <tr>'+'\n'
      HTML += '        <label id="LabelSinglecore" style="cursor: pointer" title="Display the data table for single-core cables">' + '\n'
      HTML += '            <input type="checkbox" id="Singlecorebox">Single-core cables' + '\n'
      HTML += '        </label><br>' + '\n'
      HTML += '    </tr>'+'\n'
   }

   if(phases_pipe_cable>0){
      HTML += '    <tr>'+'\n'
      HTML += '        <label id="LabelPipetype" style="cursor: pointer" title="Display the data table for pipe-type cables">' + '\n'
      HTML += '            <input type="checkbox" id="Pipetypebox">Pipe-type cables' + '\n'
      HTML += '        </label><br>' + '\n'
      HTML += '    </tr>'+'\n'
   }

   HTML += '</table>'+'\n'

   List_of_Tables.innerHTML = HTML

   // ***********************************************************************************
   // 3. set/reset previous status (display/hide and checked/unchecked)
   // ************************************************************************************ 

   // If at least one table exists, show the list of tables (checkboxes)
   if(phases_single_wire+phases_bundle+phases_single_core+phases_pipe_cable > 0 || MultilayerSoil.checked){
      document.getElementById("List_of_Tables_section").style.display = "block"  // there is some input, show section
      if(ListofTablesExistBefore){
         document.getElementById("List_of_Tables").style.display = oldListofTablesState; // reset old state shown/hidden
      }else{
         document.getElementById("List_of_Tables").style.display = "";   // newly created, show it
      }
   }else{
      document.getElementById("List_of_Tables_section").style.display = "none"  // there is no input, hide section
   }

   if(MultilayerSoil.checked){
      if(SoilLayersCheckbox != null){  // It is null if did not exist before in the list
         // here, it existed before, check if it changed
                  // before was unchecked    or  nlayers changed
         if(!previous_multilayerSoil_checked || previous_n_soil_layers != n_soil_layers){   //REVISAR AQUI
            //here, there was a change 
            document.getElementById("SoilLayersbox").checked = true;  // set checked
         }else{
            //here, there was no change 
            document.getElementById("SoilLayersbox").checked = OldMultilaterSoilState;  // reset previous status
         }
      }else{
         document.getElementById("SoilLayersbox").checked = true;  // newly created, set checked
      }
   }

   if(phases_single_wire>0){
      if(SingleWiresCheckbox != null){
         if(previousphases_single_wire != phases_single_wire){
            document.getElementById("Overheadsinglebox").checked = true; // there was a change, set checked
         }else{
            document.getElementById("Overheadsinglebox").checked = OldOverheadsinglebox; // no change reset old status
         }
      }else{
         document.getElementById("Overheadsinglebox").checked = true; // newly created, set checked
      }
   }

   if(phases_bundle>0){
      if(BundleCheckbox != null){
         if(previousphases_bundle != phases_bundle){
            document.getElementById("Overheadbundlebox").checked = true; // there was a change, set checked
         }else{
            document.getElementById("Overheadbundlebox").checked = OldBundlebox; // no change, reset old status
         }
      }else{
         document.getElementById("Overheadbundlebox").checked = true; // newly created, set checked
      }
   }

   if(phases_single_core>0){
      if(SingleCoreCheckbox != null){
         if(previousphases_single_core != phases_single_core){
            document.getElementById("Singlecorebox").checked = true; // there was a change, checked
         }else{
            document.getElementById("Singlecorebox").checked = OldSingleCorebox;  // no change, reset previous status
         }
      }else{
         document.getElementById("Singlecorebox").checked = true; // newly created, set checked
      }
   }

   if(phases_pipe_cable>0){
      if(PipeCheckbox != null){
         if(previousphases_pipe_cable != phases_pipe_cable){
            document.getElementById("Pipetypebox").checked = true; // there was a change, set checked
         }else{
            document.getElementById("Pipetypebox").checked = OldPipebox;  // no change, reset old status
         }
      }else{
         document.getElementById("Pipetypebox").checked = true; // newly created, set checked
      }
   }

}


//+ DESCRIPTION: This function sets the dhtmlx grids and the list of tables (checkboxes)
function set_DataTables(){

   if (MultilayerSoil.checked) {
      Set_SoilLayers_Grid();  //Initialize the SoilLayers GRID, first call
   }else{
      if(SoilLayers_handle_Set){
         SoilLayers_handle.destructor()
      }
      SoilLayersDataGrid.value = '';    // clean data
      SoilLayers_handle_Set = false;
   }

   if (phases_single_wire == 0) {
      if(Overhead_line_data_handle_Set){
         Overhead_line_data_handle.destructor()
      }
      OverheadLineDataGrid.value = '';    // clean old data (if exist)
      Overhead_line_data_handle_Set = false;
   }else{
      Set_Single_wire_data_Grid();  //Initialize the OVERHEAD-LINE GRID, first call
   }

   if (phases_bundle == 0) {
      if(Bundle_data_handle_Set){
         Bundle_data_handle.destructor()
      }
      BundleDataGrid.value = '';  //clean old data
      Bundle_data_handle_Set = false;
   }else{
      Set_bundle_data_Grid();  
   }

   if (phases_single_core == 0) {
      
      if(Single_core_main_handle_Set){
         Single_core_main_handle.destructor()
      }
      SingleCoreMainGrid.value = '';    // clean old data (if exist)
      Single_core_main_handle_Set = false;

      loadFromCableDatabaseButton.style.display = 'none'

      if(Single_core_data_handle_Set){
         Single_core_data_handle.destructor()
      }
      SingleCoreDataGrid.value = '';
      Single_core_data_handle_Set = false;

   }else{
      Set_Single_core_MAIN_Grid();  //Initialize the single-core MAIN GRID, first call
      Set_Single_core_DATA_Grid();  //Initialize the single-core DATA GRID, first call
   }

   if (phases_pipe_cable == 0) {
      if(Pipe_cable_main_handle_Set){
         Pipe_cable_main_handle.destructor()
      }
      PipeCableMainGrid.value = '';  // clean old data (if exist)
      Pipe_cable_main_handle_Set = false;

      if(Pipe_cable_data_handle_Set){
         Pipe_cable_data_handle.destructor();
      }
      PipeCableDataGrid.value = '';
      Pipe_cable_data_handle_Set = false;
   }else{
      Set_Pipe_cable_MAIN_Grid();  //Initialize the pipe-cable MAIN GRID, first call
      Set_Pipe_cable_DATA_Grid();  //Initialize the pipe-cable DATA GRID, first call
   }
}


//+ DESCRIPTION: Save the main selector state (line/cable/both/pul-parameters) when there is a change (onfocus)
//              The old state of selector is used in setOptions_MainTable()
function selectorFocus(){ 
   old_selectorState = SelectedType.value;
}


//+ DESCRIPTION: Following a change, check if data exists, if yes, a confirmation is asked because existing tables will be deleted.
//              Called form onchange in main selector
function changeInMainSelector(){

   if (SelectedType.value=="cable"){

      if (phases_single_wire+phases_bundle>0){
         var confirmation = confirm("Existing overhead line data will be deleted.\nContinue?"); 
         if (!confirmation){
            SelectedType.value = old_selectorState;
            return;
         }
         // clean existing overhead line data inputs, tables will be cleaned in function set_DataTables() called below
         Nphases_single_wire.value = '0';
         phases_single_wire = 0;           
         previousphases_single_wire = 0;
         Nphases_bundle.value = '0'
         phases_bundle = 0;                
         previousphases_bundle = 0;
      }

      pul_parameters_file.value = '';
      SelectedFileDisplay.innerHTML = '';

   }else if (SelectedType.value=="lines"){

      if (phases_pipe_cable+phases_single_core>0){
         var confirmation = confirm("Existing cable data will be deleted.\nContinue?"); 
         if (!confirmation){
            SelectedType.value = old_selectorState;
            return;
         }
         // clean existing cable data inputs, tables will be cleaned in function set_DataTables() called below
         Nphases_single_core.value = '0';
         phases_single_core = 0;           
         previousphases_single_core = 0;
         Nphases_pipe_cable.value = '0'
         phases_pipe_cable = 0;                
         previousphases_pipe_cable = 0;
      }

      pul_parameters_file.value = '';
      SelectedFileDisplay.innerHTML = '';

   }else if (SelectedType.value=="none" || SelectedType.value=="pul_parameters"){

      if(SelectedType.value=="none"){
         pul_parameters_file.value = '';
         SelectedFileDisplay.innerHTML = '';
      }

      if (phases_single_wire+phases_bundle+phases_pipe_cable+phases_single_core>0){
         var confirmation = confirm("Existing line/cable data will be deleted.\nContinue?"); 
         if (!confirmation){
            SelectedType.value = old_selectorState;
            return;
         }
         // clean existing overhead line data inputs, tables will be cleaned in function set_DataTables() called below
         Nphases_single_wire.value = '0';
         phases_single_wire = 0;           
         previousphases_single_wire = 0;
         Nphases_bundle.value = '0'
         phases_bundle = 0;                
         previousphases_bundle = 0;
         // clean existing cable data inputs, tables will be cleaned in function set_DataTables() called below
         Nphases_single_core.value = '0';
         phases_single_core = 0;           
         previousphases_single_core = 0;
         Nphases_pipe_cable.value = '0'
         phases_pipe_cable = 0;                
         previousphases_pipe_cable = 0;
      }
      
      if(MultilayerSoil.checked){
         var confirmation = confirm("Existing soil layers data will be deleted.\nContinue?"); 
         if (!confirmation){
            SelectedType.value = old_selectorState;
            return;
         }
      }

   }

   ShowHide_canvas();          // If zero inputs, canvas is hidden
   setOptions_MainTable();     
   set_ListOfTables();
   set_DataTables();
   ShowHide_dataTables();

}


//Description: This function hides the canvas when no lines/cables are selected.
function ShowHide_canvas(){

   var Total_phases = phases_single_wire + phases_bundle + phases_single_core + phases_pipe_cable;
   if(Total_phases>0){
      var Style = ""
   }else{
      var Style = "none"
   }

   document.getElementById("CanvasTable").style.display = Style;
   document.getElementById("DrawingOptions").style.display = Style;
}


//+ DESCRIPTION: This function displays/hides the tables of data when the user checks/unchecks the list of tables.
//              This function is called only when there is a change in the liste of tables (checkboxes)
function ShowHide_dataTables(){

   if(MultilayerSoil.checked){
      if(document.getElementById("SoilLayersbox").checked){
         document.getElementById("SoilLayersGrid").style.display = ""
      }else{
         document.getElementById("SoilLayersGrid").style.display = "none"
      }
   }else{
      document.getElementById("SoilLayersGrid").style.display = "none"
   }
   
   if(phases_single_wire>0){
      if(document.getElementById("Overheadsinglebox").checked){
         document.getElementById("OverheadLineGrid").style.display = ""
      }else{
         document.getElementById("OverheadLineGrid").style.display = "none"
      }
   }else{
      document.getElementById("OverheadLineGrid").style.display = "none"
   }

   if(phases_bundle>0){
      if(document.getElementById("Overheadbundlebox").checked){
         document.getElementById("BundleGrid").style.display = ""
      }else{
         document.getElementById("BundleGrid").style.display = "none"
      }
   }else{
      document.getElementById("BundleGrid").style.display = "none"
   }

   if(phases_single_core>0){
      if(document.getElementById("Singlecorebox").checked){
         document.getElementById("MainGridSingleCore").style.display = ""
         document.getElementById("CondInsu_SingleCore_grid").style.display = ""
         document.getElementById("loadFromCableDatabaseButton").style.display = ""
      }else{
         document.getElementById("MainGridSingleCore").style.display = "none"
         document.getElementById("CondInsu_SingleCore_grid").style.display = "none"
         document.getElementById("loadFromCableDatabaseButton").style.display = "none"
      }
   }else{
      document.getElementById("MainGridSingleCore").style.display = "none"
      document.getElementById("CondInsu_SingleCore_grid").style.display = "none"
      document.getElementById("loadFromCableDatabaseButton").style.display = "none"
   }

   if(phases_pipe_cable>0){
      if(document.getElementById("Pipetypebox").checked){
         document.getElementById("MainGridPipeCable").style.display = ""
         document.getElementById("CondInsu_PipeCable_grid").style.display = ""
      }else{
         document.getElementById("MainGridPipeCable").style.display = "none"
         document.getElementById("CondInsu_PipeCable_grid").style.display = "none"
      }
   }else{
      document.getElementById("MainGridPipeCable").style.display = "none"
      document.getElementById("CondInsu_PipeCable_grid").style.display = "none"
   }

}


// DESCRIPTION: This function opens a dialog box to select a file that contains the p.u.l. parameters 
function fileselect(){

   var isfileok = mydatafile.openDialog('Select p.u.l. parameters file','mat',designfilepath);

   if(isfileok){
      //project_folder_name.value = mydatafile.fullName(); //global
      SelectedFileDisplay.innerHTML = 'File selected: '+mydatafile;
   }

   pul_parameters_file.value = mydatafile;

}








//**************************************************************************************************
//  Section 3. Functions related to soil options
//************************************************************************************************ */


// DESCRIPTION: Test the resistivity value entered for homogeneous soil
function test_resistivity() {
   var retError = testposrealvalue(Resistivity,old_soil_resistivity);
   if(!retError){ old_soil_resistivity = Resistivity.value; }
   return retError
}


// DESCRIPTION: Test the value given to relative permeability for homogeneous soil
function test_mu() {
   var retError = testposrealvalue(soil_mu,old_soil_permeability);
   if(!retError){ old_soil_permeability = soil_mu.value; }
   return retError
}


// DESCRIPTION: Test the value given to relative permittivity for homogeneous soil
function test_epsilon() {
   var retError = testposrealvalue(soil_epsilon,old_soil_permittivity);
   if(!retError){ old_soil_permittivity = soil_epsilon.value; }
   return retError
}


//+ DESCRIPTION: Test the number of soil layers entered
// The minimum value accepted is 2 and maximum 4. For 1 layer, the multilayer option should be unchecked.
function test_NlayersSoil() {

       //*no #..# usage also rejects ##                minimum 2 layers                 limit to 4 layers
   if( !notspecialopt("SoilLayers") || parseInt(SoilLayers.value)<2 || parseInt(SoilLayers.value)>4 ){    
      SoilLayers.value = previous_n_soil_layers;  // reset to previous value
      return
   }
   //n_soil_layers is global
   var retError = testposintvalue(SoilLayers,n_soil_layers+''); 
   if(!retError){

      // The entered value was accepted, reset options and tables for the new value
      n_soil_layers = parseInt(SoilLayers.value);

      setOptions_MainTable();    // Set the options in the main table, not needed f 

      set_ListOfTables();  // Set the list of tables containing data (list of checkboxes) 
      //Must do it since a new checkbox must be added

      set_DataTables();   // Set the data tables, all DHTMLX grids
      //Here Set_bundle_data_Grid is called, TODO: no need!

      ShowHide_dataTables();        
      // show/hide data tables according to the status on the list of tables (checked/unchecked)
      // No need to do all tables

      ShowHide_canvas();   // May hide canvas in case all phases are equal to 0

      // update the value of previous number of layers
      previous_n_soil_layers = n_soil_layers;

   }

   return

}



// +DESCRIPTION: This function sets the number of rows/cables of the table containing the data of single-wire conductors
function Set_SoilLayers_Grid() {
   
   if(!SoilLayers_handle_Set){
      
      var width_columns = '115,100,130,130,130';
      var data_properties = 'None,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

      // Define labels 
      if(length_units.value == 'Metric'){
         GridLabels = "Layer,Thickness (m),Resistivity (Ohm m),Relative permittivity,Relative permeability"
      }else{         
         GridLabels = "Layer,Thickness (ft),Resistivity (Ohm mile),Relative permittivity,Relative permeability"
      }   
      
      SoilLayers_handle = CreateNColumnGrids(5,0,
                                    SoilLayers_dataGridPlace,
                                    SoilLayers_data_column,      // where we also put title
                                    SoilLayersDataGrid,          // grid content
                                    GridLabels,
                                    width_columns,
                                    data_properties);

      SoilLayers_handle.setUserData('','BlockedKeys','Insert','CtrlDelete');
   
      SoilLayers_handle_Set = true;
   }
   

   var NumberOfRows_required = n_soil_layers;
   CurrentNumberOfRows = SoilLayers_handle.getRowsNum()
       
   if( CurrentNumberOfRows != NumberOfRows_required ) {
      //Remove rows or add rows, puts 0 in added data
      resetDataGridRows(SoilLayers_handle,NumberOfRows_required,'')
   }   
   
}


// +DESCRIPTION: This function resets the rows/columns of the data table for soil layers.
// It is called when the units (Metric/English are changed)
function reset_soil_layers_grid() {
   if (MultilayerSoil.checked) {
      SoilLayers_handle_Set = false;
      SoilLayersDataGrid.value = SoilLayers_handle.serializeToCSV();
      Set_SoilLayers_Grid();
   }
}







//**************************************************************************************************
//  Section 4. Functions related to overhead line with single wires
//************************************************************************************************ */

// DESCRIPTION: This function sets the number of rows/cables of the table containing the data of single-wire conductors
function Set_Single_wire_data_Grid() {
   
   if(!Overhead_line_data_handle_Set){
      
      var metric_units = 'Metric'.localeCompare(length_units.value);
      var cond_Rdc = 'DC_resistance'.localeCompare(Conductor_characteristic.value);
      var Ncolumns_overhead = '';
      var width_columns = '';
      var data_prop = '';

      // Define labels 
      if(Midspan.checked){

         if(Hollow.checked){ // both Midspan and hollow checked

            Ncolumns_overhead = 10;     
            width_columns = '115,100,130,130,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){ 
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (m),Tower height (m),Midspan height (m),Inner radius (cm),Outer radius (cm),DC resistance (Ohm/km),Conductor relative permeability,Conductor relative permittivity"
               }else{         GridLabels = "Conductor,Phase,Horizontal position (m),Tower height (m),Midspan height (m),Inner radius (cm),Outer radius (cm),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"     }
            }else{
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),Inner radius (in),Outer radius (in),DC resistance (Ohm/mile),Conductor relative permeability,Conductor relative permittivity"
               }else{         GridLabels = "Conductor,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),Inner radius (in),Outer radius (in),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"  }
            }

         }else{ // midspan checked and hollow unchecked

            Ncolumns_overhead = 9;     
            width_columns = '115,100,130,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){ 
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (m),Tower height (m),Midspan height (m),Radius (cm),DC resistance (Ohm/km),Conductor relative permeability,Conductor relative permittivity"
               }else{         GridLabels = "Conductor,Phase,Horizontal position (m),Tower height (m),Midspan height (m),Radius (cm),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"     }
            }else{
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),Radius (in),DC resistance (Ohm/mile),Conductor relative permeability,Conductor relative permittivity"
               }else{         GridLabels = "Conductor,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),Radius (in),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"   }
            }
         }

      }else{   

         if(Hollow.checked){ // midspan unchecked and hollow checked
            
            Ncolumns_overhead = 9;     
            width_columns = '115,100,130,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){ 
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (m),Height (m),Inner radius (cm),Outer radius (cm),DC resistance (Ohm/km),Conductor relative permeability,Conductor relative permittivity"
               }else{         GridLabels = "Conductor,Phase,Horizontal position (m),Height (m),Inner radius (cm),Outer radius (cm),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"     }
            }else{
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (ft),Height (ft),Inner radius (in),Outer radius (in),DC resistance (Ohm/mile),Conductor relative permeability,Conductor relative permittivity"
               }else{         GridLabels = "Conductor,Phase,Horizontal position (ft),Height (ft),Inner radius (in),Outer radius (in),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"   }
            }

         }else{ // Both midspan and hollow unchecked

            Ncolumns_overhead = 8;    
            width_columns = '115,100,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (m),Height (m),Radius (cm),DC resistance (Ohm/km),Conductor relative permeability,Conductor relative permittivity"
               }else{         GridLabels = "Conductor,Phase,Horizontal position (m),Height (m),Radius (cm),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"    }
            }else{
               if(!cond_Rdc){ GridLabels = "Conductor,Phase,Horizontal position (ft),Height (ft),Radius (in),DC resistance (Ohm/mile),Conductor relative permeability,Conductor relative permittivity"
               }else{        GridLabels = "Conductor,Phase,Horizontal position (ft),Height (ft),Radius (in),Resistivity (Ohm m),Conductor relative permeability,Conductor relative permittivity"  }
            }
         }
      }
      
      Overhead_line_data_handle = CreateNColumnGrids(Ncolumns_overhead,0,
                                    OverheadLine_dataGridPlace,
                                    OverheadLine_data_column,      //where we also put title
                                    OverheadLineDataGrid,          // grid content
                                    GridLabels,
                                    width_columns,
                                    data_prop);

      Overhead_line_data_handle.setUserData('','BlockedKeys','Insert','CtrlDelete');
   
      Overhead_line_data_handle_Set = true;
   }
   
   var NumberOfRows = phases_single_wire;
   CurrentNumberOfRows = Overhead_line_data_handle.getRowsNum()
       
   if( CurrentNumberOfRows != NumberOfRows ) {
      //Remove rows or add rows, puts 0 in added data
      resetDataGridRows(Overhead_line_data_handle,NumberOfRows,'')
   }   

   // By default set permeability and permittivity values to 1, these values can be changed by user
   var Nrows = Overhead_line_data_handle.getRowsNum();
   var Ncols = Overhead_line_data_handle.getColumnsNum();

   if(Ncols == 8){
      var permeability_column = 6;
      var permittivity_column = 7;
   }else if(Ncols == 9){
      var permeability_column = 7;
      var permittivity_column = 8;
   }else if(Ncols == 10){
      var permeability_column = 8;
      var permittivity_column = 9;
   }

   for(var row=0; row<Nrows; row++){
      var rowId = Overhead_line_data_handle.getRowId(row)
      var ur = Overhead_line_data_handle.cells(rowId,permeability_column).getValue();
      if(ur == ""){   // only if the cell is empty, 1 is set 
         Overhead_line_data_handle.cells(rowId,permeability_column).setValue(1); 
      }
      var er = Overhead_line_data_handle.cells(rowId,permittivity_column).getValue();
      if(er == ""){   // only if the cell is empty, 1 is set 
         Overhead_line_data_handle.cells(rowId,permittivity_column).setValue(1); 
      }
   }
   
}


// +DESCRIPTION: Test the number of single-wire conductors entered.
//              Called when there is a change by entering a number or by pushing the button +1 or -1
function testNphases_single_wire() {
  
          //*no #..# usage also rejects ##                limit to 50 conductors
   if( !notspecialopt("Nphases_single_wire") || parseInt(Nphases_single_wire.value)> 150 ){    
      Nphases_single_wire.value = phases_single_wire;
      return
   }

   var retError = testpos0intvalue(Nphases_single_wire,phases_single_wire+''); 

   if(!retError){

      // The entered value was accepted, reset options and tables for the new value
      phases_single_wire = parseInt(Nphases_single_wire.value);

      setOptions_MainTable();       // Set the options in the main table 
      set_ListOfTables();           // Set the list of tables containing data (list of checkboxes)               
      set_DataTables();             // Set the data tables
      ShowHide_dataTables();        // show/hide data tables according to the status on the list of tables (checked/unchecked)
      ShowHide_canvas();            // May hide canvas in case all phases are equal to 0

      // update the value of previous phases
      previousphases_single_wire = phases_single_wire

   }

   return 

} 




//+ DESCRIPTION: This function sets the rows/columns of the data table for single-wire conductors.
function reset_single_wire_grid() {
   if (phases_single_wire > 0) {
      Overhead_line_data_handle_Set = false;
      OverheadLineDataGrid.value = Overhead_line_data_handle.serializeToCSV();
      Set_Single_wire_data_Grid();
   }
}

//+ Callback when Multilayer soil is checked or unchecked
function onMultilayerSoilChange(){
   setOptions_MainTable(); //shows all selections in the main table
   set_ListOfTables();     //onMultilayerSoilChange
   set_DataTables();
   ShowHide_dataTables();
   initializeDrawing()
}




//**************************************************************************************************
// Section 5. Functions related to overhead line with bundle of conductors
//************************************************************************************************ */

// +DESCRIPTION: This function sets the number of rows/cables of the table containing the data of bundle conductors
function Set_bundle_data_Grid() {  

   if(!Bundle_data_handle_Set){
      
      var metric_units = 'Metric'.localeCompare(length_units.value);    
      var cond_Rdc = 'DC_resistance'.localeCompare(Conductor_characteristic.value);  
      //the above will give 0 when DC Resistance is selected, strange coding
      var Ncolumns_bundle = '';
      var width_columns = '';
      var data_prop = '';

      // Define labels 
      if(Midspan.checked){

         if(Hollow.checked){ // both Midspan and hollow checked

            Ncolumns_bundle = 13; //Midspan an Hollow    
            width_columns = '100,100,130,130,130,130,130,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,IntPosNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){ 
               if(!cond_Rdc){ //+ with DC resistance data, METRIC, MIDSPAN, HOLLOW
                  GridLabels = "Bundle,Phase,Horizontal position (m),Tower height (m),Midspan height (m),"+
                               "Number of conductors,Angle (deg),"+
                               "Bundle radius (cm),Inner radius (cm),Outer radius (cm),DC resistance (Ohm/km),"+
                               "Conductor relative permeability,Conductor relative permittivity"
               }else{         //+With Resistivity, METRIC, MIDSPAN, HOLLOW     
                  GridLabels = "Bundle,Phase,Horizontal position (m),Tower height (m),Midspan height (m),"+
                               "Number of conductors,Angle (deg),"+
                               "Bundle radius (cm),Inner radius (cm),Outer radius (cm),Resistivity (Ohm m),"+
                               "Conductor relative permeability,Conductor relative permittivity"     }
            }else{
               if(!cond_Rdc){ //+ with DC Resistance, English, MIDSPAN, HOLLOW
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),"+
                               "Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Inner radius (in),Outer radius (in),DC resistance (Ohm/mile),"+
                               "Conductor relative permeability,Conductor relative permittivity"
               }else{        //+ with Resistivity, ENGLISH, MIDSPAN, HOLLOW
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),"+
                               "Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Inner radius (in),Outer radius (in),Resistivity (Ohm m),"+
                               "Conductor relative permeability,Conductor relative permittivity"   }
            }

         }else{ // midspan checked and hollow unchecked

            Ncolumns_bundle = 12;   //just midspan  
            width_columns = '100,100,130,130,130,130,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,IntPosNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){ 
               if(!cond_Rdc){ //+ with DC resistance, METRIC, Midspan
                      GridLabels = 
                      "Bundle,Phase,Horizontal position (m),Tower height (m),Midspan height (m),"+
                      "Number of conductors,Angle (deg),"+
                      "Bundle radius (cm),Conductor radius (cm),DC resistance (Ohm/km),"+
                      "Conductor relative permeability,Conductor relative permittivity"
               }else{   //+ with Resistivity, METRIC, Midspan      
                      GridLabels = 
                      "Bundle,Phase,Horizontal position (m),Tower height (m),Midspan height (m),"+
                      "Number of conductors,Angle (deg),"+
                      "Bundle radius (cm),Conductor radius (cm),Resistivity (Ohm m),"+
                      "Conductor relative permeability,Conductor relative permittivity"     }
            }else{
               if(!cond_Rdc){ //+ with DC resistance, english, Midspan
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),"+
                               "Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Conductor radius (in),DC resistance (Ohm/mile),"+
                               "Conductor relative permeability,Conductor relative permittivity"
               }else{        //+ with Resistivity, ENGLISH, Midspan   
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Tower height (ft),Midspan height (ft),"+
                               "Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Conductor radius (in),Resistivity (Ohm m),"+
                               "Conductor relative permeability,Conductor relative permittivity"   }
            }
         }

      }else{   

         if(Hollow.checked){ // midspan unchecked and hollow checked
            
            Ncolumns_bundle = 12;   //Hollow  
            width_columns = '100,100,130,130,130,130,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,IntPosNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){ 
               if(!cond_Rdc){ //+ DC resistance, Hollow, Metric
                  GridLabels = "Bundle,Phase,Horizontal position (m),Height (m),Number of conductors,Angle (deg),"+
                               "Bundle radius (cm),Inner radius (cm),Outer radius (cm),DC resistance (Ohm/km),"+  //+DC R
                               "Conductor relative permeability,Conductor relative permittivity"
               }else{        //+Resistivity, Hollow, Metric  
                  GridLabels = "Bundle,Phase,Horizontal position (m),Height (m),Number of conductors,Angle (deg),"+
                               "Bundle radius (cm),Inner radius (cm),Outer radius (cm),Resistivity (Ohm m),"+     //+ Resistivity
                               "Conductor relative permeability,Conductor relative permittivity"     
               }
            }else{
               if(!cond_Rdc){ //+ DC resistance, Hollow, English 
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Height (ft),Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Inner radius (in),Outer radius (in),DC resistance (Ohm/mile),"+  //+DC R
                               "Conductor relative permeability,Conductor relative permittivity"
               }else{         //+Resistivity, Hollow, ENglish
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Height (ft),Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Inner radius (in),Outer radius (in),Resistivity (Ohm m),"+       //+ Resistivity
                               "Conductor relative permeability,Conductor relative permittivity"   }
            }

         }else{ // Both midspan and hollow unchecked

            Ncolumns_bundle = 11;    //no hollow no midspan
            width_columns = '100,100,130,130,130,130,130,130,130,130,130';
            data_prop = 'None,IntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,IntPosNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

            if(!metric_units){
               if(!cond_Rdc){   //+ DC resistance, metric
                  GridLabels = "Bundle,Phase,Horizontal position (m),Height (m),Number of conductors,Angle (deg),"+
                               "Bundle radius (cm),Conductor radius (cm),DC resistance (Ohm/km),"+
                               "Conductor relative permeability,Conductor relative permittivity"
               }else{           //+ Resistivty, metric
                  GridLabels = "Bundle,Phase,Horizontal position (m),Height (m),Number of conductors,Angle (deg),"+
                               "Bundle radius (cm),Conductor radius (cm),Resistivity (Ohm m),"+
                               "Conductor relative permeability,Conductor relative permittivity"    }
            }else{
               if(!cond_Rdc){  //+ DC resistance, English
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Height (ft),Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Conductor radius (in),DC resistance (Ohm/mile),"+
                               "Conductor relative permeability,Conductor relative permittivity"
               }else{          //+ Resistiviy, English
                  GridLabels = "Bundle,Phase,Horizontal position (ft),Height (ft),Number of conductors,Angle (deg),"+
                               "Bundle radius (in),Conductor radius (in),Resistivity (Ohm m),"+
                               "Conductor relative permeability,Conductor relative permittivity"  
               }
            }
         }
      }   
      Bundle_data_handle = CreateNColumnGrids(Ncolumns_bundle,0,
                                    Bundle_dataGridPlace,
                                    Bundle_data_column,      //where we also put title
                                    BundleDataGrid,          // grid content
                                    GridLabels,
                                    width_columns,
                                    data_prop);

      Bundle_data_handle.setUserData('','BlockedKeys','Insert','CtrlDelete');
   
      Bundle_data_handle_Set = true;
   }
   

   //*We jump here if the grid is already available
   var NumberOfRows = phases_bundle;  //*Bundled conductors, number
   CurrentNumberOfRows = Bundle_data_handle.getRowsNum()
       
   if( CurrentNumberOfRows != NumberOfRows ) {
      //Remove rows or add rows, puts 0 in added data
      resetDataGridRows(Bundle_data_handle,NumberOfRows,'')
   }   

   // By default set permeability and permittivity values to 1, these values can be changed by user
   var Nrows = Bundle_data_handle.getRowsNum();
   var Ncols = Bundle_data_handle.getColumnsNum();
   var permeability_column = Ncols-2
   var permittivity_column = Ncols-1

   for(var row=0; row<Nrows; row++){
      var rowId = Bundle_data_handle.getRowId(row)

      var ur = Bundle_data_handle.cells(rowId,permeability_column).getValue();
      if(ur == ""){   // only if the cell is empty, 1 is set 
         Bundle_data_handle.cells(rowId,permeability_column).setValue(1); 
      }
      var er = Bundle_data_handle.cells(rowId,permittivity_column).getValue();
      if(er == ""){   // only if the cell is empty, 1 is set 
         Bundle_data_handle.cells(rowId,permittivity_column).setValue(1); 
      }
   }
   
}


// +DESCRIPTION: Test the number of bundle conductors entered
//              Called when there is a change by entering a number or by pushing the button +1 or -1
function testNphases_bundle() {

      //*no #..# usage also rejects ##                limit to 50 conductors
   if( !notspecialopt("Nphases_bundle") || parseInt(Nphases_bundle.value)>150 ){    
      Nphases_bundle.value = phases_bundle;
      //phases_bundle is global
      return
   }

   var retError = testpos0intvalue(Nphases_bundle,phases_bundle+''); 

   if(!retError){

      // The entered value was accepted, reset options and tables for the new value
      phases_bundle = parseInt(Nphases_bundle.value);

      setOptions_MainTable();       // Set the options in the main table 
      set_ListOfTables();           // Set the list of tables containing data (list of checkboxes)               
      set_DataTables();             // Set the data tables
      ShowHide_dataTables();        // show/hide data tables according to the status on the list of tables (checked/unchecked)
      ShowHide_canvas();            // May hide canvas in case all phases are equal to 0

      // update the value of previous phases
      previousphases_bundle = phases_bundle

   }

   return

} 





//***********************************************************************************************
// Section 6. Aditional options for overhead lines (database, and consider midspan and hollow)
//********************************************************************************************* */

// DESCRIPTION: This function adds/removes columns of single-wire and bundle conductors tables when the Midspan checkbox is checked/unchecked
function modify_midspan() {

   if (Hollow.checked) { 
      if (Midspan.checked){ // Hollow was checked and Midspan just checked

         if (phases_single_wire > 0) { // add column to single-wires table
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) {
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }
               }
            }
         }

         if (phases_bundle > 0) { // add column to bundles table
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) {
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }
               }
            }
         }

      }else{ // Hollow was checked and Midspan just unchecked

         if (phases_single_wire > 0) { // remove column to single-wires table
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) {
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\n';
                  }
               }
            }
         }

         if (phases_bundle > 0) { // remove column to bundles table
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) {
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\t'+ row_n[11] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\t'+ row_n[11] +'\n';
                  }
               }
            }
         }

      }

   }else{
      if (Midspan.checked){ // Hollow was unchecked and Midspan just checked

         if (phases_single_wire > 0) { // add column to single-wires table
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) {
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\n';
                  }
               }
            }
         }

         if (phases_bundle > 0) { // add column to bundles table
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) {
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\n';
                  }
               }
            }
         }

      }else{ // Hollow was unchecked and Midspan just unchecked

         if (phases_single_wire > 0) { // remove column to single-wires table
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) {
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }
               }
            }
         }

         if (phases_bundle > 0) { // remove column to bundles table
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) {
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }
               }
            }
         }
      }
   }

   Overhead_line_data_handle_Set = false;
   if (phases_single_wire > 0) {
      OverheadLineDataGrid.value = new_table_single_wires;
   }
   Set_Single_wire_data_Grid();

   Bundle_data_handle_Set = false;
   if (phases_bundle > 0) {
      BundleDataGrid.value = new_table_bundles;
   }
   Set_bundle_data_Grid();
   
}


// DESCRIPTION: This function adds/removes columns of single-wire and bundle conductors tables when the Hollow checkbox is checked/unchecked
function modify_hollow() {
   
   if (Midspan.checked) {

      if (Hollow.checked){ // Midspan was already checked and hollow just checked

         if (phases_single_wire > 0) {
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) {  // add column to single wires table
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ '' +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ '' +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }
               }
            }
         }

         if (phases_bundle > 0) {
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) {  // add column to bundle conductors table
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ '' + '\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ '' +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }
               }
            }
         }

      }else{ // Midspan was checked and hollow just unchecked

         if (phases_single_wire > 0) {
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) { // remove column to single wires table
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[5] +'\t'+ row_n[6] + '\t'+ row_n[7] +'\t'+ row_n[8] +'\n';
                  }
               }
            }
         }

         if (phases_bundle > 0) {
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) { // remove column to bundle conductors table
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\t'+ row_n[11] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\t'+ row_n[11] +'\n';
                  }
               }
            }
         }

      }

   }else{ // Midspan unchecked

      if (Hollow.checked){ // Midspan was unchecked and hollow just checked

         if (phases_single_wire > 0) {
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) {  // add column to single wires table
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ '' +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\n';
                  }
               }
            }
         }
         
         if (phases_bundle > 0) {
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) {  // add column to bundles table
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+  '' +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ '' +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\n';
                  }
               }
            }
         }

      }else{ //Midspan was unchecked and hollow just unchecked

         if (phases_single_wire > 0) {
            var new_table_single_wires = '';
            for (n=0; n<phases_single_wire; n++) { // remove column to single wires table
               row_n = Overhead_line_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_single_wires = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }else{
                     new_table_single_wires = new_table_single_wires + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\n';
                  }
               }
            }
         }
         
         if (phases_bundle > 0) {
            var new_table_bundles = '';
            for (n=0; n<phases_bundle; n++) { // remove column to single bundles table
               row_n = Bundle_data_handle.serializeToCSV();
               if(row_n.length>0){
                  row_n = row_n.split('\n');
                  row_n = row_n[n].split('\t');
                  if (n==0){
                     new_table_bundles = row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }else{
                     new_table_bundles = new_table_bundles + row_n[0] +'\t'+ row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] +'\n';
                  }
               }
            }
         }
      }

   }

   Overhead_line_data_handle_Set = false;
   if (phases_single_wire > 0) {
      OverheadLineDataGrid.value = new_table_single_wires;
   }
   Set_Single_wire_data_Grid();

   Bundle_data_handle_Set = false;
   if (phases_bundle > 0) {
      BundleDataGrid.value = new_table_bundles;
   }
   Set_bundle_data_Grid();

}


//* DESCRIPTION: Callback for the button 'Use Overhead Line Database'
// This function first creates the line data base object, then opens the corresponding database mask.
// When the user closes the database mask, the database file path and the line id selected are given as 
// return. These return variables are used to retrieve the information from the selected database .xml file
// and populate the tables and options in the mask of LineCable_Data.
function loadFromDatabaseButtonClick(){
   
   // read field units from LineCable_Data mask, used to set the units in the database mask
   var units = document.getElementById("UnitChoice").value; 

   // Create the database object
   var emtpdb = new emtpLineDatabase();
   emtpdb.ACTIVE_UNIT_SYSTEM = (units=="English") ? "English" : "Metric";
   
   // Put the database object in global values (myglobal)
   myglobal.emtpdb = emtpdb;
   // Get the list of databases available, "Line" refers to the overhead line database (it is given as default)
   myglobal.EMTPDB_FILELIST = emtpdb.listOfDatabaseXMLFile(myglobal, "Line")


   // ================= Open the database mask ===================================================================

   var ret = emtpdb.loadLineFromDatabase(myglobal); 
   //  ret contains the information about the selected line from the database

   if(ret == null) return ;// no database selected
   
   var database_filepath = ret[0];  // database file (could be the default or an user-created database)
   var            lineid = ret[1];  // id (name) of the selected line

   // ============================================================================================================


   // Before modifying, emtpy input boxes for all types of lines/cables in LineCable_Data
   document.getElementById("NbSingleWire").value = "0";
   document.getElementById("NbBundle").value = "0";
   document.getElementById("NbSinglecore").value = "0";
   document.getElementById("NbPipetype").value = "0";

   var distanceUnit = (units=="English") ? "ft":"m" ;
   var diameterUnit = (units=="English") ? "in":"cm" ;
   var bundleSpacingUnit = (units=="English") ? "in":"cm" ;
   var DCresistanceUnit = (units=="English") ? "Ohm/mi":"Ohm/km" ;
   //midspan and hollow checkboxes state
   var midspan = document.getElementById('Midspanbox').checked;
   var  hollow = document.getElementById('Hollowbox').checked;

   
   //1. Get the details of the selected line 
   var conductorIndexes = emtpdb.getConductorIndexes(database_filepath,lineid);

   //1.1 Update phase number in the mask
   nphases = conductorIndexes.length;

   //seperate index for both grids
   var bindex = -1;   // index for bundle conductors grid
   var oindex = -1;   // index for single-wire conductors grid

   // 2.0 Update parameters loop
   for (var index=0; index<conductorIndexes.length; index++){ // rows
      //adjust number of bundle or overhead
      var nbSC = document.getElementById("NbSingleWire").value;
      var nbBC = document.getElementById("NbBundle").value;
      var Grid;
      var cell;
      //chosen index
      var uindex;
      //get data
      var conductorIndex=conductorIndexes[index];
      var conductorID=emtpdb.getHangingPointParameter(database_filepath,lineid,conductorIndex,"ConductorID");
      var phaseNumber=(conductorIndex.charAt(0)=='P') ? conductorIndex.substr(1) : '0'; // If ground conductor put the phase number to 0    
      //check if conductor is a bundle or overheadline
      if (phaseNumber =="0" ){ // Ground conductor
         var bundled=emtpdb.getGroundConductorParameterUnit(database_filepath,conductorID,"Bundled");
         var DCresistance=emtpdb.getGroundConductorParameterUnit(database_filepath,conductorID,"DCResistance",DCresistanceUnit);
         var diameter=emtpdb.getGroundConductorParameterUnit(database_filepath,conductorID,"OuterDiameter",diameterUnit);
      }else{//not ground conductor
         var bundled=emtpdb.getPhaseConductorParameterUnit(database_filepath,conductorID,"Bundled");		
         var DCresistance=emtpdb.getPhaseConductorParameterUnit(database_filepath,conductorID,"DCResistance",DCresistanceUnit);
         var diameter=emtpdb.getPhaseConductorParameterUnit(database_filepath,conductorID,"OuterDiameter",diameterUnit)	;
      }
      //distance variables
      var x=emtpdb.getHangingPointParameterUnit(database_filepath,lineid,conductorIndex,"HorizontalDistance",distanceUnit);
      var y=emtpdb.getHangingPointParameterUnit(database_filepath,lineid,conductorIndex,"ConductorHeight",distanceUnit);
      //select appropriate grid
      if (bundled=="Yes"){
         
         //adjust number of bundle input
         if(nbBC==""){nbBC=0};
         //add to the number of bundled conductors
         nbBC=parseInt(nbBC);
         nbBC++;
         nbBC=nbBC+"";
         //update the HTML
         document.getElementById("NbBundle").value=nbBC;
         //update grid
         testNphases_bundle();
         Grid = Bundle_data_handle;
         bindex++
         //get apropriate index
         uindex = bindex;
      }else{
         //adjust number of singlecore input
         if(nbSC==""){nbSC=0};
         //add to the number of singlecore
         nbSC=parseInt(nbSC);
         nbSC++;
         nbSC=nbSC+"";
         //update the HTML
         document.getElementById("NbSingleWire").value=nbSC;
         //update grid
         testNphases_single_wire();
         Grid = Overhead_line_data_handle;
         oindex++
         //get appropriate index
         uindex = oindex;
      }
   
      //hide table if needed
      ShowHide_dataTables();
      //phase make this part look better;
      cell = 1;
      var phaseNumber=(conductorIndex.charAt(0)=='P') ? conductorIndex.substr(1) : '0'; // If ground conductor put the phase number to 0    
      Grid.cellByIndex(uindex,cell).setValue(phaseNumber);
      /*
      //FOR BUNDLES
      if(bundled=="Yes"){
         var numCol=10;
         if(midspan){numCol++};
         if(hollow){numCol++};
         for (var col=1;col<numCol;col++){

         }
      }
      */
      //DC resistance
      bundled=="Yes"?cell = 8:cell=5;
      if (midspan){cell++};
      if (hollow){cell++};
      Grid.cellByIndex(uindex,cell).setValue(DCresistance);
      //Outside diameter
      bundled=="Yes"?cell=7:cell=4;
      if(midspan){cell++};
      if(hollow){cell++};
      Grid.cellByIndex(uindex,cell).setValue(diameter)
      //distances
      cell=2;
      Grid.cellByIndex(uindex,cell).setValue(x);
      cell=3;
      Grid.cellByIndex(uindex,cell).setValue(y);
      //Set Conductor relative permeability and permitivity to 1
      bundled=="Yes"?cell = 9:cell=6;
      if(midspan){cell++};
      if(hollow){cell++};
      Grid.cellByIndex(uindex,cell).setValue(1);
      cell++;
      Grid.cellByIndex(uindex,cell).setValue(1);
      
      //bundles 
      if(bundled=="Yes"){
         if (phaseNumber =="0" ){ // Ground conductor
            var number= emtpdb.getGroundConductorParameter(database_filepath,conductorID,"Number")
            var spacing=emtpdb.getGroundConductorParameterUnit(database_filepath,conductorID,"Spacing",bundleSpacingUnit)
            var angularPotion=emtpdb.getGroundConductorParameter(database_filepath,conductorID,"Angle")
         }else{
            var number= emtpdb.getPhaseConductorParameter(database_filepath,conductorID,"Number")
            var spacing=emtpdb.getPhaseConductorParameterUnit(database_filepath,conductorID,"Spacing",bundleSpacingUnit)
            var angularPotion=emtpdb.getPhaseConductorParameter(database_filepath,conductorID,"Angle")
         }   
         //formula spacing to bundle radius
         var BundleRadius = Math.round(parseFloat(spacing)/(2*Math.sin(Math.PI/parseInt(number)))*1000)/1000
         //coerce
         BundleRadius += "";
         //fill table
         cell=4;
         if(midspan){cell++};
         Grid.cellByIndex(uindex,cell).setValue(number);
         cell=5;
         if(midspan){cell++};
         Grid.cellByIndex(uindex,cell).setValue(angularPotion);
         cell=6;
         if(midspan){cell++};
         Grid.cellByIndex(uindex,cell).setValue(BundleRadius)  
      }
   }

   // Update the drawing in LineCable_Data
   initializeDrawing();
   
   // Verify that all number of conductors are valid
   testNphases_bundle();
   testNphases_single_wire();
   testNphases_pipe_cable();
   testNphases_single_core();

   // Show/hide tables according to selections
   ShowHide_dataTables();  

}


// +DESCRIPTION: This function sets the rows/columns of the data table for single-wire conductors.
function reset_bundle_grid() {
   if (phases_bundle > 0) {
      Bundle_data_handle_Set = false;
      BundleDataGrid.value = Bundle_data_handle.serializeToCSV();
      Set_bundle_data_Grid();
   }
}





// **************************************************************************************************
//  Section 7. Functions related to single-core cables
// ************************************************************************************************ */

// DESCRIPTION: This function sets the number of rows/cables of the table containing the main data of single-core cables
function Set_Single_core_MAIN_Grid() {  

   if(!Single_core_main_handle_Set){
      
      var metric_units = 'Metric'.localeCompare(length_units.value);
      
      if(!metric_units){ 
         GridLabels = "Cable,Number of conductors,Horizontal position (m),Ground depth (m),Radius (cm)"
      }else{
         GridLabels = "Cable,Number of conductors,Horizontal position (ft),Ground depth (ft),Radius (in)"
      }

      num_char = 'None,IntPosNotNamed,RealOrIntNotNamed,RealOrIntNotNamed,RealOrIntPosNotNamed';

      Single_core_main_handle = CreateNColumnGrids(5,0,
                                 Single_core_mainGridPlace,
                                 Single_core_main_grid,      // where we also put title
                                 SingleCoreMainGrid,         // grid content
                                 GridLabels,
                                 '100,120,130,130,130',
                                 num_char);

      Single_core_main_handle.setUserData('','BlockedKeys','Insert','CtrlDelete');

      // Shows the single-core database button when a row is selected in the single-core main table
      Single_core_main_handle.attachEvent("onRowSelect",function(){ loadFromCableDatabaseButton.style.display='block' });

      Single_core_main_handle_Set = true;

   }

   var NumberOfRows = phases_single_core;
   
   CurrentNumberOfRows = Single_core_main_handle.getRowsNum()
       
   if( CurrentNumberOfRows != NumberOfRows ) {
      //Remove rows or add rows, puts 0 in added data
      resetDataGridRows(Single_core_main_handle,NumberOfRows,'')
   }   

}


// DESCRIPTION: This function sets the number of rows/columns of the table containing the conductors/insulators data of 
//              single-core cables when the user press enter in the second column (Number of conductors entry) of the main table
function Set_Single_core_DATA_GridFromEnter(event){
   //if the user pressed enter
   if (event.keyCode === 13) {
      Set_Single_core_DATA_Grid()
  }
}


// DESCRIPTION: the grid of single-core insulators/conductors. This function is called when the main table is changed. 
//              The number of rows in Single_core_data_handle depends on the number of conductors entered in Single_core_main_handle
function Set_Single_core_DATA_Grid() {
   
   var Ncolumns_singlecore;
   var GridCol_Width;
   var Gridprops;
   var master_string = '';
   var sum_master_string_elements = 0;
   var Nconductors_cable_n

   DataValues = Single_core_main_handle.serializeToCSV(); 
   DataValues = DataValues.split('\n'); 

   // Read the number of conductors for each single-core cable in MAIN grid
   for (cable_n=0; cable_n<phases_single_core; cable_n++) {
      idrow = Single_core_main_handle.getRowId(cable_n)  //currently tested row
      Nconductors_cable_n = Single_core_main_handle.cells(idrow,1).getValue();

      if(Nconductors_cable_n>3){ // Catch case of entering number of conductors > 3 (not allowed)
         alert('The maximum "Number of conductors" for single-core cables is 3.')
         Nconductors_cable_n = 3;
         Single_core_main_handle.cells(idrow,1).setValue('3')
      }

      if(cable_n == 0){     master_string = Nconductors_cable_n;
      }else{                master_string = master_string + ',' + Nconductors_cable_n;     } 

      if( !isNaN( parseInt(Nconductors_cable_n) ) ){
         sum_master_string_elements = sum_master_string_elements + parseInt(Nconductors_cable_n);
      }

   } 

   // Comments:
   // 'master_string' contains the number of conductors of each cable. This information is used to 
   // generate the Conductors/Insulators table to indicate the cable each conductor belong to.
   // 'sum_master_string_elements' contains the sum of conductors of all single-core cables.

   if (sum_master_string_elements > 0) { // if the column 'Number of conductors' in the main table is not empty, create the table

      var metric_units = 'Metric'.localeCompare(length_units.value);    
      
      // Define grid labels (depending on units)
      if(!metric_units){ 
         GridLabels = "Cable,Conductor,Phase,Inner radius (cm),Outer radius (cm),Conductor resistivity (Ohm m),Conductor relative permeability,Conductor relative permettivity,Insulator relative permittivity,Insulator loss factor"
      }else{
         GridLabels = "Cable,Conductor,Phase,Inner radius (in),Outer radius (in),Conductor resistivity (Ohm mile),Conductor relative permeability,Conductor relative permettivity,Insulator relative permittivity,Insulator loss factor"
      }
      Ncolumns_singlecore = 10
      GridCol_Width = '100,120,100,130,130,130,130,130,130,130';
      Gridprops = 'None,RealOrIntNotNamedNotZero,IntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosZeroNotNamed';

      
      if(Single_core_data_handle_Set){ // Single_core_data_handle_Set is a boolean variable, true if the table exists
         //Remember the old configuration of data, and print on new size
         
         var old_table = Single_core_data_handle.serializeToCSV();
         var new_table;
         old_table = old_table.split('\n')

         for(var n=0; n<old_table.length; n++) {
            row_n = old_table[n];
            row_n = row_n.split('\t')
            if(n==0){
               new_table = row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8]+'\t'+ row_n[9] + '\n';
            }else{
               new_table = new_table + row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8]+'\t'+ row_n[9] + '\n';
            }
         }
         SingleCoreDataGrid.value = new_table;

      }

      Single_core_data_handle = CreateNColumnGrids_slave(Ncolumns_singlecore,0,
                                 Single_core_dataGridPlace,
                                 Single_core_data_column,      //where we also put title
                                 SingleCoreDataGrid,          //hidden data location
                                 GridLabels,
                                 GridCol_Width,
                                 Gridprops,
                                 master_string);

      Single_core_data_handle.setUserData('','BlockedKeys','Insert','CtrlDelete');
      Single_core_data_handle_Set = true;
   }

}


// +DESCRIPTION: Test the number of single-core cables entered
//              Called when there is a change by entering a number or by pushing the button +1 or -1
function testNphases_single_core(){

           //*no #..# usage also rejects ##                limit to 50 conductors
   if( !notspecialopt("Nphases_single_core") || parseInt(Nphases_single_core.value)>150 ){    
      Nphases_single_core.value = phases_single_core;
      return
   }

   var retError = testpos0intvalue(Nphases_single_core,phases_single_core+''); 
      
   if(!retError){

      // The entered value was accepted, reset options and tables for the new value
      phases_single_core = parseInt(Nphases_single_core.value);

      setOptions_MainTable();       // Set the options in the main table 
      set_ListOfTables();           // Set the list of tables containing data (list of checkboxes)   
      set_DataTables();             // Set the data tables
      ShowHide_dataTables();        // show/hide data tables according to the status on the list of tables (checked/unchecked)
      ShowHide_canvas();            // May hide canvas in case all phases are equal to 0

      // update the value of previous phases
      previousphases_single_core = phases_single_core;

   }

   return

}       
   


// DESCRIPTION: This function sets the rows/columns of the main data table for single-core cables.
function reset_single_core_main_grid(){

   if (phases_single_core > 0) {

      if (Single_core_main_handle_Set){
         Single_core_main_handle_Set = false;
         SingleCoreMainGrid.value = Single_core_main_handle.serializeToCSV();
         Set_Single_core_MAIN_Grid();
      }

      if (Single_core_data_handle_Set){
         DataValues = Single_core_data_handle.serializeToCSV();
         
         // remove the first element of each row (inserted automatically)
         if( DataValues.length > 0 ) { //only if data exists

            DataValues = DataValues.split('\n');
            var irow = 0;
            NewDataValues = '';

            for (irow = 0; irow < DataValues.length; irow++) {
               values_irow = DataValues[irow];
               values_irow = values_irow.split('\t')
               values_irow = values_irow[1] + '\t' + values_irow[2] + '\t' + values_irow[3] + '\t' + values_irow[4] + '\t' + values_irow[5] + '\t' + values_irow[6] + '\t' + values_irow[7];

               if (irow == 0) {
                  NewDataValues = values_irow + '\n'
               }else{
                  NewDataValues = NewDataValues + values_irow + '\n'
               }
            }
         }

         SingleCoreDataGrid.value = NewDataValues;
         Set_Single_core_DATA_Grid();
      }
   }
   
}


// +DESCRIPTION: Callback for the button 'Use Cable Database'
function loadFromCableDatabaseButtonClick(){

   // Identify the selected row in the table of single core cables
   var selected_Row_Id = Single_core_main_handle.getSelectedRowId();

   if(selected_Row_Id==null){
      alert("Select a row in the 'Single-core cable main data' table")
      return
   }

   // get the index of the selected row
   var row = Single_core_main_handle.getRowIndex(selected_Row_Id);

   // read field units from LineCable_Data device, used to set the units in the database GUI
   var units = document.getElementById("UnitChoice").value;

   // Create the database object
   var emtpdb= new emtpCableDatabase();
   emtpdb.ACTIVE_UNIT_SYSTEM = (units=="English") ? "English" : "Metric";

   // put the database object in global values (myglobal)
   myglobal.emtpdb = emtpdb;
   myglobal.EMTPDB_FILELIST = emtpdb.listOfDatabaseXMLFile(myglobal, "Cable");

   // =================== Open database GUI =============================================================

   var CableSelection = emtpdb.load_SingleCore_Cable_From_Database(myglobal); 

   if(CableSelection == null) return ; // no cable selected

   var database_filepath = CableSelection[0];  // database file (could be the default or an user-created database)
   var           cableid = CableSelection[1];  // id (name) of the selected cable

   // ====================================================================================================

   // Define required units
   if(units=='Metric'){   RequiredAreaUnits = 'cm2';   RequiredDistanceUnits = 'cm';  }
   else{                  RequiredAreaUnits = 'in2';   RequiredDistanceUnits = 'in';  }

   // Define the number of conductors for the selected single-core cable from the database
   // There can be at most 4 conductors: core (always required), sheath, neutral and armor 
   NumConductors = 1;  // Start at 1 accounting for the core conductor
   
   // Define the Core outer radius. From the Cabledatabase, the core cannot be a hollow conductor, then its inner radius is always zero.
   // From the database, the size of the core is given , then it is read and used to get its radius.
   var core_area = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "CoreSize")
   var core_area_units = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "CoreSizeUnit")

   if(core_area_units != RequiredAreaUnits){ // Check if units from database coincide with the required ones, if not, do conversion
      core_area = emtpdb.getAreaValue(core_area, core_area_units, RequiredAreaUnits);
   }
   var core_radius = Math.sqrt(core_area/3.141592653589793);
   var OuterRadius = parseFloat(core_radius);  // partial outer radius (we add more depending on the existing layers)

   // Find the thickness of the insulator covering the core
   var InsulationThickness = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationThickness")
   var InsulationThicknessUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationThicknessUnit")
   if(InsulationThicknessUnits != RequiredDistanceUnits){
      InsulationThickness = emtpdb.getDistanceValue(InsulationThickness, InsulationThicknessUnits, RequiredDistanceUnits)
   }
   OuterRadius = OuterRadius + parseFloat(InsulationThickness);  // partial outer radius
   var Sheath_inner_radius = OuterRadius;  // The outer radius of the insulator equals the inner radius of the sheath (if considered)
   // If sheath is not consider, 'Sheath_inner_radius' won't be used.
   // This sheath inner radius will be overwritten in case semiconductor screens are involved, considered below.

   // Check if there is inner and/or outer semiconductor layers under/above core insulation
   var SemiconductorScreens = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationSemiconductingScreens")
   
   if(SemiconductorScreens=='Inner and Outer'){
      var InnerScreenThickness = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationInnerScreenThickness");
      var InnerScreenThicknessUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationInnerScreenThicknessUnit");
      if(InnerScreenThicknessUnits != RequiredDistanceUnits){
         InnerScreenThickness = emtpdb.getDistanceValue(InnerScreenThickness, InnerScreenThicknessUnits, RequiredDistanceUnits)
      }
      var OuterScreenThickness = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationOuterScreenThickness");
      var OuterScreenThicknessUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationOuterScreenThicknessUnit");
      if(OuterScreenThicknessUnits != RequiredDistanceUnits){
         OuterScreenThickness = emtpdb.getDistanceValue(OuterScreenThickness, OuterScreenThicknessUnits, RequiredDistanceUnits)
      }
      OuterRadius = OuterRadius + parseFloat(InnerScreenThickness) + parseFloat(OuterScreenThickness);

   }else if(SemiconductorScreens=='Inner'){
      var InnerScreenThickness = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationInnerScreenThickness");
      var InnerScreenThicknessUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationInnerScreenThicknessUnit");
      if(InnerScreenThicknessUnits != RequiredDistanceUnits){
         InnerScreenThickness = emtpdb.getDistanceValue(InnerScreenThickness, InnerScreenThicknessUnits, RequiredDistanceUnits);
      }
      OuterRadius = OuterRadius + parseFloat(InnerScreenThickness);

   }else if(SemiconductorScreens=='Outer'){
      var OuterScreenThickness = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationOuterScreenThickness");
      var OuterScreenThicknessUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationOuterScreenThicknessUnit");
      if(OuterScreenThicknessUnits != RequiredDistanceUnits){
         OuterScreenThickness = emtpdb.getDistanceValue(OuterScreenThickness, OuterScreenThicknessUnits, RequiredDistanceUnits)
      }
      OuterRadius = OuterRadius + parseInt(OuterScreenThickness);
   }

   // Find if Sheath conductor exists. If yes, get parameters
   var SheathType = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "SheathType")
   if( SheathType != "None" ){
      NumConductors++; // Add one to the counter for total number of conductors
      var SheathThickness = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "SheathThickness");
      var SheathThicknessUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "SheathThicknessUnit");
      if(SheathThicknessUnits != RequiredDistanceUnits){
         SheathThickness = emtpdb.getDistanceValue(SheathThickness, SheathThicknessUnits, RequiredDistanceUnits);
      }
      var Sheath_inner_radius = OuterRadius;
      OuterRadius = OuterRadius + parseFloat(SheathThickness);
      var Sheath_outer_radius = OuterRadius;
   }

   // Find if Neutral conductor exists. If yes, get parameters
   var NeutralType = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ConcentricNeutralType")
   if( NeutralType != "None" ){
      NumConductors++; // Add one to the counter for total number of conductors
      var NeutralWireSize = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ConcentricNeutralWireSize");
      var NeutralWireSizeUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ConcentricNeutralWireSizeUnit");
   
      if(NeutralWireSizeUnits != RequiredAreaUnits){
         NeutralWireSize = emtpdb.getAreaValue(NeutralWireSize, NeutralWireSizeUnits, RequiredAreaUnits);
      }

      var NeutralWireRadius = Math.sqrt(NeutralWireSize/3.141592653589793);

      var Neutral_inner_radius = OuterRadius;
      OuterRadius = OuterRadius + 2*parseFloat(NeutralWireRadius);
      var Neutral_outer_radius = OuterRadius;
   }

   // Find if Armor conductor exists. If yes, get parameters
   var ArmorType = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ArmorType")
   if( ArmorType != "None" ){
      NumConductors++; // Add one to the counter for total number of conductors
      var ArmorWireSize = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ArmorWireSize");
      var ArmorWireSizeUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ArmorWireSizeUnit");

      if( ArmorWireSizeUnits != RequiredAreaUnits ){
         ArmorWireSize = emtpdb.getAreaValue(ArmorWireSize, ArmorWireSizeUnits, RequiredAreaUnits);
      }

      var ArmorRadius = Math.sqrt(ArmorWireSize/3.141592653589793);      

      var Armor_inner_radius = OuterRadius;
      OuterRadius = OuterRadius + 2*ArmorRadius;
      var Armor_outer_radius = OuterRadius;
   }

   // Find Jacket outer radius, this corresponds to the total outer cable radius
   JacketMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "JacketMaterialID");
   if( JacketMaterialID != "None" ){
      var JacketThickness = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "JacketThickness");
      var JacketThicknessUnits = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "JacketThicknessUnit");
      if(JacketThicknessUnits != RequiredDistanceUnits){
         JacketThickness = emtpdb.getDistanceValue(JacketThickness, JacketThicknessUnits, RequiredDistanceUnits);
      }
      var OuterRadius = OuterRadius + parseFloat(JacketThickness);
   }

   // Populate cells in Single-core Main data grid (first table)
   Single_core_main_handle.cellByIndex(row,1).setValue(NumConductors);
   Single_core_main_handle.cellByIndex(row,4).setValue(OuterRadius.toFixed(8));

   // Define the Single-core Data grid, the number of rows depend on "NumConductors" in Main data grid 
   Set_Single_core_DATA_Grid()

   // Find the total number of single-core cables
   var Nsingle_core_cables = document.getElementById("NbSinglecore").value; 

   // Define the number of rows in the Single-core conductors/insulators data grid (second table)
   var nrows_data_grid = 0;
   for(ii=0; ii<Nsingle_core_cables; ii++){
      if( ii == row ){
         var core_row = nrows_data_grid;
      }
      var Ncond_ii_cable = Single_core_main_handle.cellByIndex(ii,1).getValue();
      if( Ncond_ii_cable != ''){
         nrows_data_grid = nrows_data_grid + parseInt(Ncond_ii_cable);
      }
   }

   // Core conductor material parameters
   var CoreMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "CoreMaterialID");
   var CoreRho = emtpdb.getConductorMaterialParameter(database_filepath, CoreMaterialID, "ElectricalResistivity")
   var CoreConductorRelativePermeability = emtpdb.getConductorMaterialParameter(database_filepath, CoreMaterialID, "RelativePermeability")
   // Insulator material covering the core conductor
   var InsulationMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "InsulationMaterialID")
   var CoreInsulatorRelativePermittivity = emtpdb.getInsulatorMaterialParameter(database_filepath, InsulationMaterialID, "RelativePermittivity")
   var CoreInsulatorLossFactor = emtpdb.getInsulatorMaterialParameter(database_filepath, InsulationMaterialID, "LossFactor")

   // Populate the conductors/insulators materials for the core conductor
   Single_core_data_handle.cellByIndex(core_row,3).setValue(0);
   Single_core_data_handle.cellByIndex(core_row,4).setValue(core_radius.toFixed(8));
   Single_core_data_handle.cellByIndex(core_row,5).setValue(CoreRho);
   Single_core_data_handle.cellByIndex(core_row,6).setValue(1); // Conductor relative permittivity (not available in database)
   Single_core_data_handle.cellByIndex(core_row,7).setValue(CoreConductorRelativePermeability);
   Single_core_data_handle.cellByIndex(core_row,8).setValue(CoreInsulatorRelativePermittivity);
   Single_core_data_handle.cellByIndex(core_row,9).setValue(CoreInsulatorLossFactor);
   var last_row = core_row;

   // If sheath conductor exists, populate the corresponding table row
   if( SheathType != "None" ){  
      var sheath_row = last_row + 1;
      // Sheath conductor material parameters
      var SheathMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "SheathMaterialID");
      var SheathRho = emtpdb.getConductorMaterialParameter(database_filepath, SheathMaterialID, "ElectricalResistivity")
      var SheathConductorRelativePermeability = emtpdb.getConductorMaterialParameter(database_filepath, SheathMaterialID, "RelativePermeability")
      // Insulator material covering the sheath
      // Note by Jesus: I think that the material covering the sheath should not be necessarily the jacket, 
      // I believe there could be another insulator layer, but database by Manu was made like that. Consider modifying
      var JacketMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "JacketMaterialID")
      var SheathInsulatorRelativePermittivity = emtpdb.getInsulatorMaterialParameter(database_filepath, JacketMaterialID, "RelativePermittivity")
      var SheathInsulatorLossFactor = emtpdb.getInsulatorMaterialParameter(database_filepath, JacketMaterialID, "LossFactor")

      // Populate the conductors/insulators materials for the sheath conductor
      Single_core_data_handle.cellByIndex(sheath_row,3).setValue(Sheath_inner_radius.toFixed(8));
      Single_core_data_handle.cellByIndex(sheath_row,4).setValue(Sheath_outer_radius.toFixed(8));
      Single_core_data_handle.cellByIndex(sheath_row,5).setValue(SheathRho);
      Single_core_data_handle.cellByIndex(sheath_row,6).setValue(1);
      Single_core_data_handle.cellByIndex(sheath_row,7).setValue(SheathConductorRelativePermeability);
      Single_core_data_handle.cellByIndex(sheath_row,8).setValue(SheathInsulatorRelativePermittivity);
      Single_core_data_handle.cellByIndex(sheath_row,9).setValue(SheathInsulatorLossFactor);
      last_row++
   }

   // If Neutral conductor exists, populate the corresponding table row
   if(NeutralType != "None"){
      var neutral_row = last_row + 1;
      // Neutral conductor material parameters
      var NeutralMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ConcentricNeutralMaterialID");
      var NeutralRho = emtpdb.getConductorMaterialParameter(database_filepath, NeutralMaterialID, "ElectricalResistivity")
      var NeutralConductorRelativePermeability = emtpdb.getConductorMaterialParameter(database_filepath, NeutralMaterialID, "RelativePermeability")
      // Insulator material covering the neutral conductor
      var JacketMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "JacketMaterialID")
      var NeutralInsulatorRelativePermittivity = emtpdb.getInsulatorMaterialParameter(database_filepath, JacketMaterialID, "RelativePermittivity")
      var NeutralInsulatorLossFactor = emtpdb.getInsulatorMaterialParameter(database_filepath, JacketMaterialID, "LossFactor")

      // Populate the conductors/insulators materials for the neutral conductor/insulator
      Single_core_data_handle.cellByIndex(neutral_row,3).setValue(Neutral_inner_radius.toFixed(8));
      Single_core_data_handle.cellByIndex(neutral_row,4).setValue(Neutral_outer_radius.toFixed(8));
      Single_core_data_handle.cellByIndex(neutral_row,5).setValue(NeutralRho);
      Single_core_data_handle.cellByIndex(neutral_row,6).setValue(1);
      Single_core_data_handle.cellByIndex(neutral_row,7).setValue(NeutralConductorRelativePermeability);
      Single_core_data_handle.cellByIndex(neutral_row,8).setValue(NeutralInsulatorRelativePermittivity);
      Single_core_data_handle.cellByIndex(neutral_row,9).setValue(NeutralInsulatorLossFactor);
      last_row++
   }

   // If Armor conductor exists, populate the corresponding table row
   if( ArmorType != "None" ){ 
      var armor_row = last_row + 1;
      // Armor conductor material parameters
      var ArmorMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "ArmorMaterialID");
      var ArmorRho = emtpdb.getConductorMaterialParameter(database_filepath, ArmorMaterialID, "ElectricalResistivity")
      var ArmorConductorRelativePermeability = emtpdb.getConductorMaterialParameter(database_filepath, ArmorMaterialID, "RelativePermeability")
      // Insulator material covering the armor. This must be the outer jacket.
      var JacketMaterialID = emtpdb.getSingleCoreCableParameter(database_filepath, cableid, "JacketMaterialID")
      var ArmorInsulatorRelativePermittivity = emtpdb.getInsulatorMaterialParameter(database_filepath, JacketMaterialID, "RelativePermittivity")
      var ArmorInsulatorLossFactor = emtpdb.getInsulatorMaterialParameter(database_filepath, JacketMaterialID, "LossFactor")

      // Populate the conductors/insulators materials for the armor conductor/insulator
      Single_core_data_handle.cellByIndex(armor_row,3).setValue(Armor_inner_radius.toFixed(8));
      Single_core_data_handle.cellByIndex(armor_row,4).setValue(Armor_outer_radius.toFixed(8));
      Single_core_data_handle.cellByIndex(armor_row,5).setValue(ArmorRho);
      Single_core_data_handle.cellByIndex(armor_row,6).setValue(1);
      Single_core_data_handle.cellByIndex(armor_row,7).setValue(ArmorConductorRelativePermeability);
      Single_core_data_handle.cellByIndex(armor_row,8).setValue(ArmorInsulatorRelativePermittivity);
      Single_core_data_handle.cellByIndex(armor_row,9).setValue(ArmorInsulatorLossFactor);
   }
   
}








// **************************************************************************************************
//  Section 8. Functions related to pipe-type cables
// ************************************************************************************************ */

// DESCRIPTION: This function sets the rows/columns of the main data table for pipe-type cables. 
function Set_Pipe_cable_MAIN_Grid() {  

   if(!Pipe_cable_main_handle_Set){

      var metric_units = 'Metric'.localeCompare(length_units.value);

      if(!metric_units){ 
         GridLabels = "Cable,Number of conductors,Horizontal position (m),Ground depth (m),Radius (cm),Pipe inner insulator relative permittivity,Pipe inner insulator loss factor"
      }else{
         GridLabels = "Cable,Number of conductors,Horizontal position (ft),Ground depth (ft),Radius (in),Pipe inner insulator relative permittivity,Pipe inner insulator loss factor"
      }

      num_char = 'None,IntPosNotNamed,RealOrIntNotNamed,RealOrIntNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed';

      Pipe_cable_main_handle = CreateNColumnGrids(8,0,
                                 Pipe_cable_mainGridPlace,
                                 Pipe_cable_main_grid,      //where we also put title
                                 PipeCableMainGrid,         // grid content
                                 GridLabels,
                                 '110,130,130,130,130,130,130',
                                 num_char);        

      Pipe_cable_main_handle.setUserData('','BlockedKeys','Insert','CtrlDelete');
 
      Pipe_cable_main_handle_Set = true;

   }

   var NumberOfRows = phases_pipe_cable;
   CurrentNumberOfRows = Pipe_cable_main_handle.getRowsNum()
       
   if( CurrentNumberOfRows != NumberOfRows ) {
      //Remove rows or add rows, puts 0 in added data
      resetDataGridRows(Pipe_cable_main_handle,NumberOfRows,'')
   }   

}


// DESCRIPTION: This function sets the number of rows/columns of the table containing the conductors/insulators data of pipe-type cables 
//              when the user press enter in the second column (Number of conductors) of the main table
function Set_Pipe_cable_DATA_GridFromEnter(event){
   //if the user pressed enter
   if (event.keyCode === 13) {
      Set_Pipe_cable_DATA_Grid()
  }
}


//* DESCRIPTION: This function sets the number of rows/columns of the table containing the conductors/insulators data of pipe-type cables
function Set_Pipe_cable_DATA_Grid() { 

   var stranded = document.getElementById("strandedbox").checked
   // Define labels 
   var master_string = '';
   var sum_master_string_elements = 0;
   var width_columns;
   var Num_columns;
   var Columns_prop;
   
   DataValues = Pipe_cable_main_handle.serializeToCSV();
   DataValues = DataValues.split('\n'); 

   //* Read the number of conductors for each single-core cable in MAIN grid
   for (cable_n=0; cable_n<phases_pipe_cable; cable_n++) {
      row_n = DataValues[cable_n];
      row_n = row_n.split('\t');
      cablenn = cable_n + 1;
      if(cable_n == 0){
         master_string = row_n[0];
         sum_master_string_elements = row_n[0];
      }else{
         master_string = master_string + ',' + row_n[0]; 
         sum_master_string_elements = sum_master_string_elements + row_n[0];
      } 
   }
 
   if (sum_master_string_elements > 0) { // if the column 'Number of conductors' in the main table is not empty, do the table

      var metric_units = 'Metric'.localeCompare(length_units.value);    

      if (stranded){

         if(!metric_units){ 
            GridLabels = "Pipe cable,Conductor,Phase,Distance from pipe center (cm),Angle (deg),Inner radius (cm),Outer radius (cm),Number of stranded conductors,Conductor resistivity (Ohm m),Conductor relative permeability,Conductor relative permettivity,Insulator relative permittivity,Insulator outer radius (cm),Insulator loss factor"
         }else{
            GridLabels = "Pipe cable,Conductor,Phase,Distance from pipe center (in),Angle (deg),Inner radius (in),Outer radius (in),Number of stranded conductors,Conductor resistivity (Ohm m),Conductor relative permeability,Conductor relative permettivity,Insulator relative permittivity,Insulator outer radius (in),Insulator loss factor"
         }
         width_columns = '100,100,100,130,130,130,130,130,130,130,130,130,130,130';
         Num_columns=14;
         Columns_prop = 'None,RealOrIntNotNamedNotZero,IntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosZeroNotNamed';
         
         if(Pipe_cable_data_handle_Set){ // Remove the first colum, the slave grid functions adds it automatically
           
            var old_table = Pipe_cable_data_handle.serializeToCSV();
            var new_table;
            old_table = old_table.split('\n');
      
            for(var n=0; n<old_table.length; n++) {
               row_n = old_table[n];
               row_n = row_n.split('\t');
               if(n==0){
                  new_table = row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] + '\t' + row_n[11] + '\t' + row_n[12] + '\n';
               }else{
                  new_table = new_table + row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] + '\t' + row_n[11] + '\t' + row_n[12]  +'\n';
               }
            }
            PipeCableDataGrid.value = new_table;
            
         }

      }else{

         if(!metric_units){ 
            GridLabels = "Pipe cable,Conductor,Phase,Distance from pipe center (cm),Angle (deg),Inner radius (cm),Outer radius (cm),Conductor resistivity (Ohm m),Conductor relative permeability,Conductor relative permettivity,Insulator relative permittivity,Insulator outer radius (cm),Insulator loss factor"
         }else{
            GridLabels = "Pipe cable,Conductor,Phase,Distance from pipe center (in),Angle (deg),Inner radius (in),Outer radius (in),Conductor resistivity (Ohm m),Conductor relative permeability,Conductor relative permettivity,Insulator relative permittivity,Insulator outer radius (in),Insulator loss factor"
         }
         width_columns = '100,115,100,130,130,130,130,130,130,130,130,130,130';
         Num_columns=13;
         Columns_prop = 'None,RealOrIntNotNamedNotZero,IntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosZeroNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosNotNamed,RealOrIntPosZeroNotNamed';

         if(Pipe_cable_data_handle_Set){ // Remove the first colum, the slave grid functions adds it automatically
            
            var old_table = Pipe_cable_data_handle.serializeToCSV();
            var new_table;
            old_table = old_table.split('\n');
      
            for(var n=0; n<old_table.length; n++) {
               row_n = old_table[n];
               row_n = row_n.split('\t');
               if(n==0){
                  new_table = row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10] + '\t' + row_n[11] +'\n';
               }else{
                  new_table = new_table + row_n[1] +'\t'+ row_n[2] +'\t'+ row_n[3] +'\t'+ row_n[4] +'\t'+ row_n[5] +'\t'+ row_n[6] +'\t'+ row_n[7] +'\t'+ row_n[8] +'\t'+ row_n[9] +'\t'+ row_n[10]  + '\t' + row_n[11] + '\n';
               }
            }
            PipeCableDataGrid.value = new_table;

         }
      }

      Pipe_cable_data_handle = CreateNColumnGrids_slave(Num_columns,0,
                                 Pipe_cable_dataGridPlace,
                                 Pipe_cable_data_column,      //where we also put title
                                 PipeCableDataGrid,          //hidden data location
                                 GridLabels,
                                 width_columns,
                                 Columns_prop,
                                 master_string);

      Pipe_cable_data_handle.setUserData('','BlockedKeys','Insert','CtrlDelete');
      Pipe_cable_data_handle_Set = true;

      //* Indicate the pipe conductor in the table: this is removed since causes error when testing table
      /*k = master_string.split(',');
      row = -1;  // account for grid header (1 row by default)
      for (cable_n=0; cable_n<phases_pipe_cable; cable_n++) {
         row = row + parseInt(k[cable_n]);
         if(!isNaN(row)){    // it can be NaN when adding more cables, then the cell with number of conductors is empty
            id_row_of_pipe_conductor = Pipe_cable_data_handle.getRowId(row)
            conductor = Pipe_cable_data_handle.cells(id_row_of_pipe_conductor,1).getValue()
            var new_value = conductor + ' (pipe)';
            Pipe_cable_data_handle.cells(id_row_of_pipe_conductor,1).setValue(new_value)
         }
      }
      */


   }

}


// DESCRIPTION: Test the number of pipe-type cables entered
function testNphases_pipe_cable() {

          //*no #..# usage also rejects ##                limit to 50 conductors
   if( !notspecialopt("Nphases_pipe_cable") || parseInt(Nphases_pipe_cable.value)>50 ){    
      Nphases_pipe_cable.value = phases_pipe_cable;
      return
   }

   var retError = testpos0intvalue(Nphases_pipe_cable,phases_pipe_cable+''); 

   if(!retError){

      // The entered value was accepted, reset options and tables for the new value
      phases_pipe_cable = parseInt(Nphases_pipe_cable.value);

      setOptions_MainTable();       // Set the options in the main table 
      set_ListOfTables();           // Set the list of tables containing data (list of checkboxes)             
      set_DataTables();             // Set the data tables
      ShowHide_dataTables();        // show/hide data tables according to the status on the list of tables (checked/unchecked)
      ShowHide_canvas();            // May hide canvas in case all phases are equal to 0

      // update the value of previous phases
      previousphases_pipe_cable = phases_pipe_cable;

   }

   return
   
}       




// DESCRIPTION: This function changes the number of rows/columns of the table containing the conductors/insulators data of pipe-type cables
function reset_pipe_cable_main_grid(){

   if (phases_pipe_cable > 0) {
      if (Pipe_cable_main_handle_Set){
         Pipe_cable_main_handle_Set = false;
         PipeCableMainGrid.value = Pipe_cable_main_handle.serializeToCSV();
         Set_Pipe_cable_MAIN_Grid();
      }

      if (Pipe_cable_data_handle_Set){
         DataValues = Pipe_cable_data_handle.serializeToCSV();
         
         // remove the first element of each row (inserted automatically)
         if( DataValues.length > 0 ) { //only if data exists

            DataValues = DataValues.split('\n');
            var irow = 0;
            NewDataValues = '';

            for (irow = 0; irow < DataValues.length; irow++) {
               values_irow = DataValues[irow];
               values_irow = values_irow.split('\t')
               values_irow = values_irow[1] + '\t' + values_irow[2] + '\t' + values_irow[3] + '\t' + values_irow[4] + '\t' + values_irow[5] + '\t' + values_irow[6] + '\t' + values_irow[7] + '\t' + values_irow[8] + '\t' + values_irow[9];

               if (irow == 0) {
                  NewDataValues = values_irow + '\n'
               }else{
                  NewDataValues = NewDataValues + values_irow + '\n'
               }
            }
      }

         PipeCableDataGrid.value = NewDataValues;
         Set_Pipe_cable_DATA_Grid();
      }
   }
   
}


// DESCRIPTION: This function modifies (adds/removes columns) the pipe-type cable conductors tables when stranded checkbox is checked/unchecked
function modify_stranded(){
   var new_table_pipetype_data = ''
   var PipetypeBOOL;
   
   var stranded = document.getElementById("strandedbox").checked
   
   try{
      var PipetypeConducTable = Pipe_cable_data_handle.serializeToCSV();
      PipetypeBOOL = true;
   }catch(err){
      PipetypeBOOL = false;
   }

   if (PipetypeBOOL){
      
      PipetypeConducTable = PipetypeConducTable.split('\n');
      NumConductors = PipetypeConducTable.length;
      if (NumConductors>0){
         for (var k=0; k<NumConductors;k++){
            var PipetypeRow = PipetypeConducTable[k].split('\t');
            if (stranded){
               new_table_pipetype_data += PipetypeRow[1] + '\t' + PipetypeRow[2] + '\t' + PipetypeRow[3] + '\t' + PipetypeRow[4] + '\t' + PipetypeRow[5] + '\t' + "1" + '\t' + PipetypeRow[6] + '\t' + PipetypeRow[7] + '\t' + PipetypeRow[8] + '\t' + PipetypeRow[9] + '\t' + PipetypeRow[10] + '\t' + PipetypeRow[11] +  '\n'; 
            }else{
               new_table_pipetype_data += PipetypeRow[1] + '\t' + PipetypeRow[2] + '\t' + PipetypeRow[3] + '\t' + PipetypeRow[4] + '\t' + PipetypeRow[5] + '\t'  + PipetypeRow[7] + '\t' + PipetypeRow[8] + '\t' + PipetypeRow[9] + '\t' + PipetypeRow[10] + '\t' + PipetypeRow[11] + '\t' + PipetypeRow[12] + '\n'; 
            }
         }
         Pipe_cable_data_handle_Set = false;
         PipeCableDataGrid.value = new_table_pipetype_data;
         Set_Pipe_cable_DATA_Grid();
      }  
   }   
}






//**********************************************************************************
// Section 9. Other functions (line length and additional table definition)
//******************************************************************************** */


// DESCRIPTION: Test line/cable length value entered
function test_length() {
   var retError = testposrealvalue(Length,old_linecable_length);
   if(!retError){ old_linecable_length = Length.value; }
   return retError
}


//+ DESCRIPTION: This function is a specific table definition for the conductors/insulators data of 
//              single-core and pipe-type cable cases. For these cases, the first column of the table
//              corresponds to the number of cable that a conductor belongs to.
//              The number of rows in this table depends on the number of conductors given to each cable, 
//              this information is passed with the master_string variable
// 
function CreateNColumnGrids_slave(Ncolumns,NextraLines,GridPlace,GridPlaceTitle,GridDataHidden,
       GridHeaderString, GridColumnSizes, GridColumnValidators,master_string) {

   master_string = master_string.split(',')
   
   //*Setup the grid
   GridSetupData=GridSetup(GridPlace.id,GridHeaderString,GridColumnSizes);
   
   GridHandle=GridSetupData.GridHandle;  //this is the handle of the grid
   
   //*Resize
   var jj;
   var ResizeString='false'; var AlignString='left'; var EditString='rowcount'
   for(jj=1;jj<=Ncolumns;jj++){
      ResizeString+=',true'
      AlignString+= ',left'
      EditString += ',ed'
   }
   GridHandle.enableResizing(ResizeString);
   
   GridHandle.setColAlign(AlignString);  //first column is for numbering
   GridHandle.setColTypes(EditString);
   //rowcount is defined in GridSetup
   
   GridHandle.setSerializationLevel(false,false,false,false,false,false);
   
   GridHandle.enableAutoHeight(true,650)
   GridHandle.enableMultiline(false);
   //+Initialize
   GridHandle.init();
   
   
   //*Width and height of the complete grid, can be done only after the init call
   GridPlace.style.width =GridSetupData.total_width;
   GridPlace.style.height=350;
   
   
   //+Send to the grid, this must have \t <0-------
   DataValues = GridDataHidden.value; //has \t
   
   //*Add the first column to existing  columns-------------------------------------------------
   // \t is the separator
   // (DataValues.length) could be 0, after split, become 1
   DataValues = DataValues.split('\n');
   var irow = 0;
   var rid
   var kpos = 0;
   NumberedDataValues = '';
   
   for (irow = 0; irow < master_string.length; irow++) {
      rid = irow + 1;
      for (icond = 0; icond < master_string[irow]; icond++){
         conductor_number = icond + 1;
         if (DataValues[kpos] == null || DataValues[kpos]==undefined){
            NumberedDataValues = NumberedDataValues + rid + '\t' + conductor_number + '\n'
            kpos = kpos +1; 
         }else{
            NumberedDataValues = NumberedDataValues + rid + '\t' + conductor_number + '\t' + DataValues[kpos] + '\n'
            kpos = kpos +1;
         }
      }
   }
   
   GridHandle.parse(NumberedDataValues, "csv");
   
   //+Only data columns must be serialized
   GridHandle.setSerializableColumns(ResizeString);
   
   //+Add extra rows
   if(NextraLines > 0 ){
   RowsNum=(GridHandle.getRowsNum());
   //gives only the data rows, not the header
   if( (RowsNum+NextraLines) < 25 ){
       var irow
       for(irow=RowsNum+1;irow<RowsNum+NextraLines;irow++){
           GridHandle.addRow(GridHandle.getUID(),"" );  //added at the end irow+',,'
       }
       RenumberRows()
   }
   }
   
   
   //+Enable validation: true means validate
   GridHandle.enableValidation(ResizeString)
   GridHandle.setColValidators(GridColumnValidators);
   //The above functions are defined in the GridSetup function
   
   //*Return to reduced size height (it does not work otherwise)
   GridPlace.style.height=250
   
   //*Tooltip
   GridPlaceTitle.title=GridSetupData.tooltip
   //GridPlace.title     =GridSetupData.tooltip

   //*IsLine/CableData grid property, does not show up in for...in loops nor keys list
   // By Philippe Gauthier
   // This property is necessary for tables created with the function CreateNColumnGrids_slave.
   // These tables are used in LineCable_Data device for the data of conductors for single-core and
   // pipe-type cables
   Object.defineProperty(GridHandle, 'isDataGrid',{
      value: true,
      enumerable: false
   })
   // end of IsLine/CableData grid property
   
   return GridHandle
   
}

