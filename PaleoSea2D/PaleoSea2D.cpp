/* Infos on the plugin:
Author and copyright: Marc Laurencelle
Last updates: July 2018
*/

#include "stdifm.h"
#include "PaleoSea2D.h"
#include <random>
#include <stdio.h>
#include <time.h>
#include <signal.h>
#include <windows.h>
#include <iostream>

#include "MarcExportFFdataNCDF4.h" //special module to export results(t) to a netCDF binary file (also by Marc L.)

/* Future development ideas:
- ?TODO+(soon): Improve BETApermanentDirichletMassConcBCs_LateAtClayTop by adding conditions e.g. aver.Kclay(t)<Kthres before applying Diri BCs effectively: this will be important when I'll put a pseudo-clay thin&permeable layer right from the beginning of the simulation...
- Review the uses of 'doexportdata' in the code. Maybe make it more independent of globXi...?
- TODO+++(soon, ~easy): Enhance the nodal exports by preparing an extended array of node indeces to export, when altexportusernodes==true, so that PuserexportnodesA (or a new array for the merge) merges the currglobtop nodes & the user-specified steady nodes indeces. Ensure, also, that nodal-group-id values are assigned to the currglobtop nodes (using reserved grpIds < 0; let's start using -1 for the whole current top).
->>> Then, do a new better analysis of vertical hydraulic gradients in the clay aquitard, to see if my empirical gradh_fw prediction formulae work well. Check also if there's an easy way to infer h_sw from h_fw + other infos.
- Better distinguish nodal and elemental properties in the code (+links between nodes and elements, if still required)
- Better types/structures for arrays, more modern C++ based: read and try from http://stackoverflow.com/questions/16137953/is-there-a-function-to-copy-an-array-in-c-c
- Try to use OnTimeStepConstraint to make time reach exactly the clayyi(t) change-point, then to force the next dt to be very small. Not very important, but would be nice. Base this simple idea on Peter's advice at: http://forum.mikepoweredbydhi.com/index.php?topic=1312.msg3230#msg3230
- IfmEditProperties, or find another way to get the Cs (max. c) ref. value!
- Automatic&Wise change of density ratios to simplify computations in zones where concentrations become very low (-->0), if allowed by the solver(t)
- Review if the "if (BETAwithsourceterms && (currclayYi > 0))" is indeed appropriate or maybe the second part is not? (in cases where "source-termed" elements are located in the rock aquifer)

DONE:
x Change memory allocation of arrays for the C++ style, using new type [size] and delete[] pointer, instead of old C style with malloc and free.
x Storage of constant distributed params inside the plugin (e.g. nodal coordinates)
x Extract all Y coordinates @PreSimulation by array=for()...=IfmGetY, to avoid repeating this operation thousand times in loops and time steps...
x Create a procedure which displays the global progress of the simulation based on final-start time, and real-clock-time since start... Display at each minute for example...
x Solve/debug the memory issue bug which occurs when the plugin is used twice by running two times a simulation in the same FEFLOW GUI session.
*/

IfmModule g_pMod;  /* Global handle related to this plugin */

#define countof(x) (sizeof(x)/sizeof(*(x)))

#pragma region IFM_Definitions
/* --- IFMREG_BEGIN --- */
/*  -- Do not edit! --  */

static IfmResult OnInitModule (IfmModule);
static void OnExitModule (IfmModule);
static IfmResult OnBeginDocument (IfmDocument);
static void OnEndDocument (IfmDocument);
static void Serialize (IfmDocument, IfmArchive);
static void OnEditDocument (IfmDocument, Widget);
static void PreSimulation (IfmDocument);
static void PostSimulation (IfmDocument);
static void PreTimeStep (IfmDocument);
static void PostTimeStep (IfmDocument);
static void PreFlowSimulation (IfmDocument);
static void PostFlowSimulation (IfmDocument);
static void PreMassSimulation (IfmDocument, int);
static void PostMassSimulation (IfmDocument, int);
static IfmBool OnTimeStepConstraint (IfmDocument, double, double*);

/*
 * Enter a short description between the quotation marks in the following lines:
 */
static const char szDesc[] = 
  "This is Marc's elaborate plugin for running postglacial paleo-hydrogeological simulations (4rd gen. -07/2018)";

#ifdef __cplusplus
extern "C"
#endif /* __cplusplus */

IfmResult RegisterModule(IfmModule pMod)
{
  if (IfmGetFeflowVersion (pMod) < IFM_REQUIRED_VERSION)
    return False;
  g_pMod = pMod;
  IfmRegisterModule (pMod, "SIMULATION", "PALEOSEA2D", "PaleoSea2D", 0x1000);
  IfmSetDescriptionString (pMod, szDesc);
  IfmSetCopyrightPath (pMod, "PaleoSea2D.txt");
  IfmSetHtmlPage (pMod, "PaleoSea2D.htm");
  IfmSetPrimarySource (pMod, "PaleoSea2D.cpp");
  IfmRegisterProc (pMod, "OnInitModule", 1, (IfmProc)OnInitModule);
  IfmRegisterProc (pMod, "OnExitModule", 1, (IfmProc)OnExitModule);
  IfmRegisterProc (pMod, "OnBeginDocument", 1, (IfmProc)OnBeginDocument);
  IfmRegisterProc (pMod, "OnEndDocument", 1, (IfmProc)OnEndDocument);
  IfmRegisterProc (pMod, "Serialize", 1, (IfmProc)Serialize);
  IfmRegisterProc (pMod, "OnEditDocument", 1, (IfmProc)OnEditDocument);
  IfmRegisterProc (pMod, "PreSimulation", 1, (IfmProc)PreSimulation);
  IfmRegisterProc (pMod, "PostSimulation", 1, (IfmProc)PostSimulation);
  IfmRegisterProc (pMod, "PreTimeStep", 1, (IfmProc)PreTimeStep);
  IfmRegisterProc (pMod, "PostTimeStep", 1, (IfmProc)PostTimeStep);
  IfmRegisterProc (pMod, "OnTimeStepConstraint", 1, (IfmProc)OnTimeStepConstraint);
  IfmRegisterProc (pMod, "PreFlowSimulation", 1, (IfmProc)PreFlowSimulation);
  IfmRegisterProc (pMod, "PostFlowSimulation", 1, (IfmProc)PostFlowSimulation);
  IfmRegisterProc (pMod, "PreMassSimulation", 1, (IfmProc)PreMassSimulation);
  IfmRegisterProc (pMod, "PostMassSimulation", 1, (IfmProc)PostMassSimulation);
  return True;
}

static void Serialize (IfmDocument pDoc, IfmArchive pArc)
{
  Cpaleosea2d::FromHandle(pDoc)->Serialize (pDoc, pArc);
}
static void OnEditDocument (IfmDocument pDoc, Widget wParent)
{
  Cpaleosea2d::FromHandle(pDoc)->OnEditDocument (pDoc, wParent);
}
static void PreSimulation (IfmDocument pDoc)
{
  Cpaleosea2d::FromHandle(pDoc)->PreSimulation (pDoc);
}
static void PostSimulation (IfmDocument pDoc)
{
  Cpaleosea2d::FromHandle(pDoc)->PostSimulation (pDoc);
}
static void PreTimeStep (IfmDocument pDoc)
{
  Cpaleosea2d::FromHandle(pDoc)->PreTimeStep (pDoc);
}
static void PostTimeStep (IfmDocument pDoc)
{
  Cpaleosea2d::FromHandle(pDoc)->PostTimeStep (pDoc);
}
static IfmBool OnTimeStepConstraint (IfmDocument pDoc, double tNow, double* dtProposed)
{
  return Cpaleosea2d::FromHandle(pDoc)->OnTimeStepConstraint (pDoc, tNow, dtProposed);
}
static void PreFlowSimulation (IfmDocument pDoc)
{
  Cpaleosea2d::FromHandle(pDoc)->PreFlowSimulation (pDoc);
}
static void PostFlowSimulation (IfmDocument pDoc)
{
  Cpaleosea2d::FromHandle(pDoc)->PostFlowSimulation (pDoc);
}
static void PreMassSimulation (IfmDocument pDoc, int iSpecies)
{
  Cpaleosea2d::FromHandle(pDoc)->PreMassSimulation (pDoc, iSpecies);
}
static void PostMassSimulation (IfmDocument pDoc, int iSpecies)
{
  Cpaleosea2d::FromHandle(pDoc)->PostMassSimulation (pDoc, iSpecies);
}

/* --- IFMREG_END --- */
#pragma endregion

/* ================= BEGIN OF PLUGIN CUSTOM HEADER ================= */

/* !!! BETA MODES under test !!! */
bool BETAwithsourceterms; //NEW FEATURE 23nov; enabled if a betaclaysrcE elem. ref. distr. is present


const bool BETApermanentDirichletMassConcBCs_LateAtClayTop = true; //SPECIAL TEST 30oct:
	/* That's the best solution we've got to deal with the changing dominance of diffusion over advection,
	   once clays appear in the model. I've made a small R code to determine the approximate threshold Y velocities
	   at which this dominance changes, from dispersion or diffusion alone (vy < threshold) to advection (vy > threshold).
	   THE CODE IS: \R scripts\Competition advection vs diffusion front advances v0b.R */

/* Abandoned beta modes (but still present in the code) */
const bool BETApermanentDirichletMassConcBCs = false; //SPECIAL TEST 26oct: Makes the simulation of convection very hard and prevents sustained convergence of the model run. Abandoned.
const bool BETAautomUseMassBCconstraints = false; //SPECIAL TEST 26oct: Conclusion: BOF! Even with Divergence Form, it does not produce stable outputs, including a very oscillating Mass Imbalance (with the amplitude growing, very early after start!). Therefore, this option shall be abandoned, as well as the possibility of the FEFLOW constraint parameter, since it is not functioning properly.
const bool DisplayBudgetComput_infos = true;

const bool EXPORT_MASS_RBUDGETs = false; //NEW HARD OPTION 21july2017 to avoid strange bugs during FF internal budget calculation for this plugin

/* Arrays for ALL MODES - (entire model domain) */
bool *Ptopnodes; //1d array[nnodes] indicating nodes at the top of the model, where most transient updates(t) are taking place by this plugin; (from user nodal distrib. 'topnodes' in older modes; static in older modes, BUT dynamically redefined in the new withclaymesh mode)
bool *Ptopelems; //1d array[nelems] indicating elements which may need update of their TR parameter by this plugin (from user elem. distrib. 'topelemszmid'; static in older modes, but IGNORED & SKIPPED in the new withclaymesh mode)
bool *Pcurrfloodednodes; //1d array[nnodes] of currently flooded nodes (updated @PreTimeStep, after isflooded)
bool *Pinflownodes; //1d array[nnodes] of currently inflowing nodes (acc. to veloc. @PostFlowSim)
double *PallnodesY; //1d array[nnodes] of Y_local coordinate for all nodes (imported @PreSimulation); coordinates do not change during the simulation: X,Y are FIXED, STEADY.
int *PglobXiElemA; //NEW BETA 1d array[nelems] of horizontal 'X' integer indeces for all elements of the model

/* Arrays for optional Data Export (compatible with ALL MODES) */
int *PglobXiraw; //NEW 1d static array loaded @PreSimulation, for sorting top nodes according to their horizontal (X) order...
int *Puexportnodes; //NEW (TODO write the comment); ~raw rounded user export nodes values; full array
int *PuserexportnodesA; //NEW (TODO write comment); filtered (small) array for the user export nodes (altern. mode)
int *Puserexportnodes_uvaluesA; //NEW; complements PuserexportnodesA with the user data value assigned to these nodes
int *PeffexportnodesA; //BETA; effective small array of export nodes (altern. mode); include user-specified nodes, optionally preceded by curr.glob.top.nodes; BEWARE: it is simply a pointer to the PuserexportnodesA if nbglobXivals == 0 !!
int *Peffexportnodes_uvaluesA; //NEW; complements PeffexportnodesA with the user data value assigned to these nodes... or the automatic expNautotopId (-999) value for the curr.glob.top.nodes; BEWARE: it is simply a pointer to the Peffexportnodes_uvaluesA if nbglobXivals == 0 !!
/* Full arrays (although partially filled / updated) used as temporary buffer for exporting some data */
double *Pnodalmassrbudgets; //mass rate budget computed values; full-range but not necessarily updated for all nodes
/* Small arrays for the current global top nodes of the active model domain (sorted according to globXiraw);
   arrays range from 0 to maxglobXival (so length = max... +1) */
int *Pcurrglobtopnodes; //array of node indeces referring to the current glob.top.nodes (not necessarily well sorted; depends on the ordering of the nodes within the FEFLOW model)
bool *Pcurrglobtop_isinflow; //current inflow status for each glob.top.node; for display only for the moment...
bool *Pprevglobtop_isinflow; //PREVIOUS inflow status for each glob.top.node; for display only for the moment...
double *Pcurrglobtop_nodalmassrbudgets; //current mass rate budget for each glob.top.node; BETA, and ALMOST UNUSED for the moment...

/* Arrays for the NEW withclaymesh MODE ONLY - (entire model domain) */
bool *PtopclayBCnodes; //NEW BETA 1d bool array[nnodes]: current top of clay with BCs
bool *PtoprockBCnodes; //NEW BETA 1d bool array[nnodes]: stable top of rock with BCs
bool *PupdateNbelowclaytop; //NEW BETA 1d bool array[nnodes]: nodes below the current top of clay aquitard AND WHICH require initialization due to recent activation of the corresponding elements
bool *PupdateNoverclaytop; //NEW BETA 1d bool array[nnodes]: nodes over (above) the current top of clay aquitard AND WHICH should be neutralized due to recent deactivation of the corresponding elements
bool *PcleanNBCprevclaytop;  //NEW BETA 1d bool array[nnodes]...
int *PclayYimain; //NEW BETA 1d array[nnodes] of nodal clayYimain indeces >= 0 only for non-sub-rows of nodes in the clay aquitard (-1 elsewhere); static, from user data
int *PclayYisub; //NEW BETA 1d array[nnodes] of nodal clayYisub indeces >= 0 only for sub-rows of nodes in the clay aquitard (-1 elsewhere); static, from user data
int *PclayXiraw; //NEW BETA 1d array[nnodes] of nodal clayXiraw indeces >= 0 for all nodes within the clay aquitard (-1 elsewhere); static, from user data
int *Pclaylayeri; //NEW BETA 1d array[nelems] of elemental claylayeri indeces >= 1 for all elements within the clay aquitard (-1 elsewhere); static, from user data
double *PclaysrcBetaCoefE; //NEW BETA 1d array[nelems] of elemental 'betaclaysrcE' values = [REVIEWED Oct.2017] [width of the current column of elements] / [total volume of the elements receiving S&C-driven gw fluxes in this column (@ globXiE) of the model mesh]
bool *PelemcurractiveE; //NEW BETA 1d array[nelems] of current active state of the elements of the entire domain (updated along with the activ/deactivation process @UpdateClayParamsAccToYi) (BOOL instead of INT values)

/* Arrays for the OLDER MODES ~MOSTLY - (entire model domain) */
double *Ptopelemszmid; //1d array of zmid values for all elements (the useful values are for topelems) (from user elem. distrib. 'topelemszmid')
double *Pfloodedduring; //1d array of flood duration values at topnodes (updated @PostTimeStep)
double *Pclayacc_nodal; //1d array of clay accum. values at topnodes (updated @PostTimeStep)
double *Pclayaccfinal_nodal; //1d array of final clay accumulation targets at topnodes (constant!)
double *Pclayacc_elemtal; //1d array of clay accum. values at topelems (updated @PostTimeStep)
double *Pclaysedrate_nodal; //1d array of clay constant sedim. rates at topnodes nodes, in meters/day (defined @PreSimulation)
double *Pclaysedrate_elemtal; //1d array of clay constant sedim. rates at topelems, in meters/day (defined @PreSimulation)
long *PelemsAtnodesBR; //1d array of indeces of the element closest to the ith node without position gap (sufficient for quad. mesh)
long *PelemsAtnodesBL; //1d array of indeces of the element closest to the ith node WITH position gap (y-gap, x-gap) (required for triang. mesh)
long *PelemsAtnodesBM; //1d array of indeces of the element closest to the ith node WITH position gap (y-gap, x±0) (required for triang. mesh)

/* New withclaymesh mode */
int *Pcurrclaytopnodes; //NEW BETA 1d dynamic array identifying the nodes at the current top of the clay aquitard; (std 0..size-1 range, where size = nbclayXivals); dynamically updated
int *Pprevclaytopnodes; //NEW BETA 1d dynamic array ... = copy of Pcurrclaytopnodes made when clayyi changes, just before updating Pcurrclaytopnodes; (std 0..size-1 range, where size = nbclayXivals)
bool *Pcurrclaytopisinflow; //NEW BETA 1d dynamic array identifying the inflowing nodes at the current top of the clay aquitard; (std 0..size-1 range, where size = nbclayXivals); dynamically updated

/* Small arrays used by the new withclaymesh mode for initializing heads for nodes inside newly activated clay elements (when withhgrad==true) */
double *Pclaytopnodes_tmp_prevtophbc; //small array storing head values of the nodes of the previous clay-top; but filled and used only when hasclayYiincreased==true (i.e. clayYi changes by increasing: currclayYi > prevclayYi)
double *Pclaytopnodes_tmp_currtopnewhbc; //small array with the head values to be assigned to the nodes of the current clay-top
double *Pclaytopnodes_tmp_prevtopY; //small array storing Y values of the nodes of the previous clay-top; (see _prevtophbc for technical details)
double *Pclaytopnodes_tmp_currtopY; //small array storing Y values of the nodes of the current clay-top; (see _prevtophbc for technical details)
/* Small arrays used by the new withclaymesh mode for initializing mass conc. of newly activated nodes at and below the new clay top */
double *Pclaytopnodes_tmp_prevtopcbc; //small array storing mass conc. values of the nodes of the previous clay-top; see the similar ..._prevtophbc for additional details
double *Pclaytopnodes_tmp_currtopnewcbc; //small array with mass conc. values to be used for assigning the updated first-type BC's to the current clay top nodes (only for the relevant nodes as identified by '..._cbc_chgneeded'[] == true)
bool   *Pclaytopnodes_tmp_currtopcbc_chgneeded; //small array to use as a condition along with Pclaytopnodes_tmp_currtopnewcbc (which was filled with the cbcval_if_outflow values in the loop; see the code)

/* Full arrays used by the NEW withclaymesh MODE ONLY, for overwritting nodal values of the entire model
   at the end of the first=bad time step where clayYi has just changed (following advice from Carlos,
   email received on October 17th, 2016) */
double *Pallnodes_heads; //full array of nodal head values to overwrite model values @PostTimeStep right after clayYi has changed
double *Pallnodes_mconcs; //full array of nodal mass conc. values to overwrite model values @PostTimeStep right after clayYi has changed
bool *Pallcurractivenodes; //identifies all nodes which are currently active i.e. within the domain of activated elements
bool *Pallcurrclayactivenodes; //similar to Pallcurractivenodes but true only if &in the clay aquitard (base included)
bool *Psteadyallrocknodes; //identifies all nodes which are within the permanently-active subdomain of rock-facies elements
bool *Psteadytoprocknodes; //identifies all nodes which bound the rock aquifer at its top, confined or not by clays. (Currently for withclaymesh mode only)

/* GLOBAL PARAMETERS FOR ALL MODES */
int nnodes; //number of nodes (@PreSim); static
int nelems; //number of elements (@PreSim); static
int ntopnodes; //current number of topnodes (@PreSim); may change dynamically in the withclaymesh mode... (BETA)
int ncurractivenodes; //current number of active nodes in the domain; in the withclaymesh MODE: variable integer calculated when the bool vector Pallcurractivenodes is updated (during UpdateClayParamsAccToYi procedure); in the OLDER MODES: set equal to nnodes (constant).
int ninflownodes; //current number of inflowing topnodes; updated by calling UpdateInflowDataForTopNodes...
int noutflownodes; //current number of inflowing topnodes; updated by calling UpdateInflowDataForTopNodes...
int nbcctopnodes; //current number of topnodes having a BC constraint (Bcc); relevant only if special option BETAautomUseMassBCconstraints is enabled; count is updated in UpdateInflowDataForTopNodes and used only for display in the Log
int ninflownodes_gtop; //NEW (usable only if doexportdata==true, since it requires globXi...)
int noutflownodes_gtop; //NEW...
int nflowchgnodes_gtop; //NEW...
int ntopelems; //number of topelems (@PreSim); static; defined only for NON-withclaymesh modes
int nuserexportnodes; //NEW... it's a count of the non-null (>=0) values in Puexportnodes
int neffexportnodes; //NEWER, BETA...

bool domasstransport = true; //to skip mass-transport related features if the problem is flow-only
bool isglobXipresent = false; 
bool doexportdata = false; //NEW optional export to a fixed-path netCDF binary data file; assigned internally @PreSim (to true only if a globxiraw nodal user data is present)
bool altexportusernodes = false; //NEW alternative mode for data export, enabled (internally) if the user provides an 'exportnodes' nodal User Data
bool altexportautoinclgtop = false; //NEW BETA...
bool meshquadtype = true; //Mesh type (triangular vs. quadrilateral); changes the way transient TR are updated; QUAD.MESH is HIGHLY RECOMMENDED, HOWEVER!
bool divergformtransp = false; //NEW problem setting value (obtained from IfmIsDivergenceFormTransport @PreSimulation)
bool wiseTRupdate = false; //A way for updating transfer rates (TR) according to the topmidz of topelems; assigned internally @PreSim, but disabled (forced false) when the withclaymesh mode is active
bool withbotclayconc = false; //Improved way of assigning the concentration at top nodes, taking the retardation effect of the clay aquitard into account!; not applicable to the new withclaymesh mode

/* GENERAL internal booleans */
bool firstTSresumed; //true if simulation was resumed (i.e. started at a time > 0); returns to false RIGHT AFTER the first time step is done (or is always false if not resuming)
bool isFirstExportedTS; //NEW; true only for the first completed time step to be exported; false afterwards; CURRENTLY does NOT consider simulation Resume...
bool initialpreTSupdate = false; //New global param. for all modes, which has a true value only in the early first time step, to drive some initializations; in the new withclaymesh mode, forces a full update of the clay nodes & elements (through a UpdateClayParamsAccToYi call); then, whatever the mode, forces initialization of in/outflow array(s) and other variables (through a UpdateInflowDataForTopNodes call)
bool relevanttimestep; //boolean flag which indicates if the current time step is relevant (true) or should be ignored (false), e.g. when writing results(t) to disk
bool firstflowsol; //NEW
bool firstmasssol; //NEW

/* SPECIFIC to the new withclaymesh mode */
bool needupdatearoundclaytop; //NEW BETA internal boolean which tells if some nodes need a status update due to recent activation/deactivation of their corresponding elements...; only used by the new withclaymesh mode
bool hasclayYichanged; //indicates if the clayYi max-layer index has changed (so that the number of active clay layers will be updated in the running model); this logical is updated within UpdateClayParamsAccToYi calls...
bool hasclayYiincreased; //indicates if the clayYi max-layer index has increased; simply a more specific alternative to 'hasclayYichanged'; this logical also is updated within UpdateClayParamsAccToYi calls...

/* SPECIFIC to the botclayconc sub-mode */
double origclayL; //reference thickness for the clay aquitard, used to prepare the botclayconc TS (ADJUST and recompile if Lref is changed in the R code!)
double endofmaxsalinitytime; //time (in days) at which the salinity starts to drop down (BETA; TODO: should be defined by the user via a time series or smt similar); HARDWIRED for now, SO MUST BE ADAPTED TO NEW TS DATA!
double maxtransftimebclaycTS; //maximum transformed time (in days) beyond which the fwconc value may be used instead of querying extrapolation of the botclayconc time series; HARDWIRED for now, SO MUST BE ADAPTED TO NEW TS DATA!

/* TIME SERIES ids */
int rslPID = -1; //rsl TS identifier
int swcPID = -1; //swconc TS identifier
int bclaycPID = -1; //botclayconc TS identifier
int maxcPID = -1; //max. concentration constant TS identifier (BETA for extracting Cs reference value despite absence of appropriate IFM function)
int mincPID = -1; //min. concentration constant TS identifier (BETA for extracting C0 reference value despite absence of appropriate IFM function)
int fwconcPID = -1; //salinity of freshwater constant TS identifier
int sedratePID = -1; //sedrate TS identifier
int claytopPID = -1; //claytop TS identifier (constant 1-step TS)
int clayYiPID = -1; //clayyi TS identifier
int seddurPID = -1; //sedimduring TS identifier (constant 1-step TS)
int rndPID = -1; //random activation TS identifier
int clayeqkvPID = -1; //NEW clayeqkv TS identifier
int clayYiFixMBCPID = -1; //NEW clayYiFixMBC identifier (i.e. clayYi fix mass BC threshold value)
int clayqdbotPID = -1; //NEW clayqdbot TS identifier

/* MODE activity flags */
bool withclayTimeDepsedrate = false; //true if a sedrate time series is found
bool withclayConstsedrate = false; //true if a target elevation for the top of clay is specified instead of transient sedrate(t)
bool withclaymesh = false; //true if a clayyi time series is found (...and other conditions are met)
bool withtimedepKclays = false; //NEW true if a clayeqkv time series is found (...and in withclaymesh mode)
bool withSOMEclayaccumANYmode; //NEW summary boolean which indicates if the current mode(s) do accumulate clay(t), or not.
bool withtransientswconc = false; //true if a swconc time series is found, and applicable; ; if true, 'swconc' TS is used for the transient salinity(t) of the flooding waters; else the constant defswconc hard-wired value is used.
bool withrndnoiseaddtoconcBCs = false; //true if a rndconc time series is found, and applicable
bool isaglobxiEloaded = false; //NEW true if a globxiE elemental user data array could be successfully loaded

/* For allowing proper stopping of simulation with CTRL+C in Console mode */
bool globvar_stopsimul = false;

/* SOME CONSTANTS FOR THIS PLUGIN */
#define expNautotopId -999 //exportnode reserved Id for the auto.curr.glob.top.nodes joined to the exportnodes...
#define daysperyear 365.2425 //Factor for conversion of years into days
#define secsperday 86400.0 //Factor for conversion of days into seconds
#define K_TILLS 8.64E-2
#define K_CLAYS 8.64E-5
const double int_tol = 1e-6; //tolerance for rounding (or flooring) of real values when converting to integer values
const double len_tol = 1e-6; //tolerance for equality testing between length values (in meters)
const double rsl_tol = 1e-4; //tolerance for equality testing between relative sea level (RSL) values (in meters: 0.1 mm)
const double conc_tol = 1e-3; //tolerance for equality testing between concentration values (in mg/L)
const double sedrate_tol = 2.5e-8; //tolerance for equality testing between sedimentation rate values (in meters/day); tol set to ~0.01 mm/year
const double time_tol = 1.0 / secsperday; //tolerance for equality testing between time values (in days); tol set to 1 second
const double densr_tol = 1e-8; //tolerance for equality testing between density ratio values (dimensionless)
const double nodata_double = std::numeric_limits<double>::quiet_NaN();
const double defrandomduration = 10.0 * daysperyear; //Default duration of the addition of random noise onto fixed concentration nodal BCs (in days)

//...including a few mode-specific constants:
const int rockYi_cstval = -2; //NEW BETA clayyi value for top rock outside of the clay zone (where topBC will apply always); for withclaymesh mode

/* SOME DEFAULT VALUES... */
double swconc_atstart; //NEW for source term mainly
double defswconc; //Default concentration for seawater (in mg/L) - set early @PreSimulation (~Hard-Wired)
double deffwconc; //Default concentration for freshwater (in mg/L) - set early @PreSimulation (~Hard-Wired)
double fwconc = -9999.0; //Constant concentration of inflowing freshwater (in mg/L)
						 /* NOTE: This parameter may be used to set minconc = fwconc, if no 'minconc' constant TS is found. Otherwise, fwconc may be different, i.e. fwconc >= minconc. */
double maxconc = -9999.0; //Maximum concentration to come in the model (in mg/L), apart from numerical oscillation errors. This is the concentration related to the max. density implicitly defined through the domain-uniform density ratio.
						  /* DEVELOPER'S NOTE: maxconc value should be defined from extraction of Problem Settings'
						  Cs ref. value (but no Ifm function found yet). So, for now, the value is either assumed
						  default (from defswconc) or specified in a constant time series named 'maxconc'
						  (extracted at time=startingtime). */
double minconc = -9999.0; //Minimum concentration (in mg/L): C0 ref. value in the Problem Settings (see comment on maxconc for more details...), related to the min. density at which the density effect is neutral. Its exceptional default value is 0.0 mg/L, but the normal assignation method is by reading 'minconc' TS constant value, if present; otherwise, minconc = fwconc (if available; should always be).
double densr = -9999.0; //density ratio (from elem #1); set @PreSimulation; defaults to 0.0 if invalid value is read in the elem. data, or if Flow-Only simulation (!domasstransport).
double CauchyKinit = K_TILLS; //K value for tills = 1e-6 m/s, in m/d @ Fluid-Transfer BC topnodes (Cauchy for flow, 3rd kind) 
double Cauchybinit = 5.0; //Thickness of tills (already present at start) for the Fluid-Transfer BC topnodes
double CauchyKclays = K_CLAYS; //K value for clays = 1e-9 m/s, in m/d @ Fluid-Transfer BC topnodes (Cauchy for flow, 3rd kind) 

double currsedimrate = 0.0; //sedimentation rate, global transient, active only if !withclayConstsedrate; in meters per day (m/d); variable(t)
double claytoptargetYval = -1.0; //final elevation of the top of the forming clay aquitard; in meters following the Y coordinate system; static
double claysedimduring = -1.0; //duration (i.e. final absolute time) of the period of constant rate sedimentation (in days) (starting from t0 = 0.0 day, even when resuming simul.); static
double simulstarttime = -1.0; //IfmGetAbsoluteSimulationTime value at simulation start (in days); static
double simulendtime_apriori = -1.0; //IfmGetFinalSimulationTime value read at simulation start (in days); static, A PRIORI final time (as the user or some simulation error might stop the run before reaching this final time)
double currKofclays = nodata_double; //current GLOBAL equivalent(average) vertical K of clays (in m/d); optionally time-dependent (if a 'clayeqkv' TS is present), else it's the constant static hard-wired value K_CLAYS; BETA-TODO Someday make it more adaptative (esp. to the LOCAL final thickness of the aquitard!)
double currSrcGlobDarcyflux = nodata_double; //BETA! current GLOBAL reference fluid Darcy fluid-flow flux (in m/d) for the single-layer elem. source terms in clays; BETA-TODO Someday make it more adaptative (esp. to the LOCAL final thickness of the aquitard!)
bool areSrcTermsAlreadyZeroed = false; //BETA! initially we do not presume that the model's source terms are zeroed. (works along with isthecurrentGlobDfluxsignif; see in the code)

//Storage of previous-time-step value for a few time series
double lastrslval = -9999.0;
double lastswconc = -9999.0;
double lastsedimrate = -9999.0;

std::string problemfilename; //without folder path; general
std::string resultsfolderpath; //folder path only ("./results/")
std::string exportfilename; //without folder path; for netCDF-format data export
std::string exportfilepath;
double lastexporttime = -9999.0; //last sim. time at which data export was made (in days)
double exportdeltatime = 1.0 * daysperyear; //minimum time-duration between data export calls (internally in days); to export at all time steps, simply assign a very small value (or even <= 0.0) to the parameter.
const bool autodefine_exportdeltatime = true; //NEW OPTION

/* SPECIFIC to the withclaymesh mode */
int prevclayYi = -9999; //previous currclayYi value stored for detecting change...
int currclayYi = -9999; //id of the topmost clay layer at current time (must be an integer)
double currclayYi_rawreal; //raw real value as extracted from the time series (for display purposes only)
double clayYi_FixMassBC_thresv; //NEW threshold value that, when currclayYi_rawreal >= ..it.., Mass-Conc. BCs apply to every clay-top node
int maxclaylayeri = -9999; //maximum id value in the claylayeri elemental user data (integer); EQUIVALENT to the number (count) of different clay-layer indeces (because ranging from 1 to maxclaylayeri)
int minclaylayeri = +9999; //minimum id value in the claylayeri elemental user data (integer); MUST BE equal to 1!
int maxclayXival = -9999; //maximum id value in the clayxiraw nodal user data (integer)
int minclayXival = +9999; //minimum id value in the clayxiraw nodal user data (integer)
int nbclayXivals = 0; //number of different clayXi values, simply based on the max-min interval bounds (nb = max - min + 1)
int maxclayYival = -9999; //maximum id value in the clayyimain nodal user data (integer)
int minclayYival = +9999; //minimum id value in the clayyimain nodal user data (integer)
int maxglobXival = -9999; //maximum id value in the globXiraw nodal user data (integer)
int minglobXival = +9999; //minimum id value in the globXiraw nodal user data (integer)
int nbglobXivals = 0; //number of different globXiraw values, simply based on the max-min interval bounds (nb = max - min + 1)
int maxglobXiEval = -9999; //maximum id value in the globXiraw ELEMENTAL user data (integer)
int minglobXiEval = +9999; //minimum id value in the globXiraw ELEMENTAL user data (integer)

/* RANDOM NUMBER GENERATOR initialization */
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dist(0, 1);

/* Data for monitoring simulation run progress */
time_t realtime_atsimstart;
time_t realtime_nextrefresh;
double rt_prev_remain_pred; //previous predicted remaining delta-time, in seconds (internal)
const time_t realtime_ProgRefrDelta = 30; //progress refresh delta-real-time, in seconds (int)

/* ==================== END OF PLUGIN CUSTOM HEADER =========================== */

static IfmResult OnInitModule (IfmModule pMod)
{
  /*
   * TODO: Add your own initialization code here ...
   */
  return True;
}

static void OnExitModule (IfmModule pMod)
{
  /*
   * TODO: Add your own code here ...
   */
}

static IfmResult OnBeginDocument (IfmDocument pDoc)
{
  if (IfmDocumentVersion (pDoc) < IFM_CURRENT_DOCUMENT_VERSION)
    return false;

  try {
    IfmDocumentSetUserData(pDoc, new Cpaleosea2d(pDoc));
  }
  catch (...) {
    return false;
  }

  return true;
}

static void OnEndDocument (IfmDocument pDoc)
{
  delete Cpaleosea2d::FromHandle(pDoc);
}

///////////////////////////////////////////////////////////////////////////
// Implementation of Cpaleosea2d

// Constructor
Cpaleosea2d::Cpaleosea2d (IfmDocument pDoc)
  : m_pDoc(pDoc)
  , editmenu_nEnum(6)
{
  /*
   * TODO: Add your own code here ...
   */
}

// Destructor
Cpaleosea2d::~Cpaleosea2d ()
{
  /*
   * TODO: Add your own code here ...
   */
}

// Obtaining class instance from document handle
Cpaleosea2d* Cpaleosea2d::FromHandle (IfmDocument pDoc)
{
  return reinterpret_cast<Cpaleosea2d*>(IfmDocumentGetUserData(pDoc));
}

// Callbacks

// Handling of Ctrl-C
static void ctrlc_handler(int sig)
{
	signal(SIGINT, SIG_DFL);
	std::cout << "\n\n*** Control-C encountered! --> Simulation Break initiated... ***\n\n" << std::endl;
	//Tech. note on endl: The only difference is that std::endl flushes the output buffer, and '\n' doesn't.
	//(source and example: http://stackoverflow.com/questions/213907/c-stdendl-vs-n)
	globvar_stopsimul = true;
	Beep(880, 200);
}
//C++ documentation related to this:
// function signal: http://www.cplusplus.com/reference/csignal/signal/
// Basic Input/Output: http://www.cplusplus.com/doc/tutorial/basic_io/


tm decomptimediff(double tdsecs, bool uptoyears)
{
	//Code based on: http://stackoverflow.com/a/19205700/3433903
	int seconds, minutes, hours;
	double days;
	double years = 0.0;
	seconds = floor(tdsecs);
	days = seconds / 86400;
	seconds = seconds % 86400;
	hours = seconds / 3600;
	hours = hours % 3600;
	minutes = seconds / 60;
	minutes = minutes % 60;
	seconds = seconds % 60;
	if (uptoyears) {
		years = floor(days / daysperyear);
		days = floor(days - years*daysperyear);
	}
	else {
		days = floor(days);
	}
	tm dtdout;
	dtdout.tm_sec = seconds;
	dtdout.tm_min = minutes;
	dtdout.tm_hour = hours;
	dtdout.tm_yday = days;
	dtdout.tm_year = years;
	return dtdout;
}

int sprintf_decomptimediff(char* const _Buffer, size_t const _BufferCount, tm dtd, bool uptodays = true)
{
	int retval;
	if (uptodays && dtd.tm_yday > 0) {
		retval = sprintf_s(_Buffer, _BufferCount, "%ddays %dh %02dm %02ds", dtd.tm_yday, dtd.tm_hour, dtd.tm_min, dtd.tm_sec);
	}
	else {
		retval = sprintf_s(_Buffer, _BufferCount, "%dh %02dm %02ds", dtd.tm_hour, dtd.tm_min, dtd.tm_sec);
	}
	return retval;
}

void UpdateSimulRunProgressInfo(IfmDocument pDoc)
{
	time_t curr_realtime;
	time(&curr_realtime);
	double curr_simtime = IfmGetAbsoluteSimulationTime(pDoc);
	double rt_elapsdif = difftime(curr_realtime, realtime_atsimstart); //in real seconds
	double st_elapsdif = curr_simtime - simulstarttime; //in simulated days
	double st_wholedif = simulendtime_apriori - simulstarttime; //in simulated days
	double pcsimdone = 100.0 * st_elapsdif / st_wholedif;
//	double rt_wholedif_pred = st_elapsdif / pcsimdone;
	double rt_remain_pred = (simulendtime_apriori - curr_simtime) / (curr_simtime - simulstarttime) * rt_elapsdif;
	bool dorefresh = curr_realtime > realtime_nextrefresh;
	if (dorefresh) {
		double rt_pred_decrease = rt_prev_remain_pred - rt_remain_pred;
		double reliab_absC = (rt_pred_decrease - realtime_ProgRefrDelta) / realtime_ProgRefrDelta;
		double reliab_relC = -rt_pred_decrease / rt_prev_remain_pred;
//		bool reliab_pred = abs(rt_remain_pred - rt_prev_remain_pred) < (2.0 * realtime_ProgRefrDelta);
		bool is_reliab_abs = abs(reliab_absC) < 1.0;
		bool is_reliab_rel = abs(reliab_relC) < 0.01;
		char txtbuff[180], txtbelap[50], txtbrem[50], txtbrelinfo[50];
		tm rtelap = decomptimediff(rt_elapsdif, false);
		tm rtrem = decomptimediff(rt_remain_pred, false);
		sprintf_decomptimediff(txtbelap, 50, rtelap, true);
		sprintf_decomptimediff(txtbrem, 50, rtrem, true);
		bool showrelrel = !is_reliab_abs && (is_reliab_rel || abs(reliab_absC) > 10.0);
		if (showrelrel) sprintf_s(txtbrelinfo, 50, "rel. %+.1f pc", 100.0*reliab_relC);
		if (!showrelrel) sprintf_s(txtbrelinfo, 50, "abs. %+.0f", round(reliab_absC));
		sprintf_s(txtbuff, 180, "Progress:\n                %.2f pc done; time elapsed: %s; time remaining: ((  ~%s  )) [%s %s].\n ", pcsimdone, txtbelap, txtbrem, is_reliab_abs || is_reliab_rel ? "RELIABLE" : "NOT reliable", txtbrelinfo);
		IfmInfo(pDoc, txtbuff);
		realtime_nextrefresh = curr_realtime + realtime_ProgRefrDelta;
		rt_prev_remain_pred = rt_remain_pred;
	}
}
void CheckPresenceOfBCCdata(IfmDocument pDoc, bool neutralize = false)
{
	int i;
	int cntflowt1bcc, cntmasst1bcc;
	cntflowt1bcc = 0;
	cntmasst1bcc = 0;
	for (i = 0; i < nnodes; i++) {
		if (IfmGetBccFlowType(pDoc, i) != IfmBCC_NONE) {
			cntflowt1bcc++;
			if(neutralize) IfmSetBccFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBCC_NONE, 0, 0, IfmMIN_BCC_TYPE);
		}
		if (domasstransport && (IfmGetBccMassType(pDoc, i) != IfmBCC_NONE)) {
			cntmasst1bcc++;
			if (neutralize) IfmSetBccMassTypeAndValueAtCurrentTime(pDoc, i, IfmBCC_NONE, 0, 0, IfmMIN_BCC_TYPE);
		}
	}

	char txtbuffer[180];
	if (cntflowt1bcc > 0) {
		sprintf_s(txtbuffer, 180, "[PreSimul] Beware: %d nodes with non-null BCC for Fluid Flow where detected%s", cntflowt1bcc, neutralize ? " and neutralized." : ". Please verify if this is adequate.");
		IfmWarning(pDoc, txtbuffer);
	}
	if (cntmasst1bcc > 0) {
		sprintf_s(txtbuffer, 180, "[PreSimul] Beware: %d nodes with non-null BCC for Mass Transport where detected%s", cntmasst1bcc, neutralize ? " and neutralized." : ". Please verify if this is adequate.");
		IfmWarning(pDoc, txtbuffer);
	}
}

const bool updateClayConstantKaswell = true; //NEW option so that the model+plugin relies no more on the initial K values of likely deactivated clay elements

//Function to call ONLY in the withclaymesh mode
void UpdateClayParamsAccToYi(IfmDocument pDoc, bool updateAll, double currtime_inUC)
{
	int i;
	char txtbuffer[180];

	int xiindex; //clayXi RELATIVE value, with minclayXival as the reference, and starting at 0; internal 1d-array index to manage special models where minclayXival > 0
	int xival; //clayXi ABSOLUTE value, as specified in the nodal user data
	int gxival; //NEW
	int gxiindex;

	/* Node related scalars and arrays */
	bool somexi, someyimain, someyisub;
	bool yimsprevactiv; //indicates if the node was included in the active part of the clay aquitard with the previous clay-top config.; by comparison of nodal yi value; note that currently active nodes are not excluded, so that this boolean array is not a sufficient condition for selecting exclusively "nodes that were previously active and are no more"...
	bool yimscurractiv; //indicates if the node are included in the active part of the clay aquitard with the current clay-top config.; by comparison of nodal yi value.
	bool prevtopNmainonly; //identifies nodes of the previous clay-top (when clayYi is changing)
	bool belowtopNmain, currtopNmain, overtopNmain, belowtopNsub, currtopNsub, overtopNsub; //Various bool arrays related to clay nodes
	bool steadytoprockN; //identifies nodes at the top of the unconfined (~outcropping) rock aquifer (i.e. no overlying clay)

	/* Element related scalars */
	bool someyiE, belowtopE, currtopE, overtopE;
	int lyival; //layer index in the clay aquitard, for current element (from base to top)

	/* Updating these important global variables */
	hasclayYichanged = currclayYi != prevclayYi;
	hasclayYiincreased = currclayYi > prevclayYi;

	//1. Updating active state of clay elements (if clayYi has changed)
	if (hasclayYichanged || updateAll) {
		for (i = 0; i < nelems; i++) {
			lyival = Pclaylayeri[i];
			someyiE = (lyival >= 1);
			belowtopE = (someyiE && (currclayYi > 1)) ? (lyival < currclayYi) : false;
			currtopE = (someyiE && (currclayYi > 0)) ? (lyival == currclayYi) : false;
			overtopE = someyiE ? (lyival > currclayYi) : false;
			//then, the 'active or not' updates (ONLY in the clay facies, where someyiE == true)
			if (belowtopE || currtopE || overtopE) {
				IfmSetMatElementActive(pDoc, i, belowtopE || currtopE ? 1 : 0);
				PelemcurractiveE[i] = belowtopE || currtopE;
				// if belowtopE --> 1 : activation of layers below current top of clays (optional: NO MORE, I think)
				// if currtopE --> 1 : activation of elements of the new top layer of clays
				// if overtopE --> 0 : deactivation of layers above current top of clays (optional: NO MORE, I think)
			}
		}
		//TRYING REMOVED: IfmUpdateElementExtents(pDoc); //BETA: Is this call adequate and useful?
	}

	//1b. Updating K of clays (at every time step, but only in some sub-modes) --- NEW, BETA!!!
	currKofclays = withtimedepKclays ? IfmInterpolatePowerValue(pDoc, clayeqkvPID, currtime_inUC) : K_CLAYS;
	if ((currclayYi > 0) && (withtimedepKclays || (updateClayConstantKaswell && (hasclayYichanged || updateAll)))) {
		for (i = 0; i < nelems; i++) {
			lyival = Pclaylayeri[i];
			someyiE = (lyival >= 1);
			belowtopE = (someyiE && (currclayYi > 1)) ? (lyival < currclayYi) : false;
			currtopE = (someyiE && (currclayYi > 0)) ? (lyival == currclayYi) : false;
			//NOT.USED: overtopE = someyiE ? (lyival > currclayYi) : false;
			if (belowtopE || currtopE) IfmSetMatConductivityValue2D(pDoc, i, currKofclays);
		}
	}

	//1c. [BETA!!!] Updating Source terms in clays (at every time step, but only in some sub-modes) --- NEW, BETA!!!
	currSrcGlobDarcyflux = BETAwithsourceterms ? IfmInterpolatePowerValue(pDoc, clayqdbotPID, currtime_inUC) : nodata_double;

	bool isthecurrentGlobDfluxsignif = !isnan(currSrcGlobDarcyflux) && (currSrcGlobDarcyflux > 1e-16); //i.e. is it >0.0 ± a small numerical error
	if (!isthecurrentGlobDfluxsignif) currSrcGlobDarcyflux = 0.0; //forcing a real zero value in such case!
	//This isthecurrentGlobDfluxsignif verification allows to optimize the later part of the simulation by skipping the
	// update if it is indeed not necessary. [BETA!]

	//REMOVED SINCE IT IS MOST LIKELY NOT NECESSARY!: if (hasclayYichanged) areSrcTermsAlreadyZeroed = false; //safety in case of newly activated clay elements when Dflux still is null (== 0.0 m/d)
	// \-> That is because the assignation for-loop is carried out on every element with a value of PclaysrcBetaCoefE even if for the elements that are not yet activated!

	bool doSrcTermsRequireAnUpdate = isthecurrentGlobDfluxsignif || !areSrcTermsAlreadyZeroed;

	bool someSrcE;
	double volumelemi; //volume of the current element (in the for loop) (in m3); or 0.0 m3 if the element is currently inactive (TODO improve some day... see comment below)
	double computedQf; //source term value for fluid
	double computedQm; //source term value for mass
	double infototFluidSrc = 0.0; //total fluid income from these distributed sources (in m3/d)
	double infototMassSrc = 0.0; //total mass income from these distributed sources (in g/d)
	int cntSrc = 0; //count of the number of elements with updated fluid & mass source terms for the current time step
	if (BETAwithsourceterms && doSrcTermsRequireAnUpdate) {
		for (i = 0; i < nelems; i++) {
			//lyival = Pclaylayeri[i];
			//someyiE = (lyival >= 1);
			someSrcE = PclaysrcBetaCoefE[i] > 0.0;
			if (someSrcE) {
				volumelemi = IfmGetElementalContent(pDoc, IfmTOTAL_VOLUME, i); //TODO: Could be optimized by reading these info only once (as is already done in the ExportFF module...); still, it should be safe even when the element is inactive (then equals 0.0 m3)
				//REPLACED (Oct.2017): computedQf = (currSrcQdarcy * PclaysrcBetaCoefE[i]) / volumelemi; //HERE ASSUMING the fluid in injected in one single element per column of the model.
				computedQf = currSrcGlobDarcyflux * PclaysrcBetaCoefE[i]; //NEW BETA formula (Oct.2017)
				computedQm = computedQf * swconc_atstart;
				IfmSetMatFlowSinkSource(pDoc, i, computedQf);
				IfmSetMatMassSinkSource(pDoc, i, computedQm);
				infototFluidSrc = infototFluidSrc + computedQf*volumelemi;
				infototMassSrc = infototMassSrc + computedQm*volumelemi;
				cntSrc++;
			}
		}
		sprintf_s(txtbuffer, 180, "[WCM @ UpdateClayParamsAccToYi] SRC TERMS: updated for n=%d elements, based on currSrcGlobDarcyflux = %g m/d.", cntSrc, currSrcGlobDarcyflux);
		IfmInfo(pDoc, txtbuffer);
		sprintf_s(txtbuffer, 180, "[WCM @ UpdateClayParamsAccToYi] SRC TERMS: Expected total gain from distributed sources: of fluid = %g m3/d; of mass = %g g/d", infototFluidSrc, infototMassSrc);
		IfmInfo(pDoc, txtbuffer);

		areSrcTermsAlreadyZeroed = !isthecurrentGlobDfluxsignif;
	}

	//TEMP BETA msg:
	if (BETAwithsourceterms && !doSrcTermsRequireAnUpdate) {
		sprintf_s(txtbuffer, 180, "[WCM @ UpdateClayParamsAccToYi] SRC TERMS: NO update was made for the current time step. In most cases this is because the Dflux is negligible (current = %g m/d).", currSrcGlobDarcyflux);
		IfmInfo(pDoc, txtbuffer);
	}

	//2. Updating nodal values of plugin internal arrays (in prep. for later update of nodal BC's and states)
	if (hasclayYichanged || updateAll) {
		ntopnodes = 0; //initialization before counting such nodes within this loop
		ncurractivenodes = 0; //initialization before counting such nodes within this loop

		//Variables for optional export of debugging info as nodal User Data
		long debugchgclayRDId = IfmGetNodalRefDistrIdByName(pDoc, "debug_chgclay");
		bool beta_debugchgclay = debugchgclayRDId >= 0;
		double debugchgclay_vali;

		for (i = 0; i < nnodes; i++) {
			xival = PclayXiraw[i];
			somexi = (xival >= 0);
			if (somexi) xiindex = xival - minclayXival;

			someyimain = (PclayYimain[i] >= 0);
			someyisub = (PclayYisub[i] >= 0);
			yimsprevactiv = (someyimain && PclayYimain[i] <= prevclayYi) || (someyisub && PclayYisub[i] <= prevclayYi);
			yimscurractiv = (someyimain && PclayYimain[i] <= currclayYi) || (someyisub && PclayYisub[i] <= currclayYi);
			prevtopNmainonly = someyimain ? (PclayYimain[i] == prevclayYi) : false;
			belowtopNmain = someyimain ? (PclayYimain[i] < currclayYi) : false;
			belowtopNsub = someyisub ? (PclayYisub[i] < currclayYi) : false;
			currtopNmain = someyimain ? (PclayYimain[i] == currclayYi) : false;
			currtopNsub = someyisub ? (PclayYisub[i] == currclayYi) : false;
			overtopNmain = someyimain ? (PclayYimain[i] > currclayYi) : false;
			overtopNsub = someyisub ? (PclayYisub[i] > currclayYi) : false;
			steadytoprockN = (PclayYimain[i] == rockYi_cstval);

			//building the array with indeces of the nodes forming the new current top of clays
			if (currtopNmain) Pprevclaytopnodes[xiindex] = Pcurrclaytopnodes[xiindex]; //temp. storage of previous top nodes indeces
			if (currtopNmain) Pcurrclaytopnodes[xiindex] = i;
			if (currtopNmain && !somexi) IfmWarning(pDoc, "[WCM @ UpdateClayParamsAccToYi] Unexpected xi value at a currtopNmain node: somexi == false. Why?");

			//arrays with interpreted transient data (will change at several steps in time, according to clayyi TS)
			PtopclayBCnodes[i] = currtopNmain;
			PtoprockBCnodes[i] = steadytoprockN;
			Ptopnodes[i] = currtopNmain | steadytoprockN; //top of either clay aquitard or exposed rock aquifer; fully redefined for all nodes here

			//OLD VERSIONS:
			//PupdateNbelowclaytop[i] = !updateAll ? ((hasclayYichanged && !hasclayYichglower) ? currtopNsub : false) : (currtopNsub | belowtopNmain | belowtopNsub);

			PupdateNbelowclaytop[i] = (currtopNsub || belowtopNmain || belowtopNsub) &&
				(updateAll || hasclayYichanged && (!yimsprevactiv && yimscurractiv));
			PupdateNoverclaytop[i] = (overtopNmain || overtopNsub) &&
				(updateAll || hasclayYichanged && (yimsprevactiv && !yimscurractiv));
			PcleanNBCprevclaytop[i] = prevtopNmainonly && !PupdateNbelowclaytop[i] && !PupdateNoverclaytop[i] &&
				hasclayYichanged;
			//so that nodes constituting the previous clay-top will require such cleaning (i.e. BCs removed) only if they have not already been treated by the two other boolean criteria...

			//Other general node classifiers
			Pallcurrclayactivenodes[i] = (currtopNmain || currtopNsub || belowtopNmain || belowtopNsub);
			Pallcurractivenodes[i] = Pallcurrclayactivenodes[i] || steadytoprockN || Psteadyallrocknodes[i];
			if (Pallcurractivenodes[i]) ncurractivenodes++;

			//Optional export for debugging (effective only if 'debug_chgclay' nodal User Data is present)
			if (beta_debugchgclay && updateAll) {
				debugchgclay_vali = PupdateNbelowclaytop[i] ? -1 : (PupdateNoverclaytop[i] ? 1 : 0);
				if (PcleanNBCprevclaytop[i]) debugchgclay_vali = -2;
				IfmSetNodalRefDistrValue(pDoc, debugchgclayRDId, i, debugchgclay_vali); //(optional) User Data update in the simulation results
			}

			//updating Pcurrglobtopnodes (if optional export is active)
			if (isglobXipresent && Ptopnodes[i]) {
				gxival = PglobXiraw[i];
				if (gxival >= 0) {
					gxiindex = gxival - minglobXival;
					Pcurrglobtopnodes[gxiindex] = i;
					if (altexportautoinclgtop) PeffexportnodesA[gxiindex] = i; //should be safe... (BETA)
				}
				else {
					sprintf_s(txtbuffer, 180, "[WCM @ UpdateClayParamsAccToYi] Unexpected globXiraw (%d < 0) at a Ptopnodes[i=%d] == true node. Why?", gxival, i);
					IfmWarning(pDoc, txtbuffer);
				}
			}

			//updating summary variables (length 1)
			if (Ptopnodes[i]) ntopnodes++;
		}
	}

	//3. Counting of nodes requiring (re-)initialization --> needupdatearoundclaytop (t|f)
	if (hasclayYichanged || updateAll) {
		int cnttot = 0; //total count of nodes requiring update of fluid&masstransp states and/or BCs
		int cntbelow = 0, cntover = 0, cntprevtop = 0;
		for (i = 0; i < nnodes; i++) {
			if (PupdateNbelowclaytop[i]) cntbelow++;
			if (PupdateNoverclaytop[i]) cntover++;
			if (PcleanNBCprevclaytop[i]) cntprevtop++;
		}
		cnttot = cntbelow + cntover + cntprevtop;
		needupdatearoundclaytop = cnttot > 0;
		if (needupdatearoundclaytop) {
			sprintf_s(txtbuffer, 180, "[WCM @ UpdateClayParamsAccToYi] %d nodes requiring state and/or BC update around current clay top (below=%d, over=%d, previous.top=%d).", cnttot, cntbelow, cntover, cntprevtop);
			IfmInfo(pDoc, txtbuffer);
		}
	}
	else needupdatearoundclaytop = false;
}

//Function to call whatever the plugin MODE; it will adapt depending on withclaymesh et al.
void UpdateInflowDataForTopNodes(IfmDocument pDoc, bool initial = false)
{
	int i;
	int xiindex; //clayXi RELATIVE value, with minclayXival as the reference, and starting at 0; internal 1d-array index to manage special models where minclayXival > 0
	int xival; //clayXi ABSOLUTE value, as specified in the nodal user data
	int gxiindex;
	bool somexi;
	if (initial) {
		//INITIAL mode:
		// >>> Defines a starting inflow state for all nodes: only clay and rock top nodes may have true values, while all non-top nodes are set to false.
		for (i = 0; i < nnodes; i++) {
			Pinflownodes[i] = Ptopnodes[i] ? (firstTSresumed ? IfmGetBcMassType(pDoc, i) == IfmBC_DIRICHLET : true) : false;
			//BETA TODO Review if line above is compatible with the new withclaymesh mode?...
		}
		//NEW [BETA]: initializing globtopnodes-aligned 1D array for the inflow state of current globtop nodes:
		if (isglobXipresent) {
			int nodei = -1; //nodal index in the FULL 1D array of nodes (i.e. different than the globtop-aligned array)
			for (gxiindex = 0; gxiindex < nbglobXivals; gxiindex++) {
				nodei = Pcurrglobtopnodes[gxiindex];
				Pcurrglobtop_isinflow[gxiindex] = Pinflownodes[nodei];
			}
			//Initialization of the counters: = indeterminate value (yet should not appear in log during the simulation)
			ninflownodes_gtop = -1; 
			noutflownodes_gtop = -1;
			nflowchgnodes_gtop = -1;
		}
	}
	if(withclaymesh && (hasclayYichanged || initial)) {
		//Mainly for each time step where number of clay layers is changed (in the withclaymesh mode, of course):
		// >>> Transfers the last inflow state from the previous top node to the current top node at this same Xi location
		// (not relevant in NON-withclaymesh modes, as in those modes the top nodes remain the same throughout the simulation)
		//Difference between the initial vs. normal modes in there:
		// initial: topnodes are not changing, so that we can get the inflow info directly from nodes Ptopnodes[i] == true and store it in the special array for current-clay-top nodes. ~Initialization of Pcurrclaytopisinflow array.
		// normal: topnodes have just changed, so that we must instead get the most recent values for the previous top nodes, to ensure continuity in in/outflow states at the top when clayYi changes. ~Update of Pcurrclaytopisinflow array.
		//Conversely, this IF is skipped when Yi has not changed, because in such cases, Pinflownodes[i] is already up-to-date and correctly aligned/assigned to the current top nodes, including the clay.
		for (i = 0; i < nnodes; i++) {
			if (Ptopnodes[i]) {
				xival = PclayXiraw[i];
				somexi = (xival >= 0); //so that top-nodes not in clay (so in rock) are autom. excluded hereafter:
				if (somexi) xiindex = xival - minclayXival;
				if (somexi && initial) Pcurrclaytopisinflow[xiindex] = Pinflownodes[i]; //initial mode
				if (somexi && !initial) Pinflownodes[i] = Pcurrclaytopisinflow[xiindex]; //normal mode
				//DON'T DO THIS WARNING, because there are also the top nodes of rock!:
				// if (!somexi) IfmWarning(pDoc, "Strangely a topnode has an undefined xival!?...");
			}
			else Pinflownodes[i] = false;
		}
	}
	//Updating the count of inflowing nodes at top nodes
	int tmpcntinfl = 0; //counter of inflow nodes, for display (from the previous flow solution)
	int tmpcntmbcc = 0; //counter of topnodes with mass BCc, for display (BETA; likely constant count?)
	for (i = 0; i < nnodes; i++) { //separate loop quite not optimal but I want it executed before IfmInfo
		if (!Ptopnodes[i]) continue;
		if (Pinflownodes[i]) tmpcntinfl++;
		if (BETAautomUseMassBCconstraints) {
			if (IfmIsBccMassSet(pDoc, i, IfmMIN_BCC_TYPE) == True) tmpcntmbcc++;
		}
	}
	ninflownodes = tmpcntinfl;
	noutflownodes = ntopnodes - ninflownodes;
	nbcctopnodes = tmpcntmbcc;

	//NEW: Updating the count of inflowing nodes at GLOBAL-top nodes
	int tmpcntchg = 0;
	tmpcntinfl = 0;
	if (isglobXipresent) {
		for (gxiindex = 0; gxiindex < nbglobXivals; gxiindex++) {
			if (Pcurrglobtop_isinflow[gxiindex]) tmpcntinfl++;
			if (Pcurrglobtop_isinflow[gxiindex] != Pprevglobtop_isinflow[gxiindex]) tmpcntchg++;
		}
		ninflownodes_gtop = tmpcntinfl;
		noutflownodes_gtop = nbglobXivals - tmpcntinfl;
		nflowchgnodes_gtop = tmpcntchg;
	}
	else {
		ninflownodes_gtop = -1;
		noutflownodes_gtop = -1;
		nflowchgnodes_gtop = -1;
	}
}


void Cpaleosea2d::Serialize (IfmDocument pDoc, IfmArchive pArc)
{
  /*
   * TODO: Add your own code here ...
   */
}

void CreateUserData_forPaleoSea2D(IfmDocument pDoc)
{
	std::string nodals[11] = { "nodefacies", "hdensadd", "Yoftoprock", "globxiraw", "clayyimain", "clayyisub", "clayxiraw", "rndnormN", "hsteady0", "cinit0noisy", "exportnodes" };
	const int nnodal = 11;

	std::string elemtals[7] = { "facies", "layerglob", "claylayeri", "toprockYatmid", "rndnormE", "contentzones", "monitzones" };
	const int nelemtal = 7;

	int i;
	int rd_id = -1; //current index for the found or not found nodal/elem. ref. distrib. (@User Data)
	bool rdexists = false;
	char txtbuff[200];

	//Nodal
	for (i = 0; i < nnodal; i++) {
		rd_id = IfmGetNodalRefDistrIdByName(pDoc, nodals[i].c_str());
		rdexists = rd_id >= 0;
		if (!rdexists) {
			IfmCreateNodalRefDistr(pDoc, nodals[i].c_str());
			sprintf_s(txtbuff, 200, "PaleoSea2D model preparation :: nodal user data named '%s' was created.", nodals[i].c_str());
			IfmInfo(pDoc, txtbuff);
		}
		else {
			sprintf_s(txtbuff, 200, "PaleoSea2D model preparation :: nodal user data named '%s' already exists.", nodals[i].c_str());
			IfmInfo(pDoc, txtbuff);
		}
	}

	//Elemental
	for (i = 0; i < nelemtal; i++) {
		rd_id = IfmGetElementalRefDistrIdByName(pDoc, elemtals[i].c_str());
		rdexists = rd_id >= 0;
		if (!rdexists) {
			IfmCreateElementalRefDistr(pDoc, elemtals[i].c_str());
			sprintf_s(txtbuff, 200, "PaleoSea2D model preparation :: elemental user data named '%s' was created.", elemtals[i].c_str());
			IfmInfo(pDoc, txtbuff);
		}
		else {
			sprintf_s(txtbuff, 200, "PaleoSea2D model preparation :: elemental user data named '%s' already exists.", elemtals[i].c_str());
			IfmInfo(pDoc, txtbuff);
		}
	}
}

void Cpaleosea2d::OnEditDocument (IfmDocument pDoc, Widget wParent)
{
	static IfmPropExtEnum enums[] = {
		{ "Task 1: CREATE all required User Data", 1 },
		{ "DO NOTHING", -1 },
		{ 0 },
	};
	static int idef = -1;
	static const char szDesc2[] = "Choice of tasks that can be run to facilitate the preparation of the FEFLOW model, which later will be simulated using the PaleoSea2D plugin.";

	IfmProperty props[] = {
		{ "Task to do", IfmPROP_ENUM, &editmenu_nEnum, &idef, szDesc2, enums },
	};

	//force default value (Do Nothing) every time for the enum property:
	editmenu_nEnum = idef;

	IfmEditProperties(pDoc, "MENU for the PaleoSea2D plugin (c) 2018 Marc Laurencelle", "", props, countof(props));
	if (editmenu_nEnum == 1) {
		IfmInfo(pDoc, "TASK #1 is starting...");
		CreateUserData_forPaleoSea2D(pDoc);
		IfmInfo(pDoc, "TASK #1 has completed successfully.");
	}
}

ncdfheader ncdf_globtopN;
const IfmContentType MyActiveContentTypes[4] = { IfmDILUTED_MASS,IfmFLUID_CONTENT,IfmVOID_VOLUME,IfmTOTAL_VOLUME };

int initialize_ncdf_exports(IfmDocument pDoc)
{
	if (!doexportdata) return 0;
	IfmInfo(pDoc, "[Export init.] Running initialize_ncdf_exports (once)...");
	int ncerr = 0;

	//BETA VERSION FOR the new User ExportNodes mode:
	if (altexportusernodes) {
		ncerr = netcdf_specify_export_nodes(pDoc, ncdf_globtopN, PeffexportnodesA, neffexportnodes, !altexportautoinclgtop, Peffexportnodes_uvaluesA, true);
	}
	else {
		//NOTE: Normally, nbglobXivals should be identical (equal) to ntopnodes.
		ncerr = netcdf_specify_export_nodes(pDoc, ncdf_globtopN, Pcurrglobtopnodes, nbglobXivals, !withclaymesh, nullptr, true);
	}
	if (ncerr) ERRpropag(ncerr); //Will quit the function if there was an error (returning the err.code)
	
	/* Explicit initializations of the ncdf structured object */
	ncdf_globtopN.nbelems_fixed = IfmGetNumberOfElements(pDoc); //must be set prior to the calls below (TODO minor: make it clearer in the code and headers!)
	ncdf_globtopN.nbnodes_fixed = IfmGetNumberOfNodes(pDoc); //same (related to above comment)
	ncdf_globtopN.putNodeFF = true;
	ncdf_globtopN.putX = true;
	ncdf_globtopN.putY = true;
	ncdf_globtopN.putYoftoprock = true;
	ncdf_globtopN.puthead = true;
	ncdf_globtopN.putconc = true && domasstransport;
	ncdf_globtopN.putDfluxcomponents = true;
	ncdf_globtopN.putmratebudg = true && domasstransport && EXPORT_MASS_RBUDGETs;
	ncdf_globtopN.PtocompmrbudgetA = domasstransport && EXPORT_MASS_RBUDGETs ? Pnodalmassrbudgets : nullptr; //BETA!
	ncdf_globtopN.steadyActiveElems = !withclaymesh;
	ncdf_globtopN.extern_transientACTIVelemA = withclaymesh ? PelemcurractiveE : nullptr;
	ncdf_globtopN.putconcquantiles = true && domasstransport; //must be set prior to the call below (TODO minor: make it clearer in the code and headers!)
	ncdf_globtopN.putmassconvectmonit = true && domasstransport;
	ncdf_globtopN.putnbnegconcs = true && domasstransport;
	ncdf_globtopN.putglobextrconcs = true && domasstransport;
	ncdf_globtopN.putvconcgrads = true && domasstransport;
	/* Calling specialized initialization functions, then */
	/* 1. for content zone (CZ) monitoring */
	if (!netcdftools_init_contentzones(pDoc, ncdf_globtopN, MyActiveContentTypes, sizeof(MyActiveContentTypes) / sizeof(MyActiveContentTypes[0]), "contentzones")) { //BETA! sizeof trick based on http://stackoverflow.com/questions/4108313/how-do-i-find-the-length-of-an-array
		IfmWarning(pDoc, "TODO-BETA: An error msg+stop should be coded here: when calling netcdftools_init_contentzones...");
	}
	/* 2a. for the shared steady data arrays (effective porosity for mass; total elem. volume) */
	netcdftools_init_steadydata(pDoc, ncdf_globtopN, false, false); //Note: now autom. adapts to .put active options!
	/* 2b. for the shared dynamic data arrays (BETA... complete this comment LATER) */
	netcdftools_init_dynamicdatabuffers(pDoc, ncdf_globtopN);
	/* 3. for monit zone (MZ) monitoring (for speed of convective fingers) */
	if (ncdf_globtopN.putmassconvectmonit) {
		ncdf_globtopN.ref_C0 = minconc;
		ncdf_globtopN.ref_CS = maxconc;
		ncdf_globtopN.DPF_normconcmin = 0.01;
		if (!netcdftools_init_monitzones(pDoc, ncdf_globtopN, "monitzones", -1, "toprockYatmid", -1)) {
			IfmWarning(pDoc, "TODO-BETA: A warning msg should be coded here: when calling netcdftools_init_monitzones without successfully loading and preparing everything...");
		}
	}

	char timestamp[100];
	//std::string timestamp_str;
	time_t timer; //normally 64-bit compatible
	struct tm tm_info; //stores time decomposed as year, month, day, hour, min, etc.
	time(&timer); //normally 64-bit compatible (calls _time64())
	errno_t err_time = localtime_s(&tm_info, &timer); //normally 64-bit compatible (calls _localtime64_s())
	if (err_time == 0) {
		strftime(timestamp, 100, "_%Y%m%d_%H%M%S", &tm_info);
	}
	else {
		clock_t now = clock(); //somewhat to a gettickcount, but more compatible
		sprintf_s(timestamp, 100, "_tc%d", now);
		IfmWarning(pDoc, "[Export init.] Could not build a stamp from local time; will use clock ~tick count instead.");
	}
	//timestamp_str = timestamp;

	problemfilename = IfmGetProblemTitle(pDoc); //somefile.fem
	exportfilename = IfmForceExtension(pDoc, problemfilename.c_str(), ""); //somefile [no extension yet!]
	resultsfolderpath = IfmGetAbsolutePath(pDoc, "../results/", IfmGetFileDirectory(pDoc, IfmGetProblemPath(pDoc)));
	exportfilepath = resultsfolderpath + exportfilename + timestamp + ".nc";

	strcpy_s(ncdf_globtopN.filepath, _MAX_PATH, exportfilepath.c_str());
	
	char txtbuff[_MAX_PATH + 50];

	sprintf_s(txtbuff, _MAX_PATH + 50, "[Export init.] Creating binary file <%s>...", exportfilepath.c_str());

	IfmWarning(pDoc, txtbuff);

	if ((ncerr = netcdf_create_file_headers(pDoc, ncdf_globtopN))) ERRpropag(ncerr);
	
	return NC_NOERR; //Returns NC_NOERR (0) if everything is okay; returns sooner in code otherwise (with a value != 0).
}


void Cpaleosea2d::PreSimulation (IfmDocument pDoc)
{
	globvar_stopsimul = false;
	signal(SIGINT, ctrlc_handler);

	IfmInfo(pDoc, "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =");
	IfmInfo(pDoc, "  PaleoSea2D plugin for simulating a marine transgression using time-dependent BC's  ");
	IfmInfo(pDoc, "                             (c) Marc Laurencelle, 2018                              ");
	IfmInfo(pDoc, "<< preSimulation begins >>");

	meshquadtype = true; //IfmGetMeshType(pDoc) == IfmMSH_QUAD ? true : false; This API function still does not work with the FEFLOW 7.1.003 SDK!
	//TODO-minor-later: Contact FEFLOW Support to resolve this issue from their side! (IfmGetMeshType)
	IfmInfo(pDoc, "[PreSimul] Internal param. meshquadtype = true is FORCED (hard-wired). The plugin is indeed incompatible with triangular meshes.");

	char txtbuffer[180];
	int i;
	double nodval, elemval;
	nnodes = IfmGetNumberOfNodes(pDoc);
	nelems = IfmGetNumberOfElements(pDoc);
	
	if (IfmGetProblemType(pDoc) != IfmTYPE_SATURATED) {
		IfmAlert(pDoc, NULL, "  OK  ", "The plugin works only for in Saturated Mode (as problem type).");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (IfmGetProblemProjection(pDoc) != IfmPROJ_VERTICAL_2D) {
		IfmAlert(pDoc, NULL, "  OK  ", "The plugin requires a 2D vertical projection (as proj. type).");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (!(IfmGetProblemClass(pDoc) == IfmPCLS_FLOW || IfmGetProblemClass(pDoc) == IfmPCLS_MASS_TRANSPORT)) {
		IfmAlert(pDoc, NULL, "  OK  ", "The plugin works only for models involving flow and/or mass transport (as problem class).");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (IfmGetTimeClass(pDoc) != IfmTCLS_UNSTEADY) {
		IfmAlert(pDoc, NULL, "  OK  ", "The plugin works only for transient problems (i.e. 'Unsteady' time class).");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	domasstransport = IfmGetProblemClass(pDoc) == IfmPCLS_MASS_TRANSPORT;
	IfmInfo(pDoc, domasstransport ? "[PreSimul] Normal plugin mode: coupled fluid flow & mass transport + density dependence." : "[PreSimul] Simplistic plugin mode: transient fluid flow only (without any density effect).");

	//NEW feature
	if (domasstransport) {
		divergformtransp = IfmIsDivergenceFormTransport(pDoc) == True;
		if (divergformtransp) {
			IfmWarning(pDoc, "* BETA! [PreSimul] Detected 'Divergence form' mode for the transport equation: >> Transient outflowing nodes will have no mass t. BC(t), and thus purely-advective, zero-dispersive, mass fluxes!");
		}
		else {
			IfmInfo(pDoc, "[PreSimul] Default 'Convective form' mode for the transport equation (i.e. transient outflowing nodes with purely-advective, zero-dispersive, mass fluxes, as commonly done).");
		}

		//BETA-TESTING unsafe feature
		if (BETApermanentDirichletMassConcBCs) {
			IfmWarning(pDoc, "* BETA! [PreSimul] BETApermanentDirichletMassConcBCs == true hard-wired in the plugin. Careful [!!!]");
			IfmAlert(pDoc, NULL, "  OK  ", "Beware: BETApermanentDirichletMassConcBCs is active in the plugin DLL!");
		}

		//BETA-TESTING unsafe feature --- TODO: This mode should be completely removed! (said on Dec. 2nd)
		if (BETAautomUseMassBCconstraints) {
			IfmWarning(pDoc, "* BETA! [PreSimul] BETAautomUseMassBCconstraints == true hard-wired in the plugin. Careful [!!!]");
			IfmWarning(pDoc, "** [!! BETAautomUseMassBCconstraints] (TRYING initial Bcc setting only with UNCHANGED conc.value(XXtXX) afterwards -- 26 Oct.)");
			IfmAlert(pDoc, NULL, "  OK  ", "Beware: BETAautomUseMassBCconstraints is active in the plugin DLL!");
		}
	}

	simulstarttime = IfmGetAbsoluteSimulationTime(pDoc);
	simulendtime_apriori = IfmGetFinalSimulationTime(pDoc);

	std::string PCname = "";
	rslPID = -1;
	swcPID = -1;
	bclaycPID = -1;
	maxcPID = -1;
	sedratePID = -1;
	claytopPID = -1;
	seddurPID = -1;
	rndPID = -1;
	clayYiPID = -1;
	clayeqkvPID = -1;
	clayqdbotPID = -1;
	i = IfmGetPowerCurve(pDoc, 0);
	while (i>0) {
		PCname = IfmGetPowerComment(pDoc, i);
		/* compulsory time series */
		if (PCname == "rsl") rslPID = i;
		/* recommended t.s. (else using a def. const. value) */
		if (PCname == "swconc") swcPID = i;
		if (PCname == "maxconc") maxcPID = i;
		if (PCname == "minconc") mincPID = i;
		if (PCname == "fwconc") fwconcPID = i;
		if (PCname == "rndconc") rndPID = i;
		/* mode-specific t.s. */
		if (PCname == "botclayconc") bclaycPID = i;
		if (PCname == "sedrate") sedratePID = i; //old, almost deprecated, method for clay sedimentation
		if (PCname == "claytop") claytopPID = i;
		if (PCname == "sedimduring") seddurPID = i;
		if (PCname == "clayyi") clayYiPID = i;
		if (PCname == "clayeqkv") clayeqkvPID = i;
		if (PCname == "clayyifixmbc") clayYiFixMBCPID = i;
		if (PCname == "clayqdbot") clayqdbotPID = i;
		//IfmInfo(pDoc, PCname.c_str());
		i = IfmGetPowerCurve(pDoc, i);
	}
	
	/* Verifying is time series requirements are met, or else simulation may be aborted. */
	if (rslPID < 0) {
		IfmAlert(pDoc, NULL, "  OK  ", "The plugin requires a time series named rsl.");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}
	if (domasstransport) {
		defswconc = 35000.0;
		deffwconc = 250.0;
		/* checking TS 'swconc' (VARIABLE or constant) */
		withtransientswconc = swcPID >= 0;
		if (withtransientswconc) {
			swconc_atstart = IfmInterpolatePowerValue(pDoc, swcPID, simulstarttime); //NEW BETA
			IfmInfo(pDoc, "[PreSimul] Solute concentration (~TDS salinity) of flooding waters will be based on 'swconc' time series.");
		}
		else {
			sprintf_s(txtbuffer, 180, "[PreSimul] Detected no time series named 'swconc'. The default constant value (%.3f g/L) will be used.", defswconc / 1000.0);
			IfmWarning(pDoc, txtbuffer);
			IfmWarning(pDoc, "[PreSimul] ...and hence, FLOODED & inflowing nodes at model surface will have C = cst.swconc whatever the time!");
		}
		/* setting maxconc (~Cs) from the TS (CONSTANT) */
		if (maxcPID >= 0) {
			//With a user-specified setting for maxconc:
			maxconc = IfmInterpolatePowerValue(pDoc, maxcPID, simulstarttime);
			//An IfmInfo writes down that value to the Log later in the code...
		}
		else {
			//In absence of user setting: maxconc is set to the default hard-wired value:
			maxconc = defswconc;
			sprintf_s(txtbuffer, 180, "[PreSimul] Detected no time series named 'maxconc'. The default constant value (%.3f g/L) will be used.", maxconc / 1000.0);
			IfmWarning(pDoc, txtbuffer);
		}
		/* setting fwconc from the TS (CONSTANT) */
		if (fwconcPID >= 0) {
			fwconc = IfmInterpolatePowerValue(pDoc, fwconcPID, simulstarttime);
			sprintf_s(txtbuffer, 180, "[PreSimul] A value of %.3f g/L will be used as the CONSTANT solute concentration for freshwater (based on 'fwconc' TS).", fwconc / 1000.0);
			IfmInfo(pDoc, txtbuffer);
		}
		else {
			fwconc = deffwconc;
			sprintf_s(txtbuffer, 180, "[PreSimul] Detected no time series named 'fwconc'. The default constant value (%.3f g/L) will be used.", fwconc / 1000.0);
			IfmWarning(pDoc, txtbuffer);
		}
		/* setting minconc (~C0) from the TS (CONSTANT) */
		if (mincPID >= 0) {
			minconc = IfmInterpolatePowerValue(pDoc, mincPID, simulstarttime);
			//An IfmInfo writes down that value to the Log later in the code...
		}
		else {
			minconc = fwconc > 0.0 ? fwconc : 0.0; //i.e. minconc = fwconc, or = 0.0 if fwconc is negative (unusual case)
			sprintf_s(txtbuffer, 180, "[PreSimul] Detected no time series named 'minconc'. The default assignation method applies: minconc = fwconc = %.3f g/L (CONSTANT).", minconc / 1000.0);
			IfmWarning(pDoc, txtbuffer);
		}
		/* checking TS 'rndconc' (VARIABLE or constant) */
		withrndnoiseaddtoconcBCs = rndPID >= 0;
		if (!withrndnoiseaddtoconcBCs) {
			sprintf_s(txtbuffer, 180, "[PreSimul] Detected no time series named 'rndconc'. A default duration of %.3f years will be used for random noise addition to fixed conc. BCs.", defrandomduration / daysperyear);
			IfmWarning(pDoc, txtbuffer);
		}
	}
	else {
		defswconc = 0.0;
		deffwconc = 0.0;
		minconc = 0.0; //though facultative...
		maxconc = 0.0; //though facultative...
		withtransientswconc = false;
		withrndnoiseaddtoconcBCs = false;
	}

	//Extraction of the uniform constant value for density ratio (for the full 2Dxz domain)
	if (domasstransport) {
		densr = IfmGetMatFlowDensityRatio(pDoc, 0); //gets the value of element #1 (assumes uniform field)
		if (densr < densr_tol) {
			densr = 0.0; //to prevent from negative values
			IfmWarning(pDoc, "[PreSimul] A density ratio <= 0 was detected (from elem #1): fw-like heads will be applied at topnodes.");
		}
		else {
			/* info on density ratio */
			sprintf_s(txtbuffer, 180, "[PreSimul] A density ratio = %f was detected (from elem #1): saltwater heads will be applied at topnodes.", densr);
			IfmInfo(pDoc, txtbuffer);
			/* info on maxconc (Cs) */
			sprintf_s(txtbuffer, 180, "[PreSimul] A max. concentration ref. value of Cs = %.3f g/L is used for this simulation.", maxconc / 1000.0);
			IfmInfo(pDoc, txtbuffer);
			/* info on minconc (C0) */
			sprintf_s(txtbuffer, 180, "[PreSimul] A min. concentration ref. value of C0 = %.3f g/L is used for this simulation.", minconc / 1000.0);
			IfmInfo(pDoc, txtbuffer);
		};
	}
	else densr = 0.0;

	//Detecting the mode for clay representation:
	withclaymesh = clayYiPID >= 0;
	withtimedepKclays = clayeqkvPID >= 0;

	//BETA-TESTING unsafe feature
	if (domasstransport) {
		if (BETApermanentDirichletMassConcBCs_LateAtClayTop && !withclaymesh) {
			IfmWarning(pDoc, "[CHG TO NON-BETA MSG] BETApermanentDirichletMassConcBCs_LateAtClayTop == true has no effect in modes other than withclaymesh.");
		}
		if (BETApermanentDirichletMassConcBCs_LateAtClayTop) {
			IfmWarning(pDoc, "[CHG TO NON-BETA MSG] BETApermanentDirichletMassConcBCs_LateAtClayTop == true hard-wired in the plugin. Be careful.");
			//IfmAlert(pDoc, NULL, "  OK  ", "Beware: BETApermanentDirichletMassConcBCs_LateAtClayTop is active in the plugin DLL!");
		}
	}

	if (!withclaymesh && domasstransport) {
		if (bclaycPID < 0) {
			withbotclayconc = false;
			sprintf_s(txtbuffer, 180, "[PreSimul] Detected no time series named 'botclayconc'. Values(t) from 'swconc' will be used directly for c(top rock), thus producing questionable results!");
			IfmWarning(pDoc, txtbuffer);
		}
		else {
			//Parameter values below are fixed in the code and should correspond to the values used for generating the TS in my R script
			origclayL = 20.0;
			endofmaxsalinitytime = 1500.0 * daysperyear;
			maxtransftimebclaycTS = 50000.0 * daysperyear;
			withbotclayconc = true;
			IfmInfo(pDoc, "[PreSimul] A 'botclayconc' TS was detected. Slow 1Dz diffusive mass transport within the clay aquitard will thus be considered when computing c(t) BCs at top nodes.");
			sprintf_s(txtbuffer, 180, "** However, be careful as related params are hard-wired in the plugin: origclayL=%.1f m, endofmaxsalinitytime=%.0f a, maxtransftimebclaycTS=%.0f a.", origclayL, endofmaxsalinitytime / daysperyear, maxtransftimebclaycTS / daysperyear);
			IfmWarning(pDoc, txtbuffer);
		}
		if (withtimedepKclays) {
			IfmWarning(pDoc, "[PreSimul] A 'clayeqkv' TS was detected, but it will be ignored as it is compatible only with withclaymesh mode.");
			withtimedepKclays = false;
		}
	}
	else withbotclayconc = false;

	//Management of top* nodal | elemental markers (User Data)
	long rdidtn = IfmGetNodalRefDistrIdByName(pDoc, "topnodes"); //ref. distr. id for 'topnodes'
	long elzmRDId = IfmGetElementalRefDistrIdByName(pDoc, "topelemszmid");
	long rdidlyiE = IfmGetElementalRefDistrIdByName(pDoc, "claylayeri");
	if (!withclaymesh) {
		if (rdidtn < 0) {
			IfmAlert(pDoc, NULL, "  OK  ", "The plugin requires a User Data nodal distrib. named topnodes.");
			IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
			return;
		}
		wiseTRupdate = elzmRDId >= 0;
		if (!wiseTRupdate) {
			IfmWarning(pDoc, "[PreSimul] Detected no User Data elemental distrib. named 'topelemszmid'. Old & flawed method will be used for linking nodal z to elements.");
		}
		if (rdidlyiE >= 0) IfmInfo(pDoc, "[PreSimul] The 'claylayeri' elem. user-data distrib. will not be used in this simulation, due to the NON-withclaymesh mode.");
	}
	else {
		/* if withclaymesh... */
		wiseTRupdate = false;
		if (rdidtn >= 0) IfmInfo(pDoc, "[PreSimul] The 'topnodes' nodal user-data distrib. will not be used in this simulation, due to the withclaymesh mode.");
		if (elzmRDId >= 0) IfmInfo(pDoc, "[PreSimul] The 'topelemszmid' nodal user-data distrib. will not be used in this simulation, due to the withclaymesh mode.");
		if (rdidlyiE < 0) {
			IfmAlert(pDoc, NULL, "  OK  ", "The withclaymesh mode requires an elemental user-data distrib. named 'claylayeri'.");
			IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
			return;
		}
	}

	//NOW COMMON TO ALL MODES (to allow export @ top nodes)
	long rdgxraw = -1; //index for the 'globXiraw' nodal ref. distrib. (@User Data); for ALL MODES!
	rdgxraw = IfmGetNodalRefDistrIdByName(pDoc, "globxiraw");
	isglobXipresent = rdgxraw >= 0;
	doexportdata = isglobXipresent;

	//BUT overwritten by the alternative export mode, if applicable:
	long rdusrexpnod = -1; //index for the 'globXiraw' nodal ref. distrib. (@User Data); for ALL MODES!
	rdusrexpnod = IfmGetNodalRefDistrIdByName(pDoc, "exportnodes");
	altexportusernodes = rdusrexpnod >= 0;
	altexportautoinclgtop = false; //default early initialization
	doexportdata = doexportdata || altexportusernodes;
	doexportdata = true; //TEMP HARD-WIRE to make sure the altexportautoinclgtop can apply even when "exportnodes" ref. distrib. is absent. TODO SOON: Restructure the code below to better determine the export=true vs export=false scenarios.

	//Management of clayyi*|clayxi* nodal | elemental markers (User Data)
	long rdcymain = -1; //index for the 'clayyimain' nodal ref. distrib. (@User Data)
	long rdcysub = -1; //index for the 'clayyisub' nodal ref. distrib. (@User Data)
	long rdcxraw = -1; //index for the 'clayxiraw' nodal ref. distrib. (@User Data)
	if (!withclaymesh) {
		IfmInfo(pDoc, "[PreSimul] Detected no time series named 'clayyi'. The mode 'withclaymesh' won't apply there.");
	}
	else {
		IfmInfo(pDoc, "[PreSimul.WCM] A 'clayyi' time series is detected. Enabling the withclaymesh (WCM) mode, which performs explicit update of active clay elements & nodes with time.");
		//IfmWarning(pDoc, "* BETA: Beware, however. The 'withclaymesh' mode still is in BETA stage of development!");
		rdcymain = IfmGetNodalRefDistrIdByName(pDoc, "clayyimain");
		rdcysub = IfmGetNodalRefDistrIdByName(pDoc, "clayyisub");
		rdcxraw = IfmGetNodalRefDistrIdByName(pDoc, "clayxiraw");
		if ((rdcymain < 0) || (rdcysub < 0) || (rdcxraw < 0)) {
			withclaymesh = false;
			IfmWarning(pDoc, "[PreSimul.WCM] Sorry: forced to DISABLE withclaymesh mode, since some of clayyimain | clayyisub | clayxiraw nodal user-data distributions are MISSING!");
		}
		else {
			//Initialization of the variable with time series value at t0:
			currclayYi_rawreal = IfmInterpolatePowerValue(pDoc, clayYiPID, simulstarttime);
			currclayYi = floor(currclayYi_rawreal + int_tol);
		}
	}

	//MARK: At this point, we now know for sure if withclaymesh is active or not.

	clayYi_FixMassBC_thresv = 9999999.9; //default, quasi-infinite value, so that by default, clayYi never is greater than this threshold value
	if (clayYiFixMBCPID >= 0) {
		if (withclaymesh && BETApermanentDirichletMassConcBCs_LateAtClayTop) {
			clayYi_FixMassBC_thresv = IfmInterpolatePowerValue(pDoc, clayYiFixMBCPID, 0.0);
			sprintf_s(txtbuffer, 180, "[PreSimul.WCM] Fixed Mass-Conc. BC's will apply to clay-top nodes as soon as clayYi >= %.3f (clayYi_FixMassBC_thresv).", clayYi_FixMassBC_thresv);
			IfmInfo(pDoc, txtbuffer);
		}
		else {
			IfmInfo(pDoc, "[PreSimul] Detected a time series named 'clayyifixmbc' but not used since it is irrelevant when withclaymesh == false.");
		}
		IfmInfo(pDoc, "** (related with hard-wired internal option BETApermanentDirichletMassConcBCs_LateAtClayTop!)");
	}
	else {
		if(withclaymesh) IfmInfo(pDoc, "[PreSimul.WCM] Detected no time series named 'clayyifixmbc'. Mass-transport BCs at clay-top nodes will continue to depend on Vy direction.");
	}

	//NEW: globxiE...
	isaglobxiEloaded = false; //safe initial default
	PglobXiElemA = nullptr; //safe initial default
	maxglobXiEval = -9999; //initialization
	minglobXiEval = +9999; //initialization
	long rdidglobxiE = IfmGetElementalRefDistrIdByName(pDoc, "globxiE"); //ref. distrib. Id for the globxiE elemental user data
	if (rdidglobxiE >= 0) {
		PglobXiElemA = new int[nelems];
		double rawgxiev; //raw globxiE value being read
		int cntbadgxie = 0;
		for (i = 0; i < nelems; i++) {
			rawgxiev = IfmGetElementalRefDistrValue(pDoc, rdidglobxiE, i);
			PglobXiElemA[i] = floor(rawgxiev + 0.5); //for numeric rounding
			if (PglobXiElemA[i] >= 0) {
				if (PglobXiElemA[i] > maxglobXiEval) maxglobXiEval = PglobXiElemA[i];
				if (PglobXiElemA[i] < minglobXiEval) minglobXiEval = PglobXiElemA[i];
			}
			else {
				PglobXiElemA[i] = -1; //to make sure all non 'valid' indeces have the the same diagnostic index value
				cntbadgxie++;
			}
		}
		if (cntbadgxie == 0) {
			isaglobxiEloaded = true;
			sprintf_s(txtbuffer, 180, "[PreSimul BETA] isaglobxiEloaded == true : could load 'globxiE' successfully (maxglobXiEval = %d; min = %d).", maxglobXiEval, minglobXiEval);
			IfmInfo(pDoc, txtbuffer);
			IfmWarning(pDoc, "[PreSimul BETA] HOWEVER, globxiE values are not used yet: for the moment, the source term coefficients are provided manually. [!!!!!]");
		}
		else {
			isaglobxiEloaded = false;
			sprintf_s(txtbuffer, 180, "[PreSimul BETA] isaglobxiEloaded == false : could NOT load 'globxiE' successfully, as there are %d elements with invalid index values!", cntbadgxie);
			IfmWarning(pDoc, txtbuffer);
			delete[] PglobXiElemA;
			PglobXiElemA = nullptr;
		}
	}

	//BETA for NEW source term simul.
	long rdclaysrcE = IfmGetElementalRefDistrIdByName(pDoc, "betaclaysrcE");
	BETAwithsourceterms = withclaymesh && clayqdbotPID >= 0 && rdclaysrcE >= 0;
	if (BETAwithsourceterms) IfmWarning(pDoc, "* BETA: [PreSimul.WCM] BETAwithsourceterms == true; NEW EARLY-BETA FEATURE!!!");
	if (BETAwithsourceterms) {
		//Specific initializations... (BETA)
		areSrcTermsAlreadyZeroed = false;
	}


	//Management of implicit clay accum. modes (repr. as transient Cauchy BC transfer rates...)
	if (!withclaymesh) {
		withclayConstsedrate = claytopPID >= 0; //has priority over the mode right below
		withclayTimeDepsedrate = sedratePID >= 0; //(lower priority mode)
		if (!withclayConstsedrate) {
			IfmInfo(pDoc, "[PreSimul] Detected no time series named 'claytop'. The 'withclayConstsedrate' cannot apply there.");
		}
		if (!withclayTimeDepsedrate) {
			IfmInfo(pDoc, "[PreSimul] Detected no time series named 'sedrate'. The 'withclayTimeDepsedrate' cannot apply there.");
		}
		if (withclayTimeDepsedrate && withclayConstsedrate) {
			IfmWarning(pDoc, "[PreSimul] Two concurrent modes detected related to clay accum. (both 'claytop' & 'sedrate' TS present). Enabling withclayConstsedrate mode (thereby ignoring the 'sedrate' TS).");
			withclayTimeDepsedrate = false;
		}
		if (!withclayTimeDepsedrate && !withclayConstsedrate) {
			IfmWarning(pDoc, "[PreSimul] Clay accum. cannot be considered in the simulation, since none of the related TS were found!");
		}
		if(withclayTimeDepsedrate) {
			IfmInfo(pDoc, "[PreSimul] Enabling the 'withclayTimeDepsedrate' mode. Cauchy b.c. transfer rates will be computed according to cumulative clay accum. at a global time-dep. rate (based on 'sedrate' TS).");
		}
		if (withclayConstsedrate) {
			IfmInfo(pDoc, "[PreSimul] Enabling the 'withclayConstsedrate' mode. Cauchy b.c. transfer rates will be computed according to regular clay accum. at locally constant sedim. rate (since 'claytop' TS detected)");
			claytoptargetYval = IfmInterpolatePowerValue(pDoc, claytopPID, 0.0);
			if (seddurPID < 0) {
				IfmAlert(pDoc, NULL, "  OK  ", "The plugin's withclayConstsedrate mode requires a time series named 'sedimduring', but none is detected. Simulation aborted.");
				IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
				withclayConstsedrate = false;
				return;
			}
			else {
				claysedimduring = IfmInterpolatePowerValue(pDoc, seddurPID, 0.0);
				sprintf_s(txtbuffer, 180, "[PreSimul] Clays will reach a top elevation of %.3f after sedimenting from 0.0 to %.3f years. (withclayConstsedrate mode)", claytoptargetYval, claysedimduring / daysperyear);
				IfmInfo(pDoc, txtbuffer);
				//withclayConstsedrate = true; ALREADY SET
			}
		}
	}
	else {
		//if withclaymesh...
		withclayTimeDepsedrate = false;
		withclayConstsedrate = false;
		if (withtimedepKclays) {
			IfmWarning(pDoc, "* BETA: [PreSimul.WCM] Enabling the withtimedepKclays sub-mode. Clay elements will have time-dependent K values as soon and as long as they are activated.");
			IfmWarning(pDoc, "** HOWEVER, this new feature uses a single time series (typically made for the max. final clay accum.).");
			IfmWarning(pDoc, "** K values will thus likely be UNDERestimated where local final thickness is lower.");
		}
	}
	withSOMEclayaccumANYmode = withclaymesh || withclayConstsedrate || withclayTimeDepsedrate;

	//Initial conditions and Model initialization
	currsedimrate = 0.0; //initializing; will be updated in PreTimeStep callbacks
	firstTSresumed = simulstarttime > time_tol;
	isFirstExportedTS = true; //TODO: Try to find a way to deal with simulation resume, and netCDF export!...
	bool loadedclayaccum = false; //internal: if true, will import clayacc User Data as initial clay thicknesses (nodal or elemental) @ simulation start; relevant ONLY IF !withclayConstsedrate
	bool loadedfloodeddur = false; //internal: if true, will import floodeddur User Data as initial flooded duration @ simulation start
	long rdidcaN = -1; //index for the 'clayaccumnodal' nodal ref. distrib. (@User Data) (in PreSimulation)
	long rdidcaE = -1; //index for the 'clayaccumelem' elemental ref. distrib. (@User Data) (in PreSimulation)
	long rdidfd = -1; //index for the 'flooddur' nodal ref. distrib. (@User Data) (in PreSimulation)

	if (firstTSresumed) {
		IfmWarning(pDoc, "[RISKY OPERATION @ PreSimul] Continuing (~resuming) simulation starting with current state and BC's (unchanged for the first iteration).");

		if (wiseTRupdate) {
			//elemental user data presence verification: 'clayaccumelem' (cumulative clay accumulated ~over this element)
			if (!withclayConstsedrate) { //since such info ain't required in the withclayConstsedrate mode
				rdidcaE = IfmGetElementalRefDistrIdByName(pDoc, "clayaccumelem");
				if (rdidcaE < 0) IfmWarning(pDoc, "[PreSimul] 'clayaccumelem' RefDistr (User Data) could not be loaded.");
				if (rdidcaE >= 0) {
					loadedclayaccum = true;
				}
			}
		}
		else {
			//nodal user data presence verification: 'clayaccumnodal' (cumulative clay accumulated ~over this node)
			rdidcaN = IfmGetNodalRefDistrIdByName(pDoc, "clayaccumnodal");
			if (rdidcaN < 0) IfmWarning(pDoc, "[PreSimul] 'clayaccumnodal' RefDistr (User Data) could not be loaded.");
			if (rdidcaN >= 0) {
				loadedclayaccum = true;
			}
		}
		
		//nodal user data presence verification: 'flooddur' (flood duration)
		rdidfd = IfmGetNodalRefDistrIdByName(pDoc, "flooddur");
		if (rdidfd < 0) IfmWarning(pDoc, "[PreSimul] 'flooddur' RefDistr (User Data) could not be loaded.");
		if (rdidfd >= 0) {
			loadedfloodeddur = true;
		}
	}
	else {
		loadedclayaccum = false;
		loadedfloodeddur = false;
	}
	//std::string InfoStr = "rslPID=" + std::to_string(rslPID);
	//IfmInfo(pDoc, InfoStr.c_str());
	//Vtopnodes.resize(nnodes);
	
	//DEBUG info: IfmInfo(pDoc, "allocating of the internal arrays...");
	
	/* Allocation of arrays: 1. COMMON TO ALL MODES */
	Ptopnodes = new bool[nnodes];
	Ptopelems = new bool[nelems];
	Pcurrfloodednodes = new bool[nnodes];
	Pinflownodes = new bool[nnodes];
	PallnodesY = new double[nnodes];

	/* Allocation of arrays: 2. specific to OLDER MODES */
	if (!withclaymesh) {
		Ptopelemszmid = new double[nelems];
		Pfloodedduring = new double[nnodes];
		Pclayacc_nodal = new double[nnodes];
		Pclayacc_elemtal = new double[nelems];
		Pclayaccfinal_nodal = new double[nnodes];
		Pclaysedrate_nodal = new double[nnodes];
		Pclaysedrate_elemtal = new double[nelems];
		PelemsAtnodesBR = new long[nnodes];
		PelemsAtnodesBL = new long[nnodes];
		PelemsAtnodesBM = new long[nnodes];
	}

	/* Allocation of arrays: 3. specific to THE NEW withclaymesh MODE */
	if (withclaymesh) {
		PclayYimain = new int[nnodes];
		PclayYisub = new int[nnodes];
		PclayXiraw = new int[nnodes];
		Pclaylayeri = new int[nelems];
		PtopclayBCnodes = new bool[nnodes];
		PtoprockBCnodes = new bool[nnodes];
		PupdateNbelowclaytop = new bool[nnodes];
		PupdateNoverclaytop = new bool[nnodes];
		PcleanNBCprevclaytop = new bool[nnodes];

		//Arrays to be completely filled at relevant time steps (for ~copy-paste to bypass clayYi-changing the t.step)
		Pallnodes_heads = new double[nnodes];
		Pallnodes_mconcs = new double[nnodes];
		//Other arrays...
		Pallcurractivenodes = new bool[nnodes];
		Pallcurrclayactivenodes = new bool[nnodes];
		Psteadyallrocknodes = new bool[nnodes];
		Psteadytoprocknodes = new bool[nnodes];
		//...for the source term feature
		if (BETAwithsourceterms) PclaysrcBetaCoefE = new double[nelems]; else PclaysrcBetaCoefE = nullptr;
		//...for the active state of elements
		PelemcurractiveE = new bool[nelems];
	}

	IfmInfo(pDoc, withSOMEclayaccumANYmode ? "[PreSimul] The model includes a clay accumulation process. Therefore: evolving upper boundary (top BC's and/or elements & nodes)." : "[PreSimul] The model does NOT include any clay accumulation process. Therefore: simple Dirichlet BCs for flow.");

	//NO MORE RELEVANT:
	/* if (!firstTSresumed) {
		IfmInfo(pDoc, domasstransport ? "[PreSimul] Initializing flow & transport related data, as well as plugin internal arrays..." : "[PreSimul] Initializing flow related data, as well as plugin internal arrays...");
	} */

	/* Importing Y coordinate data */
	for (i = 0; i < nnodes; i++) {
		PallnodesY[i] = IfmGetY(pDoc, i);
	}

	//Warns the user if BC constraints are present in the model (as it should be avoided!)
	if(domasstransport) CheckPresenceOfBCCdata(pDoc, false);

	//Preparing arrays for Data Export (optional, only if 'globxiraw' user data is present)
	//(otherwise, arrays are == nullptr, maxglobXival remains equal to -9999, and nbglobXivals equals 0)
	maxglobXival = -9999;
	minglobXival = +9999;
	nbglobXivals = 0;
	if (isglobXipresent) {
		double gxirawv;
		PglobXiraw = new int[nnodes]; //Array storing global Xi indeces...
		for (i = 0; i < nnodes; i++) {
			//global Xi indeces:
			gxirawv = IfmGetNodalRefDistrValue(pDoc, rdgxraw, i); //extracting from User Data
			PglobXiraw[i] = floor(gxirawv + 0.5); //storing into internal array
			if (PglobXiraw[i] > maxglobXival) maxglobXival = PglobXiraw[i]; //maximum value progressive detection
			if (PglobXiraw[i] < minglobXival) minglobXival = PglobXiraw[i]; //minimum value progressive detection
		}
		nbglobXivals = maxglobXival - minglobXival + 1;
		sprintf_s(txtbuffer, 180, "[Export init.] Optional Export/Aligned-Storage of simulation results is ACTIVE, with maxglobXival = %d (min = %d), hence nbglobXivals = %d.", maxglobXival, minglobXival, nbglobXivals);
		IfmInfo(pDoc, txtbuffer);

		//Then, allocating arrays of size dependent on nbglobXivals (and thus on maxglobXival-minglobXival) defined just above
		Pcurrglobtopnodes = new int[nbglobXivals];
		Pcurrglobtop_isinflow = new bool[nbglobXivals]; //(current...)
		Pprevglobtop_isinflow = new bool[nbglobXivals]; //(PREVIOUS...)
		Pcurrglobtop_nodalmassrbudgets = new double[nbglobXivals];
	}
	else {
		PglobXiraw = nullptr;
		Pcurrglobtopnodes = nullptr;
		Pcurrglobtop_isinflow = nullptr;
		Pprevglobtop_isinflow = nullptr;
		Pcurrglobtop_nodalmassrbudgets = nullptr;
		//(...and nbglobXivals will remain equal to 0)
	}

	Puexportnodes = nullptr; //initial default
	nuserexportnodes = 0;
	neffexportnodes = 0;
	if (doexportdata) {
		if (domasstransport && EXPORT_MASS_RBUDGETs) {
			Pnodalmassrbudgets = new double[nnodes]; //Array storing mass budget computed values...
			/* Common for both export modes */
			for (i = 0; i < nnodes; i++) {
				//initializing computed-mass-rate-budget internal array:
				Pnodalmassrbudgets[i] = nodata_double;
			}
		}
		else {
			Pnodalmassrbudgets = nullptr;
		}
		/* Specific to the alternative export mode */
		if(altexportusernodes) {
			double expnodv;
			int expnodintv;
			bool expnodactiv;
			Puexportnodes = new int[nnodes]; //Array storing global Xi indeces...
			nuserexportnodes = 0;
			for (i = 0; i < nnodes; i++) {
				//user export nodes indeces:
				expnodv = IfmGetNodalRefDistrValue(pDoc, rdusrexpnod, i); //extracting from User Data
				expnodintv = floor(expnodv + 0.5); //storing into internal array
				expnodactiv = expnodintv >= 0;
				Puexportnodes[i] = expnodactiv ? expnodintv : -1; //This will deal with user nodata's (-99999 or alike)
				if (expnodactiv) nuserexportnodes++;
			} //hence at the end of this for loop, nuserexportnodes is defined.
			/* Then, filling the filtered small array */
			PuserexportnodesA = new int[nuserexportnodes];
			Puserexportnodes_uvaluesA = new int[nuserexportnodes];
			int j_filling = 0;
			for (i = 0; i < nnodes; i++) {
				expnodintv = Puexportnodes[i];
				if (expnodintv < 0) continue;
				PuserexportnodesA[j_filling] = i;
				Puserexportnodes_uvaluesA[j_filling] = expnodintv;
				j_filling++;
			}
			if (j_filling != nuserexportnodes) {
				IfmWarning(pDoc, "[Export init.] Strange inequality: j_filling != nuserexportnodes. [BETA!]");
			}

			altexportautoinclgtop = withclaymesh && (nbglobXivals > 0);
			if (altexportautoinclgtop) {
				errno_t err;
				neffexportnodes = nbglobXivals + nuserexportnodes;
				PeffexportnodesA = new int[neffexportnodes];
				Peffexportnodes_uvaluesA = new int[neffexportnodes];
				err = memcpy_s(&PeffexportnodesA[nbglobXivals], neffexportnodes * sizeof(int), PuserexportnodesA, nuserexportnodes * sizeof(int));
				if (err) IfmWarning(pDoc, "Error executing memcpy_s for initializing PeffexportnodesA."); //BETA; such error should never happen.
				err = memcpy_s(&Peffexportnodes_uvaluesA[nbglobXivals], neffexportnodes * sizeof(int), Puserexportnodes_uvaluesA, nuserexportnodes * sizeof(int));
				if (err) IfmWarning(pDoc, "Error executing memcpy_s for initializing Peffexportnodes_uvaluesA."); //BETA; such error should never happen.
				for (i = 0; i < nbglobXivals; i++) {
					PeffexportnodesA[i] = -1; //temporary invalid value that will be replaced as soon as curr.glob.top.nodes are defined-initialized
					Peffexportnodes_uvaluesA[i] = expNautotopId; //initializing the Id value for all curr.glob.top.nodes
				}
			}
			else {
				neffexportnodes = nuserexportnodes;
				PeffexportnodesA = PuserexportnodesA;
				Peffexportnodes_uvaluesA = Puserexportnodes_uvaluesA;
			}

		} //BOF(OLD-COMMENT): else it's already done in the if(isglobXipresent) above...

		/* Optional auto-adjustment of the export-delta-time */
		if (autodefine_exportdeltatime) {
			double simdur_apriori_yrs = (simulendtime_apriori - simulstarttime) / daysperyear;
			double tmp_deltat_yrs;
			tmp_deltat_yrs = simdur_apriori_yrs > 1000.0 ? 1.0 :
				simdur_apriori_yrs > 100.0 ? 0.1 :
				0.05;
			exportdeltatime = tmp_deltat_yrs * daysperyear;
			sprintf_s(txtbuffer, 180, "[Export init.] Frequence of export was adjusted automatically to exportdeltatime ~ %.3f yrs.", tmp_deltat_yrs);
			IfmInfo(pDoc, txtbuffer);
		}
	}
	else {
		Pnodalmassrbudgets = nullptr;
	}

	double yimainv, yisubv, xirawv, layeryiv, srciv; //NEW BETA for withclaymesh; all are clay-related indeces
	maxclaylayeri = -9999;
	minclaylayeri = +9999;
	maxclayXival = -9999;
	minclayXival = +9999;
	nbclayXivals = 0; //just for further safety
	maxclayYival = -9999;
	minclayYival = +9999;
	if (withclaymesh) {
		//BEGIN of the withclaymesh-mode section:

		/* WCM PART 1: Nodal indeces */
		for (i = 0; i < nnodes; i++) {
			yimainv = IfmGetNodalRefDistrValue(pDoc, rdcymain, i);
			yisubv = IfmGetNodalRefDistrValue(pDoc, rdcysub, i);
			xirawv = IfmGetNodalRefDistrValue(pDoc, rdcxraw, i);

			//arrays storing input data
			PclayYimain[i] = floor(yimainv + 0.5); //floor(x+0.5) does the same thing as round(x) would, and works well with negative numbers as well
			PclayYisub[i] = floor(yisubv + 0.5); //idem round technique as above
			PclayXiraw[i] = floor(xirawv + 0.5); //idem round technique as above

			//maximum value progressive detection
			if (PclayXiraw[i] >= 0) {
				if (PclayXiraw[i] > maxclayXival) maxclayXival = PclayXiraw[i];
				if (PclayXiraw[i] < minclayXival) minclayXival = PclayXiraw[i];
			}
			if (PclayYimain[i] >= 0) {
				if (PclayYimain[i] > maxclayYival) maxclayYival = PclayYimain[i];
				if (PclayYimain[i] < minclayYival) minclayYival = PclayYimain[i];
			}

			Psteadyallrocknodes[i] = (PclayYimain[i] < 0) && (PclayYisub[i] < 0); //TODO-maybe: ...Or could be prepared using nodefacies (in which case this User Data would have to be loaded prior to this...)
			//clayYimain values should be either -1 (most often) or -2 (for the unconfined top rock nodes) in the rock subdomain;
			//clayYisub values should be -1 for every node in the rock subdomain.
			//The dual condition ensures no clay node is included.
			//Therefore, this boolean selection excludes the top-rock nodes which are also clay-base nodes.

			Psteadytoprocknodes[i] = (PclayYimain[i] == 0) || (PclayYimain[i] == rockYi_cstval);
		}
		nbclayXivals = maxclayXival - minclayXival + 1;
		sprintf_s(txtbuffer, 180, "[PreSimul.WCM] Dimensions detected for the clay 'numerical' facies (NODAL): maxclayXival = %d (min = %d), hence nbclayXivals = %d; maxclayYival = %d (min = %d)", maxclayXival, minclayXival, nbclayXivals, maxclayYival, minclayYival);
		IfmInfo(pDoc, txtbuffer);
		if (maxclayXival < 0 || minclayXival < 0 || minclayXival == +9999 || nbclayXivals <= 0) {
			IfmAlert(pDoc, NULL, "  OK  ", "The plugin detected an inconsistent range for the clayxiraw values.\n(This range should respect: 0 <= minimum <= maximum.)\n(...so that nbclayXivals >= 1)");
			IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
			return;
		}
		if (maxclayYival < 0 || minclayYival != 0) {
			IfmAlert(pDoc, NULL, "  OK  ", "The plugin detected an inconsistent range for the clayyimain values.\n(This range should respect: 0 == minimum <= maximum.)");
			IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
			return;
		}

		/* WCM PART 2: Elemental indeces */
		for (i = 0; i < nelems; i++) {
			layeryiv = IfmGetElementalRefDistrValue(pDoc, rdidlyiE, i);
			Pclaylayeri[i] = floor(layeryiv + 0.5); //idem as in node loop above (i.e. for numeric rounding)
			if (Pclaylayeri[i] >= 0) {
				//NOTE: This way, it is possible to detect non-optimal models in which minclaylayeri == 0 rather than 1.
				if (Pclaylayeri[i] > maxclaylayeri) maxclaylayeri = Pclaylayeri[i];
				if (Pclaylayeri[i] < minclaylayeri) minclaylayeri = Pclaylayeri[i];
			}
			//...and for the source terms...
			if (BETAwithsourceterms) {
				srciv = IfmGetElementalRefDistrValue(pDoc, rdclaysrcE, i);
				PclaysrcBetaCoefE[i] = !isnan(srciv) && srciv > 0.0 ? srciv : -1.0;
			}
			PelemcurractiveE[i] = IfmGetMatElementActive(pDoc, i) != 0; //That's just the PreSimul initialization of the array; values normally are overwritten right after first @PreTimeStep, in WCM.UpdateClayParamsAccToYi(...)
		}
		sprintf_s(txtbuffer, 180, "[PreSimul.WCM] Dimensions detected for the clay 'numerical' facies (ELEMENTAL): maxclaylayeri = %d (min = %d)", maxclaylayeri, minclaylayeri);
		IfmInfo(pDoc, txtbuffer);
		if (maxclaylayeri < 1 || minclaylayeri < 1) {
			IfmAlert(pDoc, NULL, "  OK  ", "The plugin detected an inconsistent range for the claylayeri values.\n(This range should respect: 1 == minimum <= maximum.)");
			IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
			return;
		}

		//Array allocation operations which required nbclayXivals and thus maxclayXival to be defined:
		Pcurrclaytopnodes = new int[nbclayXivals];
		Pprevclaytopnodes = new int[nbclayXivals];
		Pcurrclaytopisinflow = new bool[nbclayXivals];
		Pclaytopnodes_tmp_prevtophbc = new double[nbclayXivals];
		Pclaytopnodes_tmp_currtopnewhbc = new double[nbclayXivals];
		Pclaytopnodes_tmp_prevtopY = new double[nbclayXivals];
		Pclaytopnodes_tmp_currtopY = new double[nbclayXivals];
		Pclaytopnodes_tmp_prevtopcbc = new double[nbclayXivals];
		Pclaytopnodes_tmp_currtopnewcbc = new double[nbclayXivals];
		Pclaytopnodes_tmp_currtopcbc_chgneeded = new bool[nbclayXivals];

		prevclayYi = currclayYi; //so that this phase of initialization of claymesh is not interpreted as a "change" in clayYi...

		//END of the withclaymesh-mode section.
	}
	else {
		//BEGIN of the NON-withclaymesh-mode section:
		int gxival; //temporary variable for global Xi indeces...
		int gxiindex;
		ntopnodes = 0;
		for (i = 0; i < nnodes; i++) {
			nodval = IfmGetNodalRefDistrValue(pDoc, rdidtn, i);
			Ptopnodes[i] = (nodval > len_tol);
			Pfloodedduring[i] = loadedfloodeddur ? IfmGetNodalRefDistrValue(pDoc, rdidfd, i) : 0.0;
			if (Ptopnodes[i]) {
				ntopnodes++;
				Pcurrfloodednodes[i] = false; //initialization
				Pclayaccfinal_nodal[i] = (claytoptargetYval - PallnodesY[i]); //NEW, BETA!
				if (Pclayaccfinal_nodal[i] < 0.0) Pclayaccfinal_nodal[i] = 0.0; //prevents negative clay thicknesses where ground.surface.Y > target.clay.top.Y
				Pclaysedrate_nodal[i] = (claytoptargetYval - PallnodesY[i]) / claysedimduring;
				if (Pclaysedrate_nodal[i] < 0.0) Pclaysedrate_nodal[i] = 0.0; //prevents negative sed. rates where ground.surface.Y > target.clay.top.Y
				Pinflownodes[i] = firstTSresumed ? IfmGetBcMassType(pDoc, i) == IfmBC_DIRICHLET : true;
				/*      if(!firstTSresumed) { ##### This section is not relevant anymore, I think! (Nov 16th, 2015)
				//ERRONEOUS! IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_DIRICHLET, 0, PallnodesY[i]); //fw heads to be updated at first PreTimeStep
				IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_CAUCHY, 0, hbcval);
				IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_DIRICHLET, 0, fwconc);
				//IfmSetBccMassTypeAndValueAtCurrentTime(pDoc, i, 1, 0, 0, 0);
				} */

				if (!wiseTRupdate && withSOMEclayaccumANYmode) {
					Pclayacc_nodal[i] = loadedclayaccum ? IfmGetNodalRefDistrValue(pDoc, rdidcaN, i) : 0.0;
					PelemsAtnodesBR[i] = IfmFindElementAtXY(pDoc, IfmGetX(pDoc, i) + 0.001, PallnodesY[i] - 0.001, 0.0);
					PelemsAtnodesBL[i] = IfmFindElementAtXY(pDoc, IfmGetX(pDoc, i) - 0.001, PallnodesY[i] - 0.001, 0.0);
					PelemsAtnodesBM[i] = IfmFindElementAtXY(pDoc, IfmGetX(pDoc, i), PallnodesY[i] - 0.001, 0.0);
				}
				else {
					//Following arrays are not of any use when there is NO clay accum.
					Pclayacc_nodal[i] = -1.0;
					PelemsAtnodesBR[i] = -1;
					PelemsAtnodesBL[i] = -1;
					PelemsAtnodesBM[i] = -1;
				}

				//Defining Pcurrglobtopnodes once for all (as it's steady in NON-withclaymesh modes)
				if (isglobXipresent) {
					gxival = PglobXiraw[i];
					if (gxival >= 0) {
						gxiindex = gxival - minglobXival;
						Pcurrglobtopnodes[gxiindex] = i;
						Pcurrglobtop_isinflow[gxiindex] = Pinflownodes[i]; //NEW line BETA: Is it ok? >>TODO test in non-WCM mode.
						//PROGRAMMER'S NOTE: It is irrelevant here to update PeffexportnodesA[gxiindex] since this array does never include curr.glob.top.nodes in non-WCM modes.
					}
					else {
						sprintf_s(txtbuffer, 180, "[PreSimul] Unexpected globXiraw (%d < 0) at a Ptopnodes[i=%d] == true node. Why?", gxival, i);
						IfmWarning(pDoc, txtbuffer);
					}
				}
			}
			else {
				//(else: not a top node...)
				if (!wiseTRupdate) {
					Pclayacc_nodal[i] = -1.0;
					PelemsAtnodesBR[i] = -1;
					PelemsAtnodesBL[i] = -1;
					PelemsAtnodesBM[i] = -1;
				}
			}
		}
		/* //DEBUG
		sprintf_s(txtbuffer, 180, "rdidtn=%d, rdidcaN=%d, rdidfd=%d, elzmRDId=%d", rdidtn, rdidcaN, rdidfd, elzmRDId);
		IfmWarning(pDoc, txtbuffer);
		//DEBUG */
		if (wiseTRupdate) {
			IfmInfo(pDoc, "[PreSimul] Clay accumulation in the 'wiseTRupdate' elemental mode: preparing internal arrays...");
			for (i = 0; i < nelems; i++) {
				elemval = IfmGetElementalRefDistrValue(pDoc, elzmRDId, i);
				Ptopelemszmid[i] = elemval;
				Ptopelems[i] = (elemval > len_tol); //BETA (verifying against NoData would be better)
				if(Ptopelems[i]) ntopelems++;
				Pclayacc_elemtal[i] = loadedclayaccum ? IfmGetElementalRefDistrValue(pDoc, rdidcaE, i) : 0.0;
				Pclaysedrate_elemtal[i] = (claytoptargetYval - elemval) / claysedimduring;
				if (Pclaysedrate_elemtal[i] < 0.0) Pclaysedrate_elemtal[i] = 0.0; //prevents negative sed. rates where ground.surface.Y > target.clay.top.Y
			}
		}
		needupdatearoundclaytop = false; //actually never required in the older modes
		hasclayYichanged = false; //actually never required in the older modes

		//Just to make it clearer...
		Pcurrclaytopnodes = nullptr;
		Pprevclaytopnodes = nullptr;
		Pcurrclaytopisinflow = nullptr;

		//END of the NON-withclaymesh-mode section.
	}

	/* Initializing timers */
	time(&realtime_atsimstart);
	realtime_nextrefresh = realtime_atsimstart + realtime_ProgRefrDelta;
	rt_prev_remain_pred = realtime_atsimstart + daysperyear * secsperday; //some far away value : = 1 year

	initialpreTSupdate = true; //so that the first Update... function calls (for clay aquitard and/or for in/outflow data) will initialize instead of just updating
	IfmInfo(pDoc, "<< preSimulation completed >>");
}

void Cpaleosea2d::PostSimulation (IfmDocument pDoc)
{
	int i;
	if (!withclaymesh) {
		if (wiseTRupdate) {
			long rdidcaE = IfmGetElementalRefDistrIdByName(pDoc, "clayaccumelem"); //ref. distrib. index for the 'clayaccumelem' elemental ref. distrib. (@User Data) (in PostSimulation)
			if (rdidcaE<0) rdidcaE = IfmCreateElementalRefDistr(pDoc, "clayaccumelem");
			if (rdidcaE<0) {
				IfmWarning(pDoc, "[PostSimul] 'clayaccumelem' RefDistr (User Data) could not be created for some reason.");
			}
			else {
				IfmInfo(pDoc, "[PostSimul] 'clayaccumelem' RefDistr (User Data) was created or updated.");
				for (i = 0; i < nelems; i++) {
					IfmSetElementalRefDistrValue(pDoc, rdidcaE, i, Ptopelems[i] ? Pclayacc_elemtal[i] : -1.0);
				}
			}
		}
		else {
			long rdidcaN = IfmGetNodalRefDistrIdByName(pDoc, "clayaccumnodal"); //ref. distrib. index for the 'clayaccumnodal' nodal ref. distrib. (@User Data) (in PostSimulation)
			if (rdidcaN < 0) rdidcaN = IfmCreateNodalRefDistr(pDoc, "clayaccumnodal");
			if (rdidcaN < 0) {
				IfmWarning(pDoc, "[PostSimul] 'clayaccumnodal' RefDistr (User Data) could not be created for some reason.");
			}
			else {
				IfmInfo(pDoc, "[PostSimul] 'clayaccumnodal' RefDistr (User Data) was created or updated.");
				for (i = 0; i < nnodes; i++) {
					IfmSetNodalRefDistrValue(pDoc, rdidcaN, i, Ptopnodes[i] ? Pclayacc_nodal[i] : -1.0);
				}
			}
		}

		long rdidfd = IfmGetNodalRefDistrIdByName(pDoc, "flooddur"); //ref. distrib. index for the 'flooddur' nodal ref. distrib. (@User Data) (in PostSimulation)
		if (rdidfd < 0) rdidfd = IfmCreateNodalRefDistr(pDoc, "flooddur");
		if (rdidfd < 0) {
			IfmWarning(pDoc, "[PostSimul] 'flooddur' RefDistr (User Data) could not be created for some reason.");
		}
		else {
			IfmInfo(pDoc, "[PostSimul] 'flooddur' RefDistr (User Data) was created or updated.");
			for (i = 0; i < nnodes; i++) {
				IfmSetNodalRefDistrValue(pDoc, rdidfd, i, Ptopnodes[i] ? Pfloodedduring[i] : 0.0);
			}
		}
	}

	//Close export
	int ncerr = 0;
	if (doexportdata) {
		ncerr = netcdf_close_file(pDoc, ncdf_globtopN);
		if (ncerr) IfmWarning(pDoc, "[Export @ PostSimul] There was an error while trying to close the NetCDF file...");
	}

	//Deallocate arrays common to ALL MODES:
	delete[] Ptopnodes;
	delete[] Ptopelems;
	delete[] Pcurrfloodednodes;
	delete[] Pinflownodes;
	delete[] PallnodesY;

	//Deallocate arrays for the optional data export:
	if (doexportdata) {
		if(Pnodalmassrbudgets != nullptr) delete[] Pnodalmassrbudgets;
		if (!altexportusernodes) {
			delete[] PglobXiraw;
			delete[] Pcurrglobtopnodes;
			delete[] Pcurrglobtop_isinflow;
			delete[] Pprevglobtop_isinflow;
			delete[] Pcurrglobtop_nodalmassrbudgets;
		}
		else {
			delete[] Puexportnodes;
			delete[] PuserexportnodesA;
			delete[] Puserexportnodes_uvaluesA;
			if (altexportautoinclgtop) {
				delete[] PeffexportnodesA;
				delete[] Peffexportnodes_uvaluesA;
			}
		}
	}

	if(isaglobxiEloaded) delete[] PglobXiElemA;

	//Deallocate arrays specific to OLDER MODES:
	delete[] Ptopelemszmid;
	delete[] Pclayacc_nodal;
	delete[] Pclayaccfinal_nodal;
	delete[] Pclayacc_elemtal;
	delete[] Pclaysedrate_nodal;
	delete[] Pclaysedrate_elemtal;
	delete[] Pfloodedduring;
	delete[] PelemsAtnodesBR;
	delete[] PelemsAtnodesBL;
	delete[] PelemsAtnodesBM;

	if (withclaymesh) {
		//Operations for the new dynamic clay mode only:
		delete[] PclayYimain;
		delete[] PclayYisub;
		delete[] PclayXiraw;
		delete[] Pclaylayeri;
		delete[] PtopclayBCnodes;
		delete[] PtoprockBCnodes;
		delete[] PupdateNbelowclaytop;
		delete[] PupdateNoverclaytop;
		delete[] PcleanNBCprevclaytop;

		delete[] Pcurrclaytopnodes;
		delete[] Pprevclaytopnodes;
		delete[] Pcurrclaytopisinflow;

		delete[] Pclaytopnodes_tmp_prevtophbc;
		delete[] Pclaytopnodes_tmp_currtopnewhbc;
		delete[] Pclaytopnodes_tmp_prevtopY;
		delete[] Pclaytopnodes_tmp_currtopY;
		delete[] Pclaytopnodes_tmp_prevtopcbc;
		delete[] Pclaytopnodes_tmp_currtopnewcbc;
		delete[] Pclaytopnodes_tmp_currtopcbc_chgneeded;

		delete[] Pallnodes_heads;
		delete[] Pallnodes_mconcs;

		delete[] Pallcurractivenodes;
		delete[] Pallcurrclayactivenodes;
		delete[] Psteadyallrocknodes;
		delete[] Psteadytoprocknodes;

		if (BETAwithsourceterms) delete[] PclaysrcBetaCoefE;
		delete[] PelemcurractiveE;
	}

	IfmInfo(pDoc, "<< postSimulation done >> (paused, stopped or finished).");
}

//NEW BETA internal options for PreTimeStep:
const bool withhgrad = true; //internal hardwired option (for the withclaymesh mode only)
const bool nodata_conc_overclaytop = true; //now testing "true"; ... ... before, I've been testing "false" [BUT NEED FEFLOW SUPPORT!!!] since Fluid density appears not to be updated adequately after c of newly activated nodes is initialized!
const bool withcgrad = true; //NEW BETA internal hardwired option (for the withclaymesh mode only)

void Cpaleosea2d::PreTimeStep (IfmDocument pDoc)
{
	int i;
	char txtbuffer[300], simtimetxtbuff[60], sedimtxtbuff[50], ioflowtxtbuff[80], clayeqkvtxtbuff[50];

	double rslval; //RSL value at current time (extracted from rsl TS)
	double swconc; //seawater concentration global value (in mg/L) at current time (extracted from sal TS)

	double zval; //Y_local value at the node
	double hbcval, hbcvaleff; //equiv. fw head at the sea bottom if flooded, else h=z, computed per node; "eff" for effectively applied to the inside nodes when withhgrad==true
	double csurfwval; //concentration local value (in mg/L) at the "model surface", i.e. at sea/lake bottom if flooded, else that of the recharge water (=fwconc); computed per node (noise-free, original value from the time series or default constants)
	double cbcval, cbcvaleff; //concentration local value (in mg/L) at the top nodes of the model (which may include a retardation effect by the clay aquitard if enabled and present, in NON-withclaymesh modes where model top is that of the rock aquifer even in the confined parts); computed per node (noise may be added, also per node!); "eff" for effectively applied to the inside nodes when withcgrad==true
	double cbcval_if_outflow; //NEW BETA... TODO write description based on code effect

	bool addnoise; //Should noise be added to the concentration BCs at current top nodes? (normally a time-dep. optional based on time series 'rndconc') (but: always disabled if !domasstransport)

	long elemAtnodeBR, elemAtnodeBL, elemAtnodeBM;
	double computedKveq;
	double computedTR;

	bool isflooded; //is the node flooded at current time?
	double currtime; //absolute simulation time (in days) at the START of the current time step [Yes, verified 25-10-2016]

	currtime = IfmGetAbsoluteSimulationTime(pDoc); //(time at time step start)
	//Extracting time series...
	rslval = IfmInterpolatePowerValue(pDoc, rslPID, currtime);
	currsedimrate = withclayTimeDepsedrate ? IfmInterpolatePowerValue(pDoc, sedratePID, currtime) : 0.0;
	swconc = withtransientswconc ? IfmInterpolatePowerValue(pDoc, swcPID, currtime) : defswconc;
	addnoise = domasstransport ? (withrndnoiseaddtoconcBCs ? IfmInterpolatePowerValue(pDoc, rndPID, currtime)>int_tol : currtime<defrandomduration) : false;

	if (withclaymesh) {
		//Updating the geometrical state of the clay aquitard (if withclaymesh mode is active)
		currclayYi_rawreal = IfmInterpolatePowerValue(pDoc, clayYiPID, currtime);
		currclayYi = floor(currclayYi_rawreal + int_tol);
		UpdateClayParamsAccToYi(pDoc, initialpreTSupdate, currtime);
		/* As a result, Ptopnodes is updated if necessary, PupdateNbelowclaytop is non-all-false if relevant,
		   and clay elements are activated/deactivated as required. Also, currKofclays value is updated. */
		relevanttimestep = !hasclayYiincreased;
	}
	else {
		//Otherwise, in the older modes, these are set to default constant values
		currclayYi_rawreal = nodata_double;
		currclayYi = -1;
		ncurractivenodes = nnodes;
		needupdatearoundclaytop = false;
		hasclayYichanged = false;
		hasclayYiincreased = false;
		relevanttimestep = true;
	}

	/* COMMON Initialization or Update of in/outflow array and nodal counts */
	UpdateInflowDataForTopNodes(pDoc, initialpreTSupdate);	//As a result, Pinflownodes array, ninflownodes and noutflownodes, as well as Pcurrglobtop_isinflow array, are updated.

	//Verify if update of BCs is needed
	bool sedimactive = withclaymesh || (withclayConstsedrate ? (currtime < (claysedimduring + time_tol)) : (withclayTimeDepsedrate ? (currsedimrate > (0.0 + sedrate_tol)) : false));
	bool needupdate_gen = true; //BETA: Améliorer (ou retirer!?); ne devrait pas contourner l'update des cfixées selon inflow/outflow nodes! etc. //(abs(rslval-lastrslval)>len_tol) | (abs(swconc-lastswconc)>conc_tol) | sedimactive | addnoise | firstTSresumed;

	//simtimetxtbuff
	if (currtime - simulstarttime < daysperyear) {
		/* the starting, first year simulated */
		sprintf_s(simtimetxtbuff, 60, "%8.4f d", currtime);
	}
	else {
		/* then the rest... */
		sprintf_s(simtimetxtbuff, 60, "%8.2f a", currtime / daysperyear);
	}

	//sedimtxtbuff
	if (sedimactive) {
		if (withclaymesh) {
			sprintf_s(sedimtxtbuff, 50, "clayyi(t)=%6.3f (%d)", currclayYi_rawreal, currclayYi);
		}
		else if (withclayConstsedrate) {
			strcpy_s(sedimtxtbuff, 50, "loc.cst.sed.rate");
		}
		else if (withclayTimeDepsedrate) {
			sprintf_s(sedimtxtbuff, 50, "glob.sedrate(t)=%.3g m/d", currsedimrate);
		}
		else {
			strcpy_s(sedimtxtbuff, 50, "UNKNOWN sedim mode!?");
		}
	}
	else {
		strcpy_s(sedimtxtbuff, 50, "no clay sedim.");
	}
	//clayeqkvtxtbuff
	if (withclaymesh && (withtimedepKclays || updateClayConstantKaswell)) {
		if (currclayYi > 0) {
			sprintf_s(clayeqkvtxtbuff, 50, "; Kclays(t)= %.3e m/s", currKofclays / secsperday);
		}
		else {
			strcpy_s(clayeqkvtxtbuff, 50, "; no clay elems yet");
		}
	}
	else {
		strcpy_s(clayeqkvtxtbuff, 50, "");
	}
	//ioflowtxtbuff
	if (isglobXipresent) {
		sprintf_s(ioflowtxtbuff, 80, "nGtop=%d (in=%d, out=%d, chg=%d)", nbglobXivals, ninflownodes_gtop, noutflownodes_gtop, nflowchgnodes_gtop);
	}
	else {
		sprintf_s(ioflowtxtbuff, 80, "ntop=%d (in=%d, out=%d, mBCc=%d)", ntopnodes, ninflownodes, noutflownodes, nbcctopnodes);
	}
	//whole line of text to display in the Log
	sprintf_s(txtbuffer, 300, "[PreTS] %s BCs; t= %s; rsl= %.3f m, sal= %6.3f g/L; %s; %s; %s%s", needupdate_gen ? "Updating" : "No chg in", simtimetxtbuff, rslval, swconc / 1000.0, ioflowtxtbuff, addnoise ? "noisy c BCs" : "noise-free BCs", sedimtxtbuff, clayeqkvtxtbuff);
	IfmInfo(pDoc, txtbuffer);

	if (false) {
		//BETA STOP [TO BE REMOVED!]
		IfmAlert(pDoc, NULL, "  OK  ", "BETA STOPPING! (early in PreTimeStep)");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	//Variables for withbotclayconc MODE (beta)
	double currclayL; //'corrected' local clay length (thickness) at current ith top node, to avoid errors when extracting from the botclayconc TS with transformed time
	double relativtime = (currtime - endofmaxsalinitytime) / daysperyear; //relative current time (in years) beyond (>) endofmaxsalinitytime
	double transftime; //transformed time (in years) to extract from the botclayconc TS for thicknesses different to the one that was used to generate the TS values
	double transfabstime; //transformed absolute time (in days) computed from transftime right after update
	double suggbclayconcval; //suggested concentration (in mg/L) at bottom of clay (local)

	if (!firstTSresumed && initialpreTSupdate && domasstransport) {
		IfmInfo(pDoc, "* BETA: [PreTS 0] Initialization of the mass conc. state at top nodes (once, at this first time step). [NEW; does it improve early budget(t)?]");
	}
	
	if (needupdate_gen) {
		int xiindex; //clayXi RELATIVE value, with minclayXival as the reference, and starting at 0; internal 1d-array index to manage special models where minclayXival > 0
		int xival; //clayXi ABSOLUTE value, as specified in the nodal user data
		bool somexi;
		bool isinflowNi; //temp. boolean: is it or should it be an inflow at current top node[i] in the for-loop? This is based on Pinflownodes most often, except when hasclayYichanged (in withclaymesh mode) where the clay-top inflow conditions saved from the previous clay-top nodes are used instead
		bool doforceDiriBCNi; //NEW BETA temp. boolean: should the top node BC be forced to be a Dirichlet BC? (especially relevant in the withclaymesh mode)

		for (i = 0; i < nnodes; i++) {
			//Storing all nodal head anc conc. values if clayYi is changing (in withclaymesh mode only)
			if (withclaymesh && hasclayYichanged) {
				Pallnodes_heads[i] = IfmGetResultsFlowHeadValue(pDoc, i);
				Pallnodes_mconcs[i] = domasstransport ? IfmGetResultsTransportMassValue(pDoc, i) : 0.0;
			}
			if (Ptopnodes[i]) {
				//Reading horizontal in-clay node index (for the withclaymesh mode only):
				if (withclaymesh) {
					xival = PclayXiraw[i];
					somexi = (xival >= 0);
					if (somexi) xiindex = xival - minclayXival;
				}
				
				//Flooded state detection with some tolerance: > z ± 0.1 mm
				zval = PallnodesY[i];
				isflooded = rslval > (zval + rsl_tol);
				
				//New c value for the ith node (without noise), first estimated without the botclayconc time series:
				csurfwval = isflooded ? swconc : fwconc;
				//Then, a better estimate if botclayconc time series is available and if the mode requires it:
				if (!withclaymesh && withbotclayconc && (Pclayaccfinal_nodal[i] > len_tol) && (currtime > endofmaxsalinitytime)) {
					//Thickness to use (=Lextrap in R)
					currclayL = Pclayaccfinal_nodal[i];
					if (currclayL < 2.0) currclayL = 2.0; //a lower limit to L, to avoid too high transftime values; TODO: Could be moved just after preparation of Pclayaccfinal_nodal
					transftime = pow(pow(relativtime, 3.0 / 5.0)*origclayL / currclayL, 5.0 / 3.0);
					transfabstime = endofmaxsalinitytime + transftime*daysperyear;
					suggbclayconcval = transfabstime < maxtransftimebclaycTS ? IfmInterpolatePowerValue(pDoc, bclaycPID, transfabstime) : fwconc;
					cbcval = suggbclayconcval > csurfwval ? suggbclayconcval : csurfwval; //so that effective local concentration at top rock is never lower than current conc. of the flooding water
					/* if(i == 2127) { //CODE FOR VERIFYING THE BEHAVIOR OF THIS SECTION or DEBUGGING
					sprintf_s(txtbuffer, 180, "BETA DEBUG bclay: origclayL=%.3f, currclayL=%.4f m, reltime=%.3f a, transftime=%.3f a, suggc=%.3f mg/L, cbcval=%.3f mg/L. OK?",
					origclayL, currclayL, relativtime, transftime, suggbclayconcval, cbcval);
					IfmWarning (pDoc, txtbuffer);
					} */
				}
				else {
					cbcval = csurfwval;
				}

				//New equivalent freshwater head value for the ith node (with the new c value without prior to addition of noise)
				hbcval = isflooded ? rslval + ((rslval - zval) * (csurfwval - minconc) / (maxconc - minconc) * densr) : zval;
				if (withclaymesh && somexi) Pclaytopnodes_tmp_currtopnewhbc[xiindex] = hbcval;

				/* DEBUG EVOL OF transient distrib. parameters
				if(i == 7291) { //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!************************************RETIRER DEBOGAGE*******************!!!!!!!!!!!!!!!!!!!!
				sprintf_s(txtbuffer, 180, "DEBUG TMP: time= %.2f a, rsl= %.3f m, sal= %.3f mg/L, z@node= %.3f m, head@BC= %.3f m (conc@BC= %.3f mg/L; maxconc= %.3f) ***", currtime/daysperyear, rslval, swconc, zval, hbcval, csurfwval, maxconc);
				IfmInfo(pDoc, txtbuffer);
				}
				*/
				//Noise is added to the c value for the ith node (NOISY), after equiv. fw head computation
				if (addnoise) cbcval = cbcval + 0.01 * cbcval * (dist(gen) - 0.5); //random uniform noise is added on fixed conc. BC values: ± 0.5%

																					  /* DEBUG nodal assignments
																					  cbcctype = IfmGetBccMassType(pDoc, i);
																					  cbccset = IfmIsBccMassSet(pDoc, i, 0);
																					  sprintf(txtbuffer,"time= %f, rsl= %f, node= %d, zval= %f, massbcct= %d, cbccset= %d", currtime, rslval, i+1, zval, cbcctype, cbccset);
																					  IfmInfo(pDoc, txtbuffer);
																					  */

				/* ASSIGNEMENT OF FLUID FLOW related BCs (and Material properties) */
				//ACTIONS for the OLDER MODES only:
				if (!firstTSresumed && !withclaymesh) {
					if (withSOMEclayaccumANYmode) {
						//setting Cauchy (3rd kind) b.c. for flow
						IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_CAUCHY, 0, hbcval);
						//updating the transfer rates (in & out) --- nodal mode (unwise)
						if (!wiseTRupdate) {
							elemAtnodeBR = PelemsAtnodesBR[i];
							computedKveq = (Cauchybinit + Pclayacc_nodal[i]) / (Cauchybinit / CauchyKinit + Pclayacc_nodal[i] / CauchyKclays);
							computedTR = computedKveq / (Cauchybinit + Pclayacc_nodal[i]);
							if (elemAtnodeBR >= 0) {
								IfmSetMatFlowTransferIn(pDoc, elemAtnodeBR, computedTR);
								IfmSetMatFlowTransferOut(pDoc, elemAtnodeBR, computedTR);
							}
							if (!meshquadtype) {
								//THIS SUB-SECTION CODE is FAR from ideal!!
								elemAtnodeBL = PelemsAtnodesBL[i];
								if (elemAtnodeBL >= 0) {
									IfmSetMatFlowTransferIn(pDoc, elemAtnodeBL, computedTR);
									IfmSetMatFlowTransferOut(pDoc, elemAtnodeBL, computedTR);
								}
								elemAtnodeBM = PelemsAtnodesBM[i];
								if (elemAtnodeBM >= 0) {
									IfmSetMatFlowTransferIn(pDoc, elemAtnodeBM, computedTR);
									IfmSetMatFlowTransferOut(pDoc, elemAtnodeBM, computedTR);
								}
							}
						}
					}
					else {
						//setting Dirichlet (1st kind) b.c. for flow - That is the simpler case with no clay accum.
						IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_DIRICHLET, 0, hbcval);
					}
				}
				//ACTIONS for the NEW withclaymesh MODE only
				// (yes, only this mode, because older modes use Cauchy flow BC's!):
				if (!firstTSresumed && withclaymesh) {
					//setting Dirichlet (1st kind) b.c. for flow
					IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_DIRICHLET, 0, hbcval);
				}
				
				/* ASSIGNEMENT OF MASS TRANSPORT related BCs (and state values) */
				//COMMON ACTIONS for all MODES (NEW | OLD): mainly assignement of BCs to the top nodes
				isinflowNi = Pinflownodes[i];
				doforceDiriBCNi = BETApermanentDirichletMassConcBCs ||
					(withclaymesh && (true || BETApermanentDirichletMassConcBCs_LateAtClayTop) && somexi && currclayYi_rawreal >= clayYi_FixMassBC_thresv); //here, somexi means it requires a "clay-top" node
				if (!firstTSresumed && domasstransport) {
					if (BETAautomUseMassBCconstraints) {
						//DOES NOT WORK: IfmSetBccMassValueAtCurrentTime(pDoc, i, 0.0, IfmMIN_BCC_TYPE);
						if (currtime < time_tol) {
							//TRYING initial Bcc setting only with constant conc. value -- 26 Oct.
							IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_DIRICHLET, 0, cbcval);
							IfmSetBccMassTypeAndValueAtCurrentTime(pDoc, i, IfmBCC_DIRICHLET, 0, 0.0, IfmMIN_BCC_TYPE);
						} //else do NOTHING HERE!
					}
					else {
						if (isinflowNi || doforceDiriBCNi) {
							/* LOCAL FLUID INFLOW condition (based on Darcy Y velocity) */
							/* A. NEW: Initialization of the conc. state at top nodes at first time step */
							if (initialpreTSupdate && true) { //(in addition to the prior IFs (!firstTSresumed && (isinflowNi || BETApermanentDirichletMassConcBCs))...)
								//BETA BETA! : under test 28oct (Once validated, [TODO] please remove the 'true' from the IF above...)
								IfmSetResultsTransportMassValue(pDoc, i, cbcval);
								IfmSetResultsTransportMassIntermediateValue(pDoc, i, cbcval); //NEW, testing
								IfmSetResultsTransportMassPreviousTimeValue(pDoc, i, cbcval);
							}
							/* B. Setting Dirichlet (1st kind) b.c. for mass transport in case of local inflow */
							IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_DIRICHLET, 0, cbcval);
						}
						else {
							/* LOCAL FLUID OUTFLOW condition (based on Darcy Y velocity) */
							if (!divergformtransp) {
								/* (with the default 'Convective form' of the transport equation) */
								//!!!!!!!! TODO-TRYING 4th kind: Better Mass Budget than 'no BC' ?? (for Convective form) !!!!!!!!
								//Setting Well (4th kind) b.c. for mass transport in case of local outflow (~free outflow BC):
								IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_SINGLE_WELL, 0, 0.0);
								//DISABLED forTESTING: IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0.0);
							}
							else {
								/* (with the alternative 'Divergence form' of the transport equation) */
								IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0.0);
							}
						}
					}
				}
				Pcurrfloodednodes[i] = isflooded; //(updated in all modes)

				//Temporarily storing h and c data at top nodes to copy-paste these values to newly activated nodes below
				// (note that this IF runs only for topnodes; see the upper IF, further above)
				if (withclaymesh) { //TODO ~25oct: VERIFY if this section can be skipped when !hasclayYiincreased ???
					if (somexi) {
						//Preparing h and Y data necessary for performing linear interpolation (by a simple "rule of three")
						// of head values between previous and new clay-top nodes. BETA!
						// Note that this procedure is not relevant in case of decreasing clayYi, as current values of
						//  the nodes below a new-but-lower clay-top are already the better values. 
						if (hasclayYiincreased) {
							//TODO: Could use my new internal array Pallnodes_heads instead of the IFM extract method...
							// and also internal steady arrays for Y coords...
							Pclaytopnodes_tmp_prevtophbc[xiindex] = IfmGetResultsFlowHeadValue(pDoc, Pprevclaytopnodes[xiindex]);
							Pclaytopnodes_tmp_prevtopcbc[xiindex] = domasstransport ? IfmGetResultsTransportMassValue(pDoc, Pprevclaytopnodes[xiindex]) : 0.0;
							Pclaytopnodes_tmp_prevtopY[xiindex] = PallnodesY[Pprevclaytopnodes[xiindex]];
							Pclaytopnodes_tmp_currtopY[xiindex] = PallnodesY[i];
						}
						//Updating concentration ("c") values to be assigned to the current-clay-top nodes...
						// ...and to newly activated nodes just below, when clayYi has changed.
						// Note that the value used depends on the local fluid flow direction:
						//  if qY > 0 (i.e. upward, outflow), every node at this X gets the c value from the previous top for this X;
						//  if qY < 0 (i.e. downward, inflow), every node at this X gets the new c value from the swconc(t) time series.
						Pclaytopnodes_tmp_currtopcbc_chgneeded[xiindex] = withcgrad || isinflowNi || hasclayYiincreased; //NEW MODIF, BETA!
						//REPLACED BY BELOW: cbcval_if_outflow = (!isinflowNi && hasclayYiincreased) ? IfmGetResultsTransportMassValue(pDoc, Pprevclaytopnodes[xiindex]) : nodata_double;
						cbcval_if_outflow = withcgrad ? cbcval : ((!isinflowNi && hasclayYiincreased) ? Pclaytopnodes_tmp_prevtopcbc[xiindex] : nodata_double); //NEW VERSION, BETA!
						Pclaytopnodes_tmp_currtopnewcbc[xiindex] = isinflowNi ? cbcval : cbcval_if_outflow;
					}
				} //(end of this 'withclaymesh' subsection)
			} //(end of the long if topnode)
		} //(end of nodal for loop)

		//BETA section UNDER TEST: Initializing state values of the new clay-top nodes, right after clay elements were activated
		if (withclaymesh && hasclayYiincreased) {
			for (xiindex = 0; xiindex < nbclayXivals; xiindex++) {
				//NOTE: Here, xiindex acts like a clayXi value RELATIVE to the minimum clayXi value.
				xival = minclayXival + xiindex;
				i = Pcurrclaytopnodes[xiindex];
				hbcval = Pclaytopnodes_tmp_currtopnewhbc[xiindex]; //IMPORTANT!
				//initializing state values of the newly activated nodes of the new CLAY-TOP (yes, only the nodes of the new top!)
				// (conc. first, head after, so that local fluid density is correctly estimated; LIKELY NOT NECESSARY SEQUENCE...)
				if (false) {
					//Earlier approach (deprecated)
					if (Pclaytopnodes_tmp_currtopcbc_chgneeded[xiindex] && domasstransport) IfmSetResultsTransportMassValue(pDoc, i, Pclaytopnodes_tmp_currtopnewcbc[xiindex]); //TODO-someday: Verify how that works when clayyi(t0) starts with a value far from 0... because I think that currently it does not manage the "initial" assignement of state values to abruptly activated elements (first PreTimeStep right after PreSimulation initialization)...
					IfmSetResultsFlowHeadValue(pDoc, i, hbcval);
					IfmSetResultsFlowHeadPreviousTimeValue(pDoc, i, hbcval);
				}
				else {
					//Recommended approach (acc. to Carlos)
					if (Pclaytopnodes_tmp_currtopcbc_chgneeded[xiindex]) Pallnodes_mconcs[i] = Pclaytopnodes_tmp_currtopnewcbc[xiindex];
					Pallnodes_heads[i] = hbcval;
				}
			}
			xiindex = -1; //BOF: just to make sure the variable is no more a for-loop iterator...
		}

		/* For the new withclaymesh mode only:
		    overwriting state of nodes around the current clay-top, and cleaning of previous BCs */
		if(withclaymesh && needupdatearoundclaytop) {
			//TODO: Make the text 'and Mass transport' dependent on mode domasstransport = f/t
			std::string txtsubstr1;
			txtsubstr1 = domasstransport ? " and Mass transport" : "";
			sprintf_s(txtbuffer, 240, "[WCM @ PreTS] Fluid flow%s nodal state and BC values are being updated due to chg clay top (clayyi from %d to %d).", txtsubstr1.c_str(), prevclayYi, currclayYi);
			IfmInfo(pDoc, txtbuffer);
			for (i = 0; i < nnodes; i++) {
				if (PupdateNbelowclaytop[i]) {
					xival = PclayXiraw[i];
					somexi = (xival >= 0);
					if (somexi) {
						//translating the clayXi value into its corresponding index (i.e. relative) value
						xiindex = xival - minclayXival;
						//calculating hbcvaleff by linear interpolation (only if withhgrad)
						hbcval = Pclaytopnodes_tmp_currtopnewhbc[xiindex];
						if (withhgrad) {
							hbcvaleff = Pclaytopnodes_tmp_prevtophbc[xiindex] + (hbcval - Pclaytopnodes_tmp_prevtophbc[xiindex]) * (PallnodesY[i] - Pclaytopnodes_tmp_prevtopY[xiindex]) / (Pclaytopnodes_tmp_currtopY[xiindex] - Pclaytopnodes_tmp_prevtopY[xiindex]);
						}
						else hbcvaleff = hbcval;
						//calculating cbcvaleff by linear interpolation (only if withcgrad) --- BETA!!
						cbcval = Pclaytopnodes_tmp_currtopnewcbc[xiindex];
						if (withcgrad) {
							cbcvaleff = Pclaytopnodes_tmp_prevtopcbc[xiindex] + (cbcval - Pclaytopnodes_tmp_prevtopcbc[xiindex]) * (PallnodesY[i] - Pclaytopnodes_tmp_prevtopY[xiindex]) / (Pclaytopnodes_tmp_currtopY[xiindex] - Pclaytopnodes_tmp_prevtopY[xiindex]);
						}
						else cbcvaleff = cbcval;
						//NEW NOISE ADDED, UNDER TEST (line below). PROBLEM IS, THAT IS TWICE THE NOISE of the BCs...!
						if (addnoise) cbcvaleff = cbcvaleff + 0.01 * cbcvaleff * (dist(gen) - 0.5); //random uniform noise is added on fixed conc. BC values: ± 0.5%
						//state values (conc. first, head after, so that local fluid density is correctly estimated)
						if (false) {
							//Earlier approach (deprecated) (1)
							if (Pclaytopnodes_tmp_currtopcbc_chgneeded[xiindex]) IfmSetResultsTransportMassValue(pDoc, i, Pclaytopnodes_tmp_currtopnewcbc[xiindex]); //TODO-someday: Verify how that works when clayyi(t0) starts with a value far from 0... because I think that currently it does not manage the "initial" assignement of state values to abruptly activated elements (first PreTimeStep right after PreSimulation initialization)...
							IfmSetResultsFlowHeadValue(pDoc, i, hbcvaleff);
							IfmSetResultsFlowHeadPreviousTimeValue(pDoc, i, hbcvaleff);
						}
						else {
							//Recommended approach (acc. to Carlos)
							if (Pclaytopnodes_tmp_currtopcbc_chgneeded[xiindex]) Pallnodes_mconcs[i] = cbcvaleff;
							Pallnodes_heads[i] = hbcvaleff;
						}
					}
					//BC's
					IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0);
					if(domasstransport) IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0);
				}
				if (PupdateNoverclaytop[i]) {
					//state values --> NoData
					if (false) {
						//Earlier approach (deprecated)
						if (nodata_conc_overclaytop && domasstransport) IfmSetResultsTransportMassValue(pDoc, i, nodata_double);
						IfmSetResultsFlowHeadValue(pDoc, i, nodata_double);
					}
					else {
						//Recommended approach (acc. to Carlos)
						if (nodata_conc_overclaytop) Pallnodes_mconcs[i] = nodata_double;
						Pallnodes_heads[i] = nodata_double;
					}
					//BC's
					IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0);
					if(domasstransport) IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0);
				}
				if (PcleanNBCprevclaytop[i]) { //should concern only nodes below the (new) current clay-top, normally...
					//BC's
					IfmSetBcFlowTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0);
					if (domasstransport) IfmSetBcMassTypeAndValueAtCurrentTime(pDoc, i, IfmBC_NONE, 0, 0);
				}
			} //(end of that nodal for loop)
		}
		
		//Older mode with wise Transfer Rates... (so it is skipped if withclaymesh)
		if (!withclaymesh && wiseTRupdate && withSOMEclayaccumANYmode) {
			for (i = 0; i < nelems; i++) {
				if (!Ptopelems[i]) continue;
				computedKveq = (Cauchybinit + Pclayacc_elemtal[i]) / (Cauchybinit / CauchyKinit + Pclayacc_elemtal[i] / CauchyKclays);
				computedTR = computedKveq / (Cauchybinit + Pclayacc_elemtal[i]);
				IfmSetMatFlowTransferIn(pDoc, i, computedTR);
				IfmSetMatFlowTransferOut(pDoc, i, computedTR);
			}
		}

		//store current time series values as 'previous' values for the next iteration
		lastrslval = rslval;
		lastswconc = swconc;
		lastsedimrate = currsedimrate;
	}
	firstTSresumed = false;
	initialpreTSupdate = false; //so that the @PreTimeStep initializations happen only once.

	if (false && needupdatearoundclaytop && currclayYi > 0) {
		IfmAlert(pDoc, NULL, "  OK  ", "STOPPING TO SEE Values after updated claytop");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (false) {
		//BETA STOP [TO BE REMOVED!]
		IfmAlert(pDoc, NULL, "  OK  ", "BETA STOPPING! (at the end of PreTimeStep)");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if(withclaymesh) prevclayYi = currclayYi;
	firstflowsol = true;
	firstmasssol = true;
}

void OverwriteAllNodesWithArrayValues(IfmDocument pDoc, bool heads=true, bool concs=true, bool prevvals=true, bool veloczero=false)
{
	if (!(withclaymesh && hasclayYichanged)) return; //SAFETY...
	int i;
	//REMOVED: const double multforprev = 0.9999;
	for (i = 0; i < nnodes; i++) {
		if (domasstransport) {
			if (concs) IfmSetResultsTransportMassValue(pDoc, i, Pallnodes_mconcs[i]);
			if (concs && prevvals) {
				IfmSetResultsTransportMassIntermediateValue(pDoc, i, Pallnodes_mconcs[i]);
				IfmSetResultsTransportMassPreviousTimeValue(pDoc, i, Pallnodes_mconcs[i]);
			}
		}
		if (heads) IfmSetResultsFlowHeadValue(pDoc, i, Pallnodes_heads[i]);
		if (heads && prevvals) IfmSetResultsFlowHeadPreviousTimeValue(pDoc, i, Pallnodes_heads[i]);
		if (veloczero) {
			//UNDER TEST !!!
			IfmSetXVelocityValue(pDoc, i, 0.0);
			IfmSetYVelocityValue(pDoc, i, 0.0);
		}
	}
	char txtbuffer[180];
	sprintf_s(txtbuffer, 180, "[WCM @ *many*] Applying OverwriteAllNodesWithArrayValues(heads=%d, concs=%d, prevvals=%d, veloczero=%d) with plugin 'Pallnodes...' arrays", heads, concs, prevvals, veloczero);
	IfmInfo(pDoc, txtbuffer);
	/* Programmer's notes:
	I've tested many things trying to prevent the solver from bugging between Pre- and PostMassSimulation, without success.
	Assigning previous values from Pallnodes_* arrays multiplied by a non-unity coefficient (e.g. 0.9999) did nothing.
	Replacing the NaN velocity values at newly activated nodes did nothing either.
	In conclusion, it seems I cannot avoid the convergence failure for the time step when clayYi change is applied.
	That's not a real problem since I force FEFLOW to repeat (do again) this same time step once elements and nodes
	are activated and correctly initialized (as desired). Written on 2016-10-23.
		On another aspect, the prevvals option was found to be very efficient in boosting performances of the solver
	when clayYi changes. When it is disabled, previous solution values = NaN make it very hard for the solver, which
	has to do a lot of iterations to rebuild some ~valid values internally. In addition to these undesired iterations,
	it results in a dramatically decreased time step (dt), e.g. from 0.025a to 1e-10a, which has thus to be avoided!
	In short, prevvals is a great option to continue to use!
	*/
}

void ExtractYvelocAtTopToIdentifyInflowNodes(IfmDocument pDoc)
{
	int i;
	int xiindex; //clayXi RELATIVE value, with minclayXival as the reference, and starting at 0; internal 1d-array index to manage special models where minclayXival > 0
	int xival; //clayXi ABSOLUTE value, as specified in the nodal user data
	int gxival;
	int gxiindex;
	bool somexi;
	double nodYveloc;
//	double x, y, z = 0.0; //only z is initialized here
//	IfmBool validveloc;
	char warntxtbuff[180];
	if (false) {
		sprintf_s(warntxtbuff, 180, "[PostTS] IfmIsTimeStepRejected = %d", IfmIsTimeStepRejected(pDoc));
		IfmInfo(pDoc, warntxtbuff);
	}

	if (withclaymesh && hasclayYichanged) {
		//In such case, inflow info arrays are not updated, thereby keeping info of the previous clay-top (i.e. before current ongoing clayYi change).
		// Therefore, 'initial' inflow 'state' at each node of the new top will be based on Pcurrclaytopisinflow[].
		// Also note that this boolean condition is similar to verifying if(!relevanttimestep), although less strict.
		if (false) IfmInfo(pDoc, "*DEBUG-INFO* [ExtractYveloc @ PostTS] Extraction skipped because withclaymesh & hasclayYichanged");
		return;
	}

	bool velocfieldpresent = IfmIsVelocityFieldPresent(pDoc) == True;
	if (!velocfieldpresent) {
		IfmWarning(pDoc, "[ExtractYveloc @ PostTS] Strangely, no velocity field is present at that iteration. Pinflownodes won't be updated. Please make sure that's okay.");
		return;
	}

	for (i = 0; i < nnodes; i++) {
		if (!Ptopnodes[i]) continue;
/* TESTING THE SIMPLER NODAL Velocity Extraction approach : further below 
		x = IfmGetX(pDoc, i);
		y = PallnodesY[i];
		nodYveloc = IfmGetResultsYVelocityValueAtXYZ(pDoc, x, y, z, &validveloc);
*/
		nodYveloc = IfmGetResultsYVelocityValue(pDoc, i);
		//REMOVED DEBUG WARNING:
		/*if (true && needupdatearoundclaytop && currclayYi > 0) {
			sprintf_s(warntxtbuff, 180, "Y Darcy veloc. @ node %d (x,y,z)=(%.3f,%.3f,%.3f) = %g.", i + 1, x, y, z, nodYveloc);
			IfmInfo(pDoc, warntxtbuff);
		} */
/*		
		if (!validveloc) { //console warning for the positive case
			sprintf_s(warntxtbuff, 180, "Y Darcy veloc. could not be extracted @ node %d (x,y,z)=(%.3f,%.3f,%.3f).", i + 1, x, y, z);
			IfmWarning(pDoc, warntxtbuff);
		}
*/
		//First, updating Pinflownodes, sufficient for OLDER MODES, OR when withclaymesh & clayYi is NOT changing.
		if (!withclaymesh || !hasclayYichanged) {
			Pinflownodes[i] = nodYveloc <= 0.0;
		}
		//Then, updating Pcurrclaytopisinflow (for eventual use when clayYi will change...)
		// (but this update is done only when withclaymesh & clayYi is not changing, because we can't infer anything
		//  from the undefined NaN values forced into the transient model by FEFLOW when new elements are activated).
		if (withclaymesh && !hasclayYichanged) {
			xival = PclayXiraw[i];
			somexi = (xival >= 0); //normally should always be true here
			if (somexi) {
				xiindex = xival - minclayXival;
				Pcurrclaytopisinflow[xiindex] = Pinflownodes[i];
			}
		}

		//Finally, optional storage of data for eventual data export
		if (isglobXipresent) {
			//NEW, BETA!
			gxival = PglobXiraw[i];
			if (gxival >= 0) {
				gxiindex = gxival - minglobXival;
				Pprevglobtop_isinflow[gxiindex] = Pcurrglobtop_isinflow[gxiindex];
				Pcurrglobtop_isinflow[gxiindex] = Pinflownodes[i];
			}
		}
	}
}

//FOR NOW: BOF... not used
void DoComputationOfOverallBudgets(IfmDocument pDoc)
{
	IfmBudget *fb = IfmBudgetFlowCreate(pDoc);
	int i;
	double totfrb = 0.0;
	double ffluxi;
	for (i = 0; i < nnodes; i++) {
		ffluxi = IfmBudgetQueryFlowAtNode(pDoc, fb, i);
		totfrb = totfrb + ffluxi;
	}
	IfmBudgetClose(pDoc, fb);

}

//PROGRAMMER'S NOTE: This procedure is never called if EXPORT_MASS_RBUDGETs ain't true!
void ForceComputationOfBudgets(IfmDocument pDoc)
{
	int i; //flexible index (can be used for nodal indeces as well as for globXi indeces
	int *idbudgetnodes; //array of nodal indeces...
	double *computedmrbudgets;
	bool *boolvec;
	int bveclen;
	int cnt, truecnt_shouldbe;
	char boolvecname[32];

	/* ---Applied info from FEFLOW Help--- :
	   In transient applications the budget evaluation should be put in a 'PostTimeStep' callback rather than
	    in the 'PostFlowSimulation' because transient budget vectors are yet incomplete at that simulation stage
		and incorrect flux values result!. (in help for: IfmBudgetFlowCreate) */

	/* BEWARE: This procedure does not yet work well for user-export-nodes, in the sense that it may export 
	   nodata or even uninitialized mass-rate-budget values at nodes which are not common with globXiraw top nodes... */

	if (isglobXipresent) {
		/* THE NEW WAY, with global top nodes (i.e. based on Pcurrglobtopnodes array of FF node indeces) */
		// (Be careful, as it replaces what's been prepared by the IF section right before)

		strncpy_s(boolvecname, 32, "Xi-ordered global top nodes", 32);
		truecnt_shouldbe = nbglobXivals;
		cnt = truecnt_shouldbe;

		idbudgetnodes = new int[cnt];
		computedmrbudgets = new double[cnt];

		//Just for copying the current global-top node indeces to the internal temporary array 'idbudgetnodes':
		//(note: order is important here, and is kept)
		for (i = 0; i < cnt; i++) {
			idbudgetnodes[i] = Pcurrglobtopnodes[i];
		}

		if (ntopnodes != truecnt_shouldbe) {
			IfmWarning(pDoc, "[Budg @ PostTS] Plugin internal non-fatal error: Inconsistency: idbudgetnodes count is not equal to ntopnodes!");
		}
	}
	if (false) {
		/* THE OLD WAY, with full array... */

		if (false) {
			/* only the current top nodes */
			strncpy_s(boolvecname, 32, "top nodes (clay & rock)", 32);
			boolvec = Ptopnodes;
			bveclen = nnodes;
			truecnt_shouldbe = ntopnodes;
		}
		if (true) {
			/* all actives nodes */
			strncpy_s(boolvecname, 32, "active nodes (clay & rock)", 32);
			boolvec = Pallcurractivenodes;
			bveclen = nnodes;
			truecnt_shouldbe = ncurractivenodes;
		}

		idbudgetnodes = new int[truecnt_shouldbe];
		computedmrbudgets = new double[truecnt_shouldbe];
		cnt = 0;
		for (i = 0; i < bveclen; i++) {
			if (boolvec[i]) {
				idbudgetnodes[cnt] = i;
				cnt++;
				if (cnt > truecnt_shouldbe) {
					IfmWarning(pDoc, "[Budg @ PostTS] Plugin internal non-fatal error: Inconsistency: idbudgetnodes count became > truecnt_shouldbe!");
					break;
				}
			}
		}
	}

	char txtbuffer[180];
	if (cnt != truecnt_shouldbe) {
		sprintf_s(txtbuffer, 180, "[Budg @ PostTS] Inconsistent current number of budget-nodes: count(t)=%d; vs. ntrue-should-be(t)=%d.", cnt, truecnt_shouldbe);
		IfmWarning(pDoc, txtbuffer);
		if (cnt > truecnt_shouldbe) cnt = truecnt_shouldbe;
	}
	else {
		if (DisplayBudgetComput_infos) {
			sprintf_s(txtbuffer, 180, "[Budg @ PostTS] Flow & Mass transport budgets are now computed for the current %d %s.", cnt, boolvecname);
			IfmInfo(pDoc, txtbuffer);
		}
	}
	IfmBudgetCompute(pDoc, IfmPCLS_MASS_TRANSPORT, idbudgetnodes, cnt, computedmrbudgets, nullptr);
	
	//Copying back the results in the full array... and the global-top array as well (NEW, for optional data export)
	if (isglobXipresent) {
		//Here, i is the index within the filtered array.
		for (i = 0; i < cnt; i++) {
			Pnodalmassrbudgets[idbudgetnodes[i]] = computedmrbudgets[i]; //for the OLD WAYS
			Pcurrglobtop_nodalmassrbudgets[i] = computedmrbudgets[i]; //NEW WAY (BUT not programmed yet!)
		}
	}

	//Deallocation, whatever doexportdata is active or not
	delete[] idbudgetnodes;
	delete[] computedmrbudgets;

	if (DisplayBudgetComput_infos) IfmInfo(pDoc, "[Budg @ PostTS] Budget computations done for this time step.");
}

void Cpaleosea2d::PostTimeStep (IfmDocument pDoc)
{
	int i;
	char txtbuffer[180];

	double currtime = IfmGetAbsoluteSimulationTime(pDoc); //current time (in days) at the END of the current time step
	double currTS = IfmGetCurrentTimeIncrement(pDoc); //current time step length (in days)

	if (false && needupdatearoundclaytop && currclayYi > 0) {
		//BETA STOP [TO BE REMOVED!]
		IfmInfo(pDoc, "calling: PostTimeStep (and aborting right before final overwrite of h & c values)");
		IfmWarning(pDoc, "Do you see the default h and c values?");
		IfmAlert(pDoc, NULL, "  OK  ", "BETA STOPPING @ PostTimeStep");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (withclaymesh && hasclayYichanged && hasclayYiincreased) OverwriteAllNodesWithArrayValues(pDoc);

	if (!withclaymesh) {
		bool elemisflooded;
		if (wiseTRupdate) {
			for (i = 0; i < nelems; i++) {
				if (!Ptopelems[i]) continue;
				elemisflooded = lastrslval > (Ptopelemszmid[i] + rsl_tol); //flooded state detection (elemental) with some tolerance: > z ± 0.1 mm
				Pclayacc_elemtal[i] = withclayConstsedrate ? (currtime < claysedimduring ? Pclaysedrate_elemtal[i] * currtime : Pclaysedrate_elemtal[i] * claysedimduring) : Pclayacc_elemtal[i] + currsedimrate * currTS * (elemisflooded ? 1.0 : 0.0);
				//Note: lastrslval was updated during PreTimeStep so is up-to-date for the current, ending time step.
				//TODO: Make some use of elemisflooded... (like an output User Data Elem. Distrib...)
			}
		}
		for (i = 0; i < nnodes; i++) {
			if (!Ptopnodes[i]) continue;
			if (!wiseTRupdate) Pclayacc_nodal[i] = withclayConstsedrate ? (currtime < claysedimduring ? Pclaysedrate_nodal[i] * currtime : Pclaysedrate_nodal[i] * claysedimduring) : Pclayacc_nodal[i] + (Pcurrfloodednodes[i] ? currsedimrate * currTS : 0.0);
			Pfloodedduring[i] = Pfloodedduring[i] + (Pcurrfloodednodes[i] ? currTS : 0.0);
		}
	}

	//[RECENTLY MOVED, from @PostFlowSimulation; October 18th, 2016]
	//Extraction of vertical velocity Y component nodal values to infer in/outflow state of the top nodes,
	// which will be used for constraining Mass transport BC's for the next time step
	// (always performed, except when withclaymesh & hasclayYichanged; so it is for ALL MODES)
	ExtractYvelocAtTopToIdentifyInflowNodes(pDoc);

	//NEW; ok for ALL MODES
	if(domasstransport && doexportdata && EXPORT_MASS_RBUDGETs) ForceComputationOfBudgets(pDoc); //BETA IF conditions (TODO review them!)
	//also TODO: Consider moving this task into the Effective Export IF section further below, to avoid computing budgets if no such data needs to be exported at current TS.

	//NEW! Excellent & efficient way to repeat the simulation of the time step when it involved a change in active elements
	if (withclaymesh && hasclayYichanged && hasclayYiincreased) IfmSetSimulationControlFlag(pDoc, IfmCTL_REPEAT);
	//BETA messages for the programmer:
	if (!relevanttimestep) {
		//REPLACED by cst text: sprintf_s(txtbuffer, 180, "[@PostTimeStep] relevanttimestep = %d", relevanttimestep);
		IfmInfo(pDoc, "[PostTS] Time step is rejected by the plugin (e.g. to be repeated again with the updated clay mesh).");
	}
	if (withclaymesh && hasclayYichanged && !hasclayYiincreased) {
		IfmWarning(pDoc, "[BETA for the PROGRAMMER @ PostTS] Please verify if time step should be relevant==false also when clayYi decreases. (???)");
	}

	if (relevanttimestep && doexportdata) {
		//EXPORT (in test!)
		int ncerr;
		if (isFirstExportedTS) {
			ncerr = initialize_ncdf_exports(pDoc);
			lastexporttime = currtime;
			if (ncerr != NC_NOERR) {
				IfmWarning(pDoc, "[Export init.] Due to error(s) during its initialization, NetCDF export is forced DISABLED!");
				IfmWarning(pDoc, "** (Note: Related dynamically allocated objects may not be properly deallocated at simul. end.)");
				IfmWarning(pDoc, "** (The most common cause preventing creation of the output file is a MISSING 'results' FOLDER. Please verify this.)");
				//IfmWarning(pDoc, "** (MOREOVER, note that FEFLOW + PaleoSea2D will crash soon because this is not perfectly handled yet!)");
				/*
				    TODO SOON: YET STILL, IT WOULD BE IMPORTANT TO VERIFY IF THE SIMUL. WORKS OKAY EVEN WHEN doexportdata == false. Is it?
				*/
				doexportdata = false;
			}
		}
		if (doexportdata && (isFirstExportedTS || (currtime - lastexporttime) > exportdeltatime)) {
			ncerr = netcdf_writedata(pDoc, ncdf_globtopN);
			if (!isFirstExportedTS) lastexporttime = lastexporttime + exportdeltatime; //To ensure more regular export times; REPLACED = currtime;
			sprintf_s(txtbuffer, 180, "[Export] Current post-time-step results exported to the NetCDF file (at t = %.3f a) (for %zd nodes, %d CZ, and %d MZ) (gminconc = %.4g g/L) (err.code %d).", currtime / daysperyear, ncdf_globtopN.nodeindex_len, ncdf_globtopN.nbcontentzones, ncdf_globtopN.nbmonitzones, ncdf_globtopN.comp_globminc / 1000.0, ncerr);
			IfmInfo(pDoc, txtbuffer);
		}
		//Finally:
		isFirstExportedTS = false;
	}
	if (globvar_stopsimul) {
		IfmInfo(pDoc, "[CTRL-C @ PostTS] Breaking simulation at the end of this time step... (since globvar_stopsimul == true was DETECTED)");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_BREAK);
		Beep(440, 500);
	}
	else UpdateSimulRunProgressInfo(pDoc);
}

const double dtAfterYichg = 1e-6 * daysperyear;
const bool TSCactive = false; //BETA parameter for activating or not the OnTimeStepConstraint under-test callback

IfmBool Cpaleosea2d::OnTimeStepConstraint (IfmDocument pDoc, double tNow, double* dtProposed)
{
	double dtProposed0 = *dtProposed;
	bool proposed0isOk = true;
	if (TSCactive && hasclayYichanged && (dtProposed0 > dtAfterYichg)) {
		*dtProposed = dtAfterYichg;
		proposed0isOk = false;
		char txtbuffer[180];
		sprintf_s(txtbuffer, 180, "[OnTimeStepConstraint] dt for next time step was reduced to %.2e a, because clay-top has changed. (tNow=%.4e d; dt_pre=%.4e d; dt_post=%.4e d)", *dtProposed / daysperyear, tNow, dtProposed0, *dtProposed);
		IfmInfo(pDoc, txtbuffer);
	}
	return proposed0isOk;
}

void Cpaleosea2d::PreFlowSimulation (IfmDocument pDoc)
{
	if (false && needupdatearoundclaytop && currclayYi > 0) {
		//BETA STOP [TO BE REMOVED!]
		IfmInfo(pDoc, "calling: PreFlowSimulation");
		IfmWarning(pDoc, "Do you see the default h and c values?");
		IfmAlert(pDoc, NULL, "  OK  ", "BETA STOPPING @ PreFlowSimulation");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (withclaymesh && hasclayYichanged && hasclayYiincreased && firstflowsol) OverwriteAllNodesWithArrayValues(pDoc, true, false, true);
}

void Cpaleosea2d::PostFlowSimulation (IfmDocument pDoc)
{
	if (false && needupdatearoundclaytop && currclayYi > 0) {
		//BETA STOP [TO BE REMOVED!]
		IfmInfo(pDoc, "calling: PostFlowSimulation");
		IfmWarning(pDoc, "Do you see the default h and c values?");
		IfmAlert(pDoc, NULL, "  OK  ", "BETA STOPPING @ PostFlowSimulation");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (withclaymesh && hasclayYichanged && hasclayYiincreased && firstflowsol) OverwriteAllNodesWithArrayValues(pDoc, true, false, true);

	//NOTE: Darcy velocities extraction moved to @PostTimeStep

	firstflowsol = false;
}

void Cpaleosea2d::PreMassSimulation (IfmDocument pDoc, int iSpecies)
{
	if (false && needupdatearoundclaytop && currclayYi > 0) {
		//BETA STOP [TO BE REMOVED!]
		IfmInfo(pDoc, "in callback: PreMassSimulation (right after the overwritting task...)");
		IfmWarning(pDoc, "Do you see the default h and c values?");
		IfmAlert(pDoc, NULL, "  OK  ", "BETA STOPPING @ PreMassSimulation");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (withclaymesh && hasclayYichanged && hasclayYiincreased && firstmasssol) {
		OverwriteAllNodesWithArrayValues(pDoc, false, true, true, true);
		IfmWarning(pDoc, "[WCM @ PreMass] Please ignore the following convergence failure warning; it cannot be avoided when new elements are activated. (~MassSimulation)");
	}

	//BETA CODE TO REMOVE DEFINITIVELY SOON
	/* if (withclaymesh && hasclayYichanged && false) {
		double tempval = IfmGetResultsYVelocityValue(pDoc, 1799); //rock-top / clay-base: should always have defined velocities
		double tempval2 = IfmGetResultsYVelocityValue(pDoc, Pcurrclaytopnodes[0]);
		char txtbuffer[180];
		sprintf_s(txtbuffer, 180, "[PreMassSimulation TEMP] Y veloc values: at leftmost clay-top node q_Y = %g;  at chosen rock-top node q_Y = %g", tempval2, tempval);
		IfmWarning(pDoc, txtbuffer);
	} */
}

void Cpaleosea2d::PostMassSimulation (IfmDocument pDoc, int iSpecies)
{
	if (false && needupdatearoundclaytop && currclayYi > 0) {
		//BETA STOP [TO BE REMOVED!]
		IfmInfo(pDoc, "in callback: PostMassSimulation (right before overwritting mass conc. values again...)");
		IfmWarning(pDoc, "Do you see the default h and c values?");
		IfmAlert(pDoc, NULL, "  OK  ", "BETA STOPPING @ PostMassSimulation");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	}

	if (withclaymesh && hasclayYichanged && hasclayYiincreased && firstmasssol) OverwriteAllNodesWithArrayValues(pDoc, false, true, true, false);

	firstmasssol = false;
}

