// Analysis class of V0s
// OliverM 2019 Lund

#include "compInstructions.h"	// !THIS INCLUDES THE CUTS NAMESPACE!

#ifndef MYANALYSISV0_H
#define MYANALYSISV0_H

#include "TObject.h"
#include "TString.h"
#include "MyAnalysis.h"
#include "TH2D.h"


class TFile;	// forward declaration
class TList;
class TH1D;
class TH2D;
class TH3D;
class TNtuple;
class MyV0;
class MyEvent;
class MyTrack;
class MyParticle;
class TransverseSpherocity;

// Defining a namespace with constants
namespace V0consts {

	const bool MAKE_EXCLUSIVE = false; // RUN TRUE ONLY IF NSPHERO==9

	const Int_t NSPECIES = 4;//4;
	const Int_t NTYPE = 3; 
	const Int_t NMULTI = 5;//5;
	const Int_t NSPHERO = 9;//9;
	const Int_t NREGIONS = 3;
	const char* SPECIES[] = {"inc","K0s","L","Lbar"};
	const Float_t MASSES[] = {0., 0.497614, 1.11568, 1.11568};
	const Float_t XIMASS = 1.32171;
	const char* TYPE[NTYPE] = {"D","RC","MC"};
	const char* MULTI[] = {"MB","V0M","NCharged","V0M01","NCharged01"};//,"RTTrans","RTNear","RTAway"};
	//const char* MULTI[] = {"MB","V0M","NCharged","V0M","NCharged"};//temporary
	const char* PLOTS_MULTI[] = {"MB","V0M 0-10%","CL1 0-10%", "V0M 0-1%","CL1 0-1%"};//, "R_{T} Trans.","R_{T} Near","R_{T} Away"};
	const char* SPHERO[] = {"MB","Jetty20","Iso20", "Jetty10","Iso10","Jetty5","Iso5","Jetty1","Iso1"};//,"0-1","1-2","2-3","3-4","4-5"};
	const char* PLOTS_SPHERO[] = {"MB","Jetty 20%","Iso 20%", "Jetty 10%","Iso 10%","Jetty 5%","Iso %5","Jetty %1","Iso %1"};//,"0-1","1-2","2-3","3-4","4-5";
	const char* REGIONS[NREGIONS] = {"Trans","Near","Away"};
	const char* PLOTS_REGIONS[NREGIONS] = {"Trans.","Near","Away"};
	//const Int_t NPTBINS = 35;
	//const Double_t XBINS[NPTBINS+1] = { 0.00, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 
	/*const Int_t NPTBINS = 55;
	const Double_t XBINS[NPTBINS+1] = { 
		0.00, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 
		0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
    	0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 
		1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 
		3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 
		9.00, 10.00, 11.00, 12.00, 13.00, 14.00 };*/

	/*const Int_t NPTBINS = 44;		//official MB spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 3.9, 4.2, 
		4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 7.5, 8.0, 8.5, 9.2, 
		10.0, 11.0, 12.0, 13.5, 15.0 };*/

	/*const Int_t NPTBINS = 51;		//Omar K+- spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 
		0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 
		0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 
		1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 
		2.80, 3.00, 3.20, 3.40, 3.60, 3.80,	4.00, 4.50, 5.00, 
		5.50, 6.00, 6.50, 7.00, 8.00, 10.0, 20.0 };*/

	/*const Int_t NPTBINS = 54;		//official K+- spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.00, 0.05, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 
		0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 
		0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 
		1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 
		2.80, 3.00, 3.20, 3.40, 3.60, 3.80,	4.00, 4.50, 5.00, 
		5.50, 6.00, 6.50, 7.00, 8.00, 10.0, 13.0, 20.0 };
	*/

	/*const Int_t NPTBINS = 38;		//official K0s V0M spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 
		0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 
		1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 
		3.6, 3.9, 4.2, 4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 
		8.0, 10.0, 12.0 };*/

	/*const Int_t NPTBINS = 36;		//official K0s V0M spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 
		0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 
		1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 
		3.6, 3.9, 4.2, 4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 
		8.0, 10.0, 12.0 };*/


	/*const Int_t NPTBINS = 7;		//sys stability check
	const Double_t XBINS[NPTBINS+1] = {
		0.0, 0.1, 0.8, 
		1.7, 
		3.3, 
		7.0, 
		12.0 };*/

	const Int_t NPTBINS = 42;		//superset of official K0s and L V0M spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 
		0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 
		1.8, 1.9, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 2.9, 3.0, 3.3, 
		3.4, 3.6, 3.9, 4.0, 4.2, 4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 
		8.0, 10.0, 12.0 };

	/*const Int_t NPTBINS = 16;		// offi HM L spectra
	const Double_t XBINS[NPTBINS+1] = { 
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.5, 2.9, 3.4, 4.0, 5.0, 
		6.5, 8.0 };*/

	/*const Int_t NPTBINS = 16;
	const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 5.0, 
		6.5, 8.0 };*/

	/*const Int_t NPTBINS = 13;
	const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions, 1gev cut
		1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 5.0, 
		6.5, 8.0 };*/

	/*const Int_t NPTBINS = 18;
	const Double_t XBINS[NPTBINS+1] = { 
		0.00, 0.15, 0.20, 0.30, 0.40, 
		0.60, 0.80, 1.00, 1.40, 1.80, 
		2.20, 2.80, 3.40, 4.00, 5.50, 
		7.00, 9.00,	11.00, 14.00 };*/

	/*const Int_t NPTBINS2 = 28;
	const Double_t XBINS2[NPTBINS2+1] = { 
		0.00, 0.10, 0.14, 0.18, 0.25, 0.35, 
		0.45, 0.55, 0.65, 0.75, 0.85,
    	0.95, 1.10, 1.30, 1.50, 1.70, 
		1.90, 2.20, 2.60, 3.00, 3.40, 
		3.80, 4.50, 5.50, 6.50, 8.00, 
		10.00, 12.00, 14.00 };*/

	/*const Int_t NPTBINS2 = 28;
	const Double_t XBINS2[NPTBINS2+1] = { 
		0.00, 0.12, 0.16, 0.20, 0.30, 0.40, 
		0.50, 0.60, 0.70, 0.80, 0.90,
    	1.00, 1.20, 1.40, 1.60, 1.80, 
		2.00, 2.40, 2.80, 3.20, 3.60, 
		4.00, 5.00, 6.00, 7.00, 8.00, 
		10.00, 12.00, 14.00 };*/

	/*const Int_t NPTBINS2 = 18;
	const Double_t XBINS2[NPTBINS2+1] = { 
		0.00, 0.10, 0.20, 0.30, 0.40, 
		0.60, 0.80, 1.00, 1.40, 1.80, 
		2.20, 2.80, 3.40, 4.00, 5.50, 
		7.00, 9.00,	11.00, 15.00 };*/

	const Int_t NPTBINS2 = 18;
	const Double_t XBINS2[NPTBINS2+1] = { 
		0.00, 0.10, 0.20, 0.30, 0.40, 
		0.60, 0.80, 1.00, 1.40, 1.80, 
		2.20, 2.80, 3.40, 4.00, 5.00, 
		7.00, 9.00,	11.00, 15.00 };


	/*const Int_t NPTBINS2 = 44;		//official MB spectra
	const Double_t XBINS2[NPTBINS2+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 3.9, 4.2, 
		4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 7.5, 8.0, 8.5, 9.2, 
		10.0, 11.0, 12.0, 13.5, 15.0 };*/

	/*const Int_t NPTBINS2 = 55; 
	const Double_t XBINS2[NPTBINS2+1] = { 
		0.00, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 
		0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
    	0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 
		1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 
		3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 
		9.00, 10.00, 11.00, 12.00, 13.00, 14.00 };*/

	//const Int_t NRTBINS = 26;
	//const Double_t RTBINS[NRTBINS+1] = {
	//	-0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7,
	//	1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7,
	//	3.9, 4.1, 4.3, 4.5, 4.7, 4.9, 5.1	};

	const Int_t NRTBINS = 13;
	const Double_t RTBINS[NRTBINS+1] = {
		-0.1, 0.3, 0.7, 1.1, 1.5,
		1.9, 2.3, 2.7, 3.1, 3.5,
		3.9, 4.3, 4.7, 5.1	};


	/*const Int_t NRTBINS0 = 6;
	const Double_t RTBINS0[NRTBINS0+1] = {
		0., 5.0, 0., 0.5, 1., 1.5, 5.	};*/

	const Int_t NRTBINS0 = 6;
	const Double_t RTBINS0[NRTBINS0+1] = {
		0., 5.0, 0., 0.5, 1.5, 2.5, 5.	};

	/*const Int_t NRTBINS0 = 4;
	const Double_t RTBINS0[NRTBINS0+1] = {
		0., 5.0, 0., 1., 5.	};*/


	const Int_t NRTPTBINS = 3;
	const Double_t RT_PTRANGE[NRTPTBINS][2] = {
			{ 0.5 , 2.0 },
			{ 2.0 , 5.0 },
			{ 0.15 , 14.0 } };

	/*const Int_t NXIPTBINS = 7; 
  	const Double_t XIXBINS[NXIPTBINS+1] = 
  		{0.6, 1.2, 1.6, 2.2, 2.8, 3.6, 5.0, 6.5};*/

  	const Int_t NXIPTBINS = 45;		//official MB spectra
	const Double_t XIXBINS[NXIPTBINS+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 3.9, 4.2, 
		4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 7.5, 8.0, 8.5, 9.2, 
		10.0, 11.0, 12.0, 13.5, 15.0, 20.0 };

	/*const Int_t NXIPTBINS = 16;
	const Double_t XIXBINS[NXIPTBINS+1] = { 
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.5, 2.9, 3.4, 4.0, 5.0, 
		6.5, 8.0 };*/


	const Int_t PDG_IDS[] = {-999, 310, 3122, -3122};
	const Int_t PDG_IDS_DPOS[] = {-999, 211, 2212, 211};
	const Int_t PDG_IDS_DNEG[] = {-999, -211, -211, -2212};
	const char* SPECNAMES[] = {"inc.","K^{0}_{s}","#Lambda","#bar{#Lambda}"};
	//const Int_t COLOURS[7] = {kAzure-3,kOrange+8,kGreen+2,kMagenta+2, kViolet+10,kPink+10,kGreen+2};
	const Int_t COLOURS[7] = {kRed, kBlack, kGreen+2, kBlue, kViolet+10,kPink+10,kGreen+2};

	//const Int_t NEVENTTYPES = 34; //1+2+2+4+6+4+4 +2+4+4
	/*const char* EVENTTYPES[NEVENTTYPES] = {"MB pre-ES", "MB post-ES", "FHM", "MHM",
			"ISO", "JETTY", "FHM ISO", "FHM JETTY", "MHM ISO", "MHM JETTY",
			"RT", "RT int.", "RT 1-2", "RT 2-3", "RT 3-4", "RT 4-5",
			"FHM ISO MC", "FHM JETTY MC", "MHM ISO MC", "MHM JETTY MC" ,
			"MB ES rejected", "MB post-ES no V", "MB post-ES bad V", "MB post-ES good V" ,
			"FHM 0-1%", "MHM 0-1%",
			"FHM 0-1% ISO", "FHM 0-1% JETTY", "MHM 0-1% ISO", "MHM 0-1% JETTY",
			"FHM 0-1% ISO MC", "FHM 0-1% JETTY MC", "MHM 0-1% ISO MC", "MHM 0-1% JETTY MC"};*/

	const Int_t NEVENTTYPES = 10; //1+2+2+4+6+4+4 +2+4+4
	const char* EVENTTYPES[NEVENTTYPES] = {"MB pre-ES", "MB post-ES", 
	"MB ES rejected", "MB post-ES no V", "MB post-ES bad V", "MB post-ES good V" ,
	"FHM", "MHM", "FHM 0-1%", "MHM 0-1%"
	};

	const Float_t RT_DEN		= 7.449;
	const Float_t RT_DEN_MC		= 7.525;

	enum sysCuts { NoCut,														//0
		DCAdd, CPA, RadiusL, RadiusH,											//4
		FastSignal, CompMassK0s, CompMassL,						//8
		CompMassLbar, LifetimeK0s,					//12
		LifetimeL, LifetimeLbar, 												//14
		cutsV0cuts,																//15
		NSigmaTPCposPiL, NSigmaTPCposPiH, NSigmaTPCnegPiL, NSigmaTPCnegPiH,		//19
		NSigmaTPCposPrL, NSigmaTPCposPrH, NSigmaTPCnegPrL, NSigmaTPCnegPrH,		//23
		DCAPVpos, NClusterpos, NClusterFpos, 									//26
		DCAPVneg, NClusterneg, NClusterFneg,  									//29
		cutsSizeof 																//30
	};

	enum sysSources {
		sysRadiusL, sysDCAdd, sysCPA, sysFastSignal,
		sysCompMass, sysLifetime, sysNSigmaTPC, sysDCAPVpos,
		sysDCAPVneg, sysNCluster, sysNClusterF, sysSizeof
	};

	const char* SYSTS[sysSizeof] = { 
		"sysRadiusL", "sysDCAdd", "sysCPA", "sysFastSignal",
		"sysCompMass", "sysLifetime", "sysNSigmaTPC", "sysDCAPVpos",
		"sysDCAPVneg", "sysNCluster", "sysNClusterF" 
	};

	const char* PLOTS_SYSTS[sysSizeof] = { 
		"radius", "DCA_{d-d}", "cos PA", "fast signals",
		"comp. mass #sigma's", "lifetime", "TPC #sigma's", "DCA_{PV-d+}",
		"DCA_{PV-d-}", "TPC crossed rows", "TPC find. ratio" 
	};

	enum sysVars {
		loosest, loose, deflt, tight, tightest, sysVarsSizeof
	};

	const char* SYSTVAR[5] = {
		"loosest", "loose", "deflt", "tight", "tightest"
	};

	Double_t sysVar[NSPECIES-1][sysSizeof][sysVarsSizeof] = {
		{
			{0.3, 0.4, 0.5, 0.6, 0.7},
			{1.5, 1.25, 1.0, 0.75, 0.50},
			{0.95, 0.96, 0.97, 0.98, 0.99},
			{1., 1., 1., 2., 2.},
			{2.5, 3., 4., 5., 5.5},
			{0., 0., 0., 0., 0.},
			{6.5, 6., 5., 4., 3.5},
			{0.05, 0.055, 0.06, 0.07, 0.08},
			{0.05, 0.055, 0.06, 0.07, 0.08},
			{70., 70., 70., 75., 80.},
			{0.8, 0.8, 0.8, 0.95, 0.95}
		},
		{
			{0.3, 0.4, 0.5, 0.6, 0.7},
			{1.5, 1.25, 1.0, 0.75, 0.50},
			{0.993, 0.994, 0.995, 0.996, 0.997},
			{1., 1., 1., 2., 2.},
			{2.5, 3., 4., 5., 5.5},
			{35., 35., 30., 25., 25.},
			{6.5, 6., 5., 4., 3.5},
			{0.05, 0.055, 0.06, 0.07, 0.08},
			{0.05, 0.055, 0.06, 0.07, 0.08},
			{70., 70., 70., 75., 80.},
			{0.8, 0.8, 0.8, 0.95, 0.95}
		},
		{
			{0.3, 0.4, 0.5, 0.6, 0.7},
			{1.5, 1.25, 1.0, 0.75, 0.50},
			{0.993, 0.994, 0.995, 0.996, 0.997},
			{1., 1., 1., 2., 2.},
			{2.5, 3., 4., 5., 5.5},
			{35., 35., 30., 25., 25.},
			{6.5, 6., 5., 4., 3.5},
			{0.05, 0.055, 0.06, 0.07, 0.08},
			{0.05, 0.055, 0.06, 0.07, 0.08},
			{70., 70., 70., 75., 80.},
			{0.8, 0.8, 0.8, 0.95, 0.95}
		}
	};

}

class MyAnalysisV0: public MyAnalysis {

	public:
		//myAnalysisV0();
		MyAnalysisV0();	
		~MyAnalysisV0() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t CreateHistograms();
		Bool_t BorrowHistograms();

		Bool_t ProcessV0toHist(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph);
		Bool_t ProcessV0toTree(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu);
		Bool_t ProcessTrack(MyTrack &t, Int_t Type, Int_t Mu, Int_t Sph);
		Bool_t SelectEvent(MyEvent &ev, Int_t flag);
		Int_t ClassifyEvent(MyEvent &event, Int_t ntracks);
		Bool_t IsCentral(MyEvent &ev, Int_t Mu);
		Bool_t IsV0(MyV0 &v0, Int_t Sp, Int_t Type);
		Bool_t IsV0VaryCut(MyV0 &v0, Int_t Sp, Int_t Type, Int_t VarCut, Float_t VarVal);
		Bool_t IsTrans(Double_t phi1, Double_t phiTrig);
		Int_t  WhatRegion(Double_t phi1, Double_t phiTrig);
		Bool_t SelectV0Daughter(MyTrack &tr);
		Bool_t SelectParticle(MyParticle &p);
		Bool_t SelectTrack(MyTrack &tr);
		void ProcessV0SystVar(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph);


		void DoEfficiency();
		void DoEfficiencyFromTrees();
		void DoLambdaFeeddown();
		TH2D* RebinTH2(TH2D* h);

		ClassDef(MyAnalysisV0,2);

	protected:

		TList* mList;

		Bool_t mFlagMC;
		Bool_t mFlagHist;

		Double_t bugR;
		Double_t bugPt;
		TransverseSpherocity* mTS[V0consts::NTYPE][V0consts::NMULTI];
		TransverseSpherocity* mTSNorm[V0consts::NTYPE][V0consts::NMULTI];
		
		Double_t eventRt;
		Int_t nChTrans;
		Double_t phiLead;
		Double_t ptLead;

		TF1* mParMuK0s;
		TF1* mParSigK0s;
		TF1* mParMuL;
		TF1* mParSigL;

		Int_t mCutsDir[V0consts::NSPECIES][V0consts::NTYPE][V0consts::cutsSizeof+1];
		Double_t mCutsVal[V0consts::NSPECIES][V0consts::cutsSizeof+1];

		// MONITORS
		TH1D* hEventMonitor;
		TH1D* hTrackMonitor;
		TH1D* hV0Monitor;
		TH1D* hParticleMonitor;

		// EVENT INFO HISTOGRAMS
		TH1D* hEventCuts;
		TH1D* hEventVz;
		TH1D* hEventV0MCentrality;
		TH1D* hEventRefMult;
		TH2D* hEventV0MCentvRefMult;

		TH1D* hEventType;
		TH1D* hEventSpherocityV0M;
		TH1D* hEventSpherocityNCharged;
		TH1D* hEventSpherocityV0M01;
		TH1D* hEventSpherocityNCharged01;
		TH2D* hEventTSMCvRC;
		TH2D* hEventTSNormMCvRC;
		TH2D* hEventTSMCvNorm;
		TH2D* hEventTSRCvNorm;
		TH2D* hEventMultvSpheroD;
		TH2D* hEventMultvSpheroMC;

		TH2D* hLeadPhivPt;
		TH1D* hNchvLeadPt;
		TH2D* hNchvLeadPt2;
		TH1D* hNchTrans;
		TH1D* hNchTransMC;
		TH2D* hNchTransRCvMC;
		TH1D* hRt;
		TH1D* hRtMC;
		TH2D* hRtRCvMC;
		TH1D* hRt2;
		TH1D* hRt2MC;
		TH2D* hRt2RCvMC;
		TH2D* hLeadPtvNchTrans0;
		TH2D* hLeadPtvNchTrans;
		TH2D* hNchTransvSpherocityV0M;
		TH2D* hNchTransvSpherocityNCharged;
		TH2D* hTrackDPhivNchTrans;
		TH2D* hParticleDPhivNchTrans;
		TH2D* hV0DPhivNchTrans;

		// TRACK HISTOGRAMS
		TH1D* hTrackPt[V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2D* hTrackEtavPhi[V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

		// MC PARTICLE HISTOGRAMS
		TH1D* hV0Efficiency[V0consts::NSPECIES];
		TH1D* hV0EfficiencyEta[V0consts::NSPECIES];
		TH1D* hV0EfficiencyRt[V0consts::NSPECIES][V0consts::NREGIONS];
		TH1D* hV0Feeddown[V0consts::NSPECIES];
		TH1D* hV0FeeddownPDG[V0consts::NSPECIES];
		TH3D* hV0FeeddownMatrix[V0consts::NSPECIES];
		TH1D* hV0FeeddownMotherPt[V0consts::NSPECIES];
		TH3D* hV0FeeddownMatrixXi0[V0consts::NSPECIES];
		TH1D* hV0FeeddownMotherPtXi0[V0consts::NSPECIES];

		TH2D* hParticlePrimaryvPDG;
		TH2D* hProtonNchTransvPt[V0consts::NREGIONS];
		TH2D* hPionNchTransvPt[V0consts::NREGIONS];
		TH2D* hLambdaNchTransvPt[V0consts::NREGIONS];
		TH2D* hK0sNchTransvPt[V0consts::NREGIONS];

		// V0 HISTOGRAMS
		TH1D* hV0Radius;
		TH1D* hV0ProperT;
		TH1D* hV0Pt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0Eta[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0Y[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2D* hV0EtavY[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2D* hV0IMvPt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2D* hV0DpiNsigTPCvpt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2D* hV0DprNsigTPCvpt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2D* hV0KFIMvIM[V0consts::NSPECIES][V0consts::NTYPE];
		TH2D* hV0DeltaIMvPt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2D* hV0BaselineIMvPt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2D* hV0FastSignalIMvPt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2D* hV0CutIMvPt[V0consts::NSPECIES][V0consts::NTYPE][25];

		TH1D* hV0PtNoTrigger[V0consts::NSPECIES];

		// DETAILED MC CLOSURE STUDY
		TH2D* hV0IMvPtPrimary[V0consts::NSPECIES];
		TH2D* hV0IMvPtPrimaryPDG[V0consts::NSPECIES];
		TH2D* hV0IMvPtSecondary[V0consts::NSPECIES];
		TH2D* hV0IMvPtSecondaryPDG[V0consts::NSPECIES];
		TH2D* hV0IMvPtSecondaryXi[V0consts::NSPECIES];
		TH2D* hV0IMvPtBackground[V0consts::NSPECIES];
		TH2D* hK0sLowDaughtersPDG;
		TH2D* hK0sLowIMvPDG;
		TH2D* hK0sLowIMvDaughterPDGPos;
		TH2D* hK0sLowIMvDaughterPDGNeg;
		TH2D* hK0sHighDaughtersPDG;
		TH2D* hK0sHighIMvPDG;
		TH2D* hK0sHighIMvDaughterPDGPos;
		TH2D* hK0sHighIMvDaughterPDGNeg;

		// V0 RC V MC HISTOGRAMS
		TH2D* hV0PtRCvMC[V0consts::NSPECIES];
		TH2D* hV0EtaRCvMC[V0consts::NSPECIES];
		TH2D* hV0PhiRCvMC[V0consts::NSPECIES];

		// V0 NTUPLES
		TNtuple* tV0PtMCMB[V0consts::NSPECIES];
		TNtuple* tV0massRCMB[V0consts::NSPECIES];
		TNtuple* tV0PtMCRt[V0consts::NSPECIES][V0consts::NREGIONS];
		TNtuple* tV0massRt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];

		// SYSTEMATICS STUDY
		TH2D* hV0IMvRadiusL[V0consts::NSPECIES];
		TH2D* hV0IMvDCAdd[V0consts::NSPECIES];
		TH2D* hV0IMvCPA[V0consts::NSPECIES];
		TH2D* hV0IMvFastSignal[V0consts::NSPECIES];
		TH2D* hV0IMvCompMass[V0consts::NSPECIES];
		TH2D* hV0IMvLifetime[V0consts::NSPECIES];
		TH2D* hV0IMvNSigmaTPC[V0consts::NSPECIES];
		TH2D* hV0IMvDCAPVpos[V0consts::NSPECIES];
		TH2D* hV0IMvDCAPVneg[V0consts::NSPECIES];
		TH2D* hV0IMvNCluster[V0consts::NSPECIES];
		TH2D* hV0IMvNClusterF[V0consts::NSPECIES];

		TH2D* hV0IMvPtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];

};
#endif