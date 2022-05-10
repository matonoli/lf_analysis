// Analysis class of V0s
// OliverM 2019 Lund

#include "compInstructions.h"	// !THIS INCLUDES THE CUTS NAMESPACE!

#ifndef MYANALYSISV0_H
#define MYANALYSISV0_H

#include "TObject.h"
#include "TString.h"
#include "MyAnalysis.h"
#include "TH2F.h"


class TFile;	// forward declaration
class TList;
class TH1F;
class TH2F;
class TH3F;
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
	const Int_t NSPHERO = 1;//9;
	const Int_t NREGIONS = 5;
	const char* SPECIES[] = {"inc","K0s","L","Lbar"};
	const Float_t MASSES[] = {0., 0.497614, 1.11568, 1.11568};
	const Float_t XIMASS = 1.32171;
	const char* TYPE[NTYPE] = {"D","RC","MC"};
	const char* MULTI[] = {"MB","V0M","NCharged","V0M01","NCharged01"};//,"RTTrans","RTNear","RTAway"};
	//const char* MULTI[] = {"MB","V0M","NCharged","V0M","NCharged"};//temporary
	const char* PLOTS_MULTI[] = {"MB","V0M 0-10%","CL1 0-10%", "V0M 0-1%","CL1 0-1%"};//, "R_{T} Trans.","R_{T} Near","R_{T} Away"};
	const char* SPHERO[] = {"MB","Jetty20","Iso20", "Jetty10","Iso10","Jetty5","Iso5","Jetty1","Iso1"};//,"0-1","1-2","2-3","3-4","4-5"};
	const char* PLOTS_SPHERO[] = {"MB","Jetty 20%","Iso 20%", "Jetty 10%","Iso 10%","Jetty 5%","Iso %5","Jetty %1","Iso %1"};//,"0-1","1-2","2-3","3-4","4-5";
	const char* REGIONS[NREGIONS] = {"Trans","Near","Away","TransMin","TransMax"};
	const char* PLOTS_REGIONS[NREGIONS] = {"Trans.","Near","Away","Trans.-min","Trans.-max"};
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

	/*const Int_t NPTBINS = 42;		//superset of official K0s and L V0M spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 
		0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 
		1.8, 1.9, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 2.9, 3.0, 3.3, 
		3.4, 3.6, 3.9, 4.0, 4.2, 4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 
		8.0, 10.0, 12.0 };*/

	/*const Int_t NPTBINS = 16;		// offi HM L spectra
	const Double_t XBINS[NPTBINS+1] = { 
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.5, 2.9, 3.4, 4.0, 5.0, 
		6.5, 8.0 };*/

	const Int_t NPTBINS = 16;
	const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 5.0, 
		6.5, 8.0 };

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
		0., 5., 0., 0.5, 1.5, 2.5, 5.	};

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

	const Float_t RT_DEN		= 7.225;
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
		Bool_t ProcessV0toHistRT(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Reg, Int_t Nt, Int_t NtMin);
		Bool_t ProcessV0toTree(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu);
		Bool_t ProcessTrack(MyTrack &t, Int_t Type, Int_t Mu, Int_t Sph);
		Int_t ClassifyEvent(MyEvent &event, Int_t ntracks);
		Bool_t IsCentral(MyEvent &ev, Int_t Mu);
		Bool_t IsV0(MyV0 &v0, Int_t Sp, Int_t Type);
		Bool_t IsV0VaryCut(MyV0 &v0, Int_t Sp, Int_t Type, Int_t VarCut, Float_t VarVal);
		Bool_t IsTrans(Double_t phi1, Double_t phiTrig);
		Int_t IsMinOrMax(Bool_t isA, Int_t reg);
		Int_t  WhatRegion(Double_t phi1, Double_t phiTrig);
		Int_t  WhatRegionSide(Double_t phi1, Double_t phiTrig);
		Bool_t SelectV0Daughter(MyTrack &tr);
		Bool_t SelectParticle(MyParticle &p);
		Bool_t SelectTrack(MyTrack &tr);
		void ProcessV0SystVar(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph);


		void DoEfficiency();
		TH2F* RebinTH2(TH2F* h);
		TH3F* RebinTH3(TH3F* h);

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
		TH1F* hEventMonitor;
		TH1F* hTrackMonitor;
		TH1F* hV0Monitor;
		TH1F* hParticleMonitor;

		// EVENT INFO HISTOGRAMS
		TH1F* hEventCuts;
		TH1F* hEventVz;
		TH1F* hEventV0MCentrality;
		TH1F* hEventRefMult;
		TH2F* hEventV0MCentvRefMult;

		TH1F* hEventType;
		TH1F* hEventSpherocityV0M;
		TH1F* hEventSpherocityNCharged;
		TH1F* hEventSpherocityV0M01;
		TH1F* hEventSpherocityNCharged01;
		TH2F* hEventTSMCvRC;
		TH2F* hEventTSNormMCvRC;
		TH2F* hEventTSMCvNorm;
		TH2F* hEventTSRCvNorm;
		TH2F* hEventMultvSpheroD;
		TH2F* hEventMultvSpheroMC;

		TH2F* hLeadPhivPt;
		TH1F* hNchvLeadPt;
		TH2F* hNchvLeadPt2;
		TH2F* hNchMinvLeadPt2;
		TH2F* hNchMaxvLeadPt2;
		TH2F* hMeanPtvLeadPt2;
		TH2F* hMeanPtMinvLeadPt2;
		TH2F* hMeanPtMaxvLeadPt2;
		TH1F* hNchTrans;
		TH1F* hNchTrans2011;
		TH1F* hNchTransHybrid;
		TH1F* hNchTransNo2011Hybrid;
		TH1F* hNchTrans2011OrHybrid;
		TH1F* hTrackTransPhi2011;
		TH1F* hTrackTransPhiHybrid;
		TH1F* hTrackTransPhiNo2011Hybrid;
		TH1F* hTrackTransPhi2011OrHybrid;
		TH1F* hNchTransMC;
		TH2F* hNchTransRCvMC;
		TH2F* hNtvNtMin;
		TH2F* hNtvNtMax;
		TH2F* hNtMaxvNtMin;
		TH1F* hRt;
		TH1F* hRtMC;
		TH2F* hRtRCvMC;
		TH1F* hRt2;
		TH1F* hRt2MC;
		TH2F* hRt2RCvMC;
		TH2F* hLeadPtvNchTrans0;
		TH2F* hLeadPtvNchTrans;
		TH2F* hLeadPtvNchTransMin;
		TH2F* hLeadPtvNchTransMax;
		TH2F* hNchTransvSpherocityV0M;
		TH2F* hNchTransvSpherocityNCharged;
		TH2F* hV0DPhivNchTrans;
		TH2F* hTrackDPhivNchTrans;
		TH2F* hTrackDPhivNchTransMin;
		TH2F* hParticleDPhivNchTrans;
		TH2F* hK0sDPhivNchTrans;
		TH2F* hLDPhivNchTrans;
		TH2F* hLbarDPhivNchTrans;
		TH2F* hK0sDPhivNchTransMC;
		TH2F* hLDPhivNchTransMC;
		TH2F* hLbarDPhivNchTransMC;

		// TRACK HISTOGRAMS
		TH1F* hTrackPt[V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2F* hTrackEtavPhi[V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

		// MC PARTICLE HISTOGRAMS
		TH1F* hV0Efficiency[V0consts::NSPECIES];
		TH1F* hV0EfficiencyEta[V0consts::NSPECIES];
		TH1F* hV0EfficiencyRt[V0consts::NSPECIES][V0consts::NREGIONS];
		TH1F* hV0Feeddown[V0consts::NSPECIES];
		TH1F* hV0FeeddownPDG[V0consts::NSPECIES];
		TH3F* hV0FeeddownMatrix[V0consts::NSPECIES];
		TH1F* hV0FeeddownMotherPt[V0consts::NSPECIES];
		TH3F* hV0FeeddownMatrixXi0[V0consts::NSPECIES];
		TH1F* hV0FeeddownMotherPtXi0[V0consts::NSPECIES];

		TH2F* hParticlePrimaryvPDG;
		TH2F* hProtonNchTransvPt[V0consts::NREGIONS];
		TH2F* hPionNchTransvPt[V0consts::NREGIONS];
		TH2F* hLambdaNchTransvPt[V0consts::NREGIONS];
		TH2F* hK0sNchTransvPt[V0consts::NREGIONS];

		// V0 HISTOGRAMS
		TH1F* hV0Radius;
		TH1F* hV0ProperT;
		TH1F* hV0Pt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0Eta[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0Y[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2F* hV0EtavY[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2F* hV0IMvPt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2F* hV0DpiNsigTPCvpt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2F* hV0DprNsigTPCvpt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2F* hV0KFIMvIM[V0consts::NSPECIES][V0consts::NTYPE];
		TH2F* hV0DeltaIMvPt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2F* hV0BaselineIMvPt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2F* hV0FastSignalIMvPt[V0consts::NSPECIES][V0consts::NTYPE];
		TH2F* hV0CutIMvPt[V0consts::NSPECIES][V0consts::NTYPE][25];

		// V0 RT HISTOGRAMS
		TH2F* hV0PtNt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMin[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0EtaNt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PhiNt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH3F* hV0IMvPtNt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH3F* hV0IMvPtNtMin[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];

		TH1F* hV0PtNoTrigger[V0consts::NSPECIES];

		// DETAILED MC CLOSURE STUDY
		TH2F* hV0IMvPtPrimary[V0consts::NSPECIES];
		TH2F* hV0IMvPtPrimaryPDG[V0consts::NSPECIES];
		TH2F* hV0IMvPtSecondary[V0consts::NSPECIES];
		TH2F* hV0IMvPtSecondaryPDG[V0consts::NSPECIES];
		TH2F* hV0IMvPtSecondaryXi[V0consts::NSPECIES];
		TH2F* hV0IMvPtBackground[V0consts::NSPECIES];
		TH2F* hK0sLowDaughtersPDG;
		TH2F* hK0sLowIMvPDG;
		TH2F* hK0sLowIMvDaughterPDGPos;
		TH2F* hK0sLowIMvDaughterPDGNeg;
		TH2F* hK0sHighDaughtersPDG;
		TH2F* hK0sHighIMvPDG;
		TH2F* hK0sHighIMvDaughterPDGPos;
		TH2F* hK0sHighIMvDaughterPDGNeg;

		// V0 RC V MC HISTOGRAMS
		TH2F* hV0PtRCvMC[V0consts::NSPECIES];
		TH2F* hV0EtaRCvMC[V0consts::NSPECIES];
		TH2F* hV0PhiRCvMC[V0consts::NSPECIES];

		// V0 NTUPLES
		TNtuple* tV0PtMCMB[V0consts::NSPECIES];
		TNtuple* tV0massRCMB[V0consts::NSPECIES];

		// SYSTEMATICS STUDY
		TH2F* hV0IMvRadiusL[V0consts::NSPECIES];
		TH2F* hV0IMvDCAdd[V0consts::NSPECIES];
		TH2F* hV0IMvCPA[V0consts::NSPECIES];
		TH2F* hV0IMvFastSignal[V0consts::NSPECIES];
		TH2F* hV0IMvCompMass[V0consts::NSPECIES];
		TH2F* hV0IMvLifetime[V0consts::NSPECIES];
		TH2F* hV0IMvNSigmaTPC[V0consts::NSPECIES];
		TH2F* hV0IMvDCAPVpos[V0consts::NSPECIES];
		TH2F* hV0IMvDCAPVneg[V0consts::NSPECIES];
		TH2F* hV0IMvNCluster[V0consts::NSPECIES];
		TH2F* hV0IMvNClusterF[V0consts::NSPECIES];

		TH2F* hV0IMvPtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];

};
#endif