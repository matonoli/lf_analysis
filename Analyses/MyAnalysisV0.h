// Analysis class of V0s
// OliverM 2019 Lund

#ifndef __MyAnalysisV0__
#define __MyAnalysisV0__

#include "TObject.h"
#include "TString.h"
#include "../MyAnalysis.h"
#include "cuts01.h"			// namespace cuts

class TFile;	// forward declaration
class TList;
class TH1D;
class TH2D;
class TNtuple;
class MyV0;
class MyEvent;
class MyTrack;
class MyParticle;
class TransverseSpherocity;

// Defining a namespace with constants
namespace V0consts {
	const Int_t NSPECIES = 4;
	const Int_t NTYPE = 3; 
	const Int_t NMULTI = 6;
	const Int_t NSPHERO = 8;
	const Int_t NREGIONS = 3;
	const char* SPECIES[NSPECIES] = {"inc","K0s","L","Lbar"};
	const Float_t MASSES[NSPECIES] = {0., 0.497614, 1.11568, 1.11568};
	const char* TYPE[NTYPE] = {"D","RC","MC"};
	const char* MULTI[NMULTI] = {"MB","V0M","NCharged","RTTrans","RTNear","RTAway"};
	const char* PLOTS_MULTI[NMULTI] = {"MB","V0M 0-10%","CL1", "R_{T} Trans.","R_{T} Near","R_{T} Away"};
	const char* SPHERO[NSPHERO] = {"MB","Jetty","Iso","0-1","1-2","2-3","3-4","4-5"};
	const char* REGIONS[NREGIONS] = {"Trans","Near","Away"};
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

	const Int_t NPTBINS = 44;		//official MB spectra
	const Double_t XBINS[NPTBINS+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 3.9, 4.2, 
		4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 7.5, 8.0, 8.5, 9.2, 
		10.0, 11.0, 12.0, 13.5, 15.0 };

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
		0.00, 0.15, 0.20, 0.30, 0.40, 
		0.60, 0.80, 1.00, 1.40, 1.80, 
		2.20, 2.80, 3.40, 4.00, 5.50, 
		7.00, 9.00,	11.00, 14.00 };*/

	const Int_t NPTBINS2 = 44;		//official MB spectra
	const Double_t XBINS2[NPTBINS2+1] = {
		0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 3.9, 4.2, 
		4.6, 5.0, 5.4, 5.9, 6.5, 7.0, 7.5, 8.0, 8.5, 9.2, 
		10.0, 11.0, 12.0, 13.5, 15.0 };

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


	const Int_t PDG_IDS[NSPECIES] = {-999, 310, 3122, -3122};
	const char* SPECNAMES[NSPECIES] = {"inc.","K^{0}_{s}","#Lambda","#bar{#Lambda}"};
	const Int_t COLOURS[7] = {kAzure-3,kOrange+8,kGreen+2,kMagenta+2, kViolet+10,kPink+10,kGreen+2};

	const Int_t NEVENTTYPES = 24; //1+2+2+4+6+4+4
	const char* EVENTTYPES[NEVENTTYPES] = {"MB pre-ES", "MB post-ES", "FHM", "MHM",
			"ISO", "JETTY", "FHM ISO", "FHM JETTY", "MHM ISO", "MHM JETTY",
			"RT", "RT 0-1", "RT 1-2", "RT 2-3", "RT 3-4", "RT 4-5",
			"FHM ISO MC", "FHM JETTY MC", "MHM ISO MC", "MHM JETTY MC" ,
			"MB ES rejected", "MB post-ES no V", "MB post-ES bad V", "MB post-ES good V"};

	const Float_t RT_DEN		= 7.449;
	const Float_t RT_DEN_MC		= 7.525; 

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
		Int_t ClassifyEvent(MyEvent &event, TClonesArray* trackArray = 0);
		Bool_t IsCentral(MyEvent &ev, Int_t Mu);
		Bool_t IsV0(MyV0 &v0, Int_t Sp, Int_t Type);
		Bool_t IsTrans(Double_t phi1, Double_t phiTrig);
		Int_t  WhatRegion(Double_t phi1, Double_t phiTrig);
		Bool_t SelectV0Daughter(MyTrack &tr);
		Bool_t SelectParticle(MyParticle &p);
		Bool_t SelectTrack(MyTrack &tr);


		void DoEfficiency();
		void DoEfficiencyFromTrees();
		void DoLambdaFeeddown();

		ClassDef(MyAnalysisV0,1);

	protected:

		TList* mList;

		Bool_t mFlagMC;
		Bool_t mFlagHist;

		Double_t bugR;
		Double_t bugPt;
		TransverseSpherocity* mTS[V0consts::NTYPE][V0consts::NMULTI-1];
		TransverseSpherocity* mTSNorm[V0consts::NTYPE][V0consts::NMULTI-1];
		
		Double_t eventRt;
		Int_t nChTrans;
		Double_t phiLead;
		Double_t ptLead;

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
		TH2D* hEventTSMCvRC;
		TH2D* hEventTSNormMCvRC;
		TH2D* hEventTSMCvNorm;
		TH2D* hEventTSRCvNorm;

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
		TH1D* hV0PtCut[25];

		// V0 NTUPLES
		TNtuple* tV0PtMCMB[V0consts::NSPECIES];
		TNtuple* tV0massRCMB[V0consts::NSPECIES];
		TNtuple* tV0PtMCRt[V0consts::NSPECIES][V0consts::NREGIONS];
		TNtuple* tV0massRt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
};
#endif