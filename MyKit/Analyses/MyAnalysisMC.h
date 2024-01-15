// Analysis class of V0s
// OliverM 2019 Lund

#include "compInstructions.h"	// !THIS INCLUDES THE CUTS NAMESPACE!

#ifndef MYANALYSISMC_H
#define MYANALYSISMC_H

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
class UnfoldNTclass;

namespace MCconsts {
		
	const Int_t NTYPE = 3; 
	const char* TYPE[NTYPE] = {"RC","MC","D"};
			
	const Int_t NREGIONS = 5;
	const char* REGIONS[] = {"Trans","Near","Away","TransMin","TransMax"};
	const char* PLOTS_REGIONS[] = {"Trans.","Near","Away","Trans.-min","Trans.-max"};

	const Int_t NSPECIES = 6;
	const char* SPECIES[] = {"piKp","pr","XiInc","Xi","phiInc","phi"};
			
	const int NRTBINS = 5;
	const double RTBINS[NRTBINS+1] = {5.0,0.0, 0.85, 1.5, 2.5, 5.0};
			
	const Float_t RT_DEN		= 7.225;
	const Float_t RT_DEN_MC		= 7.525;
		
	// PT BINS
	const Int_t NPIPTBINS = 16;
	const Double_t PIXBINS[NPIPTBINS+1] = { // divisible by omar's pions
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 5.0, 
		6.5, 8.0 };
		
	const Int_t NXIPTBINS = 7; 
  	const Double_t XIXBINS[NXIPTBINS+1] = {
  		0.6, 1.2, 1.6, 2.2, 
  		2.8, 3.6, 5.0, 6.5 };

  	const Int_t NPHIPTBINS = 14; 
  	const Double_t PHIXBINS[NPHIPTBINS+1] = {
  		0.5, 0.7, 0.9, 1.2, 
  		1.4, 1.6, 1.8, 2.0,
  		2.2, 2.6, 3.0, 3.5,
  		4.0, 5.0, 8.0 };

  	const Double_t NPTBINS[NSPECIES] = {
  		NPIPTBINS, NPIPTBINS, NXIPTBINS, NXIPTBINS, NPHIPTBINS, NPHIPTBINS
  	};

  	const Double_t* XBINS[NSPECIES] = {
  		PIXBINS, PIXBINS, XIXBINS, XIXBINS, PHIXBINS, PHIXBINS
  	};
		
  	const Int_t NEVENTTYPES = 10; //1+2+2+4+6+4+4 +2+4+4
	const char* EVENTTYPES[NEVENTTYPES] = {"MB pre-ES", "MB post-ES", 
		"MB ES rejected", "MB post-ES no V", "MB post-ES bad V", "MB post-ES good V" ,
		"FHM", "MHM", "FHM 0-1%", "MHM 0-1%"
	};
}

class MyAnalysisMC: public MyAnalysis {

	public:



		//myAnalysisV0();
		MyAnalysisMC();	
		~MyAnalysisMC() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t CreateHistograms();
		Bool_t BorrowHistograms();

		Int_t ClassifyEvent(MyEvent &event, Int_t ntracks);
		Bool_t IsTrans(Double_t phi1, Double_t phiTrig);
		Int_t IsMinOrMax(Bool_t isA, Int_t reg);
		Int_t  WhatRegion(Double_t phi1, Double_t phiTrig);
		Int_t  WhatRegionSide(Double_t phi1, Double_t phiTrig);
		
		Bool_t SelectParticle(MyParticle &p);
		Bool_t SelectTrack(MyTrack &tr);
		Bool_t SelectV0Daughter(MyTrack &tr);
		Bool_t IsGeometricalCut(Float_t phiprime, Float_t pt);
		Double_t FlipNegativeAngle(Double_t phi);
		TH2F* FlipMatrix(TH2F* h);
		TH2F* ScaleWidthTH2(TH2F* h);
		TH2F* DivideTH2ByTH1(TH2F* h, TH1D* d);
		TH2F* ScaleTH2Rows(TH2F* h, Double_t d);

		Bool_t CalculateEfficiencies();
		Bool_t Normalise();
		Bool_t CorrectEfficiency();
		Bool_t Unfold();
		void LoadDataXi();
		void BinHistogramsIntoRT();


		void DoUnfoldingNt();
		void DoUnfolding1D();
		void DoUnfoldingNtMin();
		void DoUnfolding1DMin();
		void DoUnfoldingNtMax();
		void DoUnfolding1DMax();
		
		/*Bool_t ProcessV0toHist(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph);
		Bool_t ProcessV0toHistRT(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Reg, Int_t Nt, Int_t NtMin, Int_t NtMax);
		Bool_t ProcessV0toTree(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu);
		Bool_t ProcessTrack(MyTrack &t, Int_t Type, Int_t Mu, Int_t Sph);
*/

		ClassDef(MyAnalysisMC,2);

	protected:

		// NECESSARY
		TList* mList;
		Bool_t mFlagMC;
		Bool_t mFlagHist;
		Double_t bugR;
		Double_t bugPt;
		UnfoldNTclass* mUnf;

		// DECLARE AUXILLIARY ENUMS
		enum { RC, MC, D, sizeofTypes };
		enum { piKp, pr, XiInc, Xi, phiInc, phi, sizeofSpecies };

		// IMPORTANT EVENT GLOBALS		
		Double_t eventRt;
		Int_t nChTrans;
		Double_t phiLead;
		Double_t phiPrimeLead;
		Double_t ptLead;

		// MONITORS
		TH1D* hEventMonitor;
		TH1D* hEventType;
		TH1D* hEventVz;

		// NT HISTOGRAMS
		TH2F* hLeadPhivPt;
		TH2F* hLeadPhiPrimevPt;
		TH1F* hLeadPDG;
		TH1F* hLeadMotherPDG;
		TH1F* hLeadDCAPri;
		TH1F* hLeadDCASec;

		TH1F* hNchvLeadPt;
		TH2F* hNchvLeadPt2;
		TH2F* hNchMinvLeadPt2;
		TH2F* hNchMaxvLeadPt2;
		TH1F* hNchTrans;
		TH2F* hNtvNtMin;
		TH2F* hNtvNtMax;
		TH2F* hNtMaxvNtMin;
		TH2F* hLeadPtvNchTrans;
		TH2F* hLeadPtvNchTrans0;
		TH2F* hLeadPtvNchTransMin;
		TH2F* hLeadPtvNchTransMax;
		TH1F* hNchTransRC;
		TH1F* hNchTransMC;
		TH2F* hNchTransRCvMC;
		TH1F* hNchTransMinRC;
		TH1F* hNchTransMinMC;
		TH2F* hNchTransMinRCvMC;
		TH1F* hNchTransMaxRC;
		TH1F* hNchTransMaxMC;
		TH2F* hNchTransMaxRCvMC;
		TH1F* hNchTransMCTrigMC;
		TH1F* hNchTransMinMCTrigMC;
		TH1F* hNchTransMaxMCTrigMC;
		TH2F* hMeanPtvLeadPt2;
		TH2F* hMeanPtMinvLeadPt2;
		TH2F* hMeanPtMaxvLeadPt2;
		TH2F* hPhiDaughterRegionsPt[MCconsts::NTYPE];
		TH2F* hPhiDaughterDPhiPt[MCconsts::NTYPE];
		TH1F* hXiBachDCAXY;
		TH1F* hXiPrDCAXY;

		//PID NT HISTOGRAMS
		TH1D* hPIDPt[MCconsts::NSPECIES][MCconsts::NTYPE];
		TH1D* hPIDEffi[MCconsts::NSPECIES];
		TH2F* hPIDDPhivNchTrans[MCconsts::NSPECIES];

		//PID NT HISTOGRAMS
		TH2F* hPIDPtNt[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		TH2F* hPIDPtNtMin[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		TH2F* hPIDPtNtMax[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		
		// POST PROCESSING HISTOGRAMS
		TH1D* hPIDPtCorr[MCconsts::NSPECIES][MCconsts::NTYPE];
		TH2F* hPIDPtNtCorr[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		TH2F* hPIDPtNtMinCorr[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		TH2F* hPIDPtNtMaxCorr[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		
		TH2F* hPIDPtNtCorrUnf[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		TH2F* hPIDPtNtMinCorrUnf[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];
		TH2F* hPIDPtNtMaxCorrUnf[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS];

		// UNFOLDING HISTOGRAMS
		TH1F* hNt;
		TH1F* hNtRec;
		TH1F* hNtGen;
		TH2F* hNtRM;
		TH1F* hNtUnf;
		TH1F* hNtClosure;
		TH1F* hRtUnf;

		TH1F* hNtMin;
		TH1F* hNtMinRec;
		TH1F* hNtMinGen;
		TH2F* hNtMinRM;
		TH1F* hNtMinUnf;
		TH1F* hNtMinClosure;
		TH1F* hRtMinUnf;

		TH1F* hNtMax;
		TH1F* hNtMaxRec;
		TH1F* hNtMaxGen;
		TH2F* hNtMaxRM;
		TH1F* hNtMaxUnf;
		TH1F* hNtMaxClosure;
		TH1F* hRtMaxUnf;

		// DATA
		TH1F* hPIDPtRtCorrUnf[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS][MCconsts::NRTBINS];
		TH1F* hPIDPtRtMinCorrUnf[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS][MCconsts::NRTBINS];
		TH1F* hPIDPtRtMaxCorrUnf[MCconsts::NSPECIES][MCconsts::NTYPE][MCconsts::NREGIONS][MCconsts::NRTBINS];



};
#endif