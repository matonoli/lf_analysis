// Analysis class of V0s -- signal extraction sub-class
// OliverM 2019 Lund

#ifndef MYANALYSISV0EXTRACT_H
#define MYANALYSISV0EXTRACT_H

#include "TObject.h"
#include "TString.h"
#include "MyAnalysis.h"
#include "MyAnalysisV0.h"
#include "TF1.h"

class TFile;	// forward declaration
class TList;
class TH1F;
class TH2F;
class TNtuple;
class TTree;
class MyV0;
class MyEvent;
class MyTrack;
class MyParticle;

Bool_t gFReject;
Double_t gFLeft;
Double_t gFRight;
Double_t gfpol3(Double_t *x, Double_t *par)
{
    if (gFReject && x[0] > gFLeft && x[0] < gFRight) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

class MyAnalysisV0extract: public MyAnalysis {

	public:
		//MyAnalysisV0extract();
		MyAnalysisV0extract();	
		~MyAnalysisV0extract() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t BorrowHistograms();
		Bool_t CreateHistograms();
		Double_t* ExtractYieldSB(TH1F* hist = 0);
		Double_t* ExtractYieldSBVarySigma(Double_t nsig = 6., TH1F* hist = 0);
		Double_t* ExtractYieldFit(TH1F* hist = 0, Int_t Type = 0, Int_t MB = 0);
		Double_t* ExtractYieldFitRt(TTree* tree = 0, Int_t Type = 0);
		Double_t* ExtractYieldFitPtTree(TTree* tree = 0, Int_t Type = 0);
		void MakeExclusiveS0Bins();

		void DrawConstraints();
		void DefineSidebands();
		void TakeoverSidebands();
		void DoClosureTest(Int_t opt = 0);
		void DrawPad(Int_t Sp, Int_t Type);
		void ConvertNttoRt();

		void SetMCInputFile(const Char_t *name);
		void ProducePtSpectraFromHists();
		void ProducePtSpectraFromHistsRt();
		void ProducePtSpectraFromTrees();
		void ProduceRtSpectraFromTrees();

		ClassDef(MyAnalysisV0extract,1);

	protected:

		TString mOutName;
		TList* mList;
		TFile* mFileMC;

		Int_t iCan = 0;
		Int_t canCounter = 0;
		Int_t canCounterRt = 0;
		TCanvas* cFitsSB[V0consts::NSPECIES];
		TCanvas* cFits[(V0consts::NSPECIES-1) * V0consts::NMULTI * V0consts::NSPHERO * 2];
		TCanvas* cFitsRt[(V0consts::NSPECIES-1) * V0consts::NREGIONS * V0consts::NPTBINS * 60 * 2];
		TCanvas* cFitsPtTree[(V0consts::NSPECIES-1) * V0consts::NREGIONS * V0consts::NRTBINS0 * 2];

		Int_t nBins;
		const Double_t* xBins;
		Int_t nPads;
		Int_t nchmax;
		Int_t increm;

		//SIDEBANDS
		TH1F* hSidebandMean[V0consts::NSPECIES];
		TH1F* hSidebandSigma[V0consts::NSPECIES];
		TH1F* hSidebandSF[V0consts::NSPECIES];

		TF1* mParMuK0s = 0x0;
		TF1* mParSigK0s = 0x0;
		TF1* mParMuL = 0x0;
		TF1* mParSigL = 0x0;


		// V0 HISTOGRAMS
		//borrowed
		TH1F* hNchTrans;
		TH2F* hV0IMvPt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0Pt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		
		TNtuple* tV0massRCMB[V0consts::NSPECIES];
		TH2F* hV0IMvPtPrimary[V0consts::NSPECIES];
		TH2F* hV0IMvPtPrimaryPDG[V0consts::NSPECIES];
		TH2F* hV0IMvPtSecondary[V0consts::NSPECIES];
		TH2F* hV0IMvPtSecondaryPDG[V0consts::NSPECIES];
		TH2F* hV0IMvPtSecondaryXi[V0consts::NSPECIES];
		TH2F* hV0IMvPtBackground[V0consts::NSPECIES];

		TH3F* hV0IMvPtNt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0IMvPtRt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];


		

		//owned
		TH1F* hRtV0Yields[V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];

		TH1F* hV0PtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtRtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1F* hV0RtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];

		TH2F* hV0PtNtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMinFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];

		TH1F* hV0PtFitPrimary[V0consts::NSPECIES];
		TH1F* hV0PtFitPrimaryPDG[V0consts::NSPECIES];
		TH1F* hV0PtFitSecondary[V0consts::NSPECIES];
		TH1F* hV0PtFitSecondaryPDG[V0consts::NSPECIES];
		TH1F* hV0PtFitSecondaryXi[V0consts::NSPECIES];
		TH1F* hV0PtFitBackground[V0consts::NSPECIES];

		TH1F* hClosureTest[V0consts::NSPECIES];

		TH1F* hFitParam0[V0consts::NSPECIES];
		TH1F* hFitParam1[V0consts::NSPECIES];
		TH1F* hFitParam2[V0consts::NSPECIES];
		TH1F* hFitParam3[V0consts::NSPECIES];

		Double_t parRt0[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		Double_t parRt1[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		Double_t parRt2[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		Double_t parRt3[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];

		TF1* fPol1_0[V0consts::NSPECIES];
		TF1* fPol1_1[V0consts::NSPECIES];
		TF1* fPol1_2[V0consts::NSPECIES];
		TF1* fPol1_3[V0consts::NSPECIES];

};
#endif