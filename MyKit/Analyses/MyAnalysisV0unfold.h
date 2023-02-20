// Analysis class of V0s -- unfolding sub-class
// OliverM 2019 Lund

#ifndef MYANALYSISV0UNFOLD_H
#define MYANALYSISV0UNFOLD_H

#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "MyAnalysis.h"
#include "MyAnalysisV0.h"
#include "TF1.h"

class TFile;	// forward declaration
class TList;
class TH1F;
class TH2F;
class TNTuple;
class MyEvent;
class UnfoldNTclass;


class MyAnalysisV0unfold: public MyAnalysis {

	public:
		//MyAnalysisV0unfold();
		MyAnalysisV0unfold();	
		~MyAnalysisV0unfold() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t BorrowHistograms();
		Bool_t NormaliseMC();
		Bool_t CreateHistograms();
		Bool_t CloneHistograms();
		void SetMCInputFile(const Char_t *name);

		void DoUnfoldingNt();
		void DoUnfolding1D();
		void DoUnfoldingNtMin();
		void DoUnfolding1DMin();

		TH2F* FlipMatrix(TH2F* h);
		void ComparisonPublished(TH1F* hRT);


		ClassDef(MyAnalysisV0unfold,1);

	protected:

		TString mOutName;
		TFile* mFileMC;
		TList* mList;

		UnfoldNTclass* mUnf;

		// HISTOGRAMS

		// BORROWED
		TH2F* hV0PtNt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMin[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMinFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH1F* hNchTransMCTrigMC;

		// OWNED
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

		TH2F* hV0PtNtFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMinFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		

};
#endif