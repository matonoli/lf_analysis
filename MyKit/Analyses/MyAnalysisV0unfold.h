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
		Bool_t SumLambdas();
		void SetMCInputFile(const Char_t *name);

		void BinHistogramsIntoRT();
		void ApplySystematics();

		void DoUnfoldingNt();
		void DoUnfolding1D();
		void DoUnfoldingNtMin();
		void DoUnfolding1DMin();
		void DoUnfoldingNtMax();
		void DoUnfolding1DMax();

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
		TH2F* hV0PtNtMax[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMinFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMaxFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH1F* hNchTransMCTrigMC;

		TH1F* hV0PtRtSysSum[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSumUnc[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];

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

		TH1F* hNtMax;
		TH1F* hNtMaxRec;
		TH1F* hNtMaxGen;
		TH2F* hNtMaxRM;
		TH1F* hNtMaxUnf;
		TH1F* hNtMaxClosure;
		TH1F* hRtMaxUnf;

		TH2F* hV0PtNtFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMinFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		TH2F* hV0PtNtMaxFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];
		
		TH1F* hV0PtRtFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtMinFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtMaxFitCorrUnf[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];

		TH1F* hV0PtRtFitCorrUnfSyst[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtMinFitCorrUnfSyst[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtMaxFitCorrUnfSyst[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];

		TH1F* hV0PtRtFitCorrUnfSystUnc[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtMinFitCorrUnfSystUnc[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtMaxFitCorrUnfSystUnc[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS];

};
#endif