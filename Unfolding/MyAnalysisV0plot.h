// Analysis class of V0s -- final results calculating and plotting sub-class
// OliverM 2019 Lund

#ifndef MYANALYSISV0PLOT_H
#define MYANALYSISV0PLOT_H

#include "TObject.h"
#include "TString.h"
#include "MyAnalysis.h"
#include "MyAnalysisV0.h"

class TFile;	// forward declaration
class TList;
class TH1F;
class TH2F;
class MyV0;
class MyEvent;
class MyTrack;
class MyParticle;


class MyAnalysisV0plot: public MyAnalysis {

	public:
		//MyAnalysisV0plot();
		MyAnalysisV0plot();	
		~MyAnalysisV0plot() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t BorrowHistograms();
		Bool_t CreateHistograms();
		Bool_t CloneHistograms();

		void MakeFinalFiguresSpherocity();
		void MakeFinalFiguresRt();
		void MakeFinalFiguresEvent();
		void MakeFinalFiguresRatios();

		ClassDef(MyAnalysisV0plot,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TList* mList;

		// V0 HISTOGRAMS
		//borrowed
		TH1F* hEventSpherocityV0M;
		TH1F* hEventSpherocityNCharged;
		TH1F* hEventSpherocityV0M01;
		TH1F* hEventSpherocityNCharged01;
		TH1F* hNchTrans;
		TH1F* hRt2;
		TH2F* hNchvLeadPt2;
		TH2F* hTrackDPhivNchTrans;
		TH2F* hV0DPhivNchTrans;

		TH2F* hNchTransRCvMC;
		TH2F* hLeadPtvNchTrans0;
		TH2F* hLeadPtvNchTrans;

		TH1F* hV0Efficiency[V0consts::NSPECIES];
		TH1F* hV0EfficiencyRt[V0consts::NSPECIES][V0consts::NREGIONS];

		TH1F* hV0PtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtRtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1F* hV0PtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0Pt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hTrackPt[V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

		
		TH2F* hProtonNchTransvPt[V0consts::NREGIONS];
		TH2F* hPionNchTransvPt[V0consts::NREGIONS];
		TH2F* hLambdaNchTransvPt[V0consts::NREGIONS];
		TH2F* hK0sNchTransvPt[V0consts::NREGIONS];

		TH1F* hV0RtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		TH1F* hV0RtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];

		//owned
		TH1F* hBtoM[V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hBtoMRt[V0consts::NREGIONS+1][V0consts::NRTBINS0];
		TH1F* hV0toNchDR[2][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

};
#endif