// Analysis class of V0s -- final results calculating and plotting sub-class
// OliverM 2019 Lund

#ifndef __MyAnalysisV0plot__
#define __MyAnalysisV0plot__

#include "TObject.h"
#include "TString.h"
#include "../MyAnalysis.h"
#include "MyAnalysisV0.h"

class TFile;	// forward declaration
class TList;
class TH1D;
class TH2D;
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

		ClassDef(MyAnalysisV0plot,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TList* mList;

		// V0 HISTOGRAMS
		//borrowed
		TH1D* hEventSpherocityV0M;
		TH1D* hEventSpherocityNCharged;
		TH1D* hNchTrans;
		TH1D* hRt2;

		TH1D* hV0PtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtRtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1D* hV0PtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0Pt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hTrackPt[V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

		TH1D* hV0RtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		TH1D* hV0RtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];

		//owned
		TH1D* hBtoM[V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hBtoMRt[V0consts::NREGIONS+1][V0consts::NRTBINS0];
		TH1D* hV0toNchDR[2][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

};
#endif