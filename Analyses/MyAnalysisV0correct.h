// Analysis class of V0s -- efficiency correction sub-class
// OliverM 2019 Lund

#ifndef __MyAnalysisV0correct__
#define __MyAnalysisV0correct__

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


class MyAnalysisV0correct: public MyAnalysis {

	public:
		//MyAnalysisV0correct();
		MyAnalysisV0correct();	
		~MyAnalysisV0correct() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t BorrowHistograms();
		Bool_t CreateHistograms();
		Bool_t CloneHistograms();

		void SetMCInputFile(const Char_t *name);
		void NormaliseSpectra();
		void LoadEfficiency();
		void CorrectSpectra();

		ClassDef(MyAnalysisV0correct,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TList* mList;

		// V0 HISTOGRAMS
		//borrowed
		TH1D* hEventType;
		TH1D* hNchTrans;
		TH1D* hRt;
		TH1D* hRt2;
		TH1D* hRtV0Yields;
		TH1D* hV0PtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtRtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1D* hV0RtFit[V0consts::NSPECIES][1][1];
		TH1D* hV0Efficiency[V0consts::NSPECIES];

		//owned
		TH1D* hRtRebin;
		TH1D* hV0PtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtRtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1D* hV0RtFitCorr[V0consts::NSPECIES][1][1];

};
#endif