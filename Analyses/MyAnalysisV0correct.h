// Analysis class of V0s -- efficiency correction sub-class
// OliverM 2019 Lund

#ifndef __MyAnalysisV0correct__
#define __MyAnalysisV0correct__

#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "../MyAnalysis.h"
#include "MyAnalysisV0.h"

class TFile;	// forward declaration
class TList;
class TH1D;
class TH2D;
class TNTuple;
class MyV0;
class MyEvent;
class MyTrack;
class MyParticle;

Double_t rap_correction(Double_t* x, Double_t* par)
{
  Double_t pt = x[0];
  Double_t eta  = par[0];
  Double_t mass = par[1];
  const Double_t mt = TMath::Sqrt(pt*pt + mass*mass);
  const Double_t rap = TMath::ASinH(pt/mt*TMath::SinH(eta));
  //  return rap/eta;
  return rap/eta;
}

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
		void DoEfficiencyFromFile();
		void DoEfficiencyFromTrees();
		void CorrectSpectra();
		void StudyCuts();

		//Double_t rap_correction(Double_t* x, Double_t* par);

		ClassDef(MyAnalysisV0correct,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TList* mList;

		Double_t NormEta = 0;

		// V0 HISTOGRAMS
		//borrowed
		TH1D* hEventType;
		TH1D* hNchTrans;
		TH1D* hRt;
		TH1D* hRt2;
		TH1D* hRtV0Yields[V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		TH1D* hV0PtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtRtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1D* hV0RtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		TH1D* hV0Efficiency[V0consts::NSPECIES];
		TH1D* hV0EfficiencyRt[V0consts::NSPECIES][V0consts::NREGIONS];
		TH1D* hV0PtCut[25];
		TH1D* hV0PtCutMC[25];

		TNtuple* tV0PtMCMB[V0consts::NSPECIES];
		TNtuple* tV0massRCMB[V0consts::NSPECIES];
		TNtuple* tV0PtMCRt[V0consts::NSPECIES][V0consts::NREGIONS];
		TNtuple* tV0massRt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];

		//owned
		TH1D* hRtRebin;
		TH1D* hRt2Rebin;
		TH1D* hV0PtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtRtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1D* hV0RtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];

};
#endif