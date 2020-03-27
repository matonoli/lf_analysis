// Analysis class of V0s -- efficiency correction sub-class
// OliverM 2019 Lund

#ifndef MYANALYSISV0CORRECT_H
#define MYANALYSISV0CORRECT_H

#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "MyAnalysis.h"
#include "MyAnalysisV0.h"
#include "TF1.h"

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

Double_t
LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 *
LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e8);
  return fLevyTsallis;
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
		void SetXiSpectraFile(const Char_t *name);
		void NormaliseSpectra();
		void LoadEfficiency();
		void DoEfficiencyFromFile();
		void DoEfficiencyFromTrees();
		void CorrectForFeeddown();
		void CorrectSpectra();
		void StudyCuts();
		void DoClosureTest(Int_t opt = 0);
		void DoXCheckV0M();

		//Double_t rap_correction(Double_t* x, Double_t* par);

		ClassDef(MyAnalysisV0correct,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TFile* mFileXi = 0x0;
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
		TH1D* hV0Pt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

		TH1D* hClosureTestCorr[V0consts::NSPECIES];

		TNtuple* tV0PtMCMB[V0consts::NSPECIES];
		TNtuple* tV0massRCMB[V0consts::NSPECIES];
		TNtuple* tV0PtMCRt[V0consts::NSPECIES][V0consts::NREGIONS];
		TNtuple* tV0massRt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS];

		TH2D* hV0FeeddownMatrix[V0consts::NSPECIES];
		TH1D* hV0FeeddownMotherPt[V0consts::NSPECIES];

		//owned
		TH1D* hRtRebin;
		TH1D* hRt2Rebin;
		TH1D* hV0PtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtRtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTBINS0];
		TH1D* hV0RtFitCorr[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NREGIONS][V0consts::NRTPTBINS];
		TH1D* hV0PtFeeddown[V0consts::NSPECIES][V0consts::NMULTI];
		TH1D* hXiPt[V0consts::NSPECIES][V0consts::NMULTI];

};
#endif