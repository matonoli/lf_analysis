// Analysis class of V0s -- efficiency correction sub-class
// OliverM 2019 Lund

#ifndef MYANALYSISV0SYST_H
#define MYANALYSISV0SYST_H

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


class MyAnalysisV0syst: public MyAnalysis {

	public:
		//MyAnalysisV0syst();
		MyAnalysisV0syst();	
		~MyAnalysisV0syst() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t BorrowHistograms();
		Bool_t CreateHistograms();
		Bool_t CloneHistograms();
		void SetMCInputFile(const Char_t *name);

		void StudyRawYieldLoss();
		void ProcessRawYieldLossHist(TH2D* hist, TH1D* yieldhist, Int_t Sp, Double_t loose, Int_t opt);
		void DrawRawYieldLossHist(TH1D* da, TH1D* mc, Double_t ymax);
		void DrawVariation(Double_t cut, Int_t col, Int_t styl, TVirtualPad* can);

		void MakeEfficiencies();
		TH1D* ProcessEfficiency(TH2D* hist, TH1D* hmc, Int_t Sp);

		void MakeCorrectedYields();
		void MakeBarlowChecks();
		void MakeBarlowChecksPt();
		void MakeDeviations();
		void AddDeviations();
		TH1D* DivideAndComputeRogerBarlow(TH1D* h1, TH1D *h2);

		void CalculateSignalExSys();
		void CalculateFeeddownSys();
		void DrawMirrored(TH1D* hist, const char* opt = "");

		ClassDef(MyAnalysisV0syst,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TList* mList;

		Double_t nRBsigmas = 1.;

		// V0 HISTOGRAMS
		//borrowed
		TH2D* hV0IMvRadiusL[V0consts::NSPECIES];
		TH2D* hV0IMvDCAdd[V0consts::NSPECIES];
		TH2D* hV0IMvCPA[V0consts::NSPECIES];
		TH2D* hV0IMvFastSignal[V0consts::NSPECIES];
		TH2D* hV0IMvCompMass[V0consts::NSPECIES];
		TH2D* hV0IMvLifetime[V0consts::NSPECIES];
		TH2D* hV0IMvNSigmaTPC[V0consts::NSPECIES];
		TH2D* hV0IMvDCAPVpos[V0consts::NSPECIES];
		TH2D* hV0IMvDCAPVneg[V0consts::NSPECIES];
		TH2D* hV0IMvNCluster[V0consts::NSPECIES];
		TH2D* hV0IMvNClusterF[V0consts::NSPECIES];

		TH2D* hMCV0IMvRadiusL[V0consts::NSPECIES];
		TH2D* hMCV0IMvDCAdd[V0consts::NSPECIES];
		TH2D* hMCV0IMvCPA[V0consts::NSPECIES];
		TH2D* hMCV0IMvFastSignal[V0consts::NSPECIES];
		TH2D* hMCV0IMvCompMass[V0consts::NSPECIES];
		TH2D* hMCV0IMvLifetime[V0consts::NSPECIES];
		TH2D* hMCV0IMvNSigmaTPC[V0consts::NSPECIES];
		TH2D* hMCV0IMvDCAPVpos[V0consts::NSPECIES];
		TH2D* hMCV0IMvDCAPVneg[V0consts::NSPECIES];
		TH2D* hMCV0IMvNCluster[V0consts::NSPECIES];
		TH2D* hMCV0IMvNClusterF[V0consts::NSPECIES];

		TH2D* hV0IMvPtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH2D* hMCV0IMvPtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1D* hMCV0Pt[V0consts::NSPECIES];

		TH1D* hV0PtFeeddown[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtFeeddownXi0[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtFeeddownXiErr[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];

		//owned
		TH1D* hV0YieldvRadiusL[V0consts::NSPECIES];
		TH1D* hV0YieldvDCAdd[V0consts::NSPECIES];
		TH1D* hV0YieldvCPA[V0consts::NSPECIES];
		TH1D* hV0YieldvFastSignal[V0consts::NSPECIES];
		TH1D* hV0YieldvCompMass[V0consts::NSPECIES];
		TH1D* hV0YieldvLifetime[V0consts::NSPECIES];
		TH1D* hV0YieldvNSigmaTPC[V0consts::NSPECIES];
		TH1D* hV0YieldvDCAPVpos[V0consts::NSPECIES];
		TH1D* hV0YieldvDCAPVneg[V0consts::NSPECIES];
		TH1D* hV0YieldvNCluster[V0consts::NSPECIES];
		TH1D* hV0YieldvNClusterF[V0consts::NSPECIES];
		
		TH1D* hMCV0YieldvRadiusL[V0consts::NSPECIES];
		TH1D* hMCV0YieldvDCAdd[V0consts::NSPECIES];
		TH1D* hMCV0YieldvCPA[V0consts::NSPECIES];
		TH1D* hMCV0YieldvFastSignal[V0consts::NSPECIES];
		TH1D* hMCV0YieldvCompMass[V0consts::NSPECIES];
		TH1D* hMCV0YieldvLifetime[V0consts::NSPECIES];
		TH1D* hMCV0YieldvNSigmaTPC[V0consts::NSPECIES];
		TH1D* hMCV0YieldvDCAPVpos[V0consts::NSPECIES];
		TH1D* hMCV0YieldvDCAPVneg[V0consts::NSPECIES];
		TH1D* hMCV0YieldvNCluster[V0consts::NSPECIES];
		TH1D* hMCV0YieldvNClusterF[V0consts::NSPECIES];

		TH1D* hV0PtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1D* hV0PtSysMaxD[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof];
		TH1D* hV0EfficiencySys[V0consts::NSPECIES][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1D* hV0PtSysSum[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysSumUnc[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];

		TH1D* hV0PtSysSigExLoose[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysSigExTight[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysMaxDSigEx[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysBudget;
		TH1D* hV0PtSysEffi;
		TH1D* hV0PtSysFeeddown[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysMaxDFeeddown[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysFeeddownXiErr[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysMaxDFeeddownXiErr[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0PtSysMaxDFeeddownTotal[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];

		TH1D* hFracBudget[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hFracEffi[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hFracCuts[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hFracSigEx[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hFracFD[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];


		TH1D* hRBcheck[V0consts::NSPECIES][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1D* hRBcheckPt[V0consts::NSPECIES][V0consts::NMULTI][V0consts::sysSizeof][V0consts::sysVarsSizeof];

		TF1* mParMuK0s = 0x0;
		TF1* mParSigK0s = 0x0;
		TF1* mParMuL = 0x0;
		TF1* mParSigL = 0x0;


};
#endif