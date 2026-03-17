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
class TH1F;
class TH2F;
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
		void ProcessRawYieldLossHist(TH2F* hist, TH1F* yieldhist, Int_t Sp, Double_t loose, Int_t opt);
		void DrawRawYieldLossHist(TH1F* da, TH1F* mc, Double_t ymax);
		void DrawVariation(Double_t cut, Int_t col, Int_t styl, TVirtualPad* can);

		void MakeEfficiencies();
		TH1F* ProcessEfficiency(TH2F* hist, TH1F* hmc, Int_t Sp);

		void MakeCorrectedYields();
		void MakeBarlowChecks();
		void MakeBarlowChecksPt();
		void BinHistogramsIntoRTrec();
		Int_t GetRTBinRec(const int& binRt, bool isLowEdge);

		void MakeDeviations();
		void CalculateAndPlotMaxDeviations();
		void CalculateMaxDeviationsNt();
		void AddDeviations();
		void AddDeviationsNt();
		TH1F* DivideAndComputeRogerBarlow(TH1F* h1, TH1F *h2);

		void CalculateSignalExSys();
		void CalculateSignalExSysNt();
		void CalculateFeeddownSys();
		void DrawMirrored(TH1F* hist, const char* opt = "");
		TH1F* CalculateRatiosRBSigmas(TH1F* ha, TH1F* hd);
		TH1F* CalculateSigExMaxD(TH1F* hm, TH1F* h1, TH1F* h2, TH1F* hd);

		ClassDef(MyAnalysisV0syst,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TList* mList;

		Double_t nRBsigmas = 1.;

		// V0 HISTOGRAMS
		//borrowed
		TH2F* hV0IMvRadiusL[V0consts::NSPECIES];
		TH2F* hV0IMvDCAdd[V0consts::NSPECIES];
		TH2F* hV0IMvCPA[V0consts::NSPECIES];
		TH2F* hV0IMvFastSignal[V0consts::NSPECIES];
		TH2F* hV0IMvCompMass[V0consts::NSPECIES];
		TH2F* hV0IMvLifetime[V0consts::NSPECIES];
		TH2F* hV0IMvNSigmaTPC[V0consts::NSPECIES];
		TH2F* hV0IMvDCAPVpos[V0consts::NSPECIES];
		TH2F* hV0IMvDCAPVneg[V0consts::NSPECIES];
		TH2F* hV0IMvNCluster[V0consts::NSPECIES];
		TH2F* hV0IMvNClusterF[V0consts::NSPECIES];

		TH2F* hMCV0IMvRadiusL[V0consts::NSPECIES];
		TH2F* hMCV0IMvDCAdd[V0consts::NSPECIES];
		TH2F* hMCV0IMvCPA[V0consts::NSPECIES];
		TH2F* hMCV0IMvFastSignal[V0consts::NSPECIES];
		TH2F* hMCV0IMvCompMass[V0consts::NSPECIES];
		TH2F* hMCV0IMvLifetime[V0consts::NSPECIES];
		TH2F* hMCV0IMvNSigmaTPC[V0consts::NSPECIES];
		TH2F* hMCV0IMvDCAPVpos[V0consts::NSPECIES];
		TH2F* hMCV0IMvDCAPVneg[V0consts::NSPECIES];
		TH2F* hMCV0IMvNCluster[V0consts::NSPECIES];
		TH2F* hMCV0IMvNClusterF[V0consts::NSPECIES];

		TH2F* hV0IMvPtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH2F* hMCV0IMvPtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1F* hMCV0Pt[V0consts::NSPECIES];

		TH3F* hV0IMvPtNtSys[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof][V0consts::sysVarsSizeof];

		TH1F* hV0PtFeeddown[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtFeeddownXi0[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtFeeddownXiErr[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];

		//owned
		TH1F* hV0YieldvRadiusL[V0consts::NSPECIES];
		TH1F* hV0YieldvDCAdd[V0consts::NSPECIES];
		TH1F* hV0YieldvCPA[V0consts::NSPECIES];
		TH1F* hV0YieldvFastSignal[V0consts::NSPECIES];
		TH1F* hV0YieldvCompMass[V0consts::NSPECIES];
		TH1F* hV0YieldvLifetime[V0consts::NSPECIES];
		TH1F* hV0YieldvNSigmaTPC[V0consts::NSPECIES];
		TH1F* hV0YieldvDCAPVpos[V0consts::NSPECIES];
		TH1F* hV0YieldvDCAPVneg[V0consts::NSPECIES];
		TH1F* hV0YieldvNCluster[V0consts::NSPECIES];
		TH1F* hV0YieldvNClusterF[V0consts::NSPECIES];
		
		TH1F* hMCV0YieldvRadiusL[V0consts::NSPECIES];
		TH1F* hMCV0YieldvDCAdd[V0consts::NSPECIES];
		TH1F* hMCV0YieldvCPA[V0consts::NSPECIES];
		TH1F* hMCV0YieldvFastSignal[V0consts::NSPECIES];
		TH1F* hMCV0YieldvCompMass[V0consts::NSPECIES];
		TH1F* hMCV0YieldvLifetime[V0consts::NSPECIES];
		TH1F* hMCV0YieldvNSigmaTPC[V0consts::NSPECIES];
		TH1F* hMCV0YieldvDCAPVpos[V0consts::NSPECIES];
		TH1F* hMCV0YieldvDCAPVneg[V0consts::NSPECIES];
		TH1F* hMCV0YieldvNCluster[V0consts::NSPECIES];
		TH1F* hMCV0YieldvNClusterF[V0consts::NSPECIES];

		TH1F* hV0PtSys[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1F* hV0PtSysRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1F* hV0PtSysMaxD[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof];
		TH1F* hV0PtSysMaxDRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO][V0consts::sysSizeof];
		TH1F* hV0EfficiencySys[V0consts::NSPECIES][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1F* hV0PtSysSum[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysSumUnc[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];

		// NT SYSTEMATICS
		TH2F* hV0PtNtSys[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH2F* hV0PtNtSysRatioToHM[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH2F* hV0PtNtSysMaxD[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof];
		TH2F* hV0PtNtSysSum[V0consts::NSPECIES][V0consts::NREGIONS];
		TH2F* hV0PtNtSysSumUnc[V0consts::NSPECIES][V0consts::NREGIONS];

		// RT SYSTEMATICS
		TH1F* hV0PtRtSys[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof][V0consts::sysVarsSizeof][V0consts::NRTBINS];
		TH1F* hV0PtRtSysRatioToHM[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof][V0consts::sysVarsSizeof][V0consts::NRTBINS];
		TH1F* hV0PtRtSysMaxDRatioToHM[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof][V0consts::NRTBINS];
		TH1F* hV0PtRtSysMaxD[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::sysSizeof][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSum[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSumUnc[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSigExLoose[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSigExDeflt[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSigExTight[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysMaxDSigEx[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSigExLooseRatioToHM[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSigExDefltRatioToHM[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysSigExTightRatioToHM[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		TH1F* hV0PtRtSysMaxDSigExRatioToHM[V0consts::NSPECIES][V0consts::NREGIONS][V0consts::NRTBINS];
		


		TH1F* hV0PtSysSigExLoose[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysSigExTight[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDSigEx[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysSigExLooseRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysSigExTightRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDSigExRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysBudget;
		TH1F* hV0PtSysEffi;
		TH1F* hV0PtSysExpBiasJetty;
		TH1F* hV0PtSysExpBiasIso;
		TH1F* hV0PtSysFeeddown[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDFeeddown[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysFeeddownXiErr[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDFeeddownXiErr[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDFeeddownTotal[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysFeeddownRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDFeeddownRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysFeeddownXiErrRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDFeeddownXiErrRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hV0PtSysMaxDFeeddownTotalRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];




		TH1F* hFracBudget[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracEffi[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracExpBias[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracCuts[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracSigEx[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracFD[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];

		TH1F* hFracEffiRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracExpBiasRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracCutsRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracSigExRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];
		TH1F* hFracFDRatioToHM[V0consts::NSPECIES][V0consts::NMULTI][V0consts::NSPHERO];


		TH1F* hRBcheck[V0consts::NSPECIES][V0consts::sysSizeof][V0consts::sysVarsSizeof];
		TH1F* hRBcheckPt[V0consts::NSPECIES][V0consts::NMULTI][V0consts::sysSizeof][V0consts::sysVarsSizeof];

		TF1* mParMuK0s = 0x0;
		TF1* mParSigK0s = 0x0;
		TF1* mParMuL = 0x0;
		TF1* mParSigL = 0x0;


};
#endif