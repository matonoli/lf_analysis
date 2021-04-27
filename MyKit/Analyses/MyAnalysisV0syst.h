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
		void DrawVariation(Double_t cut, Int_t col, TVirtualPad* can);

		ClassDef(MyAnalysisV0syst,1);

	protected:

		TString mOutName;
		//TFile* mFileOut;
		TFile* mFileMC;
		TList* mList;

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
		

};
#endif