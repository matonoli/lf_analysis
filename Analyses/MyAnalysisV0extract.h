// Analysis class of V0s -- signal extraction sub-class
// OliverM 2019 Lund

#ifndef __MyAnalysisV0extract__
#define __MyAnalysisV0extract__

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


class MyAnalysisV0extract: public MyAnalysis {

	public:
		//MyAnalysisV0extract();
		MyAnalysisV0extract();	
		~MyAnalysisV0extract() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t BorrowHistograms();
		Bool_t CreateHistograms();
		Double_t* ExtractYieldFit(TH1D* hist = 0);

		ClassDef(MyAnalysisV0extract,1);

	protected:

		TString mOutName;
		TList* mList;

		Int_t iCan = 0;
		Int_t canCounter = 0;
		TCanvas* cFits[(V0consts::NSPECIES-1) * V0consts::NMULTI * V0consts::NSPHERO * 2];


		// V0 HISTOGRAMS
		//borrowed
		TH2D* hV0IMvPt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

		//owned
		TH1D* hV0PtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

};
#endif