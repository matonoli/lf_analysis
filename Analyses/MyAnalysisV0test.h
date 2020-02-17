// Analysis class of V0s -- final results calculating and testting sub-class
// OliverM 2019 Lund

#ifndef __MyAnalysisV0test__
#define __MyAnalysisV0test__

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


class MyAnalysisV0test: public MyAnalysis {

	public:
		//MyAnalysisV0test();
		MyAnalysisV0test();	
		~MyAnalysisV0test() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t CreateHistograms();
		
		ClassDef(MyAnalysisV0test,1);

	protected:

		TH1D* hTest;

};
#endif