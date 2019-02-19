// Analysis class of V0s
// OliverM 2019 Lund

#ifndef __MyAnalysisV0__
#define __MyAnalysisV0__

#include "TObject.h"
#include "MyAnalysis.h"

class AliAnalysisPIDEvent;	// forward declaration
class TClonesArray;
class TH1D;
class MyV0;

// Defining global constants here, should be changed.
enum { nBinsSpecies = 4, nBinsMulti = 3, nBinsSphero = 3};
const char* SPECIES[nBinsSpecies] = {"inc","K0s","L","Lbar"};
const char* MULTI[nBinsMulti] = {"MB","V0M","NCharged"};
const char* SPHERO[nBinsSphero] = {"MB","Jetty","Iso"};

class MyAnalysisV0: public MyAnalysis {

	public:
		//myAnalysisV0();
		MyAnalysisV0();	
		~MyAnalysisV0() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t CreateHistograms();
		Bool_t ProcessV0(MyV0 v0, Int_t Sp, Int_t Mu, Int_t Sph);

		ClassDef(MyAnalysisV0,1);

	protected:

		AliAnalysisPIDEvent* mEvent = 0;

		TClonesArray* bTracks = 0;
		TClonesArray* bV0s = 0;
		TClonesArray* bParticles = 0;

		Bool_t mFlagMC;

		// MONITORS
		TH1D* hEventMonitor;
		TH1D* hTrackMonitor;
		TH1D* hV0Monitor;
		TH1D* hParticleMonitor;

		// EVENT INFO HISTOGRAMS

		// TRACK HISTOGRAMS

		// MC PARTICLE HISTOGRAMS

		// V0 HISTOGRAMS
		TH1D* hV0Pt[nBinsSpecies][nBinsMulti][nBinsSphero];
};
#endif