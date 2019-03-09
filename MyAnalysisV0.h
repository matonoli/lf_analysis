// Analysis class of V0s
// OliverM 2019 Lund

#ifndef __MyAnalysisV0__
#define __MyAnalysisV0__

#include "TObject.h"
#include "MyAnalysis.h"
#include "cuts00.h"			// namespace cuts

class AliAnalysisPIDEvent;	// forward declaration
class TClonesArray;
class TH1D;
class TH2D;
class MyV0;
class MyEvent;
class MyTrack;
class MyParticle;
class TransverseSpherocity;

// Defining a namespace with constants
namespace V0consts {
	const Int_t NSPECIES = 4;
	const Int_t NTYPE = 3; 
	const Int_t NMULTI = 3;
	const Int_t NSPHERO = 3;
	const char* SPECIES[NSPECIES] = {"inc","K0s","L","Lbar"};
	const char* TYPE[NTYPE] = {"D","RC","MC"};
	const char* MULTI[NMULTI] = {"MB","V0M","NCharged"};
	const char* SPHERO[NSPHERO] = {"MB","Jetty","Iso"};
	const Int_t NPTBINS = 35;
	const Double_t XBINS[NPTBINS+1] = { 0.00, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 
		1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 
		3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 
		9.00, 10.00, 11.00, 12.00, 13.00, 14.00 };

	const Int_t PDG_IDS[NSPECIES] = {-999, 310, 3122, -3122};

}

class MyAnalysisV0: public MyAnalysis {

	public:
		//myAnalysisV0();
		MyAnalysisV0();	
		~MyAnalysisV0() { }
		Int_t Init();
		Int_t Make(Int_t iEv);
		Int_t Finish();
		Bool_t CreateHistograms();
		Bool_t ProcessV0(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph);
		Bool_t SelectEvent(MyEvent &ev);
		Bool_t IsCentral(MyEvent &ev, Int_t Mu);
		Bool_t IsV0(MyV0 &v0, Int_t Sp, Int_t Type);
		Bool_t SelectV0Daughter(MyTrack &tr);
		Bool_t SelectParticle(MyParticle &p);
		Bool_t SelectTrack(MyTrack &tr);
		Double_t* ExtractYieldFit(TH1D* hist = 0);

		ClassDef(MyAnalysisV0,1);

	protected:

		AliAnalysisPIDEvent* mEvent = 0;

		TClonesArray* bTracks = 0;
		TClonesArray* bV0s = 0;
		TClonesArray* bParticles = 0;

		Bool_t mFlagMC;
		Double_t bugR;
		Double_t bugPt;
		TransverseSpherocity* mTS;

		// MONITORS
		TH1D* hEventMonitor;
		TH1D* hTrackMonitor;
		TH1D* hV0Monitor;
		TH1D* hParticleMonitor;

		// EVENT INFO HISTOGRAMS
		TH1D* hEventV0MCentrality;
		TH1D* hEventRefMult;
		TH2D* hEventV0MCentvRefMult;

		TH1D* hEventSpherocity;

		// TRACK HISTOGRAMS

		// MC PARTICLE HISTOGRAMS

		// V0 HISTOGRAMS
		TH1D* hV0Pt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH1D* hV0Eta[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
		TH2D* hV0IMvPt[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];

		TH1D* hV0PtFit[V0consts::NSPECIES][V0consts::NTYPE][V0consts::NMULTI][V0consts::NSPHERO];
};
#endif