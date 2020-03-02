// MyKit analysis handler class
// OliverM 2019 Lund

#include "inputformat.h"

#ifndef __MyHandler__
#define __MyHandler__
#include "TObject.h"
#include "TClonesArray.h"
#include <AliESDEvent.h>

class MyAnalysis;	// forward declaration
class TChain;
class TFile;
class TCanvas;
class TPad;
class TVirtualPad;
class TH1D;
class TLegend;
class TDirectory;
class TString;
class TROOT;
class AliAnalysisPIDEvent;
class AliAnalysisPIDTrack;
class AliAnalysisPIDV0;
class AliAnalysisPIDParticle;

class AliESDEvent;
class AliESDtrack;
class AliESDv0;
class AliESDtrackCuts;

#if INPUTFORMAT == 1
	typedef AliAnalysisPIDEvent AnyEvent;
#elif INPUTFORMAT == 2
	typedef AliESDEvent AnyEvent;
#endif


class MyHandler: public TObject {

	public:
		MyHandler();	
		~MyHandler() { }
		void AddAnalysis(MyAnalysis* ana = 0);
		Int_t Init();
		Int_t LoadInputTree(const Char_t *inputFile, const Char_t *chainName);
		Int_t LoadInputHist(const Char_t *inputFile);
		Int_t SetDirectory(TDirectory* d) { mDir = d;};
		Int_t SetROOT(TROOT* r) { mROOT = r;};
		Int_t Make(Int_t iEv);
		Int_t Finish();
		void SetOutputName(const Char_t *name) { mOutName = TString(name);};
		
		Bool_t GetFlagMC() 				const {return mFlagMC;};
		Bool_t GetFlagHist()			const {return mFlagHist;};
		void RebinPt(Bool_t opt)		{ mRebinPt = opt;};
		Bool_t IsRebinPt()				const { return mRebinPt;};
		
		TChain* chain() 					const {return mChain;};
		TDirectory* directory() 			const {return mDir;};
		TROOT* root()	 					const {return mROOT;};
		TFile* file() 						const {return mFile;};
		TFile* filehist() 					const {return mFileHist;};
		Int_t nAnalysis()					const {return nAna;};
		Int_t getNAnalyses()				const {return nAnalyses;};
		MyAnalysis* analysis(Int_t iAna)	const {return mAnalysis[iAna];};

		void SetEvent(AliESDEvent* ev)		{mEvent = ev;};	
		
		AnyEvent* event()					const {return mEvent;};
		#if INPUTFORMAT == 1
		TClonesArray* tracks()						const {return bTracks;};
		TClonesArray* v0s()							const {return bV0s;};
		TClonesArray* particles()					const {return bParticles;};
		AliAnalysisPIDTrack* track(Int_t i)			const {return (AliAnalysisPIDTrack*)bTracks->At(i);};
		AliAnalysisPIDV0* v0(Int_t i)				const {return (AliAnalysisPIDV0*)bV0s->At(i);};
		AliAnalysisPIDParticle* particle(Int_t i)	const {return (AliAnalysisPIDParticle*)bParticles->At(i);};
		Int_t getNtracks()							const {return bTracks->GetEntriesFast();};
		Int_t getNv0s()								const {return bV0s->GetEntriesFast();};
		Int_t getNparticles()						const {return bParticles->GetEntriesFast();};
		
		#elif INPUTFORMAT == 2
		AliESDtrack* track(Int_t i)				const {return (AliESDtrack*)mEvent->GetTrack(i);};
		AliESDv0* v0(Int_t i)					const {return (AliESDv0*)mEvent->GetV0(i);};
		AliESDtrack* particle(Int_t i)			const {return (AliESDtrack*)mEvent->GetTrack(i);};
		Int_t getNtracks()						const {return mEvent->GetNumberOfTracks();};
		Int_t getNv0s()							const {return mEvent->GetNumberOfV0s();};
		Int_t getNparticles()					const {return mEvent->GetNumberOfTracks();};	//needs a fix

		void SetupTrackCuts();
		AliESDtrackCuts* trackCuts2010()		const {return mTrackCuts2010;};
		AliESDtrackCuts* trackCuts2011()		const {return mTrackCuts2011;};
		AliESDtrackCuts* trackCutsTPCOnly()		const {return mTrackCutsTPCOnly;};
		AliESDtrackCuts* trackCuts2011sys()		const {return mTrackCuts2011sys;};
		AliESDtrackCuts* trackCutsV0d()			const {return mTrackCutsV0d;};
		#endif
		

		void DrawCut(Double_t cut, Int_t direction, TVirtualPad* can);
		void MakeNiceHistogram(TH1D* h, Int_t col);
		void MakeNiceLegend(TLegend* leg, Float_t size = 0.07, Int_t columns = 2);
		void MakeRatioPlot(TH1D* hn, TH1D* hd, TCanvas* c, Double_t low, Double_t high, Double_t lowx, Double_t highx);
		void MakeZoomPlot(TH1D* h, TCanvas* c, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax);

		Double_t DeltaPhi(Double_t phi1, Double_t phi2);
		
		ClassDef(MyHandler,1);

	protected:
		
		MyAnalysis* mAnalysis[5];
		Int_t nAnalyses = 0;
		TChain* mChain;
		TFile* mFile;
		TFile* mFileHist;
		TDirectory* mDir;
		TROOT* mROOT;
		TString mOutName;

		Int_t nAna = 0;


		AnyEvent* mEvent = 0;
		TClonesArray* bTracks = 0;
		TClonesArray* bV0s = 0;
		TClonesArray* bParticles = 0;

		Bool_t mFlagMC = 0;
		Bool_t mFlagHist = 0;
		Bool_t mRebinPt;

		AliESDtrackCuts* mTrackCuts2010 = 0;
		AliESDtrackCuts* mTrackCuts2011 = 0;
		AliESDtrackCuts* mTrackCutsTPCOnly = 0;
		AliESDtrackCuts* mTrackCuts2011sys = 0;
		AliESDtrackCuts* mTrackCutsV0d = 0;


};
#endif