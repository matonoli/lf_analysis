// MyKit analysis handler class
// OliverM 2019 Lund

#ifndef __MyHandler__
#define __MyHandler__
#include "TObject.h"
#include "TClonesArray.h"

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
//class TClonesArray;

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
		MyAnalysis* analysis(Int_t iAna)	const {return mAnalysis[iAna];};
		
		AliAnalysisPIDEvent* event()		const {return mEvent;};
		TClonesArray* tracks()				const {return bTracks;};
		TClonesArray* v0s()					const {return bV0s;};
		TClonesArray* particles()			const {return bParticles;};
		AliAnalysisPIDTrack* track(Int_t i)			const {return (AliAnalysisPIDTrack*)bTracks->At(i);};
		AliAnalysisPIDV0* v0(Int_t i)				const {return (AliAnalysisPIDV0*)bV0s->At(i);};
		AliAnalysisPIDParticle* particle(Int_t i)	const {return (AliAnalysisPIDParticle*)bParticles->At(i);};
		

		void DrawCut(Double_t cut, Int_t direction, TVirtualPad* can);
		void MakeNiceHistogram(TH1D* h, Int_t col);
		void MakeNiceLegend(TLegend* leg, Float_t size = 0.07, Int_t columns = 2);
		void MakeRatioPlot(TH1D* hn, TH1D* hd, TCanvas* c, Double_t low, Double_t high);
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


		AliAnalysisPIDEvent* mEvent = 0;
		TClonesArray* bTracks = 0;
		TClonesArray* bV0s = 0;
		TClonesArray* bParticles = 0;

		Bool_t mFlagMC = 0;
		Bool_t mFlagHist = 0;
		Bool_t mRebinPt;


};
#endif