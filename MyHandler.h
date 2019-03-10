// MyKit analysis handler class
// OliverM 2019 Lund

#ifndef __MyHandler__
#define __MyHandler__
#include "TObject.h"

class MyAnalysis;	// forward declaration
class TChain;
class TFile;
class TCanvas;
class TH1D;
class TLegend;

class MyHandler: public TObject {

	public:
		MyHandler() { }	
		~MyHandler() { }
		void AddAnalysis(MyAnalysis* ana = 0);
		Int_t Init();
		Int_t LoadInput(const Char_t *inputFile, const Char_t *chainName);
		Int_t Make(Int_t iEv);
		Int_t Finish();
		
		Bool_t GetFlagMC() 		const {return mFlagMC;};
		Bool_t GetFlagHist()	const {return mFlagHist;};
		TChain* Chain() 		{return mChain;};

		void DrawCut(Double_t cut, Int_t direction, TCanvas* can);
		void MakeNiceHistogram(TH1D* h, Int_t col);
		void MakeNiceLegend(TLegend* leg, Float_t size = 0.07, Int_t columns = 2);
		
		ClassDef(MyHandler,1);

	protected:
		
		MyAnalysis* mAnalysis;  // should be changed into an array to allow multiple analyses
		TChain* mChain;
		TFile* mFile;

		Bool_t mFlagMC;
		Bool_t mFlagHist;


};
#endif