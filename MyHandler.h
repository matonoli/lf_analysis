// MyKit analysis handler class
// OliverM 2019 Lund

#ifndef __MyHandler__
#define __MyHandler__
#include "TObject.h"

class MyAnalysis;	// forward declaration
class TChain;
class TFile;

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
		TChain* Chain() 		{return mChain;};
		
		ClassDef(MyHandler,1);

	protected:
		
		MyAnalysis* mAnalysis;  // should be changed into an array to allow multiple analyses
		TChain* mChain;
		TFile* mFile;

		Bool_t mFlagMC;
		Bool_t mFlagHist;


};
#endif