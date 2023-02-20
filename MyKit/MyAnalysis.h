// MyKit analysis base abstract class
// OliverM 2019 Lund

#ifndef MYANALYSIS_H
#define MYANALYSIS_H
#include "TObject.h"

class MyHandler; // forward declaration
class TDirectoryFile;

class MyAnalysis: public TObject {

	public:
		MyAnalysis() { }	
		~MyAnalysis() { }
		Int_t SetHandler(MyHandler* h = 0);
		Int_t SetDirectory();
		Int_t TakeoverHistograms(const Char_t* sourceName);

		Bool_t IsRun()				const {return mRun;};
		void SetRun(Bool_t run)		{ mRun = run;};

		virtual Int_t Init() = 0;				// must be defined in daughters
		virtual Int_t Make(Int_t iEv) = 0;		// must be defined in daughters
		virtual Int_t Finish() = 0;				// must be defined in daughters

		TDirectoryFile* dirFile()	const {return mDirFile;};
		
		
		ClassDef(MyAnalysis,1);

	protected:
		MyHandler* mHandler;
		TDirectoryFile* mDirFile;
		Bool_t mRun = true;
};
#endif