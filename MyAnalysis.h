// MyKit analysis base abstract class
// OliverM 2019 Lund

#ifndef __MyAnalysis__
#define __MyAnalysis__
#include "TObject.h"

class MyHandler; // forward declaration

class MyAnalysis: public TObject {

	public:
		MyAnalysis() { }	
		~MyAnalysis() { }
		Int_t SetHandler(MyHandler* h = 0);
		virtual Int_t Init() = 0;				// must be defined in daughters
		virtual Int_t Make(Int_t iEv) = 0;		// must be defined in daughters
		virtual Int_t Finish() = 0;				// must be defined in daughters
		
		ClassDef(MyAnalysis,1);

	protected:
		MyHandler* mHandler;		
};
#endif