// AliAnalysisPIDTrack wrapper class
// OliverM 2019 Lund

#ifndef __MyTrack__
#define __MyTrack__
#include <AliAnalysisPIDTrack.h>
#include "TObject.h"

class AliAnalysisPIDTrack;

class MyTrack: public TObject {

	public:
		MyTrack() { }
		MyTrack(AliAnalysisPIDTrack* tr) : mAliTrack(tr) { }
		~MyTrack() { }
		Float_t GetPt() 		const { return mAliTrack->GetPt();};

		// add safety measure to getters for if mAliTrack is an invalid pointer
		
		ClassDef(MyTrack,1);

	private:
		AliAnalysisPIDTrack* mAliTrack; 
};
#endif