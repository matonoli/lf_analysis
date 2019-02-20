// AliAnalysisPIDV0 wrapper class
// OliverM 2019 Lund

#ifndef __MyV0__
#define __MyV0__
#include <AliAnalysisPIDV0.h>
#include "TObject.h"

class AliAnalysisPIDTrack; // forward declaration

class MyV0: public TObject {

	public:
		MyV0() { }
		MyV0(AliAnalysisPIDV0* v0) : mAliV0(v0) { }	
		~MyV0() { }
		Float_t GetPt() 					const { return mAliV0->GetPt();};
		Float_t GetEta() 					const { return mAliV0->GetEta();};
		AliAnalysisPIDTrack* GetTrackPos() 	const { return mAliV0->GetPosAnalysisTrack();};
		
		ClassDef(MyV0,1);

	private:
		AliAnalysisPIDV0* mAliV0;
};
#endif