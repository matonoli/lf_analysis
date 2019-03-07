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
		Float_t GetPt() 						const { return mAliTrack->GetPt();};
		Float_t GetEta() 						const { return mAliTrack->GetEta();};
		Float_t GetNSigmaPionTPC() 				const { return mAliTrack->GetNSigmaPionTPC();};
		Float_t GetNSigmaProtonTPC() 			const { return mAliTrack->GetNSigmaProtonTPC();};
		Float_t GetDCApvXY() 					const { return mAliTrack->GetImpactParameter(0);};

		// add safety measure to getters for if mAliTrack is an invalid pointer
		
		ClassDef(MyTrack,1);

	private:
		AliAnalysisPIDTrack* mAliTrack; 
};
#endif