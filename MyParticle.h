// AliAnalysisPIDParticle wrapper class
// OliverM 2019 Lund

#ifndef __MyParticle__
#define __MyParticle__
#include <AliAnalysisPIDParticle.h>
#include "TObject.h"

class AliAnalysisPIDParticle;

class MyParticle: public TObject {

	public:
		MyParticle() { }
		MyParticle(AliAnalysisPIDParticle* p) : mAliParticle(p) { }
		~MyParticle() { }
		Float_t GetPt() 						const { return mAliParticle->GetPt();};
		Float_t GetEta() 						const { return mAliParticle->GetEta();};
		Int_t GetLabel()						const { return mAliParticle->GetLabel();};
		Int_t GetPdgCode()						const { return mAliParticle->GetPdgCode();};

		// add safety measure to getters for if mAliParticle is an invalid pointer
		
		ClassDef(MyParticle,1);

	private:
		AliAnalysisPIDParticle* mAliParticle; 
};
#endif