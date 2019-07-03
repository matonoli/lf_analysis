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
		Float_t GetPx() 						const { return mAliParticle->GetPt()*TMath::Cos(mAliParticle->GetPhi());};
		Float_t GetPy() 						const { return mAliParticle->GetPt()*TMath::Sin(mAliParticle->GetPhi());};
		Float_t GetPhi() 						const { return mAliParticle->GetPhi();};
		Float_t GetEta() 						const { return mAliParticle->GetEta();};
		Float_t GetSign() 						const { return mAliParticle->GetSign();};	// needs fix
		Int_t GetLabel()						const { return mAliParticle->GetLabel();};
		Int_t GetPdgCode()						const { return mAliParticle->GetPdgCode();};
		Int_t GetMotherPdgCode()				const { return mAliParticle->GetMotherPdgCode();};
		Bool_t GetIsPrimary()					const { return mIsPrimary;};

		Bool_t SetIsPrimary(Bool_t val)			{ mIsPrimary = val;};

		// add safety measure to getters for if mAliParticle is an invalid pointer
		
		ClassDef(MyParticle,1);

	private:
		AliAnalysisPIDParticle* mAliParticle;
		Bool_t mIsPrimary = 1;
};
#endif