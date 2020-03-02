// AliAnalysisPIDParticle wrapper class
// OliverM 2019 Lund

#ifndef __MyParticle__
#define __MyParticle__
#include <AliAnalysisPIDParticle.h>

#include <AliESDtrack.h>
#include "TObject.h"

class AliAnalysisPIDParticle;

class MyParticle: public TObject {

	public:
		MyParticle() { }
		MyParticle(AliESDtrack* p) : mAliParticle(p) { }
		~MyParticle() { }
		Float_t GetPt() 						const { return 1;};//mAliParticle->GetPt();};
		Float_t GetPx() 						const { return 1;};//mAliParticle->GetPt()*TMath::Cos(mAliParticle->GetPhi());};
		Float_t GetPy() 						const { return 1;};//mAliParticle->GetPt()*TMath::Sin(mAliParticle->GetPhi());};
		Float_t GetPhi() 						const { return 1;};//mAliParticle->GetPhi();};
		Float_t GetEta() 						const { return 1;};//mAliParticle->GetEta();};
		Float_t GetY()	 						const { return 1;};//mAliParticle->GetY();};
		Int_t GetSign() const;// 						const { return mAliParticle->GetSign();};	// needs fix
		Int_t GetLabel()						const { return 1;};//mAliParticle->GetLabel();};
		Int_t GetPdgCode()						const { return 1;};//mAliParticle->GetPdgCode();};
		Int_t GetMotherPdgCode()				const { return 1;};//mAliParticle->GetMotherPdgCode();};
		Bool_t GetIsPrimary()					const { return 1;};//mIsPrimary;};

		Bool_t SetIsPrimary(Bool_t val)			{ mIsPrimary = val;};

		// add safety measure to getters for if mAliParticle is an invalid pointer
		
		ClassDef(MyParticle,1);

	private:
		AliESDtrack* mAliParticle;
		Bool_t mIsPrimary = 1;
};
#endif