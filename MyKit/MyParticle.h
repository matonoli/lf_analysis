// AliAnalysisPIDParticle wrapper class
// OliverM 2019 Lund

#include "inputformat.h"

#ifndef __MyParticle__
#define __MyParticle__
#include <AliAnalysisPIDParticle.h>

#include <TParticle.h>
#include "TObject.h"

#if INPUTFORMAT == 1
	typedef AliAnalysisPIDParticle AnyParticle;
#elif INPUTFORMAT == 2
	typedef TParticle AnyParticle;
#endif


class MyParticle: public TObject {

	public:
		MyParticle() { }
		MyParticle(AnyParticle* p) : mAliParticle(p) { }
		~MyParticle() { }
		Int_t GetSign() const;// 						const { return mAliParticle->GetSign();};	// needs fix
		Bool_t GetIsPrimary()					const { return mIsPrimary;};

		Bool_t SetIsPrimary(Bool_t val)			{ mIsPrimary = val;};
#if INPUTFORMAT == 1
		Float_t GetPt() 						const { return mAliParticle->GetPt();};
		Float_t GetPx() 						const { return mAliParticle->GetPt()*TMath::Cos(mAliParticle->GetPhi());};
		Float_t GetPy() 						const { return mAliParticle->GetPt()*TMath::Sin(mAliParticle->GetPhi());};
		Float_t GetPhi() 						const { return mAliParticle->GetPhi();};
		Float_t GetEta() 						const { return mAliParticle->GetEta();};
		Float_t GetY()	 						const { return mAliParticle->GetY();};
		Int_t GetLabel()						const { return mAliParticle->GetLabel();};
		Int_t GetPdgCode()						const { return mAliParticle->GetPdgCode();};
		Int_t GetMotherPdgCode()				const { return mAliParticle->GetMotherPdgCode();};
#elif INPUTFORMAT == 2
		Float_t GetPt() 						const { return mAliParticle->Pt();};
		Float_t GetPx() 						const { return mAliParticle->Px();};
		Float_t GetPy() 						const { return mAliParticle->Py();};
		Float_t GetPhi() 						const { return mAliParticle->Phi();};
		Float_t GetEta() 						const { return mAliParticle->Eta();};
		Float_t GetY()	 						const { return mAliParticle->Y();};
		Int_t GetLabel()						const { return mLabel;};// label is fMCEvent->Particle(label)
		void SetLabel(Int_t lab)				{ mLabel = lab;};
		Int_t GetPdgCode()						const { return mAliParticle->GetPdgCode();};
		//Int_t GetMotherPdgCode()				const { return 0;}; // needs to be accessed manually
		

#endif	

		// add safety measure to getters for if mAliParticle is an invalid pointer
		
		ClassDef(MyParticle,1);

	private:
		AnyParticle* mAliParticle;
		Bool_t mIsPrimary = 1;

		Int_t mLabel = 0;
};
#endif

/* todo
get particles from handler
add setlabel to analysis code
add extra conditions for selection - is primary
am i calculating mc spherocity in analysis code only from primaries ??

add GetMCLabel to myVO.h
check what is v0.ismcprimary

there is possibly double counting in reconstructed v0s (often there are more with same label)
*/