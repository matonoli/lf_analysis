// AliAnalysisPIDEvent wrapper class
// OliverM 2019 Lund

#ifndef __MyEvent__
#define __MyEvent__
#include <AliAnalysisPIDEvent.h>
#include "TObject.h"

class MyEvent: public TObject {

	public:
		MyEvent() { }
		MyEvent(AliAnalysisPIDEvent* ev) : mAliEvent(ev) { }	
		~MyEvent() { }
		Float_t GetZ() const { return mAliEvent->GetVertexZ();};
		Float_t GetV0MCentrality() const { return mAliEvent->GetV0Mmultiplicity();};
		Float_t GetRefMult() const { return mAliEvent->GetReferenceMultiplicity();};
		Bool_t IsGoodAliEvent() const { return mAliEvent->AcceptEvent(kFALSE,0);}; // checks flag + collision candidate
		Bool_t IsPileupFromSPD() const { return mAliEvent->IsPileup();};
		Bool_t HasVertex()		const { return mAliEvent->HasVertex();};
		Bool_t AcceptVertex() const { return mAliEvent->AcceptVertex();};
		UChar_t GetCentralityQuality() const { return mAliEvent->GetCentralityQuality();};
		Bool_t IsCollisionCandidate() const { return mAliEvent->IsCollisionCandidate();};
		
		void SetCheckFlag(Int_t flag)	{ mAliEvent->SetCheckFlag(flag);};
		Int_t GetEventFlags()	const { return mAliEvent->GetEventFlags();};
		Bool_t CheckFlag() const { return mAliEvent->CheckFlag();};
		Int_t GetCheckFlag()	const { return mAliEvent->GetCheckFlag();};
		
		ClassDef(MyEvent,1);

	private:
		AliAnalysisPIDEvent* mAliEvent;
};
#endif