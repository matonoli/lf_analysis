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
		Bool_t IsGoodAliEvent() const { return mAliEvent->AcceptEvent(1,0);}; // vertex position
		
		ClassDef(MyEvent,1);

	private:
		AliAnalysisPIDEvent* mAliEvent;
};
#endif