// AliAnalysisPIDCascadeEvent wrapper class
// OliverM 2019 Lund

#include "compInstructions.h"

#ifndef MYEVENT_H
#define MYEVENT_H
#include <AliAnalysisPIDCascadeEvent.h>
#include <AliESDEvent.h>
#include "TObject.h"

#if INPUTFORMAT == 1
	typedef AliAnalysisPIDCascadeEvent AnyEvent;
#elif INPUTFORMAT == 2
	typedef AliESDEvent AnyEvent;
	enum EventFlags_t {
	    kNotPileupInSPD = 1,
	    kNotPileupInMV = 2,
	    kNotPileupInMB = 4,
	    kINELgtZERO = 8,
	    kNoInconsistentVtx = 16,
	    kNoV0Asym = 32,
	    kVertexSelected2015pp=64,
	    kSPDandTrkVtxExists=128,
	    kPassProximityCut=256,
	    kAll = 511
  	};
#endif

class MyEvent: public TObject {


	public:
		MyEvent() : mAliEvent(0) { }
		MyEvent(AnyEvent* ev) : mAliEvent(ev) { }	
		~MyEvent() { }
#if INPUTFORMAT == 1
		Float_t GetZ()						const { return mAliEvent->GetVertexZ();};
		Float_t GetV0MCentrality()			const { return mAliEvent->GetV0Mmultiplicity();};
		Int_t GetRefMult()					const { return mAliEvent->GetReferenceMultiplicity();};
		Bool_t IsGoodAliEvent() 			const { return mAliEvent->AcceptEvent(kFALSE,0);}; // checks flag + collision candidate
		Bool_t HasVertex()					const { return mAliEvent->HasVertex();};
		
		void SetCheckFlag(Int_t flag)	{ mAliEvent->SetCheckFlag(flag);};
		Bool_t CheckFlag() 					const { return mAliEvent->CheckFlag();};
#elif INPUTFORMAT == 2
		Float_t GetZ()						const { return mAliEvent->GetPrimaryVertexTracks()->GetZ();};
		Float_t GetV0MCentrality()			{ return (mV0Mmultiplicity==0) ? CalculateV0Mmultiplicity() : mV0Mmultiplicity;};
		Float_t CalculateV0Mmultiplicity();
		Int_t GetRefMult()					{ return (mRefMult==0) ? CalculateRefMult() : mRefMult;};
		Int_t CalculateRefMult();
		Bool_t IsGoodAliEvent(); // const { return mAliEvent->AcceptEvent(kFALSE,0);}; // checks flag + collision candidate
		Bool_t HasVertex()					{ return (!mHasVertexCal) ? CalculateHasVertex() : mHasVertex;};
		Bool_t CalculateHasVertex();

		void SetCheckFlag(Int_t flag)	{ mFlagToCheck = flag;};
		Bool_t CheckFlag();
		Bool_t CalculateVertexSelection(Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
		//void CalculateFlags();
#endif
#if INPUTFORMAT == 1	
		Bool_t IsPileupFromSPD() 			const { return mAliEvent->IsPileup();};
		Bool_t AcceptVertex() 				const { return mAliEvent->AcceptVertex();};
		UChar_t GetCentralityQuality() 		const { return mAliEvent->GetCentralityQuality();};
		Bool_t IsCollisionCandidate() 		const { return mAliEvent->IsCollisionCandidate();};
		
		Int_t GetEventFlags()				const { return mAliEvent->GetEventFlags();};
		Int_t GetCheckFlag()				const { return mAliEvent->GetCheckFlag();};
		
#endif
		ClassDef(MyEvent,1);

	private:
		AnyEvent* mAliEvent;

#if INPUTFORMAT == 2
		Float_t mV0Mmultiplicity = 0;
		Float_t mRefMult = 0;
		Int_t mEventFlags = 0;
		Int_t mFlagToCheck = 0;
		Bool_t mVertexSelected2015pp = 0;

		Bool_t mHasVertex = 0;
		Bool_t mHasVertexCal = 0;
#endif		

};
#endif