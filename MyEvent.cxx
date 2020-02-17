#include "MyEvent.h"

#if INPUTFORMAT == 2
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#endif

ClassImp(MyEvent)

#if INPUTFORMAT == 2
Float_t MyEvent::CalculateV0Mmultiplicity() {
	
	
	AliMultSelection *ams = (AliMultSelection*)mAliEvent->FindListObject("MultSelection");
	if (!ams) mV0Mmultiplicity = 0;
	else mV0Mmultiplicity = ams->GetMultiplicityPercentile("V0M");
	
	return mV0Mmultiplicity;
}

Int_t MyEvent::CalculateRefMult() {
	
	mRefMult = AliESDtrackCuts::GetReferenceMultiplicity(mAliEvent, AliESDtrackCuts::kTrackletsITSTPC,0.8);
	/// Returns a negative value if no reliable estimate can be provided:
	///   -1 SPD vertex missing
	///   -2 SPD VertexerZ dispersion too large
	///   -3 Track vertex missing (not checked for kTracklets)
	///   -4 Distance between SPD and track vertices too large (not checked for kTracklets)
	return mRefMult;
}

Bool_t MyEvent::IsGoodAliEvent() {
	//emulate AcceptEvent
  	Int_t isEventSelected = (((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
  	//if (centralityQuality != 0) return kFALSE; //not needed
	if (!(isEventSelected & AliVEvent::kAny)) return 0;
  	if (!(isEventSelected & AliVEvent::kINT7)) return 0;
  	//if((fEventFlags&fgFlagToCheck)!=fgFlagToCheck) return kFALSE;
	// set flag and check it  
	return 1;
  
}

Bool_t MyEvent::CalculateHasVertex() {
	
	mHasVertex = true;
	if (mAliEvent->GetPrimaryVertexTracks()->GetNContributors() < 1) {
		if (mAliEvent->GetPrimaryVertexSPD()->GetNContributors() < 1) mHasVertex = false;
		TString vtxTyp = mAliEvent->GetPrimaryVertexSPD()->GetTitle();
    	Double_t cov[6]={0};
    	mAliEvent->GetPrimaryVertexSPD()->GetCovarianceMatrix(cov);
    	Double_t zRes = TMath::Sqrt(cov[5]);
    	if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) mHasVertex = false;
	}
	mHasVertexCal = true;
	return mHasVertex;
}
#endif