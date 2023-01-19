#include "MyEvent.h"

#if INPUTFORMAT == 2
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
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

  	// INEL>0
  	if (AliESDtrackCuts::GetReferenceMultiplicity(mAliEvent, AliESDtrackCuts::kTracklets,1.0) < 1) return 0;

  	// Incomplete DAQ
  	if (mAliEvent->IsIncompleteDAQ()) return 0;

  	// SPD vs tracklet rejection of background
  	// copied AliAnalysisUtils::IsSPDClusterVsTrackletBG()
	Int_t nClustersLayer0 = mAliEvent->GetNumberOfITSClusters(0);
	Int_t nClustersLayer1 = mAliEvent->GetNumberOfITSClusters(1);
	Int_t nTracklets      = mAliEvent->GetMultiplicity()->GetNumberOfTracklets();
	const Float_t ASPDCvsTCut(65.);
	const Float_t BSPDCvsTCut(4.);
	if (nClustersLayer0 + nClustersLayer1 > ASPDCvsTCut + nTracklets*BSPDCvsTCut) return 0;

	// Check other conditions similarly to using event flags in trees
  	if (!CheckFlag()) return 0;
  	//if((fEventFlags&fgFlagToCheck)!=fgFlagToCheck) return kFALSE;
	// set flag and check it  
	return 1;
  
}

Bool_t MyEvent::CheckFlag() {
	
	AliMultSelection *ams = (AliMultSelection*)mAliEvent->FindListObject("MultSelection");
	if (mFlagToCheck & kNotPileupInSPD &&
			//!ams->GetThisEventIsNotPileup()) return 0;
			mAliEvent->IsPileupFromSPD()) return 0;
	if (mFlagToCheck & kNotPileupInSPD &&
			mAliEvent->IsPileupFromSPDInMultBins()) return 0;
	if (mFlagToCheck & kVertexSelected2015pp ||
		mFlagToCheck & kSPDandTrkVtxExists ||
		mFlagToCheck & kPassProximityCut ) mVertexSelected2015pp = CalculateVertexSelection(1,0,1);

	if (mFlagToCheck & kVertexSelected2015pp && !mVertexSelected2015pp) return 0;


	return 1;		   
  
}

Bool_t MyEvent::CalculateVertexSelection(Bool_t checkSPDres, //enable check on vtx resolution
	Bool_t requireSPDandTrk,	//ask for both trk and SPD vertex
	Bool_t checkProximity )		//apply cut on relative position of spd and trk verteces
{

	if (!mAliEvent) return 0;
	const AliESDVertex * trkVertex = mAliEvent->GetPrimaryVertexTracks();
	const AliESDVertex * spdVertex = mAliEvent->GetPrimaryVertexSPD();
	Bool_t hasSPD = spdVertex->GetStatus();
	Bool_t hasTrk = trkVertex->GetStatus();
  
	if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;
	  
  	//reject events if none between the SPD or track verteces are available
 	//if no trk vertex, try to fall back to SPD vertex;
	if (!hasTrk) {
		if (!hasSPD) return kFALSE;
    	//on demand check the spd vertex resolution and reject if not satisfied
		if (checkSPDres) {//} && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
			if (!spdVertex) return kFALSE;
    		if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
		}    
	} else {
		if (hasSPD) {
		//if enabled check the spd vertex resolution and reject if not satisfied
    	//if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
			if (checkSPDres) {// && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
				if (!spdVertex) return kFALSE;
    			if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
			}   	
			if (checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5) return kFALSE;
    	}
  	}

  //Not needed here, done separately
  //Cut on the vertex z position
  //const AliESDVertex * vertex = esd->GetPrimaryVertex();
  //if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
  return kTRUE;
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