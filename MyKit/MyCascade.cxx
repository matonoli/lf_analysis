#include <TVector3.h>
#include <TParticle.h>

#include "AliKFParticle.h"
#include "AliExternalTrackParam.h"
#include "AliKFVertex.h"

#include "MyCascade.h"
#include "MyTrack.h"



ClassImp(MyCascade)

MyCascade::MyCascade(AnyCascade* cascade) : mAliCascade(cascade), 
	mPhicalculated(false) { 

	mPhi = 0;

}

#if INPUTFORMAT==1

Double_t MyCascade::GetPhi() {

	if (mPhicalculated) {
		return mPhi;
	} 
	else {
		TVector3 momNeg, momPos, momBach, momTot;
		momNeg.SetPtEtaPhi(this->GetV0()->GetNegAnalysisTrack()->GetPt(),this->GetV0()->GetNegAnalysisTrack()->GetEta(),this->GetV0()->GetNegAnalysisTrack()->GetPhi());
		momPos.SetPtEtaPhi(this->GetV0()->GetPosAnalysisTrack()->GetPt(),this->GetV0()->GetPosAnalysisTrack()->GetEta(),this->GetV0()->GetPosAnalysisTrack()->GetPhi());
		momBach.SetPtEtaPhi(this->GetBachTrack()->GetPt(),this->GetBachTrack()->GetEta(),this->GetBachTrack()->GetPhi());
		
		momTot = momNeg + momPos + momBach;

		mPhi = momTot.Phi();
		mPhicalculated = true;
		return mPhi; 
	}
}

#endif
