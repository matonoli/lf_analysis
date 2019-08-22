#include <TVector3.h>

#include "MyV0.h"


ClassImp(MyV0)

MyV0::MyV0(AliAnalysisPIDV0* v0) : mAliV0(v0), 
	mAPcalculated(false), mPhicalculated(false) { 

	mAP[0] = 0;
	mAP[1] = 0;
	mPhi = 0;
}

Double_t* MyV0::CalculateAP() {

	if (mAPcalculated) {
		return mAP;
	} 
	else {
		TVector3 momNeg, momPos, momTot;
		momNeg.SetPtEtaPhi(this->GetNegTrack()->GetPt(),this->GetNegTrack()->GetEta(),this->GetNegTrack()->GetPhi());
		momPos.SetPtEtaPhi(this->GetPosTrack()->GetPt(),this->GetPosTrack()->GetEta(),this->GetPosTrack()->GetPhi());
		//momTot.SetPtEtaPhi(v0->GetPt(),v0->GetEta(),v0->GetPhi());
		momTot = momNeg + momPos;

		Double_t lQlNeg = momNeg.Dot(momTot)/momTot.Mag();
		Double_t lQlPos = momPos.Dot(momTot)/momTot.Mag();

	 	mAP[0] = (lQlPos - lQlNeg)/(lQlPos + lQlNeg);
		mAP[1] = momPos.Perp(momTot);
		mAPcalculated = true;
		return mAP; 
	}
}

Double_t MyV0::GetPhi() {

	if (mPhicalculated) {
		return mPhi;
	} 
	else {
		TVector3 momNeg, momPos, momTot;
		momNeg.SetPtEtaPhi(this->GetNegTrack()->GetPt(),this->GetNegTrack()->GetEta(),this->GetNegTrack()->GetPhi());
		momPos.SetPtEtaPhi(this->GetPosTrack()->GetPt(),this->GetPosTrack()->GetEta(),this->GetPosTrack()->GetPhi());
		//momTot.SetPtEtaPhi(v0->GetPt(),v0->GetEta(),v0->GetPhi());
		momTot = momNeg + momPos;

		mPhi = momTot.Phi();
		mPhicalculated = true;
		return mPhi; 
	}
}

Double_t MyV0::CalculateY(Int_t Sp) {

	if (mYcalculated[Sp]) {
		return mY[Sp];
	} 
	else {
		if (Sp > 3) return -99;
		TLorentzVector mom;
		Double_t V0mass[] = { 0., 0.497614, 1.11568, 1.11568};
		mom.SetPtEtaPhiM(this->GetPt(),this->GetEta(),this->GetPhi(),V0mass[Sp]);

		mY[Sp] = mom.Rapidity();
		mYcalculated[Sp] = true;
		return mY[Sp]; 
	}
}
