#include <TVector3.h>
#include <TParticle.h>

#include "MyV0.h"
#include "MyTrack.h"



ClassImp(MyV0)

MyV0::MyV0(AnyV0* v0) : mAliV0(v0), 
	mAPcalculated(false), mPhicalculated(false) { 

	for (int iSp = 0; iSp < 4; iSp++) mYcalculated[iSp] = false;
	mAP[0] = 0;
	mAP[1] = 0;
	mPhi = 0;
}

#if INPUTFORMAT == 2
Double_t MyV0::GetIMK0s() {
	//mAliV0->ChangeMassHypothesis(310);			// THIS ALSO CHANGES PDG ID
	return (mAliV0->GetEffMass(2,2)-0.497614);		// Using AliEsdV0::GetEffMass(int, int) instead
													// 0: e, 1: mu-, 2: pi+, 3: K+, 4: p
}


Double_t MyV0::GetIML() {
	//mAliV0->ChangeMassHypothesis(3122);
	return (mAliV0->GetEffMass(4,2)-1.11568);
}

Double_t MyV0::GetIMLbar() {
	//mAliV0->ChangeMassHypothesis(-3122);
	return (mAliV0->GetEffMass(2,4)-1.11568);
}

Double_t* MyV0::CalculateAP() {

	mAP[0] = mAliV0->AlphaV0();
	mAP[1] = mAliV0->PtArmV0();
	return mAP;
}

Bool_t MyV0::HasFastSignal() {
  
  // logical or
  Bool_t hasFast = kFALSE;
  if ((this->GetNegTrack()->GetStatus() & AliESDtrack::kITSrefit)) hasFast = kTRUE;
  if ((this->GetPosTrack()->GetStatus() & AliESDtrack::kITSrefit)) hasFast = kTRUE;
  
  MyTrack trP(this->GetPosTrack());
  MyTrack trN(this->GetNegTrack());
  if (trP.HasTOFPID()) hasFast = kTRUE;
  if (trN.HasTOFPID()) hasFast = kTRUE;
    
  return hasFast;
}


Int_t MyV0::GetMCLabel() const {

	TParticle* pd = mHandler->mcstack()->Particle(TMath::Abs(this->GetPosTrack()->GetLabel()));
	return pd->GetFirstMother();
}

Int_t MyV0::GetMCPdgCode() {

	TParticle* pd = mHandler->mcstack()->Particle(this->GetMCLabel());
	return pd->GetPdgCode();
}

#endif


#if INPUTFORMAT==1
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
		//cout << "rap for sp " << Sp << " is " << mY[Sp] << endl;
		return mY[Sp]; 
	}
}

Bool_t MyV0::HasFastSignal() {
  
  // logical or
  Bool_t hasFast = kFALSE;
  if ((this->GetNegTrack()->GetStatus() & AliESDtrack::kITSrefit)) hasFast = kTRUE;
  if ((this->GetPosTrack()->GetStatus() & AliESDtrack::kITSrefit)) hasFast = kTRUE;
  if (this->GetNegTrack()->HasTOFPID()) hasFast = kTRUE;
  if (this->GetPosTrack()->HasTOFPID()) hasFast = kTRUE;
    
  return hasFast;
}
#endif
