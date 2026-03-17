#include <TVector3.h>
#include <TParticle.h>

#include "AliKFParticle.h"
#include "AliExternalTrackParam.h"
#include "AliKFVertex.h"

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

Double_t MyV0::GetCut(Int_t Cut) {
	
	MyTrack trP;
	MyTrack trN;
	if (Cut>cutsV0cuts) {
		trP = MyTrack(this->GetPosTrack());
		trN = MyTrack(this->GetNegTrack());
	}

	switch (Cut) {

		default: break;
		case NoCut : return 1;
		case DCAdd : return this->GetDCAdd(); break;
		case CPA : return this->GetCPA(); break;
		case RadiusL : return this->GetRadius(); break;
		case RadiusH : return this->GetRadius(); break;
		case FastSignal : return this->HasFastSignal(); break;
		case CompMassK0s : return this->GetIMK0s(); break;
		case CompMassL : return this->GetIML(); break;
		case CompMassLbar : return this->GetIMLbar(); break;
		case LifetimeK0s : return (this->GetPt() > 0 ) ? (0.497614+this->GetIMK0s())*this->GetRadius()/this->GetPt() : 0; break;
		case LifetimeL : return (this->GetPt() > 0 ) ? (1.11568+this->GetIML())*this->GetRadius()/this->GetPt() : 0; break;
		case LifetimeLbar : return (this->GetPt() > 0 ) ? (1.11568+this->GetIMLbar())*this->GetRadius()/this->GetPt() : 0; break;
		case cutsV0cuts : break;

		case NSigmaTPCposPiL : return trP.GetNSigmaPionTPC(); break;
		case NSigmaTPCposPiH : return trP.GetNSigmaPionTPC(); break;
		case NSigmaTPCposPrL : return trP.GetNSigmaProtonTPC(); break;
		case NSigmaTPCposPrH : return trP.GetNSigmaProtonTPC(); break;
		case NSigmaTPCnegPiL : return trN.GetNSigmaPionTPC(); break;
		case NSigmaTPCnegPiH : return trN.GetNSigmaPionTPC(); break;
		case NSigmaTPCnegPrL : return trN.GetNSigmaProtonTPC(); break;
		case NSigmaTPCnegPrH : return trN.GetNSigmaProtonTPC(); break;
		case DCAPVpos : return TMath::Abs(trP.GetDCApvXY()); break;
		case NClusterpos : return trP.GetTPCnc(); break;
		case NClusterFpos : return (trP.GetTPCNclsF()>0) ? trP.GetTPCnc()/trP.GetTPCNclsF() : 0; break;
		case DCAPVneg : return TMath::Abs(trN.GetDCApvXY()); break;
		case NClusterneg : return trN.GetTPCnc(); break;
		case NClusterFneg : return (trN.GetTPCNclsF()>0) ? trN.GetTPCnc()/trN.GetTPCNclsF() : 0; break;
		case cutsSizeof : break;
	}

	return 1;
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

Double_t MyV0::GetKFIMK0s() {
	
	AliKFParticle mother;
	const AliExternalTrackParam* paramP = mAliV0->GetParamP();
	const AliExternalTrackParam* paramN = mAliV0->GetParamN();

	if (paramP->GetParameter()[4] < 0)	{
		paramP = mAliV0->GetParamN();
		paramN = mAliV0->GetParamP();	}

	AliKFParticle daughterP( *(paramP), 211);
	AliKFParticle daughterN( *(paramN), -211);

	mother += daughterP;
	mother += daughterN;

	mother.SetField(mHandler->event()->GetMagneticField());
	AliKFVertex PrimaryVtxKF(*(mHandler->event()->GetPrimaryVertex()));
	mother.SetProductionVertex(PrimaryVtxKF);

	return (mother.GetMass() - 0.497614);
	
}

Double_t MyV0::GetKFIML() {
	
	AliKFParticle mother;
	const AliExternalTrackParam* paramP = mAliV0->GetParamP();
	const AliExternalTrackParam* paramN = mAliV0->GetParamN();

	if (paramP->GetParameter()[4] < 0)	{
		paramP = mAliV0->GetParamN();
		paramN = mAliV0->GetParamP();	}

	AliKFParticle daughterP( *(paramP), 2212);
	AliKFParticle daughterN( *(paramN), -211);

	mother += daughterP;
	mother += daughterN;

	mother.SetField(mHandler->event()->GetMagneticField());
	AliKFVertex PrimaryVtxKF(*(mHandler->event()->GetPrimaryVertex()));
	mother.SetProductionVertex(PrimaryVtxKF);

	return (mother.GetMass() - 1.11568);
	
}

Double_t MyV0::GetKFIMLbar() {
	
	AliKFParticle mother;
	const AliExternalTrackParam* paramP = mAliV0->GetParamP();
	const AliExternalTrackParam* paramN = mAliV0->GetParamN();

	if (paramP->GetParameter()[4] < 0)	{
		paramP = mAliV0->GetParamN();
		paramN = mAliV0->GetParamP();	}

	AliKFParticle daughterP( *(paramP), 211);
	AliKFParticle daughterN( *(paramN), -2212);

	mother += daughterP;
	mother += daughterN;

	mother.SetField(mHandler->event()->GetMagneticField());
	AliKFVertex PrimaryVtxKF(*(mHandler->event()->GetPrimaryVertex()));
	mother.SetProductionVertex(PrimaryVtxKF);

	return (mother.GetMass() - 1.11568);
	
}

Double_t* MyV0::CalculateAP() {

	mAP[0] = mAliV0->AlphaV0();
	mAP[1] = mAliV0->PtArmV0();
	return mAP;
}

Int_t MyV0::HasFastSignal() {
  
  Int_t howMany = 0;
  
  // logical or
  if ((this->GetNegTrack()->GetStatus() & AliESDtrack::kITSrefit)) howMany++;
  if ((this->GetPosTrack()->GetStatus() & AliESDtrack::kITSrefit)) howMany++;
  
  MyTrack trP(this->GetPosTrack());
  MyTrack trN(this->GetNegTrack());
  if (trP.HasTOFPID()) howMany++;
  if (trN.HasTOFPID()) howMany++;
    
  return howMany;
}


Int_t MyV0::GetMCLabel() const {

	TParticle* pd = mHandler->mcstack()->Particle(TMath::Abs(this->GetPosTrack()->GetLabel()));
	TParticle* nd = mHandler->mcstack()->Particle(TMath::Abs(this->GetNegTrack()->GetLabel()));
	if (pd->GetFirstMother() == nd->GetFirstMother() && pd->GetFirstMother() > 0)
		return pd->GetFirstMother();
	else
		return 0;
}

Int_t MyV0::GetMCPdgCode() {

	TParticle* pd = mHandler->mcstack()->Particle(this->GetMCLabel());
	return pd->GetPdgCode();
}

Int_t MyV0::GetPosTrackPdg() {

	return mHandler->mcstack()->Particle(TMath::Abs(this->GetPosTrack()->GetLabel()))->GetPdgCode();
}

Int_t MyV0::GetNegTrackPdg() {

	return mHandler->mcstack()->Particle(TMath::Abs(this->GetNegTrack()->GetLabel()))->GetPdgCode();
}

Double_t MyV0::CalculateY(Int_t Sp) {

	if (Sp==1) return mAliV0->RapK0Short();
	if (Sp==2) return mAliV0->RapLambda();
	if (Sp==3) return mAliV0->RapLambda();
	return -99.;
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

Int_t MyV0::HasFastSignal() {
  
  // logical or
  Int_t howMany = 0;
  if ((this->GetNegTrack()->GetStatus() & AliESDtrack::kITSrefit)) howMany++;
  if ((this->GetPosTrack()->GetStatus() & AliESDtrack::kITSrefit)) howMany++;
  if (this->GetNegTrack()->HasTOFPID()) howMany++;
  if (this->GetPosTrack()->HasTOFPID()) howMany++;
    
  return howMany++;
}
#endif
