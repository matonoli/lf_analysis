#include "MyParticle.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>


ClassImp(MyParticle)

Int_t MyParticle::GetSign() const
{

  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
  TParticlePDG *ppdg = dbpdg->GetParticle(this->GetPdgCode());
  if (!ppdg)
    return 0;
  return TMath::Nint(ppdg->Charge());
}

#if INPUTFORMAT == 1
Float_t MyParticle::GetY() const
{
	TLorentzVector pVec;
	TDatabasePDG *dbpdg = TDatabasePDG::Instance();
	TParticlePDG *ppdg = dbpdg->GetParticle(this->GetPdgCode());
	if (!ppdg)
		return 0;
	Float_t mass = ppdg->Mass(); 
	pVec.SetPtEtaPhiM(this->GetPt(),this->GetEta(),this->GetPhi(),mass);
	return pVec.Rapidity();	  
}
#endif