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
Float_t MyParticle::GetY(Int_t Sp) const
{
	TLorentzVector pVec;
	const Float_t MASSES[] = {0., 0.497614, 1.11568, 1.11568};
	if (Sp < 4 && Sp > -1) pVec.SetPtEtaPhiM(this->GetPt(),this->GetEta(),this->GetPhi(),MASSES[Sp]);
	return pVec.Rapidity();	  
}
#endif