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