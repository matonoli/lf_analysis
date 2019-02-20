#include <TH1.h>
#include <TChain.h>
#include <TBranch.h>

#include "MyAnalysisV0.h"
#include "MyEvent.h"
#include "MyTrack.h"
#include "MyV0.h"
#include "MyHandler.h"

//#include <AliAnalysisPIDV0.h>

using namespace V0consts;

ClassImp(MyAnalysisV0)

MyAnalysisV0::MyAnalysisV0() {

}

Int_t MyAnalysisV0::Init() {

	mFlagMC = mHandler->GetFlagMC();

	//SetCutValues();

	TH1::SetDefaultSumw2();
	CreateHistograms();

	// Link tree variables to arrays -- needs generalisation
	mHandler->Chain()->SetBranchAddress("AnalysisEvent",&mEvent);
	mHandler->Chain()->SetBranchAddress("AnalysisTrack",&bTracks);
	mHandler->Chain()->SetBranchAddress("AnalysisV0Track",&bV0s);
	if (mFlagMC) mHandler->Chain()->SetBranchAddress("AnalysisParticle",&bParticles);

	printf("Analysis %s initiated with flag MC %i \n", this->GetName(), mFlagMC);

	return 0;
}

Int_t MyAnalysisV0::Make(Int_t iEv) {
	printf("Looping in analysis %i \n", iEv);

	// EVENT INFO HISTOGRAMS
	hEventMonitor->Fill(0);
	if (!mEvent) return 1;
	MyEvent event(mEvent);
	hEventMonitor->Fill(1);

	// EVENT SELECTION AND CLASSIFICATION
	enum { multMB, V0M, NCharged };
	enum { sphMB, Jetty, Iso };

	// TRACK LOOP
	Int_t nTracks = bTracks->GetEntriesFast();
	for (int iTr = 0; iTr < nTracks; ++iTr)
	{
		if (!(AliAnalysisPIDTrack*)bTracks->At(iTr)) continue;
		MyTrack t((AliAnalysisPIDTrack*)bTracks->At(iTr));
		//printf("Pt is %f \n", t.GetPt());
	}

	// MC V0 ANALYSIS: PARTICLES LOOP
	if (mFlagMC) {
		Int_t nParticles = bParticles->GetEntriesFast();
	}

	// V0 DATA ANALYSIS: V0 CANDIDATES LOOP
	Int_t nV0s = bV0s->GetEntriesFast();
	for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
		
		hV0Monitor->Fill(0);
		if (!(AliAnalysisPIDV0*)bV0s->At(iV0)) continue;
		MyV0 v0((AliAnalysisPIDV0*)bV0s->At(iV0));

		for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
			if (IsV0(v0,iSp)) ProcessV0(v0,iSp,multMB,sphMB);		}

	}

	return 0;	
}

Bool_t MyAnalysisV0::ProcessV0(MyV0 &v0, Int_t Sp, Int_t Mu, Int_t Sph) {
	
	//printf("v0 pt is %f, args are %i , %i , %i \n", v0.GetPt(), Sp, Mu, Sph);
	hV0Pt[Sp][Mu][Sph]->Fill(v0.GetPt());
	hV0Eta[Sp][Mu][Sph]->Fill(v0.GetEta());
	return true;	
}

Bool_t MyAnalysisV0::IsV0(MyV0 &v0, Int_t Sp) {
	
	if (v0.GetEta() < cuts::V0_ETA[0]) 	return false;
	if (v0.GetEta() > cuts::V0_ETA[1]) 	return false;
	if (v0.GetPt() < cuts::V0_PT[0]) 	return false;
	if (v0.GetPt() > cuts::V0_PT[1]) 	return false;
	if (v0.GetDCAdd() > cuts::V0_DCADD) return false;
	if (v0.GetCPA() < cuts::V0_CPA) 	return false;
	if (v0.GetRadius() < cuts::V0_R[0]) return false;
	if (v0.GetRadius() > cuts::V0_R[1]) return false;

	if (!Sp) return true;
	MyTrack trP(v0.GetPosTrack()); 
	MyTrack trN(v0.GetNegTrack());

	switch (Sp) {
		default : 
			break;
		case 1 	: // K0s
			if (trP.GetNSigmaPionTPC() < cuts::K0_NSIGTPC[0]) 	return false;
			if (trP.GetNSigmaPionTPC() > cuts::K0_NSIGTPC[1]) 	return false;
			if (trN.GetNSigmaPionTPC() < cuts::K0_NSIGTPC[0]) 	return false;
			if (trN.GetNSigmaPionTPC() > cuts::K0_NSIGTPC[1]) 	return false;
			break;
		case 2 	: // L
			if (trP.GetNSigmaProtonTPC() < cuts::K0_NSIGTPC[0]) return false;
			if (trP.GetNSigmaProtonTPC() > cuts::K0_NSIGTPC[1]) return false;
			if (trN.GetNSigmaPionTPC() < cuts::K0_NSIGTPC[0]) 	return false;
			if (trN.GetNSigmaPionTPC() > cuts::K0_NSIGTPC[1]) 	return false;
			break;
		case 3 	: // Lbar
			if (trP.GetNSigmaPionTPC() < cuts::K0_NSIGTPC[0]) 	return false;
			if (trP.GetNSigmaPionTPC() > cuts::K0_NSIGTPC[1]) 	return false;
			if (trN.GetNSigmaProtonTPC() < cuts::K0_NSIGTPC[0]) return false;
			if (trN.GetNSigmaProtonTPC() > cuts::K0_NSIGTPC[1]) return false;
			break;	}

	return true;	
}

Bool_t MyAnalysisV0::CreateHistograms() {

	// MONITORS
	hEventMonitor 			= new TH1D("hEventMonitor","; Step; Entries",10,-0.5,9.5);
	hTrackMonitor 			= new TH1D("hTrackMonitor","; Step; Entries",10,-0.5,9.5);
	hV0Monitor  			= new TH1D("hV0Monitor","; Step; Entries",10,-0.5,9.5);
	hParticleMonitor 		= new TH1D("hParticleMonitor","; Step; Entries",10,-0.5,9.5);

	// EVENT INFO HISTOGRAMS

	// TRACK HISTOGRAMS

	// MC PARTICLE HISTOGRAMS

	// V0 HISTOGRAMS
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
				
		hV0Pt[iSp][iMu][iSph] 		= new TH1D(Form("hV0Pt_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); Entries", 100, 0, 10);
		hV0Eta[iSp][iMu][iSph] 		= new TH1D(Form("hV0Eta_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),
			";V0 #eta; Entries", 100, -1., 1.);

	} } }
		

}

Int_t MyAnalysisV0::Finish() {
	printf("Finishing analysis %s \n",this->GetName());
	return 0;	
}

