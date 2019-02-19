#include <TH1.h>
#include <TChain.h>
#include <TBranch.h>

#include "MyAnalysisV0.h"
#include "MyEvent.h"
#include "MyTrack.h"
#include "MyV0.h"
#include "MyHandler.h"

//#include <AliAnalysisPIDV0.h>

ClassImp(MyAnalysisV0)

MyAnalysisV0::MyAnalysisV0() {

}

Int_t MyAnalysisV0::Init() {

	mFlagMC = mHandler->GetFlagMC();

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
	MyEvent event(mEvent);

	// EVENT SELECTION AND CLASSIFICATION
	enum { multMB, V0M, NCharged };
	enum { sphMB, Jetty, Iso };

	// TRACK LOOP
	Int_t nTracks = bTracks->GetEntriesFast();
	for (int iTr = 0; iTr < nTracks; ++iTr)
	{
		Int_t monitor = 0;
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
		
		Int_t monitor = 0;
		hV0Monitor->Fill(monitor++);

		MyV0 v0((AliAnalysisPIDV0*)bV0s->At(iV0));
		for (int iSp = 0; iSp < nBinsSpecies; ++iSp)	{
			ProcessV0(v0,iSp,multMB,sphMB);		}

		hV0Pt[0][0][0]->Fill(v0.GetPt());
	}

	return 0;	
}

Bool_t MyAnalysisV0::ProcessV0(MyV0 v0, Int_t Sp, Int_t Mu, Int_t Sph) {
	
	//printf("args are %i , %i , %i \n", Sp, Mu, Sph);
	return true;	
}

Bool_t MyAnalysisV0::CreateHistograms() {

	// MONITORS
	hEventMonitor 			= new TH1D("hEventMonitor","; Step; Entries",10,-0.5,9.5);
	hTrackMonitor 			= new TH1D("hTrackMonitor","; Step; Entries",10,-0.5,9.5);
	hV0Monitor 				= new TH1D("hV0Monitor","; Step; Entries",10,-0.5,9.5);
	hParticleMonitor 		= new TH1D("hParticleMonitor","; Step; Entries",10,-0.5,9.5);

	// EVENT INFO HISTOGRAMS

	// TRACK HISTOGRAMS

	// MC PARTICLE HISTOGRAMS

	// V0 HISTOGRAMS
	for (int iSp = 0; iSp < nBinsSpecies; ++iSp)	{
	for (int iMu = 0; iMu < nBinsMulti; ++iMu)		{
	for (int iSph = 0; iSph < nBinsSphero; ++iSph)	{
				
		hV0Pt[iSp][iMu][iSph] = new TH1D(Form("hV0Pt_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),
			"V0 Pt; (GeV/#it{c}); Entries", 100, 0, 10);

	} } }
		

}

Int_t MyAnalysisV0::Finish() {
	printf("Finishing analysis \n");
	return 0;	
}

