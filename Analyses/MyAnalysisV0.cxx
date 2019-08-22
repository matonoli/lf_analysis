#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TList.h>
#include <TFile.h>
#include <TLegend.h>
#include <TNamed.h>
#include <THashList.h>
#include <TNtuple.h>

#include "MyAnalysisV0.h"
#include "../MyEvent.h"
#include "../MyTrack.h"
#include "../MyParticle.h"
#include "../MyV0.h"
#include "../MyHandler.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooPlot.h"

#include "TransverseSpherocity/TransverseSpherocity.h"

//#include <AliAnalysisPIDV0.h>

using namespace V0consts;
using namespace RooFit;

ClassImp(MyAnalysisV0)

MyAnalysisV0::MyAnalysisV0() {

}

Int_t MyAnalysisV0::Init() {

	TString dfName(this->GetName());
	dfName = Form("%s_%i",dfName.Data(),mHandler->nAnalysis());
	mDirFile = new TDirectoryFile(dfName,dfName,"",mHandler->file());
	mDirFile->cd();

	mFlagMC = mHandler->GetFlagMC();
	mFlagHist = mHandler->GetFlagHist();

	printf("Initialising analysis %s with flag MC %i and Hist %i \n", 
		this->GetName(), mFlagMC, mFlagHist);

	//if (mFlagHist) return 0;

	TH1::SetDefaultSumw2();
	if (mFlagHist)	BorrowHistograms();
	else 			CreateHistograms();

	mList = (TList*)mHandler->directory()->GetList();

	// initialising treebug checker
	bugR = 0;
	bugPt = 0;

	// spherocities
	for (Int_t iType = 0; iType < NTYPE; ++iType)	{
	for (Int_t iMu = 0; iMu < NMULTI-1; ++iMu)	{
		mTS[iType][iMu] = new TransverseSpherocity();
		mTS[iType][iMu]->SetMinMulti(10);
		mTSNorm[iType][iMu] = new TransverseSpherocity();
		mTSNorm[iType][iMu]->SetMinMulti(10);
	}	}

	return 0;
}

Int_t MyAnalysisV0::Make(Int_t iEv) {
	//printf("Looping in analysis %i \n", iEv);

	if (mFlagHist) return 0;

	// EVENT INFO HISTOGRAMS
	hEventMonitor->Fill(0);
	//cout << "mevent is " << mHandler->event() << endl;
	if (!mHandler->event()) return 1;
	MyEvent event(mHandler->event());
	hEventMonitor->Fill(1);

	// EVENT SELECTION
	hEventType->Fill(EVENTTYPES[0],1);	//preES MB
	if (!SelectEvent(event)) return 0;
	hEventType->Fill(EVENTTYPES[1],1);	//postES MB
	hEventMonitor->Fill(2);

	// BUG HOTFIX FOR AURORATREES
	MyV0 bugfix;
	if (mHandler->v0(0)) { bugfix = MyV0(mHandler->v0(0));
	if (TMath::Abs(bugfix.GetRadius()-bugR) < 0.0001
		&& TMath::Abs(bugfix.GetPt()-bugPt) < 0.0001) return 0;
	bugR = bugfix.GetRadius(); bugPt = bugfix.GetPt(); }
	hEventMonitor->Fill(3);

	// EVENT CLASSIFICATION
	enum { multMB, V0M, NCharged, RT };
	enum { sphMB, Jetty, Iso };
	enum { D, RC, MC };
	Int_t isEventCentral = 0;
	hEventV0MCentrality->Fill(event.GetV0MCentrality());
	hEventRefMult->Fill(event.GetRefMult());
	hEventV0MCentvRefMult->Fill(event.GetRefMult(),event.GetV0MCentrality());
	Bool_t isEventFHM = IsCentral(event,V0M);
	Bool_t isEventMHM = IsCentral(event,NCharged);
	if (isEventFHM)	hEventType->Fill(EVENTTYPES[2],1);
	if (isEventMHM)	hEventType->Fill(EVENTTYPES[3],1);
	
	isEventCentral = (isEventFHM || isEventMHM);
	//isEventCentral += IsCentral(event,NCharged);

	for (Int_t iType = 0; iType < NTYPE; ++iType)	{
	for (Int_t iMu = 0; iMu < NMULTI-1; ++iMu)	{
		mTS[iType][iMu]->Reset();
		mTSNorm[iType][iMu]->Reset();
	}	}
	Double_t eventTS = -99.;
	Double_t eventTSRC = -99.;
	Double_t eventTS2[6] = {-99.,-99.,-99.,-99.,-99.,-99.};
	enum { tsD, tsRC, tsMC, tsDNorm, tsRCNorm, tsMCNorm };


	// TRACK LOOP TO CONSTRUCT SPHEROCITY AND FIND LEADING
	hEventMonitor->Fill(4);
	Int_t nTracks = mHandler->tracks()->GetEntriesFast();
	Int_t nParticles = (mFlagMC) ? mHandler->particles()->GetEntriesFast() : 0;

	ptLead = -99.; 
	phiLead = -99.;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));

		//if (!SelectTrack(t)) continue;
		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;
		

		if (SelectTrack(t)) {				// primary dca cut
			if (t.GetPt() > ptLead) {
				ptLead = t.GetPt();
				phiLead = t.GetPhi();	}
		}

		//+add tpc refit 
		if (!t.IskITSrefit()) 		continue;
		if (!t.IsTPCOnlyRefit())	continue;
		//0.15 cut ?

		// for Data calc. s0 from tracks
		if (isEventCentral) {
			mTS[0][0]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[0][0]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
		if (isEventFHM) {
			mTS[0][V0M]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[0][V0M]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
		if (isEventMHM) {
			mTS[0][NCharged]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[0][NCharged]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
	}
	hLeadPhivPt->Fill(ptLead,phiLead);

	// SPHERO CALCULATION ON PARTICLE LEVEL IN MC
	if (mFlagMC && isEventCentral) {
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP));

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.GetPt() > 0.15)	{

				mTS[2][0]->AddTrack(p.GetPx(), p.GetPy());
				mTSNorm[2][0]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventFHM) mTS[2][V0M]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventFHM) mTSNorm[2][V0M]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventMHM) mTS[2][NCharged]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventMHM) mTSNorm[2][NCharged]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

			}	
		}
	}

	// EVENT SPHEROCITY CLASSIFICATION
	Bool_t isEventIso[NMULTI-1] 	= { 0, 0, 0 };
	Bool_t isEventJetty[NMULTI-1] 	= { 0, 0, 0 };
	Bool_t isEventIsoMC[NMULTI-1] 	= { 0, 0, 0 };
	Bool_t isEventJettyMC[NMULTI-1] = { 0, 0, 0 };

	if (isEventCentral) {

		eventTS2[tsD] 		= mTS[D][0]->GetTransverseSpherocityTracks();
		eventTS2[tsRC] 		= mTS[D][0]->GetTransverseSpherocityTracks();	//in MC reconstructed equals data
		eventTS2[tsMC] 		= mTS[MC][0]->GetTransverseSpherocityTracks();
		eventTS2[tsDNorm] 	= mTSNorm[D][0]->GetTransverseSpherocityTracks();
		eventTS2[tsRCNorm]	= mTSNorm[D][0]->GetTransverseSpherocityTracks();
		eventTS2[tsMCNorm]	= mTSNorm[MC][0]->GetTransverseSpherocityTracks();

		//eventTS = (mFlagMC) ? eventTS2[tsMC] : eventTS2[tsD];
		//eventTS = (mFlagMC) ? eventTS2[tsMCNorm] : eventTS2[tsDNorm];
		//eventTS = (mFlagMC) ? eventTS2[tsMCNorm] : eventTS2[tsDNorm];		// in mc study default on particle level
		eventTS = eventTS2[tsDNorm];

		isEventIso[V0M]		= (eventTS > cuts::EV_SPH_ISO[V0M-1] && eventTS < 1.) ;
		isEventJetty[V0M]	= (eventTS < cuts::EV_SPH_JETTY[V0M-1] && eventTS > 0.);
		isEventIsoMC[V0M]	= mFlagMC && (eventTS2[tsMCNorm] > cuts::EV_SPH_ISO[V0M-1] && eventTS2[tsMCNorm] < 1.) ;		// in mc study also particle level 
		isEventJettyMC[V0M]	= mFlagMC && (eventTS2[tsMCNorm] < cuts::EV_SPH_JETTY[V0M-1] && eventTS2[tsMCNorm] > 0.) ;

		isEventIso[NCharged]		= (eventTS > cuts::EV_SPH_ISO[NCharged-1] && eventTS < 1.) ;
		isEventJetty[NCharged]		= (eventTS < cuts::EV_SPH_JETTY[NCharged-1] && eventTS > 0.);
		isEventIsoMC[NCharged]		= mFlagMC && (eventTS2[tsMCNorm] > cuts::EV_SPH_ISO[NCharged-1] && eventTS2[tsMCNorm] < 1.) ;		// in mc study also particle level 
		isEventJettyMC[NCharged]	= mFlagMC && (eventTS2[tsMCNorm] < cuts::EV_SPH_JETTY[NCharged-1] && eventTS2[tsMCNorm] > 0.) ;
		
		if (isEventFHM) hEventSpherocityV0M->Fill(eventTS);
		if (isEventMHM) hEventSpherocityNCharged->Fill(eventTS);
		hEventTSMCvRC->Fill(eventTS2[tsMC],eventTS2[tsRC]);
		hEventTSNormMCvRC->Fill(eventTS2[tsMCNorm],eventTS2[tsRCNorm]);
		hEventTSMCvNorm->Fill(eventTS2[tsMC],eventTS2[tsMCNorm]);
		hEventTSRCvNorm->Fill(eventTS2[tsRC],eventTS2[tsRCNorm]);
	}

	if (isEventFHM && isEventIso[V0M])			hEventType->Fill(EVENTTYPES[6],1);
	if (isEventFHM && isEventJetty[V0M])		hEventType->Fill(EVENTTYPES[7],1);
	if (isEventMHM && isEventIso[NCharged])		hEventType->Fill(EVENTTYPES[8],1);
	if (isEventMHM && isEventJetty[NCharged])	hEventType->Fill(EVENTTYPES[9],1);
	if (isEventFHM && isEventIsoMC[V0M])		hEventType->Fill(EVENTTYPES[16],1);
	if (isEventFHM && isEventJettyMC[V0M])		hEventType->Fill(EVENTTYPES[17],1);
	if (isEventMHM && isEventIsoMC[NCharged])	hEventType->Fill(EVENTTYPES[18],1);
	if (isEventMHM && isEventJettyMC[NCharged])	hEventType->Fill(EVENTTYPES[19],1);

	// TRACK SPHEROCITY STUDY LOOP
	nChTrans = 0;
	Int_t nChTrans0 = 0;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));

		if (!SelectTrack(t)) continue;

		ProcessTrack(t,D,multMB,sphMB);
		if (mFlagMC) ProcessTrack(t,RC,multMB,sphMB);

		if (isEventFHM) {
			ProcessTrack(t,D,V0M,sphMB);
			if (mFlagMC)	ProcessTrack(t,RC,V0M,sphMB);
			if (isEventJetty[V0M])		ProcessTrack(t,D,V0M,Jetty);
			if (isEventIso[V0M])		ProcessTrack(t,D,V0M,Iso);			//
			if (isEventJettyMC[V0M])	ProcessTrack(t,RC,V0M,Jetty);		// study tracks also for true sphero, stored under RC
			if (isEventIsoMC[V0M])		ProcessTrack(t,RC,V0M,Iso);		} 

		if (isEventMHM) {
			ProcessTrack(t,D,NCharged,sphMB);
			ProcessTrack(t,RC,NCharged,sphMB);
			if (isEventJetty[NCharged])		ProcessTrack(t,D,NCharged,Jetty);
			if (isEventIso[NCharged])		ProcessTrack(t,D,NCharged,Iso);
			if (isEventJettyMC[NCharged])	ProcessTrack(t,RC,NCharged,Jetty);
			if (isEventIsoMC[NCharged])		ProcessTrack(t,RC,NCharged,Iso);		}

		/*for (int iMu = 0; iMu < isEventFHM+isEventMHM+1; ++iMu) {
		for (int iSph = 0; iSph < isEventJetty+isEventIso+1; ++iSph) {
			ProcessTrack(t,D,iMu+iMu*(!isEventFHM),iSph+iSph*(!isEventJetty));
		}	}	*/
	}

	// RT NCH CALCULATION
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));
		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;
		
		// RT DETERMINATION
		if (!IsTrans(t.GetPhi(),phiLead))	continue;
		if (!t.IskITSrefit()) continue;
		if (!t.IsTPCOnlyRefit()) continue;
		if (t.GetPt()<0.15) continue;
		
		hNchvLeadPt->Fill(ptLead);
		nChTrans0++;
		if (ptLead < 5.) continue;
		if (ptLead > 40.) continue;
		nChTrans++;
	}

	hNchvLeadPt2->Fill(ptLead,nChTrans0);
	eventRt = 0;
	if (ptLead>5. && ptLead < 40.) {
		hNchTrans->Fill(nChTrans);
		eventRt = (double)nChTrans/RT_DEN;
		hRt->Fill(eventRt);		}

	Bool_t isEventRT = false;
	enum { rt01, rt12, rt23, rt34, rt45, rtsizeof };
	Bool_t isRT[rtsizeof] = {false, false, false, false, false};
	if (ptLead>5. && ptLead<40.) {
		isEventRT = true;
		for (int iRt = 0; iRt < rtsizeof; ++iRt)	{
			isRT[iRt] = (eventRt > (Double_t)iRt 
				&& eventRt < (Double_t)(iRt+1));	}
	}

	if (isEventRT) hEventType->Fill(EVENTTYPES[10],1);
	for (int iRt = 0; iRt < rtsizeof; ++iRt)	{
		if (isRT[iRt]) hEventType->Fill(EVENTTYPES[11+iRt],1); }

	// <pT> vs RT
	hLeadPtvNchTrans0->Fill(nChTrans0,ptLead);
	if (isEventRT) hLeadPtvNchTrans->Fill(nChTrans,ptLead);
	if (isEventRT && isEventFHM) hNchTransvSpherocityV0M->Fill(eventTS,nChTrans);
	if (isEventRT && isEventMHM) hNchTransvSpherocityNCharged->Fill(eventTS,nChTrans);

	// STUDY DPHI DISTRIBUTION IN RT EVENTS
	if (isEventRT) for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));
		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;
		
		if (!t.IskITSrefit()) continue;
		if (!t.IsTPCOnlyRefit()) continue;
		if (t.GetPt()<0.15) continue;
		
		if (TMath::Abs(t.GetPt()-ptLead)>1E-5) hTrackDPhivNchTrans->Fill(nChTrans,mHandler->DeltaPhi(phiLead,t.GetPhi())); //TMath::Abs(phiLead-t.GetPhi()));

	}



	//////////////////////////////////////////////////////////////////
	// MC RT CALCULATION
	//////////////////////////////////////////////////////////////////
	Double_t ptLeadMC = -99., phiLeadMC = -99.;
	if (mFlagMC) {
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP));

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0)	{

				if (p.GetPt() > ptLeadMC) {
					ptLeadMC = p.GetPt();
					phiLeadMC = p.GetPhi();
				}
			}	
		}

		// RT
		Int_t nChTransMC = 0;
		Double_t eventRtMC = 0;
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP));

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.GetPt() > 0.15)	{

				if (!IsTrans(p.GetPhi(),phiLeadMC)) continue;
				if (ptLeadMC<5.0) continue;
				if (ptLeadMC>40.0) continue;

				nChTransMC++;
			}
		}

		if ( (ptLeadMC > 5. && ptLeadMC < 40.) && (ptLead > 5. && ptLead < 40.) )	{   // or ||
			hNchTransMC->Fill(nChTransMC);
			hNchTransRCvMC->Fill(nChTransMC,nChTrans);
			eventRtMC = (double)nChTransMC/RT_DEN_MC;
			hRtMC->Fill(eventRtMC);
			hRtRCvMC->Fill(eventRtMC,eventRt);		}

		if (ptLeadMC > 5. && ptLeadMC < 40.) for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP));

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetPt() > 0.15)	{

				if (p.GetSign() != 0 && TMath::Abs(p.GetPt()-ptLeadMC)>1E-5) hParticleDPhivNchTrans->Fill(nChTransMC,mHandler->DeltaPhi(phiLeadMC,p.GetPhi()));
				Int_t region = WhatRegion(p.GetPhi(),phiLeadMC);
				if (TMath::Abs(p.GetPdgCode())==2212) 	hProtonNchTransvPt[region]->Fill(p.GetPt(),nChTransMC);
				if (TMath::Abs(p.GetPdgCode())==211) 	hPionNchTransvPt[region]->Fill(p.GetPt(),nChTransMC);
				if (TMath::Abs(p.GetPdgCode())==3122) 	hLambdaNchTransvPt[region]->Fill(p.GetPt(),nChTransMC);
				if (TMath::Abs(p.GetPdgCode())==310) 	hK0sNchTransvPt[region]->Fill(p.GetPt(),nChTransMC);
				// b/m study in MC

			}
		}
	}
	/////////////////////////////////////////////////////////////////


	// MC V0 ANALYSIS: PARTICLES LOOP
	hEventMonitor->Fill(5);
	std::vector<Int_t> PartLabels;
	std::vector<Int_t> PartIds; 
	if (mFlagMC) {
		//Int_t nParticles = mHandler->particles()->GetEntriesFast();
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP));


			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0)	{

				hTrackPt[MC][multMB][sphMB]->Fill(p.GetPt());

				if (isEventFHM) {
					hTrackPt[MC][V0M][sphMB]->Fill(p.GetPt());
					if (isEventJettyMC[V0M])	hTrackPt[MC][V0M][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC[V0M])		hTrackPt[MC][V0M][Iso]->Fill(p.GetPt()); }

				if (isEventMHM) {
					hTrackPt[MC][NCharged][sphMB]->Fill(p.GetPt());
					if (isEventJettyMC[NCharged])	hTrackPt[MC][NCharged][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC[NCharged])		hTrackPt[MC][NCharged][Iso]->Fill(p.GetPt()); }

				if (isEventRT && IsTrans(p.GetPhi(),phiLead)) {
					hTrackPt[MC][RT][sphMB]->Fill(p.GetPt());
					for (int iRt = 0; iRt < rtsizeof; ++iRt) {
						if (isRT[iRt])	hTrackPt[MC][RT][3+iRt]->Fill(p.GetPt());	}
				}

			}
				
				/*for (int iMu = 0; iMu < isEventFHM+isEventMHM+1; ++iMu) {
				for (int iSph = 0; iSph < isEventJetty+isEventIso+1; ++iSph) {
					hTrackPt[MC][iMu+iMu*(!isEventFHM)][iSph+iSph*(!isEventJetty)]->Fill(p.GetPt());
				}	}*/

			if (!SelectParticle(p)) continue;
			PartLabels.push_back(p.GetLabel());
			PartIds.push_back(iP);

			for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
				if (p.GetPdgCode() == PDG_IDS[iSp])		{

					hV0Pt[iSp][MC][multMB][sphMB]->Fill(p.GetPt());
					tV0PtMCMB[iSp]->Fill(p.GetPt());

					if (isEventFHM) {
						hV0Pt[iSp][MC][V0M][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC[V0M])	hV0Pt[iSp][MC][V0M][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC[V0M])		hV0Pt[iSp][MC][V0M][Iso]->Fill(p.GetPt());		}

					if (isEventMHM) {
						hV0Pt[iSp][MC][NCharged][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC[NCharged])	hV0Pt[iSp][MC][NCharged][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC[NCharged])		hV0Pt[iSp][MC][NCharged][Iso]->Fill(p.GetPt());		}

					if (isEventRT)	{
							Int_t region = WhatRegion(p.GetPhi(),phiLead);

							if (iSp>0) tV0PtMCRt[iSp][region]->Fill(p.GetPt());
						}

					if (isEventRT && IsTrans(p.GetPhi(),phiLead)) {
						hV0Pt[iSp][MC][RT][sphMB]->Fill(p.GetPt());
						for (int iRt = 0; iRt < rtsizeof; ++iRt) {
							if (isRT[iRt])	hV0Pt[iSp][MC][RT][3+iRt]->Fill(p.GetPt());	}
					}

					/*for (int iMu = 0; iMu < isEventFHM+isEventMHM+1; ++iMu) {
					for (int iSph = 0; iSph < isEventJetty+isEventIso+1; ++iSph) {
						hV0Pt[iSp][MC][iMu+iMu*(!isEventFHM)][iSph+iSph*(!isEventJetty)]->Fill(p.GetPt());
					}	}*/
				}
			}

			for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
				//if (p.GetPdgCode() == PDG_IDS[iSp]){
					//hV0Pt[iSp][MC][0][0]->Fill(p.GetPt());	}
			}
		}

	}

	// V0 DATA ANALYSIS: V0 CANDIDATES LOOP
	hEventMonitor->Fill(6);
	Int_t nV0s = mHandler->v0s()->GetEntriesFast();
	for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
		
		hV0Monitor->Fill(0);
		if (!mHandler->v0(iV0)) continue;
		MyV0 v0(mHandler->v0(iV0));


		if (mFlagMC) {
			MyParticle v0mc;
			Bool_t MCfound = false;
			for (unsigned int iP = 0; iP < PartLabels.size(); ++iP)	{
				if (PartLabels[iP] == v0.GetMCLabel()) {
					v0mc = MyParticle(mHandler->particle(PartIds[iP]));
					if (v0mc.GetPdgCode() != v0mc.GetMotherPdgCode() 
							&& TMath::Abs(v0mc.GetMotherPdgCode()) != 311
							&& TMath::Abs(v0mc.GetMotherPdgCode()) > 10)	{
						v0mc.SetIsPrimary(0);	}
					MCfound = true;
					break;	}
			}

			if (MCfound) {
				for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
					if (IsV0(v0,iSp,RC)) {
						if (v0mc.GetIsPrimary()) hV0Feeddown[iSp]->Fill(v0.GetPt());
						else hV0FeeddownPDG[iSp]->Fill(v0mc.GetMotherPdgCode()); //printf("code is %i and %i \n", v0mc.GetPdgCode(), v0mc.GetMotherPdgCode());
						
						ProcessV0toHist(v0,iSp,RC,multMB,sphMB);		
						ProcessV0toTree(v0,iSp,RC,0);		

						if (isEventFHM) {
							ProcessV0toHist(v0,iSp,RC,V0M,sphMB);
							if (isEventJettyMC[V0M])	ProcessV0toHist(v0,iSp,RC,V0M,Jetty);			// study RC spectra for true spher. ev classification
							if (isEventIsoMC[V0M])		ProcessV0toHist(v0,iSp,RC,V0M,Iso);
						}
						
						if (isEventMHM) {
							ProcessV0toHist(v0,iSp,RC,NCharged,sphMB);
							if (isEventJettyMC[NCharged])	ProcessV0toHist(v0,iSp,RC,NCharged,Jetty);
							if (isEventIsoMC[NCharged])		ProcessV0toHist(v0,iSp,RC,NCharged,Iso);
						}

						if (isEventRT)	{
							Int_t region = WhatRegion(v0.GetPhi(),phiLead);
							if (iSp>0) ProcessV0toTree(v0,iSp,RC,RT+region);
						}

						if (isEventRT && IsTrans(v0.GetPhi(),phiLead)) {
							ProcessV0toHist(v0,iSp,RC,RT,sphMB);
							for (int iRt = 0; iRt < rtsizeof; ++iRt) {
								if (isRT[iRt])	ProcessV0toHist(v0,iSp,RC,RT,3+iRt);	}
						}
					}
				}
			}
		}

		//for (int iMu = 0; iMu < isEventFHM+isEventMHM+isEventRT+1; ++iMu) {
		//for (int iSph = 0; iSph < isEventJetty+isEventIso+1; ++iSph) {
		for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
			if (IsV0(v0,iSp,D)) {

				ProcessV0toHist(v0,iSp,D,multMB,sphMB);

				if (isEventFHM) {
					ProcessV0toHist(v0,iSp,D,V0M,sphMB);
					if (isEventJetty[V0M])	ProcessV0toHist(v0,iSp,D,V0M,Jetty);
					if (isEventIso[V0M])	ProcessV0toHist(v0,iSp,D,V0M,Iso); }

				if (isEventMHM) {
					ProcessV0toHist(v0,iSp,D,NCharged,sphMB);
					if (isEventJetty[NCharged])	ProcessV0toHist(v0,iSp,D,NCharged,Jetty);
					if (isEventIso[NCharged])	ProcessV0toHist(v0,iSp,D,NCharged,Iso); }
				
				if (isEventRT) {
					
					Int_t region = WhatRegion(v0.GetPhi(),phiLead);
					if (iSp>0) ProcessV0toTree(v0,iSp,D,RT+region);

					ProcessV0toHist(v0,iSp,D,RT+region,sphMB);
					for (int iRt = 0; iRt < rtsizeof; ++iRt) {
						if (isRT[iRt])	ProcessV0toHist(v0,iSp,D,RT+region,3+iRt);	}
				}

				//	ProcessV0toHist(v0,iSp,D,iMu+iMu*(!isEventFHM)*(!isEventRT),iSph+iSph*(!isEventJetty));
			}
		}		

	}

	hEventMonitor->Fill(7);
	return 0;	
}

Bool_t MyAnalysisV0::SelectEvent(MyEvent &ev) {

	//if (!ev.IsGoodAliEvent())			return false;
	//if(bV0s->GetEntriesFast() < 1)		return false;
	//if(bTracks->GetEntriesFast() < 1)		return false;
	hEventCuts->Fill(0);
	if (!ev.AcceptVertex())				return false;
	hEventCuts->Fill(1);
	if (ev.IsPileupFromSPD())			return false;
	hEventCuts->Fill(2);
	//if (!ev.CheckFlag())				return false;
	//if (!ev.IsCollisionCandidate())		return false;
	//if (ev.GetCentralityQuality() != 0)	return false;

	return true;
}


Bool_t MyAnalysisV0::IsCentral(MyEvent &ev, Int_t Mu) {

	switch (Mu) {
		default: 
			break;
		case 1: 
			if (ev.GetV0MCentrality() < 10. 
				&& ev.GetRefMult() > 9.9)		return true;
			break;

		case 2: 
			if (ev.GetRefMult() < cuts::EV_MHM[1] 
				&& ev.GetRefMult() > cuts::EV_MHM[0])		return true;
			break;
	}

	return false;
}

Bool_t MyAnalysisV0::ProcessV0toHist(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph) {
	
	//printf("v0 pt is %f, args are %i , %i , %i \n", v0.GetPt(), Sp, Mu, Sph);

	hV0Pt[Sp][Type][Mu][Sph]->Fill(v0.GetPt());
	hV0Eta[Sp][Type][Mu][Sph]->Fill(v0.GetEta());
	Double_t v0mass[] = {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	hV0IMvPt[Sp][Type][Mu][Sph]->Fill(v0.GetPt(),v0mass[Sp]);

	if (Mu>2 && Type==0 && Sph == 0 && Sp==1) {
		hV0DPhivNchTrans->Fill(nChTrans,mHandler->DeltaPhi(phiLead,v0.GetPhi()));
	}

	//if (Sp>0 && Type==0 && Mu==3 && Sph==0) {
	//	tV0massRt[Sp][0][0]->Fill(v0mass[Sp],v0.GetPt(),eventRt);	}

	return true;	
}

Bool_t MyAnalysisV0::ProcessV0toTree(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu) {
	
	Double_t v0mass[] = {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	
	if (Mu > 2) {
		Int_t Reg = Mu-3;
		tV0massRt[Sp][Type][Reg]->Fill(v0mass[Sp],v0.GetPt(),nChTrans);
	}

	if (Mu == 0 && Type==1) {
		tV0massRCMB[Sp]->Fill(v0mass[Sp],v0.GetPt());	
	}


	return true;	
}

Bool_t MyAnalysisV0::ProcessTrack(MyTrack &t, Int_t Type, Int_t Mu, Int_t Sph) {

	hTrackPt[Type][Mu][Sph]->Fill(t.GetPt());
	hTrackEtavPhi[Type][Mu][Sph]->Fill(t.GetPhi(),t.GetEta());

}

Bool_t MyAnalysisV0::IsV0(MyV0 &v0, Int_t Sp, Int_t Type) {
	
	//printf("is primary %i \n", v0.IsMCPrimary());

	if (Type>0) {
		if (v0.GetMCPdgCode() != PDG_IDS[Sp])	return false; }
		
		//if (!v0.IsMCPrimary())					return false; } // always 1 ?

	//if (Type==0 && Sp==2) if (TMath::Abs(v0.GetIML()) > 0.01) return false;

	//if (v0.GetEta() < cuts::V0_ETA[0]) 	return false;
	//if (v0.GetEta() > cuts::V0_ETA[1]) 	return false;
	if (v0.CalculateY(Sp) < cuts::V0_Y[0]) 	return false;
	if (v0.CalculateY(Sp) > cuts::V0_Y[1]) 	return false;
	if (v0.GetPt() < cuts::V0_PT[0]) 	return false;
	if (v0.GetPt() > cuts::V0_PT[1]) 	return false;
	if (v0.GetDCAdd() > cuts::V0_DCADD) return false;
	if (v0.GetCPA() < cuts::V0_CPA) 	return false;
	if (v0.GetRadius() < cuts::V0_R[0]) return false;
	if (v0.GetRadius() > cuts::V0_R[1]) return false;

	//if (!Sp) return true;
	MyTrack trP(v0.GetPosTrack()); 
	MyTrack trN(v0.GetNegTrack());
	if (!SelectV0Daughter(trP)) return false;
	if (!SelectV0Daughter(trN)) return false;

	switch (Sp) {
		default : 
			break;
		case 1 	: // K0s
			if (*(v0.CalculateAP()+1) < cuts::K0S_AP*TMath::Abs(*(v0.CalculateAP()+0))) return false;
			if (trP.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0]) 	return false;
			if (trP.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1]) 	return false;
			if (trN.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0]) 	return false;
			if (trN.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1]) 	return false;
			break;
		case 2 	: // L
			if (v0.GetPt() < cuts::L_PT[0]) 	return false;
			if (v0.GetPt() > cuts::L_PT[1]) 	return false;
			if (trP.GetNSigmaProtonTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			if (trP.GetNSigmaProtonTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			if (trN.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			if (trN.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			break;
		case 3 	: // Lbar
			if (v0.GetPt() < cuts::L_PT[0]) 	return false;
			if (v0.GetPt() > cuts::L_PT[1]) 	return false;
			if (trP.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			if (trP.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			if (trN.GetNSigmaProtonTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			if (trN.GetNSigmaProtonTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			break;	}

	return true;	
}

Bool_t MyAnalysisV0::IsTrans(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	if (TMath::Abs(dphi) < TMath::Pi()/3.) 		return false;
	if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return false;
	
	return true;
}

Int_t MyAnalysisV0::WhatRegion(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	
	if (TMath::Abs(dphi) < TMath::Pi()/3.)			return 1;
	else if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return 2;
	else return 0;
}

Bool_t MyAnalysisV0::SelectV0Daughter(MyTrack &tr) {

	if (tr.GetEta() < cuts::K0S_D_ETA[0])	return false;
	if (tr.GetEta() > cuts::K0S_D_ETA[1])	return false;
	if (TMath::Abs(tr.GetDCApvXY()) < cuts::K0S_D_DCAPVXY)	return false;

	return true;
}

Bool_t MyAnalysisV0::SelectParticle(MyParticle &p) {

	//if (p.GetEta() < cuts::V0_ETA[0]) 		return false;
	//if (p.GetEta() > cuts::V0_ETA[1]) 		return false;
	if (p.GetY() < cuts::V0_Y[0]) 		return false;
	if (p.GetY() > cuts::V0_Y[1]) 		return false;
	if (p.GetPdgCode() != PDG_IDS[1]
		&& p.GetPdgCode() != PDG_IDS[2]
		&& p.GetPdgCode() != PDG_IDS[3])	return false;
	
	//printf("code is %i and %i \n", p.GetPdgCode(), p.GetMotherPdgCode());

	return true;
}

Bool_t MyAnalysisV0::SelectTrack(MyTrack &tr) {

	if (tr.GetEta() < cuts::V0_ETA[0]) 		return false;
	if (tr.GetEta() > cuts::V0_ETA[1]) 		return false;
	if (TMath::Abs(tr.GetDCApvXY()) > 
		cuts::TR_PRIMARY_PAR[0] + 
		cuts::TR_PRIMARY_PAR[1]/TMath::Power(tr.GetPt(),cuts::TR_PRIMARY_PAR[2])) return false;

	if (!tr.IsITSTPC2011())					return false;


	return true;
}

Bool_t MyAnalysisV0::CreateHistograms() {

	// MONITORS
	hEventMonitor 			= new TH1D("hEventMonitor","; Step; Entries",10,-0.5,9.5);
	hTrackMonitor 			= new TH1D("hTrackMonitor","; Step; Entries",10,-0.5,9.5);
	hV0Monitor  			= new TH1D("hV0Monitor","; Step; Entries",10,-0.5,9.5);
	hParticleMonitor 		= new TH1D("hParticleMonitor","; Step; Entries",10,-0.5,9.5);

	// EVENT INFO HISTOGRAMS
	hEventCuts	 			= new TH1D("hEventCuts","; Step; Entries",10,-0.5,9.5);
	hEventV0MCentrality		= new TH1D("hEventV0MCentrality","; V0M Centrality; Entries",300,0,150);
	hEventRefMult			= new TH1D("hEventRefMult","; Reference multiplicity; Entries",150,0,150);
	hEventV0MCentvRefMult	= new TH2D("hEventV0MCentvRefMult","; Reference multiplicity; V0M Centrality; Entries"
		,150,0,150,300,0,150);

	hEventSpherocityV0M			= new TH1D("hEventSpherocityV0M","; S_{O}; Entries",400,-0.1,1.1);
	hEventSpherocityNCharged	= new TH1D("hEventSpherocityNCharged","; S_{O}; Entries",400,-0.1,1.1);
	hEventType					= new TH1D("hEventType",";; Counts", NEVENTTYPES, 0, NEVENTTYPES);
	for (int iBin = 1; iBin <= NEVENTTYPES; iBin++) {
		hEventType->GetXaxis()->SetBinLabel(iBin,EVENTTYPES[iBin-1]); }

	hEventTSMCvRC 			= new TH2D("hEventTSMCvRC",";S_{O} with MC particles; S_{O} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSNormMCvRC		= new TH2D("hEventTSNormMCvRC",";S_{O}^{p_{T}=1} with MC particles; S_{O}^{p_{T}=1} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSMCvNorm			= new TH2D("hEventTSMCvNorm",";S_{O} with MC particles; S_{O}^{p_{T}=1} with MC particles", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSRCvNorm			= new TH2D("hEventTSRCvNorm",";S_{O} with RC tracks; S_{O}^{p_{T}=1} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);

	hLeadPhivPt				= new TH2D("hLeadPhivPt","; p_{T} (GeV/#it{c}); #phi", 200, 0., 30., 400, -0.2, 6.4);
	hNchvLeadPt				= new TH1D("hNchvLeadPt","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 200, 0., 30.);
	hNchvLeadPt2			= new TH2D("hNchvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 90, 0., 30.,50,-0.5,49.5);
	hNchTrans				= new TH1D("hNchTrans","; N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransMC 			= new TH1D("hNchTransMC","; MC N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransRCvMC			= new TH2D("hNchTransRCvMC", ";MC N_ch [trans.]; RC N_ch [trans.]",50,-0.5,49.5,50,-0.5,49.5);
	hRt						= new TH1D("hRt","; R_{T}; Entries",4000, -0.02, 5.02);//4.975);
	hRtMC					= new TH1D("hRtMC","; MC R_{T}; Entries",4000, -0.02, 5.02);//4.975);
	hRtRCvMC 				= new TH2D("hRtRCMC","; MC R_{T}; R_{T}", 200, -0.02, 5.02, 200, -0.02, 5.02);
 	hLeadPtvNchTrans0		= new TH2D("hLeadPtvNchTrans0","; N_{ch}^{trans}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hLeadPtvNchTrans		= new TH2D("hLeadPtvNchTrans","; N_{ch}^{trans}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);

 	hNchTransvSpherocityV0M			= new TH2D("hNchTransvSpherocityV0M",";S_{O};N_{ch}^{trans}",400,-0.1,1.1,50,-0.5,49.5);
 	hNchTransvSpherocityNCharged	= new TH2D("hNchTransvSpherocityNCharged",";S_{O};N_{ch}^{trans}",400,-0.1,1.1,50,-0.5,49.5);
 	hV0DPhivNchTrans				= new TH2D("hV0DPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hTrackDPhivNchTrans				= new TH2D("hTrackDPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hParticleDPhivNchTrans			= new TH2D("hParticleDPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);

	// TRACK HISTOGRAMS
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hTrackPt[iType][iMu][iSph]			= new TH1D(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";track p_{T} (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hTrackEtavPhi[iType][iMu][iSph]		= new TH2D(Form("hTrackEtavPhi_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";track #phi; track #eta", 400, -0.2, 6.4, 400, -1., 1.);		

	} } } 

	// MC PARTICLE HISTOGRAMS
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
		hProtonNchTransvPt[iReg]		= new TH2D(Form("hProtonNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
		hPionNchTransvPt[iReg]			= new TH2D(Form("hPionNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
		hLambdaNchTransvPt[iReg]		= new TH2D(Form("hLambdaNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
		hK0sNchTransvPt[iReg]			= new TH2D(Form("hK0sNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
	}

	// V0 HISTOGRAMS
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0Pt[iSp][iType][iMu][iSph]			= new TH1D(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
			//";V0 Pt (GeV/#it{c}); Entries",								400, 0, 20);
		hV0Eta[iSp][iType][iMu][iSph]			= new TH1D(Form("hV0Eta_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 #eta; Entries", 										200, -1., 1.);
		hV0IMvPt[iSp][iType][iMu][iSph]		= new TH2D(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.1, 0.1);

		//hV0PtFit[iSp][iType][iMu][iSph]		= new TH1D(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
		//	";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		

	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0Feeddown[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Feeddown_%s",SPECIES[iSp]));
		hV0FeeddownPDG[iSp] =	new TH1D(Form("hV0FeeddownPDG_%s",SPECIES[iSp]),";PDG ID;Entries",20000,-10000,10000);
	}


	// V0 NTUPLES
	Int_t nType = (mFlagMC) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{		
	
		tV0massRt[iSp][iType][iReg] = new TNtuple(Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),"v0 rt mass tree","MassDT:lPt:lNchTrans");

	}	}	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
		tV0PtMCMB[iSp]		= new TNtuple(Form("tV0PtMCMB_%s",SPECIES[iSp]),"v0 MC MB pt tree","lPt");
		tV0massRCMB[iSp]	= new TNtuple(Form("tV0massRCMB_%s",SPECIES[iSp]),"v0 RC MB mass tree","MassDT:lPt");
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{	
			tV0PtMCRt[iSp][iReg]	= new TNtuple(Form("tV0PtMCMB_%s_%s",SPECIES[iSp],REGIONS[iReg]),"v0 MC Rt pt tree","lPt");
		}
	}


}

Bool_t MyAnalysisV0::BorrowHistograms() {

	mDirFile = (TDirectoryFile*)mHandler->filehist()->Get("MyAnalysisV0_0");	// this was smart of me i think
	//mDirFile->ls();

	// MONITORS
	hEventMonitor 				= (TH1D*)mDirFile->Get("hEventMonitor");
	hTrackMonitor 				= (TH1D*)mDirFile->Get("hTrackMonitor");
	hV0Monitor  				= (TH1D*)mDirFile->Get("hV0Monitor");
	hParticleMonitor 			= (TH1D*)mDirFile->Get("hParticleMonitor");
	hEventType					= (TH1D*)mDirFile->Get("hEventType");

	// EVENT INFO HISTOGRAMS
	hEventV0MCentrality			= (TH1D*)mDirFile->Get("hEventV0MCentrality");
	hEventRefMult				= (TH1D*)mDirFile->Get("hEventRefMult");
	hEventV0MCentvRefMult		= (TH2D*)mDirFile->Get("hEventV0MCentvRefMult");

	hEventSpherocityV0M			= (TH1D*)mDirFile->Get("hEventSpherocityV0M");
	hEventSpherocityNCharged	= (TH1D*)mDirFile->Get("hEventSpherocityNCharged");
	hEventTSMCvRC				= (TH2D*)mDirFile->Get("hEventTSMCvRC");
	hEventTSNormMCvRC			= (TH2D*)mDirFile->Get("hEventTSNormMCvRC");
	hEventTSMCvNorm				= (TH2D*)mDirFile->Get("hEventTSMCvNorm");
	hEventTSRCvNorm				= (TH2D*)mDirFile->Get("hEventTSRCvNorm");

	hLeadPhivPt					= (TH2D*)mDirFile->Get("hLeadPhivPt");
	hNchvLeadPt					= (TH1D*)mDirFile->Get("hNchvLeadPt");
	hNchvLeadPt2				= (TH2D*)mDirFile->Get("hNchvLeadPt2");
	hNchTrans					= (TH1D*)mDirFile->Get("hNchTrans");
	hRt							= (TH1D*)mDirFile->Get("hRt");
	//hLeadPtvRt					= (TH2D*)mDirFile->Get("hLeadPtvRt");

	// TRACK HISTOGRAMS
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hTrackPt[iType][iMu][iSph]	= (TH1D*)mDirFile->Get(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hTrackEtavPhi[iType][iMu][iSph]	= (TH2D*)mDirFile->Get(Form("hTrackEtavPhi_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));

	} } } 

	// MC PARTICLE HISTOGRAMS

	// V0 HISTOGRAMS
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0Pt[iSp][iType][iMu][iSph]			= (TH1D*)mDirFile->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0Eta[iSp][iType][iMu][iSph]			= (TH1D*)mDirFile->Get(Form("hV0Eta_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0IMvPt[iSp][iType][iMu][iSph]			= (TH2D*)mDirFile->Get(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		// actually perhaps only histos which need to be processed within this analysis need to get fetched


		//hV0PtFit[iSp][iType][iMu][iSph]		= new TH1D(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
		//	";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		
	} } } }	

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0Feeddown[iSp] 	= (TH1D*)mDirFile->Get(Form("hV0Feeddown_%s",SPECIES[iSp]));
		hV0FeeddownPDG[iSp] = (TH1D*)mDirFile->Get(Form("hV0FeeddownPDG_%s",SPECIES[iSp]));
	}

	Int_t nType = (mFlagMC) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{		
	
		tV0massRt[iSp][iType][iReg] = (TNtuple*)mDirFile->Get(Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));

	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
		tV0PtMCMB[iSp]		= (TNtuple*)mDirFile->Get(Form("tV0PtMCMB_%s",SPECIES[iSp]));
		tV0massRCMB[iSp]	= (TNtuple*)mDirFile->Get(Form("tV0massRCMB_%s",SPECIES[iSp]));
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{	
			tV0PtMCRt[iSp][iReg]	= (TNtuple*)mDirFile->Get(Form("tV0PtMCMB_%s_%s",SPECIES[iSp],REGIONS[iReg]));
		}
	}		

}

Int_t MyAnalysisV0::Finish() {
	
	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();


	TH1D* hLeadPt = hLeadPhivPt->ProjectionX();
	hNchvLeadPt->Divide(hLeadPt);
	

	//MAKE RT BINNING
	hRt2 = (TH1D*)hNchTrans->Clone("hRt2");
	Double_t rt_den = hNchTrans->GetMean();
	Int_t nbins = hNchTrans->GetXaxis()->GetNbins();
	Double_t rtbins[nbins+1];
	//cout << "nb " << nbins << " rt " << rtbins << endl;
	for (int iBin = 0; iBin < nbins+1; ++iBin)	{
		//cout << "a bin " << hNchTrans->GetBinLowEdge(iBin+1) << endl;
		if (rt_den>0) rtbins[iBin] = (double)hNchTrans->GetBinLowEdge(iBin+1)/rt_den;
		//cout << "b bin " << rtbins[iBin] << endl;
	}

	hRt2->SetBins(nbins,rtbins);

	//if (mFlagMC) DoEfficiency();
	if (mFlagMC) DoEfficiencyFromTrees();
	if (mFlagMC) DoLambdaFeeddown();

	return 0;	
}

void MyAnalysisV0::DoEfficiency() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0Efficiency[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Efficiency_%s",SPECIES[iSp]));

		hV0Efficiency[iSp]->SetTitle("; V0 pT (GeV/#it{c}); Efficiency");
		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Efficiency[iSp]->Divide(hV0Pt[iSp][2][0][0]);
		//hV0Efficiency[iSp]->Write();
	}

}

void MyAnalysisV0::DoEfficiencyFromTrees() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		
		// minimum bias histos
		hV0Efficiency[iSp] = new TH1D(Form("hV0Efficiency_%s",SPECIES[iSp]),"; V0 pT (GeV/#it{c}); Efficiency",NPTBINS2,XBINS2);
		TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

		tV0massRCMB[iSp]->Draw(Form("lPt>>hV0Efficiency_%s",SPECIES[iSp]),"","goff");
		tV0PtMCMB[iSp]->Draw("lPt>>hDen","","goff");

		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Efficiency[iSp]->Divide(hDen);
		delete hDen;

		// rt histos
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{		
			hV0EfficiencyRt[iSp][iReg] = new TH1D(Form("hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]),"; V0 pT (GeV/#it{c}); Efficiency",NPTBINS2,XBINS2);
			TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

			tV0massRt[iSp][1][iReg]->Draw(Form("lPt>>hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]),"","goff");
			tV0PtMCRt[iSp][iReg]->Draw("lPt>>hDen","","goff");

			hV0EfficiencyRt[iSp][iReg]->GetYaxis()->SetRangeUser(0.,0.65);
			hV0EfficiencyRt[iSp][iReg]->GetXaxis()->SetRangeUser(0.,14.0);

			hV0EfficiencyRt[iSp][iReg]->Divide(hDen);
			delete hDen;
		}
	}

}

void MyAnalysisV0::DoLambdaFeeddown() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		//hV0Feeddown[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Feeddown_%s",SPECIES[iSp]));

		hV0Feeddown[iSp]->SetTitle("; V0 pT (GeV/#it{c}); Feed-down fraction");
		hV0Feeddown[iSp]->GetYaxis()->SetRangeUser(0.,0.30);
		hV0Feeddown[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Feeddown[iSp]->Divide(hV0Pt[iSp][1][0][0]);
		//hV0Efficiency[iSp]->Write();
	}

}
