#include <iostream>

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
#include "MyEvent.h"
#include "MyTrack.h"
#include "MyParticle.h"
#include "MyV0.h"
#include "MyHandler.h"

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

#include "TransverseSpherocity.h"


using namespace std;
using namespace V0consts;
using namespace RooFit;

ClassImp(MyAnalysisV0)

MyAnalysisV0::MyAnalysisV0() {

}

Int_t MyAnalysisV0::Init() {

	TString dfName(this->GetName());
	dfName = Form("%s_%i",dfName.Data(),mHandler->nAnalysis());
	mDirFile = new TDirectoryFile(dfName,dfName,"",mHandler->file());		// removed mother dir -- unnecessary?
	mDirFile->cd();

	mFlagMC = mHandler->GetFlagMC();
	mFlagHist = mHandler->GetFlagHist();

	printf("Initialising analysis %s with flag MC %i and Hist %i \n", 
		this->GetName(), mFlagMC, mFlagHist);

	TH1::SetDefaultSumw2(1);
	TH1::AddDirectory(kTRUE);
	if (mFlagHist)	BorrowHistograms();
	else 			CreateHistograms();

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
	if (!mHandler->event()) return 1;
	MyEvent event(mHandler->event());
	hEventMonitor->Fill(1);


	
	// EVENT SELECTION
	hEventType->Fill(EVENTTYPES[0],1);	//preES MB
	hEventVz->Fill(event.GetZ());
	

	switch (ClassifyEvent(event,mHandler->getNtracks())) {
		default : break;
		case (-2) : hEventType->Fill(EVENTTYPES[20],1);		// MB, rejected by Event Selection
			return 0;
			break;
		case (-1) : hEventType->Fill(EVENTTYPES[21],1);		// MB, passed ES, no vertex
			return 0;
			break;
		case (0) : hEventType->Fill(EVENTTYPES[22],1);		// MB, passed ES, bad vertex
			return 0;
			break;
		case (1) : hEventType->Fill(EVENTTYPES[23],1);		// MB, passed ES, good vertex
			break;
	}


	hEventMonitor->Fill(2);


	// BUG HOTFIX FOR AURORATREES
	MyV0 bugfix;
	if (mHandler->v0(0)) { bugfix = MyV0(mHandler->v0(0));
	bugfix.SetHandler(mHandler);
	if (TMath::Abs(bugfix.GetRadius()-bugR) < 0.0001
		&& TMath::Abs(bugfix.GetPt()-bugPt) < 0.0001) return 0;
	bugR = bugfix.GetRadius(); bugPt = bugfix.GetPt(); }
	hEventMonitor->Fill(3);

	// EVENT CLASSIFICATION
	enum { multMB, V0M, NCharged, V0M01, NCharged01, RT };
	enum { sphMB, Jetty, Iso };
	enum { D, RC, MC };
	Int_t isEventCentral = 0;
	hEventV0MCentrality->Fill(event.GetV0MCentrality());
	hEventRefMult->Fill(event.GetRefMult());
	hEventV0MCentvRefMult->Fill(event.GetRefMult(),event.GetV0MCentrality());
	Bool_t isEventFHM 	= IsCentral(event,V0M);				// forward-rapidity high-multiplicity 0-10
	Bool_t isEventMHM 	= IsCentral(event,NCharged);		// mid-rapidity high-multiplicity 0-10
	Bool_t isEventFHM01 = IsCentral(event,V0M01);			// forward-rapidity high-multiplicity 0-1
	Bool_t isEventMHM01 = IsCentral(event,NCharged01);		// mid-rapidity high-multiplicity 0-1
	if (isEventFHM)		hEventType->Fill(EVENTTYPES[2],1);
	if (isEventMHM)		hEventType->Fill(EVENTTYPES[3],1);
	if (isEventFHM01)	hEventType->Fill(EVENTTYPES[24],1);
	if (isEventMHM01)	hEventType->Fill(EVENTTYPES[25],1);
	
	isEventCentral = (isEventFHM || isEventMHM);

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
	Int_t nTracks = mHandler->getNtracks();
	Int_t nParticles = (mFlagMC) ? mHandler->getNparticles() : 0;

	ptLead = -99.; 
	phiLead = -99.;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr)); t.SetHandler(mHandler);

		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;
		if (SelectTrack(t)) {				// strict cuts, has primary dca cut
			if (t.GetPt() > ptLead) {		// search for leading track among primaries
				ptLead = t.GetPt();
				phiLead = t.GetPhi();	}
		}

		if (!t.IskITSrefit()) 		continue;	// cuts for tracks entering spherocity
		if (!t.IsTPCOnlyRefit())	continue;
		if (t.GetPt() < 0.15)		continue;

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
		if (isEventFHM01) {
			mTS[0][V0M01]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[0][V0M01]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
		if (isEventMHM01) {
			mTS[0][NCharged01]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[0][NCharged01]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
	}
	hLeadPhivPt->Fill(ptLead,phiLead);

	// SPHERO CALCULATION ON PARTICLE LEVEL IN MC
	if (mFlagMC && isEventCentral) {
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);
			
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.GetPt() > 0.15 && p.IsPrimary())	{	// cuts for particles entering spherocity

				mTS[2][0]->AddTrack(p.GetPx(), p.GetPy());
				mTSNorm[2][0]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventFHM) mTS[2][V0M]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventFHM) mTSNorm[2][V0M]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventMHM) mTS[2][NCharged]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventMHM) mTSNorm[2][NCharged]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventFHM01) mTS[2][V0M01]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventFHM01) mTSNorm[2][V0M01]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventMHM01) mTS[2][NCharged01]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventMHM01) mTSNorm[2][NCharged01]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

			}	
		}
	}

	// EVENT SPHEROCITY CLASSIFICATION
	Bool_t isEventIso[NMULTI-1] 	= { 0, 0, 0, 0, 0 };
	Bool_t isEventJetty[NMULTI-1] 	= { 0, 0, 0, 0, 0 };
	Bool_t isEventIsoMC[NMULTI-1] 	= { 0, 0, 0, 0, 0 };
	Bool_t isEventJettyMC[NMULTI-1] = { 0, 0, 0, 0, 0 };

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
		eventTS = eventTS2[tsDNorm];		// select a default spherocity

		isEventIso[V0M]		= (eventTS > cuts::EV_SPH_ISO[V0M-1] && eventTS < 1.) ;
		isEventJetty[V0M]	= (eventTS < cuts::EV_SPH_JETTY[V0M-1] && eventTS > 0.);
		isEventIsoMC[V0M]	= mFlagMC && (eventTS2[tsMCNorm] > cuts::EV_SPH_ISO[V0M-1] && eventTS2[tsMCNorm] < 1.) ;		// in mc study also particle level 
		isEventJettyMC[V0M]	= mFlagMC && (eventTS2[tsMCNorm] < cuts::EV_SPH_JETTY[V0M-1] && eventTS2[tsMCNorm] > 0.) ;

		isEventIso[NCharged]		= (eventTS > cuts::EV_SPH_ISO[NCharged-1] && eventTS < 1.) ;
		isEventJetty[NCharged]		= (eventTS < cuts::EV_SPH_JETTY[NCharged-1] && eventTS > 0.);
		isEventIsoMC[NCharged]		= mFlagMC && (eventTS2[tsMCNorm] > cuts::EV_SPH_ISO[NCharged-1] && eventTS2[tsMCNorm] < 1.) ;		// in mc study also particle level 
		isEventJettyMC[NCharged]	= mFlagMC && (eventTS2[tsMCNorm] < cuts::EV_SPH_JETTY[NCharged-1] && eventTS2[tsMCNorm] > 0.) ;

		isEventIso[V0M01]		= (eventTS > cuts::EV_SPH_ISO[V0M01-1] && eventTS < 1.) ;
		isEventJetty[V0M01]		= (eventTS < cuts::EV_SPH_JETTY[V0M01-1] && eventTS > 0.);
		isEventIsoMC[V0M01]		= mFlagMC && (eventTS2[tsMCNorm] > cuts::EV_SPH_ISO[V0M01-1] && eventTS2[tsMCNorm] < 1.) ;		// in mc study also particle level 
		isEventJettyMC[V0M01]	= mFlagMC && (eventTS2[tsMCNorm] < cuts::EV_SPH_JETTY[V0M01-1] && eventTS2[tsMCNorm] > 0.) ;

		isEventIso[NCharged01]		= (eventTS > cuts::EV_SPH_ISO[NCharged01-1] && eventTS < 1.) ;
		isEventJetty[NCharged01]	= (eventTS < cuts::EV_SPH_JETTY[NCharged01-1] && eventTS > 0.);
		isEventIsoMC[NCharged01]	= mFlagMC && (eventTS2[tsMCNorm] > cuts::EV_SPH_ISO[NCharged01-1] && eventTS2[tsMCNorm] < 1.) ;		// in mc study also particle level 
		isEventJettyMC[NCharged01]	= mFlagMC && (eventTS2[tsMCNorm] < cuts::EV_SPH_JETTY[NCharged01-1] && eventTS2[tsMCNorm] > 0.) ;
		
		if (isEventFHM) hEventSpherocityV0M->Fill(eventTS);
		if (isEventMHM) hEventSpherocityNCharged->Fill(eventTS);
		if (isEventFHM01) hEventSpherocityV0M01->Fill(eventTS);
		if (isEventMHM01) hEventSpherocityNCharged01->Fill(eventTS);

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

	if (isEventFHM01 && isEventIso[V0M01])			hEventType->Fill(EVENTTYPES[26],1);
	if (isEventFHM01 && isEventJetty[V0M01])		hEventType->Fill(EVENTTYPES[27],1);
	if (isEventMHM01 && isEventIso[NCharged01])		hEventType->Fill(EVENTTYPES[28],1);
	if (isEventMHM01 && isEventJetty[NCharged01])	hEventType->Fill(EVENTTYPES[29],1);
	if (isEventFHM01 && isEventIsoMC[V0M01])		hEventType->Fill(EVENTTYPES[30],1);
	if (isEventFHM01 && isEventJettyMC[V0M01])		hEventType->Fill(EVENTTYPES[31],1);
	if (isEventMHM01 && isEventIsoMC[NCharged01])	hEventType->Fill(EVENTTYPES[32],1);
	if (isEventMHM01 && isEventJettyMC[NCharged01])	hEventType->Fill(EVENTTYPES[33],1);

	// TRACK DEPENDENCE ON SPHEROCITY STUDY LOOP
	nChTrans = 0;
	Int_t nChTrans0 = 0;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));	t.SetHandler(mHandler);

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

		if (isEventFHM01) {
			ProcessTrack(t,D,V0M01,sphMB);
			if (mFlagMC)				ProcessTrack(t,RC,V0M01,sphMB);
			if (isEventJetty[V0M01])	ProcessTrack(t,D,V0M01,Jetty);
			if (isEventIso[V0M01])		ProcessTrack(t,D,V0M01,Iso);			//
			if (isEventJettyMC[V0M01])	ProcessTrack(t,RC,V0M01,Jetty);		// study tracks also for true sphero, stored under RC
			if (isEventIsoMC[V0M01])	ProcessTrack(t,RC,V0M01,Iso);		} 

		if (isEventMHM01) {
			ProcessTrack(t,D,NCharged01,sphMB);
			ProcessTrack(t,RC,NCharged01,sphMB);
			if (isEventJetty[NCharged01])		ProcessTrack(t,D,NCharged01,Jetty);
			if (isEventIso[NCharged01])			ProcessTrack(t,D,NCharged01,Iso);
			if (isEventJettyMC[NCharged01])		ProcessTrack(t,RC,NCharged01,Jetty);
			if (isEventIsoMC[NCharged01])		ProcessTrack(t,RC,NCharged01,Iso);		}

	}

	// RT NCH CALCULATION
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr)); t.SetHandler(mHandler);
		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;
		
		// RT DETERMINATION
		if (!IsTrans(t.GetPhi(),phiLead))	continue;
		if (!t.IskITSrefit()) continue;		// cuts for tracks entering Nch^trans calculations
		if (!t.IsTPCOnlyRefit()) continue;
		if (t.GetPt()<0.15) continue;
		
		hNchvLeadPt->Fill(ptLead);
		nChTrans0++;
		if (ptLead < 5.) continue;
		if (ptLead > 40.) continue;
		nChTrans++;		// increment only if a leading particle was found
	}

	hNchvLeadPt2->Fill(ptLead,nChTrans0);
	eventRt = 0;
	if (ptLead>5. && ptLead < 40.) {
		hNchTrans->Fill(nChTrans);
		eventRt = (double)nChTrans/RT_DEN;
		hRt->Fill(eventRt);		}

	// CLASSIFYING EVENT IN RT BINS DURING RUNTIME (optional)
	Bool_t isEventRT = false;		
	enum { rt, rtsizeof }; //rt01, rt12, rt23, rt34, rt45, rtsizeof };
	Bool_t isRT[rtsizeof] = {false}; //, false, false, false, false};
	if (ptLead>5. && ptLead<40.) {
		isEventRT = true;
		//for (int iRt = 0; iRt < rtsizeof; ++iRt)	{
		//	isRT[iRt] = (eventRt > (Double_t)iRt 
		//		&& eventRt < (Double_t)(iRt+1));	}
		isRT[0] = isEventRT;
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
		t.SetHandler(mHandler);
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
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.IsPrimary())	{

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
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.GetPt() > 0.15 && p.IsPrimary())	{

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
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetPt() > 0.15 && p.IsPrimary())	{

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

		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);


			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.IsPrimary())	{

				hTrackPt[MC][multMB][sphMB]->Fill(p.GetPt());

				if (isEventFHM) {
					hTrackPt[MC][V0M][sphMB]->Fill(p.GetPt());
					if (isEventJettyMC[V0M])	hTrackPt[MC][V0M][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC[V0M])		hTrackPt[MC][V0M][Iso]->Fill(p.GetPt()); }

				if (isEventMHM) {
					hTrackPt[MC][NCharged][sphMB]->Fill(p.GetPt());
					if (isEventJettyMC[NCharged])	hTrackPt[MC][NCharged][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC[NCharged])		hTrackPt[MC][NCharged][Iso]->Fill(p.GetPt()); }

				if (isEventFHM01) {
					hTrackPt[MC][V0M01][sphMB]->Fill(p.GetPt());
					if (isEventJettyMC[V0M01])		hTrackPt[MC][V0M01][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC[V0M01])		hTrackPt[MC][V0M01][Iso]->Fill(p.GetPt()); }

				if (isEventMHM01) {
					hTrackPt[MC][NCharged01][sphMB]->Fill(p.GetPt());
					if (isEventJettyMC[NCharged01])		hTrackPt[MC][NCharged01][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC[NCharged01])		hTrackPt[MC][NCharged01][Iso]->Fill(p.GetPt()); }


				if (isEventRT && IsTrans(p.GetPhi(),phiLead)) {
					hTrackPt[MC][RT][sphMB]->Fill(p.GetPt());
					for (int iRt = 0; iRt < rtsizeof; ++iRt) {
						if (isRT[iRt])	hTrackPt[MC][RT][3+iRt]->Fill(p.GetPt());	}
				}

			}

			//if (p.GetPdgCode() == 310) cout << "eta " << p.GetEta() << " y " << p.GetY() << endl;
				
			if (!SelectParticle(p)) continue;		// also contains Xi's !
			
			PartLabels.push_back(p.GetLabel());		// filling particles in a vector to avoid unnecessary nesting loops
			PartIds.push_back(iP);					// only particles in good |eta| with K0s,L,Lbar,Xi pdg's

			if (!p.IsPrimary()) continue;

			if (p.GetPdgCode() == 3312) hV0FeeddownMotherPt[2]->Fill(p.GetPt());		// filling primary Xi
			if (p.GetPdgCode() == -3312) hV0FeeddownMotherPt[3]->Fill(p.GetPt());

			for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
				if (p.GetPdgCode() == PDG_IDS[iSp])		{

					// FILLING MC INFORMATION FOR V0s particles
					hV0Pt[iSp][MC][multMB][sphMB]->Fill(p.GetPt());
					hV0Eta[iSp][MC][multMB][sphMB]->Fill(p.GetEta());
					hV0Y[iSp][MC][multMB][sphMB]->Fill(p.GetY());

					tV0PtMCMB[iSp]->Fill(p.GetPt());

					if (isEventFHM) {
						hV0Pt[iSp][MC][V0M][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC[V0M])	hV0Pt[iSp][MC][V0M][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC[V0M])		hV0Pt[iSp][MC][V0M][Iso]->Fill(p.GetPt());		}

					if (isEventMHM) {
						hV0Pt[iSp][MC][NCharged][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC[NCharged])	hV0Pt[iSp][MC][NCharged][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC[NCharged])		hV0Pt[iSp][MC][NCharged][Iso]->Fill(p.GetPt());		}

					if (isEventFHM01) {
						hV0Pt[iSp][MC][V0M01][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC[V0M01])		hV0Pt[iSp][MC][V0M01][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC[V0M01])		hV0Pt[iSp][MC][V0M01][Iso]->Fill(p.GetPt());		}

					if (isEventMHM01) {
						hV0Pt[iSp][MC][NCharged01][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC[NCharged01])		hV0Pt[iSp][MC][NCharged01][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC[NCharged01])		hV0Pt[iSp][MC][NCharged01][Iso]->Fill(p.GetPt());		}


					if (isEventRT)	{
							Int_t region = WhatRegion(p.GetPhi(),phiLead);

							if (iSp>0) tV0PtMCRt[iSp][region]->Fill(p.GetPt());
						}

					if (isEventRT && IsTrans(p.GetPhi(),phiLead)) {
						hV0Pt[iSp][MC][RT][sphMB]->Fill(p.GetPt());
						for (int iRt = 0; iRt < rtsizeof; ++iRt) {
							if (isRT[iRt])	hV0Pt[iSp][MC][RT][3+iRt]->Fill(p.GetPt());	}
					}

				}
			}

		}

	}

	// V0 DATA ANALYSIS: V0 CANDIDATES LOOP
	hEventMonitor->Fill(6);
	Int_t nV0s = mHandler->getNv0s();
	for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
		
		hV0Monitor->Fill(0);
		if (!mHandler->v0(iV0)) continue;
		MyV0 v0(mHandler->v0(iV0));
		v0.SetHandler(mHandler);


		if (mFlagMC) {
			MyParticle v0mc; v0mc.SetHandler(mHandler);
			Bool_t MCfound = false;

			// LOOKING IF V0 HAS A MC PARTNER
			for (unsigned int iP = 0; iP < PartLabels.size(); ++iP)	{
				if (PartLabels[iP] == v0.GetMCLabel()) {		// perhaps ask for pdgid as well?
					
					v0mc = MyParticle(mHandler->particle(PartIds[iP])); v0mc.SetHandler(mHandler); v0mc.SetLabel(PartIds[iP]);

					//if (v0mc.GetPdgCode() != v0.GetMCPdgCode()) {
						//printf("Codes are %i and %i and pt %f \n",v0mc.GetPdgCode(),v0.GetMCPdgCode(), v0.GetPt());
						//cout << "WRONG PARTICLE-V0 ASSOCIATION" << endl;
						//return 1; }
					//}

					if (MCfound) {
						cout << "PARTICLE FOUND TWICE" << endl;
						return 1; }

					if (v0mc.GetPdgCode() == v0.GetMCPdgCode()) MCfound = true;
					// found a MC particle with the same label and PDGCode
					// MyV0::GetMCPdgCode() also requires both daughters to point to the same mother
				}
			}

			for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
				//cout << "Testing 1 iSp " << iSp << " with pdg " << v0.GetMCPdgCode() << " label " << v0.GetMCLabel() << endl;
					
				if (IsV0(v0,iSp,RC)) {


					//if (iSp==1) cout << "detected k0s with label " << v0.GetMCLabel() << endl;
					//else hV0FeeddownPDG[iSp]->Fill(v0mc.GetMotherPdgCode()); //printf("code is %i and %i \n", v0mc.GetPdgCode(), v0mc.GetMotherPdgCode());
						
						

					// ask for primary
					// get the MC label of mother
					// get the Xi
					// fill histo  l.pt and xi.pt
					//printf("primary %i and grandma pdg code %i \n", v0.IsMCPrimary(), v0.GetMCPrimaryPdgCode());
					/*if (!v0.IsMCPrimary()) {			// ask for secondary
						hV0FeeddownPDG[iSp]->Fill(v0mc.GetMotherPdgCode());
						
						if (iSp>1 && TMath::Abs(v0mc.GetMotherPdgCode()) == 3312) {	// consider Xi->L
							MyParticle priMC;
							for (unsigned int iP = 0; iP < PartLabels.size(); ++iP)	{	// find particle Xi
								if (PartLxiabels[iP] == v0.GetMCPrimaryLabel()) {
									priMC = MyParticle(mHandler->particle(PartIds[iP]));
								}
							}

							if (TMath::Abs(priMC.GetPdgCode()) != 3312) {
								cout << "GRANDMOTHER NOT XI" << endl;
								return 1;
							}

							hV0FeeddownMatrix[iSp]->Fill(priMC.GetPt(),v0.GetPt());
						}
						continue;
					}*/

					// also require daughter pdg id for association
					Bool_t properDaughters = ((v0.GetPosTrackPdg() == PDG_IDS_DPOS[iSp]) && (v0.GetNegTrackPdg() == PDG_IDS_DNEG[iSp]));

					if (!v0.IsMCPrimary() && MCfound) {
							// using secondary particles to build the feeddown matrix


							//cout << "v0 not primary " << iSp << endl;
							//if (iSp>1 && TMath::Abs(v0mc.GetMotherPdgCode()) == 3312) {
						if (iSp>1) {
							MyParticle xiMC; xiMC.SetHandler(mHandler);
						
							for (unsigned int iP = 0; iP < PartLabels.size(); ++iP)	{
								xiMC = MyParticle(mHandler->particle(PartIds[iP])); xiMC.SetHandler(mHandler); xiMC.SetLabel(PartIds[iP]);
								//cout << "found particle w id " << xiMC.GetPdgCode() << endl;
								if (xiMC.GetPdgCode() == 3312 && iSp==2 && xiMC.IsPrimary() && v0mc.GetMotherPdgCode() == 3312)	{
									hV0FeeddownMatrix[iSp]->Fill(xiMC.GetPt(),v0.GetPt());
								}
								if (xiMC.GetPdgCode() == -3312 && iSp==3 && xiMC.IsPrimary() && v0mc.GetMotherPdgCode() == -3312)	{
									hV0FeeddownMatrix[iSp]->Fill(xiMC.GetPt(),v0.GetPt());
								}
							}
						}
					}

					

					if (!v0.IsMCPrimary()) continue;	// considering only primaries for efficiency
					if (!MCfound) continue;
					if (!properDaughters) continue;

					// MC v RC histograms
					hV0PtRCvMC[iSp]->Fill(v0mc.GetPt(), v0.GetPt());
					hV0EtaRCvMC[iSp]->Fill(v0mc.GetEta(), v0.GetEta());
					hV0PhiRCvMC[iSp]->Fill(v0mc.GetPhi(), v0.GetPhi());
						

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

					if (isEventFHM01) {
						ProcessV0toHist(v0,iSp,RC,V0M01,sphMB);
						if (isEventJettyMC[V0M01])		ProcessV0toHist(v0,iSp,RC,V0M01,Jetty);			// study RC spectra for true spher. ev classification
						if (isEventIsoMC[V0M01])		ProcessV0toHist(v0,iSp,RC,V0M01,Iso);
					}
						
					if (isEventMHM01) {
						ProcessV0toHist(v0,iSp,RC,NCharged01,sphMB);
						if (isEventJettyMC[NCharged01])		ProcessV0toHist(v0,iSp,RC,NCharged01,Jetty);
						if (isEventIsoMC[NCharged01])		ProcessV0toHist(v0,iSp,RC,NCharged01,Iso);
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
		
		hV0Radius->Fill(v0.GetRadius());
		hV0ProperT->Fill(v0.GetIML()*v0.GetRadius()/v0.GetPt());
		
		// STUDY OF V0s IN DATA
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


				if (isEventFHM01) {
					ProcessV0toHist(v0,iSp,D,V0M01,sphMB);
					if (isEventJetty[V0M01])	ProcessV0toHist(v0,iSp,D,V0M01,Jetty);
					if (isEventIso[V0M01])		ProcessV0toHist(v0,iSp,D,V0M01,Iso); }

				if (isEventMHM01) {
					ProcessV0toHist(v0,iSp,D,NCharged01,sphMB);
					if (isEventJetty[NCharged01])	ProcessV0toHist(v0,iSp,D,NCharged01,Jetty);
					if (isEventIso[NCharged01])		ProcessV0toHist(v0,iSp,D,NCharged01,Iso); }
				

				if (isEventRT) {
					
					Int_t region = WhatRegion(v0.GetPhi(),phiLead);
					if (iSp>0) ProcessV0toTree(v0,iSp,D,RT+region);

					ProcessV0toHist(v0,iSp,D,RT+region,sphMB);
					for (int iRt = 0; iRt < rtsizeof; ++iRt) {
						if (isRT[iRt])	ProcessV0toHist(v0,iSp,D,RT+region,3+iRt);	}
				}

			}
		}		

	}

	hEventMonitor->Fill(7);
	return 0;	
}

Bool_t MyAnalysisV0::SelectEvent(MyEvent &ev, Int_t flag) {

	// NOT USED (CLASSIFYEVENT USED INSTEAD)

	ev.SetCheckFlag(flag);
	if (!ev.IsGoodAliEvent()) return false;

	//if (!ev.IsGoodAliEvent())			return false;
	//if(bV0s->GetEntriesFast() < 1)		return false;
	//if(bTracks->GetEntriesFast() < 1)		return false;
	//hEventCuts->Fill(0);
	//if (!ev.AcceptVertex())				return false;
	//hEventCuts->Fill(1);
	//if (ev.IsPileupFromSPD())			return false;
	//hEventCuts->Fill(2);
	//if (!ev.CheckFlag())				return false;
	//if (!ev.IsCollisionCandidate())		return false;
	//if (ev.GetCentralityQuality() != 0)	return false;

	return true;
}

Int_t MyAnalysisV0::ClassifyEvent(MyEvent &event, Int_t ntracks)
{

	// usage of flags also emulated for ESDs to ensure we're doing the same thing
	/*
	enum EventFlags_t {
	    kNotPileupInSPD = 1,
	    kNotPileupInMV = 2,
	    kNotPileupInMB = 4,
	    kINELgtZERO = 8,
	    kNoInconsistentVtx = 16,
	    kNoV0Asym = 32,
	    kVertexSelected2015pp=64,
	    kSPDandTrkVtxExists=128,
	    kPassProximityCut=256,
	    kAll = 511
  	};
  	*/

	event.SetCheckFlag(1);
	if(!event.IsGoodAliEvent()) {
		//    cout << "Pileup or similar" << endl;
		return -2;
	}

	if(!event.HasVertex()) {
		return -1;
	}

	event.SetCheckFlag(448);		// = 256 + 128 + 64
	if(!event.IsGoodAliEvent()) {
		return -1;
	}
  
	if(TMath::Abs(event.GetZ())>10) {
		return 0;
	}

	if(ntracks<1) { 
		return 1;
	}

  	// Update TOF response
  	//AliAnalysisPIDCascadeTrack* track = (AliAnalysisPIDCascadeTrack*)trackArray->At(0);
  	//track->UpdateTOFResponse(event);
  	// ^^ IS THIS NECESSARY?

  return 1;
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

		case 3: 
			if (ev.GetV0MCentrality() < 1. 
				&& ev.GetRefMult() > 9.9)		return true;
			break;

		case 4: 
			if (ev.GetRefMult() < cuts::EV_MHM01[1] 
				&& ev.GetRefMult() > cuts::EV_MHM01[0])		return true;
			break;
	}

	return false;
}

Bool_t MyAnalysisV0::ProcessV0toHist(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph) {
	
	//printf("V0 p_{T} is %f, args are %i , %i , %i \n", v0.GetPt(), Sp, Mu, Sph);

	hV0Pt[Sp][Type][Mu][Sph]->Fill(v0.GetPt());
	hV0Eta[Sp][Type][Mu][Sph]->Fill(v0.GetEta());
	hV0Y[Sp][Type][Mu][Sph]->Fill(v0.CalculateY(Sp));
	
	hV0EtavY[Sp][Type][Mu][Sph]->Fill(v0.CalculateY(Sp),v0.GetEta());
	Double_t v0mass[] 	= {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	Double_t v0massKF[] = {0., v0.GetKFIMK0s(), v0.GetKFIML(), v0.GetKFIMLbar()};
	hV0IMvPt[Sp][Type][Mu][Sph]->Fill(v0.GetPt(),v0mass[Sp]);

	if (Mu>4 && Type==0 && Sph == 0 && Sp==1) {
		hV0DPhivNchTrans->Fill(nChTrans,mHandler->DeltaPhi(phiLead,v0.GetPhi()));
	}

	MyTrack trP(v0.GetPosTrack());	trP.SetHandler(mHandler);
	MyTrack trN(v0.GetNegTrack());	trN.SetHandler(mHandler); 
	if (Sph==0 && Mu == 0) {
		if (Sp!=2) hV0DpiNsigTPCvpt[Sp][Type]->Fill(trP.GetPt(),trP.GetNSigmaPionTPC());
		if (Sp!=3) hV0DpiNsigTPCvpt[Sp][Type]->Fill(trN.GetPt(),trN.GetNSigmaPionTPC());
		if (Sp==2) hV0DprNsigTPCvpt[Sp][Type]->Fill(trP.GetPt(),trP.GetNSigmaProtonTPC());
		if (Sp==3) hV0DprNsigTPCvpt[Sp][Type]->Fill(trN.GetPt(),trN.GetNSigmaProtonTPC());

		hV0KFIMvIM[Sp][Type]->Fill(v0mass[Sp],v0massKF[Sp]);
		hV0DeltaIMvPt[Sp][Type]->Fill(v0.GetPt(),v0mass[Sp]-v0massKF[Sp]);

		// CUT STUDIES
		hV0BaselineIMvPt[Sp][Type]->Fill(v0.GetPt(),v0mass[Sp]);
		if (v0.HasFastSignal()) hV0FastSignalIMvPt[Sp][Type]->Fill(v0.GetPt(),v0mass[Sp]);
	}

	return true;	
}

Bool_t MyAnalysisV0::ProcessV0toTree(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu) {
	
	Double_t v0mass[] = {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	Double_t v0massKF[] = {0., v0.GetKFIMK0s(), v0.GetKFIML(), v0.GetKFIMLbar()};
	
	if (Mu > 4) {
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

	//Bool_t isPromising = (TMath::Abs(v0.GetIMK0s())<0.002) && (1) && (Type==mFlagMC) && (Sp == 1);
	Double_t v0mass[] = {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	Double_t v0massKF[] = {0., v0.GetKFIMK0s(), v0.GetKFIML(), v0.GetKFIMLbar()};
	Int_t cutN = 1;

	if (Type>0) {
		//cout << Sp << ": checking " << v0.GetMCPdgCode() << " vs " << PDG_IDS[Sp] << " label " << v0.GetMCLabel() << endl; 
		if (v0.GetMCPdgCode() != PDG_IDS[Sp])	return false; 
		//cout << "passed \n"; 
	}
		
		//if (!v0.IsMCPrimary())					return false; } // always 1 ?

	//if (Type==0 && Sp==2) if (TMath::Abs(v0.GetIML()) > 0.01) return false;

	if (!v0.IsOffline())	return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //1
	if (v0.GetEta() < cuts::V0_ETA[0]) 	return false;
	//if (v0.CalculateY(Sp) < cuts::V0_Y[0]) 	return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //2
	if (v0.GetEta() > cuts::V0_ETA[1]) 	return false;
	//if (v0.CalculateY(Sp) > cuts::V0_Y[1]) 	return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //3
	if (v0.GetPt() < cuts::V0_PT[0]) 		return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //4
	if (v0.GetPt() > cuts::V0_PT[1]) 		return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //5
	if (v0.GetDCAdd() > cuts::V0_DCADD) 	return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //6
	if (v0.GetCPA() < cuts::V0_CPA) 		return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //7
	if (v0.GetRadius() < cuts::V0_R[0]) 	return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //8
	if (v0.GetRadius() > cuts::V0_R[1]) 	return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //9
	if (!v0.HasFastSignal() && cuts::V0_FASTSIGNAL)		return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //10
	
	
	//if (!Sp) return true;
	MyTrack trP(v0.GetPosTrack());
	trP.SetHandler(mHandler); 
	MyTrack trN(v0.GetNegTrack());
	trN.SetHandler(mHandler);
	
	if (!SelectV0Daughter(trP)) return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //11
	if (!SelectV0Daughter(trN)) return false;
	hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //12

	
	
	switch (Sp) {
		default :
			
			break;
		case 1 	: 
			
			// K0s
			//if (*(v0.CalculateAP()+1) < cuts::K0S_AP*TMath::Abs(*(v0.CalculateAP()+0))) return false;
			//if (v0mass[Sp]*v0.GetRadius()/v0.GetPt() > cuts::K0S_TAU)	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //13
			if (TMath::Abs(v0mass[2]) < cuts::K0S_COMP_M) 			return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //14
			if (TMath::Abs(v0mass[3]) < cuts::K0S_COMP_M) 			return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //15
			//cout << "tpc " << trP.GetNSigmaPionTPC() << " " << trN.GetNSigmaPionTPC() << endl;
			if (!mFlagMC && trP.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0]) 	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //16
			if (!mFlagMC && trP.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1]) 	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //17
			if (!mFlagMC && trN.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0]) 	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //18
			if (!mFlagMC && trN.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1]) 	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //19
			
			break;
		case 2 	: // L
		
			if (v0.GetPt() < cuts::L_PT[0]) 	return false;
			if (v0.GetPt() > cuts::L_PT[1]) 	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //13
			if (v0mass[Sp]*v0.GetRadius()/v0.GetPt() > cuts::L_TAU)	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //14
			if (TMath::Abs(v0mass[1]) < cuts::L_COMP_M) 			return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //15
			if (v0.GetCPA() < cuts::L_CPA)					 		return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //16
			if (!mFlagMC && trP.GetNSigmaProtonTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //17
			if (!mFlagMC && trP.GetNSigmaProtonTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //18
			if (!mFlagMC && trN.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //19
			if (!mFlagMC && trN.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //20
			
			break;
		case 3 	: // Lbar
		
			if (v0.GetPt() < cuts::L_PT[0]) 	return false;
			if (v0.GetPt() > cuts::L_PT[1]) 	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //13
			if (v0mass[Sp]*v0.GetRadius()/v0.GetPt() > cuts::L_TAU)	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //14
			if (TMath::Abs(v0mass[1]) < cuts::L_COMP_M) 			return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //15
			//if (TMath::Abs(v0mass[2]) < cuts::K0S_COMP_M) 			return false;
			if (v0.GetCPA() < cuts::L_CPA)					 		return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //16
			if (!mFlagMC && trP.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //17
			if (!mFlagMC && trP.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //18
			if (!mFlagMC && trN.GetNSigmaProtonTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //19
			if (!mFlagMC && trN.GetNSigmaProtonTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //20
			
			break;	}

			
		hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //21

	
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
	//if (!tr.IsITSTPC2011())	return false;

	if ( cuts::V0_D_GOODTRACK && !tr.IsGoodV0daughter())	return false;

	//printf("values are %f and %f \n", tr.GetTPCnc(),tr.GetTPCnc()/tr.GetTPCNclsF());

	if (tr.GetTPCnc() < cuts::V0_D_NCR)							return false;
	if (tr.GetTPCNclsF() > 0 && (tr.GetTPCnc()/tr.GetTPCNclsF()) < cuts::V0_D_NRATIO)		return false;

	return true;
}

Bool_t MyAnalysisV0::SelectParticle(MyParticle &p) {

	if (p.GetEta() < cuts::V0_ETA[0]) 		return false;
	if (p.GetEta() > cuts::V0_ETA[1]) 		return false;
	//if (p.GetY() < cuts::V0_Y[0]) 		return false;
	//if (p.GetY() > cuts::V0_Y[1]) 		return false;
	if (p.GetPdgCode() != PDG_IDS[1]
		&& p.GetPdgCode() != PDG_IDS[2]
		&& p.GetPdgCode() != PDG_IDS[3]
		&& p.GetPdgCode() != 3312			// also Xi pm for feed-down study
		&& p.GetPdgCode() != -3312 )		return false;
	
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
	hEventVz				= new TH1D("hEventVz","; vz; Entries",400,-50,50);
	hEventV0MCentrality		= new TH1D("hEventV0MCentrality","; V0M Centrality; Entries",300,0,150);
	hEventRefMult			= new TH1D("hEventRefMult","; Reference multiplicity; Entries",150,0,150);
	hEventV0MCentvRefMult	= new TH2D("hEventV0MCentvRefMult","; Reference multiplicity; V0M Centrality; Entries"
		,150,0,150,300,0,150);

	hEventSpherocityV0M			= new TH1D("hEventSpherocityV0M","; S_{O}; Entries",400,-0.1,1.1);
	hEventSpherocityNCharged	= new TH1D("hEventSpherocityNCharged","; S_{O}; Entries",400,-0.1,1.1);
	hEventSpherocityV0M01		= new TH1D("hEventSpherocityV0M01","; S_{O}; Entries",400,-0.1,1.1);
	hEventSpherocityNCharged01	= new TH1D("hEventSpherocityNCharged01","; S_{O}; Entries",400,-0.1,1.1);
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
				
		if (iMu > 4 && (iSph < 3 && iSph)) continue;
		if (iMu < 5 && iSph > 2) continue; 
		hTrackPt[iType][iMu][iSph]			= new TH1D(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";track p_{T} (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hTrackEtavPhi[iType][iMu][iSph]		= new TH2D(Form("hTrackEtavPhi_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";track #phi; track #eta", 400, -0.2, 6.4, 400, -1., 1.);		

	} } } 

	// MC PARTICLE HISTOGRAMS
	hParticlePrimaryvPDG				= new TH2D("hParticlePrimaryvPDG", "; PDG id; Primary", 10000,-5000,5000,2,-0.5,1.5);
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
		hProtonNchTransvPt[iReg]		= new TH2D(Form("hProtonNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
		hPionNchTransvPt[iReg]			= new TH2D(Form("hPionNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
		hLambdaNchTransvPt[iReg]		= new TH2D(Form("hLambdaNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
		hK0sNchTransvPt[iReg]			= new TH2D(Form("hK0sNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS2, XBINS2, 50, -0.5, 49.5);
	}

	// V0 HISTOGRAMS
	hV0Radius		= new TH1D("hV0Radius",";V0 radius; Entries",400,0.,150.);
	hV0ProperT		= new TH1D("hV0ProperT",";V0 proper #tau; Entries",400,0.,150.);
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu > 4 && (iSph < 3 && iSph)) continue;
		if (iMu < 5 && iSph > 2) continue; 
		hV0Pt[iSp][iType][iMu][iSph]			= new TH1D(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 p_{T} (GeV/#it{c}); Entries",								NPTBINS,XBINS);
			//";V0 p_{T} (GeV/#it{c}); Entries",								400, 0, 20);
		hV0Eta[iSp][iType][iMu][iSph]			= new TH1D(Form("hV0Eta_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 #eta; Entries", 										200, -1., 1.);
		hV0Y[iSp][iType][iMu][iSph]				= new TH1D(Form("hV0Y_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 y; Entries", 										200, -1., 1.);

		hV0EtavY[iSp][iType][iMu][iSph]			= new TH2D(Form("hV0EtavY_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 #eta; V0 y; Entries", 										200, -1., 1., 200, -1., 1.);
		hV0IMvPt[iSp][iType][iMu][iSph]		= new TH2D(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.1, 0.1);


	} } } }

	Int_t nType = (mFlagMC) ? 2 : 1;
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
		hV0DpiNsigTPCvpt[iSp][iType]		= new TH2D(Form("hV0DpiNsigTPCvpt_%s_%s",SPECIES[iSp],TYPE[iType]),
			";V0 daughter p_{T} (GeV/#it{c}); V0 daughter n#sigma_{TPC}^{#pi}; Entries", 										300, 0., 15., 200, -10., 10.);
		hV0DprNsigTPCvpt[iSp][iType]		= new TH2D(Form("hV0DprNsigTPCvpt_%s_%s",SPECIES[iSp],TYPE[iType]),
			";V0 daughter p_{T} (GeV/#it{c}); V0 daughter n#sigma_{TPC}^{p}; Entries", 										300, 0., 15., 200, -10., 10.);

		hV0KFIMvIM[iSp][iType]		= new TH2D(Form("hV0KFIMvIM_%s_%s",SPECIES[iSp],TYPE[iType]),
			";V0 m (GeV/#it{c}^{2}); V0 KF m (GeV/#it{c}^{2}); Entries",		400, -0.1, 0.1, 400, -0.1, 0.1);
		hV0DeltaIMvPt[iSp][iType]	= new TH2D(Form("hV0DeltaIMvPt_%s_%s",SPECIES[iSp],TYPE[iType]),
			";V0 p_{T} (GeV/#it{c}); V0 m-m_{KF} (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 400, -0.2, 0.2);

		hV0BaselineIMvPt[iSp][iType]	= new TH2D(Form("hV0BaselineIMvPt_%s_%s",SPECIES[iSp],TYPE[iType]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Fraction removed",		NPTBINS, XBINS, 25, -0.05, 0.05);
		hV0FastSignalIMvPt[iSp][iType]	= new TH2D(Form("hV0FastSignalIMvPt_%s_%s",SPECIES[iSp],TYPE[iType]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Fraction removed",		NPTBINS, XBINS, 25, -0.05, 0.05);

	for (int iCut = 0; iCut < 25; ++iCut)	{
		
		hV0CutIMvPt[iSp][iType][iCut]			= new TH2D(Form("hV0CutIMvPt_%s_%s_%i",SPECIES[iSp],TYPE[iType],iCut),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Fraction removed",		NPTBINS, XBINS, 25, -0.05, 0.05);
	}
	
	}	}

	// V0 RC v MC HISTOGRAMS
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
		hV0PtRCvMC[iSp]			= new TH2D(Form("hV0PtRCvMC_%s",SPECIES[iSp]),
			";V0 MC p_{T} (GeV/#it{c}); V0 RC p_{T} (GeV/#it{c}); Entries",	NPTBINS,XBINS, NPTBINS,XBINS);
		hV0EtaRCvMC[iSp]			= new TH2D(Form("hV0EtaRCvMC_%s",SPECIES[iSp]),
			";V0 MC #eta; V0 RC #eta; Entries",	200, -1., 1., 200, -1., 1.);
		hV0PhiRCvMC[iSp]			= new TH2D(Form("hV0PhiRCvMC_%s",SPECIES[iSp]),
			";V0 MC #phi; V0 RC #phi; Entries",	400, -0.2, 6.4, 400, -0.2, 6.4);
	}


	// FEED-DOWN STUDY HISTOGRAMS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

		hV0Feeddown[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Feeddown_%s",SPECIES[iSp]));
		hV0FeeddownPDG[iSp] =	new TH1D(Form("hV0FeeddownPDG_%s",SPECIES[iSp]),";PDG ID;Entries",20000,-10000,10000);
		hV0FeeddownMatrix[iSp]	= new TH2D(Form("hV0FeeddownMatrix_%s",SPECIES[iSp]),";primary grandmother p_{T}; decay V0 p_{T}", NXIPTBINS, XIXBINS, NPTBINS, XBINS);
		hV0FeeddownMotherPt[iSp] 	= new TH1D(Form("hV0FeeddownMotherPt_%s",SPECIES[iSp]),";primary grandmother p_{T}; Entries", NXIPTBINS, XIXBINS);
		
	}

	// V0 CUTS STUDY HISTOGRAMS
	


	// V0 NTUPLES TO ALSO ALLOW POST-PROCESSING REBINNING
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

	if (mHandler->filehist()->Get("MyAnalysisV0_0")->ClassName() == string("TDirectoryFile")) {
		cout << "Borrowing histograms from a TDirectoryFile" << endl;
		mDirFile = (TDirectoryFile*)mHandler->filehist()->Get("MyAnalysisV0_0");}
	if (mHandler->filehist()->Get("MyAnalysisV0_0")->ClassName() == string("THashList")) {
		cout << "Borrowing histograms from a THashList" << endl;
		THashList* hashList = (THashList*)mHandler->filehist()->Get("MyAnalysisV0_0");
		while (hashList->GetEntries()) {
			mDirFile->Append(hashList->First());
			hashList->RemoveFirst();
		}
	}
	
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

	// TRACK HISTOGRAMS
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu > 4 && (iSph < 3 && iSph)) continue;
		if (iMu < 5 && iSph > 2) continue; 
		hTrackPt[iType][iMu][iSph]	= (TH1D*)mDirFile->Get(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hTrackEtavPhi[iType][iMu][iSph]	= (TH2D*)mDirFile->Get(Form("hTrackEtavPhi_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));

	} } } 

	// MC PARTICLE HISTOGRAMS

	// V0 HISTOGRAMS
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu > 4 && (iSph < 3 && iSph)) continue;
		if (iMu < 5 && iSph > 2) continue; 
		hV0Pt[iSp][iType][iMu][iSph]			= (TH1D*)mDirFile->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0Eta[iSp][iType][iMu][iSph]			= (TH1D*)mDirFile->Get(Form("hV0Eta_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0IMvPt[iSp][iType][iMu][iSph]			= (TH2D*)mDirFile->Get(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		// actually perhaps only histos which need to be processed within this analysis need to get fetched
		
	} } } }	

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

		hV0Feeddown[iSp] 	= (TH1D*)mDirFile->Get(Form("hV0Feeddown_%s",SPECIES[iSp]));
		hV0FeeddownPDG[iSp] = (TH1D*)mDirFile->Get(Form("hV0FeeddownPDG_%s",SPECIES[iSp]));
		hV0FeeddownMatrix[iSp]		= (TH2D*)mDirFile->Get(Form("hV0FeeddownMatrix_%s",SPECIES[iSp]));
		hV0FeeddownMotherPt[iSp]	= (TH1D*)mDirFile->Get(Form("hV0FeeddownMotherPt_%s",SPECIES[iSp]));

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
	mDirFile->Write();

	TH1D* hLeadPt = hLeadPhivPt->ProjectionX();
	hNchvLeadPt->Divide(hLeadPt);
	

	//MAKE RT BINNING
	hRt2 = (TH1D*)hNchTrans->Clone("hRt2");
	Double_t rt_den = hNchTrans->GetMean();
	Int_t nbins = hNchTrans->GetXaxis()->GetNbins();
	Double_t rtbins[nbins+1];

	for (int iBin = 0; iBin < nbins+1; ++iBin)	{
		//cout << "a bin " << hNchTrans->GetBinLowEdge(iBin+1) << endl;
		if (rt_den>0) rtbins[iBin] = (double)hNchTrans->GetBinLowEdge(iBin+1)/rt_den;
		//cout << "b bin " << rtbins[iBin] << endl;
	}

	hRt2->SetBins(nbins,rtbins);

	if (mFlagMC) DoEfficiency();
	//if (mFlagMC) DoEfficiencyFromTrees();
	//if (mFlagMC && !mFlagHist) DoLambdaFeeddown();
	if (mFlagMC) DoLambdaFeeddown();

	return 0;	
}

void MyAnalysisV0::DoEfficiency() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0Efficiency[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Efficiency_%s",SPECIES[iSp]));
		hV0EfficiencyEta[iSp] = (TH1D*)hV0Eta[iSp][1][0][0]->Clone(Form("hV0EfficiencyEta_%s",SPECIES[iSp]));

		hV0Efficiency[iSp]->SetTitle("; V0 p_{T} (GeV/#it{c}); Efficiency");
		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Efficiency[iSp]->Divide(hV0Pt[iSp][2][0][0]);
		hV0EfficiencyEta[iSp]->Divide(hV0Eta[iSp][2][0][0]);
	}

}

void MyAnalysisV0::DoEfficiencyFromTrees() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		
		// minimum bias histos
		hV0Efficiency[iSp] = new TH1D(Form("hV0Efficiency_%s",SPECIES[iSp]),"; V0 p_{T} (GeV/#it{c}); Efficiency",NPTBINS2,XBINS2);
		TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

		tV0massRCMB[iSp]->Draw(Form("lPt>>hV0Efficiency_%s",SPECIES[iSp]),"","goff");
		tV0PtMCMB[iSp]->Draw("lPt>>hDen","","goff");

		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Efficiency[iSp]->Divide(hDen);
		delete hDen;

		// rt histos
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{		
			hV0EfficiencyRt[iSp][iReg] = new TH1D(Form("hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]),"; V0 p_{T} (GeV/#it{c}); Efficiency",NPTBINS2,XBINS2);
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

	/*for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		//hV0Feeddown[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Feeddown_%s",SPECIES[iSp]));

		hV0Feeddown[iSp]->SetTitle("; V0 p_{T} (GeV/#it{c}); Feed-down fraction");
		hV0Feeddown[iSp]->GetYaxis()->SetRangeUser(0.,0.30);
		hV0Feeddown[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Feeddown[iSp]->Divide(hV0Pt[iSp][1][0][0]);
		//hV0Efficiency[iSp]->Write();
	}*/

	// NORMALIZING 
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{

		Int_t nCols = hV0FeeddownMatrix[iSp]->GetNbinsX();
		Int_t nRows = hV0FeeddownMatrix[iSp]->GetNbinsY();
		//cout << "1mother pt at " << hV0FeeddownMotherPt[iSp] << endl;
		//hV0FeeddownMotherPt[iSp] = hV0FeeddownMatrix[iSp]->ProjectionX(Form("hV0FeeddownMotherPt_%s",SPECIES[iSp]),0,-1);
		//cout << "2mother pt at " << hV0FeeddownMotherPt[iSp] << endl;
		
		//hV0FeeddownMotherPt[iSp]->Scale(1,"width");

		for (int iC = 1; iC < nCols+1; ++iC)	{
			Double_t integral = hV0FeeddownMotherPt[iSp]->Integral(iC,iC);//,1,nRows);
			for (int iR = 1; iR < nRows+1; ++iR) {
				
				Double_t binContent = hV0FeeddownMatrix[iSp]->GetBinContent(iC,iR);
				if (binContent>0) hV0FeeddownMatrix[iSp]->SetBinContent(iC,iR,binContent/integral);
			}
		}
	}

}
