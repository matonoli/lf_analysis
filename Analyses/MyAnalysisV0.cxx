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

	bugR = 0;
	bugPt = 0;

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

	Double_t ptLead = -99., phiLead = -99.;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));

		if (!SelectTrack(t)) continue;
		if (!t.IskITSrefit()) continue;

		if (t.GetPt() > ptLead) {
			ptLead = t.GetPt();
			phiLead = t.GetPhi();	}

		// for Data calc. s0 from tracks
		if (isEventCentral) {
			mTS[0][0]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[0][0]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
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
			}	
		}
	}

	// EVENT SPHEROCITY CLASSIFICATION
	Bool_t isEventIso = 0;
	Bool_t isEventJetty = 0;
	Bool_t isEventIsoMC = 0;
	Bool_t isEventJettyMC = 0;

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

		isEventIso		= (eventTS > cuts::EV_SPH_ISO && eventTS < 1.) ;
		isEventJetty	= (eventTS < cuts::EV_SPH_JETTY && eventTS > 0.);

		isEventIsoMC	= mFlagMC && (eventTS2[tsMCNorm] > cuts::EV_SPH_ISO && eventTS2[tsMCNorm] < 1.) ;		// in mc study also particle level 
		isEventJettyMC	= mFlagMC && (eventTS2[tsMCNorm] < cuts::EV_SPH_JETTY && eventTS2[tsMCNorm] > 0.) ;
		
		hEventSpherocity->Fill(eventTS);
		hEventTSMCvRC->Fill(eventTS2[tsMC],eventTS2[tsRC]);
		hEventTSNormMCvRC->Fill(eventTS2[tsMCNorm],eventTS2[tsRCNorm]);
		hEventTSMCvNorm->Fill(eventTS2[tsMC],eventTS2[tsMCNorm]);
		hEventTSRCvNorm->Fill(eventTS2[tsRC],eventTS2[tsRCNorm]);
	}

	if (isEventFHM && isEventIso)		hEventType->Fill(EVENTTYPES[6],1);
	if (isEventFHM && isEventJetty)		hEventType->Fill(EVENTTYPES[7],1);
	if (isEventMHM && isEventIso)		hEventType->Fill(EVENTTYPES[8],1);
	if (isEventMHM && isEventJetty)		hEventType->Fill(EVENTTYPES[9],1);
	if (isEventFHM && isEventIsoMC)		hEventType->Fill(EVENTTYPES[16],1);
	if (isEventFHM && isEventJettyMC)	hEventType->Fill(EVENTTYPES[17],1);
	if (isEventMHM && isEventIsoMC)		hEventType->Fill(EVENTTYPES[18],1);
	if (isEventMHM && isEventJettyMC)	hEventType->Fill(EVENTTYPES[19],1);

	// TRACK SPHEROCITY STUDY LOOP
	Int_t nChTrans = 0;
	Int_t nChTrans0 = 0;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));

		if (!SelectTrack(t)) continue;

		ProcessTrack(t,D,multMB,sphMB);
		ProcessTrack(t,RC,multMB,sphMB);
		if (isEventFHM) {
			ProcessTrack(t,D,V0M,sphMB);
			ProcessTrack(t,RC,V0M,sphMB);
			if (isEventJetty)	ProcessTrack(t,D,V0M,Jetty);
			if (isEventIso)		ProcessTrack(t,D,V0M,Iso);			//
			if (isEventJettyMC)	ProcessTrack(t,RC,V0M,Jetty);		// study tracks also for true sphero, stored under RC
			if (isEventIsoMC)	ProcessTrack(t,RC,V0M,Iso);		} 

		if (isEventMHM) {
			ProcessTrack(t,D,NCharged,sphMB);
			ProcessTrack(t,RC,NCharged,sphMB);
			if (isEventJetty)	ProcessTrack(t,D,NCharged,Jetty);
			if (isEventIso)		ProcessTrack(t,D,NCharged,Iso);
			if (isEventJettyMC)	ProcessTrack(t,RC,NCharged,Jetty);
			if (isEventIsoMC)	ProcessTrack(t,RC,NCharged,Iso);		}

		/*for (int iMu = 0; iMu < isEventFHM+isEventMHM+1; ++iMu) {
		for (int iSph = 0; iSph < isEventJetty+isEventIso+1; ++iSph) {
			ProcessTrack(t,D,iMu+iMu*(!isEventFHM),iSph+iSph*(!isEventJetty));
		}	}	*/

		// RT DETERMINATION
		Double_t dphi = mHandler->DeltaPhi(t.GetPhi(),phiLead);
		if (TMath::Abs(dphi) < TMath::Pi()/3.) continue;
		if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) continue;
		if (!t.IskITSrefit()) continue;
		hNchvLeadPt->Fill(ptLead);
		nChTrans0++;
		if (ptLead < 5.) continue;
		nChTrans++;
	}

	hNchvLeadPt2->Fill(ptLead,nChTrans0);
	Double_t eventRt = 0;
	if (nChTrans>0) {
		hNchTrans->Fill(nChTrans);
		eventRt = (double)nChTrans/RT_DEN;
		hRt->Fill(eventRt);		}

	Bool_t isEventRT = false;
	enum { rt01, rt12, rt23, rt34, rt45, rtsizeof };
	Bool_t isRT[rtsizeof] = {false, false, false, false, false};
	if (eventRt) {
		isEventRT = true;
		for (int iRt = 0; iRt < rtsizeof; ++iRt)	{
			isRT[iRt] = (eventRt > (Double_t)iRt 
				&& eventRt < (Double_t)(iRt+1));	}
	}

	if (isEventRT) hEventType->Fill(EVENTTYPES[10],1);
	for (int iRt = 0; iRt < rtsizeof; ++iRt)	{
		if (isRT[iRt]) hEventType->Fill(EVENTTYPES[11+iRt],1); }

	// <pT> vs RT
	if (isEventRT) hLeadPtvRt->Fill(eventRt,ptLead);



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
					if (isEventJettyMC)		hTrackPt[MC][V0M][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC)		hTrackPt[MC][V0M][Iso]->Fill(p.GetPt()); }

				if (isEventMHM) {
					hTrackPt[MC][NCharged][sphMB]->Fill(p.GetPt());
					if (isEventJettyMC)		hTrackPt[MC][NCharged][Jetty]->Fill(p.GetPt());
					if (isEventIsoMC)		hTrackPt[MC][NCharged][Iso]->Fill(p.GetPt()); }

				if (isEventRT) {
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

					if (isEventFHM) {
						hV0Pt[iSp][MC][V0M][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC)		hV0Pt[iSp][MC][V0M][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC)		hV0Pt[iSp][MC][V0M][Iso]->Fill(p.GetPt());		}

					if (isEventMHM) {
						hV0Pt[iSp][MC][NCharged][sphMB]->Fill(p.GetPt());
						if (isEventJettyMC)		hV0Pt[iSp][MC][NCharged][Jetty]->Fill(p.GetPt());
						if (isEventIsoMC)		hV0Pt[iSp][MC][NCharged][Iso]->Fill(p.GetPt());		}

					if (isEventRT) {
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
						
						ProcessV0(v0,iSp,RC,multMB,sphMB);		

						if (isEventFHM) {
							ProcessV0(v0,iSp,RC,V0M,sphMB);
							if (isEventJettyMC)	ProcessV0(v0,iSp,RC,V0M,Jetty);			// study RC spectra for true spher. ev classification
							if (isEventIsoMC)	ProcessV0(v0,iSp,RC,V0M,Iso);
						}
						
						if (isEventMHM) {
							ProcessV0(v0,iSp,RC,NCharged,sphMB);
							if (isEventJettyMC)	ProcessV0(v0,iSp,RC,NCharged,Jetty);
							if (isEventIsoMC)	ProcessV0(v0,iSp,RC,NCharged,Iso);
						}

						if (isEventRT) {
							ProcessV0(v0,iSp,RC,RT,sphMB);
							for (int iRt = 0; iRt < rtsizeof; ++iRt) {
								if (isRT[iRt])	ProcessV0(v0,iSp,RC,RT,3+iRt);	}
						}
					}
				}
			}
		}

		//for (int iMu = 0; iMu < isEventFHM+isEventMHM+isEventRT+1; ++iMu) {
		//for (int iSph = 0; iSph < isEventJetty+isEventIso+1; ++iSph) {
		for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
			if (IsV0(v0,iSp,D)) {

				ProcessV0(v0,iSp,D,multMB,sphMB);

				if (isEventFHM) {
					ProcessV0(v0,iSp,D,V0M,sphMB);
					if (isEventJetty)	ProcessV0(v0,iSp,D,V0M,Jetty);
					if (isEventIso)		ProcessV0(v0,iSp,D,V0M,Iso); }

				if (isEventMHM) {
					ProcessV0(v0,iSp,D,NCharged,sphMB);
					if (isEventJetty)	ProcessV0(v0,iSp,D,NCharged,Jetty);
					if (isEventIso)		ProcessV0(v0,iSp,D,NCharged,Iso); }

				if (isEventRT) {
					ProcessV0(v0,iSp,D,RT,sphMB);
					for (int iRt = 0; iRt < rtsizeof; ++iRt) {
						if (isRT[iRt])	ProcessV0(v0,iSp,D,RT,3+iRt);	}
				}

				//	ProcessV0(v0,iSp,D,iMu+iMu*(!isEventFHM)*(!isEventRT),iSph+iSph*(!isEventJetty));
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
	if (ev.IsPileupFromSPD())			return false;
	if (!ev.AcceptVertex())				return false;
	if (!ev.CheckFlag())				return false;
	//if (!ev.IsCollisionCandidate())		return false;
	//if (ev.GetCentralityQuality() != 0)	return false;

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

Bool_t MyAnalysisV0::ProcessV0(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph) {
	
	//printf("v0 pt is %f, args are %i , %i , %i \n", v0.GetPt(), Sp, Mu, Sph);

	hV0Pt[Sp][Type][Mu][Sph]->Fill(v0.GetPt());
	hV0Eta[Sp][Type][Mu][Sph]->Fill(v0.GetEta());
	Double_t v0mass[] = {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	hV0IMvPt[Sp][Type][Mu][Sph]->Fill(v0.GetPt(),v0mass[Sp]);

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

	if (v0.GetEta() < cuts::V0_ETA[0]) 	return false;
	if (v0.GetEta() > cuts::V0_ETA[1]) 	return false;
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

Bool_t MyAnalysisV0::SelectV0Daughter(MyTrack &tr) {

	if (tr.GetEta() < cuts::K0S_D_ETA[0])	return false;
	if (tr.GetEta() > cuts::K0S_D_ETA[1])	return false;
	if (TMath::Abs(tr.GetDCApvXY()) < cuts::K0S_D_DCAPVXY)	return false;

	return true;
}

Bool_t MyAnalysisV0::SelectParticle(MyParticle &p) {

	if (p.GetEta() < cuts::V0_ETA[0]) 		return false;
	if (p.GetEta() > cuts::V0_ETA[1]) 		return false;
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


	return true;
}

Bool_t MyAnalysisV0::CreateHistograms() {

	// MONITORS
	hEventMonitor 			= new TH1D("hEventMonitor","; Step; Entries",10,-0.5,9.5);
	hTrackMonitor 			= new TH1D("hTrackMonitor","; Step; Entries",10,-0.5,9.5);
	hV0Monitor  			= new TH1D("hV0Monitor","; Step; Entries",10,-0.5,9.5);
	hParticleMonitor 		= new TH1D("hParticleMonitor","; Step; Entries",10,-0.5,9.5);

	// EVENT INFO HISTOGRAMS
	hEventV0MCentrality		= new TH1D("hEventV0MCentrality","; V0M Centrality; Entries",300,0,150);
	hEventRefMult			= new TH1D("hEventRefMult","; Reference multiplicity; Entries",150,0,150);
	hEventV0MCentvRefMult	= new TH2D("hEventV0MCentvRefMult","; Reference multiplicity; V0M Centrality; Entries"
		,150,0,150,300,0,150);

	hEventSpherocity		= new TH1D("hEventSpherocity","; S_{O}; Entries",400,-0.1,1.1);
	hEventType				= new TH1D("hEventType",";; Counts", NEVENTTYPES, 0, NEVENTTYPES);
	for (int iBin = 1; iBin <= NEVENTTYPES; iBin++) {
		hEventType->GetXaxis()->SetBinLabel(iBin,EVENTTYPES[iBin-1]); }

	hEventTSMCvRC 			= new TH2D("hEventTSMCvRC",";S_{O} with MC particles; S_{O} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSNormMCvRC		= new TH2D("hEventTSNormMCvRC",";S_{O}^{p_{T}=1} with MC particles; S_{O}^{p_{T}=1} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSMCvNorm			= new TH2D("hEventTSMCvNorm",";S_{O} with MC particles; S_{O}^{p_{T}=1} with MC particles", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSRCvNorm			= new TH2D("hEventTSRCvNorm",";S_{O} with RC tracks; S_{O}^{p_{T}=1} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);

	hLeadPhivPt				= new TH2D("hLeadPhivPt","; p_{T} (GeV/#it{c}); #phi", 200, 0., 30., 400, -0.2, 6.4);
	hNchvLeadPt				= new TH1D("hNchvLeadPt","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 200, 0., 30.);
	hNchvLeadPt2			= new TH2D("hNchvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 90, 0., 30.,50,0.,50.);
	hNchTrans				= new TH1D("hNchTrans","; N_ch [trans.]; Entries",100, 0., 100.);
	hRt						= new TH1D("hRt","; R_{T}; Entries",100, 0., 5.);
	hLeadPtvRt				= new TH2D("hLeadPtvRt","; R_{T}; p_{T}^{leading} (GeV/#it{c})", 100, 0., 5., 200, 0., 30.);

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


}

Bool_t MyAnalysisV0::BorrowHistograms() {

	mDirFile = (TDirectoryFile*)mHandler->filehist()->Get("MyAnalysisV0_0");
	//mDirFile->ls();

	// MONITORS
	hEventMonitor 			= (TH1D*)mDirFile->Get("hEventMonitor");
	hTrackMonitor 			= (TH1D*)mDirFile->Get("hTrackMonitor");
	hV0Monitor  			= (TH1D*)mDirFile->Get("hV0Monitor");
	hParticleMonitor 		= (TH1D*)mDirFile->Get("hParticleMonitor");
	hEventType				= (TH1D*)mDirFile->Get("hEventType");

	// EVENT INFO HISTOGRAMS
	hEventV0MCentrality		= (TH1D*)mDirFile->Get("hEventV0MCentrality");
	hEventRefMult			= (TH1D*)mDirFile->Get("hEventRefMult");
	hEventV0MCentvRefMult	= (TH2D*)mDirFile->Get("hEventV0MCentvRefMult");

	hEventSpherocity		= (TH1D*)mDirFile->Get("hEventSpherocity");
	hEventTSMCvRC			= (TH2D*)mDirFile->Get("hEventTSMCvRC");
	hEventTSNormMCvRC		= (TH2D*)mDirFile->Get("hEventTSNormMCvRC");
	hEventTSMCvNorm			= (TH2D*)mDirFile->Get("hEventTSMCvNorm");
	hEventTSRCvNorm			= (TH2D*)mDirFile->Get("hEventTSRCvNorm");

	hLeadPhivPt				= (TH2D*)mDirFile->Get("hLeadPhivPt");
	hNchvLeadPt				= (TH1D*)mDirFile->Get("hNchvLeadPt");
	hNchvLeadPt2			= (TH2D*)mDirFile->Get("hNchvLeadPt2");
	hNchTrans				= (TH1D*)mDirFile->Get("hNchTrans");
	hRt						= (TH1D*)mDirFile->Get("hRt");
	hLeadPtvRt				= (TH2D*)mDirFile->Get("hLeadPtvRt");

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

		//hV0PtFit[iSp][iType][iMu][iSph]		= new TH1D(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
		//	";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		
	} } } }	

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0Feeddown[iSp] 	= (TH1D*)mDirFile->Get(Form("hV0Feeddown_%s",SPECIES[iSp]));
		hV0FeeddownPDG[iSp] = (TH1D*)mDirFile->Get(Form("hV0FeeddownPDG_%s",SPECIES[iSp]));
	}		

}

Int_t MyAnalysisV0::Finish() {
	
	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	TH1D* hLeadPt = hLeadPhivPt->ProjectionX();
	hNchvLeadPt->Divide(hLeadPt);

	if (mFlagMC) DoEfficiency();
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
