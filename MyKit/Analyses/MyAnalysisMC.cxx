#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TH3F.h>
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
#include <TString.h>
#include <TROOT.h>

#include "MyAnalysisMC.h"
#include "MyEvent.h"
#include "MyTrack.h"
#include "MyParticle.h"
#include "MyV0.h"
#include "MyCascade.h"
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

#include "UnfoldNTclass.h"

using std::cout;
using namespace MCconsts;
using namespace RooFit;

ClassImp(MyAnalysisMC)

MyAnalysisMC::MyAnalysisMC() {

}

Int_t MyAnalysisMC::Init() {

	mDirFile->cd();

	mFlagMC = mHandler->GetFlagMC();
	mFlagHist = mHandler->GetFlagHist();

	printf("Initialising analysis %s with flag MC %i and Hist %i \n", 
		this->GetName(), mFlagMC, mFlagHist);

	TH1::SetDefaultSumw2(1);
	TH1::AddDirectory(kTRUE);
	if (mFlagHist)	{
		TakeoverHistograms("MyAnalysisMC_0");
		BorrowHistograms();
		mDirFile->cd();	}
	else 			CreateHistograms();

	return 0;
}

Int_t MyAnalysisMC::Make(Int_t iEv) {
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
		case (-2) : hEventType->Fill(EVENTTYPES[2],1);		// MB, rejected by Event Selection
			return 0;
			break;
		case (-1) : hEventType->Fill(EVENTTYPES[3],1);		// MB, passed ES, no vertex
			return 0;
			break;
		case (0) : hEventType->Fill(EVENTTYPES[4],1);		// MB, passed ES, bad vertex
			return 0;
			break;
		case (1) : hEventType->Fill(EVENTTYPES[5],1);		// MB, passed ES, good vertex
			break;
	}

	hEventMonitor->Fill(2);

	// TRACK LOOP TO FIND LEADING
	Int_t nTracks = mHandler->getNtracks();
	Int_t nCascades = mHandler->getNcascades();
	Int_t nParticles = (mFlagMC) ? mHandler->getNparticles() : 0;

	ptLead = -99.; phiLead = -99.; phiPrimeLead = -99;
	Int_t pdgLead = 0;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr)); t.SetHandler(mHandler);

		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;

		// CALCULATING FIELD+CHARGE ADJUSTED FOR PHI FOR GEOMETRICAL CUTS
		Float_t phiPrime = t.GetPhi();//fmod(t.GetPhi(),2*TMath::Pi()/18.);
		if (event.GetMagneticField()<0)		phiPrime = 2*TMath::Pi()-phiPrime;
		if (t.GetSign()<0) phiPrime = 2*TMath::Pi()-phiPrime;
		phiPrime = phiPrime + TMath::Pi()/18.;
		phiPrime = fmod(phiPrime,2*TMath::Pi()/18.);

		if (SelectTrack(t) && IsGeometricalCut(phiPrime, t.GetPt()) ) {				// strict cuts, has primary dca cut
			if (t.GetPt() > ptLead) {		// search for leading track among primaries
				ptLead = t.GetPt();
				phiLead = t.GetPhi();
				phiPrimeLead = phiPrime;
				pdgLead = t.GetMCPdgCode();	}
		}
	}
	hLeadPhivPt->Fill(ptLead,phiLead);
	hLeadPhiPrimevPt->Fill(ptLead,phiPrimeLead);
	hLeadPDG->Fill(pdgLead);


	// RT NCH CALCULATION
	nChTrans = 0;
	Int_t nChTrans0 = 0;	// NT in events also wo leading
	Int_t nChTransA = 0;	// NT in the first transverse region
	Int_t nChTransB = 0;	// NT in the other transverse region
	Double_t meanPtTrans0 = 0;
	Double_t meanPtTransA = 0;
	Double_t meanPtTransB = 0;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{

		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr)); t.SetHandler(mHandler);
		
		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;
		
		// RT DETERMINATION
		if (!IsTrans(t.GetPhi(),phiLead))	continue;
		if (t.GetPt()<0.15) continue;

		// using hybrid tracks
		if (!t.IsITSTPC2011() && !t.IsITSTPC2011HybridNone()) continue;
		
		// APPLY DCA CUT TO AVOID V0 DAUGHTERS
		if (TMath::Abs(t.GetDCApvXY()) > cuts::V0_D_DCAPVXY) continue;
		
		hNchvLeadPt->Fill(ptLead);
		nChTrans0++;
		meanPtTrans0 += t.GetPt();
		if (WhatRegionSide(t.GetPhi(),phiLead)==1) {
			nChTransA++;
			meanPtTransA += t.GetPt(); }
		if (WhatRegionSide(t.GetPhi(),phiLead)==2) {
			nChTransB++;
			meanPtTransB += t.GetPt(); }

		if (ptLead < 5.) continue;
		if (ptLead > 40.) continue;
		nChTrans++;		// increment only if a leading particle was found	
	}
	
	// DETERMINE NT_MIN AND NT_MAX
	Bool_t isSideAMin = 0;			// We need to keep track for also classifying the regions for spectra.
	Int_t nChTransMin = nChTransB;	// NT in the min transverse region
	Int_t nChTransMax = nChTransA;	// NT in the max transverse region
	Double_t meanPtTransMin = meanPtTransB;	// NT in the min transverse region
	Double_t meanPtTransMax = meanPtTransA;	// NT in the max transverse region
	if (nChTransA < nChTransB) {
		isSideAMin = 1;
		nChTransMin = nChTransA;
		nChTransMax = nChTransB;	
		meanPtTransMin = meanPtTransA;
		meanPtTransMax = meanPtTransB;
	}
	if (nChTransA == nChTransB && meanPtTransA < meanPtTransB) {
		isSideAMin = 1;
		nChTransMin = nChTransA;
		nChTransMax = nChTransB;	
		meanPtTransMin = meanPtTransA;
		meanPtTransMax = meanPtTransB;	
	}

	meanPtTrans0 = (nChTrans0>0) ? meanPtTrans0/nChTrans0 : 0;
	meanPtTransMin = (nChTransMin>0) ? meanPtTransMin/nChTransMin : 0;
	meanPtTransMax = (nChTransMax>0) ? meanPtTransMax/nChTransMax : 0;

	hNchvLeadPt2->Fill(ptLead,nChTrans0);
	hNchMinvLeadPt2->Fill(ptLead,nChTransMin);
	hNchMaxvLeadPt2->Fill(ptLead,nChTransMax);
	if (meanPtTrans0>0) hMeanPtvLeadPt2->Fill(ptLead,meanPtTrans0);
	if (meanPtTransMin>0) hMeanPtMinvLeadPt2->Fill(ptLead,meanPtTransMin);
	if (meanPtTransMax>0) hMeanPtMaxvLeadPt2->Fill(ptLead,meanPtTransMax);
	eventRt = 0;
	if (ptLead>5. && ptLead < 40.) {
		hNchTrans->Fill(nChTrans);
	}

	// CLASSIFYING EVENT FOR RT ANALYSIS
	Bool_t isEventRT = false;		
	if (ptLead>5. && ptLead<40.) {
		isEventRT = true;
	}

	if (isEventRT) {
		hNtvNtMin->Fill(nChTransMin,nChTrans);
		hNtvNtMax->Fill(nChTransMax,nChTrans);
		hNtMaxvNtMin->Fill(nChTransMin,nChTransMax);
	}


	// <pT> vs RT
	hLeadPtvNchTrans0->Fill(nChTrans0,ptLead);
	hLeadPtvNchTransMin->Fill(nChTransMin,ptLead);
	hLeadPtvNchTransMax->Fill(nChTransMax,ptLead);
	if (isEventRT) hLeadPtvNchTrans->Fill(nChTrans,ptLead);


	//////////////////////////////////////////////////////////////////
	// MC RT CALCULATION
	//////////////////////////////////////////////////////////////////
	Double_t ptLeadMC = -99., phiLeadMC = -99.;
	Double_t eventRtMC = 0;
	Int_t nChTransMC = 0;
	Int_t nChTransAMC = 0;
	Int_t nChTransBMC = 0;
	Bool_t isSideAMinMC = 0;
	Int_t nChTransMinMC = 0;
	Int_t nChTransMaxMC = 0;
	Double_t meanPtTrans0MC = 0;
	Double_t meanPtTransAMC = 0;
	Double_t meanPtTransBMC = 0;
	
	if (mFlagMC) {

		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.IsPrimary() && 
					( TMath::Abs(p.GetPdgCode()) == 211 || TMath::Abs(p.GetPdgCode()) == 321 
					|| TMath::Abs(p.GetPdgCode()) == 2212 || TMath::Abs(p.GetPdgCode()) == 11 
					|| TMath::Abs(p.GetPdgCode()) == 3312 )
				)	{

				if (p.GetPt() > ptLeadMC) {
					ptLeadMC = p.GetPt();
					phiLeadMC = p.GetPhi();
				}
			}	
		}

		// RT
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.GetPt() > 0.15 && p.IsPrimary())	{

				if (!IsTrans(p.GetPhi(),phiLeadMC)) continue;

				// Exclude Xi from NTMC
				if (TMath::Abs(p.GetPdgCode()) == 3312) continue;

				nChTransMC++;

				if (WhatRegionSide(p.GetPhi(),phiLeadMC)==1) {
					nChTransAMC++;
					meanPtTransAMC += p.GetPt();
				}
				if (WhatRegionSide(p.GetPhi(),phiLeadMC)==2) {
					nChTransBMC++;
					meanPtTransBMC += p.GetPt();
				}


				if (ptLeadMC<5.0) continue;
				if (ptLeadMC>40.0) continue;

			}
		}

		isSideAMinMC = 0;			// We need to keep track for also classifying the regions for spectra.
		nChTransMinMC = nChTransBMC;	// NT in the min transverse region
		nChTransMaxMC = nChTransAMC;	// NT in the max transverse region
		if (nChTransAMC < nChTransBMC) {
			isSideAMinMC = 1;
			nChTransMinMC = nChTransAMC;
			nChTransMaxMC = nChTransBMC;	}

		if (nChTransAMC == nChTransBMC && meanPtTransAMC < meanPtTransBMC) {
			isSideAMinMC = 1;
			nChTransMinMC = nChTransAMC;
			nChTransMaxMC = nChTransBMC;	}

		if ( (ptLeadMC > 5. && ptLeadMC < 40.) && (ptLead > 5. && ptLead < 40.) )	{   // or ||
			hNchTransRC->Fill(nChTrans);
			hNchTransMC->Fill(nChTransMC);
			hNchTransRCvMC->Fill(nChTransMC,nChTrans);
			hNchTransMinRC->Fill(nChTransMin);
			hNchTransMinMC->Fill(nChTransMinMC);
			hNchTransMinRCvMC->Fill(nChTransMinMC,nChTransMin);
			hNchTransMaxRC->Fill(nChTransMax);
			hNchTransMaxMC->Fill(nChTransMaxMC);
			hNchTransMaxRCvMC->Fill(nChTransMaxMC,nChTransMax);
		}

	}

	// CLASSIFYING EVENT FOR RT ANALYSIS
	Bool_t isEventRTMC = false;		
	if (ptLeadMC>5. && ptLeadMC<40.) {
		isEventRTMC = true;
		hNchTransMCTrigMC->Fill(nChTransMC);
		hNchTransMinMCTrigMC->Fill(nChTransMinMC);
		hNchTransMaxMCTrigMC->Fill(nChTransMaxMC);
	}

	/////////////////////////////////////////////////////////////////
	// PID ANALYSIS

	// MC PID PARTICLE LOOP
	std::vector<Int_t> PartLabels;
	std::vector<Int_t> PartIds; 
	if (mFlagMC) {

		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);

			// piKp
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.IsPrimary() && p.GetSign() != 0 &&
				(TMath::Abs(p.GetPdgCode())==211 ||
				 TMath::Abs(p.GetPdgCode())==321 || 
				 TMath::Abs(p.GetPdgCode())==2212 ) ) {

				hPIDPt[piKp][MC]->Fill(p.GetPt());
				if (isEventRTMC) {
	
					Int_t region = WhatRegion(p.GetPhi(),phiLeadMC);
	
					hPIDPtNt[piKp][MC][region]->Fill(p.GetPt(),nChTransMC);
					hPIDPtNtMin[piKp][MC][region]->Fill(p.GetPt(),nChTransMinMC);
					hPIDPtNtMax[piKp][MC][region]->Fill(p.GetPt(),nChTransMaxMC);
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {

						hPIDPtNt[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
					}
				}
			}

			// p
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.IsPrimary() && p.GetSign() != 0 &&
				(TMath::Abs(p.GetPdgCode())==2212 ) ) {

				hPIDPt[pr][MC]->Fill(p.GetPt());
				if (isEventRTMC) {
	
					Int_t region = WhatRegion(p.GetPhi(),phiLeadMC);
	
					hPIDPtNt[pr][MC][region]->Fill(p.GetPt(),nChTransMC);
					hPIDPtNtMin[pr][MC][region]->Fill(p.GetPt(),nChTransMinMC);
					hPIDPtNtMax[pr][MC][region]->Fill(p.GetPt(),nChTransMaxMC);
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {

						hPIDPtNt[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
					}
				}
			}

			// Xi
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.IsPrimary() && p.GetSign() != 0 &&
				(TMath::Abs(p.GetPdgCode())==3312) ) {

				hPIDPt[XiInc][MC]->Fill(p.GetPt());
				hPIDPt[Xi][MC]->Fill(p.GetPt());
				if (isEventRTMC) {
				
					Int_t region = WhatRegion(p.GetPhi(),phiLeadMC);
	

					hPIDPtNt[XiInc][MC][region]->Fill(p.GetPt(),nChTransMC);
					hPIDPtNtMin[XiInc][MC][region]->Fill(p.GetPt(),nChTransMinMC);
					hPIDPtNtMax[XiInc][MC][region]->Fill(p.GetPt(),nChTransMaxMC);
					hPIDPtNt[Xi][MC][region]->Fill(p.GetPt(),nChTransMC);
					hPIDPtNtMin[Xi][MC][region]->Fill(p.GetPt(),nChTransMinMC);
					hPIDPtNtMax[Xi][MC][region]->Fill(p.GetPt(),nChTransMaxMC);
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {
						hPIDPtNt[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						hPIDPtNt[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
					}
				}
			}

			// inc. phi
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1] &&
				(TMath::Abs(p.GetPdgCode())==333) ) {

				hPIDPt[phiInc][MC]->Fill(p.GetPt());
				if (isEventRTMC) {
				
					Int_t region = WhatRegion(p.GetPhi(),phiLeadMC);
	
					hPIDPtNt[phiInc][MC][region]->Fill(p.GetPt(),nChTransMC);
					hPIDPtNtMin[phiInc][MC][region]->Fill(p.GetPt(),nChTransMinMC);
					hPIDPtNtMax[phiInc][MC][region]->Fill(p.GetPt(),nChTransMaxMC);
					
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {
						hPIDPtNt[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
					}
				}
			}

			// phi->K+K-
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 &&
				(TMath::Abs(p.GetPdgCode())==321) &&
				(TMath::Abs(p.GetMotherPdgCode())==333) ) {

				TLorentzVector p1Mom;
				TLorentzVector p2Mom;
				TLorentzVector phiMom;

				

				for (int iP2 = iP+1; iP2 < nParticles; ++iP2)		{
			
					if (!mHandler->particle(iP2)) continue;
					MyParticle p2(mHandler->particle(iP2)); p2.SetHandler(mHandler); p2.SetLabel(iP2);

					if (p2.GetSign() == 0) continue;
					if (p2.GetEta() < cuts::V0_ETA[0] || p2.GetEta() > cuts::V0_ETA[1]) continue;
					if (p2.GetSign() == p.GetSign()) continue;
					if (TMath::Abs(p2.GetPdgCode())!=321) continue;
					if (TMath::Abs(p2.GetMotherPdgCode())!=333) continue;

					p1Mom.SetPtEtaPhiM(p.GetPt(),p.GetEta(),p.GetPhi(),0.493677);
					p2Mom.SetPtEtaPhiM(p2.GetPt(),p2.GetEta(),p2.GetPhi(),0.493677);
					phiMom = p1Mom+p2Mom;

					if (TMath::Abs(phiMom.M()-1.0194550)>0.01) continue;
					if (phiMom.Eta() < cuts::V0_ETA[0] || phiMom.Eta() > cuts::V0_ETA[1]) continue;

					hPIDPt[phi][MC]->Fill(phiMom.Pt());
			
					if (!isEventRTMC) continue;
				
					Int_t region = WhatRegion(FlipNegativeAngle(phiMom.Phi()),phiLeadMC);
					Int_t regiond1 = WhatRegion(p.GetPhi(),phiLeadMC);
					Int_t regiond2 = WhatRegion(p2.GetPhi(),phiLeadMC);

					hPhiDaughterDPhiPt[MC]->Fill(phiMom.Pt(),mHandler->DeltaPhi(FlipNegativeAngle(phiMom.Phi()),p.GetPhi()));

					hPhiDaughterRegionsPt[MC]->Fill(phiMom.Pt(),
							(region == regiond1 && region == regiond2) ? 0 :
							( (region == regiond1 || region == regiond2) ? 1 : 2) );

					hPIDPtNt[phi][MC][region]->Fill(phiMom.Pt(),nChTransMC);
					hPIDPtNtMin[phi][MC][region]->Fill(phiMom.Pt(),nChTransMinMC);
					hPIDPtNtMax[phi][MC][region]->Fill(phiMom.Pt(),nChTransMaxMC);
					
					Int_t regionSide = WhatRegionSide(FlipNegativeAngle(phiMom.Phi()),phiLeadMC);
					if (!region) {
						hPIDPtNt[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMC);
						hPIDPtNtMin[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMinMC);
						hPIDPtNtMax[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMaxMC);
					}
				}
			}
		}
	}

	// RC PID LOOP -- PRIMARY TRACKS
	for (int iTr = 0; iTr < nTracks; ++iTr)		{

		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));	t.SetHandler(mHandler);

		if (mFlagMC) {
			if (t.GetEta() > cuts::V0_ETA[0] && t.GetEta() < cuts::V0_ETA[1]
				&& t.GetSign() != 0 && t.IsMCPrimary() 
				&& (t.IsITSTPC2011() || t.IsITSTPC2011HybridNone()) 
				&& TMath::Abs(t.GetDCApvXY()) < cuts::V0_D_DCAPVXY)	{

				if (TMath::Abs(t.GetMCPdgCode()) != 211 &&
					TMath::Abs(t.GetMCPdgCode()) != 321 &&
					TMath::Abs(t.GetMCPdgCode()) != 2212) continue;
				
				hPIDPt[piKp][RC]->Fill(t.GetPt());

				if (isEventRT)	{

					hPIDDPhivNchTrans[piKp]->Fill(nChTrans,mHandler->DeltaPhi(phiLead,t.GetPhi()));

					Int_t region = WhatRegion(t.GetPhi(),phiLead);
					hPIDPtNt[piKp][RC][region]->Fill(t.GetPt(),nChTrans);
					hPIDPtNtMin[piKp][RC][region]->Fill(t.GetPt(),nChTransMin);
					hPIDPtNtMax[piKp][RC][region]->Fill(t.GetPt(),nChTransMax);
					
					Int_t regionSide = WhatRegionSide(t.GetPhi(),phiLead);
					if (!region) {
						hPIDPtNt[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTrans);
						hPIDPtNtMin[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMin);
						hPIDPtNtMax[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMax);
					}
				}

				if (TMath::Abs(t.GetMCPdgCode()) != 2212) continue;
				
				hPIDPt[pr][RC]->Fill(t.GetPt());

				if (isEventRT)	{

					hPIDDPhivNchTrans[pr]->Fill(nChTrans,mHandler->DeltaPhi(phiLead,t.GetPhi()));

					Int_t region = WhatRegion(t.GetPhi(),phiLead);
					hPIDPtNt[pr][RC][region]->Fill(t.GetPt(),nChTrans);
					hPIDPtNtMin[pr][RC][region]->Fill(t.GetPt(),nChTransMin);
					hPIDPtNtMax[pr][RC][region]->Fill(t.GetPt(),nChTransMax);
					
					Int_t regionSide = WhatRegionSide(t.GetPhi(),phiLead);
					if (!region) {
						hPIDPtNt[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTrans);
						hPIDPtNtMin[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMin);
						hPIDPtNtMax[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMax);
					}
				}

			}
		}
	}

	// RC PID LOOP -- CASCADES
	for (int iCa = 0; iCa < nCascades; ++iCa)		{

		if (!mHandler->cascade(iCa)) continue;
		MyCascade cas(mHandler->cascade(iCa));	cas.SetHandler(mHandler);

		if (mFlagMC) {

			if (cas.GetEta() > cuts::V0_ETA[0] && cas.GetEta() < cuts::V0_ETA[1] 
				&& (TMath::Abs(cas.GetMCPdgCode()) == 3312) ) {

				// ALL XI
				hPIDPt[XiInc][RC]->Fill(cas.GetPt());
				hXiBachDCAXY->Fill(TMath::Abs(cas.GetBachDCApvXY()) );
				hXiPrDCAXY->Fill( cas.GetMCPdgCode() > 0 ? TMath::Abs(cas.GetV0posDCApvXY()) : TMath::Abs(cas.GetV0negDCApvXY()) );

				if (isEventRT)	{

					hPIDDPhivNchTrans[XiInc]->Fill(nChTrans,mHandler->DeltaPhi(phiLead,cas.GetPhi()));

					Int_t region = WhatRegion(cas.GetPhi(),phiLead);
					hPIDPtNt[XiInc][RC][region]->Fill(cas.GetPt(),nChTrans);
					hPIDPtNtMin[XiInc][RC][region]->Fill(cas.GetPt(),nChTransMin);
					hPIDPtNtMax[XiInc][RC][region]->Fill(cas.GetPt(),nChTransMax);
					
					Int_t regionSide = WhatRegionSide(cas.GetPhi(),phiLead);
					if (!region) {
						hPIDPtNt[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTrans);
						hPIDPtNtMin[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMin);
						hPIDPtNtMax[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMax);
					}
				}

				// ONLY XI W DCA CUT	
				if (TMath::Abs(cas.GetBachDCApvXY()) > cuts::V0_D_DCAPVXY
				&& TMath::Abs(cas.GetV0posDCApvXY()) > cuts::V0_D_DCAPVXY
				&& TMath::Abs(cas.GetV0negDCApvXY()) > cuts::V0_D_DCAPVXY ) {
				 
					hPIDPt[Xi][RC]->Fill(cas.GetPt());

					if (isEventRT)	{

						hPIDDPhivNchTrans[Xi]->Fill(nChTrans,mHandler->DeltaPhi(phiLead,cas.GetPhi()));

						Int_t region = WhatRegion(cas.GetPhi(),phiLead);
						hPIDPtNt[Xi][RC][region]->Fill(cas.GetPt(),nChTrans);
						hPIDPtNtMin[Xi][RC][region]->Fill(cas.GetPt(),nChTransMin);
						hPIDPtNtMax[Xi][RC][region]->Fill(cas.GetPt(),nChTransMax);
						
						Int_t regionSide = WhatRegionSide(cas.GetPhi(),phiLead);
						if (!region) {
							hPIDPtNt[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTrans);
							hPIDPtNtMin[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMin);
							hPIDPtNtMax[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMax);
						}
					}
				}
			}
		}
	}

	// RC PID LOOP -- PHI
	for (int iTr = 0; iTr < nTracks; ++iTr)		{

		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));	t.SetHandler(mHandler);

		if (mFlagMC) {
			if (t.GetEta() > cuts::V0_ETA[0] && t.GetEta() < cuts::V0_ETA[1]
				&& t.GetSign() != 0 && t.IsMCPrimary() 
				&& (t.IsITSTPC2011() || t.IsITSTPC2011HybridNone()) 
				&& TMath::Abs(t.GetDCApvXY()) < cuts::V0_D_DCAPVXY )	{

				if (TMath::Abs(t.GetMCPdgCode()) != 321 ) continue;
				if (TMath::Abs(t.GetMCMotherPdgCode()) != 333 ) continue;

				TLorentzVector p1Mom;
				TLorentzVector p2Mom;
				TLorentzVector phiMom;

				for (int iTr2 = iTr+1; iTr2 < nTracks; ++iTr2)		{

					if (!mHandler->track(iTr2)) continue;
					MyTrack t2(mHandler->track(iTr2));	t2.SetHandler(mHandler);

					if (t2.GetEta() > cuts::V0_ETA[0] && t2.GetEta() < cuts::V0_ETA[1]
					&& t2.GetSign() != 0 && t2.IsMCPrimary() 
					&& (t2.IsITSTPC2011() || t2.IsITSTPC2011HybridNone()) 
					&& TMath::Abs(t2.GetDCApvXY()) < cuts::V0_D_DCAPVXY )	{

						if (t2.GetMCPdgCode() != -1*t.GetMCPdgCode() ) continue;
						if (t2.GetMCMotherLabel() != t.GetMCMotherLabel() ) continue;
						if (TMath::Abs(t2.GetMCMotherPdgCode()) != 333 ) continue;

						p1Mom.SetPtEtaPhiM(t.GetPt(),t.GetEta(),t.GetPhi(),0.493677);
						p2Mom.SetPtEtaPhiM(t2.GetPt(),t2.GetEta(),t2.GetPhi(),0.493677);
						phiMom = p1Mom+p2Mom;

						if (phiMom.Eta() < cuts::V0_ETA[0] || phiMom.Eta() > cuts::V0_ETA[1]) continue;

						hPIDPt[phiInc][RC]->Fill(phiMom.Pt());
						hPIDPt[phi][RC]->Fill(phiMom.Pt());
						
						if (!isEventRT) continue;

						hPIDDPhivNchTrans[phi]->Fill(nChTrans,mHandler->DeltaPhi(phiLead,phiMom.Phi()));
					
						Int_t region = WhatRegion(FlipNegativeAngle(phiMom.Phi()),phiLead);
						Int_t regiond1 = WhatRegion(t.GetPhi(),phiLead);
						Int_t regiond2 = WhatRegion(t2.GetPhi(),phiLead);

						hPhiDaughterDPhiPt[RC]->Fill(phiMom.Pt(),mHandler->DeltaPhi(FlipNegativeAngle(phiMom.Phi()),t.GetPhi()));

						hPhiDaughterRegionsPt[RC]->Fill(phiMom.Pt(),
							(region == regiond1 && region == regiond2) ? 0 :
							( (region == regiond1 || region == regiond2) ? 1 : 2) );

						hPIDPtNt[phiInc][RC][region]->Fill(phiMom.Pt(),nChTrans);
						hPIDPtNtMin[phiInc][RC][region]->Fill(phiMom.Pt(),nChTransMin);
						hPIDPtNtMax[phiInc][RC][region]->Fill(phiMom.Pt(),nChTransMax);
						hPIDPtNt[phi][RC][region]->Fill(phiMom.Pt(),nChTrans);
						hPIDPtNtMin[phi][RC][region]->Fill(phiMom.Pt(),nChTransMin);
						hPIDPtNtMax[phi][RC][region]->Fill(phiMom.Pt(),nChTransMax);
						
						Int_t regionSide = WhatRegionSide(FlipNegativeAngle(phiMom.Phi()),phiLead);
						if (!region) {
							hPIDPtNt[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTrans);
							hPIDPtNtMin[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMin);
							hPIDPtNtMax[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMax);
							hPIDPtNt[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTrans);
							hPIDPtNtMin[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMin);
							hPIDPtNtMax[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMax);
						}
					}
				}
			}
		}
	}


	hEventMonitor->Fill(7);
	return 0;	
}

Int_t MyAnalysisMC::ClassifyEvent(MyEvent &event, Int_t ntracks)
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

  	/* event.IsGoodAliEvent() does:
	if (!fIsCollisionCandidate) return kFALSE;
	if (fCentralityQuality != 0) return kFALSE;
	if (!(fIsEventSelected & AliVEvent::kINT7)) return kFALSE;
	if((fEventFlags&fgFlagToCheck)!=fgFlagToCheck) return kFALSE;
	return kTRUE;
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

  return 1;
}


Bool_t MyAnalysisMC::IsTrans(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	if (TMath::Abs(dphi) < TMath::Pi()/3.) 		return false;
	if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return false;
	
	return true;
}

Int_t MyAnalysisMC::WhatRegion(Double_t phi1, Double_t phiTrig) {
	// 0: transverse, 1: toward, 2: away
	
	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	
	if (TMath::Abs(dphi) < TMath::Pi()/3.)			return 1;
	else if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return 2;
	else return 0;
}

Int_t MyAnalysisMC::WhatRegionSide(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	if (!IsTrans(phi1,phiTrig)) return 0;
	
	if (dphi>0) return 1; // region A
	if (dphi<0) return 2; // region B
	return 0;
}

Int_t MyAnalysisMC::IsMinOrMax(Bool_t isA, Int_t reg) {

	Int_t res = 4;
	if (isA && reg == 1)	res = 3;
	if (!isA && reg == 1)	res = 4;
	if (isA && reg == 2)	res = 4;
	if (!isA && reg == 2)	res = 3;
	return res;
}

Bool_t MyAnalysisMC::SelectParticle(MyParticle &p) {

	

	return true;
}

Bool_t MyAnalysisMC::SelectTrack(MyTrack &tr) {

	if (tr.GetEta() < cuts::V0_ETA[0]) 		return false;
	if (tr.GetEta() > cuts::V0_ETA[1]) 		return false;
	if (TMath::Abs(tr.GetDCApvXY()) > 
		cuts::TR_PRIMARY_PAR[0] + 
		cuts::TR_PRIMARY_PAR[1]/TMath::Power(tr.GetPt(),cuts::TR_PRIMARY_PAR[2])) return false;

	if (!tr.IsITSTPC2011() && !tr.IsITSTPC2011HybridNone())					return false; // using hybrid tracks
	
	return true;
}

Bool_t MyAnalysisMC::IsGeometricalCut(Float_t phiprime, Float_t pt) {

	bool ret = true;
	if (pt>0) if ( phiprime<(0.12/pt + TMath::Pi()/18. + 0.035) &&
				phiprime>(0.1/pt/pt + TMath::Pi()/18. - 0.025) ) ret = false;
	
	return ret;
}

Double_t MyAnalysisMC::FlipNegativeAngle(Double_t phi) {
	
	return (phi>0.) ? phi : 2.*TMath::Pi()+phi;
}

Bool_t MyAnalysisMC::CreateHistograms() {

	Int_t nType = (mFlagMC) ? NTYPE : 1;

	// MONITORS
	hEventMonitor 			= new TH1D("hEventMonitor","; Step; Entries",10,-0.5,9.5);
	hEventVz				= new TH1D("hEventVz","; vz; Entries",400,-50,50);

	// EVENT INFO HISTOGRAMS
	hEventType					= new TH1D("hEventType",";; Counts", NEVENTTYPES, 0, NEVENTTYPES);
	for (int iBin = 1; iBin <= NEVENTTYPES; iBin++) {
		hEventType->GetXaxis()->SetBinLabel(iBin,EVENTTYPES[iBin-1]); }
	
	hLeadPhivPt				= new TH2F("hLeadPhivPt","; p_{T} (GeV/#it{c}); #phi", 200, 0., 30., 400, -0.2, 6.4);
	hLeadPhiPrimevPt		= new TH2F("hLeadPhiPrimevPt","; p_{T} (GeV/#it{c}); #phi", 200, 0., 30., 400, -0.2, 6.4);
	hLeadPDG				= new TH1F("hLeadPDG",";PDG ID;Entries",20000,-10000,10000);
	hNchvLeadPt				= new TH1F("hNchvLeadPt","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 200, 0., 30.);
	hNchvLeadPt2			= new TH2F("hNchvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 90, 0., 30.,50,-0.5,49.5);
	hNchMinvLeadPt2			= new TH2F("hNchMinvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.,min]", 90, 0., 30.,50,-0.5,49.5);
	hNchMaxvLeadPt2			= new TH2F("hNchMaxvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.,max]", 90, 0., 30.,50,-0.5,49.5);
	hMeanPtvLeadPt2			= new TH2F("hMeanPtvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); <pT> [trans.]", 90, 0., 30.,50,0.,2.5);
	hMeanPtMinvLeadPt2		= new TH2F("hMeanPtMinvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); <pT> [trans.,min]", 90, 0., 30.,50,0.,2.5);
	hMeanPtMaxvLeadPt2		= new TH2F("hMeanPtMaxvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); <pT> [trans.,max]", 90, 0., 30.,50,0.,2.5);
	hNchTrans				= new TH1F("hNchTrans","; N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransRC 			= new TH1F("hNchTransRC","; RC N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransMC 			= new TH1F("hNchTransMC","; MC N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransMCTrigMC		= new TH1F("hNchTransMCTrigMC","; MC N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransRCvMC			= new TH2F("hNchTransRCvMC", ";MC N_ch [trans.]; RC N_ch [trans.]",50,-0.5,49.5,50,-0.5,49.5);
	hNchTransMinRC 			= new TH1F("hNchTransMinRC","; RC N_ch [trans.,min]; Entries",50, -0.5, 49.5);
	hNchTransMinMC 			= new TH1F("hNchTransMinMC","; MC N_ch [trans.,min]; Entries",50, -0.5, 49.5);
	hNchTransMinMCTrigMC		= new TH1F("hNchTransMinMCTrigMC","; MC N_ch [trans.,min]; Entries",50, -0.5, 49.5);
	hNchTransMinRCvMC		= new TH2F("hNchTransMinRCvMC", ";MC N_ch [trans.,min]; RC N_ch [trans.,min]",50,-0.5,49.5,50,-0.5,49.5);
	hNchTransMaxRC 			= new TH1F("hNchTransMaxRC","; RC N_ch [trans.,max]; Entries",50, -0.5, 49.5);
	hNchTransMaxMC 			= new TH1F("hNchTransMaxMC","; MC N_ch [trans.,max]; Entries",50, -0.5, 49.5);
	hNchTransMaxMCTrigMC		= new TH1F("hNchTransMaxMCTrigMC","; MC N_ch [trans.,max]; Entries",50, -0.5, 49.5);
	hNchTransMaxRCvMC		= new TH2F("hNchTransMaxRCvMC", ";MC N_ch [trans.,max]; RC N_ch [trans.,max]",50,-0.5,49.5,50,-0.5,49.5);
	hNtvNtMin				= new TH2F("hNtvNtMin","; N_ch [trans.,min]; N_ch [trans.]",50, -0.5, 49.5,50, -0.5, 49.5);
	hNtvNtMax				= new TH2F("hNtvNtMax","; N_ch [trans.,max]; N_ch [trans.]",50, -0.5, 49.5,50, -0.5, 49.5);
	hNtMaxvNtMin			= new TH2F("hNtMaxvNtMin","; N_ch [trans.,min]; N_ch [trans.,max]",50, -0.5, 49.5,50, -0.5, 49.5);
	hLeadPtvNchTrans		= new TH2F("hLeadPtvNchTrans","; N_{ch}^{trans}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hLeadPtvNchTrans0		= new TH2F("hLeadPtvNchTrans0","; N_{ch}^{trans}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hLeadPtvNchTransMin		= new TH2F("hLeadPtvNchTransMin","; N_{ch}^{trans,min}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hLeadPtvNchTransMax		= new TH2F("hLeadPtvNchTransMax","; N_{ch}^{trans,max}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);

 	hPhiDaughterRegionsPt[MC]   = new TH2F("hPhiDaughterRegionsPt_MC",";p_{T} (GeV/#it{c}); n. of differences between m. and d. regions", NPTBINS[phi], XBINS[phi], 3,-0.5,2.5);
 	hPhiDaughterRegionsPt[RC]	= new TH2F("hPhiDaughterRegionsPt_RC",";p_{T} (GeV/#it{c}); n. of differences between m. and d. regions", NPTBINS[phi], XBINS[phi], 3,-0.5,2.5);

 	hPhiDaughterDPhiPt[MC]	= new TH2F("hPhiDaughterDPhiPt_MC",";p_{T} (GeV/#it{c}); n. of differences between m. and d. regions", NPTBINS[phi], XBINS[phi], 300, -3.2, 3.2);
 	hPhiDaughterDPhiPt[RC]	= new TH2F("hPhiDaughterDPhiPt_RC",";p_{T} (GeV/#it{c}); n. of differences between m. and d. regions", NPTBINS[phi], XBINS[phi], 300, -3.2, 3.2);

 	hXiBachDCAXY		= new TH1F("hXiBachDCAXY",";DCA_{xy}^{PV}; Entries",	400, 0.0, 1.0);
 	hXiPrDCAXY			= new TH1F("hXiPrDCAXY",";DCA_{xy}^{PV}; Entries",	400, 0.0, 1.0);

 	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{

 		hPIDDPhivNchTrans[iSp]			= new TH2F(TString::Format("hPIDDPhivNchTrans_%s",SPECIES[iSp]),"; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
	
	for (int iType = 0; iType < nType; ++iType)		{
				
		hPIDPt[iSp][iType]			= new TH1D(TString::Format("hPIDPt_%s_%s",SPECIES[iSp],TYPE[iType]),
			";p_{T} (GeV/#it{c}); Entries",	NPTBINS[iSp],XBINS[iSp]);


	} }  

 	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
				
		hPIDPtNt[iSp][iType][iReg]			= new TH2F(TString::Format("hPIDPtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [trans.]; Entries",	NPTBINS[iSp],XBINS[iSp], 50, -0.5, 49.5);
		hPIDPtNtMin[iSp][iType][iReg]			= new TH2F(TString::Format("hPIDPtNtMin_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [trans.,min]; Entries",	NPTBINS[iSp],XBINS[iSp], 50, -0.5, 49.5);
		hPIDPtNtMax[iSp][iType][iReg]			= new TH2F(TString::Format("hPIDPtNtMax_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [trans.,max]; Entries",	NPTBINS[iSp],XBINS[iSp], 50, -0.5, 49.5);
		
	} } } 


	printf("Histograms in analysis %s created \n", 
		this->GetName());

	return true;

}

Bool_t MyAnalysisMC::BorrowHistograms() {

	// Retrieve histograms that are used in Finish()
	// When running on the same file

	Int_t nType = (mFlagMC) ? NTYPE : 1;

	hNchTransMCTrigMC	= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMCTrigMC");

	hNt 			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTrans")->Clone("hNt");
	hNtRec 			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransRC")->Clone("hNtRec");
	hNtGen			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMC")->Clone("hNtGen");
	hNtRM			= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hNchTransRCvMC")->Clone("hNtRM");

	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{

		hPIDEffi[iSp]		= (TH1D*)mHandler->analysis(0)->dirFile()->Get(TString::Format("hPIDEffi_%s",SPECIES[iSp]));

	for (int iType = 0; iType < nType; ++iType)		{
				
		hPIDPt[iSp][iType]	= (TH1D*)mHandler->analysis(0)->dirFile()->Get(TString::Format("hPIDPt_%s_%s",SPECIES[iSp],TYPE[iType]));

	} }  

 	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
				
		hPIDPtNt[iSp][iType][iReg]		= (TH2F*)mHandler->analysis(0)->dirFile()->Get(TString::Format("hPIDPtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hPIDPtNtMin[iSp][iType][iReg]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get(TString::Format("hPIDPtNtMin_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hPIDPtNtMax[iSp][iType][iReg]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get(TString::Format("hPIDPtNtMax_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		
	} } } 


	return true;

}

Int_t MyAnalysisMC::Finish() {
	
	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	CalculateEfficiencies();
	CorrectEfficiency();
	Normalise();
	Unfold();

	printf("Analysis %s finished.\n",this->GetName());
	return 0;	
}

Bool_t MyAnalysisMC::CalculateEfficiencies() {

	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{

		hPIDEffi[iSp] = (TH1D*)hPIDPt[iSp][RC]->Clone(TString::Format("hPIDEffi_%s",SPECIES[iSp]));
		hPIDEffi[iSp]->Divide(hPIDEffi[iSp],hPIDPt[iSp][MC],1.,1.,"B");

		hPIDEffi[iSp]->Write();

	}

	return true;

}

Bool_t MyAnalysisMC::Normalise() {

	// MC
	for (int iSp = 0; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < NTYPE; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{

		printf("Normalising histogram %s by N_T distribution \n", hPIDPtNtCorr[iSp][iType][iReg]->GetName());

		for (int iNt = 1; iNt < hPIDPtNtCorr[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {

			Double_t NormEv = 0;			
			NormEv = (iType==MC) ? hNchTransMCTrigMC->Integral(1,50) : hNt->Integral(1,50);
			TH2F* htmp = hPIDPtNtCorr[iSp][iType][iReg];

			TH1D* hpt = (TH1D*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
		}
	}	}	}

	return true;
}

Bool_t MyAnalysisMC::CorrectEfficiency() {

	// SCALE BIN WIDTH
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
		
		hPIDPtCorr[iSp][iType] = 	
			(TH1D*)hPIDPt[iSp][iType]->Clone(TString::Format("hPIDPtCorr_%s_%s",SPECIES[iSp],TYPE[iType]));
		hPIDPtCorr[iSp][iType]->Scale(1.,"width");
	
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hPIDPtNtCorr[iSp][iType][iReg] = 	
			(TH2F*)hPIDPtNt[iSp][iType][iReg]->Clone(TString::Format("hPIDPtNtCorr_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hPIDPtNtMinCorr[iSp][iType][iReg] = 	
			(TH2F*)hPIDPtNtMin[iSp][iType][iReg]->Clone(TString::Format("hPIDPtNtMinCorr_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hPIDPtNtMaxCorr[iSp][iType][iReg] = 	
			(TH2F*)hPIDPtNtMax[iSp][iType][iReg]->Clone(TString::Format("hPIDPtNtMaxCorr_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		

		hPIDPtNtCorr[iSp][iType][iReg]			= ScaleWidthTH2(hPIDPtNtCorr[iSp][iType][iReg]);
		hPIDPtNtMinCorr[iSp][iType][iReg]		= ScaleWidthTH2(hPIDPtNtMinCorr[iSp][iType][iReg]);
		hPIDPtNtMaxCorr[iSp][iType][iReg]		= ScaleWidthTH2(hPIDPtNtMaxCorr[iSp][iType][iReg]);

	}	}	}


	// CORRECT FOR EFFICIENCY
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{

		hPIDPtCorr[iSp][RC]->Divide(hPIDEffi[iSp]);

	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	
		hPIDPtNtCorr[iSp][RC][iReg] 		= DivideTH2ByTH1(hPIDPtNtCorr[iSp][RC][iReg],hPIDEffi[iSp]);
		hPIDPtNtMinCorr[iSp][RC][iReg] 		= DivideTH2ByTH1(hPIDPtNtMinCorr[iSp][RC][iReg],hPIDEffi[iSp]);
		hPIDPtNtMaxCorr[iSp][RC][iReg] 		= DivideTH2ByTH1(hPIDPtNtMaxCorr[iSp][RC][iReg],hPIDEffi[iSp]);

	}	}


	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
		
		hPIDPtCorr[iSp][iType]->Write();

	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hPIDPtNtCorr[iSp][iType][iReg]->Write();
		hPIDPtNtMinCorr[iSp][iType][iReg]->Write();
		hPIDPtNtMaxCorr[iSp][iType][iReg]->Write();

	}	}	}

	return true;

}

Bool_t MyAnalysisMC::Unfold() {

	mUnf = new UnfoldNTclass();
	cout << "Unfolding class created " << mUnf << endl;

	DoUnfoldingNt();
	DoUnfolding1D();

	return true;

}

TH2F* MyAnalysisMC::FlipMatrix(TH2F* h) {

	TH2F* htmp = (TH2F*)h->Clone(Form("%s_flip",h->GetName()));
	for (int i = 1; i < h->GetNbinsX()+1; i++) {
	for (int j = 1; j < h->GetNbinsY()+1; j++) {
			htmp->SetBinContent(i,j,h->GetBinContent(j,i));
			htmp->SetBinError(i,j,h->GetBinError(j,i));
	}	}
	htmp->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(h->GetXaxis()->GetTitle());
	delete h;
	return htmp; 
}

TH2F* MyAnalysisMC::ScaleWidthTH2(TH2F* h) {

	if (!h) {
		printf("No histogram to scale \n");
		return 0x0; }

	for (int iX = 1; iX < h->GetNbinsX()+1; iX++) {
		for (int iY = 1; iY < h->GetNbinsY()+1; iY++) {	
			float binwidth = h->GetXaxis()->GetBinWidth(iX);
			float binc = h->GetBinContent(iX,iY);
			float bine = h->GetBinError(iX,iY);
			if (binwidth>0) h->SetBinContent(iX,iY,binc/binwidth);
			if (binwidth>0) h->SetBinError(iX,iY,bine/binwidth);
		}
	}

	return h;
}

TH2F* MyAnalysisMC::DivideTH2ByTH1(TH2F* h, TH1D* d) {

	if (!h || !d) {
		printf("Input parameters empty \n");
		return 0x0; }

	for (int iY = 1; iY < h->GetNbinsY()+1; iY++) {

		TH1F* hx = (TH1F*)h->ProjectionX("",iY,iY);
		hx->Divide(d);

		for (int iX = 1; iX < h->GetNbinsX()+1; iX++) {
			h->SetBinContent(iX,iY,hx->GetBinContent(iX));
			h->SetBinError(iX,iY,hx->GetBinError(iX));
		}

		delete hx;
	}

	return h;
}

void MyAnalysisMC::DoUnfoldingNt() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding";
	const bool eRM = false;
	const int NumberOfIters = 20;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtRM = FlipMatrix(hNtRM);
	// ROWWISE NORMALISATION IS PERFORMED IN mUnf
	//cout << "Histograms: " << hNtRec << " " << hNt << " " << hNtGen << " " << hNtRM << "\n";
	
	if (eRM) mUnf->ExtrapolateRM(hNtRM);
	mUnf->SetnIter(NumberOfIters);
	mUnf->SetError(hNtRec);
	mUnf->SetError(hNtGen);
	mUnf->SaveSolutionNT(kTRUE);
	
	if (mHandler->GetFlagMC()) mUnf->Setup(hNtRec,hNtGen,hNtRM);
	else mUnf->Setup(hNt,hNtGen,hNtRM);
	
	mUnf->Unfold();
	mUnf->V2H();
	
	printf(" - Unfolding region : %s\n",mUnf->GetRegion());

	hNtUnf = (TH1F*)(mUnf->GetUnfoldedDistH())->Clone("hNtUnf");
	if(mHandler->GetFlagMC())	hNtClosure = (TH1F*)(mUnf->GetClosureH())->Clone("hNtClosure");
	
	TH1F* hNT = (TH1F*)hNtUnf->Clone("_hNT");
	hRtUnf = (TH1F*)mUnf->RebinNT2RT(hNtUnf, kTRUE);
	hRtUnf->Scale(1.0/hRtUnf->Integral());
	
	//! Relative statistical uncertainty
	//! NT distribution with final bins to plot
	TH1F* hRelStatUnc = (TH1F*)hNtUnf->Clone("hRelStatUnc_Nt");
	hRelStatUnc->Reset();

	TH1F* hNch = (TH1F*)hNtUnf->Clone("hNT");
	hNch->Reset();

	for(int bin = 1; bin <= hRelStatUnc->GetNbinsX(); bin++){
		double yield = hNtUnf->GetBinContent(bin);
		double error = hNtUnf->GetBinError(bin);
		if( yield > 0. ) hRelStatUnc->SetBinContent(bin, error / yield);
	}

	lOut->Add(hNtRM);
	
	if (mHandler->GetFlagMC()) 	lOut->Add(hNtClosure);
	//else 	ComparisonPublished(hRtUnf);

	lOut->Add(hNtUnf);
	lOut->Add(hRtUnf);
	lOut->Add(hRelStatUnc);

	TGraph* gChi2 = nullptr;
	if (mHandler->GetFlagMC())	{
		mUnf->DrawNchClosure(hNtGen,hNtUnf,-0.5,30.0,"#it{N}_{T}","Unfolded/True",dOut);
		TVectorD iterX = (TVectorD)mUnf->GetChi2Vectors(kTRUE);
		TVectorD Chi2 = (TVectorD)mUnf->GetChi2Vectors(kFALSE);
		gChi2 = new TGraph(iterX,Chi2);
		gChi2->SetName("gChi2");
		gChi2->SetMarkerStyle(8);
		lOut->Add(gChi2);
	}

	TFile* fUnfOut;
	if (mHandler->GetFlagMC())	fUnfOut = new TFile(Form("./%s/1D_newClass_mc.root",dOut),"RECREATE");
	else	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");


	fUnfOut->cd();
	lOut->Write();
	(mUnf->GetObjArray())->Write();
	fUnfOut->Close();
	delete fUnfOut;

	mDirFile->cd();
	lOut->Write();
}

void MyAnalysisMC::DoUnfolding1D() {

	mHandler->root()->SetBatch(kTRUE);

	const char* dOut = "results_unfolding";
	const bool eRM = false;
	const int NumberOfIters = 20;
	const char* Regions[NREGIONS] = {"Trans1D","Toward","Away","TransMin1D","TransMax1D"};

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	Int_t iType = RC;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(7+iSp);
		obj->SetMCAnalysis(mHandler->GetFlagMC());
		obj->SaveSolutionNT(kFALSE);

		printf(" - Unfolding species %s vs NT in region : %s\n",SPECIES[iSp],obj->GetRegion());
		cout << "RM has " << hNtRM->GetNbinsX() << " x " << hNtRM->GetNbinsY() << endl;
		cout << "hMC has " << hPIDPtNtCorr[iSp][iType][iReg]->GetNbinsX() << endl;
		obj->UnfoldV02D(hNtRM,hPIDPtNtCorr[iSp][iType][iReg],hPIDPtNtCorr[iSp][iType][iReg]);
		
		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << hPIDPtNtCorr[iSp][iType][iReg] << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		// first is Gen, second is Unf

		// APPLY PROPER ERRORS!
		hPIDPtNtCorrUnf[iSp][iType][iReg] = (TH2F*)Arr->FindObject("_hPtvsNch")->Clone(Form("hPIDPtNtCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		for (int iX = 1; iX < hPIDPtNtCorrUnf[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hPIDPtNtCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float bine = hPIDPtNtCorrUnf[iSp][iType][iReg]->GetBinError(iX,iY);
				if (bine>0) hPIDPtNtCorrUnf[iSp][iType][iReg]->SetBinError(iX,iY,bine);
			}
		}

		if (mHandler->GetFlagMC())	obj->GetMCclosureinRTBins(hPIDPtNtCorr[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	

		// NORMALISE
		/*for (int iNt = 1; iNt < hV0PtNtFitCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {
			Double_t NormEv = 0;
			NormEv = hNtUnf->Integral(1,50) > 0 ? (double)hNtUnf->GetBinContent(iNt) / hNtUnf->Integral(1,50) : 1;
			//cout << "Normalisingggg by " << NormEv << endl;
			TH2F* htmp = hV0PtNtFitCorrUnf[iSp][iType][iReg];

			TH1F* hpt = (TH1F*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
			delete hpt;
		}*/

		lOut->Add(hPIDPtNtCorrUnf[iSp][iType][iReg]);
		fUnfOut1D->cd();
		TDirectory* dir = fUnfOut1D->mkdir(Form("%s_%s",SPECIES[iSp],REGIONS[iReg]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	}

	lOut->Write();
	fUnfOut1D->Close();
	delete fUnfOut1D;
	mDirFile->cd();
	lOut->Write();

	mHandler->root()->SetBatch(kFALSE);
}



