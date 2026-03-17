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
#include "TSystem.h"

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
	Int_t nV0s = mHandler->getNv0s();
	Int_t nCascades = mHandler->getNcascades();
	Int_t nParticles = (mFlagMC) ? mHandler->getNparticles() : 0;

	ptLead = -99.; phiLead = -99.; phiPrimeLead = -99;
	Int_t pdgLead = 0;
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr)); t.SetHandler(mHandler);

		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;

		// OPTIONAL, REQUIRE PDG
		if (TMath::Abs(t.GetMCPdgCode()) != 211 && TMath::Abs(t.GetMCPdgCode()) != 321
			&& TMath::Abs(t.GetMCPdgCode()) != 2212) continue;
		if (t.GetPt()>5.) {
			hLeadMotherPDG->Fill(t.GetMCMotherPdgCode());
			if (t.IsMCPrimary()) hLeadDCAPri->Fill(TMath::Abs(t.GetDCApvXY()));
			else hLeadDCASec->Fill(TMath::Abs(t.GetDCApvXY()));
		}
		if (!t.IsMCPrimary()) continue;

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
	if (ptLead>5. && ptLead < 40.) hLeadPDG->Fill(pdgLead);


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
		if (TMath::Abs(t.GetMCPdgCode()) != 211 && TMath::Abs(t.GetMCPdgCode()) != 321
			&& TMath::Abs(t.GetMCPdgCode()) != 2212) continue;
		
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
					|| TMath::Abs(p.GetPdgCode()) == 2212 )
				)	{

					/*|| TMath::Abs(p.GetPdgCode()) == 11 
					|| TMath::Abs(p.GetPdgCode()) == 3312 || TMath::Abs(p.GetPdgCode()) == 13
					|| TMath::Abs(p.GetPdgCode()) == 3222 || TMath::Abs(p.GetPdgCode()) == 3112
					|| TMath::Abs(p.GetPdgCode()) == 0 )*/

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

				// OPTIONAL -- ONLY PIKP FOR NT
				if (TMath::Abs(p.GetPdgCode()) != 211 && TMath::Abs(p.GetPdgCode()) != 321
					&& TMath::Abs(p.GetPdgCode()) != 2212 ) continue;

				// Exclude Xi from NTMC
				//if (TMath::Abs(p.GetPdgCode()) == 3312) continue;

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
	if (ptLead>5. && ptLead < 40.  && ptLeadMC>5. && ptLeadMC<40.) {// THE LAST TWO CAN BE REMOVED
		hNchTrans->Fill(nChTrans);
	}

	// CLASSIFYING EVENT FOR RT ANALYSIS
	Bool_t isEventRT = false;		
	if (ptLead>5. && ptLead<40.  && ptLeadMC>5. && ptLeadMC<40.) {// THE LAST TWO CAN BE REMOVED
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

	Bool_t isEventRTMC = false;		
	if (ptLeadMC>5. && ptLeadMC<40.  && ptLead>5. && ptLead<40.) {  // THE LAST TWO CAN BE REMOVED
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
					if (isEventRT) {
						hPIDPtNtRCvMC[piKp][MC][region]->Fill(p.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[piKp][MC][region]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[piKp][MC][region]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
					}
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {

						hPIDPtNt[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						if (isEventRT) {
							hPIDPtNtRCvMC[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[piKp][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
						}
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
					if (isEventRT) {
						hPIDPtNtRCvMC[pr][MC][region]->Fill(p.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[pr][MC][region]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[pr][MC][region]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
					}
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {

						hPIDPtNt[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						if (isEventRT) {
							hPIDPtNtRCvMC[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[pr][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
						}
					}
				}
			}

			// k0s
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.IsPrimary() &&
				(TMath::Abs(p.GetPdgCode())==310 ) ) {

				hPIDPt[K0s][MC]->Fill(p.GetPt());
				if (isEventRTMC) {
	
					Int_t region = WhatRegion(p.GetPhi(),phiLeadMC);
	
					hPIDPtNt[K0s][MC][region]->Fill(p.GetPt(),nChTransMC);
					hPIDPtNtMin[K0s][MC][region]->Fill(p.GetPt(),nChTransMinMC);
					hPIDPtNtMax[K0s][MC][region]->Fill(p.GetPt(),nChTransMaxMC);
					if (isEventRT) {
						hPIDPtNtRCvMC[K0s][MC][region]->Fill(p.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[K0s][MC][region]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[K0s][MC][region]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
					}
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {

						hPIDPtNt[K0s][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[K0s][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[K0s][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						if (isEventRT) {
							hPIDPtNtRCvMC[K0s][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[K0s][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[K0s][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
						}
					}
				}
			}

			// lambdas
			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.IsPrimary() &&
				(TMath::Abs(p.GetPdgCode())==3122 ) ) {

				hPIDPt[L][MC]->Fill(p.GetPt());
				if (isEventRTMC) {
	
					Int_t region = WhatRegion(p.GetPhi(),phiLeadMC);
	
					hPIDPtNt[L][MC][region]->Fill(p.GetPt(),nChTransMC);
					hPIDPtNtMin[L][MC][region]->Fill(p.GetPt(),nChTransMinMC);
					hPIDPtNtMax[L][MC][region]->Fill(p.GetPt(),nChTransMaxMC);
					if (isEventRT) {
						hPIDPtNtRCvMC[L][MC][region]->Fill(p.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[L][MC][region]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[L][MC][region]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
					}
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {

						hPIDPtNt[L][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[L][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[L][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						if (isEventRT) {
							hPIDPtNtRCvMC[L][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[L][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[L][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
						}
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
					if (isEventRT) {
						hPIDPtNtRCvMC[XiInc][MC][region]->Fill(p.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[XiInc][MC][region]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[XiInc][MC][region]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
						hPIDPtNtRCvMC[Xi][MC][region]->Fill(p.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[Xi][MC][region]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[Xi][MC][region]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
					}
	
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {
						hPIDPtNt[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						hPIDPtNt[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						if (isEventRT) {
							hPIDPtNtRCvMC[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[XiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
							hPIDPtNtRCvMC[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[Xi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
						}
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
					if (isEventRT) {
						hPIDPtNtRCvMC[phiInc][MC][region]->Fill(p.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[phiInc][MC][region]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[phiInc][MC][region]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
					}
					
					Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLeadMC);
					if (!region) {
						hPIDPtNt[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC);
						hPIDPtNtMin[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						hPIDPtNtMax[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC);
						if (isEventRT) {
							hPIDPtNtRCvMC[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[phiInc][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(p.GetPt(),nChTransMaxMC,nChTransMax);
						}
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
					if (isEventRT) {
						hPIDPtNtRCvMC[phi][MC][region]->Fill(phiMom.Pt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[phi][MC][region]->Fill(phiMom.Pt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[phi][MC][region]->Fill(phiMom.Pt(),nChTransMaxMC,nChTransMax);
					}
					
					Int_t regionSide = WhatRegionSide(FlipNegativeAngle(phiMom.Phi()),phiLeadMC);
					if (!region) {
						hPIDPtNt[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMC);
						hPIDPtNtMin[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMinMC);
						hPIDPtNtMax[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMaxMC);
						if (isEventRT) {
							hPIDPtNtRCvMC[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[phi][MC][IsMinOrMax(isSideAMinMC,regionSide)]->Fill(phiMom.Pt(),nChTransMaxMC,nChTransMax);
						}
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

					if (isEventRTMC) {
						hPIDPtNtRCvMC[piKp][RC][region]->Fill(t.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[piKp][RC][region]->Fill(t.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[piKp][RC][region]->Fill(t.GetPt(),nChTransMaxMC,nChTransMax);
					}
					
					Int_t regionSide = WhatRegionSide(t.GetPhi(),phiLead);
					if (!region) {
						hPIDPtNt[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTrans);
						hPIDPtNtMin[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMin);
						hPIDPtNtMax[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMax);
						if (isEventRTMC) {
							hPIDPtNtRCvMC[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[piKp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMaxMC,nChTransMax);
						}
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

					if (isEventRTMC) {
						hPIDPtNtRCvMC[pr][RC][region]->Fill(t.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[pr][RC][region]->Fill(t.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[pr][RC][region]->Fill(t.GetPt(),nChTransMaxMC,nChTransMax);
					}
					
					Int_t regionSide = WhatRegionSide(t.GetPhi(),phiLead);
					if (!region) {
						hPIDPtNt[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTrans);
						hPIDPtNtMin[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMin);
						hPIDPtNtMax[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMax);
						if (isEventRTMC) {
							hPIDPtNtRCvMC[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[pr][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(t.GetPt(),nChTransMaxMC,nChTransMax);
						}
					}
				}

			}
		}
	}

	// RC PID LOOP -- V0s
	for (int iV0 = 0; iV0 < nV0s; ++iV0)		{

		if (!mHandler->v0(iV0)) continue;
		MyV0 v0(mHandler->v0(iV0));	v0.SetHandler(mHandler);

		if (mFlagMC) {

			if (v0.GetEta() > cuts::V0_ETA[0] && v0.GetEta() < cuts::V0_ETA[1] 
				&& (TMath::Abs(v0.GetMCPdgCode()) == 310 || TMath::Abs(v0.GetMCPdgCode()) == 3122) && v0.IsMCPrimary() ) {

				// ONLY V0 W DCA CUT	
				if ( TMath::Abs(v0.GetPosTrack()->GetImpactParameter(0)) > cuts::V0_D_DCAPVXY
				&& TMath::Abs(v0.GetNegTrack()->GetImpactParameter(0)) > cuts::V0_D_DCAPVXY ) {

					Int_t sp = (TMath::Abs(v0.GetMCPdgCode()) == 310) ? K0s : L; 
				 
					hPIDPt[sp][RC]->Fill(v0.GetPt());

					if (isEventRT)	{

						hPIDDPhivNchTrans[sp]->Fill(nChTrans,mHandler->DeltaPhi(phiLead,v0.GetPhi()));

						Int_t region = WhatRegion(v0.GetPhi(),phiLead);
						hPIDPtNt[sp][RC][region]->Fill(v0.GetPt(),nChTrans);
						hPIDPtNtMin[sp][RC][region]->Fill(v0.GetPt(),nChTransMin);
						hPIDPtNtMax[sp][RC][region]->Fill(v0.GetPt(),nChTransMax);
						if (isEventRTMC) {
							hPIDPtNtRCvMC[sp][RC][region]->Fill(v0.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[sp][RC][region]->Fill(v0.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[sp][RC][region]->Fill(v0.GetPt(),nChTransMaxMC,nChTransMax);
						}	
						
						Int_t regionSide = WhatRegionSide(v0.GetPhi(),phiLead);
						if (!region) {
							hPIDPtNt[sp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(v0.GetPt(),nChTrans);
							hPIDPtNtMin[sp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(v0.GetPt(),nChTransMin);
							hPIDPtNtMax[sp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(v0.GetPt(),nChTransMax);
							if (isEventRTMC) {
								hPIDPtNtRCvMC[sp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(v0.GetPt(),nChTransMC,nChTrans);
								hPIDPtNtMinRCvMC[sp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(v0.GetPt(),nChTransMinMC,nChTransMin);
								hPIDPtNtMaxRCvMC[sp][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(v0.GetPt(),nChTransMaxMC,nChTransMax);
							}
						}
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
				&& (TMath::Abs(cas.GetMCPdgCode()) == 3312) && cas.IsMCPrimary() ) {

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
					if (isEventRTMC) {
						hPIDPtNtRCvMC[XiInc][RC][region]->Fill(cas.GetPt(),nChTransMC,nChTrans);
						hPIDPtNtMinRCvMC[XiInc][RC][region]->Fill(cas.GetPt(),nChTransMinMC,nChTransMin);
						hPIDPtNtMaxRCvMC[XiInc][RC][region]->Fill(cas.GetPt(),nChTransMaxMC,nChTransMax);
					}
					
					Int_t regionSide = WhatRegionSide(cas.GetPhi(),phiLead);
					if (!region) {
						hPIDPtNt[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTrans);
						hPIDPtNtMin[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMin);
						hPIDPtNtMax[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMax);
						if (isEventRTMC) {
							hPIDPtNtRCvMC[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[XiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMaxMC,nChTransMax);
						}
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
						if (isEventRTMC) {
							hPIDPtNtRCvMC[Xi][RC][region]->Fill(cas.GetPt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[Xi][RC][region]->Fill(cas.GetPt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[Xi][RC][region]->Fill(cas.GetPt(),nChTransMaxMC,nChTransMax);
						}	
						
						Int_t regionSide = WhatRegionSide(cas.GetPhi(),phiLead);
						if (!region) {
							hPIDPtNt[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTrans);
							hPIDPtNtMin[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMin);
							hPIDPtNtMax[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMax);
							if (isEventRTMC) {
								hPIDPtNtRCvMC[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMC,nChTrans);
								hPIDPtNtMinRCvMC[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMinMC,nChTransMin);
								hPIDPtNtMaxRCvMC[Xi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(cas.GetPt(),nChTransMaxMC,nChTransMax);
							}
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
						if (isEventRTMC) {
							hPIDPtNtRCvMC[phiInc][RC][region]->Fill(phiMom.Pt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[phiInc][RC][region]->Fill(phiMom.Pt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[phiInc][RC][region]->Fill(phiMom.Pt(),nChTransMaxMC,nChTransMax);
							hPIDPtNtRCvMC[phi][RC][region]->Fill(phiMom.Pt(),nChTransMC,nChTrans);
							hPIDPtNtMinRCvMC[phi][RC][region]->Fill(phiMom.Pt(),nChTransMinMC,nChTransMin);
							hPIDPtNtMaxRCvMC[phi][RC][region]->Fill(phiMom.Pt(),nChTransMaxMC,nChTransMax);
						}
						
						Int_t regionSide = WhatRegionSide(FlipNegativeAngle(phiMom.Phi()),phiLead);
						if (!region) {
							hPIDPtNt[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTrans);
							hPIDPtNtMin[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMin);
							hPIDPtNtMax[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMax);
							hPIDPtNt[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTrans);
							hPIDPtNtMin[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMin);
							hPIDPtNtMax[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMax);
							if (isEventRTMC) {
								hPIDPtNtRCvMC[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMC,nChTrans);
								hPIDPtNtMinRCvMC[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMinMC,nChTransMin);
								hPIDPtNtMaxRCvMC[phiInc][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMaxMC,nChTransMax);
								hPIDPtNtRCvMC[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMC,nChTrans);
								hPIDPtNtMinRCvMC[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMinMC,nChTransMin);
								hPIDPtNtMaxRCvMC[phi][RC][IsMinOrMax(isSideAMin,regionSide)]->Fill(phiMom.Pt(),nChTransMaxMC,nChTransMax);
							}
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

	Int_t nType = NTYPE-1;

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
	hLeadMotherPDG			= new TH1F("hLeadMotherPDG",";PDG ID;Entries",20000,-10000,10000);
	hLeadDCAPri				= new TH1F("hLeadDCAPri",";DCA_{xy}^{PV}; Entries",	400, 0.0, 1.0);
	hLeadDCASec				= new TH1F("hLeadDCASec",";DCA_{xy}^{PV}; Entries",	400, 0.0, 1.0);
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

	const int NNTBINS = 50;
    double NTBINS[NNTBINS + 1];
    for(int i = 0; i < NNTBINS+1; ++i)	NTBINS[i] = (double)i - 0.5;
    
 	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
				
		hPIDPtNt[iSp][iType][iReg]			= new TH2F(TString::Format("hPIDPtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [trans.]; Entries",	NPTBINS[iSp],XBINS[iSp], 50, -0.5, 49.5);
		hPIDPtNtMin[iSp][iType][iReg]			= new TH2F(TString::Format("hPIDPtNtMin_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [trans.,min]; Entries",	NPTBINS[iSp],XBINS[iSp], 50, -0.5, 49.5);
		hPIDPtNtMax[iSp][iType][iReg]			= new TH2F(TString::Format("hPIDPtNtMax_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [trans.,max]; Entries",	NPTBINS[iSp],XBINS[iSp], 50, -0.5, 49.5);
		
		hPIDPtNtRCvMC[iSp][iType][iReg]			= new TH3F(TString::Format("hPIDPtNtRCvMC_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [true,trans.]; N_{ch} [meas,trans.]; Entries",	NPTBINS[iSp], XBINS[iSp], NNTBINS, NTBINS, NNTBINS, NTBINS);
		hPIDPtNtMinRCvMC[iSp][iType][iReg]			= new TH3F(TString::Format("hPIDPtNtMinRCvMC_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [true,trans.,min]; N_{ch} [meas,trans.,min]; Entries",	NPTBINS[iSp],XBINS[iSp], NNTBINS, NTBINS, NNTBINS, NTBINS);
		hPIDPtNtMaxRCvMC[iSp][iType][iReg]			= new TH3F(TString::Format("hPIDPtNtMaxRCvMC_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";p_{T} (GeV/#it{c}); N_{ch} [true,trans.,max]; N_{ch} [meas,trans.,max]; Entries",	NPTBINS[iSp],XBINS[iSp], NNTBINS, NTBINS, NNTBINS, NTBINS);


	} } } 


	printf("Histograms in analysis %s created \n", 
		this->GetName());

	return true;

}

Bool_t MyAnalysisMC::BorrowHistograms() {

	// Retrieve histograms that are used in Finish()
	// When running on the same file

	Int_t nType = NTYPE-1;

	hNchTransMCTrigMC	= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMCTrigMC");

	hNt 			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTrans")->Clone("hNt");
	hNtRec 			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransRC")->Clone("hNtRec");
	hNtGen			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMC")->Clone("hNtGen");
	hNtRM			= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hNchTransRCvMC")->Clone("hNtRM");

	hNtMin			= (TH1F*)((TH2F*)mHandler->analysis(0)->dirFile()->Get("hNtvNtMin"))->ProjectionX("hNtMin");
	hNtMinRec		= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMinRC")->Clone("hNtMinRec");
	hNtMinGen		= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMinMC")->Clone("hNtMinGen");
	hNtMinRM		= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMinRCvMC")->Clone("hNtMinRM");
	
	hNtMax			= (TH1F*)((TH2F*)mHandler->analysis(0)->dirFile()->Get("hNtvNtMax"))->ProjectionX("hNtMax");
	hNtMaxRec		= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMaxRC")->Clone("hNtMaxRec");
	hNtMaxGen		= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMaxMC")->Clone("hNtMaxGen");
	hNtMaxRM		= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMaxRCvMC")->Clone("hNtMaxRM");

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

	if (mFlagHist)	{
		CalculateEfficiencies();
		CorrectEfficiency();
		Normalise();
		LoadDataXi();
		LoadDataPhi();
		RebinPhiForClosure();
		Unfold();
		BinHistogramsIntoRT();
	}

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

TH2F* MyAnalysisMC::RebinTH2(TH2F* h, Int_t nbins, const Double_t xbins[]) {

	if (!h) return 0x0;
	if (h->GetNbinsX() != nbins) {
			
		TAxis *xaxis = h->GetXaxis(); 
		TAxis *yaxis = h->GetYaxis(); 
		TH2F *htmp = new TH2F("htmp",
			TString::Format(";%s;%s;%s",xaxis->GetTitle(),yaxis->GetTitle(),h->GetZaxis()->GetTitle()),
			nbins, xbins, yaxis->GetNbins(), yaxis->GetBinLowEdge(1), yaxis->GetBinLowEdge(h->GetYaxis()->GetNbins()+1));
		for (int j=1; j<=yaxis->GetNbins();j++)	{ 
		for (int i=1; i<=xaxis->GetNbins();i++)	{ 
			htmp->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h->GetBinContent(i,j)); 
		}	}	
			
		h = (TH2F*)htmp->Clone(h->GetName());
		delete htmp;
			
	}
	return h;

}

Bool_t MyAnalysisMC::Normalise() {

	// MC
	for (int iSp = 0; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < NTYPE-1; ++iType)			{
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

			htmp = hPIDPtNtMinCorr[iSp][iType][iReg];
			hpt = (TH1D*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

			htmp = hPIDPtNtMaxCorr[iSp][iType][iReg];
			hpt = (TH1D*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
		}


		//if (iSp == Xi)	hPIDPtNtCorr[iSp][iType][iReg] = RebinTH2(hPIDPtNtCorr[iSp][iType][iReg],NXIPTBINS,XIXBINS);



	}	}	}

	return true;
}

Bool_t MyAnalysisMC::CorrectEfficiency() {

	// SCALE BIN WIDTH
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE-1; ++iType)		{
		
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
	for (int iType = 0; iType < NTYPE-1; ++iType)		{
		
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

	delete mUnf; mUnf = new UnfoldNTclass();
	DoUnfoldingNtMin();
	DoUnfolding1DMin();

	delete mUnf; mUnf = new UnfoldNTclass();
	DoUnfoldingNtMax();
	DoUnfolding1DMax();

	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
	}

	return true;

}

TH2F* MyAnalysisMC::FlipMatrix(TH2F* h, Int_t limit = 0) {

    cout << "Flipping " << h << " " << h->GetName() << endl;

    const int nBinsX = h->GetNbinsX();
    const int nBinsY = h->GetNbinsY();

    // Dynamically allocate arrays for bin edges
    double* xBinEdges = new double[nBinsX + 1];
    double* yBinEdges = new double[nBinsY + 1];

    // Extract the bin edges from the x-axis
    for (int i = 0; i <= nBinsX; i++) {
        xBinEdges[i] = h->GetXaxis()->GetBinLowEdge(i + 1);
    }

    // Extract the bin edges from the y-axis
    for (int j = 0; j <= nBinsY; j++) {
        yBinEdges[j] = h->GetYaxis()->GetBinLowEdge(j + 1);
    }

    // Create a new histogram with the flipped axis bins
    TH2F* htmp = new TH2F(Form("%s_flip", h->GetName()), h->GetTitle(),
                          nBinsY, yBinEdges, nBinsX, xBinEdges);

    // Fill the new histogram with flipped content
    for (int i = 1; i <= nBinsY; i++) {
        for (int j = 1; j <= nBinsX; j++) {
            htmp->SetBinContent(i, j, h->GetBinContent(j, i) > limit ? h->GetBinContent(j, i) : 0);
            htmp->SetBinError(i, j, h->GetBinContent(j, i) > limit ? h->GetBinError(j, i) : 0);
        }
    }

    // Cleanup: Delete dynamically allocated arrays
    delete[] xBinEdges;
    delete[] yBinEdges;
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

TH2F* MyAnalysisMC::ScaleTH2Rows(TH2F* h, Double_t d) {

	if (!h) {
		printf("Input parameters empty \n");
		return 0x0; }

	for (int iY = 1; iY < h->GetNbinsY()+1; iY++) {

		TH1F* hx = (TH1F*)h->ProjectionX("",iY,iY);
		hx->Scale(d);

		for (int iX = 1; iX < h->GetNbinsX()+1; iX++) {
			h->SetBinContent(iX,iY,hx->GetBinContent(iX));
			h->SetBinError(iX,iY,hx->GetBinError(iX));
		}

		delete hx;
	}

	return h;
}

void MyAnalysisMC::LoadDataXi() {

	TFile* fileXi = new TFile("../files/xi_results_oli_jan_30.root","READ");
	const char* xiNames[NREGIONS] = { "Transverse", "Towards", "Away", "Transverse", "Transverse"};  

	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hPIDPtNtCorr[Xi][D][iReg] = 0x0;
		hPIDPtNtMinCorr[Xi][D][iReg] = 0x0;
		hPIDPtNtMaxCorr[Xi][D][iReg] = 0x0;
		
		hPIDPtNtCorr[Xi][D][iReg] = (TH2F*)fileXi->Get(Form("hXiNtVsPt%s",xiNames[iReg]))->Clone(TString::Format("hPIDPtNtCorr_%s_%s_%s",SPECIES[Xi],TYPE[D],REGIONS[iReg]));
		hPIDPtNtMinCorr[Xi][D][iReg] = (TH2F*)fileXi->Get(Form("hXiNtMinVsPt%s",xiNames[iReg]))->Clone(TString::Format("hPIDPtNtMinCorr_%s_%s_%s",SPECIES[Xi],TYPE[D],REGIONS[iReg]));
		hPIDPtNtMaxCorr[Xi][D][iReg] = (TH2F*)fileXi->Get(Form("hXiNtMaxVsPt%s",xiNames[iReg]))->Clone(TString::Format("hPIDPtNtMaxCorr_%s_%s_%s",SPECIES[Xi],TYPE[D],REGIONS[iReg]));

		cout << "iReg " << hPIDPtNtCorr[Xi][D][iReg] << " " << hPIDPtNtMinCorr[Xi][D][iReg] << " " << hPIDPtNtMaxCorr[Xi][D][iReg] << endl;
		cout << "iReg " << hPIDPtNtCorr[Xi][D][iReg]->GetName() << " " << hPIDPtNtMinCorr[Xi][D][iReg]->GetName() << " " << hPIDPtNtMaxCorr[Xi][D][iReg]->GetName() << endl;
		hPIDPtNtCorr[Xi][D][iReg]->SetDirectory(mDirFile);
		hPIDPtNtMinCorr[Xi][D][iReg]->SetDirectory(mDirFile);
		hPIDPtNtMaxCorr[Xi][D][iReg]->SetDirectory(mDirFile);

	}

	delete fileXi;
	mDirFile->cd();


	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		// This isn't a correct way of flipping (results in equal bins) but works with the unfolding class
		hPIDPtNtCorr[Xi][D][iReg] = (TH2F*)FlipMatrix(hPIDPtNtCorr[Xi][D][iReg]);
		hPIDPtNtMinCorr[Xi][D][iReg] = (TH2F*)FlipMatrix(hPIDPtNtMinCorr[Xi][D][iReg]);
		hPIDPtNtMaxCorr[Xi][D][iReg] = (TH2F*)FlipMatrix(hPIDPtNtMaxCorr[Xi][D][iReg]);

		hPIDPtNtCorr[Xi][D][iReg]->Write();
		hPIDPtNtMinCorr[Xi][D][iReg]->Write();
		hPIDPtNtMaxCorr[Xi][D][iReg]->Write();
	}

	// UNNORMALISE XI DATA
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	Int_t iSp = Xi;
	Int_t iType = D;

		for (int iNt = 1; iNt < 51; ++iNt)	{
			Double_t NormScaleNt = hNt->Integral(iNt,iNt)/hNt->Integral(1,51);
			if (iReg==3) NormScaleNt = hNtMin->Integral(iNt,iNt)/hNtMin->Integral(1,51);
			if (iReg==4) NormScaleNt = hNtMax->Integral(iNt,iNt)/hNtMax->Integral(1,51);

			cout << "iReg " << iReg << " " << hPIDPtNtCorr[Xi][D][iReg] << " " << hPIDPtNtMinCorr[Xi][D][iReg] << " " << hPIDPtNtMaxCorr[Xi][D][iReg] << endl;

			for (int iX = 1; iX < hPIDPtNtCorr[Xi][D][iReg]->GetNbinsX()+1; iX++) {
				Double_t binc = hPIDPtNtCorr[Xi][D][iReg]->GetBinContent(iX,iNt);
				Double_t bine = hPIDPtNtCorr[Xi][D][iReg]->GetBinError(iX,iNt);
				hPIDPtNtCorr[Xi][D][iReg]->SetBinContent(iX,iNt,binc*NormScaleNt);
				hPIDPtNtCorr[Xi][D][iReg]->SetBinError(iX,iNt,bine*NormScaleNt);

				if (iReg==3) {
					Double_t binc = hPIDPtNtMinCorr[Xi][D][iReg]->GetBinContent(iX,iNt);
					Double_t bine = hPIDPtNtMinCorr[Xi][D][iReg]->GetBinError(iX,iNt);
					hPIDPtNtMinCorr[Xi][D][iReg]->SetBinContent(iX,iNt,binc*NormScaleNt);
					hPIDPtNtMinCorr[Xi][D][iReg]->SetBinError(iX,iNt,bine*NormScaleNt);
				}
				if (iReg==4) {
				Double_t binc = hPIDPtNtMaxCorr[Xi][D][iReg]->GetBinContent(iX,iNt);
				Double_t bine = hPIDPtNtMaxCorr[Xi][D][iReg]->GetBinError(iX,iNt);
				hPIDPtNtMaxCorr[Xi][D][iReg]->SetBinContent(iX,iNt,binc*NormScaleNt);
				hPIDPtNtMaxCorr[Xi][D][iReg]->SetBinError(iX,iNt,bine*NormScaleNt);
				}

			}
		}
	}
}


void MyAnalysisMC::LoadDataPhi() {

	TFile* filePhi = new TFile("../files/CompressedHistograms.root","READ");
	const char* phiNames[NREGIONS] = { "Transverse", "Toward", "Away", "Transverse", "Transverse"};  

	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hPIDPtNtCorr[phi][D][iReg] = 0x0;
		hPIDPtNtMinCorr[phi][D][iReg] = 0x0;
		hPIDPtNtMaxCorr[phi][D][iReg] = 0x0;
		
		hPIDPtNtCorr[phi][D][iReg] = (TH2F*)filePhi->Get(Form("hPhi_Vs_NT_%s",phiNames[iReg]))->Clone(TString::Format("hPIDPtNtCorr_%s_%s_%s",SPECIES[phi],TYPE[D],REGIONS[iReg]));
		//hPIDPtNtMinCorr[phi][D][iReg] = (TH2F*)fileXi->Get(Form("hXiNtMinVsPt%s",xiNames[iReg]))->Clone(TString::Format("hPIDPtNtMinCorr_%s_%s_%s",SPECIES[Xi],TYPE[D],REGIONS[iReg]));
		//hPIDPtNtMaxCorr[phi][D][iReg] = (TH2F*)fileXi->Get(Form("hXiNtMaxVsPt%s",xiNames[iReg]))->Clone(TString::Format("hPIDPtNtMaxCorr_%s_%s_%s",SPECIES[Xi],TYPE[D],REGIONS[iReg]));

		hPIDPtNtCorr[Xi][D][iReg]->SetDirectory(mDirFile);
		//hPIDPtNtMinCorr[Xi][D][iReg]->SetDirectory(mDirFile);
		//hPIDPtNtMaxCorr[Xi][D][iReg]->SetDirectory(mDirFile);

	}

	hPhiPtNtBinning = (TH2F*)hPIDPtNtCorr[phi][D][0]->Clone("hPhiPtNtBinning");

	delete filePhi;
	mDirFile->cd();


	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		//hPIDPtNtCorr[phi][D][iReg] = (TH2F*)FlipMatrix(hPIDPtNtCorr[phi][D][iReg]);
		//hPIDPtNtMinCorr[Xi][D][iReg] = (TH2F*)FlipMatrix(hPIDPtNtMinCorr[Xi][D][iReg]);
		//hPIDPtNtMaxCorr[Xi][D][iReg] = (TH2F*)FlipMatrix(hPIDPtNtMaxCorr[Xi][D][iReg]);

		hPIDPtNtCorr[phi][D][iReg]->Write();
		//hPIDPtNtMinCorr[Xi][D][iReg]->Write();
		//hPIDPtNtMaxCorr[Xi][D][iReg]->Write();
	}

	// UNNORMALISE PHI DATA
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	Int_t iSp = phi;
	Int_t iType = D;

	TH1F* hNtRebin = RebinHistogram1D(hNt,hPIDPtNtCorr[phi][D][iReg]);
	hNtRebin->Write();

		for (int iNt = 1; iNt < hNtRebin->GetNbinsX()+1; ++iNt)	{
			Double_t NormScaleNt = hNtRebin->Integral(iNt,iNt)/hNtRebin->Integral(1,51);
			if (iReg==3) NormScaleNt = hNtMin->Integral(iNt,iNt)/hNtMin->Integral(1,51);
			if (iReg==4) NormScaleNt = hNtMax->Integral(iNt,iNt)/hNtMax->Integral(1,51);

			//cout << "iReg " << iReg << " " << hPIDPtNtCorr[Xi][D][iReg] << " " << hPIDPtNtMinCorr[Xi][D][iReg] << " " << hPIDPtNtMaxCorr[Xi][D][iReg] << endl;

			for (int iX = 1; iX < hPIDPtNtCorr[phi][D][iReg]->GetNbinsX()+1; iX++) {
				Double_t binc = hPIDPtNtCorr[phi][D][iReg]->GetBinContent(iX,iNt);
				Double_t bine = hPIDPtNtCorr[phi][D][iReg]->GetBinError(iX,iNt);
				hPIDPtNtCorr[phi][D][iReg]->SetBinContent(iX,iNt,binc*NormScaleNt);
				hPIDPtNtCorr[phi][D][iReg]->SetBinError(iX,iNt,bine*NormScaleNt);

				/*
				if (iReg==3) {
					Double_t binc = hPIDPtNtMinCorr[Xi][D][iReg]->GetBinContent(iX,iNt);
					Double_t bine = hPIDPtNtMinCorr[Xi][D][iReg]->GetBinError(iX,iNt);
					hPIDPtNtMinCorr[Xi][D][iReg]->SetBinContent(iX,iNt,binc*NormScaleNt);
					hPIDPtNtMinCorr[Xi][D][iReg]->SetBinError(iX,iNt,bine*NormScaleNt);
				}
				if (iReg==4) {
				Double_t binc = hPIDPtNtMaxCorr[Xi][D][iReg]->GetBinContent(iX,iNt);
				Double_t bine = hPIDPtNtMaxCorr[Xi][D][iReg]->GetBinError(iX,iNt);
				hPIDPtNtMaxCorr[Xi][D][iReg]->SetBinContent(iX,iNt,binc*NormScaleNt);
				hPIDPtNtMaxCorr[Xi][D][iReg]->SetBinError(iX,iNt,bine*NormScaleNt);
				}*/

			}
		}

		// Now we have phi as a function of RT classes but in NT binning
		TH2F* hRebinned = new TH2F(Form("%s_rebin",hPIDPtNtCorr[phi][D][iReg]->GetName()),hPIDPtNtCorr[phi][D][iReg]->GetTitle(),
			NPHIPTBINS, PHIXBINS, 50, -0.5, 49.5);
		for (int iX = 1; iX < hPIDPtNtCorr[phi][D][iReg]->GetNbinsX()+1; iX++) {
		for (int iY = 1; iY < hPIDPtNtCorr[phi][D][iReg]->GetNbinsY()+1; iY++) {
			hRebinned->SetBinContent(iX,iY,hPIDPtNtCorr[phi][D][iReg]->GetBinContent(iX,iY));
			hRebinned->SetBinError(iX,iY,hPIDPtNtCorr[phi][D][iReg]->GetBinError(iX,iY));
		}	}
		hPIDPtNtCorr[phi][D][iReg] = hRebinned;
		hPIDPtNtCorr[phi][D][iReg]->Write();
	}
}

void MyAnalysisMC::RebinPhiForClosure() {


	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	//for (int iType = 0; iType < NTYPE-1; ++iType)		{ // full rebinning of matrix
	for (int iType = 0; iType < 1; ++iType)		{		// or only partial in n_meas

		// First store NT projections for reweighting of the unfolding matrix
		hNtPID[phi][iType][iReg] = (TH1F*)hPIDPtNtCorr[phi][iType][iReg]->ProjectionY(Form("hNtPID_%i_%i_%i",phi,iType,iReg),1,51);
		hNtPID[phi][iType][iReg]->Write();

		// Rebinning NT into coarser according to the NT binning given by histogram from data
		hPIDPtNtCorr[phi][iType][iReg] = RebinHistogram2DOnlyNT(hPIDPtNtCorr[phi][iType][iReg],hPhiPtNtBinning);

		// Moving back to histogram with equal bins so that it works w the unf. class
		TH2F* hRebinned = new TH2F(Form("%s_rebin",hPIDPtNtCorr[phi][iType][iReg]->GetName()),hPIDPtNtCorr[phi][iType][iReg]->GetTitle(),
			NPHIPTBINS, PHIXBINS, 50, -0.5, 49.5);
		for (int iX = 1; iX < hPIDPtNtCorr[phi][iType][iReg]->GetNbinsX()+1; iX++) {
		for (int iY = 1; iY < hPIDPtNtCorr[phi][iType][iReg]->GetNbinsY()+1; iY++) {
			hRebinned->SetBinContent(iX,iY,hPIDPtNtCorr[phi][iType][iReg]->GetBinContent(iX,iY));
			hRebinned->SetBinError(iX,iY,hPIDPtNtCorr[phi][iType][iReg]->GetBinError(iX,iY));
		}	}
		hPIDPtNtCorr[phi][iType][iReg] = hRebinned;
		hPIDPtNtCorr[phi][iType][iReg]->Write();


	}	}

}

TH2F* MyAnalysisMC::RebinHistogram2D(TH2F* hrm, TH2F* hptnt) {

  	const int maxBins = 100; // Maximum expected number of bins
    double newBinEdges[maxBins + 1]; // +1 for the upper edge of the last bin
    int numBins = hptnt->GetYaxis()->GetNbins();

    // Check to not exceed maxBins
    if (numBins > maxBins) {
        // Handle error: Too many bins
        return nullptr;
    }

    newBinEdges[0] = -0.5;
    // Extract bin edges from the y-axis of hptnt
    for (int iBin = 2; iBin <= numBins + 1; ++iBin) { // Include the upper edge of the last bin
        newBinEdges[iBin - 1] = 0.5 + hptnt->GetYaxis()->GetBinLowEdge(iBin);
    }

    // Create a new histogram with these bin edges for both x and y axes
    TH2F* hRebinned = new TH2F(Form("%s_rebin",hrm->GetName()), hrm->GetTitle(),
                               numBins, newBinEdges, // x-axis bins
                               numBins, newBinEdges); // y-axis bins

    // Fill the new histogram with contents from hrm, adjusting for the new binning
    for (int xBin = 1; xBin <= hrm->GetNbinsX(); ++xBin) {
        for (int yBin = 1; yBin <= hrm->GetNbinsY(); ++yBin) {

            double content = hrm->GetBinContent(xBin, yBin);
            double xCenter = hrm->GetXaxis()->GetBinCenter(xBin);
            double yCenter = hrm->GetYaxis()->GetBinCenter(yBin);
            if (content != 0) {
                // Fill the histogram using the center of the original bins
                hRebinned->Fill(xCenter, yCenter, content);
            }
        }
    }

    return hRebinned;
}

TH2F* MyAnalysisMC::RebinHistogram2DOnlyNT(TH2F* hrm, TH2F* hptnt) {

  	const int maxBins = 100; // Maximum expected number of bins
    double newBinEdges[maxBins + 1]; // +1 for the upper edge of the last bin
    int numBins = hptnt->GetYaxis()->GetNbins();

    // Check to not exceed maxBins
    if (numBins > maxBins) {
        // Handle error: Too many bins
        return nullptr;
    }

    newBinEdges[0] = -0.5;
    // Extract bin edges from the y-axis of hptnt
    for (int iBin = 2; iBin <= numBins + 1; ++iBin) { // Include the upper edge of the last bin
        newBinEdges[iBin - 1] = 0.5 + hptnt->GetYaxis()->GetBinLowEdge(iBin);
    }

    // Create a new histogram with these bin edges for both x and y axes
    TH2F* hRebinned = new TH2F(Form("%s_rebin",hrm->GetName()), hrm->GetTitle(),
                               NPHIPTBINS,PHIXBINS, // x-axis bins
                               numBins, newBinEdges); // y-axis bins

    // Fill the new histogram with contents from hrm, adjusting for the new binning
    for (int xBin = 1; xBin <= hrm->GetNbinsX(); ++xBin) {
        for (int yBin = 1; yBin <= hrm->GetNbinsY(); ++yBin) {

            double content = hrm->GetBinContent(xBin, yBin);
            double xCenter = hrm->GetXaxis()->GetBinCenter(xBin);
            double yCenter = hrm->GetYaxis()->GetBinCenter(yBin);
            if (content != 0) {
                // Fill the histogram using the center of the original bins
                hRebinned->Fill(xCenter, yCenter, content);
            }
        }
    }

    return hRebinned;
}

TH1F* MyAnalysisMC::RebinHistogram1D(TH1F* hrm, TH2F* hptnt) {

  	const int maxBins = 200; // Maximum expected number of bins
    double newBinEdges[maxBins + 1]; // +1 for the upper edge of the last bin
    int numBins = hptnt->GetYaxis()->GetNbins();

    // Check to not exceed maxBins
    if (numBins > maxBins) {
        // Handle error: Too many bins
        return nullptr;
    }

    newBinEdges[0] = -0.5;
    // Extract bin edges from the y-axis of hptnt
    for (int iBin = 2; iBin <= numBins + 1; ++iBin) { // Include the upper edge of the last bin
        newBinEdges[iBin - 1] = 0.5 + hptnt->GetYaxis()->GetBinLowEdge(iBin);
    }

    // Create a new histogram with these bin edges for both x and y axes
    TH1F* hRebinned = new TH1F(Form("%s_rebin",hrm->GetName()), hrm->GetTitle(),
                               numBins, newBinEdges);

    // Fill the new histogram with contents from hrm, adjusting for the new binning
    for (int xBin = 1; xBin <= hrm->GetNbinsX(); ++xBin) {

        double content = hrm->GetBinContent(xBin);
        double xCenter = hrm->GetXaxis()->GetBinCenter(xBin);
        if (content != 0) {
                // Fill the histogram using the center of the original bins
                hRebinned->Fill(xCenter, content);
        }
    }

    return hRebinned;
}


/*void MyAnalysisMC::DoUnfoldingNtRebin() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding_rebin";
	const bool eRM = false;
	const int NumberOfIters = 3;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtRM = FlipMatrix(hNtRM);
	TH2F* hNtRMRebin = RebinHistogram2D(hNtRM,hPIDPtNtCorr[phi][D][0]);
	TH1F* hNtRecRebin = RebinHistogram1D(hNtRec,hPIDPtNtCorr[phi][D][0]);
	TH1F* hNtGenRebin = RebinHistogram1D(hNtGen,hPIDPtNtCorr[phi][D][0]);
	TH1F* hNtRebin = RebinHistogram1D(hNt,hPIDPtNtCorr[phi][D][0]);
	// ROWWISE NORMALISATION IS PERFORMED IN mUnf
	//cout << "Histograms: " << hNtRec << " " << hNt << " " << hNtGen << " " << hNtRM << "\n";
	
	if (eRM) mUnf->ExtrapolateRM(hNtRMRebin);
	mUnf->SetnIter(NumberOfIters);
	mUnf->SetError(hNtRecRebin);
	mUnf->SetError(hNtGenRebin);
	mUnf->SaveSolutionNT(kTRUE);
	
	if (mHandler->GetFlagMC()) mUnf->Setup(hNtRecRebin,hNtGenRebin,hNtRMRebin);
	else mUnf->Setup(hNtRebin,hNtGenRebin,hNtRMRebin);
	
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

	lOut->Add(hNtRMRebin);
	
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

	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");
	fUnfOut->cd();
	lOut->Write();
	(mUnf->GetObjArray())->Write();
	fUnfOut->Close();
	delete fUnfOut;	

	mDirFile->cd();
	lOut->Write();
}*/


void MyAnalysisMC::DoUnfoldingNt() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding";
	const bool eRM = false;
	const int NumberOfIters = 5;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtRM = FlipMatrix(hNtRM,20);
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

	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");
	fUnfOut->cd();
	lOut->Write();
	(mUnf->GetObjArray())->Write();
	fUnfOut->Close();
	delete fUnfOut;	

	mDirFile->cd();
	lOut->Write();
}

void MyAnalysisMC::DoUnfoldingNtMin() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding_min";
	const bool eRM = false;
	const int NumberOfIters = 5;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtMinRM = FlipMatrix(hNtMinRM,20);
	// ROWWISE NORMALISATION IS PERFORMED IN mUnf
	//cout << "Histograms: " << hNtRec << " " << hNt << " " << hNtGen << " " << hNtRM << "\n";
	
	if (eRM) mUnf->ExtrapolateRM(hNtMinRM);
	mUnf->SetnIter(NumberOfIters);
	mUnf->SetError(hNtMinRec);
	mUnf->SetError(hNtMinGen);
	mUnf->SaveSolutionNT(kTRUE);
	
	if (mHandler->GetFlagMC()) mUnf->Setup(hNtMinRec,hNtMinGen,hNtMinRM);
	else mUnf->Setup(hNtMin,hNtMinGen,hNtMinRM);
	
	mUnf->Unfold();
	mUnf->V2H();
	
	printf(" - Unfolding region : %s_min\n",mUnf->GetRegion());

	hNtMinUnf = (TH1F*)(mUnf->GetUnfoldedDistH())->Clone("hNtMinUnf");
	if(mHandler->GetFlagMC())	hNtMinClosure = (TH1F*)(mUnf->GetClosureH())->Clone("hNtMinClosure");
	
	TH1F* hNTMin = (TH1F*)hNtMinUnf->Clone("_hNTMin");
	hRtMinUnf = (TH1F*)mUnf->RebinNT2RT(hNtMinUnf, kTRUE);
	hRtMinUnf->Scale(1.0/hRtMinUnf->Integral());
	
	//! Relative statistical uncertainty
	//! NT distribution with final bins to plot
	TH1F* hRelStatUnc = (TH1F*)hNtUnf->Clone("hRelStatUnc_NtMin");
	hRelStatUnc->Reset();

	TH1F* hNch = (TH1F*)hNtMinUnf->Clone("hNTMin");
	hNch->Reset();

	for(int bin = 1; bin <= hRelStatUnc->GetNbinsX(); bin++){
		double yield = hNtMinUnf->GetBinContent(bin);
		double error = hNtMinUnf->GetBinError(bin);
		if( yield > 0. ) hRelStatUnc->SetBinContent(bin, error / yield);
	}

	lOut->Add(hNtMinRM);
	
	if (mHandler->GetFlagMC()) 	lOut->Add(hNtMinClosure);
	//else 	ComparisonPublished(hRtUnf);

	lOut->Add(hNtMinUnf);
	lOut->Add(hRtMinUnf);
	lOut->Add(hRelStatUnc);

	TGraph* gChi2 = nullptr;
	if (mHandler->GetFlagMC())	{
		mUnf->DrawNchClosure(hNtMinGen,hNtMinUnf,-0.5,15.5,"#it{N}_{T,min.}","Unfolded/True",dOut);
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

	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");
	fUnfOut->cd();
	lOut->Write();
	(mUnf->GetObjArray())->Write();
	fUnfOut->Close();
	delete fUnfOut;

	//gSystem->Exec(Form("cp ./%s/1D_newClass_mc.root ./%s/1D_newClass_data.root",dOut));	

	mDirFile->cd();
	lOut->Write();
}

void MyAnalysisMC::DoUnfoldingNtMax() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding_max";
	const bool eRM = false;
	const int NumberOfIters = 5;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtMaxRM = FlipMatrix(hNtMaxRM,20);
	// ROWWISE NORMALISATION IS PERFORMED IN mUnf
	//cout << "Histograms: " << hNtRec << " " << hNt << " " << hNtGen << " " << hNtRM << "\n";
	
	if (eRM) mUnf->ExtrapolateRM(hNtMaxRM);
	mUnf->SetnIter(NumberOfIters);
	mUnf->SetError(hNtMaxRec);
	mUnf->SetError(hNtMaxGen);
	mUnf->SaveSolutionNT(kTRUE);
	
	if (mHandler->GetFlagMC()) mUnf->Setup(hNtMaxRec,hNtMaxGen,hNtMaxRM);
	else mUnf->Setup(hNtMax,hNtMaxGen,hNtMaxRM);
	
	mUnf->Unfold();
	mUnf->V2H();
	
	printf(" - Unfolding region : %s_max\n",mUnf->GetRegion());

	hNtMaxUnf = (TH1F*)(mUnf->GetUnfoldedDistH())->Clone("hNtMaxUnf");
	if(mHandler->GetFlagMC())	hNtMaxClosure = (TH1F*)(mUnf->GetClosureH())->Clone("hNtMaxClosure");
	
	TH1F* hNTMax = (TH1F*)hNtMaxUnf->Clone("_hNTMax");
	hRtMaxUnf = (TH1F*)mUnf->RebinNT2RT(hNtMaxUnf, kTRUE);
	hRtMaxUnf->Scale(1.0/hRtMaxUnf->Integral());
	
	//! Relative statistical uncertainty
	//! NT distribution with final bins to plot
	TH1F* hRelStatUnc = (TH1F*)hNtUnf->Clone("hRelStatUnc_NtMax");
	hRelStatUnc->Reset();

	TH1F* hNch = (TH1F*)hNtMaxUnf->Clone("hNTMax");
	hNch->Reset();

	for(int bin = 1; bin <= hRelStatUnc->GetNbinsX(); bin++){
		double yield = hNtMaxUnf->GetBinContent(bin);
		double error = hNtMaxUnf->GetBinError(bin);
		if( yield > 0. ) hRelStatUnc->SetBinContent(bin, error / yield);
	}

	lOut->Add(hNtMaxRM);
	
	if (mHandler->GetFlagMC()) 	lOut->Add(hNtMaxClosure);
	//else 	ComparisonPublished(hRtUnf);

	lOut->Add(hNtMaxUnf);
	lOut->Add(hRtMaxUnf);
	lOut->Add(hRelStatUnc);

	TGraph* gChi2 = nullptr;
	if (mHandler->GetFlagMC())	{
		mUnf->DrawNchClosure(hNtMaxGen,hNtMaxUnf,-0.5,25.5,"#it{N}_{T,max.}","Unfolded/True",dOut);
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

	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");
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
	const int NumberOfIters = 5;
	const char* Regions[NREGIONS] = {"Trans1D","Toward","Away","TransMin1D","TransMax1D"};

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	Int_t typesUnf[] = {RC, D};
	size_t sizeoftypesUnf = sizeof(typesUnf)/sizeof(typesUnf[0]);

	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 0; iReg < 5; ++iReg)		{
	for (int itype = 0; itype < sizeoftypesUnf; ++itype)		{
	Int_t iType = typesUnf[itype];
	if (iType == D && iSp != Xi && iSp != phi) continue;
	if (iType == D && iSp == phi && iReg > 2) continue;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(7+iSp);
		obj->SetMCAnalysis(iType == RC ? 1 : 0);
		obj->SaveSolutionNT(kTRUE);

		printf(" - Unfolding species %s vs NT in region : %s\n",SPECIES[iSp],obj->GetRegion());
		cout << "RM has " << hNtRM->GetNbinsX() << " x " << hNtRM->GetNbinsY() << endl;
		cout << "hMC has " << hPIDPtNtCorr[iSp][iType][iReg]->GetNbinsY() << endl;

		// Calling with the last argument makes it rebin the unfolding matrix according to its Y-axis
		obj->UnfoldV02D(hNtRM,hPIDPtNtCorr[iSp][iType][iReg],hPIDPtNtCorr[iSp][iType][iReg], 
			//(iType ==D && iSp == phi) ? hPhiPtNtBinning : nullptr ); // normal
			( (iType == RC||iType == D) && iSp == phi) ? hPhiPtNtBinning : nullptr, // only for testing closures
			//( (iType == RC||iType == D) && iSp == phi) ? (TH1F*)hPIDPtNtCorr[iSp][MC][iReg]->ProjectionY(Form("py%i_%i_%i",iSp,iType,iReg),1,51) : nullptr); // we need a reweighting distribution for rebinning
			( (iType == RC||iType == D) && iSp == phi) ? hNtPID[phi][RC][iReg] : nullptr); // we need a reweighting distribution for rebinning

		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << hPIDPtNtCorr[iSp][iType][iReg] << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		// first is Gen, second is Unf

		// APPLY PROPER ERRORS!
		hPIDPtNtCorrUnf[iSp][iType][iReg] = (TH2F*)Arr->FindObject("_hPtvsNch")->Clone(Form("hPIDPtNtCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		for (int iX = 1; iX < hPIDPtNtCorrUnf[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hPIDPtNtCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float bine = hPIDPtNtCorr[iSp][iType][iReg]->GetBinError(iX,iY);
				if (bine>0) hPIDPtNtCorrUnf[iSp][iType][iReg]->SetBinError(iX,iY,bine);
			}
		}

		if (iType==RC)	obj->GetMCclosureinRTBins(hPIDPtNtCorr[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	//normal
		//if (iSp == phi && iType==RC)	obj->GetMCclosureinRTBinsPhi(hPIDPtNtCorr[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	//testing phi closure in coarser binning
		//else if (iType==RC)	obj->GetMCclosureinRTBins(hPIDPtNtCorr[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	//normal

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
		TDirectory* dir = fUnfOut1D->mkdir(Form("%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	}	}

	lOut->Write();
	fUnfOut1D->Close();
	delete fUnfOut1D;
	mDirFile->cd();
	lOut->Write();

	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisMC::DoUnfolding1DMin() {

	mHandler->root()->SetBatch(kTRUE);

	const char* dOut = "results_unfolding_min";
	const bool eRM = false;
	const int NumberOfIters = 5;
	const char* Regions[4] = {"Trans1D","Toward","Away","TransMin1D"};
	

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};
	Int_t typesUnf[] = {RC, D};
	size_t sizeoftypesUnf = sizeof(typesUnf)/sizeof(typesUnf[0]);
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 0; iReg < NREGIONS-1; ++iReg)		{
	for (int itype = 0; itype < sizeoftypesUnf; ++itype)		{
	Int_t iType = typesUnf[itype];
	if (iType == D && iSp != Xi ) continue;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(7+iSp);
		obj->SetMCAnalysis(iType == RC ? 1 : 0);
		obj->SaveSolutionNT(kTRUE);

		printf(" - Unfolding species %s in region : %s\n",SPECIES[iSp],obj->GetRegion());
		cout << "RM has " << hNtMinRM->GetNbinsX() << " x " << hNtMinRM->GetNbinsY() << endl;
		cout << "hMC has " << hPIDPtNtMinCorr[iSp][iType][iReg]->GetNbinsX() << endl;
		obj->UnfoldV02D(hNtMinRM,hPIDPtNtMinCorr[iSp][iType][iReg],hPIDPtNtMinCorr[iSp][iType][iReg],nullptr,nullptr);
		
		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << hPIDPtNtMinCorr[iSp][iType][iReg] << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		// first is Gen, second is Unf

		// APPLY PROPER ERRORS!
		hPIDPtNtMinCorrUnf[iSp][iType][iReg] = (TH2F*)Arr->FindObject("_hPtvsNch")->Clone(Form("hPIDPtNtMinCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		for (int iX = 1; iX < hPIDPtNtMinCorrUnf[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hPIDPtNtMinCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float bine = hPIDPtNtMinCorr[iSp][iType][iReg]->GetBinError(iX,iY);
				if (bine>0) hPIDPtNtMinCorrUnf[iSp][iType][iReg]->SetBinError(iX,iY,bine);
			}
		}

		if (iType==RC)	obj->GetMCclosureinRTMinBins(hPIDPtNtMinCorr[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	

		// NORMALISE
		/*for (int iNt = 1; iNt < hPIDPtNtMinCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {
			Double_t NormEv = 0;
			NormEv = hNtMinUnf->Integral(1,50) > 0 ? (double)hNtMinUnf->GetBinContent(iNt) / hNtMinUnf->Integral(1,50) : 1;
			//cout << "Normalisingggg by " << NormEv << endl;
			TH2F* htmp = hPIDPtNtMinCorrUnf[iSp][iType][iReg];

			TH1F* hpt = (TH1F*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
			delete hpt;
		}*/

		lOut->Add(hPIDPtNtMinCorrUnf[iSp][iType][iReg]);
		fUnfOut1D->cd();
		TDirectory* dir = fUnfOut1D->mkdir(Form("%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	}	}

	lOut->Write();
	fUnfOut1D->Close();
	delete fUnfOut1D;
	mDirFile->cd();
	lOut->Write();

	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisMC::DoUnfolding1DMax() {

	
	mHandler->root()->SetBatch(kTRUE);

	const char* dOut = "results_unfolding_max";
	const bool eRM = false;
	const int NumberOfIters = 5;
	const char* Regions[5] = {"Trans1D","Toward","Away","TransMin1D","TransMax1D"};
	

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};
	Int_t typesUnf[] = {RC, D};
	size_t sizeoftypesUnf = sizeof(typesUnf)/sizeof(typesUnf[0]);
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	if (iReg == 3) continue;
	for (int itype = 0; itype < sizeoftypesUnf; ++itype)		{
	Int_t iType = typesUnf[itype];
	if (iType == D && iSp != Xi ) continue;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(7+iSp);
		obj->SetMCAnalysis(iType == RC ? 1 : 0);
		obj->SaveSolutionNT(kTRUE);

		printf(" - Unfolding species %s in region : %s\n",SPECIES[iSp],obj->GetRegion());
		cout << "RM has " << hNtMaxRM->GetNbinsX() << " x " << hNtMaxRM->GetNbinsY() << endl;
		cout << "hMC has " << hPIDPtNtMaxCorr[iSp][iType][iReg]->GetNbinsX() << endl;
		obj->UnfoldV02D(hNtMaxRM,hPIDPtNtMaxCorr[iSp][iType][iReg],hPIDPtNtMaxCorr[iSp][iType][iReg],nullptr,nullptr);
		
		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << hPIDPtNtMaxCorr[iSp][iType][iReg] << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		// first is Gen, second is Unf

		// APPLY PROPER ERRORS!
		hPIDPtNtMaxCorrUnf[iSp][iType][iReg] = (TH2F*)Arr->FindObject("_hPtvsNch")->Clone(Form("hPIDPtNtMaxCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		for (int iX = 1; iX < hPIDPtNtMaxCorrUnf[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hPIDPtNtMaxCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float bine = hPIDPtNtMaxCorr[iSp][iType][iReg]->GetBinError(iX,iY);
				if (bine>0) hPIDPtNtMaxCorrUnf[iSp][iType][iReg]->SetBinError(iX,iY,bine);
			}
		}

		if (iType==RC)	obj->GetMCclosureinRTMaxBins(hPIDPtNtMaxCorr[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	

		// NORMALISE
		/*for (int iNt = 1; iNt < hPIDPtNtMaxCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {
			Double_t NormEv = 0;
			NormEv = hNtMaxUnf->Integral(1,50) > 0 ? (double)hNtMaxUnf->GetBinContent(iNt) / hNtMaxUnf->Integral(1,50) : 1;
			//cout << "Normalisingggg by " << NormEv << endl;
			TH2F* htmp = hPIDPtNtMaxCorrUnf[iSp][iType][iReg];

			TH1F* hpt = (TH1F*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
			delete hpt;
		}*/

		lOut->Add(hPIDPtNtMaxCorrUnf[iSp][iType][iReg]);
		fUnfOut1D->cd();
		TDirectory* dir = fUnfOut1D->mkdir(Form("%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	}	}

	lOut->Write();
	fUnfOut1D->Close();
	delete fUnfOut1D;
	mDirFile->cd();
	lOut->Write();

	mHandler->root()->SetBatch(kFALSE);
}


void MyAnalysisMC::BinHistogramsIntoRT() {

	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{
	Int_t iSp = Xi;
	Int_t iType = D;

	hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] = 0x0;
	hPIDPtRtMinCorrUnf[iSp][iType][iReg][iRtBin] = 0x0;
	hPIDPtRtMaxCorrUnf[iSp][iType][iReg][iRtBin] = 0x0;

		// REBIN
		int lowedge = mUnf->GetRTBin(iRtBin,kTRUE);
		int upedge  = mUnf->GetRTBin(iRtBin,kFALSE);
		Double_t NormScaleNt = hNtUnf->Integral(lowedge,upedge)/hNtUnf->Integral(1,51);

		cout << "iReg "  << iReg << " iRtBin " << iRtBin << " " << hPIDPtNtCorrUnf[iSp][iType][iReg] << " " << hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] << endl;
		if (hPIDPtNtCorrUnf[iSp][iType][iReg]) {
			hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]
				= (TH1F*)hPIDPtNtCorrUnf[iSp][iType][iReg]->ProjectionX(Form("hPIDPtRtCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),lowedge,upedge);
			hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]->Scale(1./NormScaleNt);
		}
		cout << "iReg "  << iReg << " iRtBin " << iRtBin << " " << hPIDPtNtCorrUnf[iSp][iType][iReg] << " " << hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] << endl;


		lowedge = mUnf->GetRTMinBin(iRtBin,kTRUE);
		upedge  = mUnf->GetRTMinBin(iRtBin,kFALSE);
		Double_t NormScaleNtMin = hNtMinUnf->Integral(lowedge,upedge)/hNtMinUnf->Integral(1,51);

		
		if (iReg != 0 && iReg != 4) if (hPIDPtNtMinCorrUnf[iSp][iType][iReg]) {
			hPIDPtRtMinCorrUnf[iSp][iType][iReg][iRtBin]
				= (TH1F*)hPIDPtNtMinCorrUnf[iSp][iType][iReg]->ProjectionX(Form("hPIDPtRtMinCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),lowedge,upedge);
			hPIDPtRtMinCorrUnf[iSp][iType][iReg][iRtBin]->Scale(1./NormScaleNtMin);
		}
		
		lowedge = mUnf->GetRTMaxBin(iRtBin,kTRUE);
		upedge  = mUnf->GetRTMaxBin(iRtBin,kFALSE);
		Double_t NormScaleNtMax = hNtMaxUnf->Integral(lowedge,upedge)/hNtMaxUnf->Integral(1,51);

		if (iReg != 0 && iReg != 3) if (hPIDPtNtMaxCorrUnf[iSp][iType][iReg]) {
			hPIDPtRtMaxCorrUnf[iSp][iType][iReg][iRtBin]
				= (TH1F*)hPIDPtNtMaxCorrUnf[iSp][iType][iReg]->ProjectionX(Form("hPIDPtRtMaxCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),lowedge,upedge);
			hPIDPtRtMaxCorrUnf[iSp][iType][iReg][iRtBin]->Scale(1./NormScaleNtMax);
		}

		cout << "Writing " << iReg << " " << hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] << " " << hPIDPtRtMinCorrUnf[iSp][iType][iReg][iRtBin] << " " << hPIDPtRtMaxCorrUnf[iSp][iType][iReg][iRtBin] << endl;
	if (hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]) hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]->Write();
	if (iReg != 0 && iReg != 4)	if (hPIDPtRtMinCorrUnf[iSp][iType][iReg][iRtBin]) hPIDPtRtMinCorrUnf[iSp][iType][iReg][iRtBin]->Write();
	if (iReg != 0 && iReg != 3) if (hPIDPtRtMaxCorrUnf[iSp][iType][iReg][iRtBin]) hPIDPtRtMaxCorrUnf[iSp][iType][iReg][iRtBin]->Write();

	}	}

	TH1F* hNtUnfRebin = RebinHistogram1D(hNtUnf,hPhiPtNtBinning);

	for (int iReg = 0; iReg < 3; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{
	Int_t iSp = phi;
	Int_t iType = D;

	hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] = 0x0;

		// REBIN
		int lowedge = GetRTBinPhi(iRtBin,kTRUE);
		int upedge  = GetRTBinPhi(iRtBin,kFALSE);
		Double_t NormScaleNt = hNtUnfRebin->Integral(lowedge,upedge)/hNtUnfRebin->Integral(1,7);

		cout << "iReg "  << iReg << " iRtBin " << iRtBin << " " << hPIDPtNtCorrUnf[iSp][iType][iReg] << " " << hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] << endl;
		if (hPIDPtNtCorrUnf[iSp][iType][iReg]) {
			hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]
				= (TH1F*)hPIDPtNtCorrUnf[iSp][iType][iReg]->ProjectionX(Form("hPIDPtRtCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),lowedge,upedge);
			hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]->Scale(1./NormScaleNt);
		}
		cout << "iReg "  << iReg << " iRtBin " << iRtBin << " " << hPIDPtNtCorrUnf[iSp][iType][iReg] << " " << hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] << endl;


	cout << "Writing " << iReg << " " << hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin] << " " << hPIDPtRtMinCorrUnf[iSp][iType][iReg][iRtBin] << " " << hPIDPtRtMaxCorrUnf[iSp][iType][iReg][iRtBin] << endl;
	if (hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]) hPIDPtRtCorrUnf[iSp][iType][iReg][iRtBin]->Write();
	
	}	}

}

Int_t MyAnalysisMC::GetRTBinPhi(const int& binRt, bool isLowEdge)	{

	int binNch = -1;

	//! <NT> = 7. 43
	
	if( binRt == 0 ){
		if(isLowEdge) binNch = 1;
		else binNch = 50;
	}
	else if( binRt == 1 ){ //! From NT = 0 to NT = 6
		if(isLowEdge) binNch = 1;
		else binNch = 1;//else binNch = 4;
	}
	else if( binRt == 2 ){ //! From NT = 6 to NT = 11
		if(isLowEdge) binNch = 2;
		else binNch = 3;
	}
	else if( binRt == 3 ){//! From NT = 12 to NT = 18
		if(isLowEdge) binNch = 4;
		else binNch = 5;
	}
	else if( binRt == 4 ){//! From NT = 19 to NT = 36
		if(isLowEdge) binNch = 6;
		else binNch = 6;
	}
	else{
		if(isLowEdge) binNch = 39;
		else binNch = 50;
	}

	return binNch;

}


/*void MyAnalysisMC::CalculateUnfoldingSysNt() {

	TFile* fClosuresNt = new TFile("../files/unfolding/results_unfolding/2D_newClass_mc.root","READ");
	TFile* fClosuresNtMin = new TFile("../files/unfolding/results_unfolding_min/2D_newClass_mc.root","READ");
	TFile* fClosuresNtMax = new TFile("../files/unfolding/results_unfolding_max/2D_newClass_mc.root","READ");
	const char* Regions[NREGIONS] = {"Trans1D","Toward","Away","TransMin1D","TransMax1D"};

	// RT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)	{

		TDirectoryFile* dir;
		if (iReg==3) dir = (TDirectoryFile*)fClosuresNtMin->Get(Form("%s_%s",SPECIES[iSp],REGIONS[iReg]));
		else if (iReg==4) dir = (TDirectoryFile*)fClosuresNtMax->Get(Form("%s_%s",SPECIES[iSp],REGIONS[iReg]));
		else dir = (TDirectoryFile*)fClosuresNt->Get(Form("%s_%s",SPECIES[iSp],REGIONS[iReg]));

		// DIVIDE BY RT INTEGRATED
		TH1D* hClosureD = (TH1D*)dir->Get(Form("hUnToGen_%s_%s_%i",SPECIES[iSp],Regions[iReg],0));
		for (int iRtBin = 1; iRtBin < NRTBINS; ++iRtBin)	{
		
			TH1D* hClosure = (TH1D*)dir->Get(Form("hUnToGen_%s_%s_%i",SPECIES[iSp],Regions[iReg],iRtBin));
			hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->Divide(hClosure,hClosureD,1.,1.,"");	

			// CALCULATE UNCERTAINTIES
			for( Int_t i=1; i<hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetNbinsX()+1; i++){ 
				Double_t err1 = hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinError(i);
				Double_t errHM = hClosureD->GetBinError(i);
				Double_t errR = TMath::Sqrt(TMath::Abs(err1*err1 - errHM*errHM));
				errR = (hClosureD->GetBinContent(i) > 0) ? errR / hClosureD->GetBinContent(i) : 0;
				hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->SetBinError(i,errR);
			}

			// CALCULATE MAX DEVIATIONS
			for (int iBin = 1; iBin<hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetNbinsX()+1;iBin++) {
				Double_t maxD = 0;
				Double_t varD = TMath::Abs(hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinContent(iBin)-1.);
				if (varD>maxD && varD>0.*hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinError(iBin)) maxD=varD;
						
				Double_t model = (iReg==1||iReg==2)? 0.005*0.005 : 0.02*0.02; // model dependence systematic						
				if (hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinContent(iBin)>0) hV0PtRtSysMaxDUnfolding[iSp][iReg][iRtBin]->SetBinContent(iBin,TMath::Sqrt(model+maxD*maxD));
				if (hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinContent(iBin)>0) hV0PtRtSysMaxDUnfolding[iSp][iReg][iRtBin]->SetBinError(iBin,hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinError(iBin));
				if (hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinContent(iBin)>0) hV0PtRtSysMaxDUnfoldingRatioToHM[iSp][iReg][iRtBin]->SetBinContent(iBin,TMath::Sqrt(maxD*maxD));
				if (hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinContent(iBin)>0) hV0PtRtSysMaxDUnfoldingRatioToHM[iSp][iReg][iRtBin]->SetBinError(iBin,hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetBinError(iBin));
			}

			// SMOOTHEN BY FITTING   (TERNARIES USED FOR TREATING UNSTABLE OUTLIERS)
			TF1* fpol = new TF1("fpol",(iSp==1&&iReg==4&&iRtBin==4)?"pol1":"pol2",0.4,8.0);
			hV0PtRtSysMaxDUnfolding[iSp][iReg][iRtBin]->Fit(fpol,(iSp==2&&(iReg==4||iReg==0)&&iRtBin==4)?"W0Q":"0Q");
			for (int iBin = 1; iBin<hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetNbinsX()+1;iBin++) {
				hV0PtRtSysMaxDUnfoldingSmooth[iSp][iReg][iRtBin]->SetBinContent(iBin,
					fpol->Eval(hV0PtRtSysMaxDUnfoldingSmooth[iSp][iReg][iRtBin]->GetBinCenter(iBin)));
			}
			hV0PtRtSysMaxDUnfoldingRatioToHM[iSp][iReg][iRtBin]->Fit(fpol,(iSp==2&&(iReg==4||iReg==0)&&iRtBin==4)?"W0Q":"0Q");
			for (int iBin = 1; iBin<hV0PtRtSysUnfolding[iSp][iReg][iRtBin]->GetNbinsX()+1;iBin++) {
				hV0PtRtSysMaxDUnfoldingSmoothRatioToHM[iSp][iReg][iRtBin]->SetBinContent(iBin,
					fpol->Eval(hV0PtRtSysMaxDUnfoldingSmoothRatioToHM[iSp][iReg][iRtBin]->GetBinCenter(iBin)));
			}

		}

	}	}

	// CORRELATION OF L TO K RATIOS
	for (int iReg = 0; iReg < NREGIONS; ++iReg)	{
	for (int iRtBin = 1; iRtBin < NRTBINS; ++iRtBin)	{
		
		hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->Divide(hV0PtRtSysUnfolding[2][iReg][iRtBin],hV0PtRtSysUnfolding[1][iReg][iRtBin],1.,1.,"");	

		// CALCULATE UNCERTAINTIES
		for( Int_t i=1; i<hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetNbinsX()+1; i++){ 
			Double_t err1 = hV0PtRtSysUnfolding[2][iReg][iRtBin]->GetBinError(i);
			Double_t errHM = hV0PtRtSysUnfolding[1][iReg][iRtBin]->GetBinError(i);
			Double_t errR = TMath::Sqrt(TMath::Abs(err1*err1 - errHM*errHM));
			errR = (hV0PtRtSysUnfolding[1][iReg][iRtBin]->GetBinContent(i) > 0) ? errR / hV0PtRtSysUnfolding[1][iReg][iRtBin]->GetBinContent(i) : 0;
			hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->SetBinError(i,errR);
		}

			// CALCULATE MAX DEVIATIONS
		for (int iBin = 1; iBin<hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetNbinsX()+1;iBin++) {
			Double_t maxD = 0;
			Double_t varD = TMath::Abs(hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetBinContent(iBin)-1.);
			if (varD>maxD && varD>0.*hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetBinError(iBin)) maxD=varD;
						
			Double_t model = 0;
			if (hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetBinContent(iBin)>0) hV0PtRtSysMaxDUnfoldingLtoK0s[iReg][iRtBin]->SetBinContent(iBin,TMath::Sqrt(model+maxD*maxD));
			if (hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetBinContent(iBin)>0) hV0PtRtSysMaxDUnfoldingLtoK0s[iReg][iRtBin]->SetBinError(iBin,hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetBinError(iBin));
			
		}

			// SMOOTHEN BY FITTING   (TERNARIES USED FOR TREATING UNSTABLE OUTLIERS)
			TF1* fpol = new TF1("fpol",((iReg==2||iReg==4)&&iRtBin==4)?"pol1":"pol2",0.4,8.0);
			hV0PtRtSysMaxDUnfoldingLtoK0s[iReg][iRtBin]->Fit(fpol,(iReg==0&&iRtBin==4)?"W0Q":"0Q");
			for (int iBin = 1; iBin<hV0PtRtSysUnfoldingLtoK0s[iReg][iRtBin]->GetNbinsX()+1;iBin++) {
				hV0PtRtSysMaxDUnfoldingSmoothLtoK0s[iReg][iRtBin]->SetBinContent(iBin,
					fpol->Eval(hV0PtRtSysMaxDUnfoldingSmoothLtoK0s[iReg][iRtBin]->GetBinCenter(iBin)));
			}

		}

	}


	fClosuresNt->Close();
	fClosuresNtMin->Close();
	fClosuresNtMax->Close();
	mDirFile->cd();

}*/