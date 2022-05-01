#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TH3D.h>
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


using namespace std;		// Don't cringe.
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

	// Initialising treebug checker
	bugR = 0;
	bugPt = 0;

	// Initialize spherocities
	for (Int_t iType = 0; iType < NTYPE; ++iType)	{
	for (Int_t iMu = 0; iMu < NMULTI; ++iMu)	{
		mTS[iType][iMu] = new TransverseSpherocity();
		mTS[iType][iMu]->SetMinMulti(10);
		mTSNorm[iType][iMu] = new TransverseSpherocity();
		mTSNorm[iType][iMu]->SetMinMulti(10);
	}	}

	// Initialize inv.mass sigmas for competing mass rejection
	mParSigK0s	= new TF1("funcSigK0s","[0]+[1]*x+[2]/x",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParSigK0s->SetParameters(cuts::K0S_PARSIG[0],cuts::K0S_PARSIG[1],cuts::K0S_PARSIG[2]);
	mParMuK0s	= new TF1("funcMuK0s","(x<=1.6)*([0]+[1]*x+[2]*x*x)+(x>1.6)*[3]",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParMuK0s->SetParameters(cuts::K0S_PARMU[0],cuts::K0S_PARMU[1],cuts::K0S_PARMU[2],cuts::K0S_PARMU[3]);
	
	mParSigL	= new TF1("funcSigL","[0]+[1]*x+[2]/x",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParSigL->SetParameters(cuts::L_PARSIG[0],cuts::L_PARSIG[1],cuts::L_PARSIG[2]);
	mParMuL	= new TF1("funcMuL","(x<=1.9)*([0]+[1]*x+[2]*x*x)+(x>1.9)*([3]+[4]*x)",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParMuL->SetParameters(cuts::L_PARMU[0],cuts::L_PARMU[1],cuts::L_PARMU[2],cuts::L_PARMU[3],cuts::L_PARMU[4]);
	
	mParSigK0s->Write();
	mParSigL->Write();
	mParMuK0s->Write();
	mParMuL->Write();
	//mParSigL	= cuts::L_PARSIG;
	//mParMuL		= cuts::L_PARMU;


	// SET CUT DEFAULT VALUES AND CUT DIRECTIONS FOR SYST. STUDIES
	for (Int_t iType = 0; iType < NTYPE; ++iType)	{
	for (Int_t iSp = 1; iSp < NSPECIES; ++iSp)	{
	
		mCutsDir[iSp][iType][DCAdd] = 1; 		mCutsVal[iSp][DCAdd] = cuts::V0_DCADD;
		mCutsDir[iSp][iType][CPA] = -1;			mCutsVal[iSp][CPA] = cuts::V0_CPA;
		mCutsDir[iSp][iType][RadiusL] = -1;		mCutsVal[iSp][RadiusL] = cuts::V0_R[0];
		mCutsDir[iSp][iType][RadiusH] = 1; 		mCutsVal[iSp][RadiusH] = cuts::V0_R[1];
		mCutsDir[iSp][iType][FastSignal] = -1; 	mCutsVal[iSp][FastSignal] = cuts::V0_FASTSIGNAL;
		mCutsDir[iSp][iType][cutsV0cuts] = 0; 	mCutsVal[iSp][cutsV0cuts] = 0;

		mCutsDir[iSp][iType][DCAPVpos] = -1; 	mCutsVal[iSp][DCAPVpos] = cuts::V0_D_DCAPVXY;
		mCutsDir[iSp][iType][NClusterpos] = -1; mCutsVal[iSp][NClusterpos] = cuts::V0_D_NCR;
		mCutsDir[iSp][iType][NClusterFpos] = -1; mCutsVal[iSp][NClusterFpos] = cuts::V0_D_NRATIO;
		mCutsDir[iSp][iType][DCAPVneg] = -1; 	mCutsVal[iSp][DCAPVneg] = cuts::V0_D_DCAPVXY;
		mCutsDir[iSp][iType][NClusterneg] = -1; mCutsVal[iSp][NClusterneg] = cuts::V0_D_NCR;;
		mCutsDir[iSp][iType][NClusterFneg] = -1; mCutsVal[iSp][NClusterFneg] = cuts::V0_D_NRATIO;
		mCutsDir[iSp][iType][cutsSizeof] = 0;	mCutsVal[iSp][cutsSizeof] = 0;

		switch (iSp) {
			case 1 : 
		mCutsDir[iSp][iType][CompMassK0s] = 0; mCutsVal[iSp][CompMassK0s] = 0;
		mCutsDir[iSp][iType][CompMassL] = 1; 	mCutsVal[iSp][CompMassL] = cuts::V0_COMP_NSIG;
		mCutsDir[iSp][iType][CompMassLbar] = 1; mCutsVal[iSp][CompMassLbar] = cuts::V0_COMP_NSIG;
		mCutsDir[iSp][iType][LifetimeK0s] = 0; 	mCutsVal[iSp][LifetimeK0s] = 0;
		mCutsDir[iSp][iType][LifetimeL] = 0; 	mCutsVal[iSp][LifetimeL] = 0;
		mCutsDir[iSp][iType][LifetimeLbar] = 0; mCutsVal[iSp][LifetimeLbar] = 0;
		mCutsDir[iSp][iType][NSigmaTPCposPiL] = (!mFlagMC) ? -1 : 0; mCutsVal[iSp][NSigmaTPCposPiL] = cuts::K0S_D_NSIGTPC[0];
		mCutsDir[iSp][iType][NSigmaTPCposPiH] = (!mFlagMC) ? 1 : 0; mCutsVal[iSp][NSigmaTPCposPiH] = cuts::K0S_D_NSIGTPC[1];
		mCutsDir[iSp][iType][NSigmaTPCposPrL] = 0; mCutsVal[iSp][NSigmaTPCposPrL] = 0;
		mCutsDir[iSp][iType][NSigmaTPCposPrH] = 0; mCutsVal[iSp][NSigmaTPCposPrH] = 0;
		mCutsDir[iSp][iType][NSigmaTPCnegPiL] = (!mFlagMC) ? -1 : 0; mCutsVal[iSp][NSigmaTPCnegPiL] = cuts::K0S_D_NSIGTPC[0];
		mCutsDir[iSp][iType][NSigmaTPCnegPiH] = (!mFlagMC) ? 1 : 0; mCutsVal[iSp][NSigmaTPCnegPiH] = cuts::K0S_D_NSIGTPC[1];
		mCutsDir[iSp][iType][NSigmaTPCnegPrL] = 0; mCutsVal[iSp][NSigmaTPCnegPrL] = 0;
		mCutsDir[iSp][iType][NSigmaTPCnegPrH] = 0; mCutsVal[iSp][NSigmaTPCnegPrH] = 0;
			break;

			case 2 :
		mCutsDir[iSp][iType][CPA] = -1;			mCutsVal[iSp][CPA] = cuts::L_CPA;			
		mCutsDir[iSp][iType][CompMassK0s] = 1; mCutsVal[iSp][CompMassK0s] = cuts::V0_COMP_NSIG;
		mCutsDir[iSp][iType][CompMassL] = 0; 	mCutsVal[iSp][CompMassL] = 0;
		mCutsDir[iSp][iType][CompMassLbar] = 0; mCutsVal[iSp][CompMassLbar] = 0;
		mCutsDir[iSp][iType][LifetimeK0s] = 0; 	mCutsVal[iSp][LifetimeK0s] = 0;
		mCutsDir[iSp][iType][LifetimeL] = 1; 	mCutsVal[iSp][LifetimeL] = cuts::L_TAU;
		mCutsDir[iSp][iType][LifetimeLbar] = 0; mCutsVal[iSp][LifetimeLbar] = 0;
		mCutsDir[iSp][iType][NSigmaTPCposPiL] = 0; mCutsVal[iSp][NSigmaTPCposPiL] = 0;
		mCutsDir[iSp][iType][NSigmaTPCposPiH] = 0; mCutsVal[iSp][NSigmaTPCposPiH] = 0;
		mCutsDir[iSp][iType][NSigmaTPCposPrL] = (!mFlagMC) ? -1 : 0; mCutsVal[iSp][NSigmaTPCposPrL] = cuts::K0S_D_NSIGTPC[0];
		mCutsDir[iSp][iType][NSigmaTPCposPrH] = (!mFlagMC) ? 1 : 0; mCutsVal[iSp][NSigmaTPCposPrH] = cuts::K0S_D_NSIGTPC[1];
		mCutsDir[iSp][iType][NSigmaTPCnegPiL] = (!mFlagMC) ? -1 : 0; mCutsVal[iSp][NSigmaTPCnegPiL] = cuts::K0S_D_NSIGTPC[0];
		mCutsDir[iSp][iType][NSigmaTPCnegPiH] = (!mFlagMC) ? 1 : 0; mCutsVal[iSp][NSigmaTPCnegPiH] = cuts::K0S_D_NSIGTPC[1];
		mCutsDir[iSp][iType][NSigmaTPCnegPrL] = 0; mCutsVal[iSp][NSigmaTPCnegPrL] = 0;
		mCutsDir[iSp][iType][NSigmaTPCnegPrH] = 0; mCutsVal[iSp][NSigmaTPCnegPrH] = 0;
			break;
		
			case 3 :
		mCutsDir[iSp][iType][CPA] = -1;			mCutsVal[iSp][CPA] = cuts::L_CPA;
		mCutsDir[iSp][iType][CompMassK0s] = 1; mCutsVal[iSp][CompMassK0s] = cuts::V0_COMP_NSIG;
		mCutsDir[iSp][iType][CompMassL] = 0; 	mCutsVal[iSp][CompMassL] = 0;	
		mCutsDir[iSp][iType][CompMassLbar] = 0; mCutsVal[iSp][CompMassLbar] = 0;
		mCutsDir[iSp][iType][LifetimeK0s] = 0; 	mCutsVal[iSp][LifetimeK0s] = 0;
		mCutsDir[iSp][iType][LifetimeL] = 0; 	mCutsVal[iSp][LifetimeL] = 0;
		mCutsDir[iSp][iType][LifetimeLbar] = 1; mCutsVal[iSp][LifetimeLbar] = cuts::L_TAU;
		mCutsDir[iSp][iType][NSigmaTPCposPiL] = (!mFlagMC) ? -1 : 0; mCutsVal[iSp][NSigmaTPCposPiL] = cuts::K0S_D_NSIGTPC[0];
		mCutsDir[iSp][iType][NSigmaTPCposPiH] = (!mFlagMC) ? 1 : 0; mCutsVal[iSp][NSigmaTPCposPiH] = cuts::K0S_D_NSIGTPC[1];
		mCutsDir[iSp][iType][NSigmaTPCposPrL] = 0; mCutsVal[iSp][NSigmaTPCposPrL] = 0;
		mCutsDir[iSp][iType][NSigmaTPCposPrH] = 0; mCutsVal[iSp][NSigmaTPCposPrH] = 0;
		mCutsDir[iSp][iType][NSigmaTPCnegPiL] = 0; mCutsVal[iSp][NSigmaTPCnegPiL] = 0;
		mCutsDir[iSp][iType][NSigmaTPCnegPiH] = 0; mCutsVal[iSp][NSigmaTPCnegPiH] = 0;
		mCutsDir[iSp][iType][NSigmaTPCnegPrL] = (!mFlagMC) ? -1 : 0; mCutsVal[iSp][NSigmaTPCnegPrL] = cuts::K0S_D_NSIGTPC[0];
		mCutsDir[iSp][iType][NSigmaTPCnegPrH] = (!mFlagMC) ? 1 : 0; mCutsVal[iSp][NSigmaTPCnegPrH] = cuts::K0S_D_NSIGTPC[1];
			break;

		}

	}	}

	return 0;
}

Int_t MyAnalysisV0::Make(Int_t iEv) {
	//printf("Looping in analysis %i \n", iEv);

	if (mFlagHist) return 0;

	enum { multMB, V0M, NCharged, V0M01, NCharged01, sizeofMult, RT };
	enum { sphMB, Jetty20, Iso20, Jetty10, Iso10, Jetty5, Iso5, Jetty1, Iso1};
	enum { D, RC, MC };

	// EVENT INFO HISTOGRAMS
	hEventMonitor->Fill(0);
	if (!mHandler->event()) return 1;
	MyEvent event(mHandler->event());
	hEventMonitor->Fill(1);


	
	// EVENT SELECTION
	hEventType->Fill(EVENTTYPES[0],1);	//preES MB
	hEventVz->Fill(event.GetZ());

	// SIGNAL LOSS STUDY
	if (mFlagMC) {
		Bool_t isINEL = false;
		for (int iP = 0; iP < mHandler->getNparticles(); ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);
			
			if (p.GetSign() != 0 && TMath::Abs(p.GetEta()) < 1.0) {
				isINEL = true;
				break;	}
		}

		if (isINEL && TMath::Abs(event.GetZ()) < 10.)	{
			for (int iP = 0; iP < mHandler->getNparticles(); ++iP)		{
				
				if (!mHandler->particle(iP)) continue;
				MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);
							
				if (!SelectParticle(p))		 			continue;
				for (int iSp = 1; iSp < NSPECIES; ++iSp)	{	
					if (p.GetPdgCode() != PDG_IDS[iSp]) continue;
					hV0PtNoTrigger[iSp]->Fill(p.GetPt());
				}	
			}
		}
	}

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


	// BUG HOTFIX FOR AURORATREES, unnecessary
	if (0) {
		MyV0 bugfix;
		if (mHandler->v0(0)) { bugfix = MyV0(mHandler->v0(0));
		bugfix.SetHandler(mHandler);
		if (TMath::Abs(bugfix.GetRadius()-bugR) < 0.0001
			&& TMath::Abs(bugfix.GetPt()-bugPt) < 0.0001) return 0;
		bugR = bugfix.GetRadius(); bugPt = bugfix.GetPt(); }
	}
	hEventMonitor->Fill(3);

	// EVENT CLASSIFICATION
	Int_t isEventCentral = 0;
	hEventV0MCentrality->Fill(event.GetV0MCentrality());
	hEventRefMult->Fill(event.GetRefMult());
	hEventV0MCentvRefMult->Fill(event.GetRefMult(),event.GetV0MCentrality());
	Bool_t isEventFHM 	= (NMULTI>1) && IsCentral(event,V0M);				// forward-rapidity high-multiplicity 0-10
	Bool_t isEventMHM 	= (NMULTI>2) && IsCentral(event,NCharged);		// mid-rapidity high-multiplicity 0-10
	Bool_t isEventFHM01 = (NMULTI>3) && IsCentral(event,V0M01);			// forward-rapidity high-multiplicity 0-1
	Bool_t isEventMHM01 = (NMULTI>4) && IsCentral(event,NCharged01);		// mid-rapidity high-multiplicity 0-1
	if (isEventFHM)		hEventType->Fill(EVENTTYPES[6],1);
	if (isEventMHM)		hEventType->Fill(EVENTTYPES[7],1);
	if (isEventFHM01)	hEventType->Fill(EVENTTYPES[8],1);
	if (isEventMHM01)	hEventType->Fill(EVENTTYPES[9],1);
	
	isEventCentral = (isEventFHM || isEventMHM);
	Bool_t isEventMult[sizeofMult] = {1, isEventFHM, isEventMHM, isEventFHM01, isEventMHM01};

	// INITIALIZE SPHEROCITIES
	Double_t eventTS[NTYPE][sizeofMult];
	Bool_t isEventSphero[NTYPE][sizeofMult][NSPHERO] = {};
	for (Int_t iType = 0; iType < NTYPE; ++iType)	{
	for (Int_t iMu = 0; iMu < NMULTI; ++iMu)	{
		mTS[iType][iMu]->Reset();			// pT-weighted sphero
		mTSNorm[iType][iMu]->Reset();		// pT=1 sphero
		eventTS[iType][iMu] = -1.;
		for (Int_t iSph = 0; iSph < NSPHERO; iSph++) isEventSphero[iType][iMu][iSph] = (iSph==0) ? 1 : 0;
	}	}


	// TRACK LOOP TO CONSTRUCT SPHEROCITY AND FIND LEADING
	hEventMonitor->Fill(4);
	Int_t nTracks = mHandler->getNtracks();
	Int_t nParticles = (mFlagMC) ? mHandler->getNparticles() : 0;

	ptLead = -99.; phiLead = -99.;
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
			mTS[D][multMB]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[D][multMB]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt());
			// fill event counting tree here too 
		}

		if (isEventFHM) {
			mTS[D][V0M]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[D][V0M]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
		if (isEventMHM) {
			mTS[D][NCharged]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[D][NCharged]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
		if (isEventFHM01) {
			mTS[D][V0M01]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[D][V0M01]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
		}
		if (isEventMHM01) {
			mTS[D][NCharged01]->AddTrack(t.GetPx(), t.GetPy());
			mTSNorm[D][NCharged01]->AddTrack(t.GetPx()/t.GetPt(), t.GetPy()/t.GetPt()); 
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

				mTS[MC][multMB]->AddTrack(p.GetPx(), p.GetPy());
				mTSNorm[MC][multMB]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventFHM) mTS[MC][V0M]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventFHM) mTSNorm[MC][V0M]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventMHM) mTS[MC][NCharged]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventMHM) mTSNorm[MC][NCharged]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventFHM01) mTS[MC][V0M01]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventFHM01) mTSNorm[MC][V0M01]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

				if (isEventMHM01) mTS[MC][NCharged01]->AddTrack(p.GetPx(), p.GetPy());
				if (isEventMHM01) mTSNorm[MC][NCharged01]->AddTrack(p.GetPx()/p.GetPt(), p.GetPy()/p.GetPt());

			}	
		}
	}

	if (isEventCentral) {

		// DETERMINING SPHEROCITY CLASS OF THE EVENT
		const Float_t* cutsJetty[]	= {0x0, cuts::EV_SPH_V0M_JETTY, cuts::EV_SPH_NCharged_JETTY,
									cuts::EV_SPH_V0M01_JETTY, cuts::EV_SPH_NCharged01_JETTY};
		const Float_t* cutsIso[]	= {0x0, cuts::EV_SPH_V0M_ISO, cuts::EV_SPH_NCharged_ISO,
									cuts::EV_SPH_V0M01_ISO, cuts::EV_SPH_NCharged01_ISO};

		for (Int_t iType = 0; iType < NTYPE; ++iType)	{		// if MC, also do MC and RC
		for (Int_t iMu = 1; iMu < NMULTI; ++iMu)	{
			
			// CHOOSE PT=1 SPHEROCITY AS DEFAULT
			eventTS[iType][iMu] = mTSNorm[iType][iMu]->GetTransverseSpherocityTracks();
			for (Int_t iSph = 1; iSph < NSPHERO; ++iSph)	{
				if (iSph%2)		isEventSphero[iType][iMu][iSph] = (eventTS[iType][iMu] < cutsJetty[iMu][(iSph-1)/2] && eventTS[iType][iMu] > 0.);
				else 			isEventSphero[iType][iMu][iSph] = (eventTS[iType][iMu] > cutsIso[iMu][(iSph-1)/2] && eventTS[iType][iMu] < 1.);
				isEventSphero[iType][iMu][iSph] = isEventSphero[iType][iMu][iSph] && (NSPHERO>iSph);
			}
		}	}

		// FILLING HISTOGRAMS
		if (NMULTI>4) {
			if (isEventFHM)		hEventSpherocityV0M->Fill(eventTS[D][V0M]);
			if (isEventMHM)		hEventSpherocityNCharged->Fill(eventTS[D][NCharged]);
			if (isEventFHM01)	hEventSpherocityV0M01->Fill(eventTS[D][V0M01]);
			if (isEventMHM01)	hEventSpherocityNCharged01->Fill(eventTS[D][NCharged01]);
			hEventTSMCvRC->Fill(mTS[MC][V0M]->GetTransverseSpherocityTracks(),mTS[D][V0M]->GetTransverseSpherocityTracks());
			hEventTSNormMCvRC->Fill(eventTS[MC][V0M],eventTS[D][V0M]);
			hEventTSMCvNorm->Fill(mTS[MC][V0M]->GetTransverseSpherocityTracks(),eventTS[MC][V0M]);
			hEventTSRCvNorm->Fill(mTS[D][V0M]->GetTransverseSpherocityTracks(),eventTS[D][V0M]);
		}
	}

	// COUNT EVENT CLASS IN EVENT COUNTERS FOR NORMALIZATIONS
	for (Int_t iMu = 0; iMu < sizeofMult; ++iMu)	{
	for (Int_t iSph = 0; iSph < NSPHERO; ++iSph)	{		

		if (isEventMult[iMu] && isEventSphero[D][iMu][iSph])	hEventMultvSpheroD->Fill(SPHERO[iSph], MULTI[iMu], 1);
		if (isEventMult[iMu] && isEventSphero[MC][iMu][iSph])	hEventMultvSpheroMC->Fill(SPHERO[iSph], MULTI[iMu], 1);

	}	}

	// TRACK DEPENDENCE ON SPHEROCITY STUDY LOOP
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));	t.SetHandler(mHandler);

		if (!SelectTrack(t)) continue;

		ProcessTrack(t,D,multMB,sphMB);
		if (mFlagMC) ProcessTrack(t,RC,multMB,sphMB);

		if (isEventFHM) {
			ProcessTrack(t,D,V0M,sphMB);
			if (mFlagMC)	ProcessTrack(t,RC,V0M,sphMB);
			for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][V0M][iSph]) ProcessTrack(t,D,V0M,iSph);
			if (mFlagMC) for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M][iSph]) ProcessTrack(t,RC,V0M,iSph);
			/*if (isEventJetty[V0M])		ProcessTrack(t,D,V0M,Jetty);
			if (isEventIso[V0M])		ProcessTrack(t,D,V0M,Iso);			//
			if (isEventJettyMC[V0M])	ProcessTrack(t,RC,V0M,Jetty);		// study tracks also for true sphero, stored under RC
			if (isEventIsoMC[V0M])		ProcessTrack(t,RC,V0M,Iso);*/		
		} 

		if (isEventMHM) {
			ProcessTrack(t,D,NCharged,sphMB);
			if (mFlagMC) ProcessTrack(t,RC,NCharged,sphMB);
			for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][NCharged][iSph]) ProcessTrack(t,D,NCharged,iSph);
			if (mFlagMC) for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged][iSph]) ProcessTrack(t,RC,NCharged,iSph);
		}

		if (isEventFHM01) {
			ProcessTrack(t,D,V0M01,sphMB);
			if (mFlagMC) ProcessTrack(t,RC,V0M01,sphMB);
			for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][V0M01][iSph]) ProcessTrack(t,D,V0M01,iSph);
			if (mFlagMC) for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M01][iSph]) ProcessTrack(t,RC,V0M01,iSph);		
		} 

		if (isEventMHM01) {
			ProcessTrack(t,D,NCharged01,sphMB);
			if (mFlagMC) ProcessTrack(t,RC,NCharged01,sphMB);
			for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][NCharged01][iSph]) ProcessTrack(t,D,NCharged01,iSph);
			if (mFlagMC) for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged01][iSph]) ProcessTrack(t,RC,NCharged01,iSph);
		}

	}


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
		if (!t.IskITSrefit()) continue;		// cuts for tracks entering Nch^trans calculations
		if (!t.IsTPCOnlyRefit()) continue;
		if (t.GetPt()<0.15) continue;
		
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
		eventRt = (double)nChTrans/RT_DEN;
		hRt->Fill(eventRt);		}

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
	if (NMULTI>2) {
		if (isEventRT && isEventFHM) hNchTransvSpherocityV0M->Fill(eventTS[D][V0M],nChTrans);
		if (isEventRT && isEventMHM) hNchTransvSpherocityNCharged->Fill(eventTS[D][NCharged],nChTrans);
	}
	// STUDY DPHI DISTRIBUTION IN RT EVENTS
	if (isEventRT) for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));
		t.SetHandler(mHandler);
		if (t.GetEta() < cuts::V0_ETA[0] || t.GetEta() > cuts::V0_ETA[1] ) continue;
		
		if (!t.IskITSrefit()) continue;
		if (!t.IsTPCOnlyRefit()) continue;
		if (t.GetPt()<0.15) continue;
		
		if (TMath::Abs(t.GetPt()-ptLead)>1E-5) {
			hTrackDPhivNchTrans->Fill(nChTrans,mHandler->DeltaPhi(phiLead,t.GetPhi())); //TMath::Abs(phiLead-t.GetPhi()));
			hTrackDPhivNchTransMin->Fill(nChTransMin,mHandler->DeltaPhi(phiLead,t.GetPhi())); //TMath::Abs(phiLead-t.GetPhi()));
		}

	}




	//////////////////////////////////////////////////////////////////
	// MC RT CALCULATION
	//////////////////////////////////////////////////////////////////
	Double_t ptLeadMC = -99., phiLeadMC = -99.;
	Double_t eventRtMC = 0;
	Int_t nChTransMC = 0;
	Int_t nChTransAMC = 0;
	Int_t nChTransBMC = 0;
	Bool_t isSideAMinMC = 0;			// We need to keep track for also classifying the regions for spectra.
	Int_t nChTransMinMC = nChTransBMC;	// NT in the min transverse region
	Int_t nChTransMaxMC = nChTransAMC;	// NT in the max transverse region

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
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP)); p.SetHandler(mHandler); p.SetLabel(iP);

			if (p.GetEta() > cuts::V0_ETA[0] && p.GetEta() < cuts::V0_ETA[1]
				&& p.GetSign() != 0 && p.GetPt() > 0.15 && p.IsPrimary())	{

				if (!IsTrans(p.GetPhi(),phiLeadMC)) continue;
				nChTransMC++;

				if (WhatRegionSide(p.GetPhi(),phiLeadMC)==1) nChTransAMC++;
				if (WhatRegionSide(p.GetPhi(),phiLeadMC)==2) nChTransBMC++;


				if (ptLeadMC<5.0) continue;
				if (ptLeadMC>40.0) continue;

			}
		}

		if (nChTransAMC < nChTransBMC) {
			isSideAMinMC = 1;
			nChTransMinMC = nChTransAMC;
			nChTransMaxMC = nChTransBMC;	}

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
					for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M][iSph]) hTrackPt[MC][V0M][iSph]->Fill(p.GetPt());
				}

				if (isEventMHM) {
					hTrackPt[MC][NCharged][sphMB]->Fill(p.GetPt());
					for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged][iSph]) hTrackPt[MC][NCharged][iSph]->Fill(p.GetPt());
				}

				if (isEventFHM01) {
					hTrackPt[MC][V0M01][sphMB]->Fill(p.GetPt());
					for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M01][iSph]) hTrackPt[MC][V0M01][iSph]->Fill(p.GetPt());
				}

				if (isEventMHM01) {
					hTrackPt[MC][NCharged01][sphMB]->Fill(p.GetPt());
					for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged01][iSph]) hTrackPt[MC][NCharged01][iSph]->Fill(p.GetPt());
				}

			}

			//if (p.GetPdgCode() == 310) cout << "eta " << p.GetEta() << " y " << p.GetY() << endl;
				
			if (!SelectParticle(p)) continue;		// also contains Xi's !
			
			PartLabels.push_back(p.GetLabel());		// filling particles in a vector to avoid unnecessary nesting loops
			PartIds.push_back(iP);					// only particles in good |eta| with K0s,L,Lbar,Xi pdg's

			if (!p.IsPrimary()) continue;

			if (p.GetPdgCode() == 3312) hV0FeeddownMotherPt[2]->Fill(p.GetPt());		// filling primary Xi
			if (p.GetPdgCode() == -3312) hV0FeeddownMotherPt[3]->Fill(p.GetPt());
			if (p.GetPdgCode() == 3322) hV0FeeddownMotherPtXi0[2]->Fill(p.GetPt());		// filling primary Xi
			if (p.GetPdgCode() == -3322) hV0FeeddownMotherPtXi0[3]->Fill(p.GetPt());

			for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
				if (p.GetPdgCode() == PDG_IDS[iSp])		{

					// FILLING MC INFORMATION FOR V0s particles
					hV0Pt[iSp][MC][multMB][sphMB]->Fill(p.GetPt());
					hV0Eta[iSp][MC][multMB][sphMB]->Fill(p.GetEta());
					hV0Y[iSp][MC][multMB][sphMB]->Fill(p.GetY());

					tV0PtMCMB[iSp]->Fill(p.GetPt());

					if (isEventRT && iSp==1) hK0sDPhivNchTransMC->Fill(nChTransMC,mHandler->DeltaPhi(phiLead,p.GetPhi()));
					if (isEventRT && iSp==2) hLDPhivNchTransMC->Fill(nChTransMC,mHandler->DeltaPhi(phiLead,p.GetPhi()));
					if (isEventRT && iSp==3) hLbarDPhivNchTransMC->Fill(nChTransMC,mHandler->DeltaPhi(phiLead,p.GetPhi()));


					if (isEventFHM) {
						hV0Pt[iSp][MC][V0M][sphMB]->Fill(p.GetPt());
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M][iSph]) hV0Pt[iSp][MC][V0M][iSph]->Fill(p.GetPt());
					}

					if (isEventMHM) {
						hV0Pt[iSp][MC][NCharged][sphMB]->Fill(p.GetPt());
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged][iSph]) hV0Pt[iSp][MC][NCharged][iSph]->Fill(p.GetPt());
					}

					if (isEventFHM01) {
						hV0Pt[iSp][MC][V0M01][sphMB]->Fill(p.GetPt());
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M01][iSph]) hV0Pt[iSp][MC][V0M01][iSph]->Fill(p.GetPt());
					}

					if (isEventMHM01) {
						hV0Pt[iSp][MC][NCharged01][sphMB]->Fill(p.GetPt());
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged01][iSph]) hV0Pt[iSp][MC][NCharged01][iSph]->Fill(p.GetPt());
					}

					if (isEventRT)	{
						Int_t region = WhatRegion(p.GetPhi(),phiLead);
						hV0PtNt[iSp][MC][region]->Fill(p.GetPt(),nChTransMC);
						hV0PtNtMin[iSp][MC][region]->Fill(p.GetPt(),nChTransMinMC);
						
						Int_t regionSide = WhatRegionSide(p.GetPhi(),phiLead);
						if (!region) {
							hV0PtNt[iSp][MC][IsMinOrMax(isSideAMin,regionSide)]->Fill(p.GetPt(),nChTransMC);
							hV0PtNtMin[iSp][MC][IsMinOrMax(isSideAMin,regionSide)]->Fill(p.GetPt(),nChTransMinMC);
						}
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

					//cout << "bbb " << v0mc.GetPdgCode() << " " << v0.GetMCPdgCode() << endl;
					if (v0mc.GetPdgCode() == v0.GetMCPdgCode()) MCfound = true;
					// found a MC particle with the same label and PDGCode
					// MyV0::GetMCPdgCode() also requires both daughters to point to the same mother
				}
			}


			for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
				//cout << "Testing 1 iSp " << iSp << " with pdg " << v0.GetMCPdgCode() << " label " << v0.GetMCLabel() << endl;
				
				// SYSTEMATIC CUT VARIATIONS FOR EFF. CALCULATIONS
				if (mFlagMC) {
					Bool_t properDaughters = ((v0.GetPosTrackPdg() == PDG_IDS_DPOS[iSp]) && (v0.GetNegTrackPdg() == PDG_IDS_DNEG[iSp]));
					if (!properDaughters) continue;
					if (!v0.IsMCPrimary()) continue;	// considering only primaries for efficiency
					if (!MCfound) continue;
					ProcessV0SystVar(v0,iSp,RC,multMB,sphMB);	

				}
			}

			for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

				if (IsV0(v0,iSp,RC)) {
					
					//if (iSp==1) cout << "detected k0s with label " << v0.GetMCLabel() << endl;
					//else 
					//hV0FeeddownPDG[iSp]->Fill(v0mc.GetMotherPdgCode()); //
					//printf("code is %i and %i \n", v0mc.GetPdgCode(), v0mc.GetMotherPdgCode());
					
					Bool_t properDaughters = ((v0.GetPosTrackPdg() == PDG_IDS_DPOS[iSp]) && (v0.GetNegTrackPdg() == PDG_IDS_DNEG[iSp]));
					//if (iSp == 1 && !MCfound) cout << "mc code " << v0.GetMCPdgCode() << " pt " << v0.GetPt() << " primary " << v0.IsMCPrimary() << " daughters " << properDaughters << " granma " << v0.GetMCPrimaryPdgCode() << endl;
					if (!properDaughters) continue;
						
					/**MyParticle v0mc2; v0mc2.SetHandler(mHandler);
					for (unsigned int iP = 0; iP < PartLabels.size(); ++iP)	{
						v0mc2 = MyParticle(mHandler->particle(PartIds[iP])); v0mc2.SetHandler(mHandler); v0mc2.SetLabel(PartIds[iP]);
								

					// ask for primary
					// get the MC label of mother
					// get the Xi
					// fill histo  l.pt and xi.pt
					if (!v0.IsMCPrimary()) {
						printf("primary %i and grandma pdg code %i and v0 codes mc %i and rc %i \n", v0.IsMCPrimary(), v0.GetMCPrimaryPdgCode(),v0mc2.GetPdgCode(),v0.GetMCPdgCode());
						printf("labels are mc %i and rc %i ",v0mc2.GetLabel(),v0.GetMCPrimaryLabel());
						cout << iSp << endl;
						cout << properDaughters << endl;
						cout << "mc f " << MCfound << endl;
					}
					}
					*/

					if (!v0.IsMCPrimary()) {
							// using secondary particles to build the feeddown matrix
							// not requiring matched Particle (we store only primaries)

						if (iSp>1) {
							MyParticle xiMC; xiMC.SetHandler(mHandler);
							Double_t v0mass[] 	= {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
						
							for (unsigned int iP = 0; iP < PartLabels.size(); ++iP)	{
								xiMC = MyParticle(mHandler->particle(PartIds[iP])); xiMC.SetHandler(mHandler); xiMC.SetLabel(PartIds[iP]);
								
								if (xiMC.GetPdgCode() == 3312 && iSp==2 && xiMC.IsPrimary() && 
									v0.GetMCPrimaryPdgCode() == 3312 && (v0.GetMCPrimaryLabel() == xiMC.GetLabel()))	{
									hV0FeeddownMatrix[iSp]->Fill(xiMC.GetPt(),v0.GetPt(),v0mass[iSp]);
									hV0FeeddownMatrixXi0[iSp]->Fill(xiMC.GetPt(),v0.GetPt(),v0mass[iSp]);
								}
								if (xiMC.GetPdgCode() == -3312 && iSp==3 && xiMC.IsPrimary() && 
									v0.GetMCPrimaryPdgCode() == -3312 && (v0.GetMCPrimaryLabel() == xiMC.GetLabel()))	{
									hV0FeeddownMatrix[iSp]->Fill(xiMC.GetPt(),v0.GetPt(),v0mass[iSp]);
									hV0FeeddownMatrixXi0[iSp]->Fill(xiMC.GetPt(),v0.GetPt(),v0mass[iSp]);
								}

								if (xiMC.GetPdgCode() == 3322 && iSp==2 && xiMC.IsPrimary() && 
									v0.GetMCPrimaryPdgCode() == 3322 && (v0.GetMCPrimaryLabel() == xiMC.GetLabel()))	{
									hV0FeeddownMatrixXi0[iSp]->Fill(xiMC.GetPt(),v0.GetPt(),v0mass[iSp]);
								}
								if (xiMC.GetPdgCode() == -3322 && iSp==3 && xiMC.IsPrimary() && 
									v0.GetMCPrimaryPdgCode() == -3322 && (v0.GetMCPrimaryLabel() == xiMC.GetLabel()))	{
									hV0FeeddownMatrixXi0[iSp]->Fill(xiMC.GetPt(),v0.GetPt(),v0mass[iSp]);
								}
							}
						}
					}
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

					/*if (!v0.IsMCPrimary() && MCfound) {
							// using secondary particles to build the feeddown matrix


							cout << "v0 not primary " << iSp << endl;
							//if (iSp>1 && TMath::Abs(v0mc.GetMotherPdgCode()) == 3312) {
						if (iSp>1) {
							MyParticle xiMC; xiMC.SetHandler(mHandler);
						
							for (unsigned int iP = 0; iP < PartLabels.size(); ++iP)	{
								xiMC = MyParticle(mHandler->particle(PartIds[iP])); xiMC.SetHandler(mHandler); xiMC.SetLabel(PartIds[iP]);
								cout << "found particle w id " << xiMC.GetPdgCode() << endl;
								if (xiMC.GetPdgCode() == 3312 && iSp==2 && xiMC.IsPrimary() && v0mc.GetMotherPdgCode() == 3312)	{
									hV0FeeddownMatrix[iSp]->Fill(xiMC.GetPt(),v0.GetPt());
								}
								if (xiMC.GetPdgCode() == -3312 && iSp==3 && xiMC.IsPrimary() && v0mc.GetMotherPdgCode() == -3312)	{
									hV0FeeddownMatrix[iSp]->Fill(xiMC.GetPt(),v0.GetPt());
								}
							}
						}
					}*/


					if (!v0.IsMCPrimary()) continue;	// considering only primaries for efficiency
					if (!MCfound) continue;

					// MC v RC histograms
					hV0PtRCvMC[iSp]->Fill(v0mc.GetPt(), v0.GetPt());
					hV0EtaRCvMC[iSp]->Fill(v0mc.GetEta(), v0.GetEta());
					hV0PhiRCvMC[iSp]->Fill(v0mc.GetPhi(), v0.GetPhi());
						

					ProcessV0toHist(v0,iSp,RC,multMB,sphMB);		
					ProcessV0toTree(v0,iSp,RC,0);		

					if (isEventFHM) {
						ProcessV0toHist(v0,iSp,RC,V0M,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M][iSph]) ProcessV0toHist(v0,iSp,RC,V0M,iSph);
						// study RC spectra for true spher. ev classification
					}
						
					if (isEventMHM) {
						ProcessV0toHist(v0,iSp,RC,NCharged,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged][iSph]) ProcessV0toHist(v0,iSp,RC,NCharged,iSph);
					}

					if (isEventFHM01) {
						ProcessV0toHist(v0,iSp,RC,V0M01,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][V0M01][iSph]) ProcessV0toHist(v0,iSp,RC,V0M01,iSph);
					}
						
					if (isEventMHM01) {
						ProcessV0toHist(v0,iSp,RC,NCharged01,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[MC][NCharged01][iSph]) ProcessV0toHist(v0,iSp,RC,NCharged01,iSph);
					}

					if (isEventRT) {
						Int_t region = WhatRegion(v0.GetPhi(),phiLead);
						ProcessV0toHistRT(v0,iSp,RC,region,nChTrans,nChTransMin);
						Int_t regionSide = WhatRegionSide(v0.GetPhi(),phiLead);
						if (!region) ProcessV0toHistRT(v0,iSp,RC,IsMinOrMax(isSideAMin,regionSide),nChTrans,nChTransMin);
					}
				}
			}
		}	

		// in MC this part is only analyzed if RC selection was successful ?
		hV0Radius->Fill(v0.GetRadius());
		hV0ProperT->Fill(v0.GetIML()*v0.GetRadius()/v0.GetPt());
			
			// STUDY OF V0s IN DATA
			for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

				//printf("isV0 %i isV0vary %i \n", IsV0(v0,iSp,D), IsV0VaryCut(v0,iSp,D,RadiusL,cuts::V0_R[0]));
				if (IsV0(v0,iSp,D) != IsV0VaryCut(v0,iSp,D,NClusterpos,70.)) {
					printf("Houston we got a problem %i lt %f \n", iSp, v0.GetRadius());
					printf("%i and %i \n", IsV0(v0,iSp,D), IsV0VaryCut(v0,iSp,D,NClusterpos,70.));
				}
				Double_t v0mass[] 	= {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
				if (1) {
					for (Int_t iBin = 1; iBin < hV0IMvRadiusL[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,RadiusL,hV0IMvRadiusL[iSp]->GetXaxis()->GetBinCenter(iBin))) hV0IMvRadiusL[iSp]->Fill(hV0IMvRadiusL[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
					for (Int_t iBin = 1; iBin < hV0IMvDCAdd[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,DCAdd,hV0IMvDCAdd[iSp]->GetXaxis()->GetBinCenter(iBin))) hV0IMvDCAdd[iSp]->Fill(hV0IMvDCAdd[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
					for (Int_t iBin = 1; iBin < hV0IMvCPA[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,CPA,hV0IMvCPA[iSp]->GetXaxis()->GetBinCenter(iBin))) hV0IMvCPA[iSp]->Fill(hV0IMvCPA[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
					for (Int_t iBin = 1; iBin < hV0IMvFastSignal[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,FastSignal,hV0IMvFastSignal[iSp]->GetXaxis()->GetBinLowEdge(iBin))) hV0IMvFastSignal[iSp]->Fill(hV0IMvFastSignal[iSp]->GetXaxis()->GetBinLowEdge(iBin),v0mass[iSp]);
					
					for (Int_t iBin = 1; iBin < hV0IMvCompMass[iSp]->GetNbinsX()+1; iBin++)
						if ( (iSp==1 && IsV0VaryCut(v0,iSp,D,CompMassL,hV0IMvCompMass[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,CompMassLbar,hV0IMvCompMass[iSp]->GetXaxis()->GetBinCenter(iBin)) )
						|| (iSp>1 && IsV0VaryCut(v0,iSp,D,CompMassK0s,hV0IMvCompMass[iSp]->GetXaxis()->GetBinCenter(iBin)) ) )
						 hV0IMvCompMass[iSp]->Fill(hV0IMvCompMass[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
					
					if (iSp>1) for (Int_t iBin = 1; iBin < hV0IMvLifetime[iSp]->GetNbinsX()+1; iBin++)
						if ( (iSp==2 && IsV0VaryCut(v0,iSp,D,LifetimeL,hV0IMvLifetime[iSp]->GetXaxis()->GetBinCenter(iBin)))
						|| (iSp==3 && IsV0VaryCut(v0,iSp,D,LifetimeLbar,hV0IMvLifetime[iSp]->GetXaxis()->GetBinCenter(iBin)) ) ) 
						 hV0IMvLifetime[iSp]->Fill(hV0IMvLifetime[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
					
					for (Int_t iBin = 1; iBin < hV0IMvNSigmaTPC[iSp]->GetNbinsX()+1; iBin++)
						if ( (iSp==1 && IsV0VaryCut(v0,iSp,D,NSigmaTPCposPiL,-1.*hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCposPiH,hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCnegPiL,-1.*hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCnegPiH,hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin)) )
						|| (iSp==2 && IsV0VaryCut(v0,iSp,D,NSigmaTPCposPrL,-1.*hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCposPrH,hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCnegPiL,-1.*hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCnegPiH,hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin)) )
						|| (iSp==3 && IsV0VaryCut(v0,iSp,D,NSigmaTPCposPiL,-1.*hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCposPiH,hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCnegPrL,-1.*hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin))
									&& IsV0VaryCut(v0,iSp,D,NSigmaTPCnegPrH,hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin)) )	)
						 hV0IMvNSigmaTPC[iSp]->Fill(hV0IMvNSigmaTPC[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);

					for (Int_t iBin = 1; iBin < hV0IMvDCAPVpos[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,DCAPVpos,hV0IMvDCAPVpos[iSp]->GetXaxis()->GetBinCenter(iBin))) hV0IMvDCAPVpos[iSp]->Fill(hV0IMvDCAPVpos[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
					for (Int_t iBin = 1; iBin < hV0IMvDCAPVneg[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,DCAPVneg,hV0IMvDCAPVneg[iSp]->GetXaxis()->GetBinCenter(iBin))) hV0IMvDCAPVneg[iSp]->Fill(hV0IMvDCAPVneg[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
					
					for (Int_t iBin = 1; iBin < hV0IMvNCluster[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,NClusterpos,hV0IMvNCluster[iSp]->GetXaxis()->GetBinLowEdge(iBin))
						&& IsV0VaryCut(v0,iSp,D,NClusterneg,hV0IMvNCluster[iSp]->GetXaxis()->GetBinLowEdge(iBin)) ) 
							hV0IMvNCluster[iSp]->Fill(hV0IMvNCluster[iSp]->GetXaxis()->GetBinLowEdge(iBin),v0mass[iSp]);
					
					for (Int_t iBin = 1; iBin < hV0IMvNClusterF[iSp]->GetNbinsX()+1; iBin++)
						if (IsV0VaryCut(v0,iSp,D,NClusterFpos,hV0IMvNClusterF[iSp]->GetXaxis()->GetBinCenter(iBin))
						&& IsV0VaryCut(v0,iSp,D,NClusterFneg,hV0IMvNClusterF[iSp]->GetXaxis()->GetBinCenter(iBin)) ) 
							hV0IMvNClusterF[iSp]->Fill(hV0IMvNClusterF[iSp]->GetXaxis()->GetBinCenter(iBin),v0mass[iSp]);
				}

				//if (IsV0(v0,iSp,D) != IsV0VaryCut(v0,iSp,D,DCAPVpos,sysVar[iSp-1][sysDCAPVpos][deflt]))
				//printf("iSp %i isV0 %i isV0vary %i w cut %f and dca %f VS cut %f and dca %f \n", iSp, IsV0(v0,iSp,D), IsV0VaryCut(v0,iSp,D,DCAPVpos,sysVar[iSp-1][sysDCAPVpos][deflt]),
				//		sysVar[iSp-1][sysDCAPVpos][deflt], v0.GetCut(DCAPVpos), cuts::V0_D_DCAPVXY, TMath::Abs(v0.GetPosTrack()->GetImpactParameter(0)));

				// SYSTEMATIC CUT VARIATIONS
				if (!mFlagMC) {

					if (isEventFHM) {
						ProcessV0SystVar(v0,iSp,D,V0M,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][V0M][iSph]) ProcessV0SystVar(v0,iSp,D,V0M,iSph);
					}

					if (isEventMHM) {
						ProcessV0SystVar(v0,iSp,D,NCharged,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][NCharged][iSph]) ProcessV0SystVar(v0,iSp,D,NCharged,iSph);
					}


					if (isEventFHM01) {
						ProcessV0SystVar(v0,iSp,D,V0M01,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][V0M01][iSph]) ProcessV0SystVar(v0,iSp,D,V0M01,iSph);
					}

					if (isEventMHM01) {
						ProcessV0SystVar(v0,iSp,D,NCharged01,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][NCharged01][iSph]) ProcessV0SystVar(v0,iSp,D,NCharged01,iSph);
					}
				} 

				if (IsV0(v0,iSp,D)) {

					//if (iSp == 1) if (v0.GetMCPdgCode() == 310) continue;
					//if (iSp > 1) if (!v0.IsMCPrimary()) continue;
					ProcessV0toHist(v0,iSp,D,multMB,sphMB);
					if (isEventRT && iSp==1) hK0sDPhivNchTrans->Fill(nChTrans,mHandler->DeltaPhi(phiLead,v0.GetPhi()));
					if (isEventRT && iSp==2) hLDPhivNchTrans->Fill(nChTrans,mHandler->DeltaPhi(phiLead,v0.GetPhi()));
					if (isEventRT && iSp==3) hLbarDPhivNchTrans->Fill(nChTrans,mHandler->DeltaPhi(phiLead,v0.GetPhi()));

					if (mFlagMC) {
						if (v0.IsMCPrimary())
							hV0IMvPtPrimary[iSp]->Fill(v0.GetPt(),v0mass[iSp]);
						if (v0.IsMCPrimary() && v0.GetMCPdgCode() == PDG_IDS[iSp])
							hV0IMvPtPrimaryPDG[iSp]->Fill(v0.GetPt(),v0mass[iSp]);
						if (!v0.IsMCPrimary())
							hV0IMvPtSecondary[iSp]->Fill(v0.GetPt(),v0mass[iSp]);
						if (!v0.IsMCPrimary() && v0.GetMCPdgCode() == PDG_IDS[iSp])
							hV0IMvPtSecondaryPDG[iSp]->Fill(v0.GetPt(),v0mass[iSp]);
						if (!v0.IsMCPrimary() && TMath::Abs(v0.GetMCPrimaryPdgCode()) == 3312)
							hV0IMvPtSecondaryXi[iSp]->Fill(v0.GetPt(),v0mass[iSp]);
						if (v0.GetMCPdgCode() != PDG_IDS[iSp])
							hV0IMvPtBackground[iSp]->Fill(v0.GetPt(),v0mass[iSp]);

						if (iSp == 1 && v0.GetPt() < 0.6 && v0.GetMCPdgCode() != PDG_IDS[iSp]) {
							if (TMath::Abs(v0mass[iSp])<0.01) hK0sLowDaughtersPDG->Fill(v0.GetPosTrackPdg(),v0.GetNegTrackPdg());
							hK0sLowIMvPDG->Fill(v0.GetMCPdgCode(),v0mass[iSp]);
							hK0sLowIMvDaughterPDGPos->Fill(v0.GetPosTrackPdg(),v0mass[iSp]);
							hK0sLowIMvDaughterPDGNeg->Fill(v0.GetNegTrackPdg(),v0mass[iSp]);
						}
						if (iSp == 1 && v0.GetPt() > 3. && v0.GetMCPdgCode() != PDG_IDS[iSp]) {
							if (TMath::Abs(v0mass[iSp])<0.01) hK0sHighDaughtersPDG->Fill(v0.GetPosTrackPdg(),v0.GetNegTrackPdg());
							hK0sHighIMvPDG->Fill(v0.GetMCPdgCode(),v0mass[iSp]);
							hK0sHighIMvDaughterPDGPos->Fill(v0.GetPosTrackPdg(),v0mass[iSp]);
							hK0sHighIMvDaughterPDGNeg->Fill(v0.GetNegTrackPdg(),v0mass[iSp]);
						}
					}

					if (isEventFHM) {
						ProcessV0toHist(v0,iSp,D,V0M,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][V0M][iSph]) ProcessV0toHist(v0,iSp,D,V0M,iSph);
					}

					if (isEventMHM) {
						ProcessV0toHist(v0,iSp,D,NCharged,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][NCharged][iSph]) ProcessV0toHist(v0,iSp,D,NCharged,iSph);
					}


					if (isEventFHM01) {
						ProcessV0toHist(v0,iSp,D,V0M01,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][V0M01][iSph]) ProcessV0toHist(v0,iSp,D,V0M01,iSph);
					}

					if (isEventMHM01) {
						ProcessV0toHist(v0,iSp,D,NCharged01,sphMB);
						for (Int_t iSph = 1; iSph < NSPHERO; iSph++) if (isEventSphero[D][NCharged01][iSph]) ProcessV0toHist(v0,iSp,D,NCharged01,iSph);
					}
					

					if (isEventRT) {
						
						Int_t region = WhatRegion(v0.GetPhi(),phiLead);
						ProcessV0toHistRT(v0,iSp,D,region,nChTrans,nChTransMin);
						
						Int_t regionSide = WhatRegionSide(v0.GetPhi(),phiLead);
						if (!region) ProcessV0toHistRT(v0,iSp,D,IsMinOrMax(isSideAMin,regionSide),nChTrans,nChTransMin);
						
					}

				}
			}		

	}

	hEventMonitor->Fill(7);
	return 0;	
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
				&& ev.GetRefMult() > 0.)		return true;
				//&& ev.GetRefMult() > 9.9)		return true;
			break;

		case 2: 
			if (ev.GetRefMult() < cuts::EV_MHM[1] 
				&& ev.GetRefMult() > cuts::EV_MHM[0])		return true;
			break;

		case 3: 
			if (ev.GetV0MCentrality() < 1. 
				&& ev.GetRefMult() > 0.)		return true;
				//&& ev.GetRefMult() > 9.9)		return true;
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

Bool_t MyAnalysisV0::ProcessV0toHistRT(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Reg, Int_t Nt, Int_t NtMin) {
	
	//printf("V0 p_{T} is %f, args are %i , %i , %i \n", v0.GetPt(), Sp, Mu, Sph);

	hV0PtNt[Sp][Type][Reg]->Fill(v0.GetPt(),Nt);
	hV0EtaNt[Sp][Type][Reg]->Fill(v0.GetEta(),Nt);
	hV0PhiNt[Sp][Type][Reg]->Fill(v0.GetPhi(),Nt);
	hV0PtNtMin[Sp][Type][Reg]->Fill(v0.GetPt(),NtMin);
	
	if (Type==2) return true;
	Double_t v0mass[] 	= {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	hV0IMvPtNt[Sp][Type][Reg]->Fill(v0.GetPt(),v0mass[Sp],Nt);
	hV0IMvPtNtMin[Sp][Type][Reg]->Fill(v0.GetPt(),v0mass[Sp],NtMin);

	return true;	
}

Bool_t MyAnalysisV0::ProcessV0toTree(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu) {
	
	Double_t v0mass[] = {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	Double_t v0massKF[] = {0., v0.GetKFIMK0s(), v0.GetKFIML(), v0.GetKFIMLbar()};
	
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
	if (v0.HasFastSignal() < cuts::V0_FASTSIGNAL)		return false;
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
			//if (TMath::Abs(v0mass[2]) < cuts::K0S_COMP_M) 			return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //14
			//if (TMath::Abs(v0mass[3]) < cuts::K0S_COMP_M) 			return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //15

			if (v0mass[2] < mParMuL->Eval(v0.GetPt()) + cuts::V0_COMP_NSIG*mParSigL->Eval(v0.GetPt())
				&& v0mass[2] > mParMuL->Eval(v0.GetPt()) - cuts::V0_COMP_NSIG*mParSigL->Eval(v0.GetPt())) return false;
			if (v0mass[3] < mParMuL->Eval(v0.GetPt()) + cuts::V0_COMP_NSIG*mParSigL->Eval(v0.GetPt())
				&& v0mass[3] > mParMuL->Eval(v0.GetPt()) - cuts::V0_COMP_NSIG*mParSigL->Eval(v0.GetPt())) return false;

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
			if ((MASSES[Sp]+v0mass[Sp])*v0.GetRadius()/v0.GetPt() > cuts::L_TAU)	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //14
			
			//if (TMath::Abs(v0mass[1]) < cuts::L_COMP_M) 			return false;
			if (v0mass[1] < mParMuK0s->Eval(v0.GetPt()) + cuts::V0_COMP_NSIG*mParSigK0s->Eval(v0.GetPt())
				&& v0mass[1] > mParMuK0s->Eval(v0.GetPt()) - cuts::V0_COMP_NSIG*mParSigK0s->Eval(v0.GetPt())) return false;
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
			if ((MASSES[Sp]+v0mass[Sp])*v0.GetRadius()/v0.GetPt() > cuts::L_TAU)	return false;
			hV0CutIMvPt[Sp][Type][cutN++]->Fill(v0.GetPt(),v0mass[Sp]); //14
			//if (TMath::Abs(v0mass[1]) < cuts::L_COMP_M) 			return false;
			if (v0mass[1] < mParMuK0s->Eval(v0.GetPt()) + cuts::V0_COMP_NSIG*mParSigK0s->Eval(v0.GetPt())
				&& v0mass[1] > mParMuK0s->Eval(v0.GetPt()) - cuts::V0_COMP_NSIG*mParSigK0s->Eval(v0.GetPt())) return false;
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

Bool_t MyAnalysisV0::IsV0VaryCut(MyV0 &v0, Int_t Sp, Int_t Type, Int_t VarCut, Float_t VarVal) {

	// APPLY BASIC CUTS
	
	if (Type>0) if (v0.GetMCPdgCode() != PDG_IDS[Sp])	return false; 
	if (!v0.IsOffline())	return false;
	if (v0.GetEta() < cuts::V0_ETA[0]) 	return false;
	if (v0.GetEta() > cuts::V0_ETA[1]) 	return false;
	if (v0.GetPt() < cuts::V0_PT[0]) 	return false;
	if (v0.GetPt() > cuts::V0_PT[1]) 	return false;
	if (Sp>1) {
		if (v0.GetPt() < cuts::L_PT[0]) 	return false;
		if (v0.GetPt() > cuts::L_PT[1]) 	return false;
	}
	
	MyTrack trP(v0.GetPosTrack());
	trP.SetHandler(mHandler); 
	MyTrack trN(v0.GetNegTrack());
	trN.SetHandler(mHandler);
	
	if (trP.GetEta() < cuts::K0S_D_ETA[0])	return false;
	if (trP.GetEta() > cuts::K0S_D_ETA[1])	return false;
	if ( cuts::V0_D_GOODTRACK && !trP.IsGoodV0daughter())	return false;
	if (trN.GetEta() < cuts::K0S_D_ETA[0])	return false;
	if (trN.GetEta() > cuts::K0S_D_ETA[1])	return false;
	if ( cuts::V0_D_GOODTRACK && !trN.IsGoodV0daughter())	return false;

	// APPLY ALL OTHER CUTS FIRST EXCEPT VARIED
	for (int iCut = 0; iCut < cutsSizeof; ++iCut)	{

		//cout << "Gothereee " << iCut << endl;
		if (iCut == VarCut) continue;

		//if (iCut == 7) cout << Sp << " " << Type << " " << mCutsDir[Sp][Type][iCut] << " " << v0.GetCut(iCut) << " " << mCutsVal[Sp][iCut] << " " << mParMuK0s->Eval(v0.GetPt()) + mCutsVal[Sp][iCut]*mParSigK0s->Eval(v0.GetPt()) << endl;

		if (iCut == CompMassK0s) {
			if (mCutsDir[Sp][Type][iCut]!=0) if ( (v0.GetIMK0s() > mParMuK0s->Eval(v0.GetPt()) - mCutsVal[Sp][iCut]*mParSigK0s->Eval(v0.GetPt()))
			&&  (v0.GetIMK0s() < mParMuK0s->Eval(v0.GetPt()) + mCutsVal[Sp][iCut]*mParSigK0s->Eval(v0.GetPt())) ) return false;
			continue;	}
		if (iCut == CompMassL && VarCut != CompMassLbar) {
			if (mCutsDir[Sp][Type][iCut]!=0) if ( (v0.GetIML() > mParMuL->Eval(v0.GetPt()) - mCutsVal[Sp][iCut]*mParSigL->Eval(v0.GetPt()) )
			&&  (v0.GetIML() < mParMuL->Eval(v0.GetPt()) + mCutsVal[Sp][iCut]*mParSigL->Eval(v0.GetPt()) ) )	return false;
			continue;	}
		if (iCut == CompMassLbar && VarCut != CompMassL) {
			if (mCutsDir[Sp][Type][iCut]!=0) if ( (v0.GetIMLbar() > mParMuL->Eval(v0.GetPt()) - mCutsVal[Sp][iCut]*mParSigL->Eval(v0.GetPt()) )
			&&  (v0.GetIMLbar() < mParMuL->Eval(v0.GetPt()) + mCutsVal[Sp][iCut]*mParSigL->Eval(v0.GetPt()) ) )	return false;
			continue;	}

		//if (iCut == 12 || (iCut > 15 && iCut < 33)) continue;
		//if (iCut == 17) cout << Sp << " " << Type << " " << mCutsDir[Sp][Type][iCut] << " " << v0.GetCut(iCut) << " " << mCutsVal[Sp][iCut] << endl;
		
		if (mCutsDir[Sp][Type][iCut]>0) if (v0.GetCut(iCut) > mCutsVal[Sp][iCut]) return false;
		if (mCutsDir[Sp][Type][iCut]<0) if (v0.GetCut(iCut) < mCutsVal[Sp][iCut]) return false; 
		
	}

	

	// APPLY VARIED CUT
	if (VarCut == CompMassK0s) {
		if ( (v0.GetIMK0s() > mParMuK0s->Eval(v0.GetPt()) - VarVal*mParSigK0s->Eval(v0.GetPt()) ) 
			&& (v0.GetIMK0s() < mParMuK0s->Eval(v0.GetPt()) + VarVal*mParSigK0s->Eval(v0.GetPt()) ) ) return false;
		return true;	}
	if (VarCut == CompMassL) {
		if ( (v0.GetIML() > mParMuL->Eval(v0.GetPt()) - VarVal*mParSigL->Eval(v0.GetPt()) ) 
			&& (v0.GetIML() < mParMuL->Eval(v0.GetPt()) + VarVal*mParSigL->Eval(v0.GetPt()) ) ) return false;
		return true;	}
	if (VarCut == CompMassLbar) {
		if ( (v0.GetIMLbar() > mParMuL->Eval(v0.GetPt()) - VarVal*mParSigL->Eval(v0.GetPt()) ) 
			&& (v0.GetIMLbar() < mParMuL->Eval(v0.GetPt()) + VarVal*mParSigL->Eval(v0.GetPt()) ) ) return false;
		return true;	}
				
	
	if (mCutsDir[Sp][Type][VarCut]>0) if (v0.GetCut(VarCut) > VarVal) return false;
	if (mCutsDir[Sp][Type][VarCut]<0) if (v0.GetCut(VarCut) < VarVal) return false;

	return true;

}

void MyAnalysisV0::ProcessV0SystVar(MyV0 &v0, Int_t Sp, Int_t Type, Int_t Mu, Int_t Sph) {

	Double_t v0mass[] 	= {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};

	for (Int_t iV = 0; iV < sysVarsSizeof; iV++) {
	
		if (IsV0VaryCut(v0,Sp,Type,RadiusL,sysVar[Sp-1][sysRadiusL][iV])) hV0IMvPtSys[Sp][Mu][Sph][sysRadiusL][iV]->Fill(v0.GetPt(),v0mass[Sp]);
		if (IsV0VaryCut(v0,Sp,Type,DCAdd,sysVar[Sp-1][sysDCAdd][iV])) hV0IMvPtSys[Sp][Mu][Sph][sysDCAdd][iV]->Fill(v0.GetPt(),v0mass[Sp]);
		if (IsV0VaryCut(v0,Sp,Type,CPA,sysVar[Sp-1][sysCPA][iV])) hV0IMvPtSys[Sp][Mu][Sph][sysCPA][iV]->Fill(v0.GetPt(),v0mass[Sp]);
		if (IsV0VaryCut(v0,Sp,Type,FastSignal,sysVar[Sp-1][sysFastSignal][iV])) hV0IMvPtSys[Sp][Mu][Sph][sysFastSignal][iV]->Fill(v0.GetPt(),v0mass[Sp]);

		if ( (Sp==1 && IsV0VaryCut(v0,Sp,Type,CompMassL,sysVar[Sp-1][sysCompMass][iV])
					 && IsV0VaryCut(v0,Sp,Type,CompMassLbar,sysVar[Sp-1][sysCompMass][iV]) )
			|| (Sp>1 && IsV0VaryCut(v0,Sp,Type,CompMassK0s,sysVar[Sp-1][sysCompMass][iV]) ) )
			hV0IMvPtSys[Sp][Mu][Sph][sysCompMass][iV]->Fill(v0.GetPt(),v0mass[Sp]);			

		if ( (Sp==2 && IsV0VaryCut(v0,Sp,Type,LifetimeL,sysVar[Sp-1][sysLifetime][iV]) )
			|| (Sp==3 && IsV0VaryCut(v0,Sp,Type,LifetimeLbar,sysVar[Sp-1][sysLifetime][iV]) ) )
			hV0IMvPtSys[Sp][Mu][Sph][sysLifetime][iV]->Fill(v0.GetPt(),v0mass[Sp]);

		if ( (Sp==1 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCposPiL,-1.*sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCposPiH,sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCnegPiL,-1.*sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCnegPiH,sysVar[Sp-1][sysNSigmaTPC][iV])		)
		  || (Sp==2 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCposPrL,-1.*sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCposPrH,sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCnegPiL,-1.*sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCnegPiH,sysVar[Sp-1][sysNSigmaTPC][iV])		)
		  || (Sp==3 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCposPiL,-1.*sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCposPiH,sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCnegPrL,-1.*sysVar[Sp-1][sysNSigmaTPC][iV])
					 && IsV0VaryCut(v0,Sp,Type,NSigmaTPCnegPrH,sysVar[Sp-1][sysNSigmaTPC][iV])		)	)
			hV0IMvPtSys[Sp][Mu][Sph][sysNSigmaTPC][iV]->Fill(v0.GetPt(),v0mass[Sp]);

		if (IsV0VaryCut(v0,Sp,Type,DCAPVpos,sysVar[Sp-1][sysDCAPVpos][iV])) hV0IMvPtSys[Sp][Mu][Sph][sysDCAPVpos][iV]->Fill(v0.GetPt(),v0mass[Sp]);
		if (IsV0VaryCut(v0,Sp,Type,DCAPVneg,sysVar[Sp-1][sysDCAPVneg][iV])) hV0IMvPtSys[Sp][Mu][Sph][sysDCAPVneg][iV]->Fill(v0.GetPt(),v0mass[Sp]);

		if (IsV0VaryCut(v0,Sp,Type,NClusterpos,sysVar[Sp-1][sysNCluster][iV])
			&& IsV0VaryCut(v0,Sp,Type,NClusterneg,sysVar[Sp-1][sysNCluster][iV]) ) 
			hV0IMvPtSys[Sp][Mu][Sph][sysNCluster][iV]->Fill(v0.GetPt(),v0mass[Sp]);

		if (IsV0VaryCut(v0,Sp,Type,NClusterFpos,sysVar[Sp-1][sysNClusterF][iV])
			&& IsV0VaryCut(v0,Sp,Type,NClusterFneg,sysVar[Sp-1][sysNClusterF][iV]) ) 
			hV0IMvPtSys[Sp][Mu][Sph][sysNClusterF][iV]->Fill(v0.GetPt(),v0mass[Sp]);
	}

}

Bool_t MyAnalysisV0::IsTrans(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	if (TMath::Abs(dphi) < TMath::Pi()/3.) 		return false;
	if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return false;
	
	return true;
}

Int_t MyAnalysisV0::WhatRegion(Double_t phi1, Double_t phiTrig) {
	// 0: transverse, 1: toward, 2: away
	
	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	
	if (TMath::Abs(dphi) < TMath::Pi()/3.)			return 1;
	else if (TMath::Abs(dphi) > 2.*TMath::Pi()/3.) 	return 2;
	else return 0;
}

Int_t MyAnalysisV0::WhatRegionSide(Double_t phi1, Double_t phiTrig) {

	Double_t dphi = mHandler->DeltaPhi(phi1,phiTrig);
	if (!IsTrans(phi1,phiTrig)) return 0;
	
	if (dphi>0) return 1; // region A
	if (dphi<0) return 2; // region B
	return 0;
}

Int_t MyAnalysisV0::IsMinOrMax(Bool_t isA, Int_t reg) {

	Int_t res = 4;
	if (isA && reg == 1)	res = 3;
	if (!isA && reg == 1)	res = 4;
	if (isA && reg == 2)	res = 4;
	if (!isA && reg == 2)	res = 3;
	return res;
}

Bool_t MyAnalysisV0::SelectV0Daughter(MyTrack &tr) {

	if (tr.GetEta() < cuts::K0S_D_ETA[0])	return false;
	if (tr.GetEta() > cuts::K0S_D_ETA[1])	return false;
	if (TMath::Abs(tr.GetDCApvXY()) < cuts::V0_D_DCAPVXY)	return false;
	//if (!tr.IsITSTPC2011())	return false;

	if ( cuts::V0_D_GOODTRACK && !tr.IsGoodV0daughter())	return false;

	//printf("values are %f and %f \n", tr.GetTPCnc(),tr.GetTPCnc()/tr.GetTPCNclsF());

	if (tr.GetTPCnc() < cuts::V0_D_NCR)							return false;
	if (tr.GetTPCNclsF() > 0 && (tr.GetTPCnc()/tr.GetTPCNclsF()) < cuts::V0_D_NRATIO)		return false;

	return true;
}

Bool_t MyAnalysisV0::SelectParticle(MyParticle &p) {

	if (p.GetPdgCode() != PDG_IDS[1]
		&& p.GetPdgCode() != PDG_IDS[2]
		&& p.GetPdgCode() != PDG_IDS[3]
		&& p.GetPdgCode() != 3312			// also Xi pm for feed-down study
		&& p.GetPdgCode() != -3312 
		&& p.GetPdgCode() != 3322			// also Xi0 for feed-down study
		&& p.GetPdgCode() != -3322 )		return false;

	if (p.GetPdgCode() == -3312 || p.GetPdgCode() == 3312
			|| p.GetPdgCode() == -3322 || p.GetPdgCode() == 3322 ) {
		if (p.GetY() < cuts::V0_ETA[0]) 		return false;
		if (p.GetY() > cuts::V0_ETA[1]) 		return false;
	} else {
		if (p.GetEta() < cuts::V0_ETA[0]) 		return false;
		if (p.GetEta() > cuts::V0_ETA[1]) 		return false;
	}
	//if (p.GetY() < cuts::V0_Y[0]) 		return false;
	//if (p.GetY() > cuts::V0_Y[1]) 		return false;
	
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

	Int_t nType = (mFlagMC) ? 2 : 1;

	// MONITORS
	hEventMonitor 			= new TH1D("hEventMonitor","; Step; Entries",10,-0.5,9.5);
	hTrackMonitor 			= new TH1D("hTrackMonitor","; Step; Entries",10,-0.5,9.5);
	hV0Monitor  			= new TH1D("hV0Monitor","; Step; Entries",10,-0.5,9.5);
	hParticleMonitor 		= new TH1D("hParticleMonitor","; Step; Entries",10,-0.5,9.5);

	cout << "Created " << hEventMonitor << endl;

	// EVENT INFO HISTOGRAMS
	hEventCuts	 			= new TH1D("hEventCuts","; Step; Entries",10,-0.5,9.5);
	hEventVz				= new TH1D("hEventVz","; vz; Entries",400,-50,50);
	hEventV0MCentrality		= new TH1D("hEventV0MCentrality","; V0M Centrality; Entries",300,0,150);
	hEventRefMult			= new TH1D("hEventRefMult","; Reference multiplicity; Entries",150,0,150);
	hEventV0MCentvRefMult	= new TH2D("hEventV0MCentvRefMult","; Reference multiplicity; V0M Centrality; Entries"
		,150,0,150,300,0,150);

	hEventType					= new TH1D("hEventType",";; Counts", NEVENTTYPES, 0, NEVENTTYPES);
	for (int iBin = 1; iBin <= NEVENTTYPES; iBin++) {
		hEventType->GetXaxis()->SetBinLabel(iBin,EVENTTYPES[iBin-1]); }
	hEventSpherocityV0M			= new TH1D("hEventSpherocityV0M","; S_{O}; Entries",400,-0.1,1.1);
	hEventSpherocityNCharged	= new TH1D("hEventSpherocityNCharged","; S_{O}; Entries",400,-0.1,1.1);
	hEventSpherocityV0M01		= new TH1D("hEventSpherocityV0M01","; S_{O}; Entries",400,-0.1,1.1);
	hEventSpherocityNCharged01	= new TH1D("hEventSpherocityNCharged01","; S_{O}; Entries",400,-0.1,1.1);
	hEventTSMCvRC 			= new TH2D("hEventTSMCvRC",";S_{O} with MC particles; S_{O} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSNormMCvRC		= new TH2D("hEventTSNormMCvRC",";S_{O}^{p_{T}=1} with MC particles; S_{O}^{p_{T}=1} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSMCvNorm			= new TH2D("hEventTSMCvNorm",";S_{O} with MC particles; S_{O}^{p_{T}=1} with MC particles", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventTSRCvNorm			= new TH2D("hEventTSRCvNorm",";S_{O} with RC tracks; S_{O}^{p_{T}=1} with RC tracks", 440, -1.1, 1.1, 440, -1.1, 1.1);
	hEventMultvSpheroD		= new TH2D("hEventMultvSpheroD","", NSPHERO, 0, NSPHERO, NMULTI, 0, NMULTI);
	for (int iBin = 1; iBin <= NSPHERO; iBin++) hEventMultvSpheroD->GetXaxis()->SetBinLabel(iBin,SPHERO[iBin-1]);
	for (int iBin = 1; iBin <= NMULTI; iBin++) hEventMultvSpheroD->GetYaxis()->SetBinLabel(iBin,MULTI[iBin-1]);
	hEventMultvSpheroMC		= new TH2D("hEventMultvSpheroMC","", NSPHERO, 0, NSPHERO, NMULTI, 0, NMULTI);
	for (int iBin = 1; iBin <= NSPHERO; iBin++) hEventMultvSpheroMC->GetXaxis()->SetBinLabel(iBin,SPHERO[iBin-1]);
	for (int iBin = 1; iBin <= NMULTI; iBin++) hEventMultvSpheroMC->GetYaxis()->SetBinLabel(iBin,MULTI[iBin-1]);		


	hLeadPhivPt				= new TH2D("hLeadPhivPt","; p_{T} (GeV/#it{c}); #phi", 200, 0., 30., 400, -0.2, 6.4);
	hNchvLeadPt				= new TH1D("hNchvLeadPt","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 200, 0., 30.);
	hNchvLeadPt2			= new TH2D("hNchvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.]", 90, 0., 30.,50,-0.5,49.5);
	hNchMinvLeadPt2			= new TH2D("hNchMinvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.,min]", 90, 0., 30.,50,-0.5,49.5);
	hNchMaxvLeadPt2			= new TH2D("hNchMaxvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); N_{ch} [trans.,max]", 90, 0., 30.,50,-0.5,49.5);
	hMeanPtvLeadPt2			= new TH2D("hMeanPtvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); <pT> [trans.]", 90, 0., 30.,50,0.,2.5);
	hMeanPtMinvLeadPt2		= new TH2D("hMeanPtMinvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); <pT> [trans.,min]", 90, 0., 30.,50,0.,2.5);
	hMeanPtMaxvLeadPt2		= new TH2D("hMeanPtMaxvLeadPt2","; p_{T}^{leading} (GeV/#it{c}); <pT> [trans.,max]", 90, 0., 30.,50,0.,2.5);
	hNchTrans				= new TH1D("hNchTrans","; N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransMC 			= new TH1D("hNchTransMC","; MC N_ch [trans.]; Entries",50, -0.5, 49.5);
	hNchTransRCvMC			= new TH2D("hNchTransRCvMC", ";MC N_ch [trans.]; RC N_ch [trans.]",50,-0.5,49.5,50,-0.5,49.5);
	hNtvNtMin				= new TH2D("hNtvNtMin","; N_ch [trans.,min]; N_ch [trans.]",50, -0.5, 49.5,50, -0.5, 49.5);
	hNtvNtMax				= new TH2D("hNtvNtMax","; N_ch [trans.,max]; N_ch [trans.]",50, -0.5, 49.5,50, -0.5, 49.5);
	hNtMaxvNtMin			= new TH2D("hNtMaxvNtMin","; N_ch [trans.,max]; N_ch [trans.,min]",50, -0.5, 49.5,50, -0.5, 49.5);
	hRt						= new TH1D("hRt","; R_{T}; Entries",4000, -0.02, 5.02);//4.975);
	hRtMC					= new TH1D("hRtMC","; MC R_{T}; Entries",4000, -0.02, 5.02);//4.975);
	hRtRCvMC 				= new TH2D("hRtRCMC","; MC R_{T}; R_{T}", 200, -0.02, 5.02, 200, -0.02, 5.02);
 	hLeadPtvNchTrans0		= new TH2D("hLeadPtvNchTrans0","; N_{ch}^{trans}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hLeadPtvNchTrans		= new TH2D("hLeadPtvNchTrans","; N_{ch}^{trans}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hLeadPtvNchTransMin		= new TH2D("hLeadPtvNchTransMin","; N_{ch}^{trans,min}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hLeadPtvNchTransMax		= new TH2D("hLeadPtvNchTransMax","; N_{ch}^{trans,max}; p_{T}^{leading} (GeV/#it{c})", 50, -0.5, 49.5, 200, 0., 40.);
 	hNchTransvSpherocityV0M			= new TH2D("hNchTransvSpherocityV0M",";S_{O};N_{ch}^{trans}",400,-0.1,1.1,50,-0.5,49.5);
 	hNchTransvSpherocityNCharged	= new TH2D("hNchTransvSpherocityNCharged",";S_{O};N_{ch}^{trans}",400,-0.1,1.1,50,-0.5,49.5);
 	hV0DPhivNchTrans				= new TH2D("hV0DPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hTrackDPhivNchTrans				= new TH2D("hTrackDPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hTrackDPhivNchTransMin			= new TH2D("hTrackDPhivNchTransMin","; N_{ch}^{trans,min}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hParticleDPhivNchTrans			= new TH2D("hParticleDPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hK0sDPhivNchTrans				= new TH2D("hK0sDPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hLDPhivNchTrans					= new TH2D("hLDPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hLbarDPhivNchTrans				= new TH2D("hLbarDPhivNchTrans","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hK0sDPhivNchTransMC				= new TH2D("hK0sDPhivNchTransMC","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hLDPhivNchTransMC				= new TH2D("hLDPhivNchTransMC","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);
 	hLbarDPhivNchTransMC			= new TH2D("hLbarDPhivNchTransMC","; N_{ch}^{trans}; #phi - #phi^{lead}", 50, -0.5, 49.5, 300, -3.2, 3.2);

	// TRACK HISTOGRAMS
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
		
		if (iMu == 0 && iSph > 0) continue;		
		
		hTrackPt[iType][iMu][iSph]			= new TH1D(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";track p_{T} (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hTrackEtavPhi[iType][iMu][iSph]		= new TH2D(Form("hTrackEtavPhi_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";track #phi; track #eta", 400, -0.2, 6.4, 400, -1., 1.);		

	} } } 

	// MC PARTICLE HISTOGRAMS
	hParticlePrimaryvPDG				= new TH2D("hParticlePrimaryvPDG", "; PDG id; Primary", 10000,-5000,5000,2,-0.5,1.5);
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
		hProtonNchTransvPt[iReg]		= new TH2D(Form("hProtonNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS, XBINS, 50, -0.5, 49.5);
		hPionNchTransvPt[iReg]			= new TH2D(Form("hPionNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS, XBINS, 50, -0.5, 49.5);
		hLambdaNchTransvPt[iReg]		= new TH2D(Form("hLambdaNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS, XBINS, 50, -0.5, 49.5);
		hK0sNchTransvPt[iReg]			= new TH2D(Form("hK0sNchTransvPt_%s",REGIONS[iReg]),"; p_{T} (GeV/#it{c}); N_{ch}^{trans}", NPTBINS, XBINS, 50, -0.5, 49.5);
	}

	// V0 HISTOGRAMS
	hV0Radius		= new TH1D("hV0Radius",";V0 radius; Entries",400,0.,150.);
	hV0ProperT		= new TH1D("hV0ProperT",";V0 proper #tau; Entries",400,0.,150.);
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		if (iMu == 0 && iSph > 0) continue;		

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
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);

	} } } }

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

	// V0 RT HISTOGRAMS
	const int NBINS_IM = 1000;
	double XBINS_IM[NBINS_IM+1];// = {};
	double binStepIM = (0.2 + 0.2)/(double)NBINS_IM;
	for (int i=0; i <= NBINS_IM; i++) XBINS_IM[i] = -0.2 + (double)i*binStepIM;
	const int NBINS_NT = 60;
	double XBINS_NT[NBINS_NT+1];// = {};
	double binStepNT = (0.5 + 59.5)/(double)NBINS_NT;
	for (int i=0; i <= NBINS_NT; i++) XBINS_NT[i] = -0.5 + (double)i*binStepNT;
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
				

		hV0PtNt[iSp][iType][iReg]			= new TH2D(Form("hV0PtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 p_{T} (GeV/#it{c}); N_{ch} [trans.]; Entries",	NPTBINS,XBINS, 60, -0.5, 59.5);
		hV0PtNtMin[iSp][iType][iReg]			= new TH2D(Form("hV0PtNtMin_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 p_{T} (GeV/#it{c}); N_{ch} [trans.,min]; Entries",	NPTBINS,XBINS, 60, -0.5, 59.5);
		
		hV0EtaNt[iSp][iType][iReg]			= new TH2D(Form("hV0EtaNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 #eta; N_{ch} [trans.]; Entries", 		200, -1., 1., 60, -0.5, 59.5);
		hV0PhiNt[iSp][iType][iReg]			= new TH2D(Form("hV0PhiNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 #eta; N_{ch} [trans.]; Entries", 		200, -3.2, 3.2, 60, -0.5, 59.5);

		hV0IMvPtNt[iSp][iType][iReg]		= new TH3D(Form("hV0IMvPtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); N_{ch} [trans.]; Entries",		NPTBINS, XBINS, NBINS_IM, XBINS_IM, NBINS_NT, XBINS_NT);
		hV0IMvPtNtMin[iSp][iType][iReg]		= new TH3D(Form("hV0IMvPtNtMin_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); N_{ch} [trans.,min]; Entries",		NPTBINS, XBINS, NBINS_IM, XBINS_IM, NBINS_NT, XBINS_NT);

	} } } 


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
	const int NBINS_FDMASS = 1000;
	double XBINS_FDMASS[NBINS_FDMASS+1];// = {};
	double binStep = (0.2 + 0.2)/(double)NBINS_FDMASS;
	for (int i=0; i <= NBINS_FDMASS; i++) XBINS_FDMASS[i] = -0.2 + (double)i*binStep;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

		hV0Feeddown[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Feeddown_%s",SPECIES[iSp]));
		hV0FeeddownPDG[iSp] =	new TH1D(Form("hV0FeeddownPDG_%s",SPECIES[iSp]),";PDG ID;Entries",20000,-10000,10000);
		hV0FeeddownMatrix[iSp]	= new TH3D(Form("hV0FeeddownMatrix_%s",SPECIES[iSp]),";primary grandmother p_{T}; decay V0 p_{T}; decay V0 mass", NXIPTBINS, XIXBINS, NPTBINS, XBINS, NBINS_FDMASS, XBINS_FDMASS);
		hV0FeeddownMotherPt[iSp] 	= new TH1D(Form("hV0FeeddownMotherPt_%s",SPECIES[iSp]),";primary grandmother p_{T}; Entries", NXIPTBINS, XIXBINS);
		hV0FeeddownMatrixXi0[iSp]	= new TH3D(Form("hV0FeeddownMatrixXi0_%s",SPECIES[iSp]),";primary grandmother p_{T}; decay V0 p_{T}; decay V0 mass", NXIPTBINS, XIXBINS, NPTBINS, XBINS, NBINS_FDMASS, XBINS_FDMASS);
		hV0FeeddownMotherPtXi0[iSp] 	= new TH1D(Form("hV0FeeddownMotherPtXi0_%s",SPECIES[iSp]),";primary grandmother p_{T}; Entries", NXIPTBINS, XIXBINS);
		
	}

	// V0 CUTS STUDY HISTOGRAMS

	// SIGNAL LOSS STUDY
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
		hV0PtNoTrigger[iSp]			= new TH1D(Form("hV0PtNoTrigger_%s",SPECIES[iSp]),
			";V0 p_{T} (GeV/#it{c}); Entries",								NPTBINS,XBINS);
	}



	// DETAILED MC CLOSURE STUDY
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
		hV0IMvPtPrimary[iSp]		= new TH2D(Form("hV0IMvPtPrimary_%s",SPECIES[iSp]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtPrimaryPDG[iSp]		= new TH2D(Form("hV0IMvPtPrimaryPDG_%s",SPECIES[iSp]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSecondary[iSp]		= new TH2D(Form("hV0IMvPtSecondary_%s",SPECIES[iSp]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSecondaryPDG[iSp]		= new TH2D(Form("hV0IMvPtSecondaryPDG_%s",SPECIES[iSp]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSecondaryXi[iSp]		= new TH2D(Form("hV0IMvPtSecondaryXi_%s",SPECIES[iSp]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtBackground[iSp]		= new TH2D(Form("hV0IMvPtBackground_%s",SPECIES[iSp]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
	}
	hK0sLowDaughtersPDG			= new TH2D(Form("hK0sLowDaughtersPDG"),";daughter+ PDG ID; daughter- PDG ID; Entries",
		6000,-3000,3000, 6000,-3000,3000);
	hK0sLowIMvPDG 				= new TH2D(Form("hK0sLowIMvPDG"),";PDG ID; V0 m (GeV/#it{c}^{2}); Entries",
		20000,-10000,10000, 200, -0.1, 0.1);
	hK0sLowIMvDaughterPDGPos	= new TH2D(Form("hK0sLowIMvDaughterPDGPos"),";daughter+ PDG ID; V0 m (GeV/#it{c}^{2}); Entries",
		20000,-10000,10000, 200, -0.1, 0.1);
	hK0sLowIMvDaughterPDGNeg	= new TH2D(Form("hK0sLowIMvDaughterPDGNeg"),";daughter- PDG ID; V0 m (GeV/#it{c}^{2}); Entries",
		20000,-10000,10000, 200, -0.1, 0.1);

	hK0sHighDaughtersPDG		= new TH2D(Form("hK0sHighDaughtersPDG"),";daughter+ PDG ID; daughter- PDG ID; Entries",
		6000,-3000,3000, 6000,-3000,3000);
	hK0sHighIMvPDG 				= new TH2D(Form("hK0sHighIMvPDG"),";PDG ID; V0 m (GeV/#it{c}^{2}); Entries",
		20000,-10000,10000, 200, -0.1, 0.1);
	hK0sHighIMvDaughterPDGPos	= new TH2D(Form("hK0sHighIMvDaughterPDGPos"),";daughter+ PDG ID; V0 m (GeV/#it{c}^{2}); Entries",
		20000,-10000,10000, 200, -0.1, 0.1);
	hK0sHighIMvDaughterPDGNeg	= new TH2D(Form("hK0sHighIMvDaughterPDGNeg"),";daughter- PDG ID; V0 m (GeV/#it{c}^{2}); Entries",
		20000,-10000,10000, 200, -0.1, 0.1);


	


	// V0 NTUPLES TO ALSO ALLOW POST-PROCESSING REBINNING
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
		tV0PtMCMB[iSp]		= new TNtuple(Form("tV0PtMCMB_%s",SPECIES[iSp]),"v0 MC MB pt tree","lPt");
		tV0massRCMB[iSp]	= new TNtuple(Form("tV0massRCMB_%s",SPECIES[iSp]),"v0 RC MB mass tree","MassDT:lPt");
	}


	// SYSTEMATICS STUDY
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
		hV0IMvRadiusL[iSp]	= new TH2D(Form("hV0IMvRadiusL_%s",SPECIES[iSp]),";radius (cm); V0 m (GeV/#it{c}^{2}); Entries",
			100, 0., 1., 1000, -0.2, 0.2);
		hV0IMvDCAdd[iSp]	= new TH2D(Form("hV0IMvDCAdd_%s",SPECIES[iSp]),";DCA of daughters (cm); V0 m (GeV/#it{c}^{2}); Entries",
			100, 0., 2., 1000, -0.2, 0.2);
		hV0IMvCPA[iSp]	= (iSp==1) ? new TH2D(Form("hV0IMvCPA_%s",SPECIES[iSp]),";cos(PA); V0 m (GeV/#it{c}^{2}); Entries",
			100, 0.95, 1., 1000, -0.2, 0.2)
								 :	 new TH2D(Form("hV0IMvCPA_%s",SPECIES[iSp]),";cos(PA); V0 m (GeV/#it{c}^{2}); Entries",
			100, 0.99, 1., 1000, -0.2, 0.2);
		hV0IMvFastSignal[iSp]	= new TH2D(Form("hV0IMvFastSignal_%s",SPECIES[iSp]),";#daughters w. fast signal; V0 m (GeV/#it{c}^{2}); Entries",
			3, 1., 4., 1000, -0.2, 0.2);
		hV0IMvCompMass[iSp]	= new TH2D(Form("hV0IMvCompMass_%s",SPECIES[iSp]),";n#sigma for comp. mass rejection; V0 m (GeV/#it{c}^{2}); Entries",
			100, 2., 6., 1000, -0.2, 0.2);
		if (iSp>1) hV0IMvLifetime[iSp]	= new TH2D(Form("hV0IMvLifetime_%s",SPECIES[iSp]),";Lifetime estimate from XY plane; V0 m (GeV/#it{c}^{2}); Entries",
			100, 10., 50., 1000, -0.2, 0.2);
		hV0IMvNSigmaTPC[iSp]	= new TH2D(Form("hV0IMvNSigmaTPC_%s",SPECIES[iSp]),";n#sigma_{TPC} for daughters; V0 m (GeV/#it{c}^{2}); Entries",
			100, 1., 7., 1000, -0.2, 0.2);
		hV0IMvDCAPVpos[iSp]	= new TH2D(Form("hV0IMvDCAPVpos_%s",SPECIES[iSp]),";DCA PV-pos.daughter (cm); V0 m (GeV/#it{c}^{2}); Entries",
			100, 0.04, 0.1, 1000, -0.2, 0.2);
		hV0IMvDCAPVneg[iSp]	= new TH2D(Form("hV0IMvDCAPVneg_%s",SPECIES[iSp]),";DCA PV-neg.daughter (cm); V0 m (GeV/#it{c}^{2}); Entries",
			100, 0.04, 0.1, 1000, -0.2, 0.2);
		hV0IMvNCluster[iSp]	= new TH2D(Form("hV0IMvNCluster_%s",SPECIES[iSp]),";# crossed rows; V0 m (GeV/#it{c}^{2}); Entries",
			25, 65, 90, 1000, -0.2, 0.2);
		hV0IMvNClusterF[iSp]	= new TH2D(Form("hV0IMvNClusterF_%s",SPECIES[iSp]),";# crossed rows / # findable; V0 m (GeV/#it{c}^{2}); Entries",
			25, 0.75, 1., 1000, -0.2, 0.2);

	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)			{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{


		if (!mFlagMC) if (iMu==0 && iSph>0) continue;
		if (mFlagMC) if (iMu!=0) continue;		// efficiency is varied only in MB
		if (mFlagMC) if (iSph!=0) continue;


		hV0IMvPtSys[iSp][iMu][iSph][sysRadiusL][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysRadiusL],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAdd][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAdd],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysCPA][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysCPA],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysFastSignal][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysFastSignal],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysCompMass][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysCompMass],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysLifetime][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysLifetime],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysNSigmaTPC][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNSigmaTPC],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVpos][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAPVpos],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVneg][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAPVneg],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysNCluster][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNCluster],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);
		hV0IMvPtSys[iSp][iMu][iSph][sysNClusterF][iVar] = new TH2D(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNClusterF],SYSTVAR[iVar]),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.2, 0.2);

	}	}	}	}

	printf("Histograms in analysis %s created \n", 
		this->GetName());

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
				
		//if (iMu > 4 && (iSph < 3 && iSph)) continue;
		//if (iMu < 5 && iSph > 2) continue; 
		hTrackPt[iType][iMu][iSph]	= (TH1D*)mDirFile->Get(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hTrackEtavPhi[iType][iMu][iSph]	= (TH2D*)mDirFile->Get(Form("hTrackEtavPhi_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));

	} } } 

	// MC PARTICLE HISTOGRAMS

	// V0 HISTOGRAMS
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
		//if (iMu > 4 && (iSph < 3 && iSph)) continue;
		//if (iMu < 5 && iSph > 2) continue; 
		hV0Pt[iSp][iType][iMu][iSph]			= (TH1D*)mDirFile->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0Eta[iSp][iType][iMu][iSph]			= (TH1D*)mDirFile->Get(Form("hV0Eta_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0IMvPt[iSp][iType][iMu][iSph]			= (TH2D*)mDirFile->Get(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		// actually perhaps only histos which need to be processed within this analysis need to get fetched
		
	} } } }	

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

		hV0Feeddown[iSp] 	= (TH1D*)mDirFile->Get(Form("hV0Feeddown_%s",SPECIES[iSp]));
		hV0FeeddownPDG[iSp] = (TH1D*)mDirFile->Get(Form("hV0FeeddownPDG_%s",SPECIES[iSp]));
		hV0FeeddownMatrix[iSp]		= (TH3D*)mDirFile->Get(Form("hV0FeeddownMatrix_%s",SPECIES[iSp]));
		hV0FeeddownMotherPt[iSp]	= (TH1D*)mDirFile->Get(Form("hV0FeeddownMotherPt_%s",SPECIES[iSp]));
		hV0FeeddownMatrixXi0[iSp]		= (TH3D*)mDirFile->Get(Form("hV0FeeddownMatrixXi0_%s",SPECIES[iSp]));
		hV0FeeddownMotherPtXi0[iSp]	= (TH1D*)mDirFile->Get(Form("hV0FeeddownMotherPtXi0_%s",SPECIES[iSp]));

		hV0PtNoTrigger[iSp]	= (TH1D*)mDirFile->Get(Form("hV0PtNoTrigger_%s",SPECIES[iSp]));
		cout << "is " << hV0PtNoTrigger[iSp] << endl;
	}

	Int_t nType = (mFlagMC) ? 2 : 1;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
		tV0PtMCMB[iSp]		= (TNtuple*)mDirFile->Get(Form("tV0PtMCMB_%s",SPECIES[iSp]));
		tV0massRCMB[iSp]	= (TNtuple*)mDirFile->Get(Form("tV0massRCMB_%s",SPECIES[iSp]));
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

	//hRt2->SetBins(nbins,rtbins);	// files cannot be hadded afterwards



	if (mFlagMC) DoEfficiency();

	if (mFlagMC) {
		for (int iSp = 1; iSp < NSPECIES; ++iSp)		{	
			hV0PtNoTrigger[iSp]->Divide(hV0Pt[iSp][2][0][0]);
		}
	}

	printf("Analysis %s finished.\n",this->GetName());
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


TH2D* MyAnalysisV0::RebinTH2(TH2D* h) {

	if (!h) return 0x0;
	if (h->GetNbinsX() != NPTBINS) {
			
		TAxis *xaxis = h->GetXaxis(); 
		TAxis *yaxis = h->GetYaxis(); 
		TH2D *htmp = new TH2D("htmp",
			Form(";%s;%s;%s",xaxis->GetTitle(),yaxis->GetTitle(),h->GetZaxis()->GetTitle()),
			NPTBINS, XBINS, yaxis->GetNbins(), yaxis->GetBinLowEdge(1), yaxis->GetBinLowEdge(h->GetYaxis()->GetNbins()+1));
		for (int j=1; j<=yaxis->GetNbins();j++)	{ 
		for (int i=1; i<=xaxis->GetNbins();i++)	{ 
			htmp->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h->GetBinContent(i,j)); 
		}	}	
			
		h = (TH2D*)htmp->Clone(h->GetName());
		delete htmp;
			
	}
	return h;

}

