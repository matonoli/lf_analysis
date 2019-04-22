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

	if (mFlagHist) return 0;

	TH1::SetDefaultSumw2();
	CreateHistograms();

	mList = (TList*)mHandler->directory()->GetList();

	bugR = 0;
	bugPt = 0;

	mTS = new TransverseSpherocity();
	mTS->SetMinMulti(10);

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
	if (!SelectEvent(event)) return 0;
	hEventMonitor->Fill(2);

	// BUG HOTFIX FOR AURORATREES
	MyV0 bugfix;
	if (mHandler->v0(0)) { bugfix = MyV0(mHandler->v0(0));
	if (TMath::Abs(bugfix.GetRadius()-bugR) < 0.0001
		&& TMath::Abs(bugfix.GetPt()-bugPt) < 0.0001) return 0;
	bugR = bugfix.GetRadius(); bugPt = bugfix.GetPt(); }
	hEventMonitor->Fill(3);

	// EVENT CLASSIFICATION
	enum { multMB, V0M, NCharged };
	enum { sphMB, Jetty, Iso };
	enum { D, RC, MC };
	Int_t isEventCentral = 0;
	hEventV0MCentrality->Fill(event.GetV0MCentrality());
	hEventRefMult->Fill(event.GetRefMult());
	hEventV0MCentvRefMult->Fill(event.GetRefMult(),event.GetV0MCentrality());
	isEventCentral += IsCentral(event,V0M);
	//isEventCentral += IsCentral(event,NCharged);

	mTS->Reset();
	Double_t eventTS = -99.;

	// TRACK LOOP
	hEventMonitor->Fill(4);
	Int_t nTracks = mHandler->tracks()->GetEntriesFast();
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!mHandler->track(iTr)) continue;
		MyTrack t(mHandler->track(iTr));

		if (!SelectTrack(t)) continue;

		if (isEventCentral) mTS->AddTrack(t.GetPx(), t.GetPy());		
	}

	// EVENT SPHEROCITY CLASSIFICATION
	Bool_t isEventIso = 0;
	Bool_t isEventJetty = 0;
	if (isEventCentral) {
		eventTS = mTS->GetTransverseSpherocityTracks();
		isEventIso		= (eventTS > cuts::EV_SPH_ISO && eventTS < 1.) ;
		isEventJetty	= (eventTS < cuts::EV_SPH_JETTY && eventTS > 0.);
		hEventSpherocity->Fill(eventTS);	}

	// MC V0 ANALYSIS: PARTICLES LOOP
	hEventMonitor->Fill(5);
	std::vector<Int_t> PartLabels;
	std::vector<Int_t> PartIds; 
	if (mFlagMC) {
		Int_t nParticles = mHandler->particles()->GetEntriesFast();
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!mHandler->particle(iP)) continue;
			MyParticle p(mHandler->particle(iP));

			if (!SelectParticle(p)) continue;
			PartLabels.push_back(p.GetLabel());
			PartIds.push_back(iP);

			for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
				if (p.GetPdgCode() == PDG_IDS[iSp]){
					hV0Pt[iSp][MC][0][0]->Fill(p.GetPt());	}
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
					MCfound = true;
					break;	}
			}
			if (MCfound) {
				for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
					if (IsV0(v0,iSp,RC)) ProcessV0(v0,iSp,RC,multMB,sphMB);		}
			}
		}

		for (int iMu = 0; iMu < isEventCentral+1; ++iMu) {
		for (int iSph = 0; iSph < isEventJetty+isEventIso+1; ++iSph) {
		for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
			if (IsV0(v0,iSp,D)) ProcessV0(v0,iSp,D,iMu,iSph+iSph*isEventIso);
		}	}	}

	}

	hEventMonitor->Fill(7);
	return 0;	
}

Bool_t MyAnalysisV0::SelectEvent(MyEvent &ev) {

	if (!ev.IsGoodAliEvent())			return false;
	//if(bV0s->GetEntriesFast() < 1)		return false;
	//if(bTracks->GetEntriesFast() < 1)		return false;
}


Bool_t MyAnalysisV0::IsCentral(MyEvent &ev, Int_t Mu) {

	switch (Mu) {
		default: 
			break;
		case 1: 
			if (ev.GetV0MCentrality() < 10. 
				&& ev.GetRefMult() > 9.9)		return true;
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

Bool_t MyAnalysisV0::IsV0(MyV0 &v0, Int_t Sp, Int_t Type) {
	
	if (Type==1) if (v0.GetMCPdgCode() != PDG_IDS[Sp]) return false;

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
			if (trP.GetNSigmaProtonTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			if (trP.GetNSigmaProtonTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			if (trN.GetNSigmaPionTPC() < cuts::K0S_D_NSIGTPC[0])	return false;
			if (trN.GetNSigmaPionTPC() > cuts::K0S_D_NSIGTPC[1])	return false;
			break;
		case 3 	: // Lbar
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

	// TRACK HISTOGRAMS

	// MC PARTICLE HISTOGRAMS

	// V0 HISTOGRAMS
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < NTYPE; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)			{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
				
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

}

Int_t MyAnalysisV0::Finish() {
	
	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	return 0;	
}

void MyAnalysisV0::DoEfficiency() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0Efficiency[iSp] = (TH1D*)hV0Pt[iSp][1][0][0]->Clone(Form("hV0Efficiency_%s",SPECIES[iSp]));

		hV0Efficiency[iSp]->SetTitle("; V0 pT (GeV/#it{c}); Efficiency");
		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Efficiency[iSp]->Divide(hV0Pt[iSp][2][0][0]);
		hV0Efficiency[iSp]->Write();
	}

}

/*void MyAnalysisV0::MakeFinalFigures() {

	mFilePlots = new TFile("plots_"+mOutName,"RECREATE");
	mFilePlots->cd();

	enum { LEFT, RIGHT };

	// EVENT INFO

	// SPHEROCITY
	Double_t quantileValues[4] = {0.0, 0.2, 0.8, 1.0};
	Double_t quantileCuts[4];
	hEventSpherocity->GetXaxis()->SetRangeUser(0.0,1.0);
	//hEventSpherocity->ClearUnderflowAndOverflow();
	hEventSpherocity->GetQuantiles(4,quantileCuts,quantileValues);


	TCanvas* cSpherocity = new TCanvas("cSpherocity","",1000,800);
	mHandler->MakeNiceHistogram(hEventSpherocity,kBlack);
	hEventSpherocity->Draw();
	cSpherocity->Update();
	mHandler->DrawCut(quantileCuts[1],LEFT,cSpherocity);
	mHandler->DrawCut(quantileCuts[2],RIGHT,cSpherocity);
	TLegend *leg1 = new TLegend(0.075,0.7,0.5,0.88);
	mHandler->MakeNiceLegend(leg1,0.050,1);
	leg1->AddEntry((TObject*)0,Form("Jetty: <%4.3f", quantileCuts[1]),"");
	leg1->AddEntry((TObject*)0,Form("Iso: >%4.3f", quantileCuts[2]),"");
	leg1->Draw();
	cSpherocity->Write();
	cSpherocity->SaveAs("plots/spherocity.png");

	// PT SPECTRA
	TCanvas* cPt[4];
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		for (int iSph = 0; iSph < NSPHERO; ++iSph)	{

			mHandler->MakeNiceHistogram(hV0PtFit[iSp][0][0][0],COLOURS[0]);
			mHandler->MakeNiceHistogram(hV0PtFit[iSp][0][1][iSph],COLOURS[1+iSph]);
		}

		cPt[iSp] = new TCanvas(Form("cPt_%s",SPECIES[iSp]),"",1000,800);
		cPt[iSp]->SetLogy(1);
		hV0PtFit[iSp][0][0][0]->GetYaxis()->SetRangeUser(0.1,10.*hV0PtFit[iSp][0][0][0]->GetMaximum());
		hV0PtFit[iSp][0][0][0]->Draw();
		cPt[iSp]->Update();
		hV0PtFit[iSp][0][1][0]->Draw("same");
		hV0PtFit[iSp][0][1][1]->Draw("same");
		hV0PtFit[iSp][0][1][2]->Draw("same");

		TLegend* legPt = new TLegend(0.55,0.55,0.85,0.85);
		mHandler->MakeNiceLegend(legPt,0.04,1);
		legPt->AddEntry((TObject*)0,Form("%s   |#eta| < 0.8", SPECNAMES[iSp]),"");
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"","");
		legPt->AddEntry(hV0PtFit[iSp][0][0][0],"MB","pl");
		legPt->AddEntry(hV0PtFit[iSp][0][1][0],"FHM","pl");
		legPt->AddEntry(hV0PtFit[iSp][0][1][1],"FHM Jetty","pl");
		legPt->AddEntry(hV0PtFit[iSp][0][1][2],"FHM Iso","pl");
		legPt->Draw();

		cPt[iSp]->Write();
		cPt[iSp]->SaveAs(Form("plots/pt_%s.png",SPECIES[iSp]));
	}

	// B/M RATIO
	TH1D* hBtoM[NMULTI][NSPHERO];
	for (int iMu = 0; iMu < 2; ++iMu) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
	
		hBtoM[iMu][iSph] = (TH1D*)hV0PtFit[2][0][iMu][iSph]->Clone(Form("hBtoM_%s_%s",MULTI[iMu],SPHERO[iSph]));
		hBtoM[iMu][iSph]->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / 2K^{0}_{s}");
		hBtoM[iMu][iSph]->Add(hV0PtFit[3][0][iMu][iSph]);
		hBtoM[iMu][iSph]->Divide(hBtoM[iMu][iSph],hV0PtFit[1][0][iMu][iSph],1.,2.,"");

	}	}
	TCanvas* cBtoM[2];
	TLegend* legBtoM[2];
	cBtoM[0] = new TCanvas("cBtoM_Mult","",1000,800);
	legBtoM[0] = new TLegend(0.55,0.55,0.85,0.85);
	mHandler->MakeNiceLegend(legBtoM[0],0.04,1);
	mHandler->MakeNiceHistogram(hBtoM[0][0],COLOURS[0]);
	mHandler->MakeNiceHistogram(hBtoM[1][0],COLOURS[1]);
	hBtoM[0][0]->GetYaxis()->SetRangeUser(0.,1.0);
	hBtoM[0][0]->Draw();
	cBtoM[0]->Update();
	hBtoM[1][0]->Draw("same");

	legBtoM[0]->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
	legBtoM[0]->AddEntry((TObject*)0,"","");
	legBtoM[0]->AddEntry(hBtoM[0][0],"MB","pl");
	legBtoM[0]->AddEntry(hBtoM[1][0],"FHM","pl");
	legBtoM[0]->Draw();

	cBtoM[0]->Write();
	cBtoM[0]->SaveAs("plots/btom_mu.png");

	cBtoM[1] = new TCanvas("cBtoM_Sph","",1000,800);
	legBtoM[1] = new TLegend(0.55,0.55,0.85,0.85);
	mHandler->MakeNiceLegend(legBtoM[1],0.04,1);
	mHandler->MakeNiceHistogram(hBtoM[1][1],COLOURS[2]);
	mHandler->MakeNiceHistogram(hBtoM[1][2],COLOURS[3]);
	hBtoM[1][1]->GetYaxis()->SetRangeUser(0.,1.0);
	hBtoM[1][1]->Draw();
	cBtoM[1]->Update();
	hBtoM[1][2]->Draw("same");

	legBtoM[1]->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
	legBtoM[1]->AddEntry((TObject*)0,"","");
	legBtoM[1]->AddEntry(hBtoM[1][1],"FHM Jetty","pl");
	legBtoM[1]->AddEntry(hBtoM[1][2],"FHM Iso","pl");
	legBtoM[1]->Draw();

	cBtoM[1]->Write();
	cBtoM[1]->SaveAs("plots/btom_sph.png");


}*/
