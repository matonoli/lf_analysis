#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>

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

//#include <AliAnalysisPIDV0.h>

using namespace V0consts;
using namespace RooFit;

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
	//printf("Looping in analysis %i \n", iEv);

	// EVENT INFO HISTOGRAMS
	hEventMonitor->Fill(0);
	if (!mEvent) return 1;
	MyEvent event(mEvent);
	hEventMonitor->Fill(1);

	// EVENT SELECTION AND CLASSIFICATION
	if (!SelectEvent(event)) return 0;
	hEventMonitor->Fill(2);
	enum { multMB, V0M, NCharged };
	enum { sphMB, Jetty, Iso };
	Int_t isEventCentral = 0;
	Int_t isEventSphero = 0;
	isEventCentral += IsCentral(event,V0M);
	isEventCentral += IsCentral(event,NCharged);

	// TRACK LOOP
	hEventMonitor->Fill(3);
	Int_t nTracks = bTracks->GetEntriesFast();
	for (int iTr = 0; iTr < nTracks; ++iTr)		{
		if (!(AliAnalysisPIDTrack*)bTracks->At(iTr)) continue;
		MyTrack t((AliAnalysisPIDTrack*)bTracks->At(iTr));
		
	}

	// MC V0 ANALYSIS: PARTICLES LOOP
	hEventMonitor->Fill(4);
	std::vector<Int_t> Particles; // Pairs: labels at even, id's at odds
	if (mFlagMC) {
		Int_t nParticles = bParticles->GetEntriesFast();
		for (int iP = 0; iP < nParticles; ++iP)		{
			
			if (!(AliAnalysisPIDParticle*)bParticles->At(iP)) continue;
			MyParticle p((AliAnalysisPIDParticle*)bParticles->At(iP));

			if (!SelectParticle(p)) continue;
			Particles.push_back(p.GetLabel());
			Particles.push_back(p.GetPdgCode());
		}

	}

	// V0 DATA ANALYSIS: V0 CANDIDATES LOOP
	hEventMonitor->Fill(5);
	Int_t nV0s = bV0s->GetEntriesFast();
	for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
		
		hV0Monitor->Fill(0);
		if (!(AliAnalysisPIDV0*)bV0s->At(iV0)) continue;
		MyV0 v0((AliAnalysisPIDV0*)bV0s->At(iV0));

		//printf("vo ap is %f and %f \n", *(v0.CalculateAP()+0), *(v0.CalculateAP()+1));

		for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
			if (IsV0(v0,iSp)) ProcessV0(v0,iSp,multMB,sphMB);		}

	}

	hEventMonitor->Fill(6);
	return 0;	
}

Bool_t MyAnalysisV0::SelectEvent(MyEvent &ev) {

	if (!ev.IsGoodAliEvent())			return false;
}


Bool_t MyAnalysisV0::IsCentral(MyEvent &ev, Int_t Mu) {

	return false;
}

Bool_t MyAnalysisV0::ProcessV0(MyV0 &v0, Int_t Sp, Int_t Mu, Int_t Sph) {
	
	//printf("v0 pt is %f, args are %i , %i , %i \n", v0.GetPt(), Sp, Mu, Sph);
	hV0Pt[Sp][Mu][Sph]->Fill(v0.GetPt());
	hV0Eta[Sp][Mu][Sph]->Fill(v0.GetEta());
	Double_t v0mass[] = {0., v0.GetIMK0s(), v0.GetIML(), v0.GetIMLbar()};
	hV0IMvPt[Sp][Mu][Sph]->Fill(v0.GetPt(),v0mass[Sp]);

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

	//if (!Sp) return true;
	MyTrack trP(v0.GetPosTrack()); 
	MyTrack trN(v0.GetNegTrack());
	if (!SelectV0Daughter(trP)) return false;
	if (!SelectV0Daughter(trN)) return false;

	switch (Sp) {
		default : 
			break;
		case 1 	: // K0s
			if (*(v0.CalculateAP()+1) < cuts::K0S_AP*fabs(*(v0.CalculateAP()+0))) return false;
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
	for (int iSp = 0; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
				
		hV0Pt[iSp][iMu][iSph]			= new TH1D(Form("hV0Pt_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); Entries",								400, 0, 20);
		hV0Eta[iSp][iMu][iSph]			= new TH1D(Form("hV0Eta_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),
			";V0 #eta; Entries", 										200, -1., 1.);
		hV0IMvPt[iSp][iMu][iSph]		= new TH2D(Form("hV0IMvPt_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, 1000, -0.1, 0.1);

		hV0PtFit[iSp][iMu][iSph]		= new TH1D(Form("hV0PtFit_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		

	} } }
		

}

Int_t MyAnalysisV0::Finish() {
	printf("Finishing analysis %s \n",this->GetName());

	// EXTRACT YIELDS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iMu = 0; iMu < 1; ++iMu)		{
	for (int iSph = 0; iSph < 1; ++iSph)	{
		
		Double_t* yield = 0;
		for (int iBin = 1; iBin < NPTBINS+1; ++iBin)	{
			yield = ExtractYieldFit(hV0IMvPt[iSp][iMu][iSph]->ProjectionY(Form("Mass bin %i", iBin),iBin,iBin));
			hV0PtFit[iSp][iMu][iSph]->SetBinContent(iBin,*(yield+0));
			hV0PtFit[iSp][iMu][iSph]->SetBinError(iBin,*(yield+1));
		}
		hV0PtFit[iSp][iMu][iSph]->Scale(1,"width");
		hV0PtFit[iSp][iMu][iSph]->SetMarkerColor(1);
		hV0PtFit[iSp][iMu][iSph]->SetMarkerStyle(20);
		hV0PtFit[iSp][iMu][iSph]->SetLineColor(2);
	} } }	

	return 0;	
}

Double_t* MyAnalysisV0::ExtractYieldFit(TH1D* hist) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	hist->Rebin(8);
	RooDataHist DT_set("DT_set","DT_hist",MassDT,Import(*hist)); 

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.001,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e05);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.1);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e05);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-200,200);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e05);
		
	RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
	RooFitResult* fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-999),Verbose(false));
		
	RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));
	//printf("Errors are %f and %f, total is %f or %f wrt to %f \n", nGaus1.getError(), nGaus2.getError(), nGaus1.getError()+nGaus2.getError(),sqrt(nGaus1.getError()*nGaus1.getError()+nGaus2.getError()*nGaus2.getError()),nGaus.getPropagatedError(*fR));
	
	/*TCanvas* can1 = new TCanvas();
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	fTotal.plotOn(plot1,LineWidth(1),LineColor(kRed));
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.35*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();*/


	/*cFits[canCounter%6]->cd(1+canCounter/6);
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	fTotal.plotOn(plot1,LineWidth(1),LineColor(kRed));
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.35*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();
	TLegend *leg1 = new TLegend(0.075,0.7,0.5,0.88);
	myLegendSetUp(leg1,0.065,1);
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",xBins[canCounter/6],xBins[1+canCounter/6])," ");
	leg1->AddEntry((TObject*)0,cNames[canCounter%6]+Form(" , #chi^{2}/ndf = %4.2f",plot1->chiSquare())," ");
	leg1->Draw();*/
	
	val[0] = nGaus.getVal();
	//printf("STATUS: int from fit is %f \n", val[0]);
	val[1] = nGaus.getPropagatedError(*fR);
	//canCounter++;
	
	return val;
}