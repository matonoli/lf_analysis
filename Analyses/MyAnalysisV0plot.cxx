#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TList.h>
#include <TFile.h>
#include <TLegend.h>
#include <TNamed.h>
#include <THashList.h>
#include <TProfile.h>

#include "MyAnalysisV0plot.h"
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

//#include <AliAnalysisPIDV0.h>
using namespace V0consts;
using namespace RooFit;

ClassImp(MyAnalysisV0plot)

MyAnalysisV0plot::MyAnalysisV0plot() {

}

Int_t MyAnalysisV0plot::Init() {

	TString dfName(this->GetName());
	dfName = Form("%s_%i",dfName.Data(),mHandler->nAnalysis());
	mDirFile = new TDirectoryFile(dfName,dfName,"",mHandler->file());
	mDirFile->cd();

	printf("Initialising analysis %s  \n", 
		this->GetName());

	TH1::SetDefaultSumw2(1);
	CreateHistograms();
	


	mList = (TList*)mHandler->directory()->GetList();

	return 0;
}

Int_t MyAnalysisV0plot::Make(Int_t iEv) {
	return 0;
}


Bool_t MyAnalysisV0plot::BorrowHistograms() {

	hEventSpherocityV0M			= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hEventSpherocityV0M");
	hEventSpherocityNCharged	= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hEventSpherocityNCharged");

	hTrackDPhivNchTrans			= (TH2D*)mHandler->analysis(0)->dirFile()->Get("hTrackDPhivNchTrans");
	hV0DPhivNchTrans			= (TH2D*)mHandler->analysis(0)->dirFile()->Get("hV0DPhivNchTrans");

	hNchTrans		= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hNchTrans");
	hRt2			= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hRt2");
	hNchvLeadPt2	= (TH2D*)mHandler->analysis(0)->dirFile()->Get("hNchvLeadPt2");

	hNchTransRCvMC 		= (TH2D*)mHandler->analysis(0)->dirFile()->Get("hNchTransRCvMC");
	hLeadPtvNchTrans 	= (TH2D*)mHandler->analysis(0)->dirFile()->Get("hLeadPtvNchTrans");
	hLeadPtvNchTrans0 	= (TH2D*)mHandler->analysis(0)->dirFile()->Get("hLeadPtvNchTrans0");


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0Efficiency[iSp] = (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0Efficiency_%s",SPECIES[iSp]));
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
		hV0EfficiencyRt[iSp][iReg] = (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]));
	} }

	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
		hProtonNchTransvPt[iReg] 	= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hProtonNchTransvPt_%s",REGIONS[iReg]));
		hPionNchTransvPt[iReg] 		= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hPionNchTransvPt_%s",REGIONS[iReg]));
		hLambdaNchTransvPt[iReg] 	= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hLambdaNchTransvPt_%s",REGIONS[iReg]));
		hK0sNchTransvPt[iReg] 		= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hK0sNchTransvPt_%s",REGIONS[iReg]));

	}

	Int_t nType = (mHandler->GetFlagMC()) ? NTYPE : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0PtFitCorr[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtFitCorr_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0PtFit[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0Pt[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));

		//cout << " " << Form("hV0PtFitCorr_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]) << endl;
		//cout << "h is " << hV0PtFitCorr[iSp][iType][iMu][iSph] << endl;
		//TH1D* htmp = (TH1D*)hV0PtFitCorr[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName(),XBINS2);
		//cout << "htmp is " << htmp << " name " << htmp->GetName() << endl;
		//hV0PtFitCorr[iSp][iType][iMu][iSph] = htmp;
		//cout << "h again is " << hV0PtFitCorr[iSp][iType][iMu][iSph] << endl;
		//TString tmpName = hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName();
		//
		//hV0PtFitCorr[iSp][iType][iMu][iSph] = (TH1D*)htmp->Clone(tmpName.Data());
		//if (hV0PtFitCorr[iSp][iType][iMu][iSph])
		//hV0PtFitCorr[iSp][iType][iMu][iSph] = (TH1D*)hV0PtFitCorr[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName(),XBINS2);
		//if (hV0PtFit[iSp][iType][iMu][iSph])
		//hV0PtFit[iSp][iType][iMu][iSph]		= (TH1D*)hV0PtFit[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,hV0PtFit[iSp][iType][iMu][iSph]->GetName(),XBINS2);
		
		if (hV0Pt[iSp][iType][iMu][iSph] && mHandler->IsRebinPt()) {
		hV0Pt[iSp][iType][iMu][iSph]	= (TH1D*)hV0Pt[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,hV0Pt[iSp][iType][iMu][iSph]->GetName(),XBINS2);	}
		//hV0PtFit[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,"",XBINS2);
		//hV0Pt[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,"",XBINS2);

	} } } }

	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hTrackPt[iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));

		if (hTrackPt[iType][iMu][iSph] && mHandler->IsRebinPt())
		hTrackPt[iType][iMu][iSph]		= (TH1D*)hTrackPt[iType][iMu][iSph]->Rebin(NPTBINS2,hTrackPt[iType][iMu][iSph]->GetName(),XBINS2); 

	} } }


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{

		hV0RtFitCorr[iSp][iType][iReg][iPtBin]	= (TH1D*)mHandler->analysis(2)->dirFile()->Get(Form("hV0RtFitCorr_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iPtBin));
		hV0RtFit[iSp][iType][iReg][iPtBin]		= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hV0RtFit_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iPtBin));
	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
			
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin] 
			= (TH1D*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtRtFitCorr_%s_%s_%s_%1.1f-%1.1f",SPECIES[iSp],TYPE[iType],REGIONS[iReg],RTBINS0[iRtBin],RTBINS0[iRtBin+1]) );


	} } } }	


}

Bool_t MyAnalysisV0plot::CreateHistograms() {
	

}


Bool_t MyAnalysisV0plot::CloneHistograms() {

	//Int_t type = (mHandler->GetFlagMC()) ? 2 : 0;

	for (int iMu = 0; iMu < NMULTI; ++iMu) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
		
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hBtoM[iMu][iSph] =			(TH1D*)hV0PtFitCorr[2][0][iMu][iSph]->Clone(Form("hBtoM_%s_%s",MULTI[iMu],SPHERO[iSph]));
		//hV0toNchDR[0][iMu][iSph] =	(TH1D*)hV0PtFit[1][0][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s",SPECIES[1],MULTI[iMu],SPHERO[iSph]));
		//hV0toNchDR[1][iMu][iSph] =	(TH1D*)hV0PtFit[2][0][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s",SPECIES[2],MULTI[iMu],SPHERO[iSph]));
		
	}	}

	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
		hBtoMRt[iReg][iRtBin]	= (TH1D*)hV0PtRtFitCorr[2][0][iReg][iRtBin]->Clone(Form("hBtoMRt_%s_%1.1f-%1.1f",REGIONS[iReg],RTBINS0[iRtBin],RTBINS0[iRtBin+1]));
	}	}

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iMu = 0; iMu < NMULTI; ++iMu) {
	for (int iType = 0; iType < nType; ++iType) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 

		hV0toNchDR[0][iType][iMu][iSph] =	(TH1D*)hV0PtFit[1][iType][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s_%s",SPECIES[1],TYPE[iType],MULTI[iMu],SPHERO[iSph]));	//k0s
		hV0toNchDR[1][iType][iMu][iSph] =	(TH1D*)hV0PtFit[2][iType][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s_%s",SPECIES[2],TYPE[iType],MULTI[iMu],SPHERO[iSph]));	//l+lbar
	}	}	}

	// if flagMC show particle level ratios too
	if (mHandler->GetFlagMC()) {

		for (int iMu = 0; iMu < NMULTI; ++iMu) {
		for (int iSph = 0; iSph < NSPHERO; ++iSph) {
			if (iMu > 2 && (iSph < 3 && iSph)) continue;
			if (iMu < 3 && iSph > 2) continue; 

			hV0toNchDR[0][2][iMu][iSph] =	(TH1D*)hV0Pt[1][2][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s_%s",SPECIES[1],TYPE[2],MULTI[iMu],SPHERO[iSph]));
			hV0toNchDR[1][2][iMu][iSph] =	(TH1D*)hV0Pt[2][2][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s_%s",SPECIES[2],TYPE[2],MULTI[iMu],SPHERO[iSph]));
		}	}			// so for mc we do particle level v0 only --> add also track lvl
	}

	//mDirFile->ls();

}

Int_t MyAnalysisV0plot::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();


	BorrowHistograms();
	CloneHistograms();
	MakeFinalFiguresSpherocity();
	//MakeFinalFiguresEvent();
	//MakeFinalFiguresRt();

	return 0;	
}

void MyAnalysisV0plot::MakeFinalFiguresEvent() {

		// MEAN NCH VS LEADPT
	{	
		TCanvas* cNchvPt = new TCanvas("cNchvPt","",1000,900);
		
		TProfile* hNchvPt = (TProfile*)hNchvLeadPt2->ProfileX();
		
		mHandler->MakeNiceHistogram((TH1D*)hNchvPt,kBlack);

		hNchvPt->GetYaxis()->SetTitle("<N_{ch}^{trans}>");// / (#Deltay #Delta#phi)");
		Double_t NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
		//hNchvPt->Scale(1./NormEta);
		Double_t NormPhi = 2./3.*TMath::Pi();
		//hNchvPt->Scale(1./NormPhi);

		hNchvPt->GetYaxis()->SetRangeUser(-0.01,2.*hNchvPt->GetMaximum());
		hNchvPt->Draw();

		cNchvPt->Update();
		mHandler->DrawCut(5.,2,cNchvPt);

		TLegend *leg1 = new TLegend(0.45,0.72,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		leg1->AddEntry((TObject*)0,"#bf{Transverse}","");
		
		cNchvPt->Update();
		leg1->Draw();

		cNchvPt->Write();
		cNchvPt->SaveAs("plots/rtNchvPt.png");
	}

	// RT DISTRIBUTION
	{
		Double_t rt_den = hNchTrans->GetMean();
		TCanvas* cRtDistro = new TCanvas("cRtDistro","",1000,900);
		mHandler->MakeNiceHistogram(hRt2,kBlack);
		cRtDistro->SetLogy();

		hRt2->GetXaxis()->SetRangeUser(-0.01,6.01);
		hRt2->GetYaxis()->SetRangeUser(1.,30.*hRt2->GetMaximum());
		hRt2->GetXaxis()->SetTitle("R_{T}");
		hRt2->Draw("");
		cRtDistro->Update();
		TLegend *leg1 = new TLegend(0.45,0.60,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		leg1->AddEntry((TObject*)0,"5.0 < p_{T}^{lead} < 40.0 (GeV/#it{c})","");
		leg1->AddEntry(hRt2,"R_{T} = N_{ch}^{trans} / <N_{ch}^{trans}>","pl");
		leg1->AddEntry((TObject*)0,Form("<N_{ch}^{trans}> = %1.3f",rt_den),"");
		cRtDistro->Update();
		leg1->Draw();

		cRtDistro->Write();
		cRtDistro->SaveAs("plots/rtdistro.png");
	}

	// RT MC V RC
	/*{
		Double_t rt_den = hNchTrans->GetMean();
		Double_t rt_denMC = hNchTransMC->GetMean();
		TCanvas* cRtMCRC = new TCanvas("cRtMCRC","",1000,900);
		//mHandler->MakeNiceHistogram(hRt2,kBlack);
		cRtDistro->SetLogz();

		hRt2->GetXaxis()->SetRangeUser(-0.01,6.01);
		hRt2->GetYaxis()->SetRangeUser(1.,30.*hRt2->GetMaximum());
		hRt2->GetXaxis()->SetTitle("R_{T}");
		hRt2->Draw("");
		cRtDistro->Update();
		TLegend *leg1 = new TLegend(0.45,0.60,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		leg1->AddEntry((TObject*)0,"5.0 < p_{T}^{lead} < 40.0 (GeV/#it{c})","");
		leg1->AddEntry(hRt2,"R_{T} = N_{ch}^{trans} / <N_{ch}^{trans}>","pl");
		leg1->AddEntry((TObject*)0,Form("<N_{ch}^{trans}> = %1.3f",rt_den),"");
		cRtDistro->Update();
		leg1->Draw();

		cRtDistro->Write();
		cRtDistro->SaveAs("plots/rtdistro.png");
	}*/

	{
		Double_t rt_den = hNchTrans->GetMean();
		TCanvas* cRtDistro = new TCanvas("cRtDistro","",1000,900);
		mHandler->MakeNiceHistogram(hRt2,kBlack);
		cRtDistro->SetLogy();

		hRt2->GetXaxis()->SetRangeUser(-0.01,6.01);
		hRt2->GetYaxis()->SetRangeUser(1.,30.*hRt2->GetMaximum());
		hRt2->GetXaxis()->SetTitle("R_{T}");
		hRt2->Draw("");
		cRtDistro->Update();
		TLegend *leg1 = new TLegend(0.45,0.60,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		leg1->AddEntry((TObject*)0,"5.0 < p_{T}^{lead} < 40.0 (GeV/#it{c})","");
		leg1->AddEntry(hRt2,"R_{T} = N_{ch}^{trans} / <N_{ch}^{trans}>","pl");
		leg1->AddEntry((TObject*)0,Form("<N_{ch}^{trans}> = %1.3f",rt_den),"");
		cRtDistro->Update();
		leg1->Draw();

		cRtDistro->Write();
		cRtDistro->SaveAs("plots/rtdistro.png");
	}

	/*{
		TCanvas* cEffi[NSPECIES];
		for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
			cEffi[iSp] = new TCanvas(Form("cEffi_%s",SPECIES[iSp]),"",1000,900);
			
			mHandler->MakeNiceHistogram(hV0Efficiency[iSp],kBlack);
			hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(-0.01,0.8);
			hV0Efficiency[iSp]->Draw();
			for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
				mHandler->MakeNiceHistogram(hV0EfficiencyRt[iSp][iReg],COLOURS[iReg]);
				hV0EfficiencyRt[iSp][iReg]->Draw("same");
			}

			cEffi[iSp]->Update();
			TLegend *leg1 = new TLegend(0.45,0.60,0.85,0.85);
			mHandler->MakeNiceLegend(leg1,0.037,1);
			leg1->AddEntry((TObject*)0,Form("#bf{%s} pp #sqrt{s} = 13 TeV",SPECNAMES[iSp]),"");
			leg1->AddEntry((TObject*)0,"","");
			leg1->AddEntry(hV0Efficiency[iSp],"MB","pl");
			leg1->AddEntry(hV0EfficiencyRt[iSp][0],"R_{T} Trans","pl");
			leg1->AddEntry(hV0EfficiencyRt[iSp][1],"R_{T} Near","pl");
			leg1->AddEntry(hV0EfficiencyRt[iSp][2],"R_{T} Away","pl");
			
			cEffi[iSp]->Update();
			leg1->Draw();

			cEffi[iSp]->Write();
			cEffi[iSp]->SaveAs(Form("plots/cEffi_%s.png",SPECIES[iSp]));
		}
	}*/

	/*{
		TCanvas* cDphi = new TCanvas("cDphi","",1000,900);
		TH1D* hDphi0 = (TH1D*)hV0DPhivNchTrans->ProjectionY("hd0",1,4);
		TH1D* hDphi1 = (TH1D*)hV0DPhivNchTrans->ProjectionY("hd1",7,14);
		TH1D* hDphi2 = (TH1D*)hV0DPhivNchTrans->ProjectionY("hd2",16,21);
		mHandler->MakeNiceHistogram(hDphi0,COLOURS[0]);
		mHandler->MakeNiceHistogram(hDphi1,COLOURS[1]);
		mHandler->MakeNiceHistogram(hDphi2,COLOURS[2]);
		//cDphi->SetLogy();

		hDphi0->Rebin(8); hDphi0->Scale(1./hDphi0->Integral(0,-1));
		hDphi1->Rebin(8); hDphi1->Scale(1./hDphi1->Integral(0,-1));
		hDphi2->Rebin(8); hDphi2->Scale(1./hDphi2->Integral(0,-1));

		hDphi0->GetYaxis()->SetTitle("a.u.");
		hDphi0->Draw("");
		hDphi1->Draw("same");
		hDphi2->Draw("same");
		cDphi->Update();
		
		TLegend *leg1 = new TLegend(0.64,0.60,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry(hDphi0,"low R_{T}","pl");
		leg1->AddEntry(hDphi1,"average R_{T}","pl");
		leg1->AddEntry(hDphi2,"high R_{T}","pl");
		
		cDphi->Update();
		leg1->Draw();

		cDphi->Write();
		cDphi->SaveAs("plots/rtdphi.png");
	}

	{
		TCanvas* cPpi = new TCanvas("cPpi","",1000,900);
		Double_t rt_den = hNchTrans->GetMean();
		TLegend *leg1 = new TLegend(0.13,0.50,0.50,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV  Pythia","");
		leg1->AddEntry((TObject*)0,"5.0 < p_{T}^{lead} < 40.0 (GeV/#it{c})","");
		leg1->AddEntry((TObject*)0,"","");
		{		
		TH1D* hPpi0 	= (TH1D*)hProtonNchTransvPt[0]->ProjectionX("hPpi0",1,-1);
		TH1D* h 		= (TH1D*)hPionNchTransvPt[0]->ProjectionX("h",1,-1);
		hPpi0->Divide(hPpi0,h,1.,1.,"");
		TH1D* hPpi1 	= (TH1D*)hProtonNchTransvPt[0]->ProjectionX("hPpi1",1,2+TMath::FloorNint(1.*rt_den));
		h 				= (TH1D*)hPionNchTransvPt[0]->ProjectionX("h",1,2+TMath::FloorNint(1.*rt_den));
		hPpi1->Divide(hPpi1,h,1.,1.,"");
		TH1D* hPpi2 	= (TH1D*)hProtonNchTransvPt[0]->ProjectionX("hPpi2",3+TMath::FloorNint(1.*rt_den),50);
		h 				= (TH1D*)hPionNchTransvPt[0]->ProjectionX("h",3+TMath::FloorNint(1.*rt_den),50);
		hPpi2->Divide(hPpi2,h,1.,1.,"");

		mHandler->MakeNiceHistogram(hPpi0,COLOURS[0]);
		mHandler->MakeNiceHistogram(hPpi1,COLOURS[0]);
		mHandler->MakeNiceHistogram(hPpi2,COLOURS[0]);

		hPpi0->SetLineWidth(3);
		hPpi1->SetLineWidth(2);hPpi1->SetLineStyle(2);
		hPpi2->SetLineWidth(2);hPpi2->SetLineStyle(9);

		leg1->AddEntry(hPpi0,"#bf{Trans} R_{T} inc.","l");
		leg1->AddEntry(hPpi1,"#bf{Trans} R_{T} < 1","l");
		leg1->AddEntry(hPpi2,"#bf{Trans} R_{T} > 1","l");

		hPpi0->GetXaxis()->SetRangeUser(0.,7.);
		hPpi0->GetYaxis()->SetRangeUser(0.,0.55);
		hPpi0->GetYaxis()->SetTitle("p / #pi");
		hPpi0->Draw("hist l");
		hPpi1->Draw("hist l same");
		hPpi2->Draw("hist l same");}

		{		
		TH1D* hPpi0 	= (TH1D*)hProtonNchTransvPt[1]->ProjectionX("hPpi01",1,-1);
		TH1D* h 		= (TH1D*)hPionNchTransvPt[1]->ProjectionX("h",1,-1);
		hPpi0->Divide(hPpi0,h,1.,1.,"");
		TH1D* hPpi1 	= (TH1D*)hProtonNchTransvPt[1]->ProjectionX("hPpi11",1,2+TMath::FloorNint(1.*rt_den));
		h 				= (TH1D*)hPionNchTransvPt[1]->ProjectionX("h",1,2+TMath::FloorNint(1.*rt_den));
		hPpi1->Divide(hPpi1,h,1.,1.,"");
		TH1D* hPpi2 	= (TH1D*)hProtonNchTransvPt[1]->ProjectionX("hPpi21",3+TMath::FloorNint(1.*rt_den),50);
		h 				= (TH1D*)hPionNchTransvPt[1]->ProjectionX("h",3+TMath::FloorNint(1.*rt_den),50);
		hPpi2->Divide(hPpi2,h,1.,1.,"");

		mHandler->MakeNiceHistogram(hPpi0,COLOURS[1]);
		mHandler->MakeNiceHistogram(hPpi1,COLOURS[1]);
		mHandler->MakeNiceHistogram(hPpi2,COLOURS[1]);

		hPpi0->SetLineWidth(3);
		hPpi1->SetLineWidth(2);hPpi1->SetLineStyle(2);
		hPpi2->SetLineWidth(2);hPpi2->SetLineStyle(9);


		leg1->AddEntry(hPpi0,"#bf{Near} R_{T} inc.","l");
		leg1->AddEntry(hPpi1,"#bf{Near} R_{T} < 1","l");
		leg1->AddEntry(hPpi2,"#bf{Near} R_{T} > 1","l");

		hPpi0->Draw("hist l same");
		hPpi1->Draw("hist l same");
		hPpi2->Draw("hist l same");}
		leg1->Draw();
		cPpi->Write();
		cPpi->SaveAs("plots/cppi.png");

		/*{		
		TH1D* hPpi0 	= (TH1D*)hProtonNchTransvPt[2]->ProjectionX("hPpi02",1,-1);
		TH1D* h 		= (TH1D*)hPionNchTransvPt[2]->ProjectionX("h",1,-1);
		hPpi0->Divide(hPpi0,h,1.,1.,"");
		TH1D* hPpi1 	= (TH1D*)hProtonNchTransvPt[2]->ProjectionX("hPpi12",1,2+TMath::FloorNint(rt_den));
		h 				= (TH1D*)hPionNchTransvPt[2]->ProjectionX("h",1,2+TMath::FloorNint(rt_den));
		hPpi1->Divide(hPpi1,h,1.,1.,"");
		TH1D* hPpi2 	= (TH1D*)hProtonNchTransvPt[2]->ProjectionX("hPpi22",3+TMath::FloorNint(rt_den),50);
		h 				= (TH1D*)hPionNchTransvPt[2]->ProjectionX("h",3+TMath::FloorNint(rt_den),50);
		hPpi2->Divide(hPpi2,h,1.,1.,"");

		mHandler->MakeNiceHistogram(hPpi0,COLOURS[2]);
		mHandler->MakeNiceHistogram(hPpi1,COLOURS[2]);
		mHandler->MakeNiceHistogram(hPpi2,COLOURS[2]);

		hPpi0->SetLineWidth(3);
		hPpi1->SetLineWidth(2);hPpi1->SetLineStyle(2);
		hPpi2->SetLineWidth(2);hPpi2->SetLineStyle(9);

		hPpi0->Draw("hist l same");
		hPpi1->Draw("hist l same");
		hPpi2->Draw("hist l same");}*/


		
		/*mHandler->MakeNiceHistogram(hDphi0,COLOURS[0]);
		mHandler->MakeNiceHistogram(hDphi1,COLOURS[1]);
		mHandler->MakeNiceHistogram(hDphi2,COLOURS[2]);
		//cDphi->SetLogy();

		hDphi0->Rebin(8); hDphi0->Scale(1./hDphi0->Integral(0,-1));
		hDphi1->Rebin(8); hDphi1->Scale(1./hDphi1->Integral(0,-1));
		hDphi2->Rebin(8); hDphi2->Scale(1./hDphi2->Integral(0,-1));

		hDphi0->GetYaxis()->SetTitle("a.u.");
		hDphi0->Draw("");
		hDphi1->Draw("same");
		hDphi2->Draw("same");
		cDphi->Update();
		
		TLegend *leg1 = new TLegend(0.64,0.60,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry(hDphi0,"low R_{T}","pl");
		leg1->AddEntry(hDphi1,"average R_{T}","pl");
		leg1->AddEntry(hDphi2,"high R_{T}","pl");
		
		cDphi->Update();
		leg1->Draw();

		cDphi->Write();
		cDphi->SaveAs("plots/rtdphi.png");*/
	

	if (mHandler->GetFlagMC() ) {
		TCanvas* cRtRCvMC = new TCanvas("cRtRCvMC","",900,900);

		hNchTransRCvMC->GetXaxis()->SetTitle("generated N_{ch}^{trans}");
		hNchTransRCvMC->GetYaxis()->SetTitle("reconstructed N_{ch}^{trans}");
		hNchTransRCvMC->GetYaxis()->SetTitleSize(0.045);
		hNchTransRCvMC->GetXaxis()->SetTitleSize(0.045);

		hNchTransRCvMC->SetStats(0);
		cRtRCvMC->SetLogz();
		hNchTransRCvMC->Draw("colz");

		cRtRCvMC->Write();
		cRtRCvMC->SaveAs("plots/rtrcmc.png");
	}

	{
		TCanvas* cPtLvRt = new TCanvas("cPtLvRt","",1000,900);

		hLeadPtvNchTrans0->GetYaxis()->SetTitleSize(0.045);
		hLeadPtvNchTrans0->GetXaxis()->SetTitleSize(0.045);

		TProfile* pfx1 = hLeadPtvNchTrans0->ProfileX();
		mHandler->MakeNiceHistogram((TH1D*)pfx1,kRed);

		hLeadPtvNchTrans0->SetStats(0);
		cPtLvRt->SetLogz();
		hLeadPtvNchTrans0->Draw("colz");
		pfx1->Draw("same");

		cPtLvRt->Write();
		cPtLvRt->SaveAs("plots/ptlvrt.png");
	}

	{
		TCanvas* cPtLvRt = new TCanvas("cPtLvRt1","",1000,900);

		hLeadPtvNchTrans->GetYaxis()->SetTitleSize(0.045);
		hLeadPtvNchTrans->GetXaxis()->SetTitleSize(0.045);

		TProfile* pfx1 = hLeadPtvNchTrans->ProfileX();
		mHandler->MakeNiceHistogram((TH1D*)pfx1,kRed);

		hLeadPtvNchTrans->SetStats(0);
		cPtLvRt->SetLogz();
		hLeadPtvNchTrans->Draw("colz");
		pfx1->Draw("same");

		cPtLvRt->Write();
		cPtLvRt->SaveAs("plots/ptlvrt1.png");
	}



	

}

void MyAnalysisV0plot::MakeFinalFiguresRt() {

	const char* isMC = (mHandler->GetFlagMC()) ? "MC^{blind}_{rec}" : "Data";

	// RT PT SPECTRA for different regions
	TCanvas* cPtRt[4][3];
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	
		mHandler->MakeNiceHistogram(hV0PtRtFitCorr[iSp][0][iReg][0],kBlack);
		hV0PtRtFitCorr[iSp][0][iReg][0]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");

		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			mHandler->MakeNiceHistogram(hV0PtRtFitCorr[iSp][0][iReg][iRtBin],COLOURS[iRtBin-2]);
			hV0PtRtFitCorr[iSp][0][iReg][iRtBin]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");		
		}

		cPtRt[iSp][iReg] = new TCanvas(Form("cPtRt2_%s_%s",SPECIES[iSp],REGIONS[iReg]),"",1000,900);
		cPtRt[iSp][iReg]->SetLogy(1);
		Double_t lowerRange = 0.1*hV0PtRtFitCorr[iSp][0][iReg][2]->GetBinContent(hV0PtRtFitCorr[iSp][0][iReg][2]->FindLastBinAbove());
		hV0PtRtFitCorr[iSp][0][iReg][0]->GetYaxis()->SetRangeUser(lowerRange,10.*hV0PtRtFitCorr[iSp][0][iReg][NRTBINS0-1]->GetMaximum());
		hV0PtRtFitCorr[iSp][0][iReg][0]->Draw();
		cPtRt[iSp][iReg]->Update();
		//hV0PtFitCorr[iSp][0][3][0]->Draw("same");
		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			hV0PtRtFitCorr[iSp][0][iReg][iRtBin]->Draw("same");
		}
			//hV0PtFitCorr[iSp][0][3][6]->Draw("same");
		//hV0PtFitCorr[iSp][0][3][7]->Draw("same");
		TLegend* legPt = new TLegend(0.65,0.45,0.85,0.88);
		mHandler->MakeNiceLegend(legPt,0.04,1);
			
		legPt->AddEntry((TObject*)0,Form("#bf{%s}   |#eta| < 0.8, %s", SPECNAMES[iSp], isMC),"");
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"","");
		legPt->AddEntry((TObject*)0,Form("#bf{%s}", REGIONS[iReg]),"");
		legPt->AddEntry((TObject*)0,"","");
		legPt->AddEntry(hV0PtRtFitCorr[iSp][0][iReg][0],"R_{T} inc.","pl");
		//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][0],Form("%s (any)",PLOTS_MULTI[3]),"pl");
		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			legPt->AddEntry(hV0PtRtFitCorr[iSp][0][iReg][iRtBin],Form("%1.1f < R_{T} < %1.1f",RTBINS0[iRtBin],RTBINS0[iRtBin+1]),"pl");
		}
			//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][6],Form("%s %s",PLOTS_MULTI[3],SPHERO[6]),"pl");
		//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][7],Form("%s %s",PLOTS_MULTI[3],SPHERO[7]),"pl");
		legPt->Draw();
		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			mHandler->MakeRatioPlot(hV0PtRtFitCorr[iSp][0][iReg][iRtBin],hV0PtRtFitCorr[iSp][0][iReg][0],
			cPtRt[iSp][iReg], -0.1,4.2);
		}
			
		cPtRt[iSp][iReg]->Write();
		cPtRt[iSp][iReg]->SaveAs(Form("plots/ptrt_%s_%s.png",SPECIES[iSp],REGIONS[iReg]));

	}	}

	// B/M Ratio
	TCanvas* cBtoM[NREGIONS+1];
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
		if (iRtBin==1) continue;

		hBtoMRt[iReg][iRtBin]->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / 2K^{0}_{s}");
		hBtoMRt[iReg][iRtBin]->Add(hV0PtRtFitCorr[3][0][iReg][iRtBin]);
		hBtoMRt[iReg][iRtBin]->Divide(hBtoMRt[iReg][iRtBin],hV0PtRtFitCorr[1][0][iReg][iRtBin],1.,2.,"");

	}	}
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
		
		cBtoM[iReg] = new TCanvas(Form("cBtoMRt_%s",REGIONS[iReg]),"",1000,800);
		
		hBtoMRt[iReg][0]->GetYaxis()->SetRangeUser(0.,1.0);
		mHandler->MakeNiceHistogram(hBtoMRt[iReg][0],kBlack);
		hBtoMRt[iReg][0]->Draw();
		cBtoM[iReg]->Update();

		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			
			mHandler->MakeNiceHistogram(hBtoMRt[iReg][iRtBin],COLOURS[iRtBin-2]);
			hBtoMRt[iReg][iRtBin]->Draw("same");

		}

		TLegend* legBtoM = new TLegend(0.60,0.49,0.85,0.85);
		mHandler->MakeNiceLegend(legBtoM,0.04,1);
		legBtoM->AddEntry((TObject*)0,Form("|#eta| < 0.8, %s", isMC),"");
		legBtoM->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legBtoM->AddEntry((TObject*)0,"","");
		legBtoM->AddEntry((TObject*)0,Form("#bf{%s}", REGIONS[iReg]),"");
		legBtoM->AddEntry((TObject*)0,"","");
		legBtoM->AddEntry(hBtoMRt[iReg][0],"R_{T} inc.","pl");
		//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][0],Form("%s (any)",PLOTS_MULTI[3]),"pl");
		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			legBtoM->AddEntry(hBtoMRt[iReg][iRtBin],Form("%1.1f < R_{T} < %1.1f",RTBINS0[iRtBin],RTBINS0[iRtBin+1]),"pl");
		}
		legBtoM->Draw();

		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			mHandler->MakeRatioPlot(hBtoMRt[iReg][iRtBin],hBtoMRt[iReg][0],
			cBtoM[iReg], 0.5,1.7);
		}

		cBtoM[iReg]->Write();
		cBtoM[iReg]->SaveAs(Form("plots/btomrt_%s.png",REGIONS[iReg]));

	}

	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
		if (iRtBin==1) continue;

		hBtoMRt[3][iRtBin] = (TH1D*)hV0PtRtFitCorr[2][0][1][iRtBin]->Clone(Form("hBtoMRt_Hard_%1.1f-%1.1f",RTBINS0[iRtBin],RTBINS0[iRtBin+1]));
		hBtoMRt[3][iRtBin]->Add(hV0PtRtFitCorr[3][0][1][iRtBin]);
		hBtoMRt[3][iRtBin]->Add(hV0PtRtFitCorr[2][0][2][iRtBin]);
		hBtoMRt[3][iRtBin]->Add(hV0PtRtFitCorr[3][0][2][iRtBin]);
		hBtoMRt[3][iRtBin]->Add(hV0PtRtFitCorr[2][0][0][iRtBin],-2.);
		hBtoMRt[3][iRtBin]->Add(hV0PtRtFitCorr[3][0][0][iRtBin],-2.);
		TH1D* hDen = (TH1D*)hV0PtRtFitCorr[1][0][1][iRtBin]->Clone("hDen");
		hDen->Add(hV0PtRtFitCorr[1][0][2][iRtBin]);
		hDen->Add(hV0PtRtFitCorr[1][0][0][iRtBin],-2.);
		hBtoMRt[3][iRtBin]->Divide(hBtoMRt[3][iRtBin],hDen,1.,2.,"");
		delete hDen;
	}
	cBtoM[3] = new TCanvas(Form("cBtoMRt_Hard"),"",1000,800);
		
		hBtoMRt[3][0]->GetYaxis()->SetRangeUser(0.,1.0);
		hBtoMRt[3][0]->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / 2K^{0}_{s}");
		mHandler->MakeNiceHistogram(hBtoMRt[3][0],kBlack);
		hBtoMRt[3][0]->Draw();
		cBtoM[3]->Update();

		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			
			mHandler->MakeNiceHistogram(hBtoMRt[3][iRtBin],COLOURS[iRtBin-2]);
			hBtoMRt[3][iRtBin]->Draw("same");

		}

		TLegend* legBtoM = new TLegend(0.60,0.49,0.85,0.85);
		mHandler->MakeNiceLegend(legBtoM,0.04,1);
		legBtoM->AddEntry((TObject*)0,Form("|#eta| < 0.8, %s", isMC),"");
		legBtoM->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legBtoM->AddEntry((TObject*)0,"","");
		legBtoM->AddEntry((TObject*)0,Form("#bf{Near+Away - 2xTrans}"),"");
		legBtoM->AddEntry((TObject*)0,"","");
		legBtoM->AddEntry(hBtoMRt[3][0],"R_{T} inc.","pl");
		//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][0],Form("%s (any)",PLOTS_MULTI[3]),"pl");
		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			legBtoM->AddEntry(hBtoMRt[3][iRtBin],Form("%1.1f < R_{T} < %1.1f",RTBINS0[iRtBin],RTBINS0[iRtBin+1]),"pl");
		}
		legBtoM->Draw();

		for (int iRtBin = 2; iRtBin < NRTBINS0; ++iRtBin)		{
			mHandler->MakeRatioPlot(hBtoMRt[3][iRtBin],hBtoMRt[3][0],
			cBtoM[3], -2.5,2.5);
		}

		cBtoM[3]->Write();
		cBtoM[3]->SaveAs(Form("plots/btomrt_Hard.png"));


	{TCanvas* cBtoMRegs = new TCanvas("cBtoMRegs","",1000,900);
		
		mHandler->MakeNiceHistogram(hBtoMRt[0][0],COLOURS[0]);
		mHandler->MakeNiceHistogram(hBtoMRt[1][0],COLOURS[1]);
		mHandler->MakeNiceHistogram(hBtoMRt[2][0],COLOURS[2]);		

		hBtoMRt[0][0]->Draw();
		hBtoMRt[1][0]->Draw("same");
		hBtoMRt[2][0]->Draw("same");
		
		TLegend *leg1 = new TLegend(0.45,0.55,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry(hBtoMRt[0][0],"#bf{Trans} R_{T} inc.","pl");
		leg1->AddEntry(hBtoMRt[1][0],"#bf{Near} R_{T} inc.","pl");
		leg1->AddEntry(hBtoMRt[2][0],"#bf{Away} R_{T} inc.","pl");
		//
		leg1->Draw();

		cBtoMRegs->Write();
		cBtoMRegs->SaveAs("plots/rtregs.png");
	}


	// SELF-NORMALISED YIELDS VS RT

	TCanvas* cYieldsRt[4];
	TF1* funcLin = new TF1("funcLin","x",-1.,10.);
	funcLin->SetLineColor(kBlack);
	TH1D* dummy1 = new TH1D(); mHandler->MakeNiceHistogram(dummy1,kBlack);
	TH1D* dummy2 = new TH1D(); mHandler->MakeNiceHistogram(dummy2,kBlack);
	dummy2->SetMarkerStyle(24);
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{

		cYieldsRt[iSp] = new TCanvas(Form("cYieldsRt_%s",SPECIES[iSp]),"",1000,900);
		hV0RtFitCorr[iSp][0][0][0]->GetYaxis()->SetTitle("#V0 / <V0>");
		hV0RtFitCorr[iSp][0][0][0]->GetXaxis()->SetRangeUser(-0.01,5.1);
		hV0RtFitCorr[iSp][0][0][0]->GetYaxis()->SetRangeUser(-0.01,12.01);
		hV0RtFitCorr[iSp][0][0][0]->Draw();
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
		for (int iPtBin = 0; iPtBin < NRTPTBINS-1; ++iPtBin)		{
			mHandler->MakeNiceHistogram(hV0RtFitCorr[iSp][0][iReg][iPtBin],COLOURS[iReg]);
			if (iPtBin==1) hV0RtFitCorr[iSp][0][iReg][iPtBin]->SetMarkerStyle(24); 
			hV0RtFitCorr[iSp][0][iReg][iPtBin]->Draw("same");
		} }

		TLegend* legRtYield = new TLegend(0.15,0.45,0.38,0.88);
		mHandler->MakeNiceLegend(legRtYield,0.04,1);
			
		legRtYield->AddEntry((TObject*)0,Form("#bf{%s}  |#eta| < 0.8, %s", SPECNAMES[iSp], isMC),"");
		legRtYield->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legRtYield->AddEntry((TObject*)0,"","");
		legRtYield->AddEntry(hV0RtFitCorr[iSp][0][0][0]," #bf{Trans}","l");
		legRtYield->AddEntry(hV0RtFitCorr[iSp][0][1][0]," #bf{Near}","l");
		legRtYield->AddEntry(hV0RtFitCorr[iSp][0][2][0]," #bf{Away}","l");
		legRtYield->AddEntry(dummy1,Form("%1.1f < p_{T} < %1.1f (GeV/#it{c})",RT_PTRANGE[0][0],RT_PTRANGE[0][1]),"p");
		legRtYield->AddEntry(dummy2,Form("%1.1f < p_{T} < %1.1f (GeV/#it{c})",RT_PTRANGE[1][0],RT_PTRANGE[1][1]),"p");
		legRtYield->AddEntry(funcLin,"y=x","l");
		legRtYield->Draw();

		funcLin->Draw("same");

		cYieldsRt[iSp]->Write();
		cYieldsRt[iSp]->SaveAs(Form("plots/rtyield_%s.png",SPECIES[iSp]));
	}

	// R_T YIELDS
	{
		Double_t rt_den = hNchTrans->GetMean();
		Int_t nbins = hV0RtFitCorr[1][0][0][2]->GetNbinsX();
		Double_t rtbins[nbins+1];
		for (int iBin = 0; iBin < nbins+1; ++iBin)	{
			rtbins[iBin] = (double)hV0RtFit[1][0][0][2]->GetBinLowEdge(iBin+1)/rt_den;	}
		
		for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
			hV0RtFit[1][0][iReg][2]->Scale(1.,"width");
			mHandler->MakeNiceHistogram(hV0RtFit[1][0][iReg][2],COLOURS[iReg]);
			hV0RtFit[1][0][iReg][2]->SetBins(nbins,rtbins);
			hV0RtFit[1][0][iReg][2]->GetXaxis()->SetRangeUser(rtbins[0],5.1);
			hV0RtFit[1][0][iReg][2]->GetXaxis()->SetTitle("R_{T}");
			//hV0RtFit[1][0][iReg][2]->GetYaxis()->SetRangeUser(0.,10.1);
		}

		TCanvas* cRtRawYield = new TCanvas("cRtRawYield","",1000,900);
		//mHandler->MakeNiceHistogram(hRt2,kBlack);
		//cRtRawYield->SetLogy();

		/*hRt2->GetXaxis()->SetRangeUser(-0.01,6.01);
		hRt2->GetYaxis()->SetRangeUser(10.,20.*hRt2->GetMaximum());
		hRt2->GetXaxis()->SetTitle("R_{T}");*/
		//hRt2->Draw("");
		cRtRawYield->Update();
		hV0RtFit[1][0][1][2]->GetYaxis()->SetRangeUser(-0.01,1.2*hV0RtFit[1][0][1][2]->GetMaximum());
		hV0RtFit[1][0][1][2]->Draw("same");
		hV0RtFit[1][0][0][2]->Draw("same");
		hV0RtFit[1][0][2][2]->Draw("same");
		cRtRawYield->Update();
		mHandler->DrawCut(0.5,0,cRtRawYield);
		mHandler->DrawCut(1.,0,cRtRawYield);
		mHandler->DrawCut(1.5,0,cRtRawYield);
		cRtRawYield->Update();
		TLegend *leg1 = new TLegend(0.45,0.55,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry((TObject*)0,"5.0 < p_{T}^{lead} < 40.0 (GeV/#it{c})","");
		
		leg1->AddEntry((TObject*)0,"","");
		leg1->AddEntry(hV0RtFit[1][0][0][2],"K_{s}^{0} raw yield #bf{Trans}","pl");
		leg1->AddEntry(hV0RtFit[1][0][1][2],"K_{s}^{0} raw yield #bf{Near}","pl");
		leg1->AddEntry(hV0RtFit[1][0][2][2],"K_{s}^{0} raw yield #bf{Away}","pl");
		//
		cRtRawYield->Update();
		leg1->Draw();

		cRtRawYield->Write();
		cRtRawYield->SaveAs("plots/rtrawyields.png");
	}

}

void MyAnalysisV0plot::MakeFinalFiguresSpherocity() {

	//mFilePlots = new TFile("plots_"+mOutName,"RECREATE");
	//mFilePlots->cd();

	//mHandler->root()->SetBatch(kTRUE);

	enum { LEFT, RIGHT };
	const char* isMC = (mHandler->GetFlagMC()) ? "MC^{blind}_{rec}" : "Data";

	// EVENT INFO

	// SPHEROCITY
	Double_t quantileValues[4] = {0.0, 0.2, 0.8, 1.0};
	Double_t quantileCuts[4];
	hEventSpherocityV0M->GetXaxis()->SetRangeUser(0.0,1.0);
	//hEventSpherocity->ClearUnderflowAndOverflow();
	hEventSpherocityV0M->GetQuantiles(4,quantileCuts,quantileValues);

	TCanvas* cSpherocityV0M = new TCanvas("cSpherocityV0M","",1000,800);
	mHandler->MakeNiceHistogram(hEventSpherocityV0M,kBlack);
	hEventSpherocityV0M->Draw();
	cSpherocityV0M->Update();
	mHandler->DrawCut(quantileCuts[1],LEFT,cSpherocityV0M);
	mHandler->DrawCut(quantileCuts[2],RIGHT,cSpherocityV0M);
	TLegend *leg1 = new TLegend(0.075,0.7,0.5,0.88);
	mHandler->MakeNiceLegend(leg1,0.050,1);
	leg1->AddEntry((TObject*)0,"V0M 0-10%","");
	leg1->AddEntry((TObject*)0,Form("Jetty: <%4.3f", quantileCuts[1]),"");
	leg1->AddEntry((TObject*)0,Form("Iso: >%4.3f", quantileCuts[2]),"");
	leg1->Draw();
	cSpherocityV0M->Write();
	cSpherocityV0M->SaveAs("plots/spherocity_v0m.png");

	hEventSpherocityNCharged->GetXaxis()->SetRangeUser(0.0,1.0);
	//hEventSpherocity->ClearUnderflowAndOverflow();
	hEventSpherocityNCharged->GetQuantiles(4,quantileCuts,quantileValues);

	TCanvas* cSpherocityNCharged = new TCanvas("cSpherocityNCharged","",1000,800);
	mHandler->MakeNiceHistogram(hEventSpherocityNCharged,kBlack);
	hEventSpherocityNCharged->Draw();
	cSpherocityNCharged->Update();
	mHandler->DrawCut(quantileCuts[1],LEFT,cSpherocityNCharged);
	mHandler->DrawCut(quantileCuts[2],RIGHT,cSpherocityNCharged);
	TLegend *leg2 = new TLegend(0.075,0.7,0.5,0.88);
	mHandler->MakeNiceLegend(leg2,0.050,1);
	leg2->AddEntry((TObject*)0,"NCharged 0-10%","");
	leg2->AddEntry((TObject*)0,Form("Jetty: <%4.3f", quantileCuts[1]),"");
	leg2->AddEntry((TObject*)0,Form("Iso: >%4.3f", quantileCuts[2]),"");
	leg2->Draw();
	cSpherocityNCharged->Write();
	cSpherocityNCharged->SaveAs("plots/spherocity_NCharged.png");

	// PT SPECTRA
	TCanvas* cPt[4][2];
	for (Int_t iSp = 1; iSp < NSPECIES; ++iSp)	{				
	
		mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][0][0],COLOURS[0]);
		hV0PtFitCorr[iSp][0][0][0]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");

		for (Int_t iMu = 1; iMu < 3; ++iMu)	{
			
			for (Int_t iSph = 0; iSph < 3; ++iSph)	{

				//mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][iMu][iSph],COLOURS[iMu+iSph+(iSph>0)-(iMu>1&&iSph>0)]);
				mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][iMu][iSph],COLOURS[1+iSph]);
				hV0PtFitCorr[iSp][0][iMu][iSph]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");		

			}

			cPt[iSp][iMu-1] = new TCanvas(Form("cPt_%s_%s",SPECIES[iSp],MULTI[iMu]),"",1000,900);
			cPt[iSp][iMu-1]->SetLogy(1);
			Double_t lowerRange = 0.1*hV0PtFitCorr[iSp][0][0][0]->GetBinContent(hV0PtFitCorr[iSp][0][0][0]->FindLastBinAbove());
			if (iMu==1) hV0PtFitCorr[iSp][0][0][0]->GetYaxis()->SetRangeUser(lowerRange,10.*hV0PtFitCorr[iSp][0][0][0]->GetMaximum());
			hV0PtFitCorr[iSp][0][0][0]->Draw();
			cPt[iSp][iMu-1]->Update();
			hV0PtFitCorr[iSp][0][iMu][0]->Draw("same");
			hV0PtFitCorr[iSp][0][iMu][1]->Draw("same");
			hV0PtFitCorr[iSp][0][iMu][2]->Draw("same");

			TLegend* legPt = new TLegend(0.56,0.54,0.85,0.88);
			mHandler->MakeNiceLegend(legPt,0.04,1);
			
			legPt->AddEntry((TObject*)0,Form("%s   |#eta| < 0.8, %s", SPECNAMES[iSp], isMC),"");
			legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
			legPt->AddEntry((TObject*)0,"","");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][0][0],"MB","pl");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][iMu][0],PLOTS_MULTI[iMu],"pl");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][iMu][1],Form("%s %s",PLOTS_MULTI[iMu],SPHERO[1]),"pl");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][iMu][2],Form("%s %s",PLOTS_MULTI[iMu],SPHERO[2]),"pl");
			legPt->Draw();

			mHandler->MakeRatioPlot(hV0PtFitCorr[iSp][0][iMu][1],hV0PtFitCorr[iSp][0][iMu][0],
				cPt[iSp][iMu-1], 0.4,1.6);
			mHandler->MakeRatioPlot(hV0PtFitCorr[iSp][0][iMu][2],hV0PtFitCorr[iSp][0][iMu][0],
				cPt[iSp][iMu-1], 0.4,1.6);

			cPt[iSp][iMu-1]->Write();
			cPt[iSp][iMu-1]->SaveAs(Form("plots/pt_%s_%s.png",SPECIES[iSp],MULTI[iMu]));
		}
	}

	return;

	// RT PT SPECTRA
	TCanvas* cPtRt[4];
	for (int iSp = 1; iSp < NSPECIES; ++iSp)
	{
		mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][3][0],kRed);
		hV0PtFitCorr[iSp][0][3][0]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");

		for (int iRt = 0; iRt < 5; ++iRt)		{
			mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][3][3+iRt],COLOURS[5-iRt]);
			hV0PtFitCorr[iSp][0][3][3+iRt]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");		
		}

		cPtRt[iSp] = new TCanvas(Form("cPtRt_%s",SPECIES[iSp]),"",1000,900);
			cPtRt[iSp]->SetLogy(1);
			Double_t lowerRange = 0.1*hV0PtFitCorr[iSp][0][0][0]->GetBinContent(hV0PtFitCorr[iSp][0][0][0]->FindLastBinAbove());
			hV0PtFitCorr[iSp][0][0][0]->Draw();
			cPtRt[iSp]->Update();
			//hV0PtFitCorr[iSp][0][3][0]->Draw("same");
			hV0PtFitCorr[iSp][0][3][3]->Draw("same");
			hV0PtFitCorr[iSp][0][3][4]->Draw("same");
			hV0PtFitCorr[iSp][0][3][5]->Draw("same");
			//hV0PtFitCorr[iSp][0][3][6]->Draw("same");
			//hV0PtFitCorr[iSp][0][3][7]->Draw("same");

			TLegend* legPt = new TLegend(0.65,0.56,0.85,0.88);
			mHandler->MakeNiceLegend(legPt,0.04,1);
			
			legPt->AddEntry((TObject*)0,Form("%s   |#eta| < 0.8, %s", SPECNAMES[iSp], isMC),"");
			legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
			legPt->AddEntry((TObject*)0,"","");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][0][0],"MB","pl");
			//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][0],Form("%s (any)",PLOTS_MULTI[3]),"pl");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][3][3],Form("%s %s",PLOTS_MULTI[3],SPHERO[3]),"pl");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][3][4],Form("%s %s",PLOTS_MULTI[3],SPHERO[4]),"pl");
			legPt->AddEntry(hV0PtFitCorr[iSp][0][3][5],Form("%s %s",PLOTS_MULTI[3],SPHERO[5]),"pl");
			//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][6],Form("%s %s",PLOTS_MULTI[3],SPHERO[6]),"pl");
			//legPt->AddEntry(hV0PtFitCorr[iSp][0][3][7],Form("%s %s",PLOTS_MULTI[3],SPHERO[7]),"pl");
			legPt->Draw();

			mHandler->MakeRatioPlot(hV0PtFitCorr[iSp][0][3][3],hV0PtFitCorr[iSp][0][3][4],
				cPtRt[iSp], 0.1,3.2);
			mHandler->MakeRatioPlot(hV0PtFitCorr[iSp][0][3][4],hV0PtFitCorr[iSp][0][3][4],
				cPtRt[iSp], 0.1,3.2);
			mHandler->MakeRatioPlot(hV0PtFitCorr[iSp][0][3][5],hV0PtFitCorr[iSp][0][3][4],
				cPtRt[iSp], 0.1,3.2);

			cPtRt[iSp]->Write();
			cPtRt[iSp]->SaveAs(Form("plots/pt_%s_%s.png",SPECIES[iSp],MULTI[3]));

	}

	// B/M RATIO
	//TH1D* hBtoM[NMULTI][NSPHERO];
	for (int iMu = 0; iMu < NMULTI; ++iMu) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		
		hBtoM[iMu][iSph]->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / 2K^{0}_{s}");
		hBtoM[iMu][iSph]->Add(hV0PtFitCorr[3][0][iMu][iSph]);
		hBtoM[iMu][iSph]->Divide(hBtoM[iMu][iSph],hV0PtFitCorr[1][0][iMu][iSph],1.,2.,"");

	}	}
	
	TCanvas* cBtoM[4];
	TLegend* legBtoM[4];
	cBtoM[0] = new TCanvas("cBtoM_Mult","",1000,800);
	legBtoM[0] = new TLegend(0.55,0.55,0.85,0.85);
	mHandler->MakeNiceLegend(legBtoM[0],0.04,1);
	mHandler->MakeNiceHistogram(hBtoM[0][0],COLOURS[0]);
	mHandler->MakeNiceHistogram(hBtoM[1][0],COLOURS[1]);
	mHandler->MakeNiceHistogram(hBtoM[2][0],COLOURS[2]);
	hBtoM[0][0]->GetYaxis()->SetRangeUser(0.,1.0);
	hBtoM[0][0]->Draw();
	cBtoM[0]->Update();
	hBtoM[1][0]->Draw("same");
	hBtoM[2][0]->Draw("same");

	legBtoM[0]->AddEntry((TObject*)0,Form("pp #sqrt{s} = 13 TeV"),"");
	legBtoM[0]->AddEntry((TObject*)0,isMC,"");
	legBtoM[0]->AddEntry((TObject*)0,"","");
	legBtoM[0]->AddEntry(hBtoM[0][0],PLOTS_MULTI[0],"pl");
	legBtoM[0]->AddEntry(hBtoM[1][0],PLOTS_MULTI[1],"pl");
	legBtoM[0]->AddEntry(hBtoM[2][0],PLOTS_MULTI[2],"pl");
	legBtoM[0]->Draw();

	cBtoM[0]->Write();
	cBtoM[0]->SaveAs("plots/btom_mu.png");

	cBtoM[1] = new TCanvas("cBtoM_V0M_Sph","",1000,800);
	legBtoM[1] = new TLegend(0.55,0.55,0.85,0.85);
	mHandler->MakeNiceLegend(legBtoM[1],0.04,1);
	mHandler->MakeNiceHistogram(hBtoM[1][1],COLOURS[3]);
	mHandler->MakeNiceHistogram(hBtoM[1][2],COLOURS[4]);
	hBtoM[1][1]->GetYaxis()->SetRangeUser(0.,1.0);
	hBtoM[1][1]->Draw();
	cBtoM[1]->Update();
	hBtoM[1][2]->Draw("same");

	legBtoM[1]->AddEntry((TObject*)0,Form("pp #sqrt{s} = 13 TeV"),"");
	legBtoM[1]->AddEntry((TObject*)0,isMC,"");
	legBtoM[1]->AddEntry((TObject*)0,"","");
	legBtoM[1]->AddEntry(hBtoM[1][1],Form("%s %s",PLOTS_MULTI[1],SPHERO[1]),"pl");
	legBtoM[1]->AddEntry(hBtoM[1][2],Form("%s %s",PLOTS_MULTI[1],SPHERO[2]),"pl");
	legBtoM[1]->Draw();

	cBtoM[1]->Write();
	cBtoM[1]->SaveAs("plots/btom_v0m_sph.png");

	cBtoM[2] = new TCanvas("cBtoM_Nch_Sph","",1000,800);
	legBtoM[2] = new TLegend(0.55,0.55,0.85,0.85);
	mHandler->MakeNiceLegend(legBtoM[2],0.04,1);
	mHandler->MakeNiceHistogram(hBtoM[2][1],COLOURS[3]);
	mHandler->MakeNiceHistogram(hBtoM[2][2],COLOURS[4]);
	hBtoM[2][1]->GetYaxis()->SetRangeUser(0.,1.0);
	hBtoM[2][1]->Draw();
	cBtoM[2]->Update();
	hBtoM[2][2]->Draw("same");

	legBtoM[2]->AddEntry((TObject*)0,Form("pp #sqrt{s} = 13 TeV"),"");
	legBtoM[2]->AddEntry((TObject*)0,isMC,"");
	legBtoM[2]->AddEntry((TObject*)0,"","");
	legBtoM[2]->AddEntry(hBtoM[1][1],Form("%s %s",PLOTS_MULTI[2],SPHERO[1]),"pl");
	legBtoM[2]->AddEntry(hBtoM[1][2],Form("%s %s",PLOTS_MULTI[2],SPHERO[2]),"pl");
	legBtoM[2]->Draw();

	cBtoM[2]->Write();
	cBtoM[2]->SaveAs("plots/btom_nch_sph.png");

	cBtoM[3] = new TCanvas("cBtoM_Rt","",1000,800);
	legBtoM[3] = new TLegend(0.55,0.55,0.85,0.85);
	mHandler->MakeNiceLegend(legBtoM[3],0.04,1);
	mHandler->MakeNiceHistogram(hBtoM[3][3],COLOURS[3]);
	mHandler->MakeNiceHistogram(hBtoM[3][4],COLOURS[4]);
	mHandler->MakeNiceHistogram(hBtoM[3][5],COLOURS[5]);
	hBtoM[3][3]->GetYaxis()->SetRangeUser(0.,1.0);
	hBtoM[3][3]->Draw();
	cBtoM[3]->Update();
	hBtoM[3][4]->Draw("same");
	hBtoM[3][5]->Draw("same");

	legBtoM[3]->AddEntry((TObject*)0,Form("pp #sqrt{s} = 13 TeV"),"");
	legBtoM[3]->AddEntry((TObject*)0,isMC,"");
	legBtoM[3]->AddEntry((TObject*)0,"","");
	legBtoM[3]->AddEntry(hBtoM[3][3],Form("%s %s",PLOTS_MULTI[3],SPHERO[3]),"pl");
	legBtoM[3]->AddEntry(hBtoM[3][4],Form("%s %s",PLOTS_MULTI[3],SPHERO[4]),"pl");
	legBtoM[3]->AddEntry(hBtoM[3][5],Form("%s %s",PLOTS_MULTI[3],SPHERO[5]),"pl");
	legBtoM[3]->Draw();

	cBtoM[3]->Write();
	cBtoM[3]->SaveAs("plots/btom_rt.png");


	// DOUBLE RATIO CALC
	Int_t nType = (mHandler->GetFlagMC()) ? 3 : 1;
	for (int iSp2 = 0; iSp2 < 2; ++iSp2)		{
	for (int iType = 0; iType < nType; ++iType)	{
	for (int iSph = 0; iSph < 3; ++iSph)		{
			for (int iMu = 0; iMu < 3; ++iMu)		{

				
				if (iSp2==1) {
					if (iType!=2)	hV0toNchDR[iSp2][iType][iMu][iSph]->Add(hV0PtFit[2][iType][iMu][iSph]);
					if (iType==2)	hV0toNchDR[iSp2][iType][iMu][iSph]->Add(hV0Pt[2][iType][iMu][iSph]);	}	//adding Lbar
				
				
				if (iType==1) hV0toNchDR[iSp2][iType][iMu][iSph]->Divide(hTrackPt[1][iMu][iSph]);
				if (iType!=1) hV0toNchDR[iSp2][iType][iMu][iSph]->Divide(hTrackPt[iType][iMu][iSph]);			//dividing by track spectra of same class -> bug for particle level! (?)
				
				if (iType>0) hV0toNchDR[iSp2][iType][iMu][iSph]->SetLineWidth(2);
			}	
	}	}	}	// now we have ratios of v0 to pi(charged tracks) in mu and sph classes

	for (int iSp2 = 0; iSp2 < 2; ++iSp2)		{
	for (int iType = 0; iType < nType; ++iType)	{
	for (int iSph = 0; iSph < 3; ++iSph)		{
		hV0toNchDR[iSp2][iType][1][iSph]->Divide(hV0toNchDR[iSp2][iType][0][0]);	// dividing v0/track of v0m events by v0/track mb
		hV0toNchDR[iSp2][iType][2][iSph]->Divide(hV0toNchDR[iSp2][iType][0][0]);
	}	}	}	// iMu!=0 are now (v0/pi)_hm  / (v0/pi)_mb double ratios


	// DOUBLE RATIO PLOTTING
	TCanvas* cDR[2][2];
	TLegend* legDR[2][2];
	TString drawOpt = (mHandler->GetFlagMC()) ? "hist l" : "";
	for (int iSp2 = 0; iSp2 < 2; ++iSp2)	{
	for (int iMu = 1; iMu < 3; ++iMu)	{
		
		cDR[iSp2][iMu-1] = new TCanvas(Form("cDR_%s_%s",SPECIES[1+iSp2],MULTI[iMu]),"",1000,800);
		//cDR[iSp2][iMu-1]->SetLogy();
		hV0toNchDR[iSp2][0][iMu][0]->GetYaxis()->SetRangeUser(0.2,4);
		hV0toNchDR[iSp2][0][iMu][0]->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
		hV0toNchDR[iSp2][0][iMu][0]->GetYaxis()->SetTitle("ratio");
		mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][0][iMu][0],COLOURS[0]);
		mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][0][iMu][1],COLOURS[3]);
		mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][0][iMu][2],COLOURS[4]);
	
		hV0toNchDR[iSp2][0][iMu][0]->Draw();
		hV0toNchDR[iSp2][0][iMu][1]->Draw("same");
		hV0toNchDR[iSp2][0][iMu][2]->Draw("same");

		if (mHandler->GetFlagMC()) {
			mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][1][iMu][0],COLOURS[0]);
			mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][1][iMu][1],COLOURS[3]);
			mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][1][iMu][2],COLOURS[4]);
			mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][2][iMu][0],COLOURS[0]);
			mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][2][iMu][1],COLOURS[3]);
			mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][2][iMu][2],COLOURS[4]);
			hV0toNchDR[iSp2][1][iMu][0]->SetLineWidth(3);
			hV0toNchDR[iSp2][1][iMu][0]->SetLineStyle(2);
			hV0toNchDR[iSp2][1][iMu][1]->SetLineWidth(3);
			hV0toNchDR[iSp2][1][iMu][1]->SetLineStyle(2);
			hV0toNchDR[iSp2][1][iMu][2]->SetLineWidth(3);
			hV0toNchDR[iSp2][1][iMu][2]->SetLineStyle(2);

			hV0toNchDR[iSp2][1][iMu][0]->Draw("hist l same");
			hV0toNchDR[iSp2][1][iMu][1]->Draw("hist l same");
			hV0toNchDR[iSp2][1][iMu][2]->Draw("hist l same");
			hV0toNchDR[iSp2][2][iMu][0]->Draw("hist l same");
			hV0toNchDR[iSp2][2][iMu][1]->Draw("hist l same");
			hV0toNchDR[iSp2][2][iMu][2]->Draw("hist l same");
		}

		legDR[iSp2][iMu-1] = (!mHandler->GetFlagMC()) ? new TLegend(0.55,0.55,0.85,0.85) : new TLegend(0.55,0.49,0.85,0.85);
		if (!mHandler->GetFlagMC()) mHandler->MakeNiceLegend(legDR[iSp2][iMu-1],0.04,1);
		else 						mHandler->MakeNiceLegend(legDR[iSp2][iMu-1],0.027,1);
		legDR[iSp2][iMu-1]->AddEntry((TObject*)0,Form("pp #sqrt{s} = 13 TeV"),"");
		legDR[iSp2][iMu-1]->AddEntry((TObject*)0,Form("%s  %s", SPECNAMES[1+iSp2], PLOTS_MULTI[iMu]),"");
		legDR[iSp2][iMu-1]->AddEntry((TObject*)0,"","");
		TString labelOpt = (mHandler->GetFlagMC()) ? ", RC V0, track S_{0}" : ""; 
		legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][0][iMu][0],"HM"+labelOpt,"pl");
		legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][0][iMu][1],"HM Jetty"+labelOpt,"pl");
		legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][0][iMu][2],"HM Iso"+labelOpt,"pl");
		if (mHandler->GetFlagMC()) {
			legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][1][iMu][0],"HM, RC V0, particle S_{0}","l");
			legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][1][iMu][1],"HM Jetty, RC V0, particle S_{0}","l");
			legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][1][iMu][2],"HM Iso, RC V0, particle S_{0}","l");
			legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][2][iMu][0],"HM, MC V0, particle S_{0}","l");
			legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][2][iMu][1],"HM Jetty, MC V0, particle S_{0}","l");
			legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][2][iMu][2],"HM Iso, MC V0, particle S_{0}","l");
		}
		legDR[iSp2][iMu-1]->Draw();

		mHandler->MakeZoomPlot(hV0toNchDR[iSp2][0][iMu][0],cDR[iSp2][iMu-1],0.,2.6,0.4,1.3);
		hV0toNchDR[iSp2][0][iMu][0]->DrawCopy();
		hV0toNchDR[iSp2][0][iMu][1]->DrawCopy("same");
		hV0toNchDR[iSp2][0][iMu][2]->DrawCopy("same");
		if (mHandler->GetFlagMC()) {
			hV0toNchDR[iSp2][1][iMu][0]->DrawCopy("hist l same");
			hV0toNchDR[iSp2][1][iMu][1]->DrawCopy("hist l same");
			hV0toNchDR[iSp2][1][iMu][2]->DrawCopy("hist l same");
			hV0toNchDR[iSp2][2][iMu][0]->DrawCopy("hist l same");
			hV0toNchDR[iSp2][2][iMu][1]->DrawCopy("hist l same");
			hV0toNchDR[iSp2][2][iMu][2]->DrawCopy("hist l same"); 
		}
		cDR[iSp2][iMu-1]->cd();
		hV0toNchDR[iSp2][0][iMu][0]->GetXaxis()->SetRangeUser(0.,14.);
		hV0toNchDR[iSp2][0][iMu][0]->GetYaxis()->SetRangeUser(0.2,4);

		cDR[iSp2][iMu-1]->Write();
		cDR[iSp2][iMu-1]->SaveAs(Form("plots/dr_%s_%s.png",SPECIES[1+iSp2],MULTI[iMu]));

	}	}


	// self-normalised yields vs rt
	/*TCanvas* cRtYields[4];
	TF1* funcLin = new TF1("funcLin","x",RTBINS[0],RTBINS[NRTBINS]);
	for (int iSp = 1; iSp < NSPECIES; ++iSp) {
		mHandler->MakeNiceHistogram(hV0RtFitCorr[iSp][0][0],COLOURS[5]);
		hV0RtFitCorr[iSp][0][0]->GetYaxis()->SetTitle("N_{V0} / <N_{V0}>");
		cRtYields[iSp] = new TCanvas(Form("cRtYields_%s_%s_%s",SPECIES[iSp],TYPE[0],MULTI[3]),"",1000,800);

		hV0RtFitCorr[iSp][0][0]->GetYaxis()->SetRangeUser(-0.1,10.);
		hV0RtFitCorr[iSp][0][0]->Draw();

		TLegend* legRtYield = new TLegend(0.15,0.65,0.38,0.88);
		mHandler->MakeNiceLegend(legRtYield,0.04,1);
			
		legRtYield->AddEntry((TObject*)0,Form("|#eta| < 0.8, %s", isMC),"");
		legRtYield->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legRtYield->AddEntry((TObject*)0,"","");
		legRtYield->AddEntry(hV0RtFitCorr[iSp][0][0],Form("%s", SPECNAMES[iSp]),"pl");
		legRtYield->Draw();

		funcLin->Draw("same");

		cRtYields[iSp]->Write();
		cRtYields[iSp]->SaveAs(Form("plots/rtyield_%s.png",SPECIES[iSp]));
	}*/

	//mHandler->root()->SetBatch(kFALSE);

}