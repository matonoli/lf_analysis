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
		
		hV0RtFitCorr[iSp][0][0] 
			= (TH1D*)mHandler->analysis(2)->dirFile()->Get(Form("hV0RtFitCorr_%s_%s_%s",SPECIES[iSp],TYPE[0],MULTI[3]));
	}

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
	//MakeFinalFiguresSpherocity();
	MakeFinalFiguresRt();

	return 0;	
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
			mHandler->MakeNiceHistogram(hV0PtRtFitCorr[iSp][0][iReg][iRtBin],COLOURS[5-iRtBin]);
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
	TCanvas* cBtoM[NREGIONS];
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
			
			mHandler->MakeNiceHistogram(hBtoMRt[iReg][iRtBin],COLOURS[5-iRtBin]);
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

		cBtoM[iReg]->Write();
		cBtoM[iReg]->SaveAs(Form("plots/btomrt_%s.png",REGIONS[iReg]));

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

				mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][iMu][iSph],COLOURS[iMu+iSph+(iSph>0)-(iMu>1&&iSph>0)]);
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
				cPt[iSp][iMu-1], 0.1,3.2);
			mHandler->MakeRatioPlot(hV0PtFitCorr[iSp][0][iMu][2],hV0PtFitCorr[iSp][0][iMu][0],
				cPt[iSp][iMu-1], 0.1,3.2);

			cPt[iSp][iMu-1]->Write();
			cPt[iSp][iMu-1]->SaveAs(Form("plots/pt_%s_%s.png",SPECIES[iSp],MULTI[iMu]));
		}
	}

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
	TCanvas* cRtYields[4];
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
	}

	//mHandler->root()->SetBatch(kFALSE);

}