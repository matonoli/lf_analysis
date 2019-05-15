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

	hEventSpherocity = (TH1D*)mHandler->analysis(0)->dirFile()->Get("hEventSpherocity");

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < 1; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		hV0PtFitCorr[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtFitCorr_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		hV0PtFit[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));

	} } } }

	for (int iType = 0; iType < 1; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		hTrackPt[iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hTrackPt_%s_%s_%s",TYPE[iType],MULTI[iMu],SPHERO[iSph]));

	} } }

}

Bool_t MyAnalysisV0plot::CreateHistograms() {
	

}


Bool_t MyAnalysisV0plot::CloneHistograms() {

	for (int iMu = 0; iMu < NMULTI; ++iMu) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
		
		hBtoM[iMu][iSph] =			(TH1D*)hV0PtFitCorr[2][0][iMu][iSph]->Clone(Form("hBtoM_%s_%s",MULTI[iMu],SPHERO[iSph]));
		hV0toNchDR[0][iMu][iSph] =	(TH1D*)hV0PtFit[1][0][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s",SPECIES[1],MULTI[iMu],SPHERO[iSph]));
		hV0toNchDR[1][iMu][iSph] =	(TH1D*)hV0PtFit[2][0][iMu][iSph]->Clone(Form("hV0toNchDR_%s_%s_%s",SPECIES[2],MULTI[iMu],SPHERO[iSph]));
		
	}	}

}

Int_t MyAnalysisV0plot::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();


	BorrowHistograms();
	CloneHistograms();
	MakeFinalFigures();

	return 0;	
}

void MyAnalysisV0plot::MakeFinalFigures() {

	//mFilePlots = new TFile("plots_"+mOutName,"RECREATE");
	//mFilePlots->cd();

	enum { LEFT, RIGHT };
	const char* isMC = (mHandler->GetFlagMC()) ? "MC^{blind}_{rec}" : "Data";

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
	TCanvas* cPt[4][2];
	for (Int_t iSp = 1; iSp < NSPECIES; ++iSp)	{				
	
		mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][0][0],COLOURS[0]);
		hV0PtFitCorr[iSp][0][0][0]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");

		for (Int_t iMu = 1; iMu < NMULTI; ++iMu)	{
			
			for (Int_t iSph = 0; iSph < NSPHERO; ++iSph)	{

				mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][iMu][iSph],COLOURS[iMu+iSph+(iSph>0)-(iMu>1&&iSph>0)]);
				hV0PtFitCorr[iSp][0][iMu][iSph]->GetYaxis()->SetTitle("1/N_{ev} N/(#Deltay #Deltap_{T})  ((GeV/#it{c})^{-1})");		

			}

			cPt[iSp][iMu-1] = new TCanvas(Form("cPt_%s_%s",SPECIES[iSp],MULTI[iMu]),"",1000,800);
			cPt[iSp][iMu-1]->SetLogy(1);
			Double_t lowerRange = 0.1*hV0PtFitCorr[iSp][0][0][0]->GetBinContent(hV0PtFitCorr[iSp][0][0][0]->FindLastBinAbove());
			hV0PtFitCorr[iSp][0][0][0]->GetYaxis()->SetRangeUser(lowerRange,10.*hV0PtFitCorr[iSp][0][0][0]->GetMaximum());
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

			cPt[iSp][iMu-1]->Write();
			cPt[iSp][iMu-1]->SaveAs(Form("plots/pt_%s_%s.png",SPECIES[iSp],MULTI[iMu]));
		}
	}

	// B/M RATIO
	//TH1D* hBtoM[NMULTI][NSPHERO];
	for (int iMu = 0; iMu < NMULTI; ++iMu) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
		
		hBtoM[iMu][iSph]->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / 2K^{0}_{s}");
		hBtoM[iMu][iSph]->Add(hV0PtFitCorr[3][0][iMu][iSph]);
		hBtoM[iMu][iSph]->Divide(hBtoM[iMu][iSph],hV0PtFitCorr[1][0][iMu][iSph],1.,2.,"");

	}	}
	
	TCanvas* cBtoM[3];
	TLegend* legBtoM[3];
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


	// Double-ratio
	for (int iSp2 = 0; iSp2 < 2; ++iSp2)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			for (int iMu = 0; iMu < NMULTI; ++iMu)		{

				if (iSp2==1) hV0toNchDR[iSp2][iMu][iSph]->Add(hV0PtFit[2][0][iMu][iSph]);
				hV0toNchDR[iSp2][iMu][iSph]->Divide(hTrackPt[0][iMu][iSph]);
			}	
	}	}	// now we have ratios of v0 to pi(charged tracks)

	for (int iSp2 = 0; iSp2 < 2; ++iSp2)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		hV0toNchDR[iSp2][1][iSph]->Divide(hV0toNchDR[iSp2][0][iSph]);
		hV0toNchDR[iSp2][2][iSph]->Divide(hV0toNchDR[iSp2][0][iSph]);
	}	}	// iMu!=0 are now (v0/pi)_hm  / (v0/pi)_mb double ratios

	TCanvas* cDR[2][NMULTI-1];
	TLegend* legDR[2][NMULTI-1];
	for (int iSp2 = 0; iSp2 < 2; ++iSp2)	{
	for (int iMu = 1; iMu < NMULTI; ++iMu)	{
		//hV0toNchDR[iSp2][iMu][iSph]
		cDR[iSp2][iMu-1] = new TCanvas(Form("cDR_%s_%s",SPECIES[1+iSp2],MULTI[iMu]),"",1000,800);
		cout << "iSp2 "<< iSp2 << " iMu " << iMu << " histo " <<  hV0toNchDR[iSp2][iMu][0] << endl;
		hV0toNchDR[iSp2][iMu][0]->GetYaxis()->SetRangeUser(0.5,4.);
		hV0toNchDR[iSp2][iMu][0]->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
		hV0toNchDR[iSp2][iMu][0]->GetYaxis()->SetTitle("ratio");
		mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][iMu][0],COLOURS[0]);
		mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][iMu][1],COLOURS[3]);
		mHandler->MakeNiceHistogram(hV0toNchDR[iSp2][iMu][2],COLOURS[4]);
		hV0toNchDR[iSp2][iMu][0]->Draw();
		hV0toNchDR[iSp2][iMu][1]->Draw("same");
		hV0toNchDR[iSp2][iMu][2]->Draw("same");

		legDR[iSp2][iMu-1] = new TLegend(0.55,0.55,0.85,0.85);
		mHandler->MakeNiceLegend(legDR[iSp2][iMu-1],0.04,1);
		legDR[iSp2][iMu-1]->AddEntry((TObject*)0,Form("pp #sqrt{s} = 13 TeV"),"");
		legDR[iSp2][iMu-1]->AddEntry((TObject*)0,Form("%s  %s", SPECNAMES[iSp2], PLOTS_MULTI[iMu]),"");
		legDR[iSp2][iMu-1]->AddEntry((TObject*)0,"","");
		legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][iMu][0],"HM","pl");
		legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][iMu][1],"HM Jetty","pl");
		legDR[iSp2][iMu-1]->AddEntry(hV0toNchDR[iSp2][iMu][2],"HM Iso","pl");
		legDR[iSp2][iMu-1]->Draw();

	}	}

}