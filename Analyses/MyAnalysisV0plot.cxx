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

	} } } }

}

Bool_t MyAnalysisV0plot::CreateHistograms() {
	

}


Bool_t MyAnalysisV0plot::CloneHistograms() {

	for (int iMu = 0; iMu < 2; ++iMu) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
		
		hBtoM[iMu][iSph] = (TH1D*)hV0PtFitCorr[2][0][iMu][iSph]->Clone(Form("hBtoM_%s_%s",MULTI[iMu],SPHERO[iSph]));

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

			mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][0][0],COLOURS[0]);
			mHandler->MakeNiceHistogram(hV0PtFitCorr[iSp][0][1][iSph],COLOURS[1+iSph]);
		}

		cPt[iSp] = new TCanvas(Form("cPt_%s",SPECIES[iSp]),"",1000,800);
		cPt[iSp]->SetLogy(1);
		hV0PtFitCorr[iSp][0][0][0]->GetYaxis()->SetRangeUser(0.1,10.*hV0PtFitCorr[iSp][0][0][0]->GetMaximum());
		hV0PtFitCorr[iSp][0][0][0]->Draw();
		cPt[iSp]->Update();
		hV0PtFitCorr[iSp][0][1][0]->Draw("same");
		hV0PtFitCorr[iSp][0][1][1]->Draw("same");
		hV0PtFitCorr[iSp][0][1][2]->Draw("same");

		TLegend* legPt = new TLegend(0.55,0.55,0.85,0.85);
		mHandler->MakeNiceLegend(legPt,0.04,1);
		legPt->AddEntry((TObject*)0,Form("%s   |#eta| < 0.8", SPECNAMES[iSp]),"");
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"","");
		legPt->AddEntry(hV0PtFitCorr[iSp][0][0][0],"MB","pl");
		legPt->AddEntry(hV0PtFitCorr[iSp][0][1][0],"FHM","pl");
		legPt->AddEntry(hV0PtFitCorr[iSp][0][1][1],"FHM Jetty","pl");
		legPt->AddEntry(hV0PtFitCorr[iSp][0][1][2],"FHM Iso","pl");
		legPt->Draw();

		cPt[iSp]->Write();
		cPt[iSp]->SaveAs(Form("plots/pt_%s.png",SPECIES[iSp]));
	}

	// B/M RATIO
	//TH1D* hBtoM[NMULTI][NSPHERO];
	for (int iMu = 0; iMu < 2; ++iMu) {
	for (int iSph = 0; iSph < NSPHERO; ++iSph) {
	
		
		hBtoM[iMu][iSph]->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / 2K^{0}_{s}");
		hBtoM[iMu][iSph]->Add(hV0PtFitCorr[3][0][iMu][iSph]);
		hBtoM[iMu][iSph]->Divide(hBtoM[iMu][iSph],hV0PtFitCorr[1][0][iMu][iSph],1.,2.,"");

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


}