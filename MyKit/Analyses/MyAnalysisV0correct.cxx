#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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
#include <TNtuple.h>
#include <TCutG.h>
#include <TCut.h>

#include "MyAnalysisV0correct.h"
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
using namespace std;

ClassImp(MyAnalysisV0correct)



MyAnalysisV0correct::MyAnalysisV0correct() {

}

Int_t MyAnalysisV0correct::Init() {

	TString dfName(this->GetName());
	dfName = Form("%s_%i",dfName.Data(),mHandler->nAnalysis());
	mDirFile = new TDirectoryFile(dfName,dfName,"",mHandler->file());
	mDirFile->cd();

	printf("Initialising analysis %s  \n", 
		this->GetName());

	TH1::SetDefaultSumw2(1);
	BorrowHistograms();
	CreateHistograms();
	


	mList = (TList*)mHandler->directory()->GetList();

	return 0;
}

Int_t MyAnalysisV0correct::Make(Int_t iEv) {
	return 0;
}


Bool_t MyAnalysisV0correct::BorrowHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	hEventType			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hEventType");
	hEventMultvSpheroD 	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hEventMultvSpheroD");
	hEventMultvSpheroMC	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hEventMultvSpheroMC");
	hNchTrans			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTrans");
	hRt					= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hRt");
	hRt2				= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hRt2");
	
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{
		hRtV0Yields[iType][iReg][iPtBin]	
			= (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hRtV0Yields_%s_%s_%i",TYPE[iType],REGIONS[iReg],iPtBin));
	} } }

	
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		if (iMu == 0 && iSph > 0) continue;	
		//if (iMu > 4 && (iSph < 3 && iSph)) continue;
		//if (iMu < 5 && iSph > 2) continue; 
		hV0PtFit[iSp][iType][iMu][iSph] 
			= (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));

	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNtFit[iSp][iType][iReg]		
			= (TH2F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtNtFit_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		
	}	}	}

	if (mHandler->GetFlagMC()) {
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	Int_t iType = 2;
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNt[iSp][iType][iReg]		
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0PtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
	}	}	
	}

	if (mHandler->GetFlagMC()) {
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		Int_t iType = 2; Int_t iMu = 0; Int_t iSph = 0;
		hV0Pt[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		
		if (hV0Pt[iSp][iType][iMu][iSph]->GetNbinsX() != NPTBINS) 
			hV0Pt[iSp][iType][iMu][iSph] = (TH1D*)hV0Pt[iSp][iType][iMu][iSph]->Rebin(NPTBINS,hV0Pt[iSp][iType][iMu][iSph]->GetName(),XBINS);
	
	}	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{
		hV0RtFit[iSp][iType][iReg][iPtBin] 
			= (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0RtFit_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iPtBin));
	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
			
		hV0PtRtFit[iSp][iType][iReg][iRtBin] 
			= (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtRtFit_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin) );


	} } } }		

	
	hNchTrans					= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTrans");
	hNchTransMC					= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMC");
	hNchTransMCTrigMC			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMCTrigMC");
	hNchTransRC					= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTransRC");
	hNchTransRCvMC				= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hNchTransRCvMC");
	hNchTransMinRCvMC				= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hNchTransMinRCvMC");
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 2; iType < NTYPE; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
				
		hV0PtNt[iSp][iType][iReg] = (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0PtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hV0PtNt[iSp][iType][iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0PtNt[iSp][iType][iReg]);

	}	}	}

	if (mHandler->GetFlagMC()) {

		hPiPtMC			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hPiPtMC");
		hPiPtRC			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hPiPtRC");
		hKpmPtMC			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hKpmPtMC");
		hKpmPtRC			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hKpmPtRC");
		for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
					
			hPiPtNtRC[iReg] = (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hPiPtNtRC_%s",REGIONS[iReg]));
			hPiPtNtRC[iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hPiPtNtRC[iReg]);
			hPiPtNtMC[iReg] = (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hPiPtNtMC_%s",REGIONS[iReg]));
			hPiPtNtMC[iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hPiPtNtMC[iReg]);

			hKpmPtNtRC[iReg] = (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hKpmPtNtRC_%s",REGIONS[iReg]));
			hKpmPtNtRC[iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hKpmPtNtRC[iReg]);
			hKpmPtNtMC[iReg] = (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hKpmPtNtMC_%s",REGIONS[iReg]));
			hKpmPtNtMC[iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hKpmPtNtMC[iReg]);
		}
	}	



}

Bool_t MyAnalysisV0correct::CreateHistograms() {
	
}


Bool_t MyAnalysisV0correct::CloneHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		if (iMu == 0 && iSph > 0) continue;		
		//if (iMu > 4 && (iSph < 3 && iSph)) continue;
		//if (iMu < 5 && iSph > 2) continue; 
		hV0PtFitCorr[iSp][iType][iMu][iSph]	= (TH1F*)hV0PtFit[iSp][iType][iMu][iSph]->Clone(
			Form("hV0PtFitCorr_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]) );

		//if (hV0PtFitCorr[iSp][iType][iMu][iSph] && mHandler->IsRebinPt()) {
		//hV0PtFitCorr[iSp][iType][iMu][iSph] = (TH1F*)hV0PtFitCorr[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName(),XBINS2); }
		
	} } } }



	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNtFitCorr[iSp][iType][iReg]		= (TH2F*)hV0PtNtFit[iSp][iType][iReg]->Clone(
			Form("hV0PtNtFitCorr_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg])	);

	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
			
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin] = (TH1F*)hV0PtRtFit[iSp][iType][iReg][iRtBin]->Clone(
			Form("hV0PtRtFitCorr_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin) );


	} } } }

	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		if (!iMu && iSph) continue;
		hV0PtFeeddown[iSp][iMu][iSph]	= (TH1F*)hV0PtFitCorr[iSp][0][iMu][0]->Clone(
			Form("hV0PtFeeddown_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtFeeddownXi0[iSp][iMu][iSph]	= (TH1F*)hV0PtFitCorr[iSp][0][iMu][0]->Clone(
			Form("hV0PtFeeddownXi0_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtFeeddownXiErr[iSp][iMu][iSph]	= (TH1F*)hV0PtFitCorr[iSp][0][iMu][0]->Clone(
			Form("hV0PtFeeddownXiErr_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
	}	}	}

	if (mHandler->GetFlagMC()) {
		for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
			int iMu = 0;
			
			hClosureTestFDvSec[iSp] = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitSecondaryXi_%s",SPECIES[iSp]))->Clone(Form("hClosureTestFDvSec_%s",SPECIES[iSp]));
			hClosureTestFDvSecPDG[iSp] = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitSecondaryPDG_%s",SPECIES[iSp]))->Clone(Form("hClosureTestFDvSecPDG_%s",SPECIES[iSp]));

			hClosureTestPDG[iSp] = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitPrimaryPDG_%s",SPECIES[iSp]))->Clone(Form("hClosureTestPDG_%s",SPECIES[iSp]));
			hClosureTestPDG[iSp]->Add((TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitSecondaryPDG_%s",SPECIES[iSp])));
			cout << "AWBEALWKEBAWLEBAL " <<  hClosureTestPDG[iSp] << endl;
		}	

		for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
			int iMu = 0;
			
			hClosureTestEffi[iSp] = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitPrimary_%s",SPECIES[iSp]))->Clone(Form("hClosureTestEffi_%s",SPECIES[iSp]));
			hClosureTestEffiPDG[iSp] = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitPrimaryPDG_%s",SPECIES[iSp]))->Clone(Form("hClosureTestEffiPDG_%s",SPECIES[iSp]));
		
		}	
	}

}

Int_t MyAnalysisV0correct::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	CloneHistograms();
	hNchTrans->Write();
	hNchTransMC->Write();
	hNchTransMCTrigMC->Write();
	hNchTransRC->Write();
	hNchTransRCvMC->Write();
	hNchTransMinRCvMC->Write();
	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 2; iType < NTYPE; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
		
		// SCALE BY BIN WIDTH
		for (int iX = 1; iX < hV0PtNt[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNt[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float binwidth = hV0PtNt[iSp][iType][iReg]->GetXaxis()->GetBinWidth(iX);
				float binc = hV0PtNt[iSp][iType][iReg]->GetBinContent(iX,iY);
				float bine = hV0PtNt[iSp][iType][iReg]->GetBinError(iX,iY);
				if (binwidth>0) hV0PtNt[iSp][iType][iReg]->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) hV0PtNt[iSp][iType][iReg]->SetBinError(iX,iY,bine/binwidth);
			}
		}
		//hV0PtNt[iSp][iType][iReg]->Write();

	}	}	}


	// CORRECT PI AND KPM
	if (mHandler->GetFlagMC()) {
		hPiPtRC->Divide(hPiPtMC);
		hKpmPtRC->Divide(hKpmPtMC);
		TH2F* htmp;

		for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
		
			// SCALE BY BIN WIDTH
			for (int iX = 1; iX < hPiPtNtRC[iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hPiPtNtRC[iReg]->GetNbinsY()+1; iY++) {	

				htmp = hPiPtNtRC[iReg];
				float binwidth = htmp->GetXaxis()->GetBinWidth(iX);
				float binc = htmp->GetBinContent(iX,iY);
				float bine = htmp->GetBinError(iX,iY);
				if (binwidth>0) htmp->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) htmp->SetBinError(iX,iY,bine/binwidth);

				htmp = hPiPtNtMC[iReg];
				binwidth = htmp->GetXaxis()->GetBinWidth(iX);
				binc = htmp->GetBinContent(iX,iY);
				bine = htmp->GetBinError(iX,iY);
				if (binwidth>0) htmp->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) htmp->SetBinError(iX,iY,bine/binwidth);

				htmp = hKpmPtNtRC[iReg];
				binwidth = htmp->GetXaxis()->GetBinWidth(iX);
				binc = htmp->GetBinContent(iX,iY);
				bine = htmp->GetBinError(iX,iY);
				if (binwidth>0) htmp->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) htmp->SetBinError(iX,iY,bine/binwidth);

				htmp = hKpmPtNtMC[iReg];
				binwidth = htmp->GetXaxis()->GetBinWidth(iX);
				binc = htmp->GetBinContent(iX,iY);
				bine = htmp->GetBinError(iX,iY);
				if (binwidth>0) htmp->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) htmp->SetBinError(iX,iY,bine/binwidth);
			}	}

			// CORRECT FOR EFFI
			for (int iX = 1; iX < hPiPtNtRC[iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hPiPtNtRC[iReg]->GetNbinsY()+1; iY++) {
				
				hPiPtNtRC[iReg]->SetBinContent(iX,iY, hPiPtRC->GetBinContent(iX)>0 ? 
					(float)hPiPtNtRC[iReg]->GetBinContent(iX,iY)/hPiPtRC->GetBinContent(iX) : 0);
				hKpmPtNtRC[iReg]->SetBinContent(iX,iY, hKpmPtRC->GetBinContent(iX)>0 ? 
					(float)hKpmPtNtRC[iReg]->GetBinContent(iX,iY)/hKpmPtRC->GetBinContent(iX) : 0);

				// errors ?
			}	}
		}
	}

	
	if (MAKE_EXCLUSIVE) MakeExclusiveS0Bins(); 
	NormaliseSpectra();
	//CorrectForFeeddown();

	//LoadEfficiency();
	DoEfficiencyFromFile();
	CorrectSpectra();

	//if (!mHandler->GetFlagMC()) StudyCuts();

	//if (NMULTI>2) DoXCheckV0M();
	//if (mHandler->GetFlagMC()) DoClosureTest(0);


	CreateOutputFile("k0s_spherocity.root",1);

	printf("mb k0s spectrum final \n");
	cout << hV0PtFitCorr[1][0][0][0]->GetBinContent(30) << " at " << hV0PtFitCorr[1][0][0][0]->GetBinLowEdge(30) << endl;

	return 0;
	//draw shit, move to func later
	{
	//fast det
		TH2F* hV0BaselineIMvPt[NSPECIES][NTYPE];
		TH2F* hV0FastSignalIMvPt[NSPECIES][NTYPE];
		TH2F* hV0CutIMvPt[NSPECIES][NTYPE][25];
		TDirectoryFile* dirFile1 = new TDirectoryFile("mcFile2","mcFile2","",mHandler->file());
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("TDirectoryFile")) {
		cout << "Doing efficiency from a TDirectoryFile" << endl;
		dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");}
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("THashList")) {
		cout << "Doing efficiency from a THashList" << endl;
		THashList* hashList = (THashList*)mFileMC->Get("MyAnalysisV0_0");
		while (hashList->GetEntries()) {
			dirFile1->Append(hashList->First());
			hashList->RemoveFirst();
		}
	}
		int iSp = 1; TH1F* hd; TH1F* hmc;
		hV0FastSignalIMvPt[iSp][0]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0FastSignalIMvPt_K0s_D");
		hV0BaselineIMvPt[iSp][0]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0BaselineIMvPt_K0s_D");
		hV0FastSignalIMvPt[iSp][0]->Divide(hV0BaselineIMvPt[iSp][0]);
		hV0FastSignalIMvPt[iSp][1]	= (TH2F*)dirFile1->Get("hV0FastSignalIMvPt_K0s_D");
		hV0BaselineIMvPt[iSp][1]	= (TH2F*)dirFile1->Get("hV0BaselineIMvPt_K0s_D");
		hV0FastSignalIMvPt[iSp][1]->Divide(hV0BaselineIMvPt[iSp][1]);
		TCanvas* cn1 = new TCanvas("cn1","",1200,800); cn1->Divide(2,2,0.00005,0.00005);
		cn1->cd(1); hV0FastSignalIMvPt[iSp][0]->SetStats(0); hV0FastSignalIMvPt[iSp][0]->Draw("colz");
		cn1->cd(2); hV0FastSignalIMvPt[iSp][1]->SetStats(0); hV0FastSignalIMvPt[iSp][1]->Draw("colz");
		hd 	= (TH1F*)hV0FastSignalIMvPt[iSp][0]->ProjectionX("hd",13,13);
		hmc = (TH1F*)hV0FastSignalIMvPt[iSp][1]->ProjectionX("hmc",13,13);
		mHandler->MakeNiceHistogram(hd,kRed); mHandler->MakeNiceHistogram(hmc,kRed);
		hd->GetYaxis()->SetRangeUser(0.,1.2); hmc->GetYaxis()->SetRangeUser(0.,1.2);
		cn1->cd(3); hd->Draw();
		cn1->cd(4); hmc->Draw();
		cn1->Draw(); cn1->SaveAs("plots/pileup.png");

		hV0CutIMvPt[iSp][0][19]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_19");
		hV0CutIMvPt[iSp][0][15]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_15");
		hV0CutIMvPt[iSp][0][19]->Divide(hV0CutIMvPt[iSp][0][15]);
		hV0CutIMvPt[iSp][1][19]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_19");
		hV0CutIMvPt[iSp][1][15]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_15");
		hV0CutIMvPt[iSp][1][19]->Divide(hV0CutIMvPt[iSp][1][15]);
		TCanvas* cn2 = new TCanvas("cn2","",1200,800); cn2->Divide(2,2,0.00005,0.00005);
		cn2->cd(1); hV0CutIMvPt[iSp][0][19]->SetStats(0); hV0CutIMvPt[iSp][0][19]->Draw("colz");
		cn2->cd(2); hV0CutIMvPt[iSp][1][19]->SetStats(0); hV0CutIMvPt[iSp][1][19]->Draw("colz");
		hd 	= (TH1F*)hV0CutIMvPt[iSp][0][19]->ProjectionX("hd",13,13);
		hmc = (TH1F*)hV0CutIMvPt[iSp][1][19]->ProjectionX("hmc",13,13);
		mHandler->MakeNiceHistogram(hd,kRed); mHandler->MakeNiceHistogram(hmc,kRed);
		hd->GetYaxis()->SetRangeUser(0.,1.2); hmc->GetYaxis()->SetRangeUser(0.,1.2);
		cn2->cd(3); hd->Draw();
		cn2->cd(4); hmc->Draw();
		cn2->Draw(); cn2->SaveAs("plots/nsigma.png");

		hV0CutIMvPt[iSp][0][15]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_15");
		hV0CutIMvPt[iSp][0][13]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_13");
		hV0CutIMvPt[iSp][0][15]->Divide(hV0CutIMvPt[iSp][0][13]);
		hV0CutIMvPt[iSp][1][15]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_15");
		hV0CutIMvPt[iSp][1][13]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_13");
		hV0CutIMvPt[iSp][1][15]->Divide(hV0CutIMvPt[iSp][1][13]);
		TCanvas* cn5 = new TCanvas("cn5","",1200,800); cn5->Divide(2,2,0.00005,0.00005);
		cn5->cd(1); hV0CutIMvPt[iSp][0][15]->SetStats(0); hV0CutIMvPt[iSp][0][15]->Draw("colz");
		cn5->cd(2); hV0CutIMvPt[iSp][1][15]->SetStats(0); hV0CutIMvPt[iSp][1][15]->Draw("colz");
		hd 	= (TH1F*)hV0CutIMvPt[iSp][0][15]->ProjectionX("hd",13,13);
		hmc = (TH1F*)hV0CutIMvPt[iSp][1][15]->ProjectionX("hmc",13,13);
		mHandler->MakeNiceHistogram(hd,kRed); mHandler->MakeNiceHistogram(hmc,kRed);
		hd->GetYaxis()->SetRangeUser(0.,1.2); hmc->GetYaxis()->SetRangeUser(0.,1.2);
		cn5->cd(3); hd->Draw();
		cn5->cd(4); hmc->Draw();
		cn5->Draw(); cn5->SaveAs("plots/compmass.png");

		hV0CutIMvPt[iSp][0][13]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_13");
		hV0CutIMvPt[iSp][0][12]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_12");
		hV0CutIMvPt[iSp][0][13]->Divide(hV0CutIMvPt[iSp][0][12]);
		hV0CutIMvPt[iSp][1][13]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_13");
		hV0CutIMvPt[iSp][1][12]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_12");
		hV0CutIMvPt[iSp][1][13]->Divide(hV0CutIMvPt[iSp][1][12]);
		TCanvas* cn3 = new TCanvas("cn3","",1200,800); cn3->Divide(2,2,0.00005,0.00005);
		cn3->cd(1); hV0CutIMvPt[iSp][0][13]->SetStats(0); hV0CutIMvPt[iSp][0][13]->Draw("colz");
		cn3->cd(2); hV0CutIMvPt[iSp][1][13]->SetStats(0); hV0CutIMvPt[iSp][1][13]->Draw("colz");
		hd 	= (TH1F*)hV0CutIMvPt[iSp][0][13]->ProjectionX("hd",13,13);
		hmc = (TH1F*)hV0CutIMvPt[iSp][1][13]->ProjectionX("hmc",13,13);
		mHandler->MakeNiceHistogram(hd,kRed); mHandler->MakeNiceHistogram(hmc,kRed);
		hd->GetYaxis()->SetRangeUser(0.,1.2); hmc->GetYaxis()->SetRangeUser(0.,1.2);
		cn3->cd(3); hd->Draw();
		cn3->cd(4); hmc->Draw();
		cn3->Draw(); cn3->SaveAs("plots/lifetime.png");

		hV0CutIMvPt[iSp][0][9]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_9");
		hV0CutIMvPt[iSp][0][8]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get("hV0CutIMvPt_K0s_D_8");
		hV0CutIMvPt[iSp][0][9]->Divide(hV0CutIMvPt[iSp][0][8]);
		hV0CutIMvPt[iSp][1][9]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_9");
		hV0CutIMvPt[iSp][1][8]	= (TH2F*)dirFile1->Get("hV0CutIMvPt_K0s_D_8");
		hV0CutIMvPt[iSp][1][9]->Divide(hV0CutIMvPt[iSp][1][8]);
		TCanvas* cn4 = new TCanvas("cn4","",1200,800); cn4->Divide(2,2,0.00005,0.00005);
		cn4->cd(1); hV0CutIMvPt[iSp][0][9]->SetStats(0); hV0CutIMvPt[iSp][0][9]->Draw("colz");
		cn4->cd(2); hV0CutIMvPt[iSp][1][9]->SetStats(0); hV0CutIMvPt[iSp][1][9]->Draw("colz");
		hd 	= (TH1F*)hV0CutIMvPt[iSp][0][9]->ProjectionX("hd",13,13);
		hmc = (TH1F*)hV0CutIMvPt[iSp][1][9]->ProjectionX("hmc",13,13);
		mHandler->MakeNiceHistogram(hd,kRed); mHandler->MakeNiceHistogram(hmc,kRed);
		hd->GetYaxis()->SetRangeUser(0.,1.2); hmc->GetYaxis()->SetRangeUser(0.,1.2);
		cn4->cd(3); hd->Draw();
		cn4->cd(4); hmc->Draw();
		cn4->Draw(); cn4->SaveAs("plots/radiusUP.png");
		
	}


	return 0;	
}

void MyAnalysisV0correct::SetMCInputFile(const Char_t *name) {

	TString fileName = TString(name);
	if (fileName.Data() != "") {
		mFileMC = new TFile(fileName,"READ");
		printf("MC File %s loaded in. \n", fileName.Data()); }
	else {
		printf("No MC file loaded.");
	}
}

void MyAnalysisV0correct::SetXiSpectraFile(const Char_t *name) {

	TString fileName = TString(name);
	if (fileName.Data() != "") {
		//mFileXi = new TFile(fileName,"READ");
		mFileXi = new TFile("../official/xi_HM_spectra_sep_9_2019.root","READ");
		printf("Xi spectrum file %s loaded in. \n", fileName.Data()); }
	else {
		printf("No Xi file loaded.");
	}

}


void MyAnalysisV0correct::MakeExclusiveS0Bins() {

	// modifies the event count matrix so that exclusive bins are normalised correctly
	// creates exclusive spherocity bins, i.e 0-1, 1-5, instead of 0-1, 0-5

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	enum { sphMB, Jetty20, Iso20, Jetty10, Iso10, Jetty5, Iso5, Jetty1, Iso1};

	TH2F* hM;

	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 1; iMu < NMULTI; ++iMu)		{

		hM = (iType==1) ? hEventMultvSpheroMC : hEventMultvSpheroD;
		
		hM->SetBinContent(1+sphMB,1+iMu,hM->GetBinContent(1+sphMB,1+iMu)-hM->GetBinContent(1+Jetty20,1+iMu));
		hM->SetBinContent(1+sphMB,1+iMu,hM->GetBinContent(1+sphMB,1+iMu)-hM->GetBinContent(1+Iso20,1+iMu));

		hM->SetBinContent(1+Jetty20,1+iMu,hM->GetBinContent(1+Jetty20,1+iMu)-hM->GetBinContent(1+Jetty10,1+iMu));
		hM->SetBinContent(1+Jetty10,1+iMu,hM->GetBinContent(1+Jetty10,1+iMu)-hM->GetBinContent(1+Jetty5,1+iMu));
		hM->SetBinContent(1+Jetty5,1+iMu,hM->GetBinContent(1+Jetty5,1+iMu)-hM->GetBinContent(1+Jetty1,1+iMu));

		hM->SetBinContent(1+Iso20,1+iMu,hM->GetBinContent(1+Iso20,1+iMu)-hM->GetBinContent(1+Iso10,1+iMu));
		hM->SetBinContent(1+Iso10,1+iMu,hM->GetBinContent(1+Iso10,1+iMu)-hM->GetBinContent(1+Iso5,1+iMu));
		hM->SetBinContent(1+Iso5,1+iMu,hM->GetBinContent(1+Iso5,1+iMu)-hM->GetBinContent(1+Iso1,1+iMu));

	}	}

	hM->Write();

}

void MyAnalysisV0correct::CorrectForFeeddown() {

	NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
	Double_t NormEv = hEventType->GetBinContent(6);	// MB
	if (NormEv>0) {
		NormEv += NormEv * hEventType->GetBinContent(4) * 1./(hEventType->GetBinContent(5) + hEventType->GetBinContent(6));
	}

	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)		{
		hXiPt[iSp][iMu][iSph] = 0x0;
	}	}	}

	const int NBINS_FDMASS = 50;
	double XBINS_FDMASS[NBINS_FDMASS+1];// = {};
	double binStep = (0.05 + 0.05)/(double)NBINS_FDMASS;
	for (int i=0; i <= NBINS_FDMASS; i++) XBINS_FDMASS[i] = -0.05 + (double)i*binStep;

	TCutG* cutg;
	TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
		hV0FeeddownMatrix[iSp]		= (TH3F*)dirFile1->Get(Form("hV0FeeddownMatrix_%s",SPECIES[iSp]) );
		hV0FeeddownMotherPt[iSp]	= (TH1F*)dirFile1->Get(Form("hV0FeeddownMotherPt_%s",SPECIES[iSp]) );
		hV0FeeddownMatrixXi0[iSp]	= (TH3F*)dirFile1->Get(Form("hV0FeeddownMatrixXi0_%s",SPECIES[iSp]) );
		hV0FeeddownMotherPtXi0[iSp]	= (TH1F*)dirFile1->Get(Form("hV0FeeddownMotherPtXi0_%s",SPECIES[iSp]) );
		if (!hV0FeeddownMatrix[iSp]) {
			printf("Feed-down matrix not found, no correction performed.\n");
			return;
		}

		if (hV0FeeddownMatrix[iSp]->GetNbinsX() != NPTBINS || hV0FeeddownMatrix[iSp]->GetNbinsZ() != NBINS_FDMASS) {

			cout << "Rebinning FD matrix to right dimensions." << endl;
			
			TH3F *htmp	= new TH3F(Form("hV0FeeddownMatrix_%s",SPECIES[iSp]),";primary grandmother p_{T}; decay V0 p_{T}; mass p_{T}", NXIPTBINS, XIXBINS, NPTBINS, XBINS, NBINS_FDMASS, XBINS_FDMASS);
			TH3F *htmp2	= new TH3F(Form("hV0FeeddownMatrixXi0_%s",SPECIES[iSp]),";primary grandmother p_{T}; decay V0 p_{T}; mass p_{T}", NXIPTBINS, XIXBINS, NPTBINS, XBINS, NBINS_FDMASS, XBINS_FDMASS);
			TAxis *xaxis = hV0FeeddownMatrix[iSp]->GetXaxis(); 
			TAxis *yaxis = hV0FeeddownMatrix[iSp]->GetYaxis();
			TAxis *zaxis = hV0FeeddownMatrix[iSp]->GetZaxis(); 
			for (int k=1; k<=zaxis->GetNbins();k++)	{ 
			for (int j=1; j<=yaxis->GetNbins();j++)	{ 
			for (int i=1; i<=xaxis->GetNbins();i++)	{ 
				htmp->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),zaxis->GetBinCenter(k),hV0FeeddownMatrix[iSp]->GetBinContent(i,j,k)); 
				htmp2->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),zaxis->GetBinCenter(k),hV0FeeddownMatrixXi0[iSp]->GetBinContent(i,j,k)); 
			}	}	}
			
			hV0FeeddownMatrix[iSp] = (TH3F*)htmp->Clone(hV0FeeddownMatrix[iSp]->GetName());
			hV0FeeddownMatrixXi0[iSp] = (TH3F*)htmp2->Clone(hV0FeeddownMatrixXi0[iSp]->GetName());
			delete htmp;
			delete htmp2;
			
		}

		// NORMALIZING 
		Int_t nCols = hV0FeeddownMatrix[iSp]->GetNbinsX();
		Int_t nRows = hV0FeeddownMatrix[iSp]->GetNbinsY();
		Int_t nMass = hV0FeeddownMatrix[iSp]->GetNbinsZ();


		/*hV0FeeddownMotherPt[iSp]->Scale(1.,"width");
		hV0FeeddownMotherPt[iSp]->Scale(1./NormEv);
		hV0FeeddownMotherPt[iSp]->Scale(1./NormEta);
		hV0FeeddownMotherPt[iSp]->Scale(1.,"width");
		hV0FeeddownMotherPt[iSp]->Scale(1./NormEv);
		hV0FeeddownMotherPt[iSp]->Scale(1./NormEta);*/

		for (int iC = 1; iC < nCols+1; ++iC)	{
			//Double_t integral = hV0FeeddownMatrix[iSp]->ProjectionY("",iC,iC)->Integral(0,999);
			Double_t integral = hV0FeeddownMotherPt[iSp]->Integral(iC,iC);//,1,nRows);
			Double_t integralXi0 = hV0FeeddownMotherPtXi0[iSp]->Integral(iC,iC);//,1,nRows);
			for (int iR = 1; iR < nRows+1; ++iR)	{
			for (int iM = 1; iM < nMass+1; ++iM)	{
				
				Double_t binContent = hV0FeeddownMatrix[iSp]->GetBinContent(iC,iR,iM);
				if (binContent>0) hV0FeeddownMatrix[iSp]->SetBinContent(iC,iR,iM,binContent/integral);

				binContent = hV0FeeddownMatrixXi0[iSp]->GetBinContent(iC,iR,iM);
				if (binContent>0) hV0FeeddownMatrixXi0[iSp]->SetBinContent(iC,iR,iM,binContent/integral); //both methods should be normalized by charged xi only
			}	}
		}

		cutg = new TCutG("CUTG",11);
		cutg->SetPoint(0,1.10641,0.395916);
		cutg->SetPoint(1,4.12183,2.70515);
		cutg->SetPoint(2,8.05547,4.68511);
		cutg->SetPoint(3,15.1932,9.98092);
		cutg->SetPoint(4,15.0045,15.1427);
		cutg->SetPoint(5,9.03038,10.792);
		cutg->SetPoint(6,4.89268,5.94943);
		cutg->SetPoint(7,4.05381,4.56584);
		cutg->SetPoint(8,1.61654,1.82252);
		cutg->SetPoint(9,0.346886,0.395916);
		cutg->SetPoint(10,1.10641,0.395916);

		hV0FeeddownMatrix[iSp]->Write();
		hV0FeeddownMatrixXi0[iSp]->Write();
	}

	TH2F* hFMatrix = (TH2F*)hV0FeeddownMatrix[2]->Project3D("xy")->Clone("hFMatrix");
	/*for (int iMotherBin = 1; iMotherBin < hV0FeeddownMatrix[2]->GetNbinsX()+1; ++iMotherBin)	{
		for (int iPtBin = 1; iPtBin < hV0FeeddownMatrix[2]->GetNbinsY()+1; ++iPtBin)	{
				
				Double_t cX = hV0FeeddownMatrix[2]->GetXaxis()->GetBinCenter(iMotherBin);
				Double_t cY = hV0FeeddownMatrix[2]->GetYaxis()->GetBinCenter(iPtBin);
				hFMatrix->SetBinContent(iMotherBin,iPtBin,0);
				//if (!cutg->IsInside(cX,cY)) continue;
				hFMatrix->SetBinContent(iMotherBin,iPtBin,hV0FeeddownMatrix[2]->GetBinContent(iMotherBin,iPtBin));
	}	}*/

	//mFileXi = (mHandler->GetFlagMC()) ? 0x0 : new TFile("../official/xi_HM_spectra_sep_9_2019.root","READ");
	mFileXi = (mHandler->GetFlagMC()) ? 0x0 : new TFile("../official/xi_results_So_V0M_aug_6.root","READ");

	if (!mFileXi) {
		printf("No Xi file loaded in, using MC Xi spectra instead. \n");
		hXiPt[2][0][0] = hV0FeeddownMotherPt[2];
		hXiPt[3][0][0] = hV0FeeddownMotherPt[3];

		//hXiPt[2][0]->Rebin()

		NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
		Double_t NormEv = hEventType->GetBinContent(6);	// MB
		if (NormEv>0) {
			NormEv += NormEv * hEventType->GetBinContent(4) * 1./(hEventType->GetBinContent(5) + hEventType->GetBinContent(6));
		}
		hXiPt[2][0][0]->Scale(1.,"width");
		hXiPt[2][0][0]->Scale(1./NormEv);
		hXiPt[2][0][0]->Scale(1./NormEta);
		hXiPt[3][0][0]->Scale(1.,"width");
		hXiPt[3][0][0]->Scale(1./NormEv);
		hXiPt[3][0][0]->Scale(1./NormEta);

	} else {

		const char* xiNames[9] = {"hXi_Ref_Spectra_stat", 
		"hXi_Jetty20_Spectra_stat", "hXi_Iso20_Spectra_stat", "hXi_Jetty10_Spectra_stat", "hXi_Iso10_Spectra_stat",
		"hXi_Jetty5_Spectra_stat", "hXi_Iso5_Spectra_stat", "hXi_Jetty1_Spectra_stat", "hXi_Iso1_Spectra_stat"};

		for (int iS = 0; iS < NSPHERO; iS++) {
			if (NMULTI<2) continue;
			hXiPt[2][1][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
			hXiPt[3][1][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
		}

		mFileXi = new TFile("../official/xi_results_So_CL1_aug_6.root","READ");
		for (int iS = 0; iS < NSPHERO; iS++) {
			if (NMULTI<3) continue;
			hXiPt[2][2][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
			hXiPt[3][2][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
		}

		mFileXi = new TFile("../official/xi_results_So_V0M_top1_aug_6.root","READ");
		for (int iS = 0; iS < NSPHERO; iS++) {
			if (NMULTI<4) continue;
			hXiPt[2][3][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
			hXiPt[3][3][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
		}

		mFileXi = new TFile("../official/xi_results_So_CL1_top1_aug_6.root","READ");
		for (int iS = 0; iS < NSPHERO; iS++) {
			if (NMULTI<5) continue;
			hXiPt[2][4][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
			hXiPt[3][4][iS] = (TH1F*)mFileXi->Get(xiNames[iS]);
		}


		TFile* fileXIMB = new TFile("../official/xi_MB_spectra_sep_9_2019_sum_rebin.root","READ");
		hXiPt[2][0][0] = (TH1F*)fileXIMB->Get("hXiSumSpectrum_MB");
		hXiPt[3][0][0] = (TH1F*)fileXIMB->Get("hXiSumSpectrum_MB");

		mDirFile->cd();
	}

	// SUBTRACT FOR EXCLUSIVE S0 SPECTRA
	enum { sphMB, Jetty20, Iso20, Jetty10, Iso10, Jetty5, Iso5, Jetty1, Iso1};
	if (MAKE_EXCLUSIVE) {

		for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
		for (int iMu = 1; iMu < NMULTI; ++iMu)		{

			hXiPt[iSp][iMu][sphMB]->Scale(100.);
			hXiPt[iSp][iMu][sphMB]->Add(hXiPt[iSp][iMu][Jetty20],-20.);
			hXiPt[iSp][iMu][sphMB]->Add(hXiPt[iSp][iMu][Iso20],-20.);
			hXiPt[iSp][iMu][sphMB]->Scale(1./60.);

			hXiPt[iSp][iMu][Jetty20]->Scale(20.);
			hXiPt[iSp][iMu][Jetty20]->Add(hXiPt[iSp][iMu][Jetty10],-10.);
			hXiPt[iSp][iMu][Jetty20]->Scale(1./10.);

			hXiPt[iSp][iMu][Jetty10]->Scale(10.);
			hXiPt[iSp][iMu][Jetty10]->Add(hXiPt[iSp][iMu][Jetty5],-5.);
			hXiPt[iSp][iMu][Jetty10]->Scale(1./5.);

			hXiPt[iSp][iMu][Jetty5]->Scale(5.);
			hXiPt[iSp][iMu][Jetty5]->Add(hXiPt[iSp][iMu][Jetty1],-1.);
			hXiPt[iSp][iMu][Jetty5]->Scale(1./4.);


			hXiPt[iSp][iMu][Iso20]->Scale(20.);
			hXiPt[iSp][iMu][Iso20]->Add(hXiPt[iSp][iMu][Iso10],-10.);
			hXiPt[iSp][iMu][Iso20]->Scale(1./10.);

			hXiPt[iSp][iMu][Iso10]->Scale(10.);
			hXiPt[iSp][iMu][Iso10]->Add(hXiPt[iSp][iMu][Iso5],-5.);
			hXiPt[iSp][iMu][Iso10]->Scale(1./5.);

			hXiPt[iSp][iMu][Iso5]->Scale(5.);
			hXiPt[iSp][iMu][Iso5]->Add(hXiPt[iSp][iMu][Iso1],-1.);
			hXiPt[iSp][iMu][Iso5]->Scale(1./4.);

		}	}
	}

	// INTERPOLATE XI SPECTRA
	Double_t xixibarFactor = (mFileXi) ? 2. : 1.;
	TF1* funcLT = LevyTsallis("LT",XIMASS);
	TF1* funcLThi = LevyTsallis("LThi",XIMASS);
	TF1* funcLTlo = LevyTsallis("LTlo",XIMASS);
	if (!mFileXi) {
		funcLT->SetParameter(1,15);
		funcLT->SetParameter(2,0.5);
		funcLT->SetParameter(3,0.05); }

	TF1* mParMuL = (TF1*)mHandler->analysis(1)->dirFile()->Get("funcMuL");
	TF1* mParSigL = (TF1*)mHandler->analysis(1)->dirFile()->Get("funcSigL");
	mParMuL->Write(); mParSigL->Write();
	Double_t nSig = 5.;


	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	
		if (!hXiPt[iSp][iMu][iSph]) continue;

		cout << iSp << " " << iMu << " " << iSph << " " << hXiPt[iSp][iMu][iSph] << endl;
		//if (mFileXi) hXiPt[iSp][iMu]->Scale(0.5);
		hXiPt[iSp][iMu][iSph]->Fit(funcLT,"");
		cout << hXiPt[iSp][iMu][iSph]->GetName() << " successfully fitted" << endl;
		funcLThi->SetParameter(3,funcLT->GetParameter(3)+funcLT->GetParError(3));
		funcLThi->SetParameter(1,funcLT->GetParameter(1)-funcLT->GetParError(1));
		funcLThi->SetParameter(2,funcLT->GetParameter(2)-funcLT->GetParError(2));
		funcLTlo->SetParameter(3,funcLT->GetParameter(3)-funcLT->GetParError(3));
		funcLTlo->SetParameter(1,funcLT->GetParameter(1)+funcLT->GetParError(1));
		funcLTlo->SetParameter(2,funcLT->GetParameter(2)+funcLT->GetParError(2));
		//funcLThi->Write();
		//funcLTlo->Write();

		hXiPt[iSp][iMu][iSph]->Write();

		if (hV0PtFeeddown[iSp][iMu][iSph]->GetNbinsX() != hV0FeeddownMatrix[iSp]->GetNbinsY()) {
			printf("FD matrix has bad dimensions!\n");
			return;
		}

		for (int iBin = 1; iBin < hV0PtFeeddown[iSp][iMu][iSph]->GetNbinsX()+1; ++iBin)	{
			Double_t sum = 0;
			Double_t sumErr = 0;
			Double_t sumXi0 = 0;
			Double_t sumErrXi0 = 0;
			Double_t sumXiErr = 0;
			Double_t sumErrXiErr = 0;
			Double_t pT = hV0PtFeeddown[iSp][iMu][iSph]->GetYaxis()->GetBinCenter(iBin);

			for (int iMass = 1; iMass < hV0FeeddownMatrix[iSp]->GetNbinsZ()+1; ++iMass)	{

				Double_t sign = 1.;
				Double_t mass = hV0FeeddownMatrix[iSp]->GetZaxis()->GetBinCenter(iMass);
				if (mass > mParMuL->Eval(pT) + nSig * mParSigL->Eval(pT)) sign = -1.;
				if (mass < mParMuL->Eval(pT) - nSig * mParSigL->Eval(pT)) sign = -1.; 

				//if (iBin == 5 && iMass%25==0) cout << "sign " << sign << " mass " << mass << " interval " << mParMuL->Eval(pT) + nSig * mParSigL->Eval(pT) << endl;

				for (int iMotherBin = 1; iMotherBin < hV0FeeddownMatrix[iSp]->GetNbinsX()+1; ++iMotherBin)	{
				
					Double_t cX = hV0FeeddownMatrix[iSp]->GetXaxis()->GetBinCenter(iMotherBin);
					Double_t cY = hV0FeeddownMatrix[iSp]->GetYaxis()->GetBinCenter(iBin);
					//if (!cutg->IsInside(cX,cY)) continue;

					Double_t left 	= hV0FeeddownMatrix[iSp]->GetXaxis()->GetBinLowEdge(iMotherBin);
					Double_t right 	= hV0FeeddownMatrix[iSp]->GetXaxis()->GetBinLowEdge(iMotherBin+1);
					
					if (mFileXi) {
						sum += sign *
						hV0FeeddownMatrix[iSp]->GetBinContent(iMotherBin,iBin,iMass) *
						funcLT->Integral(left,right) / (xixibarFactor);//*(right-left)); //(2.*(right-left)); //2 because xi,xibar
						//hV0FeeddownMotherPt[iSp]->Integral(iMotherBin,iMotherBin)*(hV0FeeddownMotherPt[iSp]->GetBinLowEdge(iMotherBin+1) - hV0FeeddownMotherPt[iSp]->GetBinLowEdge(iMotherBin)) / (xixibarFactor);
						sumErr +=
						hV0FeeddownMatrix[iSp]->GetBinContent(iMotherBin,iBin,iMass) *
						//(hV0FeeddownMatrix[iSp]->GetBinError(iMotherBin,iBin) *
						//hV0FeeddownMatrix[iSp]->GetBinError(iMotherBin,iBin)) +
						(funcLT->IntegralError(left,right) *
						funcLT->IntegralError(left,right));

						sumXi0 += sign *
						hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) *
						funcLT->Integral(left,right) / (xixibarFactor);//*(right-left)); //(2.*(right-left)); //2 because xi,xibar
						//hV0FeeddownMotherPtXi0[iSp]->Integral(iMotherBin,iMotherBin)*(hV0FeeddownMotherPt[iSp]->GetBinLowEdge(iMotherBin+1) - hV0FeeddownMotherPt[iSp]->GetBinLowEdge(iMotherBin)) / (xixibarFactor);

						//if (iBin == 5 && iMass%25==0) cout << "mom pt " << left << " mass " << iMass << " adding " << sign *hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) *funcLT->Integral(left,right) / (xixibarFactor) << endl;
						//if (iBin == 5 && iMass%25==0) cout << "to total " << sumXi0 << endl;
						sumErrXi0 +=
						hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) *
						//(hV0FeeddownMatrixXi0[iSp]->GetBinError(iMotherBin,iBin) *
						//hV0FeeddownMatrixXi0[iSp]->GetBinError(iMotherBin,iBin)) +
						(funcLT->IntegralError(left,right) *
						funcLT->IntegralError(left,right));

						//cout << "AAAAAAAAAA " << funcLT->Integral(left,right) << " VS " << funcLT->IntegralError(left,right) << endl;

						sumXiErr += sign *
						hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) *
						( funcLThi->Integral(left,right) ) / (xixibarFactor);//*(right-left)); //(2.*(right-left)); //2 because xi,xibar
						//hV0FeeddownMotherPt[iSp]->Integral(iMotherBin,iMotherBin)*(hV0FeeddownMotherPt[iSp]->GetBinLowEdge(iMotherBin+1) - hV0FeeddownMotherPt[iSp]->GetBinLowEdge(iMotherBin)) / (xixibarFactor);
						sumErrXiErr +=
						hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) *
						//(hV0FeeddownMatrix[iSp]->GetBinError(iMotherBin,iBin) *
						//hV0FeeddownMatrix[iSp]->GetBinError(iMotherBin,iBin)) +
						(funcLT->IntegralError(left,right) *
						funcLT->IntegralError(left,right));

					} else {
						sum += sign *
						hV0FeeddownMatrix[iSp]->GetBinContent(iMotherBin,iBin,iMass) * 
						(right-left) * hXiPt[iSp][iMu][iSph]->Integral(iMotherBin,iMotherBin) / (xixibarFactor);
						//	funcLT->Integral(left,right) / (xixibarFactor);
						sumXi0 += sign *
						hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) * 
						//1.;
						//	hXiPt[iSp][iMu]->Integral(iMotherBin,iMotherBin) / (xixibarFactor);
						(right-left) *	hXiPt[iSp][iMu][iSph]->Integral(iMotherBin,iMotherBin) / (xixibarFactor);


						//if (iBin == 5 && iMass%25==0) cout << "mom pt " << left << " mass " << iMass << " adding " << sign * hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) * (right-left) *	hXiPt[iSp][iMu]->Integral(iMotherBin,iMotherBin) / (xixibarFactor) << endl;
						//if (iBin == 5 && iMass%25==0) cout << "to total " << sumXi0 << endl;

						sumXiErr += sign *
						hV0FeeddownMatrixXi0[iSp]->GetBinContent(iMotherBin,iBin,iMass) * 
						(right-left) *	hXiPt[iSp][iMu][iSph]->Integral(iMotherBin,iMotherBin) / (xixibarFactor);
						//funcLT->Integral(left,right) / (xixibarFactor);
						//sum +=
						//hV0FeeddownMatrix[iSp]->GetBinContent(iMotherBin,iBin) *
						//funcLT->Integral(left,right) / (xixibarFactor);//*(right-left)); //(2.*(right-left)); //2 because xi,xibar
						//printf("f %f vs h %f \n", funcLT->Integral(left,right), hXiPt[iSp][iMu]->Integral(iMotherBin,iMotherBin));
					}

				} 
			}

			hV0PtFeeddown[iSp][iMu][iSph]->SetBinContent(iBin, 2.*sum);	//2 also because xi0
			hV0PtFeeddown[iSp][iMu][iSph]->SetBinError(iBin, sumErr);


			hV0PtFeeddownXi0[iSp][iMu][iSph]->SetBinContent(iBin, 1.*sumXi0);	
			hV0PtFeeddownXi0[iSp][iMu][iSph]->SetBinError(iBin, sumErrXi0);


			hV0PtFeeddownXiErr[iSp][iMu][iSph]->SetBinContent(iBin, 1.*sumXiErr);	
			hV0PtFeeddownXiErr[iSp][iMu][iSph]->SetBinError(iBin, sumErrXiErr);
		}

		hV0PtFeeddown[iSp][iMu][iSph]->Scale(1,"width");
		hV0PtFeeddownXi0[iSp][iMu][iSph]->Scale(1,"width");
		hV0PtFeeddownXiErr[iSp][iMu][iSph]->Scale(1,"width");

		TH1F* htmp = (TH1F*)hV0PtFitCorr[iSp][0][iMu][iSph]->Clone("htmp");
		
		//for (int iSph = 0; iSph < NSPHERO; iSph++)  {
		//	if (iMu==0 && iSph>0) continue;
			cout << hV0PtFitCorr[iSp][0][iMu][iSph] << endl;
			hV0PtFitCorr[iSp][0][iMu][iSph]->Add(hV0PtFeeddownXi0[iSp][iMu][iSph],-1.);
		//}

		if (!iMu && mHandler->GetFlagMC() && !iSph) {
			//hClosureTestFDvSec[iSp] 	= (TH1F*)hV0PtFeeddown[iSp][iMu]->Clone(Form("hClosureTestFDvSec_%s",SPECIES[iSp]));
			//hClosureTestFDvSecPDG[iSp] 	= (TH1F*)hV0PtFeeddown[iSp][iMu]->Clone(Form("hClosureTestFDvSecPDG_%s",SPECIES[iSp]));
			//TH1F* htmpFD = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitSecondary_%s",SPECIES[iSp]));
			//TH1F* htmpFDPDG = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFitSecondaryPDG_%s",SPECIES[iSp]));
			hClosureTestFDvSec[iSp]->Divide(hV0PtFeeddown[iSp][iMu][iSph]);
			hClosureTestFDvSec[iSp]->Scale(2.0);
			hClosureTestFDvSecPDG[iSp]->Divide(hV0PtFeeddownXi0[iSp][iMu][iSph]);
			//hClosureTestFDvSec[iSp]->Divide(htmp);
			//hClosureTestFDvSecPDG[iSp]->Divide(htmp);
			
		}

		hV0PtFeeddownContr[iSp][iMu][iSph] = (TH1F*)hV0PtFeeddownXi0[iSp][iMu][iSph]->Clone(Form("hV0PtFeeddownContr_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		//hV0PtFeeddownContr[iSp][iMu]->Scale(1./1.6);

		hV0PtFeeddown[iSp][iMu][iSph]->Divide(htmp);

		hV0PtFeeddownXi0[iSp][iMu][iSph]->Divide(htmp);		
		hV0PtFeeddownXiErr[iSp][iMu][iSph]->Divide(htmp);		
		delete htmp;

	}	}	}


	funcLT->Write();

	cout << "bc " << hV0PtFeeddown[2][1][0]->GetBinContent(6) << " / " << hV0PtFitCorr[2][0][1][0]->GetBinContent(6) << endl;

	{
		TCanvas* cFDm = new TCanvas("cFDm","",1000,1000);
		cFDm->SetLogz();
		cFDm->SetRightMargin(0.13);
		hV0FeeddownMatrix[2]->SetStats(0);
		hV0FeeddownMatrix[2]->Project3D("xy")->Draw("colz");

		TLegend* legPt = new TLegend(0.50,.17,0.85,0.30);
		mHandler->MakeNiceLegend(legPt,0.04,1);	
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"#Xi^{-} #rightarrow #Lambda","");
		legPt->Draw();

		cFDm->SaveAs("plots/l_fdm.png");
	}
	{
		TCanvas* cFDm = new TCanvas("cFDm2","",1000,1000);
		cFDm->SetLogz();
		cFDm->SetRightMargin(0.13);
		hV0FeeddownMatrix[3]->SetStats(0);
		hV0FeeddownMatrix[3]->Project3D("xy")->Draw("colz");

		TLegend* legPt = new TLegend(0.50,.17,0.85,0.30);
		mHandler->MakeNiceLegend(legPt,0.04,1);	
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"#Xi^{+} #rightarrow #bar{#Lambda}","");
		legPt->Draw();

		cFDm->SaveAs("plots/lbar_fdm.png");
	}
	
	for (Int_t iSp=2; iSp < 4; iSp++)	{
		
		TCanvas* cFDmb = new TCanvas(Form("cFDmb_%s",SPECIES[iSp]),"",1000,1200);
		cFDmb->SetGridy();cFDmb->SetGridx();
		mHandler->MakeNiceHistogram(hV0PtFeeddown[iSp][0][0],kBlack); hV0PtFeeddown[iSp][0][0]->SetMarkerStyle(20);
		mHandler->MakeNiceHistogram(hV0PtFeeddownXi0[iSp][0][0],kBlack); hV0PtFeeddownXi0[iSp][0][0]->SetMarkerStyle(24);
		hV0PtFeeddown[iSp][0][0]->GetYaxis()->SetRangeUser(0.,0.4);
		hV0PtFeeddown[iSp][0][0]->GetYaxis()->SetTitle("fraction removed");
		hV0PtFeeddown[iSp][0][0]->Draw();
		hV0PtFeeddownXi0[iSp][0][0]->Draw("same");

		TLegend* legPt = new TLegend(0.52,.63,0.85,0.88);
		mHandler->MakeNiceLegend(legPt,0.04,1);
		
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"","");
		legPt->AddEntry((TObject*)0,Form("%s MB",SPECNAMES[iSp]),"");
		legPt->AddEntry(hV0PtFeeddown[iSp][0][0],"Double charged","pl");
		legPt->AddEntry(hV0PtFeeddownXi0[iSp][0][0],"MC ratio","pl");
		legPt->Draw();

		mHandler->MakeRatioPlot(hV0PtFeeddown[iSp][0][0],hV0PtFeeddownXi0[iSp][0][0],cFDmb,0.6,1.4,hV0PtFeeddownXi0[iSp][0][0]->GetBinLowEdge(1),hV0PtFeeddownXi0[iSp][0][0]->GetBinLowEdge(NPTBINS+1));

		cFDmb->Update();
		cFDmb->SaveAs(Form("plots/%s_fd_mb.png",SPECIES[iSp]));
	}
	{
		TCanvas* cFDhm = new TCanvas("cFDhm","",1000,1000);
		cFDhm->SetGridy();cFDhm->SetGridx();
		mHandler->MakeNiceHistogram(hV0PtFeeddown[2][1][0],kBlack); hV0PtFeeddown[2][1][0]->SetMarkerStyle(20);
		mHandler->MakeNiceHistogram(hV0PtFeeddown[3][1][0],kBlack); hV0PtFeeddown[3][1][0]->SetMarkerStyle(24);
		hV0PtFeeddown[2][1][0]->GetYaxis()->SetRangeUser(0.,0.4);
		hV0PtFeeddown[2][1][0]->GetYaxis()->SetTitle("fraction removed");
		hV0PtFeeddown[2][1][0]->Draw();
		hV0PtFeeddown[3][1][0]->Draw("same");

		TLegend* legPt = new TLegend(0.52,.63,0.85,0.88);
		mHandler->MakeNiceLegend(legPt,0.04,1);
		
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"","");
		legPt->AddEntry((TObject*)0,"V0M I-III","");
		legPt->AddEntry(hV0PtFeeddown[2][1][0],"#Lambda","pl");
		legPt->AddEntry(hV0PtFeeddown[3][1][0],"#bar{#Lambda}","pl");
		legPt->Draw();

		cFDhm->SaveAs("plots/l_fd_hm.png");
	}

}

void MyAnalysisV0correct::NormaliseSpectra() {

	/*for (int iBin = 0; iBin < NEVENTTYPES+2; ++iBin)
	{
		printf("Event type %i containing %f events\n", iBin, hEventType->GetBinContent(iBin));
	}*/

	NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
	//NormEta = (cuts::V0_Y[1] - cuts::V0_Y[0]);
	printf("Normalising all histograms by dY %f \n", NormEta);

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		if (iMu == 0 && iSph > 0) continue;	
		//if (iMu > 4 && (iSph < 3 && iSph)) continue;
		//if (iMu < 5 && iSph > 2) continue; 
		Double_t NormEv = hEventType->GetBinContent(6);	// MB
		if (iMu || iSph) {
			if (iType==1) NormEv = hEventMultvSpheroMC->GetBinContent(1+iSph,1+iMu);
			else NormEv = hEventMultvSpheroD->GetBinContent(1+iSph,1+iMu);

			//if (iMu == 3) NormEv = hEventMultvSpheroD->GetBinContent(1+iSph,1+1);
			//if (iMu == 4) NormEv = hEventMultvSpheroD->GetBinContent(1+iSph,1+2);
		}

		//if (iMu == 5 || iMu == 6 || iMu == 7)	{	NormEv = hEventType->GetBinContent(11);		// RT
		//	if (iSph == 3)	NormEv = hEventType->GetBinContent(12);		// RT 0-1
		//	if (iSph == 4)	NormEv = hEventType->GetBinContent(13);		// RT 1-2
		//	if (iSph == 5)	NormEv = hEventType->GetBinContent(14);		// RT 2-3
		//	if (iSph == 6)	NormEv = hEventType->GetBinContent(15);		// RT 3-4
		//	if (iSph == 7)	NormEv = hEventType->GetBinContent(16);	}	// RT 4-5
		//if (iType == 1) {
		//	if (iMu == 1) {
		//		if (iSph == 1 )	NormEv = hEventType->GetBinContent(18);		// FHM JET MC
		//		if (iSph == 2 )	NormEv = hEventType->GetBinContent(17);	}	// FHM ISO MC
		//	if (iMu == 2) {
		//		if (iSph == 1 )	NormEv = hEventType->GetBinContent(20);		// MHM JET MC
		//		if (iSph == 2 )	NormEv = hEventType->GetBinContent(19);	}	// MHM ISO MC
		//}
		
		if (NormEv>0) {
			if (iMu == 0) NormEv += NormEv * hEventType->GetBinContent(4) * 1./(hEventType->GetBinContent(5) + hEventType->GetBinContent(6));
		}

		if (NormEv == 0) NormEv = 1;

		printf("Normalising histogram %s by event count %f \n", hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName(), NormEv);
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEv);
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEta);

		if (mHandler->GetFlagMC() && !iMu && !iSph && !iType && iSp>1) {
			hClosureTestFDvSec[iSp]->Scale(1./NormEv);
			hClosureTestFDvSec[iSp]->Scale(1./NormEta);
			hClosureTestFDvSecPDG[iSp]->Scale(1./NormEv);
			hClosureTestFDvSecPDG[iSp]->Scale(1./NormEta);
			hClosureTestPDG[iSp]->Scale(1./NormEv);
			hClosureTestPDG[iSp]->Scale(1./NormEta);
		}
		if (mHandler->GetFlagMC() && !iMu && !iSph && !iType && iSp>0) {
			hClosureTestEffi[iSp]->Scale(1./NormEv);
			hClosureTestEffi[iSp]->Scale(1./NormEta);
			hClosureTestEffiPDG[iSp]->Scale(1./NormEv);
			hClosureTestEffiPDG[iSp]->Scale(1./NormEta);
		}
		
	} } } }

	printf("mb k0s spectrum from fit \n");
	cout << hV0PtFit[1][0][0][0]->GetBinContent(30) << endl;
	printf("mb k0s spectrum from norm \n");
	cout << hV0PtFitCorr[1][0][0][0]->GetBinContent(30) << endl;


	Double_t ntedges[NRTBINS0+1];
	for (int i=0; i<NRTBINS0+1; i++) ntedges[i] = RTBINS0[i]*RT_DEN;
	Int_t binedges[NRTBINS0+1];
	for (int i=0; i<NRTBINS0+1; i++) binedges[i] = hNchTrans->FindBin(ntedges[i]);
	for (int i=0; i<NRTBINS0+1; i++) cout << ntedges[i] << " -- " << binedges[i] << endl;


	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{

		printf("Normalising histogram %s by N_T distribution \n", hV0PtNtFitCorr[iSp][iType][iReg]->GetName());

		for (int iNt = 1; iNt < hV0PtNtFitCorr[1][0][0]->GetNbinsY()+1; iNt++) {

			Double_t NormEv = 0;
			//NormEv = hNchTrans->GetBinContent(iNt);
			NormEv = (iType==2) ? hNchTransMCTrigMC->Integral(1,50) : hNchTrans->Integral(1,50);
			TH2F* htmp = (iType==2) ? hV0PtNt[iSp][iType][iReg] : hV0PtNtFitCorr[iSp][iType][iReg];

			TH1F* hpt = (TH1F*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

		}

	hV0PtNt[iSp][iType][iReg]->Write();

	}	}	}

	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{

		printf("Normalising histograms pi and Kpm in region %s yn", REGIONS[iReg]);

		for (int iNt = 1; iNt < hPiPtNtRC[0]->GetNbinsY()+1; iNt++) {

			Double_t NormEv = 0;
			NormEv = hNchTrans->Integral(1,50);

			TH1F* hpt = (TH1F*)hPiPtNtRC[iReg]->ProjectionX(Form("%s_%i",hPiPtNtRC[iReg]->GetName(),iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hPiPtNtRC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hPiPtNtRC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

			/*hpt = (TH1F*)hPiPtNtMinRC[iReg]->ProjectionX("",iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hPiPtNtMinRC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hPiPtNtMinRC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}*/

			hpt = (TH1F*)hKpmPtNtRC[iReg]->ProjectionX(Form("%s_%i",hKpmPtNtRC[iReg]->GetName(),iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hKpmPtNtRC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hKpmPtNtRC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

			/*hpt = (TH1F*)hKpmPtNtMinRC[iReg]->ProjectionX("",iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hKpmPtNtMinRC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hKpmPtNtMinRC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}*/

			NormEv = hNchTransMCTrigMC->Integral(1,50);

			hpt = (TH1F*)hPiPtNtMC[iReg]->ProjectionX(Form("%s_%i",hPiPtNtMC[iReg]->GetName(),iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hPiPtNtMC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hPiPtNtMC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

			/*hpt = (TH1F*)hPiPtNtMinMC[iReg]->ProjectionX("",iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hPiPtNtMinMC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hPiPtNtMinMC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}*/

			hpt = (TH1F*)hKpmPtNtMC[iReg]->ProjectionX(Form("%s_%i",hKpmPtNtMC[iReg]->GetName(),iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hKpmPtNtMC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hKpmPtNtMC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

			/*hpt = (TH1F*)hKpmPtNtMinMC[iReg]->ProjectionX("",iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < hPiPtNtRC[0]->GetNbinsX()+1; iBin++) {
				hKpmPtNtMinMC[iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hKpmPtNtMinMC[iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}*/

		}		

		hPiPtNtRC[iReg]->Write();
		hPiPtNtMC[iReg]->Write();
		hKpmPtNtRC[iReg]->Write();
		hKpmPtNtMC[iReg]->Write();

	}


	//hNchTrans->Write();

	// normalise pt spectra (rt)
	/*
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{

		Double_t NormEv = 0;
		NormEv = hNchTrans->Integral(binedges[iRtBin],binedges[iRtBin+1]);
		//cout << "bla " << RTBINS0[iRtBin] << " " << RTBINS0[iRtBin]*rt_den <<
		//" bla " << RTBINS0[iRtBin+1] << " " << RTBINS0[iRtBin+1]*rt_den << endl;
		//for (int ib = 1; ib < hNchTrans->GetNbinsX(); ++ib)
		//{
			//cout << "aaaa " << 
		//	if (hNchTrans->GetBinCenter(ib) < RTBINS0[iRtBin]*rt_den-1E-05) continue;
		//	if (hNchTrans->GetBinCenter(ib) > RTBINS0[iRtBin+1]*rt_den) continue;
		//	NormEv += hNchTrans->GetBinContent(ib);
			//cout << "ib " << ib << " bins " << hNchTrans->GetBinCenter(ib) << " sum " << sum << endl;
		//}

		//if (NormEv>0) {
		//	NormEv += NormEv * hEventType->GetBinContent(4) * 1./(hEventType->GetBinContent(5) + hEventType->GetBinContent(6));
		//}

		printf("Normalising histogram %s by event count %f \n", hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->GetName(), NormEv);
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Scale(1./NormEv);
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Scale(1./NormEta);


		//Int_t leftbin = hNchTrans->FindBin(RTBINS0[iRtBin]*rt_den);
		//Int_t rightbin = hNchTrans->FindBin(RTBINS0[iRtBin+1]*rt_den);
		//Double_t NormEv	= hNchTrans->Integral(leftbin,rightbin-1);

		//cout << iSp << " " << iType << " " << iReg << " " << iRtBin << endl;
		//cout << "bla " << RTBINS0[iRtBin] << " " << RTBINS0[iRtBin]*rt_den <<
		//" bla " << RTBINS0[iRtBin+1] << " " << RTBINS0[iRtBin+1]*rt_den << " " << leftbin << " " << rightbin << endl;
		//for (int ib = 0; ib < rightbin-leftbin+1; ++ib)
		//{
		//	cout << "irtbin " << iRtBin << " lb " << leftbin << " rb " << rightbin << " bins " << hNchTrans->GetXaxis()->GetBinCenter(leftbin+ib) << endl;
		//}

		//Int_t leftbin = hRt2->FindBin(RTBINS0[iRtBin]);
		//Int_t rightbin = hRt2->FindBin(RTBINS0[iRtBin+1]);
		//Double_t NormEv	= hRt2->Integral(leftbin,rightbin-1);

	}	}	}	}
	*/
	
}

void MyAnalysisV0correct::LoadEfficiency() {

	if (!mFileMC) {
		printf("No MC file loaded! Efficiency correction not performed.\n");
		return;
	}

	TDirectoryFile* dirFile1 = new TDirectoryFile("mcFile","mcFile","",mHandler->file());
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("TDirectoryFile")) {
		cout << "Doing efficiency from a TDirectoryFile" << endl;
		dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");}
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("THashList")) {
		cout << "Doing efficiency from a THashList" << endl;
		THashList* hashList = (THashList*)mFileMC->Get("MyAnalysisV0_0");
		while (hashList->GetEntries()) {
			dirFile1->Append(hashList->First());
			hashList->RemoveFirst();
		}
	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	
		hV0Efficiency[iSp] = (TH1D*)dirFile1->Get(Form("hV0Efficiency_%s",SPECIES[iSp]));
		hV0Efficiency[iSp]->Write();
	}

}

void MyAnalysisV0correct::DoEfficiencyFromFile() {

	if (!mFileMC) {
		printf("No MC file loaded! Efficiency correction not performed.\n");
		return;
	}

	TDirectoryFile* dirFile1 = new TDirectoryFile("mcFile","mcFile","",mHandler->file());
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("TDirectoryFile")) {
		cout << "Doing efficiency from a TDirectoryFile" << endl;
		dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");}
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("THashList")) {
		cout << "Doing efficiency from a THashList" << endl;
		THashList* hashList = (THashList*)mFileMC->Get("MyAnalysisV0_0");
		while (hashList->GetEntries()) {
			dirFile1->Append(hashList->First());
			hashList->RemoveFirst();
		}
	}

	//TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < 2; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{		
	
		tV0massRt[iSp][iType][iReg] = (TNtuple*)dirFile1->Get(Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));

	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
		tV0PtMCMB[iSp]		= (TNtuple*)dirFile1->Get(Form("tV0PtMCMB_%s",SPECIES[iSp]));
		tV0massRCMB[iSp]	= (TNtuple*)dirFile1->Get(Form("tV0massRCMB_%s",SPECIES[iSp]));
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{	
			tV0PtMCRt[iSp][iReg]	= (TNtuple*)dirFile1->Get(Form("tV0PtMCMB_%s_%s",SPECIES[iSp],REGIONS[iReg]));
		}
	}

	
	DoEfficiencyFromTrees();

}

void MyAnalysisV0correct::DoEfficiencyFromTrees() {

	Float_t nSig = 4.;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		
		// minimum bias histos
		hV0Efficiency[iSp] = new TH1D(Form("hV0Efficiency_%s",SPECIES[iSp]),"; V0 pT (GeV/#it{c}); Efficiency",NPTBINS,XBINS); //NPTBINS2, XBINS2);
		TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

		//TString masscut = (iSp == 1 ) ? Form("MassDT < 0.03 && MassDT > -0.03") : Form("MassDT < 0.0075 && MassDT > -0.0075");
		TCut masscutL; TCut masscutR;
		TString mu; TString sig;
		if (iSp==1) {
			mu = Form("(lPt<=1.6)*(%f+%f*lPt+%f*lPt*lPt)+(lPt>1.6)*%f",cuts::K0S_PARMU[0],cuts::K0S_PARMU[1],cuts::K0S_PARMU[2],cuts::K0S_PARMU[3]);
			sig = Form("%f+%f*lPt+%f/lPt",cuts::K0S_PARSIG[0],cuts::K0S_PARSIG[1],cuts::K0S_PARSIG[2]);
		} else {
			mu = Form("(lPt<=1.9)*(%f+%f*lPt+%f*lPt*lPt)+(lPt>1.9)*(%f+%f*lPt)",cuts::L_PARMU[0],cuts::L_PARMU[1],cuts::L_PARMU[2],cuts::L_PARMU[3],cuts::L_PARMU[4]);
			sig = Form("%f+%f*lPt+%f/lPt",cuts::L_PARSIG[0],cuts::L_PARSIG[1],cuts::L_PARSIG[2]);
		}
		masscutR = Form("MassDT < %s + %f*(%s)", mu.Data(), nSig, sig.Data());
		masscutL = Form("MassDT > %s - %f*(%s)", mu.Data(), nSig, sig.Data());	

		tV0massRCMB[iSp]->Draw(Form("lPt>>hV0Efficiency_%s",SPECIES[iSp]),masscutL+masscutR,"goff");
		tV0PtMCMB[iSp]->Draw("lPt>>hDen","","goff");

		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,20.0);

		hV0Efficiency[iSp]->Divide(hV0Efficiency[iSp],hDen,1.,1.,"B");
		delete hDen;

	}
	

	{
		TCanvas* cEffi = new TCanvas("cEffi","",2700,900);
		cEffi->Divide(3,1,0.0005,0.0005);
		for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
			cEffi->cd(iSp);
			mHandler->MakeNiceHistogram(hV0Efficiency[iSp],COLOURS[iSp]);
			hV0Efficiency[iSp]->Draw("same");
		}
		TLegend* leg1 = new TLegend(0.19,0.13,0.9,0.25);//cFits[canCounter/NPTBINS]->BuildLegend();
		mHandler->MakeNiceLegend(leg1, 0.065, 3);
		for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
			leg1->AddEntry(hV0Efficiency[iSp],Form("%s",SPECNAMES[iSp]),"pl"); }
		leg1->Draw();
		cEffi->Write();
	}
	

	/*{
		TCanvas* cEffi[NSPECIES];
		for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
			cEffi[iSp] = new TCanvas(Form("cEffiReg_%s",SPECIES[iSp]),"",1000,900);
			
			mHandler->MakeNiceHistogram(hV0Efficiency[iSp],kBlue);
			hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(-0.01,0.8);
			hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.0,5.0);
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
			leg1->AddEntry(hV0EfficiencyRt[iSp][0],"R_{T} Trans.","pl");
			leg1->AddEntry(hV0EfficiencyRt[iSp][1],"R_{T} Near","pl");
			leg1->AddEntry(hV0EfficiencyRt[iSp][2],"R_{T} Away","pl");
			
			cEffi[iSp]->Update();
			leg1->Draw();

			cEffi[iSp]->Write();
			cEffi[iSp]->SaveAs(Form("plots/cEffiReg_%s.png",SPECIES[iSp]));
		}
	}*/

}


void MyAnalysisV0correct::CorrectSpectra() {

	Double_t MBtrigEff = 0.7448;
	

	printf("mb k0s spectrum before corr \n");
	//cout << hV0PtFitCorr[2][0][1][0]->GetBinContent(30) << endl;

	TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		
		hV0PtNoTrigger[iSp]
			= (TH1F*)dirFile1->Get(Form("hV0PtNoTrigger_%s",SPECIES[iSp]));
		
		if (hV0PtNoTrigger[iSp]->GetNbinsX() != NPTBINS)	{
			hV0PtNoTrigger[iSp] = (TH1F*)hV0PtNoTrigger[iSp]->Rebin(NPTBINS,hV0PtNoTrigger[iSp]->GetName(),XBINS);
		}

		hV0PtSignalLoss[iSp] 
			= (TH1F*)dirFile1->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[2],MULTI[0],SPHERO[0]));
		if (hV0PtSignalLoss[iSp]->GetNbinsX() != NPTBINS) 
			hV0PtSignalLoss[iSp] = (TH1F*)hV0PtSignalLoss[iSp]->Rebin(NPTBINS,hV0PtSignalLoss[iSp]->GetName(),XBINS);
	
		hV0PtSignalLoss[iSp]->Divide(hV0PtNoTrigger[iSp]);
	}
	mDirFile->cd();

	


	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	TF1* funcRapCorrection = new TF1("funcRapCorrection",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	//for (int iMu = 0; iMu < 3; ++iMu)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{


		if (iMu == 0 && iSph > 0) continue;	
		//if (iMu > 4 && (iSph < 3 && iSph)) continue;
		//if (iMu < 5 && iSph > 2) continue; 
		//hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1,"width");
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Divide(hV0Efficiency[iSp]);

		funcRapCorrection->SetParameters(0.5*NormEta,MASSES[iSp]);
		
		//cout << "at 1.0 correcting by " << funcRapCorrection->Eval(1.0) << endl;

		hV0PtFitCorr[iSp][iType][iMu][iSph]->Divide(funcRapCorrection,1.);

		if (iMu == 0) {
			hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(MBtrigEff);
			//hV0PtFitCorr[iSp][iType][iMu][iSph]->Divide(hV0PtSignalLoss[iSp]);
		}



	}	}	}	}
	funcRapCorrection->Write();


	printf("mb k0s eff \n");
	//cout << hV0Efficiency[2]->GetBinContent(30) << endl;

	printf("mb k0s spectrum after corr \n");
	//cout << hV0PtFitCorr[2][0][1][0]->GetBinContent(30) << endl;

	TF1* funcRapCorrection2 = new TF1("funcRapCorrection2",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
		funcRapCorrection2->SetParameters(0.5*NormEta,MASSES[iSp]);
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{


		for (int iNt = 1; iNt < hV0PtNtFitCorr[1][0][0]->GetNbinsY()+1; iNt++) {

			TH1F* hpt = (TH1F*)hV0PtNtFitCorr[iSp][iType][iReg]->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			hpt->Divide(hV0Efficiency[iSp]);
			//hpt->Divide(funcRapCorrection2,1.);
			//hpt->Scale(MBtrigEff);

			for (int iBin = 1; iBin < hV0PtNtFitCorr[1][0][0]->GetNbinsX()+1; iBin++) {
				hV0PtNtFitCorr[iSp][iType][iReg]->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				hV0PtNtFitCorr[iSp][iType][iReg]->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

		}

	}	}	}


	/*TF1* funcRapCorrection2 = new TF1("funcRapCorrection2",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Divide(hV0Efficiency[iSp]);

		funcRapCorrection2->SetParameters(0.5*NormEta,MASSES[iSp]);
		
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Divide(funcRapCorrection2,1.);

		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Scale(MBtrigEff);
	}	}	}	}*/

	{
		TCanvas* cSL = new TCanvas("cSL","",1000,1000);
		cSL->SetGridy();cSL->SetGridx();
		mHandler->MakeNiceHistogram(hV0PtSignalLoss[1],kBlack);
		mHandler->MakeNiceHistogram(hV0PtSignalLoss[2],kBlue);
		mHandler->MakeNiceHistogram(hV0PtSignalLoss[3],kGreen+2);
		
		hV0PtSignalLoss[1]->GetYaxis()->SetRangeUser(0.95,1.05);
		hV0PtSignalLoss[1]->GetYaxis()->SetTitle("signal loss due to trigger");
		hV0PtSignalLoss[1]->Draw();
		hV0PtSignalLoss[2]->Draw("same");
		hV0PtSignalLoss[3]->Draw("same");

		TLegend* legPt = new TLegend(0.52,.63,0.85,0.88);
		mHandler->MakeNiceLegend(legPt,0.04,1);
		
		legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
		legPt->AddEntry((TObject*)0,"","");
		legPt->AddEntry(hV0PtSignalLoss[1],"K^{0}_{S}","pl");
		legPt->AddEntry(hV0PtSignalLoss[2],"#Lambda","pl");
		legPt->AddEntry(hV0PtSignalLoss[3],"#bar{#Lambda}","pl");
		legPt->Draw();

		cSL->SaveAs("plots/sl.png");
	}



}

void MyAnalysisV0correct::StudyCuts() {

	for (int iCut = 0; iCut < 25; ++iCut)	{
		hV0PtCut[iCut]
			= (TH1F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0PtCut_%i",iCut) );
		
	}

	if (!mFileMC) {
		printf("No MC file loaded!.\n");
		return;
		} else {
		TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");
		for (int iCut = 0; iCut < 25; ++iCut)	{
			TH1F* htmp
				= (TH1F*)dirFile1->Get(Form("hV0PtCut_%i",iCut) );
			hV0PtCutMC[iCut] = (TH1F*)htmp->Clone(Form("hV0PtCutMC_%i",iCut));
			
			}
	}

	TCanvas* cCuts = new TCanvas("cCuts","",2800,2000);
	cCuts->Divide(6,5,0.00005,0.00005);

	for (int iCut = 24; iCut > 0; iCut--)	{
		//cout << "iC " << iCut << " " << hV0PtCut[iCut] << " " << hV0PtCut[iCut-1] << " " << hV0PtCutMC[iCut] << " " << hV0PtCutMC[iCut-1] << endl;
		hV0PtCut[iCut]->Divide(hV0PtCut[iCut-1]);
		hV0PtCutMC[iCut]->Divide(hV0PtCutMC[iCut-1]);
		mHandler->MakeNiceHistogram(hV0PtCut[iCut],kRed);
		mHandler->MakeNiceHistogram(hV0PtCutMC[iCut],kBlue);
		hV0PtCut[iCut]->SetMarkerSize(0.6);
		hV0PtCutMC[iCut]->SetMarkerSize(0.6);
		cCuts->cd(iCut);
		hV0PtCut[iCut]->Draw();
		hV0PtCutMC[iCut]->Draw("same");

		TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
			mHandler->MakeNiceLegend(leg1, 0.10, 1.);
			leg1->AddEntry((TObject*)0,Form("iCut %i / iCut %i",iCut, iCut-1)," ");
			leg1->AddEntry(hV0PtCut[24],Form("DATA"),"pl");
			leg1->AddEntry(hV0PtCutMC[24],Form("MC rec"),"pl");
			leg1->Draw();
	}

	cCuts->Write();

	for (int iCut = 0; iCut < 25; ++iCut)	{
		hV0PtCut[iCut]->Write();
		hV0PtCutMC[iCut]->Write();
	}

}

/*Double_t MyAnalysisV0correct::rap_correction(Double_t* x, Double_t* par)
{
  Double_t pt = x[0];
  Double_t eta  = par[0];
  Double_t mass = par[1];
  const Double_t mt = TMath::Sqrt(pt*pt + mass*mass);
  const Double_t rap = TMath::ASinH(pt/mt*TMath::SinH(eta));
  //  return rap/eta;
  return rap/eta;
} */

void MyAnalysisV0correct::DoClosureTest(Int_t opt) {

	// DIVIDING BLINDLY REC. MC DATA (CORRECTED) W MC GENERATED PARTICLES

	mHandler->root()->SetBatch(kTRUE);
	Double_t NormEv = hEventType->GetBinContent(6);	// MB
	
	if (NormEv>0) {
			NormEv += NormEv * hEventType->GetBinContent(4) * 1./(hEventType->GetBinContent(5) + hEventType->GetBinContent(6));
		}
	
	Double_t MBtrigEff = 0.7448;
	TF1* funcRapCorrection = new TF1("funcRapCorrection",rap_correction,XBINS[0],XBINS[NPTBINS],2);

	for (Int_t iSp = 1; iSp < NSPECIES; iSp++)		{
		Int_t iMu = 0; Int_t iSph = 0;	

		TH1D* hDen = (TH1D*)hV0Pt[iSp][2][iMu][iSph]->Clone(Form("hDen"));
		hDen->Scale(1./NormEv);
		hDen->Scale(1./NormEta);
		hDen->Scale(1.,"width");

		hClosureTestEffi[iSp]->Divide(hV0Efficiency[iSp]);
		hClosureTestEffiPDG[iSp]->Divide(hV0Efficiency[iSp]);
		hClosureTestEffi[iSp]->Divide(hDen);
		hClosureTestEffiPDG[iSp]->Divide(hDen);

		if (iSp>1) {
			TH1F* hcltmp = (TH1F*)hClosureTestPDG[iSp]->Clone("hcltmp");
		hcltmp->Multiply(hV0PtFeeddown[iSp][0][0]);
		hClosureTestPDG[iSp]->Add(hV0PtFeeddownContr[iSp][0][0],-1.);
		delete hcltmp;
		hClosureTestPDG[iSp]->Divide(hV0Efficiency[iSp]);
		hClosureTestPDG[iSp]->Divide(hDen);
		}

		hClosureTestCorr[iSp]	= (TH1F*)hV0PtFitCorr[iSp][0][iMu][iSph]->Clone(Form("hClosureTestCorr_%s",SPECIES[iSp]));
		hDen->Scale(MBtrigEff);
		funcRapCorrection->SetParameters(0.5*NormEta,MASSES[iSp]);
		hDen->Divide(funcRapCorrection,1.);
		mHandler->MakeNiceHistogram(hClosureTestCorr[iSp],kBlack);
		hClosureTestCorr[iSp]->GetYaxis()->SetTitle("blind rec. / MC generated");
		hClosureTestCorr[iSp]->GetYaxis()->SetRangeUser(0.7,1.3);
		hClosureTestCorr[iSp]->Divide(hDen);
		TCanvas* cClosureCorr = new TCanvas("cClosureCorr","",900,900);
		if (iSp == 3) hClosureTestCorr[3]->SetMarkerStyle(24);
		hClosureTestCorr[iSp]->Draw();
		if (iSp == 3) hClosureTestCorr[2]->Draw("same");
		TLegend* leg1 = new TLegend(0.171,0.67,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
		mHandler->MakeNiceLegend(leg1, 0.05, 1.);
		if (iSp == 3) leg1->AddEntry(hClosureTestCorr[2],Form("%s",SPECNAMES[2]),"pl");
		leg1->AddEntry(hClosureTestCorr[iSp],Form("%s",SPECNAMES[iSp]),"pl");
		leg1->Draw();
		cClosureCorr->SaveAs(Form("plots/closureCorr_%s.png",SPECIES[iSp]));
		delete hDen;
		
	}



	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0correct::DoXCheckV0M() {

	
	//mHandler->root()->SetBatch(kTRUE);
	
	TFile* fileV0M = new TFile("../official/SpectraVsMultiplicityLambdaSum.root","READ");
	TH1F* hOffiL01 = (TH1F*)fileV0M->Get("hPtLambdaSumStatOnly_V0M_00000to00100-epsPart-epsEv-Corrected");
	TH1F* hOffiL15 = (TH1F*)fileV0M->Get("hPtLambdaSumStatOnly_V0M_00100to00500-epsPart-epsEv-Corrected"); 
	TH1F* hOffiL = (TH1F*)fileV0M->Get("hPtLambdaSumStatOnly_V0M_00500to01000-epsPart-epsEv-Corrected"); 
	TH1F* hOffiLMB = (TH1F*)fileV0M->Get("hPtLambdaSumStatOnly_V0M_00000to10000-epsPart-epsEv-Corrected"); 
	
	TFile* fileMB 	= new TFile("../official/k0s_spectra.root","READ");
	TH1F* hOffiKMB 	= (TH1F*)fileMB->Get("fHistStatErrOnly"); 
	cout << "blaaadawda " << hOffiKMB << endl;

	TFile* fileMBeff 	= new TFile("../official/k0s_effi.root","READ");
	TH1F* hOffiKMBeff 	= (TH1F*)fileMBeff->Get("fHistPureEfficiency"); 
	cout << "blaaadawdaaa " << hOffiKMBeff << endl;

	TFile* fileMBraw 	= new TFile("../official/k0s_raw.root","READ");
	TH1F* hOffiKMBraw 	= (TH1F*)fileMBraw->Get("fHistPtK0ShortRaw"); 
	

	mDirFile->cd();

	hOffiL->Scale(5.); hOffiL->Add(hOffiL15,4.); hOffiL->Add(hOffiL01,1.);
	hOffiL->Scale(1./10);

	TH1F* hLLbarMB = (TH1F*)hV0PtFitCorr[3][0][0][0]->Clone("hLLbarMB");
	hLLbarMB->Add(hV0PtFitCorr[2][0][0][0]);
	hLLbarMB->Scale(1./0.7448);
	//hLLbarMB->Scale(2.);
	TH1F* hLLbarV0M = (TH1F*)hV0PtFitCorr[3][0][1][0]->Clone("hLLbarV0M"); 
	hLLbarV0M->Add(hV0PtFitCorr[2][0][1][0]);
	//hLLbarV0M->Scale(2.);
	//hLLbarV0M->Scale(1./0.7448);

	TH1F* hLLbarV0M01 = (TH1F*)hV0PtFitCorr[3][0][3][0]->Clone("hLLbarV0M01"); 
	hLLbarV0M01->Add(hV0PtFitCorr[2][0][3][0]);

	TF1* funcRapCorrectionL = new TF1("funcRapCorrectionL",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	funcRapCorrectionL->SetParameters(0.8,MASSES[2]);
	//hLLbarMB->Divide(funcRapCorrectionL,1);
	//hLLbarV0M->Divide(funcRapCorrectionL,1);


	//hLLbarV0M->Divide(hOffiL);
	//hLLbarMB->Divide(hOffiLMB);

cout << "blaaadawdaaa " << endl;
	
	if (1) {

		{

	mHandler->MakeNiceHistogram(hLLbarMB,kBlack);
	mHandler->MakeNiceHistogram(hLLbarV0M,kRed);
	mHandler->MakeNiceHistogram(hLLbarV0M01,kBlue);
	mHandler->MakeNiceHistogram((TH1F*)hOffiLMB,kGray+3);
	mHandler->MakeNiceHistogram((TH1F*)hOffiL,kRed+1);
	mHandler->MakeNiceHistogram((TH1F*)hOffiL01,kBlue+1);
	hOffiLMB->SetMarkerStyle(21); hOffiL->SetMarkerStyle(21); hOffiL01->SetMarkerStyle(21);

	
	//hMB->Divide(hOffiMB); hV0M1->Divide(hOffiV0M1); hV0M10->Divide(hOffiV0M10);
	TCanvas* cOffiL = new TCanvas("cOffiL","",1000,1000);
	//hMB->GetYaxis()->SetTitle("2K^{0}_{S} / K^{#pm}");
	hLLbarMB->GetYaxis()->SetRangeUser(1e-5,3.001);
	hLLbarMB->GetXaxis()->SetRangeUser(0.2,10.);
	cOffiL->SetGridx(); cOffiL->SetGridy(); cOffiL->SetLogy();
	hLLbarMB->Draw();
	hLLbarV0M->Draw("same");
	hLLbarV0M01->Draw("same");
	hOffiLMB->Draw("same");
	hOffiL->Draw("same");
	hOffiL01->Draw("same");

	cOffiL->cd();
	TLegend* legPt = new TLegend(0.39,0.63,0.88,0.88);
			mHandler->MakeNiceLegend(legPt,0.04,2);
			
			legPt->AddEntry((TObject*)0,Form("|#eta| < 0.8"),"");
			legPt->AddEntry((TObject*)0,"pp #sqrt{s} = 13 TeV","");
			legPt->AddEntry((TObject*)0,"","");
			legPt->AddEntry((TObject*)0,"#Lambda + #bar{#Lambda}","");
			legPt->AddEntry(hLLbarMB,"MB (V0M 0-100%)","pl");
			legPt->AddEntry(hOffiLMB,"LHC15f (offi.)","pl");
			legPt->AddEntry(hLLbarV0M01,"V0M I","pl");
			legPt->AddEntry(hOffiL01,"LHC15f (offi.)","pl");
			legPt->AddEntry(hLLbarV0M,"V0M I-III","pl");
			legPt->AddEntry(hOffiL,"LHC15f (offi.)","pl");
			legPt->Draw();


	mDirFile->cd();
	mHandler->MakeRatioPlot(hLLbarMB,(TH1F*)hOffiLMB,cOffiL,0.75,1.25,0.2,10.);
	mHandler->MakeRatioPlot(hLLbarV0M,(TH1F*)hOffiL,cOffiL,0.75,1.25,0.2,10.);
	mHandler->MakeRatioPlot(hLLbarV0M01,(TH1F*)hOffiL01,cOffiL,0.75,1.25,0.2,10.);

	/*mHandler->MakeZoomPlot(hMB,cOffiK0,0.2,1.501,0.599,1.401);
	hMB->DrawCopy();
	hV0M1->DrawCopy("same");
	hV0M10->DrawCopy("same");

	cOffiK0->cd();
	hMB->GetYaxis()->SetRangeUser(0.399,1.901);
	hMB->GetXaxis()->SetRangeUser(0.2,20.);*/

	cOffiL->SaveAs("plots/LtoL_MBandV0M_offi.png");
	cOffiL->SaveAs("plots/LtoL_MBandV0M_offi.pdf");

	}



		TCanvas* cXcheck = new TCanvas("cXcheck","",900,900);
		cXcheck->SetLogy();
		mHandler->MakeNiceHistogram((TH1F*)hOffiL,kGreen+2);
		mHandler->MakeNiceHistogram((TH1F*)hOffiLMB,kGreen+2);
		hOffiLMB->SetMarkerStyle(21);
		mHandler->MakeNiceHistogram(hLLbarV0M,kRed);
		mHandler->MakeNiceHistogram(hLLbarMB,kBlue);
		hLLbarMB->SetMarkerStyle(21);
		hOffiL->Draw();
		hLLbarV0M->Draw("same");
		hLLbarMB->Draw("same");
		hOffiLMB->Draw("same");
		mHandler->MakeRatioPlot(hLLbarMB,(TH1F*)hOffiLMB,cXcheck,0.6,1.4,0.4,8.);
		mHandler->MakeRatioPlot(hLLbarV0M,(TH1F*)hOffiL,cXcheck,0.6,1.4,0.4,8.);
	}
	cout << "blaaadawdaaa " << endl;

	TH1F* hK0sMB = (TH1F*)hV0PtFitCorr[1][0][0][0]->Clone("hK0sMB");
	TH1F* hK0sMBeff = (TH1F*)hV0Efficiency[1]->Clone("hK0sMBeff");
	TH1F* hK0sMBraw = (TH1F*)hV0PtFit[1][0][0][0]->Clone("hK0sMBraw");
	hK0sMBraw->Scale(1./NormEta);

	cout << "blaaadawdaaa " << endl;

	//Double_t NormEv = hEventType->GetBinContent(24);
	//NormEv += NormEv * hEventType->GetBinContent(22) * 1./(hEventType->GetBinContent(23) + hEventType->GetBinContent(24));
	//hK0sMBraw->Scale(1./NormEv);

	cout << "badawda " << hK0sMB << endl;
	cout << "badawda222 " << hK0sMBeff << endl;


	printf("atghh mb k0s eff \n");
	cout << hV0Efficiency[1]->GetBinContent(30) << endl;
	cout << hK0sMBeff->GetBinContent(30) << endl;

	TF1* funcRapCorrection3 = new TF1("funcRapCorrection3",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	funcRapCorrection3->SetParameters(0.8,MASSES[1]);
	//hOffiKMBeff->Divide(funcRapCorrection3,1);

	//hK0sMB->Divide(hOffiKMB);
	//hOffiKMBeff->Divide(hK0sMBeff);
	//hK0sMBraw->Divide(hOffiKMBraw);


	TCanvas* cXcheck2 = new TCanvas("cXcheck2","",900,900);
	cXcheck2->SetGridy();
	cXcheck2->SetLogy();
	//hOffiLMB->SetMarkerStyle(21);
	mHandler->MakeNiceHistogram(hK0sMB,kRed);
	mHandler->MakeNiceHistogram((TH1F*)hOffiKMB,kRed+1);
	hOffiKMB->SetMarkerStyle(21);
	//mHandler->MakeNiceHistogram(hK0sMBeff,kBlue);
	//mHandler->MakeNiceHistogram(hK0sMBraw,kBlue);
	//hLLbarMB->SetMarkerStyle(21);
	//hK0sMB->GetYaxis()->SetRangeUser(0.,2.);
	hK0sMB->GetYaxis()->SetTitle("this analysis / official MB analysis");
	hK0sMB->GetXaxis()->SetRangeUser(0.2,10.);
	
	hK0sMB->Draw();
	hOffiKMB->Draw("same");

	mHandler->MakeRatioPlotInterp(hK0sMB,(TH1F*)hOffiKMB,cXcheck2,0.799,1.201,0.2,10.);
	//hOffiKMB->Draw("same");
	
	//hK0sMBeff->Draw("same");
	//hOffiKMBeff->Draw("same");
	//hK0sMBraw->Draw("same");
	/*{TLegend *leg1 = new TLegend(0.45,0.72,0.85,0.85);
		mHandler->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry(hK0sMB,"corrected MB","pl");
		leg1->AddEntry(hOffiKMBeff,"1/efficiency","pl");
		leg1->AddEntry(hK0sMBraw,"raw yields","pl");
		leg1->Draw();
	}*/
	cXcheck2->SaveAs("plots/xcheckK0s.png");	
	//mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0correct::CreateOutputFile(const Char_t *name, Int_t Sp) {

	TString fileName = TString(name);
	TFile* outFile = new TFile(fileName,"RECREATE");

	for (Int_t iMu = 0; iMu < NMULTI; iMu++) {
	for (Int_t iSph = 0; iSph < NSPHERO; iSph++) {
		if (iMu == 0 && iSph != 0) continue;
		hV0PtFitCorr[Sp][0][iMu][iSph]->Write();
	}	}

	printf("Output file %s created. \n", fileName.Data()); 
	outFile->Close();
	mDirFile->cd();

}