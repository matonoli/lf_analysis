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

#include "MyAnalysisV0correct.h"
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

	hEventType = (TH1D*)mHandler->analysis(0)->dirFile()->Get("hEventType");

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0PtFit[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		
	} } } }

}

Bool_t MyAnalysisV0correct::CreateHistograms() {
	
}


Bool_t MyAnalysisV0correct::CloneHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0PtFitCorr[iSp][iType][iMu][iSph]	= (TH1D*)hV0PtFit[iSp][iType][iMu][iSph]->Clone(
			Form("hV0PtFitCorr_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]) );

		
	} } } }

}

Int_t MyAnalysisV0correct::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	CloneHistograms();
	NormaliseSpectra();
	LoadEfficiency();
	CorrectSpectra();

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

void MyAnalysisV0correct::NormaliseSpectra() {

	/*for (int iBin = 0; iBin < NEVENTTYPES+2; ++iBin)
	{
		printf("Event type %i containing %f events\n", iBin, hEventType->GetBinContent(iBin));
	}*/

	Double_t NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
	printf("Normalising all histograms by dEta %f \n", NormEta);

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		Double_t NormEv = hEventType->GetBinContent(2);	// MB
		if (iMu == 1)	{	NormEv = hEventType->GetBinContent(3);		// FHM
			if (iSph == 1)	NormEv = hEventType->GetBinContent(8);		// FHM JET
			if (iSph == 2)	NormEv = hEventType->GetBinContent(7);	}	// FHM ISO
		if (iMu == 2)	{	NormEv = hEventType->GetBinContent(4);		// MHM
			if (iSph == 1)	NormEv = hEventType->GetBinContent(10);		// MHM JET
			if (iSph == 2)	NormEv = hEventType->GetBinContent(9);	}	// MHM ISO
		if (iMu == 3)	{	NormEv = hEventType->GetBinContent(11);		// RT
			if (iSph == 3)	NormEv = hEventType->GetBinContent(12);		// RT 0-1
			if (iSph == 4)	NormEv = hEventType->GetBinContent(13);		// RT 1-2
			if (iSph == 5)	NormEv = hEventType->GetBinContent(14);		// RT 2-3
			if (iSph == 6)	NormEv = hEventType->GetBinContent(15);		// RT 3-4
			if (iSph == 7)	NormEv = hEventType->GetBinContent(16);	}	// RT 4-5
		if (iType == 1) {
			if (iMu == 1) {
				if (iSph == 1 )	NormEv = hEventType->GetBinContent(18);		// FHM JET MC
				if (iSph == 2 )	NormEv = hEventType->GetBinContent(17);	}	// FHM ISO MC
			if (iMu == 2) {
				if (iSph == 1 )	NormEv = hEventType->GetBinContent(20);		// MHM JET MC
				if (iSph == 2 )	NormEv = hEventType->GetBinContent(19);	}	// MHM ISO MC
		}
		

		if (NormEv == 0) NormEv = 1;

		printf("Normalising histogram %s by event count %f \n", hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName(), NormEv);
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEv);
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEta);
		
	} } } }
	
}

void MyAnalysisV0correct::LoadEfficiency() {

	if (!mFileMC) {
		printf("No MC file loaded! Efficiency correction not performed.\n");
		return;
	}

	TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	
		hV0Efficiency[iSp] = (TH1D*)dirFile1->Get(Form("hV0Efficiency_%s",SPECIES[iSp]));
		hV0Efficiency[iSp]->Write();
	}

}

void MyAnalysisV0correct::CorrectSpectra() {
	
	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	//for (int iMu = 0; iMu < 3; ++iMu)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{

		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		//hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1,"width");
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Divide(hV0Efficiency[iSp]);

	}	}	}	}

}