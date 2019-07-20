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
#include <TNtuple.h>

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

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	hEventType	= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hEventType");
	hNchTrans	= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hNchTrans");
	hRt			= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hRt");
	hRt2		= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hRt2");
	
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{
		hRtV0Yields[iType][iReg][iPtBin]	
			= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hRtV0Yields_%s_%s_%i",TYPE[iType],REGIONS[iReg],iPtBin));
	} } }

	
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0PtFit[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		
	} } } }


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{
		hV0RtFit[iSp][iType][iReg][iPtBin] 
			= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hV0RtFit_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iPtBin));
	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
			
		hV0PtRtFit[iSp][iType][iReg][iRtBin] 
			= (TH1D*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtRtFit_%s_%s_%s_%1.1f-%1.1f",SPECIES[iSp],TYPE[iType],REGIONS[iReg],RTBINS0[iRtBin],RTBINS0[iRtBin+1]) );


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

		//if (hV0PtFitCorr[iSp][iType][iMu][iSph] && mHandler->IsRebinPt()) {
		//hV0PtFitCorr[iSp][iType][iMu][iSph] = (TH1D*)hV0PtFitCorr[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName(),XBINS2); }
		
	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{

		hV0RtFitCorr[iSp][iType][iReg][iPtBin] = (TH1D*)hV0RtFit[iSp][iType][iReg][iPtBin]->Clone(
			Form("hV0RtFitCorr_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iPtBin) );
	}	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
			
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin] = (TH1D*)hV0PtRtFit[iSp][iType][iReg][iRtBin]->Clone(
			Form("hV0PtRtFitCorr_%s_%s_%s_%1.1f-%1.1f",SPECIES[iSp],TYPE[iType],REGIONS[iReg],RTBINS0[iRtBin],RTBINS0[iRtBin+1]) );


	} } } }

}

Int_t MyAnalysisV0correct::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	CloneHistograms();
	NormaliseSpectra();
	//LoadEfficiency();
	DoEfficiencyFromFile();
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
		if (iMu == 3 || iMu == 4 || iMu == 5)	{	NormEv = hEventType->GetBinContent(11);		// RT
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
		//hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEv);
		//hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEta);
		
	} } } }


	// normalise rt spectra
	Int_t increm = 2;
	Double_t rt_den = hNchTrans->GetMean();
	hRtRebin		= (TH1D*)hRt->Rebin(NRTBINS,"hRtRebin",RTBINS);
	hRt2Rebin		= (TH1D*)hV0RtFitCorr[1][0][0][0]->Clone("hRt2Rebin");
	Int_t nbins = hRt2Rebin->GetNbinsX();
	Double_t rtbins[nbins+1];

	cout << "nbins is " << nbins << endl;

	for (int iB = 0; iB < nbins; ++iB)	{
		Int_t binA = hNchTrans->FindBin(hRt2Rebin->GetBinLowEdge(iB));
		Int_t binB = hNchTrans->FindBin(hRt2Rebin->GetBinLowEdge(iB+1));
		cout << "bin a " << binA << " w center " << hRt2Rebin->GetBinLowEdge(iB) << endl;
		cout << "bin b " << binA << " w center " << hRt2Rebin->GetBinLowEdge(iB+1) << endl;
		Double_t integr = 0;
		for (int i = binA; i < binB; ++i)
		{
			cout << "in bin " << i << " integrating " << integr << " + " << hNchTrans->GetBinContent(i) << endl;
			integr += hNchTrans->GetBinContent(i);
		}
		Double_t error;
		Double_t content 	= hNchTrans->IntegralAndError(binA,binB-1,error,"");
		cout << "integral is " << integr << " vs " << content << endl;
		hRt2Rebin->SetBinContent(iB,content);
		hRt2Rebin->SetBinError(iB,error);
	}
	//hNchTransRebin	= (TH1D*)hNchTrans->Rebin(increm);

	//im dumb
	cout << "a " << hV0RtFitCorr[1][0][0][0]->Integral(hV0RtFitCorr[1][0][0][0]->FindBin(30.),hV0RtFitCorr[1][0][0][0]->FindBin(30.)) << endl;
	cout << "b " << hRt2Rebin->Integral(hRt2Rebin->FindBin(30.),hRt2Rebin->FindBin(30.)) << endl;
	cout << "c " << hRtV0Yields[0][0][0]->GetBinContent(1+1) << endl;
	cout << "d " << hRt2Rebin->Integral(1,-1) << endl;

	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
		
		hV0RtFitCorr[iSp][iType][iReg][iPtBin]->Divide(hRt2Rebin);
		// now we have #v0s/per event(rtbin)

		//cout << "rt " << hRtV0Yields << endl;
		Double_t meanV0 = hRtV0Yields[iType][iReg][iPtBin]->GetBinContent(1+iSp);
		meanV0 = meanV0/hRt2Rebin->Integral(1,-1);
		//cout << "2 meanV0 " << meanV0 << endl;
		if (meanV0>0) hV0RtFitCorr[iSp][iType][iReg][iPtBin]->Scale(1./meanV0);

		for (int iBin = 0; iBin < nbins+1; ++iBin)	{
		rtbins[iBin] = (double)hV0RtFitCorr[iSp][iType][iReg][iPtBin]->GetBinLowEdge(iBin+1)/rt_den;	}

		hV0RtFitCorr[iSp][iType][iReg][iPtBin]->SetBins(nbins,rtbins);
		hV0RtFitCorr[iSp][iType][iReg][iPtBin]->GetXaxis()->SetRangeUser(rtbins[0],5.1);
		hV0RtFitCorr[iSp][iType][iReg][iPtBin]->GetYaxis()->SetRangeUser(0.,10.1);

	}	}	}	}





	// normalise pt spectra (rt) created from trees
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{

		Double_t NormEv = 0;
		//cout << "bla " << RTBINS0[iRtBin] << " " << RTBINS0[iRtBin]*rt_den <<
		//" bla " << RTBINS0[iRtBin+1] << " " << RTBINS0[iRtBin+1]*rt_den << endl;
		for (int ib = 1; ib < hNchTrans->GetNbinsX(); ++ib)
		{
			//cout << "aaaa " << 
			if (hNchTrans->GetBinCenter(ib) < RTBINS0[iRtBin]*rt_den-1E-05) continue;
			if (hNchTrans->GetBinCenter(ib) > RTBINS0[iRtBin+1]*rt_den) continue;
			NormEv += hNchTrans->GetBinContent(ib);
			//cout << "ib " << ib << " bins " << hNchTrans->GetBinCenter(ib) << " sum " << sum << endl;
		}

		if (NormEv == 0) NormEv = 1;
		printf("Normalising histogram %s by event count %f \n", hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->GetName(), NormEv);
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Scale(1./NormEv);
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Scale(1./NormEta);


		/*Int_t leftbin = hNchTrans->FindBin(RTBINS0[iRtBin]*rt_den);
		Int_t rightbin = hNchTrans->FindBin(RTBINS0[iRtBin+1]*rt_den);
		Double_t NormEv	= hNchTrans->Integral(leftbin,rightbin-1);

		cout << iSp << " " << iType << " " << iReg << " " << iRtBin << endl;
		cout << "bla " << RTBINS0[iRtBin] << " " << RTBINS0[iRtBin]*rt_den <<
		" bla " << RTBINS0[iRtBin+1] << " " << RTBINS0[iRtBin+1]*rt_den << " " << leftbin << " " << rightbin << endl;
		for (int ib = 0; ib < rightbin-leftbin+1; ++ib)
		{
			cout << "irtbin " << iRtBin << " lb " << leftbin << " rb " << rightbin << " bins " << hNchTrans->GetXaxis()->GetBinCenter(leftbin+ib) << endl;
		}*/
		/*Int_t leftbin = hRt2->FindBin(RTBINS0[iRtBin]);
		Int_t rightbin = hRt2->FindBin(RTBINS0[iRtBin+1]);
		Double_t NormEv	= hRt2->Integral(leftbin,rightbin-1);*/

	}	}	}	}

	
}

void MyAnalysisV0correct::LoadEfficiency() {

	if (!mFileMC) {
		printf("No MC file loaded! Efficiency correction not performed.\n");
		return;
	}

	TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	
		hV0Efficiency[iSp] = (TH1D*)dirFile1->Get(Form("hV0Efficiency_%s",SPECIES[iSp]));

		if (hV0Efficiency[iSp] && mHandler->IsRebinPt()) {

			Int_t binSize = TMath::Nint((Double_t)NPTBINS/NPTBINS2);
			hV0Efficiency[iSp] = (TH1D*)hV0Efficiency[iSp]->Rebin(NPTBINS2,hV0Efficiency[iSp]->GetName(),XBINS2);
			hV0Efficiency[iSp]->Scale(1./binSize);		 }

		hV0Efficiency[iSp]->Write();
	}

}

void MyAnalysisV0correct::DoEfficiencyFromFile() {

	if (!mFileMC) {
		printf("No MC file loaded! Efficiency correction not performed.\n");
		return;
	}

	TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");

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

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		
		// minimum bias histos
		hV0Efficiency[iSp] = new TH1D(Form("hV0Efficiency_%s",SPECIES[iSp]),"; V0 pT (GeV/#it{c}); Efficiency",NPTBINS2,XBINS2);
		TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

		tV0massRCMB[iSp]->Draw(Form("lPt>>hV0Efficiency_%s",SPECIES[iSp]),"","goff");
		tV0PtMCMB[iSp]->Draw("lPt>>hDen","","goff");

		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,14.0);

		hV0Efficiency[iSp]->Divide(hDen);
		delete hDen;

		// rt histos
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{		
			hV0EfficiencyRt[iSp][iReg] = new TH1D(Form("hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]),"; V0 pT (GeV/#it{c}); Efficiency",NPTBINS2,XBINS2);
			TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

			tV0massRt[iSp][1][iReg]->Draw(Form("lPt>>hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]),"","goff");
			tV0PtMCRt[iSp][iReg]->Draw("lPt>>hDen","","goff");

			hV0EfficiencyRt[iSp][iReg]->GetYaxis()->SetRangeUser(0.,0.65);
			hV0EfficiencyRt[iSp][iReg]->GetXaxis()->SetRangeUser(0.,14.0);

			hV0EfficiencyRt[iSp][iReg]->Divide(hDen);
			delete hDen;
		}
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


	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Divide(hV0Efficiency[iSp]);
	}	}	}	}


}