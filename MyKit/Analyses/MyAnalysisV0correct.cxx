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
#include <TCutG.h>

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
		Int_t iType = 2; Int_t iMu = 0; Int_t iSph = 0;
		hV0Pt[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
			
	}


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

	
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
		
		hV0FeeddownMotherPt[iSp]	= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0FeeddownMotherPt_%s",SPECIES[iSp]) );
		
	}		// should be loaded from external files instead



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


	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < 2; ++iMu)		{
		hV0PtFeeddown[iSp][iMu]	= (TH1D*)hV0PtFitCorr[iSp][0][iMu][0]->Clone(
			Form("hV0PtFeeddown_%s_%s",SPECIES[iSp],MULTI[iMu]));
	}	}

}

Int_t MyAnalysisV0correct::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	CloneHistograms();
	
	NormaliseSpectra();
	CorrectForFeeddown();

	LoadEfficiency();
	//DoEfficiencyFromFile();
	CorrectSpectra();

	//if (!mHandler->GetFlagMC()) StudyCuts();

	DoXCheckV0M();
	if (mHandler->GetFlagMC()) DoClosureTest(0);


	printf("mb k0s spectrum final \n");
	cout << hV0PtFitCorr[1][0][0][0]->GetBinContent(30) << " at " << hV0PtFitCorr[1][0][0][0]->GetBinLowEdge(30) << endl;

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

void MyAnalysisV0correct::CorrectForFeeddown() {

	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
		hXiPt[iSp][iMu] = 0x0;
	}	}

	TCutG* cutg;
	TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
		hV0FeeddownMatrix[iSp]	= (TH2D*)dirFile1->Get(Form("hV0FeeddownMatrix_%s",SPECIES[iSp]) );
		if (!hV0FeeddownMatrix[iSp]) {
			printf("Feed-down matrix not found, no correction performed.\n");
			return;
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
	}

	TH2D* hFMatrix = (TH2D*)hV0FeeddownMatrix[2]->Clone("hFMatrix");
	for (int iMotherBin = 1; iMotherBin < hV0FeeddownMatrix[2]->GetNbinsX()+1; ++iMotherBin)	{
		for (int iPtBin = 1; iPtBin < hV0FeeddownMatrix[2]->GetNbinsY()+1; ++iPtBin)	{
				
				Double_t cX = hV0FeeddownMatrix[2]->GetXaxis()->GetBinCenter(iMotherBin);
				Double_t cY = hV0FeeddownMatrix[2]->GetYaxis()->GetBinCenter(iPtBin);
				hFMatrix->SetBinContent(iMotherBin,iPtBin,0);
				if (!cutg->IsInside(cX,cY)) continue;
				hFMatrix->SetBinContent(iMotherBin,iPtBin,hV0FeeddownMatrix[2]->GetBinContent(iMotherBin,iPtBin));
	}	}

	mFileXi = new TFile("../official/xi_HM_spectra_sep_9_2019.root","READ");
	if (!mFileXi) {
		printf("No Xi file loaded in, using MC Xi spectra instead. \n");
		hXiPt[2][0] = hV0FeeddownMotherPt[2];
		hXiPt[3][0] = hV0FeeddownMotherPt[3];
	} else {

		hXiPt[2][1] = (TH1D*)mFileXi->Get("hHMSpectrum_HM");
		hXiPt[3][1] = (TH1D*)mFileXi->Get("hHMSpectrum_HM");

		TFile* fileXIMB = new TFile("../official/xi_MB_spectra_sep_9_2019_sum_rebin.root","READ");
		hXiPt[2][0] = (TH1D*)fileXIMB->Get("hXiSumSpectrum_MB");
		hXiPt[3][0] = (TH1D*)fileXIMB->Get("hXiSumSpectrum_MB");

		if (0) {/////// change this back
		hXiPt[2][0] = (TH1D*)hV0FeeddownMotherPt[2]->Clone("hXiPt_L_MB");
		hXiPt[3][0] = (TH1D*)hV0FeeddownMotherPt[3]->Clone("hXiPt_Lbar_MB");
		// these need to be normalised
		NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
		Double_t NormEv = hEventType->GetBinContent(24);	// MB
		if (NormEv>0) {
			NormEv += NormEv * hEventType->GetBinContent(22) * 1./(hEventType->GetBinContent(23) + hEventType->GetBinContent(24));
		}
		if (NormEv == 0) NormEv = 1;

		printf("Normalising histogram %s by event count %f \n", hXiPt[2][0]->GetName(), NormEv);
		hXiPt[2][0]->Scale(2.,"width");	//in file xi is xi+xibar
		hXiPt[2][0]->Scale(1./NormEv);
		hXiPt[2][0]->Scale(1./NormEta);
		printf("Normalising histogram %s by event count %f \n", hXiPt[3][0]->GetName(), NormEv);
		hXiPt[3][0]->Scale(2.,"width");
		hXiPt[3][0]->Scale(1./NormEv);
		hXiPt[3][0]->Scale(1./NormEta);
		}//////

		mDirFile->cd();
	}

	// INTERPOLATE XI SPECTRA
	TF1* funcLT = LevyTsallis("LT",XIMASS);
	if (!mFileXi) funcLT->SetParameter(3,100000);
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	
		if (!hXiPt[iSp][iMu]) continue;
		cout << iSp << " " << iMu << " " << hXiPt[iSp][iMu] << endl;
		//if (mFileXi) hXiPt[iSp][iMu]->Scale(0.5);
		hXiPt[iSp][iMu]->Fit(funcLT,"I");

		hXiPt[iSp][iMu]->Write();

		if (hV0PtFeeddown[iSp][iMu]->GetNbinsX() != hV0FeeddownMatrix[iSp]->GetNbinsY()) {
			printf("FD matrix has bad dimensions!\n");
			return;
		}

		for (int iBin = 1; iBin < hV0PtFeeddown[iSp][iMu]->GetNbinsX()+1; ++iBin)	{
			Double_t sum = 0;
			Double_t sumErr = 0;
			for (int iMotherBin = 1; iMotherBin < hV0FeeddownMatrix[iSp]->GetNbinsX()+1; ++iMotherBin)	{
				
				Double_t cX = hV0FeeddownMatrix[iSp]->GetXaxis()->GetBinCenter(iMotherBin);
				Double_t cY = hV0FeeddownMatrix[iSp]->GetYaxis()->GetBinCenter(iBin);
				if (!cutg->IsInside(cX,cY)) continue;

				Double_t left 	= hV0FeeddownMatrix[iSp]->GetXaxis()->GetBinLowEdge(iMotherBin);
				Double_t right 	= hV0FeeddownMatrix[iSp]->GetXaxis()->GetBinLowEdge(iMotherBin+1);
				if (iBin>7 && iBin <10) printf("bin %i range %f - %f adding %f times %f -> sum before %f \n", iBin, left, right, 
					funcLT->Integral(left,right), hV0FeeddownMatrix[iSp]->GetBinContent(iMotherBin,iBin), sum);
				sum +=
					hV0FeeddownMatrix[iSp]->GetBinContent(iMotherBin,iBin) *
					funcLT->Integral(left,right) / (2.*(right-left)); //(2.*(right-left)); //2 because xi,xibar

				sumErr +=
					(hV0FeeddownMatrix[iSp]->GetBinError(iMotherBin,iBin) *
					hV0FeeddownMatrix[iSp]->GetBinError(iMotherBin,iBin)) +
					(funcLT->IntegralError(left,right) *
					funcLT->IntegralError(left,right));	
			}
			hV0PtFeeddown[iSp][iMu]->SetBinContent(iBin, 2*sum);
			hV0PtFeeddown[iSp][iMu]->SetBinError(iBin, 0);//TMath::Sqrt(sumErr));
		}
		hV0PtFeeddown[iSp][iMu]->Scale(1,"width");

		TH1D* htmp = (TH1D*)hV0PtFitCorr[iSp][0][iMu][0]->Clone("htmp");
		hV0PtFitCorr[iSp][0][iMu][0]->Add(hV0PtFeeddown[iSp][iMu],-1.);
		hV0PtFeeddown[iSp][iMu]->Divide(htmp);
		delete htmp;

	}	}


	funcLT->Write();

	cout << "bc " << hV0PtFeeddown[2][1]->GetBinContent(6) << " / " << hV0PtFitCorr[2][0][1][0]->GetBinContent(6) << endl;


}

void MyAnalysisV0correct::NormaliseSpectra() {

	/*for (int iBin = 0; iBin < NEVENTTYPES+2; ++iBin)
	{
		printf("Event type %i containing %f events\n", iBin, hEventType->GetBinContent(iBin));
	}*/

	NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
	//Double_t NormEta = (cuts::V0_Y[1] - cuts::V0_Y[0]);
	printf("Normalising all histograms by dY %f \n", NormEta);

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		Double_t NormEv = hEventType->GetBinContent(24);	// MB
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
		
		if (NormEv>0) {
			NormEv += NormEv * hEventType->GetBinContent(22) * 1./(hEventType->GetBinContent(23) + hEventType->GetBinContent(24));
		}

		if (NormEv == 0) NormEv = 1;

		printf("Normalising histogram %s by event count %f \n", hV0PtFitCorr[iSp][iType][iMu][iSph]->GetName(), NormEv);
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEv);
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1./NormEta);
		
	} } } }

	printf("mb k0s spectrum from fit \n");
	cout << hV0PtFit[1][0][0][0]->GetBinContent(30) << endl;
	printf("mb k0s spectrum from norm \n");
	cout << hV0PtFitCorr[1][0][0][0]->GetBinContent(30) << endl;


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
		//cout << "bin a " << binA << " w center " << hRt2Rebin->GetBinLowEdge(iB) << endl;
		//cout << "bin b " << binA << " w center " << hRt2Rebin->GetBinLowEdge(iB+1) << endl;
		Double_t integr = 0;
		for (int i = binA; i < binB; ++i)
		{
			//cout << "in bin " << i << " integrating " << integr << " + " << hNchTrans->GetBinContent(i) << endl;
			integr += hNchTrans->GetBinContent(i);
		}
		Double_t error;
		Double_t content 	= hNchTrans->IntegralAndError(binA,binB-1,error,"");
		//cout << "integral is " << integr << " vs " << content << endl;
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
		else { 
			//NormEv += NormEv * hEventType->GetBinContent(22) * 1./(hEventType->GetBinContent(23) + hEventType->GetBinContent(24));
		}

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

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		
		// minimum bias histos
		hV0Efficiency[iSp] = new TH1D(Form("hV0Efficiency_%s",SPECIES[iSp]),"; V0 pT (GeV/#it{c}); Efficiency",NPTBINS,XBINS); //NPTBINS2, XBINS2);
		TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

		tV0massRCMB[iSp]->Draw(Form("lPt>>hV0Efficiency_%s",SPECIES[iSp]),"","goff");
		tV0PtMCMB[iSp]->Draw("lPt>>hDen","","goff");

		hV0Efficiency[iSp]->GetYaxis()->SetRangeUser(0.,0.65);
		hV0Efficiency[iSp]->GetXaxis()->SetRangeUser(0.,15.0);

		hV0Efficiency[iSp]->Divide(hDen);
		delete hDen;



		// rt histos
		for (int iReg = 0; iReg < NREGIONS; ++iReg)		{		
			hV0EfficiencyRt[iSp][iReg] = new TH1D(Form("hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]),"; V0 pT (GeV/#it{c}); Efficiency",NPTBINS,XBINS);
			TH1D* hDen = (TH1D*)hV0Efficiency[iSp]->Clone("hDen"); // denominator with same binning

			tV0massRt[iSp][1][iReg]->Draw(Form("lPt>>hV0EfficiencyRt_%s_%s",SPECIES[iSp],REGIONS[iReg]),"","goff");
			tV0PtMCRt[iSp][iReg]->Draw("lPt>>hDen","","goff");

			//hV0EfficiencyRt[iSp][iReg]->GetYaxis()->SetRangeUser(0.,0.65);
			hV0EfficiencyRt[iSp][iReg]->GetXaxis()->SetRangeUser(0.,15.0);

			hV0EfficiencyRt[iSp][iReg]->Divide(hDen);
			delete hDen;
		}
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
		leg1->AddEntry(hV0Efficiency[1],Form("%s",SPECNAMES[1]),"pl");
		leg1->AddEntry(hV0Efficiency[2],Form("%s",SPECNAMES[2]),"pl");
		leg1->AddEntry(hV0Efficiency[3],Form("%s",SPECNAMES[3]),"pl");
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
	cout << hV0PtFitCorr[2][0][1][0]->GetBinContent(30) << endl;

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	TF1* funcRapCorrection = new TF1("funcRapCorrection",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	//for (int iMu = 0; iMu < 3; ++iMu)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{

		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		//hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(1,"width");
		hV0PtFitCorr[iSp][iType][iMu][iSph]->Divide(hV0Efficiency[iSp]);

		funcRapCorrection->SetParameters(NormEta,MASSES[iSp]);
		
		//cout << "at 1.0 correcting by " << funcRapCorrection->Eval(1.0) << endl;

		hV0PtFitCorr[iSp][iType][iMu][iSph]->Divide(funcRapCorrection,1.);

		hV0PtFitCorr[iSp][iType][iMu][iSph]->Scale(MBtrigEff);



	}	}	}	}
	funcRapCorrection->Write();


	printf("mb k0s eff \n");
	cout << hV0Efficiency[2]->GetBinContent(30) << endl;

	printf("mb k0s spectrum after corr \n");
	cout << hV0PtFitCorr[2][0][1][0]->GetBinContent(30) << endl;

	TF1* funcRapCorrection2 = new TF1("funcRapCorrection2",rap_correction,XBINS2[0],XBINS2[NPTBINS2],2);
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Divide(hV0Efficiency[iSp]);

		funcRapCorrection2->SetParameters(NormEta,MASSES[iSp]);
		
		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Divide(funcRapCorrection2,1.);

		hV0PtRtFitCorr[iSp][iType][iReg][iRtBin]->Scale(MBtrigEff);
	}	}	}	}


}

void MyAnalysisV0correct::StudyCuts() {

	for (int iCut = 0; iCut < 25; ++iCut)	{
		hV0PtCut[iCut]
			= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0PtCut_%i",iCut) );
		
	}

	if (!mFileMC) {
		printf("No MC file loaded!.\n");
		return;
		} else {
		TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");
		for (int iCut = 0; iCut < 25; ++iCut)	{
			TH1D* htmp
				= (TH1D*)dirFile1->Get(Form("hV0PtCut_%i",iCut) );
			hV0PtCutMC[iCut] = (TH1D*)htmp->Clone(Form("hV0PtCutMC_%i",iCut));
			
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
	Double_t NormEv = hEventType->GetBinContent(24);	// MB
	if (NormEv>0) {
			NormEv += NormEv * hEventType->GetBinContent(22) * 1./(hEventType->GetBinContent(23) + hEventType->GetBinContent(24));
		}
	Double_t MBtrigEff = 0.7448;
	TF1* funcRapCorrection = new TF1("funcRapCorrection",rap_correction,XBINS[0],XBINS[NPTBINS],2);


	for (Int_t iSp = 1; iSp < NSPECIES; iSp++)		{
		Int_t iMu = 0; Int_t iSph = 0;	
		hClosureTestCorr[iSp]	= (TH1D*)hV0PtFitCorr[iSp][0][iMu][iSph]->Clone(Form("hClosureTestCorr_%s",SPECIES[iSp]));
		cout << " " << hV0Pt[iSp][2][iMu][iSph] << endl;
		TH1D* hDen = (TH1D*)hV0Pt[iSp][2][iMu][iSph]->Clone(Form("hDen"));
		hDen->Scale(1./NormEv);
		hDen->Scale(1./NormEta);
		hDen->Scale(MBtrigEff);
		hDen->Scale(1.,"width");
		funcRapCorrection->SetParameters(NormEta,MASSES[iSp]);
		hDen->Divide(funcRapCorrection,1.);

		mHandler->MakeNiceHistogram(hClosureTestCorr[iSp],kBlack);
		hClosureTestCorr[iSp]->GetYaxis()->SetTitle("blind rec. / MC generated");
		hClosureTestCorr[iSp]->GetYaxis()->SetRangeUser(0.7,1.3);
		hClosureTestCorr[iSp]->Divide(hDen);
		TCanvas* cClosureCorr = new TCanvas("cClosureCorr","",900,900);
		hClosureTestCorr[iSp]->Draw();
		TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
		mHandler->MakeNiceLegend(leg1, 0.05, 1.);
		leg1->AddEntry((TObject*)0,Form("%s",SPECNAMES[iSp])," ");
		leg1->Draw();
		cClosureCorr->SaveAs(Form("tmp/closureCorr_%s.png",SPECIES[iSp]));
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

	TH1D* hLLbarMB = (TH1D*)hV0PtFitCorr[2][0][0][0]->Clone("hLLbarMB");
	hLLbarMB->Add(hV0PtFitCorr[3][0][0][0]);
	hLLbarMB->Scale(1./0.7448);
	TH1D* hLLbarV0M = (TH1D*)hV0PtFitCorr[2][0][1][0]->Clone("hLLbarV0M"); 
	hLLbarV0M->Add(hV0PtFitCorr[3][0][1][0]);
	hLLbarV0M->Scale(1./0.7448);

	TF1* funcRapCorrectionL = new TF1("funcRapCorrectionL",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	funcRapCorrectionL->SetParameters(1.6,MASSES[2]);
	//hLLbarMB->Divide(funcRapCorrectionL,1);
	//hLLbarV0M->Divide(funcRapCorrectionL,1);


	//hLLbarV0M->Divide(hOffiL);
	//hLLbarMB->Divide(hOffiLMB);

	
	if (1) {
		TCanvas* cXcheck = new TCanvas("cXcheck","",900,900);
		cXcheck->SetLogy();
		mHandler->MakeNiceHistogram((TH1D*)hOffiL,kGreen+2);
		mHandler->MakeNiceHistogram((TH1D*)hOffiLMB,kGreen+2);
		hOffiLMB->SetMarkerStyle(21);
		mHandler->MakeNiceHistogram(hLLbarV0M,kRed);
		mHandler->MakeNiceHistogram(hLLbarMB,kBlue);
		hLLbarMB->SetMarkerStyle(21);
		hOffiL->Draw();
		hLLbarV0M->Draw("same");
		hLLbarMB->Draw("same");
		hOffiLMB->Draw("same");
	}

	TH1D* hK0sMB = (TH1D*)hV0PtFitCorr[1][0][0][0]->Clone("hK0sMB");
	TH1D* hK0sMBeff = (TH1D*)hV0Efficiency[1]->Clone("hK0sMBeff");
	TH1D* hK0sMBraw = (TH1D*)hV0PtFit[1][0][0][0]->Clone("hK0sMBraw");
	hK0sMBraw->Scale(1./1.6);

	Double_t NormEv = hEventType->GetBinContent(24);
	NormEv += NormEv * hEventType->GetBinContent(22) * 1./(hEventType->GetBinContent(23) + hEventType->GetBinContent(24));
	hK0sMBraw->Scale(1./NormEv);

	cout << "badawda " << hK0sMB << endl;
	cout << "badawda222 " << hK0sMBeff << endl;


	printf("atghh mb k0s eff \n");
	cout << hV0Efficiency[1]->GetBinContent(30) << endl;
	cout << hK0sMBeff->GetBinContent(30) << endl;

	TF1* funcRapCorrection3 = new TF1("funcRapCorrection3",rap_correction,XBINS[0],XBINS[NPTBINS],2);
	funcRapCorrection3->SetParameters(1.6,MASSES[1]);
	hOffiKMBeff->Divide(funcRapCorrection3,1);

	hK0sMB->Divide(hOffiKMB);
	hOffiKMBeff->Divide(hK0sMBeff);
	hK0sMBraw->Divide(hOffiKMBraw);

	TCanvas* cXcheck2 = new TCanvas("cXcheck2","",900,900);
	//cXcheck2->SetLogy();
	//hOffiLMB->SetMarkerStyle(21);
	mHandler->MakeNiceHistogram(hK0sMB,kRed);
	mHandler->MakeNiceHistogram(hK0sMBeff,kBlue);
	mHandler->MakeNiceHistogram(hK0sMBraw,kBlue);
	//hLLbarMB->SetMarkerStyle(21);
	hK0sMB->Draw();
	//hOffiKMB->Draw("same");
	
	//hK0sMBeff->Draw("same");
	hOffiKMBeff->Draw("same");
	hK0sMBraw->Draw("same");
		
	//mHandler->root()->SetBatch(kFALSE);
}