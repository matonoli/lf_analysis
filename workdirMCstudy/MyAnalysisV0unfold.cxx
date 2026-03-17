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
#include <TCut.h>
#include <TGaxis.h>

#include "MyAnalysisV0unfold.h"
#include "MyEvent.h"
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

#include "UnfoldNTclass.h"

//#include <AliAnalysisPIDV0.h>
using namespace V0consts;
//using namespace RooFit;
using namespace std;

ClassImp(MyAnalysisV0unfold)



MyAnalysisV0unfold::MyAnalysisV0unfold() {

}

Int_t MyAnalysisV0unfold::Init() {

	mDirFile->cd();

	printf("Initialising analysis %s  \n", 
		this->GetName());

	TH1::SetDefaultSumw2(1);
	BorrowHistograms();
	NormaliseMC();
	CreateHistograms();
	


	mList = (TList*)mHandler->directory()->GetList();

	return 0;
}

Int_t MyAnalysisV0unfold::Make(Int_t iEv) {
	return 0;
}


Bool_t MyAnalysisV0unfold::BorrowHistograms() {

	printf("Borrowing histograms for analysis %s  \n", 
		this->GetName());
	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNtFitCorr[iSp][iType][iReg]		
			= (TH2F*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtNtFitCorr_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hV0PtNtMinFitCorr[iSp][iType][iReg]		
			= (TH2F*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtNtMinFitCorr_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hV0PtNtMaxFitCorr[iSp][iType][iReg]		
			= (TH2F*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtNtMaxFitCorr_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		
	}	}	}

	hNchTransMCTrigMC		= (TH1F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMCTrigMC");

	for (int iSp = 0; iSp < NSPECIES; ++iSp)		{
	for (int iType = 2; iType < NTYPE; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
				
		hV0PtNt[iSp][iType][iReg]
			= (TH2F*)mFileMC->Get(Form("MyAnalysisV0_0/hV0PtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hV0PtNt[iSp][iType][iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0PtNt[iSp][iType][iReg]);
		
		hV0PtNtMin[iSp][iType][iReg]
			= (TH2F*)mFileMC->Get(Form("MyAnalysisV0_0/hV0PtNtMin_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hV0PtNtMin[iSp][iType][iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0PtNtMin[iSp][iType][iReg]);

		hV0PtNtMax[iSp][iType][iReg]
			= (TH2F*)mFileMC->Get(Form("MyAnalysisV0_0/hV0PtNtMax_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hV0PtNtMax[iSp][iType][iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0PtNtMax[iSp][iType][iReg]);

	}	}	}
	

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

		hV0PtRtSysSum[iSp][iReg][iRtBin]		
			= (TH1F*)mHandler->analysis(3)->dirFile()->Get(Form("hV0PtRtSysSum_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSumUnc[iSp][iReg][iRtBin]		
			= (TH1F*)mHandler->analysis(3)->dirFile()->Get(Form("hV0PtRtSysSumUnc_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));

		cout << hV0PtRtSysSum[iSp][iReg][iRtBin] << " " << hV0PtRtSysSumUnc[iSp][iReg][iRtBin]	<< endl;

	}	}	}

	return true;
}

Bool_t MyAnalysisV0unfold::CreateHistograms() {

	printf("Creating histograms for analysis %s  \n", 
		this->GetName());

	
	return true;
}

Bool_t MyAnalysisV0unfold::NormaliseMC() {

	Double_t NormEta = (cuts::V0_ETA[1] - cuts::V0_ETA[0]);
	printf("Normalising all histograms by dEta %f \n", NormEta);

	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 2; iType < NTYPE; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{

		printf("Normalising histogram %s by N_T distribution \n", hV0PtNt[iSp][iType][iReg]->GetName());

		for (int iNt = 1; iNt < hV0PtNt[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {

			Double_t NormEv = 0;
			
			NormEv = hNchTransMCTrigMC->Integral(1,50);
			TH2F* htmp = hV0PtNt[iSp][iType][iReg];

			TH1D* hpt = (TH1D*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

			htmp = hV0PtNtMin[iSp][iType][iReg];

			hpt = (TH1D*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

			htmp = hV0PtNtMax[iSp][iType][iReg];

			hpt = (TH1D*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEta>0) hpt->Scale(1./NormEta);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );

			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}

		}

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

		for (int iX = 1; iX < hV0PtNtMin[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNtMin[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float binwidth = hV0PtNtMin[iSp][iType][iReg]->GetXaxis()->GetBinWidth(iX);
				float binc = hV0PtNtMin[iSp][iType][iReg]->GetBinContent(iX,iY);
				float bine = hV0PtNtMin[iSp][iType][iReg]->GetBinError(iX,iY);
				if (binwidth>0) hV0PtNtMin[iSp][iType][iReg]->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) hV0PtNtMin[iSp][iType][iReg]->SetBinError(iX,iY,bine/binwidth);
			}
		}

		for (int iX = 1; iX < hV0PtNtMax[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNtMax[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float binwidth = hV0PtNtMax[iSp][iType][iReg]->GetXaxis()->GetBinWidth(iX);
				float binc = hV0PtNtMax[iSp][iType][iReg]->GetBinContent(iX,iY);
				float bine = hV0PtNtMax[iSp][iType][iReg]->GetBinError(iX,iY);
				if (binwidth>0) hV0PtNtMax[iSp][iType][iReg]->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) hV0PtNtMax[iSp][iType][iReg]->SetBinError(iX,iY,bine/binwidth);
			}
		}

		hV0PtNt[iSp][iType][iReg]->Write();
		hV0PtNtMin[iSp][iType][iReg]->Write();
		hV0PtNtMax[iSp][iType][iReg]->Write();

	}	}	}
	
	return true;
}


Bool_t MyAnalysisV0unfold::CloneHistograms() {

	printf("Cloning histograms for analysis %s  \n", 
		this->GetName());
	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	hNt 			= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTrans")->Clone("hNt");
	hNtRec 			= (TH1F*)mFileMC->Get("MyAnalysisV0_0/hNchTransRC")->Clone("hNtRec");
	hNtGen			= (TH1F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMC")->Clone("hNtGen");
	hNtRM			= (TH2F*)mFileMC->Get("MyAnalysisV0_0/hNchTransRCvMC")->Clone("hNtRM");
	//hNtUnf 			= (TH1F*)hNt->Clone("hNtUnf");
	//hNtClosure 		= (TH1F*)hNt->Clone("hNtClosure");

	hNtMin 			= (TH1F*)((TH2F*)mHandler->analysis(0)->dirFile()->Get("hNtvNtMin"))->ProjectionX("hNtMin");
	hNtMinRec 		= (TH1F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMinRC")->Clone("hNtMinRec");
	hNtMinGen		= (TH1F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMinMC")->Clone("hNtMinGen");
	hNtMinRM		= (TH2F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMinRCvMC")->Clone("hNtMinRM");
	//hNtMinUnf 		= (TH1F*)hNtMin->Clone("hNtMinUnf");
	//hNtMinClosure 	= (TH1F*)hNtMin->Clone("hNtMinClosure");

	hNtMax 			= (TH1F*)((TH2F*)mHandler->analysis(0)->dirFile()->Get("hNtvNtMax"))->ProjectionX("hNtMax");
	hNtMaxRec 		= (TH1F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMaxRC")->Clone("hNtMaxRec");
	hNtMaxGen		= (TH1F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMaxMC")->Clone("hNtMaxGen");
	hNtMaxRM		= (TH2F*)mFileMC->Get("MyAnalysisV0_0/hNchTransMaxRCvMC")->Clone("hNtMaxRM");

	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNtFitCorrUnf[iSp][iType][iReg] = 0x0;
		hV0PtNtMinFitCorrUnf[iSp][iType][iReg] = 0x0;
		hV0PtNtMaxFitCorrUnf[iSp][iType][iReg] = 0x0;

		/*hV0PtRtFitCorrUnf[iSp][iType][iReg][iRtBin]		
		= (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]))->Clone(
			Form("hV0PtRtFitCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin)	);

		hV0PtRtMinFitCorrUnf[iSp][iType][iReg][iRtBin]		
		= (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]))->Clone(
			Form("hV0PtRtMinFitCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin)	);

		hV0PtRtMaxFitCorrUnf[iSp][iType][iReg][iRtBin]		
		= (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]))->Clone(
			Form("hV0PtRtMaxFitCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin)	);
		*/

	}	}	}


	return true;
}

Bool_t MyAnalysisV0unfold::SumLambdas() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNtFitCorr[2][iType][iReg]->Add(hV0PtNtFitCorr[3][iType][iReg],1.);
		hV0PtNtMinFitCorr[2][iType][iReg]->Add(hV0PtNtMinFitCorr[3][iType][iReg],1.);
		hV0PtNtMaxFitCorr[2][iType][iReg]->Add(hV0PtNtMaxFitCorr[3][iType][iReg],1.);


	}	}

	for (int iType = 2; iType < NTYPE; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
				
		hV0PtNt[2][iType][iReg]->Add(hV0PtNt[3][iType][iReg],1.);
		hV0PtNtMin[2][iType][iReg]->Add(hV0PtNtMin[3][iType][iReg],1.);
		hV0PtNtMax[2][iType][iReg]->Add(hV0PtNtMax[3][iType][iReg],1.);

	}	}


	return true;
}

void MyAnalysisV0unfold::SetMCInputFile(const Char_t *name) {

	TString fileName = TString(name);
	if (fileName.Data() != "") {
		mFileMC = new TFile(fileName,"READ");
		printf("MC File %s loaded in. \n", fileName.Data()); }
	else {
		printf("No MC file loaded.");
	}
}

Int_t MyAnalysisV0unfold::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	CloneHistograms();
	SumLambdas();

	mUnf = new UnfoldNTclass();
	cout << "Unfolding class created " << mUnf << endl;

	DoUnfoldingNt();
	DoUnfolding1D();

	delete mUnf; mUnf = new UnfoldNTclass();
	DoUnfoldingNtMin();
	DoUnfolding1DMin();

	delete mUnf; mUnf = new UnfoldNTclass();
	DoUnfoldingNtMax();
	DoUnfolding1DMax();

	BinHistogramsIntoRT();
	if (hV0PtRtSysSum[1][0][1]) ApplySystematics();
	
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
	}
		
	return 0;	
}

void MyAnalysisV0unfold::BinHistogramsIntoRT() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

		// REBIN
		int lowedge = mUnf->GetRTBin(iRtBin,kTRUE);
		int upedge  = mUnf->GetRTBin(iRtBin,kFALSE);
		Double_t NormScaleNt = hNtUnf->Integral(lowedge,upedge)/hNtUnf->Integral(1,51);

		cout << "BLEBLE " << hV0PtNtFitCorrUnf[iSp][iType][iReg] << " " << hV0PtNtMinFitCorrUnf[iSp][iType][iReg] << " " << hV0PtNtMaxFitCorrUnf[iSp][iType][iReg] << endl;

		if (hV0PtNtFitCorrUnf[iSp][iType][iReg]) {
			hV0PtRtFitCorrUnf[iSp][iType][iReg][iRtBin]
				= (TH1F*)hV0PtNtFitCorrUnf[iSp][iType][iReg]->ProjectionX(Form("hV0PtRtFitCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),lowedge,upedge);
			hV0PtRtFitCorrUnf[iSp][iType][iReg][iRtBin]->Scale(1./NormScaleNt);
		}


		lowedge = mUnf->GetRTMinBin(iRtBin,kTRUE);
		upedge  = mUnf->GetRTMinBin(iRtBin,kFALSE);
		Double_t NormScaleNtMin = hNtMinUnf->Integral(lowedge,upedge)/hNtMinUnf->Integral(1,51);

		if (hV0PtNtMinFitCorrUnf[iSp][iType][iReg]) {
			hV0PtRtMinFitCorrUnf[iSp][iType][iReg][iRtBin]
				= (TH1F*)hV0PtNtMinFitCorrUnf[iSp][iType][iReg]->ProjectionX(Form("hV0PtRtMinFitCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),lowedge,upedge);
			hV0PtRtMinFitCorrUnf[iSp][iType][iReg][iRtBin]->Scale(1./NormScaleNtMin);
		}
		
		lowedge = mUnf->GetRTMaxBin(iRtBin,kTRUE);
		upedge  = mUnf->GetRTMaxBin(iRtBin,kFALSE);
		Double_t NormScaleNtMax = hNtMaxUnf->Integral(lowedge,upedge)/hNtMaxUnf->Integral(1,51);

		if (hV0PtNtMaxFitCorrUnf[iSp][iType][iReg]) {
			hV0PtRtMaxFitCorrUnf[iSp][iType][iReg][iRtBin]
				= (TH1F*)hV0PtNtMaxFitCorrUnf[iSp][iType][iReg]->ProjectionX(Form("hV0PtRtMaxFitCorrUnf_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),lowedge,upedge);
			hV0PtRtMaxFitCorrUnf[iSp][iType][iReg][iRtBin]->Scale(1./NormScaleNtMax);
		}

	}	}	}	}

}

void MyAnalysisV0unfold::ApplySystematics() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	// SUM L+LBAR UNCERTAINTIES FIRST
	for (int iReg = 0; iReg < 3; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{


		for (Int_t iPt = 1; iPt < hV0PtRtSysSum[2][iReg][iRtBin]->GetNbinsX()+1; iPt++) {

			Double_t binc = hV0PtRtSysSum[2][iReg][iRtBin]->GetBinContent(iPt)*hV0PtRtSysSum[2][iReg][iRtBin]->GetBinContent(iPt);
			binc += hV0PtRtSysSum[3][iReg][iRtBin]->GetBinContent(iPt)*hV0PtRtSysSum[3][iReg][iRtBin]->GetBinContent(iPt);
			binc = TMath::Sqrt(binc)/2; // addition of uncorrelated uncertainties			
			hV0PtRtSysSum[2][iReg][iRtBin]->SetBinContent(iPt,binc);

			// and ratios to RT-int.
			binc = hV0PtRtSysSumUnc[2][iReg][iRtBin]->GetBinContent(iPt)*hV0PtRtSysSumUnc[2][iReg][iRtBin]->GetBinContent(iPt);
			binc += hV0PtRtSysSumUnc[3][iReg][iRtBin]->GetBinContent(iPt)*hV0PtRtSysSumUnc[3][iReg][iRtBin]->GetBinContent(iPt);
			binc = TMath::Sqrt(binc)/2; // addition of uncorrelated uncertainties			
			hV0PtRtSysSumUnc[2][iReg][iRtBin]->SetBinContent(iPt,binc);
		}

	}	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

		hV0PtRtFitCorrUnfSyst[iSp][iType][iReg][iRtBin] = 0x0;
		hV0PtRtMinFitCorrUnfSyst[iSp][iType][iReg][iRtBin] = 0x0;
		hV0PtRtMaxFitCorrUnfSyst[iSp][iType][iReg][iRtBin] = 0x0;
		hV0PtRtFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin] = 0x0;
		hV0PtRtMinFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin] = 0x0;
		hV0PtRtMaxFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin] = 0x0;


		hV0PtRtFitCorrUnfSyst[iSp][iType][iReg][iRtBin]
			= (TH1F*)hV0PtRtFitCorrUnf[iSp][iType][iReg][iRtBin]->Clone(
				Form("hV0PtRtFitCorrUnfSyst_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin));
		hV0PtRtFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]
			= (TH1F*)hV0PtRtFitCorrUnf[iSp][iType][iReg][iRtBin]->Clone(
				Form("hV0PtRtFitCorrUnfSystUnc_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin));

		if (iReg != 0 && iReg != 4) {
		hV0PtRtMinFitCorrUnfSyst[iSp][iType][iReg][iRtBin]
			= (TH1F*)hV0PtRtMinFitCorrUnf[iSp][iType][iReg][iRtBin]->Clone(
				Form("hV0PtRtMinFitCorrUnfSyst_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin));
		hV0PtRtMinFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]
			= (TH1F*)hV0PtRtMinFitCorrUnf[iSp][iType][iReg][iRtBin]->Clone(
				Form("hV0PtRtMinFitCorrUnfSystUnc_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin));			
		}

		if (iReg != 0 && iReg != 3) {
		hV0PtRtMaxFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]
			= (TH1F*)hV0PtRtMaxFitCorrUnf[iSp][iType][iReg][iRtBin]->Clone(
				Form("hV0PtRtMaxFitCorrUnfSystUnc_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin));
		hV0PtRtMaxFitCorrUnfSyst[iSp][iType][iReg][iRtBin]
			= (TH1F*)hV0PtRtMaxFitCorrUnf[iSp][iType][iReg][iRtBin]->Clone(
				Form("hV0PtRtMaxFitCorrUnfSyst_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin));
		}


		for (Int_t iPt = 1; iPt < hV0PtRtFitCorrUnfSyst[iSp][iType][iReg][iRtBin]->GetNbinsX()+1; iPt++) {
			hV0PtRtFitCorrUnfSyst[iSp][iType][iReg][iRtBin]->SetBinError(iPt,
				hV0PtRtFitCorrUnfSyst[iSp][iType][iReg][iRtBin]->GetBinContent(iPt)
				*hV0PtRtSysSum[iSp][iReg<3?iReg:0][iRtBin]->GetBinContent(iPt) );

			if (hV0PtRtMinFitCorrUnfSyst[iSp][iType][iReg][iRtBin]) hV0PtRtMinFitCorrUnfSyst[iSp][iType][iReg][iRtBin]->SetBinError(iPt,
				hV0PtRtMinFitCorrUnfSyst[iSp][iType][iReg][iRtBin]->GetBinContent(iPt)
				*hV0PtRtSysSum[iSp][iReg<3?iReg:0][iRtBin]->GetBinContent(iPt) );

			if (hV0PtRtMaxFitCorrUnfSyst[iSp][iType][iReg][iRtBin]) hV0PtRtMaxFitCorrUnfSyst[iSp][iType][iReg][iRtBin]->SetBinError(iPt,
				hV0PtRtMaxFitCorrUnfSyst[iSp][iType][iReg][iRtBin]->GetBinContent(iPt)
				*hV0PtRtSysSum[iSp][iReg<3?iReg:0][iRtBin]->GetBinContent(iPt) );

			hV0PtRtFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]->SetBinError(iPt,
				hV0PtRtFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]->GetBinContent(iPt)
				*hV0PtRtSysSumUnc[iSp][iReg<3?iReg:0][iRtBin]->GetBinContent(iPt) );

			if (hV0PtRtMinFitCorrUnfSyst[iSp][iType][iReg][iRtBin]) hV0PtRtMinFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]->SetBinError(iPt,
				hV0PtRtMinFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]->GetBinContent(iPt)
				*hV0PtRtSysSumUnc[iSp][iReg<3?iReg:0][iRtBin]->GetBinContent(iPt) );

			if (hV0PtRtMaxFitCorrUnfSyst[iSp][iType][iReg][iRtBin]) hV0PtRtMaxFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]->SetBinError(iPt,
				hV0PtRtMaxFitCorrUnfSystUnc[iSp][iType][iReg][iRtBin]->GetBinContent(iPt)
				*hV0PtRtSysSumUnc[iSp][iReg<3?iReg:0][iRtBin]->GetBinContent(iPt) );
		}

	}	}	}	}

}


void MyAnalysisV0unfold::DoUnfoldingNt() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding";
	const bool eRM = false;
	const int NumberOfIters = 5;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtRM = FlipMatrix(hNtRM);
	// ROWWISE NORMALISATION IS PERFORMED IN mUnf
	//cout << "Histograms: " << hNtRec << " " << hNt << " " << hNtGen << " " << hNtRM << "\n";
	
	if (eRM) mUnf->ExtrapolateRM(hNtRM);
	mUnf->SetnIter(NumberOfIters);
	mUnf->SetError(hNtRec);
	mUnf->SetError(hNtGen);
	mUnf->SaveSolutionNT(kTRUE);
	
	if (mHandler->GetFlagMC()) mUnf->Setup(hNtRec,hNtGen,hNtRM);
	else mUnf->Setup(hNt,hNtGen,hNtRM);
	
	mUnf->Unfold();
	mUnf->V2H();
	
	printf(" - Unfolding region : %s\n",mUnf->GetRegion());

	hNtUnf = (TH1F*)(mUnf->GetUnfoldedDistH())->Clone("hNtUnf");
	if(mHandler->GetFlagMC())	hNtClosure = (TH1F*)(mUnf->GetClosureH())->Clone("hNtClosure");
	
	TH1F* hNT = (TH1F*)hNtUnf->Clone("_hNT");
	hRtUnf = (TH1F*)mUnf->RebinNT2RT(hNtUnf, kTRUE);
	hRtUnf->Scale(1.0/hRtUnf->Integral());
	
	//! Relative statistical uncertainty
	//! NT distribution with final bins to plot
	TH1F* hRelStatUnc = (TH1F*)hNtUnf->Clone("hRelStatUnc_Nt");
	hRelStatUnc->Reset();

	TH1F* hNch = (TH1F*)hNtUnf->Clone("hNT");
	hNch->Reset();

	for(int bin = 1; bin <= hRelStatUnc->GetNbinsX(); bin++){
		double yield = hNtUnf->GetBinContent(bin);
		double error = hNtUnf->GetBinError(bin);
		if( yield > 0. ) hRelStatUnc->SetBinContent(bin, error / yield);
	}

	lOut->Add(hNtRM);
	
	if (mHandler->GetFlagMC()) 	lOut->Add(hNtClosure);
	//else 	ComparisonPublished(hRtUnf);

	lOut->Add(hNtUnf);
	lOut->Add(hRtUnf);
	lOut->Add(hRelStatUnc);

	TGraph* gChi2 = nullptr;
	if (mHandler->GetFlagMC())	{
		mUnf->DrawNchClosure(hNtGen,hNtUnf,-0.5,30.0,"#it{N}_{T}","Unfolded/True",dOut);
		TVectorD iterX = (TVectorD)mUnf->GetChi2Vectors(kTRUE);
		TVectorD Chi2 = (TVectorD)mUnf->GetChi2Vectors(kFALSE);
		gChi2 = new TGraph(iterX,Chi2);
		gChi2->SetName("gChi2");
		gChi2->SetMarkerStyle(8);
		lOut->Add(gChi2);
	}

	TFile* fUnfOut;
	if (mHandler->GetFlagMC())	fUnfOut = new TFile(Form("./%s/1D_newClass_mc.root",dOut),"RECREATE");
	else	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");


	fUnfOut->cd();
	lOut->Write();
	(mUnf->GetObjArray())->Write();
	fUnfOut->Close();
	delete fUnfOut;

	mDirFile->cd();
	lOut->Write();
}

void MyAnalysisV0unfold::DoUnfoldingNtMin() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding_min";
	const bool eRM = false;
	const int NumberOfIters = 5;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtMinRM = FlipMatrix(hNtMinRM);
	// ROWWISE NORMALISATION IS PERFORMED IN mUnf
	//cout << "Histograms: " << hNtRec << " " << hNt << " " << hNtGen << " " << hNtRM << "\n";
	
	if (eRM) mUnf->ExtrapolateRM(hNtMinRM);
	mUnf->SetnIter(NumberOfIters);
	mUnf->SetError(hNtMinRec);
	mUnf->SetError(hNtMinGen);
	mUnf->SaveSolutionNT(kTRUE);
	
	if (mHandler->GetFlagMC()) mUnf->Setup(hNtMinRec,hNtMinGen,hNtMinRM);
	else mUnf->Setup(hNtMin,hNtMinGen,hNtMinRM);
	
	mUnf->Unfold();
	mUnf->V2H();
	
	printf(" - Unfolding region : %s_min\n",mUnf->GetRegion());

	hNtMinUnf = (TH1F*)(mUnf->GetUnfoldedDistH())->Clone("hNtMinUnf");
	if(mHandler->GetFlagMC())	hNtMinClosure = (TH1F*)(mUnf->GetClosureH())->Clone("hNtMinClosure");
	
	TH1F* hNTMin = (TH1F*)hNtMinUnf->Clone("_hNTMin");
	hRtMinUnf = (TH1F*)mUnf->RebinNT2RT(hNtMinUnf, kTRUE);
	hRtMinUnf->Scale(1.0/hRtMinUnf->Integral());
	
	//! Relative statistical uncertainty
	//! NT distribution with final bins to plot
	TH1F* hRelStatUnc = (TH1F*)hNtUnf->Clone("hRelStatUnc_NtMin");
	hRelStatUnc->Reset();

	TH1F* hNch = (TH1F*)hNtMinUnf->Clone("hNTMin");
	hNch->Reset();

	for(int bin = 1; bin <= hRelStatUnc->GetNbinsX(); bin++){
		double yield = hNtMinUnf->GetBinContent(bin);
		double error = hNtMinUnf->GetBinError(bin);
		if( yield > 0. ) hRelStatUnc->SetBinContent(bin, error / yield);
	}

	lOut->Add(hNtMinRM);
	
	if (mHandler->GetFlagMC()) 	lOut->Add(hNtMinClosure);
	//else 	ComparisonPublished(hRtUnf);

	lOut->Add(hNtMinUnf);
	lOut->Add(hRtMinUnf);
	lOut->Add(hRelStatUnc);

	TGraph* gChi2 = nullptr;
	if (mHandler->GetFlagMC())	{
		mUnf->DrawNchClosure(hNtMinGen,hNtMinUnf,-0.5,15.5,"#it{N}_{T,min.}","Unfolded/True",dOut);
		TVectorD iterX = (TVectorD)mUnf->GetChi2Vectors(kTRUE);
		TVectorD Chi2 = (TVectorD)mUnf->GetChi2Vectors(kFALSE);
		gChi2 = new TGraph(iterX,Chi2);
		gChi2->SetName("gChi2");
		gChi2->SetMarkerStyle(8);
		lOut->Add(gChi2);
	}

	TFile* fUnfOut;
	if (mHandler->GetFlagMC())	fUnfOut = new TFile(Form("./%s/1D_newClass_mc.root",dOut),"RECREATE");
	else	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");


	fUnfOut->cd();
	lOut->Write();
	(mUnf->GetObjArray())->Write();
	fUnfOut->Close();
	delete fUnfOut;

	mDirFile->cd();
	lOut->Write();
}

void MyAnalysisV0unfold::DoUnfoldingNtMax() {

	TList* lOut = new TList();
	lOut->SetOwner();

	// Output plots
	const char* dOut = "results_unfolding_max";
	const bool eRM = false;
	const int NumberOfIters = 5;

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	hNtMaxRM = FlipMatrix(hNtMaxRM);
	// ROWWISE NORMALISATION IS PERFORMED IN mUnf
	//cout << "Histograms: " << hNtRec << " " << hNt << " " << hNtGen << " " << hNtRM << "\n";
	
	if (eRM) mUnf->ExtrapolateRM(hNtMaxRM);
	mUnf->SetnIter(NumberOfIters);
	mUnf->SetError(hNtMaxRec);
	mUnf->SetError(hNtMaxGen);
	mUnf->SaveSolutionNT(kTRUE);
	
	if (mHandler->GetFlagMC()) mUnf->Setup(hNtMaxRec,hNtMaxGen,hNtMaxRM);
	else mUnf->Setup(hNtMax,hNtMaxGen,hNtMaxRM);
	
	mUnf->Unfold();
	mUnf->V2H();
	
	printf(" - Unfolding region : %s_max\n",mUnf->GetRegion());

	hNtMaxUnf = (TH1F*)(mUnf->GetUnfoldedDistH())->Clone("hNtMaxUnf");
	if(mHandler->GetFlagMC())	hNtMaxClosure = (TH1F*)(mUnf->GetClosureH())->Clone("hNtMaxClosure");
	
	TH1F* hNTMax = (TH1F*)hNtMaxUnf->Clone("_hNTMax");
	hRtMaxUnf = (TH1F*)mUnf->RebinNT2RT(hNtMaxUnf, kTRUE);
	hRtMaxUnf->Scale(1.0/hRtMaxUnf->Integral());
	
	//! Relative statistical uncertainty
	//! NT distribution with final bins to plot
	TH1F* hRelStatUnc = (TH1F*)hNtUnf->Clone("hRelStatUnc_NtMax");
	hRelStatUnc->Reset();

	TH1F* hNch = (TH1F*)hNtMaxUnf->Clone("hNTMax");
	hNch->Reset();

	for(int bin = 1; bin <= hRelStatUnc->GetNbinsX(); bin++){
		double yield = hNtMaxUnf->GetBinContent(bin);
		double error = hNtMaxUnf->GetBinError(bin);
		if( yield > 0. ) hRelStatUnc->SetBinContent(bin, error / yield);
	}

	lOut->Add(hNtMaxRM);
	
	if (mHandler->GetFlagMC()) 	lOut->Add(hNtMaxClosure);
	//else 	ComparisonPublished(hRtUnf);

	lOut->Add(hNtMaxUnf);
	lOut->Add(hRtMaxUnf);
	lOut->Add(hRelStatUnc);

	TGraph* gChi2 = nullptr;
	if (mHandler->GetFlagMC())	{
		mUnf->DrawNchClosure(hNtMaxGen,hNtMaxUnf,-0.5,25.5,"#it{N}_{T,max.}","Unfolded/True",dOut);
		TVectorD iterX = (TVectorD)mUnf->GetChi2Vectors(kTRUE);
		TVectorD Chi2 = (TVectorD)mUnf->GetChi2Vectors(kFALSE);
		gChi2 = new TGraph(iterX,Chi2);
		gChi2->SetName("gChi2");
		gChi2->SetMarkerStyle(8);
		lOut->Add(gChi2);
	}

	TFile* fUnfOut;
	if (mHandler->GetFlagMC())	fUnfOut = new TFile(Form("./%s/1D_newClass_mc.root",dOut),"RECREATE");
	else	fUnfOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");


	fUnfOut->cd();
	lOut->Write();
	(mUnf->GetObjArray())->Write();
	fUnfOut->Close();
	delete fUnfOut;

	mDirFile->cd();
	lOut->Write();
}

void MyAnalysisV0unfold::DoUnfolding1D() {

	mHandler->root()->SetBatch(kTRUE);

	enum { D, RC, MC };
	const char* dOut = "results_unfolding";
	const bool eRM = false;
	const int NumberOfIters = 5;
	const char* Regions[NREGIONS] = {"Trans1D","Toward","Away","TransMin1D","TransMax1D"};

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	Int_t iType = mHandler->GetFlagMC() ? RC : D;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(3+iSp);
		obj->SetMCAnalysis(mHandler->GetFlagMC());
		obj->SaveSolutionNT(kFALSE);

		printf(" - Unfolding species %s vs NT in region : %s\n",SPECIES[iSp],obj->GetRegion());
		cout << "RM has " << hNtRM->GetNbinsX() << " x " << hNtRM->GetNbinsY() << endl;
		cout << "hMC has " << hV0PtNtFitCorr[iSp][iType][iReg]->GetNbinsX() << endl;
		obj->UnfoldV02D(hNtRM,hV0PtNtFitCorr[iSp][iType][iReg],hV0PtNtFitCorr[iSp][iType][iReg]);
		
		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << hV0PtNtFitCorr[iSp][iType][iReg] << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		// first is Gen, second is Unf

		// APPLY PROPER ERRORS!
		hV0PtNtFitCorrUnf[iSp][iType][iReg] = (TH2F*)Arr->FindObject("_hPtvsNch")->Clone(Form("hV0PtNtFitCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		for (int iX = 1; iX < hV0PtNtFitCorrUnf[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNtFitCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float bine = hV0PtNtFitCorr[iSp][iType][iReg]->GetBinError(iX,iY);
				if (bine>0) hV0PtNtFitCorrUnf[iSp][iType][iReg]->SetBinError(iX,iY,bine);
			}
		}

		if (mHandler->GetFlagMC())	obj->GetMCclosureinRTBins(hV0PtNt[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	

		// NORMALISE
		/*for (int iNt = 1; iNt < hV0PtNtFitCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {
			Double_t NormEv = 0;
			NormEv = hNtUnf->Integral(1,50) > 0 ? (double)hNtUnf->GetBinContent(iNt) / hNtUnf->Integral(1,50) : 1;
			//cout << "Normalisingggg by " << NormEv << endl;
			TH2F* htmp = hV0PtNtFitCorrUnf[iSp][iType][iReg];

			TH1F* hpt = (TH1F*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
			delete hpt;
		}*/

		lOut->Add(hV0PtNtFitCorrUnf[iSp][iType][iReg]);
		fUnfOut1D->cd();
		TDirectory* dir = fUnfOut1D->mkdir(Form("%s_%s",SPECIES[iSp],REGIONS[iReg]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	}

	lOut->Write();
	fUnfOut1D->Close();
	delete fUnfOut1D;
	mDirFile->cd();
	lOut->Write();

	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0unfold::DoUnfolding1DMin() {

	mHandler->root()->SetBatch(kTRUE);

	enum { D, RC, MC };
	const char* dOut = "results_unfolding_min";
	const bool eRM = false;
	const int NumberOfIters = 5;
	const char* Regions[4] = {"Trans1D","Toward","Away","TransMin1D"};
	

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 1; iReg < 4; ++iReg)		{
	Int_t iType = mHandler->GetFlagMC() ? RC : D;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(3+iSp);
		obj->SetMCAnalysis(mHandler->GetFlagMC());
		obj->SaveSolutionNT(kTRUE);

		printf(" - Unfolding species %s in region : %s\n",SPECIES[iSp],obj->GetRegion());
		cout << "RM has " << hNtMinRM->GetNbinsX() << " x " << hNtMinRM->GetNbinsY() << endl;
		cout << "hMC has " << hV0PtNtMinFitCorr[iSp][iType][iReg]->GetNbinsX() << endl;
		obj->UnfoldV02D(hNtMinRM,hV0PtNtMinFitCorr[iSp][iType][iReg],hV0PtNtMinFitCorr[iSp][iType][iReg]);
		
		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << hV0PtNtMinFitCorr[iSp][iType][iReg] << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		// first is Gen, second is Unf

		// APPLY PROPER ERRORS!
		hV0PtNtMinFitCorrUnf[iSp][iType][iReg] = (TH2F*)Arr->FindObject("_hPtvsNch")->Clone(Form("hV0PtNtMinFitCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		for (int iX = 1; iX < hV0PtNtMinFitCorrUnf[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNtMinFitCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float bine = hV0PtNtMinFitCorr[iSp][iType][iReg]->GetBinError(iX,iY);
				if (bine>0) hV0PtNtMinFitCorrUnf[iSp][iType][iReg]->SetBinError(iX,iY,bine);
			}
		}

		if (mHandler->GetFlagMC())	obj->GetMCclosureinRTMinBins(hV0PtNtMin[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	

		// NORMALISE
		/*for (int iNt = 1; iNt < hV0PtNtMinFitCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {
			Double_t NormEv = 0;
			NormEv = hNtMinUnf->Integral(1,50) > 0 ? (double)hNtMinUnf->GetBinContent(iNt) / hNtMinUnf->Integral(1,50) : 1;
			//cout << "Normalisingggg by " << NormEv << endl;
			TH2F* htmp = hV0PtNtMinFitCorrUnf[iSp][iType][iReg];

			TH1F* hpt = (TH1F*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
			delete hpt;
		}*/

		lOut->Add(hV0PtNtMinFitCorrUnf[iSp][iType][iReg]);
		fUnfOut1D->cd();
		TDirectory* dir = fUnfOut1D->mkdir(Form("%s_%s",SPECIES[iSp],REGIONS[iReg]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	}

	lOut->Write();
	fUnfOut1D->Close();
	delete fUnfOut1D;
	mDirFile->cd();
	lOut->Write();

	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0unfold::DoUnfolding1DMax() {

	mHandler->root()->SetBatch(kTRUE);

	enum { D, RC, MC };
	const char* dOut = "results_unfolding_max";
	const bool eRM = false;
	const int NumberOfIters = 5;
	const char* Regions[5] = {"Trans1D","Toward","Away","TransMin1D","TransMax1D"};
	

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 1; iReg < 5; ++iReg)		{
	if (iReg == 3) continue;
	Int_t iType = mHandler->GetFlagMC() ? RC : D;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(3+iSp);
		obj->SetMCAnalysis(mHandler->GetFlagMC());
		obj->SaveSolutionNT(kTRUE);

		printf(" - Unfolding species %s in region : %s\n",SPECIES[iSp],obj->GetRegion());
		cout << "RM has " << hNtMaxRM->GetNbinsX() << " x " << hNtMaxRM->GetNbinsY() << endl;
		cout << "hMC has " << hV0PtNtMaxFitCorr[iSp][iType][iReg]->GetNbinsX() << endl;
		obj->UnfoldV02D(hNtMaxRM,hV0PtNtMaxFitCorr[iSp][iType][iReg],hV0PtNtMaxFitCorr[iSp][iType][iReg]);
		
		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << hV0PtNtMaxFitCorr[iSp][iType][iReg] << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		// first is Gen, second is Unf

		// APPLY PROPER ERRORS!
		hV0PtNtMaxFitCorrUnf[iSp][iType][iReg] = (TH2F*)Arr->FindObject("_hPtvsNch")->Clone(Form("hV0PtNtMaxFitCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		for (int iX = 1; iX < hV0PtNtMaxFitCorrUnf[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNtMaxFitCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float bine = hV0PtNtMaxFitCorr[iSp][iType][iReg]->GetBinError(iX,iY);
				if (bine>0) hV0PtNtMaxFitCorrUnf[iSp][iType][iReg]->SetBinError(iX,iY,bine);
			}
		}

		if (mHandler->GetFlagMC())	obj->GetMCclosureinRTMaxBins(hV0PtNtMax[iSp][MC][iReg],(TH2F*)Arr->FindObject("_hPtvsNch"));	

		// NORMALISE
		/*for (int iNt = 1; iNt < hV0PtNtMinFitCorrUnf[iSp][iType][iReg]->GetNbinsY()+1; iNt++) {
			Double_t NormEv = 0;
			NormEv = hNtMinUnf->Integral(1,50) > 0 ? (double)hNtMinUnf->GetBinContent(iNt) / hNtMinUnf->Integral(1,50) : 1;
			//cout << "Normalisingggg by " << NormEv << endl;
			TH2F* htmp = hV0PtNtMinFitCorrUnf[iSp][iType][iReg];

			TH1F* hpt = (TH1F*)htmp->ProjectionX(Form("pt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iNt),iNt,iNt);
			if (NormEv>0) hpt->Scale( (iNt<51) ? 1./NormEv : 0 );
			for (int iBin = 1; iBin < htmp->GetNbinsX()+1; iBin++) {
				htmp->SetBinContent(iBin,iNt,hpt->GetBinContent(iBin));
				htmp->SetBinError(iBin,iNt,hpt->GetBinError(iBin));
			}
			delete hpt;
		}*/

		lOut->Add(hV0PtNtMaxFitCorrUnf[iSp][iType][iReg]);
		fUnfOut1D->cd();
		TDirectory* dir = fUnfOut1D->mkdir(Form("%s_%s",SPECIES[iSp],REGIONS[iReg]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	}

	lOut->Write();
	fUnfOut1D->Close();
	delete fUnfOut1D;
	mDirFile->cd();
	lOut->Write();

	mHandler->root()->SetBatch(kFALSE);
}


TH2F* MyAnalysisV0unfold::FlipMatrix(TH2F* h) {

	TH2F* htmp = (TH2F*)h->Clone(Form("%s_flip",h->GetName()));
	for (int i = 1; i < h->GetNbinsX()+1; i++) {
	for (int j = 1; j < h->GetNbinsY()+1; j++) {
			htmp->SetBinContent(i,j,h->GetBinContent(j,i));
			htmp->SetBinError(i,j,h->GetBinError(j,i));
	}	}
	htmp->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(h->GetXaxis()->GetTitle());
	delete h;
	return htmp; 
}

void MyAnalysisV0unfold::ComparisonPublished(TH1F* hRT)	{

	TFile* fIn = new TFile("./results/HEPData-ins1762350-v1-Probability_distribution_as_function_of_R_T.root","READ");
	TDirectory* dIn = (TDirectory*)fIn->Get("Probability distribution as function of R_T");
	if (!dIn) printf("Not found: Probability distribution as function of R_T\n");

	TH1F* hstat = (TH1F*)dIn->Get("Hist1D_y1");
	TH1F* he = (TH1F*)dIn->Get("Hist1D_y1_e1");

	for(int bin = 1; bin <= hstat->GetNbinsX(); bin++)
		hstat->SetBinError(bin,he->GetBinError(bin));

	for(int bin = 1; bin <= hRT->GetNbinsX(); bin++)
		hRT->SetBinContent(bin,hRT->GetBinContent(bin)/hRT->GetBinWidth(bin));

	hRT->Divide(hRT,hstat,1.0,1.0,"");

}