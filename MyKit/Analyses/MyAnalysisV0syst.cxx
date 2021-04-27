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

#include "MyAnalysisV0extract.h"
#include "MyAnalysisV0syst.h"
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

ClassImp(MyAnalysisV0syst)



MyAnalysisV0syst::MyAnalysisV0syst() {

}

Int_t MyAnalysisV0syst::Init() {

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

Int_t MyAnalysisV0syst::Make(Int_t iEv) {
	return 0;
}


Bool_t MyAnalysisV0syst::BorrowHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	//for (int iMu = 0; iMu < NMULTI; ++iMu)		{
		hV0IMvRadiusL[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvRadiusL_%s",SPECIES[iSp]));
		hV0IMvDCAdd[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvDCAdd_%s",SPECIES[iSp]));
		hV0IMvCPA[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvCPA_%s",SPECIES[iSp]));
		hV0IMvFastSignal[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvFastSignal_%s",SPECIES[iSp]));
		hV0IMvCompMass[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvCompMass_%s",SPECIES[iSp]));
		if (iSp>1) hV0IMvLifetime[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvLifetime_%s",SPECIES[iSp]));
		hV0IMvNSigmaTPC[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvNSigmaTPC_%s",SPECIES[iSp]));
		hV0IMvDCAPVpos[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvDCAPVpos_%s",SPECIES[iSp]));
		hV0IMvDCAPVneg[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvDCAPVneg_%s",SPECIES[iSp]));
		hV0IMvNCluster[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvNCluster_%s",SPECIES[iSp]));
		hV0IMvNClusterF[iSp]			
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvNClusterF_%s",SPECIES[iSp]));
	}

	// LOAD FROM MC
	if (!mFileMC) {
		printf("No MC file loaded! Comparisons of systematics with MC not performed.\n");
		return 0;
	}

	TDirectoryFile* dirFile1 = new TDirectoryFile("mcFile","mcFile","",mHandler->file());
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("TDirectoryFile")) {
		cout << "Doing systematics from a TDirectoryFile" << endl;
		dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");}
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("THashList")) {
		cout << "Doing systematics from a THashList" << endl;
		THashList* hashList = (THashList*)mFileMC->Get("MyAnalysisV0_0");
		while (hashList->GetEntries()) {
			dirFile1->Append(hashList->First());
			hashList->RemoveFirst();
		}
	}
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	//for (int iMu = 0; iMu < NMULTI; ++iMu)		{
		hMCV0IMvRadiusL[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvRadiusL_%s",SPECIES[iSp]));
		hMCV0IMvDCAdd[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvDCAdd_%s",SPECIES[iSp]));
		hMCV0IMvCPA[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvCPA_%s",SPECIES[iSp]));
		hMCV0IMvFastSignal[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvFastSignal_%s",SPECIES[iSp]));
		hMCV0IMvCompMass[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvCompMass_%s",SPECIES[iSp]));
		if (iSp>1) hMCV0IMvLifetime[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvLifetime_%s",SPECIES[iSp]));
		hMCV0IMvNSigmaTPC[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvNSigmaTPC_%s",SPECIES[iSp]));
		hMCV0IMvDCAPVpos[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvDCAPVpos_%s",SPECIES[iSp]));
		hMCV0IMvDCAPVneg[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvDCAPVneg_%s",SPECIES[iSp]));
		hMCV0IMvNCluster[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvNCluster_%s",SPECIES[iSp]));
		hMCV0IMvNClusterF[iSp]			
			= (TH2D*)dirFile1->Get(Form("hV0IMvNClusterF_%s",SPECIES[iSp]));
	} 

}

Bool_t MyAnalysisV0syst::CreateHistograms() {

	// RAW YIELD FRACTION LOSS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0YieldvRadiusL[iSp] = (TH1D*)hV0IMvRadiusL[iSp]->ProjectionX(
			Form("hV0YieldvRadiusL_%s",SPECIES[iSp]),0,-1);
		hV0YieldvDCAdd[iSp] = (TH1D*)hV0IMvDCAdd[iSp]->ProjectionX(
			Form("hV0YieldvDCAdd_%s",SPECIES[iSp]),0,-1);
		hV0YieldvCPA[iSp] = (TH1D*)hV0IMvCPA[iSp]->ProjectionX(
			Form("hV0YieldvCPA_%s",SPECIES[iSp]),0,-1);
		hV0YieldvFastSignal[iSp] = (TH1D*)hV0IMvFastSignal[iSp]->ProjectionX(
			Form("hV0YieldvFastSignal_%s",SPECIES[iSp]),0,-1);
		hV0YieldvCompMass[iSp] = (TH1D*)hV0IMvCompMass[iSp]->ProjectionX(
			Form("hV0YieldvCompMass_%s",SPECIES[iSp]),0,-1);
		if (iSp>1) hV0YieldvLifetime[iSp] = (TH1D*)hV0IMvLifetime[iSp]->ProjectionX(
			Form("hV0YieldvLifetime_%s",SPECIES[iSp]),0,-1);
		hV0YieldvNSigmaTPC[iSp] = (TH1D*)hV0IMvNSigmaTPC[iSp]->ProjectionX(
			Form("hV0YieldvNSigmaTPC_%s",SPECIES[iSp]),0,-1);
		hV0YieldvDCAPVpos[iSp] = (TH1D*)hV0IMvDCAPVpos[iSp]->ProjectionX(
			Form("hV0YieldvDCAPVpos_%s",SPECIES[iSp]),0,-1);
		hV0YieldvDCAPVneg[iSp] = (TH1D*)hV0IMvDCAPVneg[iSp]->ProjectionX(
			Form("hV0YieldvDCAPVneg_%s",SPECIES[iSp]),0,-1);
		hV0YieldvNCluster[iSp] = (TH1D*)hV0IMvNCluster[iSp]->ProjectionX(
			Form("hV0YieldvNCluster_%s",SPECIES[iSp]),0,-1);
		hV0YieldvNClusterF[iSp] = (TH1D*)hV0IMvNClusterF[iSp]->ProjectionX(
			Form("hV0YieldvNClusterF_%s",SPECIES[iSp]),0,-1);
	}

	// ROGER BARLOW CHECKS


	// MC

	if (!mFileMC) {
		printf("No MC file loaded! Comparisons of systematics with MC not performed.\n");
		return 0;
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hMCV0YieldvRadiusL[iSp] = (TH1D*)hMCV0IMvRadiusL[iSp]->ProjectionX(
			Form("hMCV0YieldvRadiusL_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvDCAdd[iSp] = (TH1D*)hMCV0IMvDCAdd[iSp]->ProjectionX(
			Form("hMCV0YieldvDCAdd_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvCPA[iSp] = (TH1D*)hMCV0IMvCPA[iSp]->ProjectionX(
			Form("hMCV0YieldvCPA_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvFastSignal[iSp] = (TH1D*)hMCV0IMvFastSignal[iSp]->ProjectionX(
			Form("hMCV0YieldvFastSignal_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvCompMass[iSp] = (TH1D*)hMCV0IMvCompMass[iSp]->ProjectionX(
			Form("hMCV0YieldvCompMass_%s",SPECIES[iSp]),0,-1);
		if (iSp>1) hMCV0YieldvLifetime[iSp] = (TH1D*)hMCV0IMvLifetime[iSp]->ProjectionX(
			Form("hMCV0YieldvLifetime_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvNSigmaTPC[iSp] = (TH1D*)hMCV0IMvNSigmaTPC[iSp]->ProjectionX(
			Form("hMCV0YieldvNSigmaTPC_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvDCAPVpos[iSp] = (TH1D*)hMCV0IMvDCAPVpos[iSp]->ProjectionX(
			Form("hMCV0YieldvDCAPVpos_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvDCAPVneg[iSp] = (TH1D*)hMCV0IMvDCAPVneg[iSp]->ProjectionX(
			Form("hMCV0YieldvDCAPVneg_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvNCluster[iSp] = (TH1D*)hMCV0IMvNCluster[iSp]->ProjectionX(
			Form("hMCV0YieldvNCluster_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvNClusterF[iSp] = (TH1D*)hMCV0IMvNClusterF[iSp]->ProjectionX(
			Form("hMCV0YieldvNClusterF_%s",SPECIES[iSp]),0,-1);
	}
	
}


Bool_t MyAnalysisV0syst::CloneHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;


}

void MyAnalysisV0syst::SetMCInputFile(const Char_t *name) {

	TString fileName = TString(name);
	if (fileName.Data() != "") {
		mFileMC = new TFile(fileName,"READ");
		printf("MC File %s loaded in. \n", fileName.Data()); }
	else {
		printf("No MC file loaded.");
	}
}

Int_t MyAnalysisV0syst::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	CloneHistograms();
	
	StudyRawYieldLoss();
	
	return 0;	
}

void MyAnalysisV0syst::StudyRawYieldLoss() {

	enum { rising, sinking };
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		ProcessRawYieldLossHist(hV0IMvRadiusL[iSp],hV0YieldvRadiusL[iSp],iSp,0.3,sinking);
		ProcessRawYieldLossHist(hV0IMvDCAdd[iSp],hV0YieldvDCAdd[iSp],iSp,1.5,rising);
		ProcessRawYieldLossHist(hV0IMvCPA[iSp],hV0YieldvCPA[iSp],iSp,(iSp>1)?0.993:0.95,sinking);
		ProcessRawYieldLossHist(hV0IMvFastSignal[iSp],hV0YieldvFastSignal[iSp],iSp,1.,sinking);
		ProcessRawYieldLossHist(hV0IMvCompMass[iSp],hV0YieldvCompMass[iSp],iSp,2.5,sinking);
		if (iSp>1) ProcessRawYieldLossHist(hV0IMvLifetime[iSp],hV0YieldvLifetime[iSp],iSp,40.,rising);
		ProcessRawYieldLossHist(hV0IMvNSigmaTPC[iSp],hV0YieldvNSigmaTPC[iSp],iSp,6.5,rising);
		ProcessRawYieldLossHist(hV0IMvDCAPVpos[iSp],hV0YieldvDCAPVpos[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hV0IMvDCAPVneg[iSp],hV0YieldvDCAPVneg[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hV0IMvNCluster[iSp],hV0YieldvNCluster[iSp],iSp,70.,sinking);
		ProcessRawYieldLossHist(hV0IMvNClusterF[iSp],hV0YieldvNClusterF[iSp],iSp,0.8,sinking);
	}

	if (!mFileMC) {
		printf("No MC file loaded! Comparisons of systematics with MC not performed.\n");
		return;
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		ProcessRawYieldLossHist(hMCV0IMvRadiusL[iSp],hMCV0YieldvRadiusL[iSp],iSp,0.3,sinking);
		ProcessRawYieldLossHist(hMCV0IMvDCAdd[iSp],hMCV0YieldvDCAdd[iSp],iSp,1.5,rising);
		ProcessRawYieldLossHist(hMCV0IMvCPA[iSp],hMCV0YieldvCPA[iSp],iSp,(iSp>1)?0.993:0.95,sinking);
		ProcessRawYieldLossHist(hMCV0IMvFastSignal[iSp],hMCV0YieldvFastSignal[iSp],iSp,1.,sinking);
		ProcessRawYieldLossHist(hMCV0IMvCompMass[iSp],hMCV0YieldvCompMass[iSp],iSp,2.5,sinking);
		if (iSp>1) ProcessRawYieldLossHist(hMCV0IMvLifetime[iSp],hMCV0YieldvLifetime[iSp],iSp,40.,rising);
		ProcessRawYieldLossHist(hMCV0IMvNSigmaTPC[iSp],hMCV0YieldvNSigmaTPC[iSp],iSp,6.5,rising);
		ProcessRawYieldLossHist(hMCV0IMvDCAPVpos[iSp],hMCV0YieldvDCAPVpos[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hMCV0IMvDCAPVneg[iSp],hMCV0YieldvDCAPVneg[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hMCV0IMvNCluster[iSp],hMCV0YieldvNCluster[iSp],iSp,70.,sinking);
		ProcessRawYieldLossHist(hMCV0IMvNClusterF[iSp],hMCV0YieldvNClusterF[iSp],iSp,0.8,sinking);
	}

	mHandler->root()->SetBatch(kTRUE);
	TCanvas* cRYL[NSPECIES]; 
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		cRYL[iSp] = new TCanvas(Form("cRYL_%s",SPECIES[iSp]), "", 2800, 2000);
		cRYL[iSp]->Divide(4,3,5e-5,5e-5);
		Int_t padC = 2; cRYL[iSp]->cd(padC);

		DrawRawYieldLossHist(hV0YieldvRadiusL[iSp],hMCV0YieldvRadiusL[iSp],(iSp>1)?0.05:0.15); 
		//DrawVariation(0.3,kRed,cRYL[iSp]->GetPad(padC)); DrawVariation(0.4,kBlue,cRYL[iSp]->GetPad(padC)); DrawVariation(0.5,kBlack,cRYL[iSp]->GetPad(padC)); DrawVariation(0.6,kGreen+2,cRYL[iSp]->GetPad(padC));
		cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvDCAdd[iSp],hMCV0YieldvDCAdd[iSp],1.0); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvCPA[iSp],hMCV0YieldvCPA[iSp],0.2); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvFastSignal[iSp],hMCV0YieldvFastSignal[iSp],0.85); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvCompMass[iSp],hMCV0YieldvCompMass[iSp],(iSp>1)?0.20:0.15); cRYL[iSp]->cd(padC++);
		if (iSp>1) DrawRawYieldLossHist(hV0YieldvLifetime[iSp],hMCV0YieldvLifetime[iSp],(iSp>1)?1.:0.15); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvNSigmaTPC[iSp],hMCV0YieldvNSigmaTPC[iSp],0.6); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvDCAPVpos[iSp],hMCV0YieldvDCAPVpos[iSp],(iSp==2)?0.15:0.04); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvDCAPVneg[iSp],hMCV0YieldvDCAPVneg[iSp],(iSp==3)?0.15:0.04); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvNCluster[iSp],hMCV0YieldvNCluster[iSp],0.03); cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvNClusterF[iSp],hMCV0YieldvNClusterF[iSp],0.25); cRYL[iSp]->cd(padC++);
		

		cRYL[iSp]->Write();
	}

	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0syst::ProcessRawYieldLossHist(TH2D* hist, TH1D* yieldhist, Int_t Sp, Double_t loose, Int_t opt) {

	Double_t* yieldLoose =
		((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSB((TH1D*)hist->ProjectionY(Form("iSp%i_iBin0",Sp),
			hist->GetXaxis()->FindBin(loose),hist->GetXaxis()->FindBin(loose)));
	//cout << *(yieldLoose) << " +- " << *(yieldLoose+1) << endl;
	Double_t cntL = *(yieldLoose);

	for (int iBin = 1; iBin < hist->GetNbinsX()+1; ++iBin)	{
		Double_t* yield =
			((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSB((TH1D*)hist->ProjectionY(Form("iSp%i_iBin0",Sp),
			iBin,iBin));

		Double_t cnt = *(yield);
		Double_t fractionY = ( cntL < cnt || cntL == 0) ? 0. :
			(1. - cnt/cntL);
		
		yieldhist->SetBinContent(iBin,fractionY);
		yieldhist->SetBinError(iBin,0);			 
	}

	// Remove junk from the file memory
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
		}
}

void MyAnalysisV0syst::DrawRawYieldLossHist(TH1D* da, TH1D* mc, Double_t ymax) {

	mHandler->MakeNiceHistogram(da,kBlack);
	mHandler->MakeNiceHistogram(mc,kRed);
	da->SetMarkerSize(0.9); mc->SetMarkerSize(0.9);
	da->GetYaxis()->SetTitle("Signal loss fraction w.r.t. loosest cut");

	da->GetYaxis()->SetRangeUser(0.,ymax);
	da->Draw("p");
	mc->Draw("p same");

}

void MyAnalysisV0syst::DrawVariation(Double_t cut, Int_t col, TVirtualPad* can) {

	Double_t x[2] = {cut, cut};
	Double_t y[2];
	y[0] = can->GetUymin(); y[1] = can->GetUymax();

	TGraph* gcut = new TGraph(2, x, y);
	gcut->SetLineWidth(2);
	gcut->SetLineColor(col);
	gcut->SetLineStyle(2);

	gcut->Draw("same");
}


