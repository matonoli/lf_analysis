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

		hV0PtNt[iSp][iType][iReg]->Write();
		hV0PtNtMin[iSp][iType][iReg]->Write();

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


	/*for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNtFitCorrUnf[iSp][iType][iReg]		= (TH2F*)hV0PtNtFitCorr[iSp][iType][iReg]->Clone(
			Form("hV0PtNtFitCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg])	);

		hV0PtNtMinFitCorrUnf[iSp][iType][iReg]		= (TH2F*)hV0PtNtMinFitCorr[iSp][iType][iReg]->Clone(
			Form("hV0PtNtMinFitCorrUnf_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg])	);

	}	}	}*/


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

	mUnf = new UnfoldNTclass();
	cout << "Unfolding class created " << mUnf << endl;

	DoUnfoldingNt();
	DoUnfolding1D();

	delete mUnf; mUnf = new UnfoldNTclass();
	DoUnfoldingNtMin();
	DoUnfolding1DMin();
	
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
	}
		
	return 0;	
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

void MyAnalysisV0unfold::DoUnfolding1D() {

	mHandler->root()->SetBatch(kTRUE);

	enum { D, RC, MC };
	const char* dOut = "results_unfolding";
	const bool eRM = false;
	const int NumberOfIters = 5;
	const char* Regions[3] = {"Trans1D","Toward","Away"};

	TFile* fUnfOut1D = (mHandler->GetFlagMC()) ? new TFile(Form("./%s/2D_newClass_mc.root",dOut),"RECREATE")
		: new TFile(Form("./%s/2D_newClass_data.root",dOut),"RECREATE");
	TList* lOut = new TList();
	lOut->SetOwner();

	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iReg = 0; iReg < 3; ++iReg)		{
	Int_t iType = mHandler->GetFlagMC() ? RC : D;

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[iReg]);
		obj->SetnIter(NumberOfIters);
		
		obj->SetPid(SPECIES[iSp]);
		obj->SetPidIdx(3+iSp);
		obj->SetMCAnalysis(mHandler->GetFlagMC());
		obj->SaveSolutionNT(kFALSE);

		printf(" - Unfolding species %s in region : %s\n",SPECIES[iSp],obj->GetRegion());
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
	const char* Regions[4] = {"Trans1D","Toward","Away","Trans1D"};
	

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