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
#include <TTree.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TF1.h>

#include "MyAnalysisV0extract.h"
#include "MyEvent.h"
#include "MyTrack.h"
#include "MyParticle.h"
#include "MyV0.h"
#include "MyHandler.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooChi2Var.h"

using namespace V0consts;
using namespace RooFit;

ClassImp(MyAnalysisV0extract)

MyAnalysisV0extract::MyAnalysisV0extract() {

}

Int_t MyAnalysisV0extract::Init() {

	mDirFile->cd();

	printf("Initialising analysis %s  \n", 
		this->GetName());

	cout << "ex mdir is " << mDirFile << endl;

	TH1::SetDefaultSumw2(1);
	CreateHistograms();
	BorrowHistograms();
	


	mList = (TList*)mHandler->directory()->GetList();

	return 0;
}

Int_t MyAnalysisV0extract::Make(Int_t iEv) {
	return 0;
}


Bool_t MyAnalysisV0extract::BorrowHistograms() {

	hNchTrans	= (TH1F*)mHandler->analysis(0)->dirFile()->Get("hNchTrans");

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
		
		if (iMu == 0 && iSph > 0) continue;	

		hV0IMvPt[iSp][iType][iMu][iSph] 
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));
		
		hV0IMvPt[iSp][iType][iMu][iSph] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPt[iSp][iType][iMu][iSph]);
	

	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
		
		hV0IMvPtNt[iSp][iType][iReg] 
			= (TH3F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtNt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		hV0IMvPtNtMin[iSp][iType][iReg] 
			= (TH3F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtNtMin_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
	
		hV0IMvPtNt[iSp][iType][iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH3(hV0IMvPtNt[iSp][iType][iReg]);
		hV0IMvPtNtMin[iSp][iType][iReg] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH3(hV0IMvPtNtMin[iSp][iType][iReg]);


	} } } 

	if (mHandler->GetFlagMC()) {
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		Int_t iType = 1, iMu = 0, iSph = 0;	

		hV0Pt[iSp][iType][iMu][iSph] 
			= (TH1D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0Pt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));

		if (hV0Pt[iSp][iType][iMu][iSph]->GetNbinsX() != NPTBINS)	{
			
			hV0Pt[iSp][iType][iMu][iSph] = (TH1D*)hV0Pt[iSp][iType][iMu][iSph]->Rebin(NPTBINS,hV0Pt[iSp][iType][iMu][iSph]->GetName(),XBINS);
		}
		

		hV0IMvPtPrimary[iSp]		= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtPrimary_%s",SPECIES[iSp]));
		hV0IMvPtPrimaryPDG[iSp]		= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtPrimaryPDG_%s",SPECIES[iSp]));
		hV0IMvPtSecondary[iSp]		= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSecondary_%s",SPECIES[iSp]));
		hV0IMvPtSecondaryPDG[iSp]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSecondaryPDG_%s",SPECIES[iSp]));
		hV0IMvPtSecondaryXi[iSp]	= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSecondaryXi_%s",SPECIES[iSp]));
		hV0IMvPtBackground[iSp]		= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtBackground_%s",SPECIES[iSp]));

		hV0IMvPtPrimary[iSp]		= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtPrimary[iSp]);
		hV0IMvPtPrimaryPDG[iSp]		= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtPrimaryPDG[iSp]);
		hV0IMvPtSecondary[iSp]		= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSecondary[iSp]);
		hV0IMvPtSecondaryPDG[iSp]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSecondaryPDG[iSp]);
		hV0IMvPtSecondaryXi[iSp]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSecondaryXi[iSp]);
		hV0IMvPtBackground[iSp]		= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtBackground[iSp]);		


	}	}

	return true;		

}

Bool_t MyAnalysisV0extract::CreateHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		hV0PtFit[iSp][iType][iMu][iSph]		= new TH1F(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		

	} } } }


	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		hV0PtNtFit[iSp][iType][iReg]		= new TH2F(Form("hV0PtNtFit_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 p_{T} (GeV/#it{c}); N_{ch} [trans.]; Entries",	NPTBINS,XBINS, 60, -0.5, 59.5);

		hV0PtNtMinFit[iSp][iType][iReg]		= new TH2F(Form("hV0PtNtMinFit_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]),
			";V0 p_{T} (GeV/#it{c}); N_{ch} [trans.,min]; Entries",	NPTBINS,XBINS, 60, -0.5, 59.5);

	}	}	}

	cout << "N OF BINS IS -------------------------" << endl;
	cout << NPTBINS << endl;

	if (mHandler->GetFlagMC()) {
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

		hV0PtFitPrimary[iSp]		= new TH1F(Form("hV0PtFitPrimary_%s",SPECIES[iSp]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hV0PtFitPrimaryPDG[iSp]		= new TH1F(Form("hV0PtFitPrimaryPDG_%s",SPECIES[iSp]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hV0PtFitSecondary[iSp]		= new TH1F(Form("hV0PtFitSecondary_%s",SPECIES[iSp]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hV0PtFitSecondaryPDG[iSp]	= new TH1F(Form("hV0PtFitSecondaryPDG_%s",SPECIES[iSp]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hV0PtFitSecondaryXi[iSp]	= new TH1F(Form("hV0PtFitSecondaryXi_%s",SPECIES[iSp]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		hV0PtFitBackground[iSp]		= new TH1F(Form("hV0PtFitBackground_%s",SPECIES[iSp]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);

	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

		hSidebandMean[iSp]	= new TH1F(Form("hSidebandMean_%s", SPECIES[iSp]),"; p_{T} bin; SB #mu;",NPTBINS,XBINS);
		hSidebandSigma[iSp]	= new TH1F(Form("hSidebandSigma_%s", SPECIES[iSp]),"; p_{T} bin; SB #sigma;",NPTBINS,XBINS);
		hSidebandSF[iSp]	= new TH1F(Form("hSidebandSF_%s", SPECIES[iSp]),"; p_{T} bin; SB scale factor;",NPTBINS,XBINS);
		/*
		hFitParam0[iSp]	= new TH1F(Form("hFitParam0_%s", SPECIES[iSp]),"; p_{T} bin; gaus1 #mu;",NPTBINS,XBINS);
		hFitParam1[iSp]	= new TH1F(Form("hFitParam1_%s", SPECIES[iSp]),"; p_{T} bin; gaus1 #sigma;",NPTBINS,XBINS);
		hFitParam2[iSp]	= new TH1F(Form("hFitParam2_%s", SPECIES[iSp]),"; p_{T} bin; gaus2 #sigma;",NPTBINS,XBINS);
		hFitParam3[iSp]	= new TH1F(Form("hFitParam3_%s", SPECIES[iSp]),"; p_{T} bin; gaus1 fraction;",NPTBINS,XBINS);
		*/
	}

	const int NBINS_IM = 125;
	double XBINS_IM[NBINS_IM+1];// = {};
	double binStepIM = (0.2 + 0.2)/(double)NBINS_IM;
	for (int i=0; i <= NBINS_IM; i++) XBINS_IM[i] = -0.2 + (double)i*binStepIM;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{

		hV0IMvPtRt[iSp][iType][iReg][iRtBin]		= new TH2F(Form("hV0IMvPtRt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),
			";V0 p_{T} (GeV/#it{c}); V0 m (GeV/#it{c}^{2}); Entries",		NPTBINS, XBINS, NBINS_IM, XBINS_IM);

		hV0PtRtFit[iSp][iType][iReg][iRtBin]		= new TH1F(Form("hV0PtRtFit_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);

	} } } }

	return true;

}


Int_t MyAnalysisV0extract::Finish() {

	printf("Finishing analysis %s \n",this->GetName());

	mDirFile->cd();

	//StudyIMShapeRC();

	//GetTemplates();

	//DefineSidebands();
	TakeoverSidebands();
	ProducePtSpectraFromHists();
	//if (MAKE_EXCLUSIVE) MakeExclusiveS0Bins();
	
	//ConvertNttoRt();
	ProducePtSpectraFromHistsRt();

	if (mHandler->GetFlagMC()) DoClosureTest(0);
	//DrawPad(0,0);

	//DrawConstraints();

	return 0;	
}

void MyAnalysisV0extract::ConvertNttoRt() {

	Double_t ntedges[NRTBINS0+1];
	for (int i=0; i<NRTBINS0+1; i++) ntedges[i] = RTBINS0[i]*RT_DEN;
	Int_t binedges[NRTBINS0+1];
	for (int i=0; i<NRTBINS0+1; i++) binedges[i] = hV0IMvPtNt[1][0][0]->GetZaxis()->FindBin(ntedges[i]);
	for (int i=0; i<NRTBINS0+1; i++) cout << ntedges[i] << " -- " << binedges[i] << endl;

	TAxis *xaxis = hV0IMvPtNt[1][0][0]->GetXaxis(); 
	TAxis *yaxis = hV0IMvPtNt[1][0][0]->GetYaxis(); 
	TAxis *zaxis = hV0IMvPtNt[1][0][0]->GetZaxis();


	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{

		for (int i=1; i<=xaxis->GetNbins();i++)	{
			for (int j=1; j<=yaxis->GetNbins();j++)	{
				int bincount = 0; 
				for (int k=binedges[iRtBin]; k<=binedges[iRtBin+1];k++)	{
					bincount += hV0IMvPtNt[iSp][iType][iReg]->GetBinContent(i,j,k);
				}
				hV0IMvPtRt[iSp][iType][iReg][iRtBin]->SetBinContent(i,j,bincount);
				hV0IMvPtRt[iSp][iType][iReg][iRtBin]->SetBinError(i,j,TMath::Sqrt(bincount));
			}
		}
		//cout << hV0IMvPtNt[iSp][iType][iReg] << endl;
		//hV0IMvPtNt[iSp][iType][iReg]->GetZaxis()->SetRange(0,99);
		//hV0IMvPtNt[iSp][iType][iReg]->GetZaxis()->SetRangeUser(ntedges[iRtBin],ntedges[iRtBin+1]);
		
		//hV0IMvPtNt[iSp][iType][iReg]->GetZaxis()->SetRange(binedges[iRtBin],binedges[iRtBin+1]);
		//hV0IMvPtRt[iSp][iType][iReg][iRtBin] = (TH2F*)hV0IMvPtNt[iSp][iType][iReg]->Project3D("yx");//,hV0IMvPtNt[iSp][iType][iReg]->GetZaxis()->FindBin(ntedges[iRtBin]),hV0IMvPtNt[iSp][iType][iReg]->GetZaxis()->FindBin(ntedges[iRtBin+1]));
		
		//hV0IMvPtRt[iSp][iType][iReg][iRtBin]->SetName(Form("hV0IMvPtRt_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iRtBin));

	}	}	}	}
} 

void MyAnalysisV0extract::DrawConstraints() { 

	TCanvas* cParams = new TCanvas("cParams","",1000,900);
	mHandler->MakeNiceHistogram(hFitParam1[1],1);
	mHandler->MakeNiceHistogram(hFitParam2[1],2);
	hFitParam1[1]->GetYaxis()->SetRangeUser(0.0,0.03);
	hFitParam1[1]->Draw();
	hFitParam2[1]->Draw("same");
}

void MyAnalysisV0extract::DefineSidebands() {

	mHandler->root()->SetBatch(kTRUE);
	nPads = TMath::Nint(TMath::Sqrt(NPTBINS));

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		int iType = 0; int iMu = 0; int iSph = 0;

		Float_t fitMin = -0.05, fitMax = 0.05;
		cFitsSB[iSp] = new TCanvas(Form("cFitsSB_iSp%i",iSp),"",2800,2000);
		cFitsSB[iSp]->Divide(nPads+1,nPads-1,0.00005,0.00005);


		for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1)	{

			TH1F* hist = (TH1F*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, iBin),
				iBin,iBin);
			

			TH1F* histrc;
			if (mHandler->GetFlagMC()) {
				//histrc = (TH1F*)hV0IMvPt[iSp][1][iMu][iSph]->ProjectionY(
				//Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, 1, iMu, iSph, iBin),
				//iBin,iBin);
				//histrc->Rebin(10);
			}

			Int_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)));// == 0);
			if (empty<10) continue;
			Double_t hmax = hist->GetMaximum();

			RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
			hist->Rebin(10);
			RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

			RooRealVar pGaus1A("pGaus1A","Mean 1",-0.0019,0.0019);
			RooRealVar pGaus1B("pGaus1B","Sigma 1",0.001,0.00001,0.015);
			RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
			RooRealVar nGaus1("nGaus1","N_{Gaus1}",0.5,0,9e08);
		
			//RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
			RooRealVar pGaus2B("pGaus2B","Sigma 2",0.008,0.0002,0.025);
			RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
			RooRealVar nGaus2("nGaus2","N_{Gaus2}",0.5,0,9e08);
		
			RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-9,9);
			RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-9,9);
			RooRealVar pPolBgC("pPolBgC","Pol. par. C",0,-9,9);
			RooChebychev fPolBg = RooChebychev("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA));
			RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,9e08);

			RooAddPdf fGaus("fGaus","fGaus",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
			RooRealVar nGaus("nGaus","N_{Gaus}",0.5*hmax,0,9e09);

			RooAddPdf fTotal = RooAddPdf("fTotal","fTotal",RooArgList(fGaus,fPolBg),RooArgList(nGaus,nPolBg));

			RooFitResult* fR = 0; 
			printf("Trying to fit iSp %i bin %i \n", iSp, iBin);
			//if (iSp == 1 && empty>10) fR = fTotal.chi2FitTo(DT_set,Save(),PrintLevel(-1));
			//if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));


			// GENERATE RMS
			RooDataSet* histSigma = fGaus.generate(MassDT,100000);
			Double_t Sigma = histSigma->sigma(MassDT);
			Double_t Mean = pGaus1A.getVal();

			// Using a very simplified fixed SB for lambdas
			if (iSp > 1) {
			//if (iSp > 0) {
				Mean = 0.;
				Sigma = 0.0018;
			}


			// SAVE PARS
			hSidebandMean[iSp]->SetBinContent(iBin,Mean);
			hSidebandMean[iSp]->SetBinError(iBin,0.0001+pGaus1A.getError());
			hSidebandSigma[iSp]->SetBinContent(iBin,Sigma);
			Double_t SigmaError = TMath::Sqrt(pGaus1B.getError()*pGaus1B.getError()+pGaus2B.getError()*pGaus2B.getError());
			hSidebandSigma[iSp]->SetBinError(iBin,0.0002+SigmaError);

			// NOW FIT BACKGROUND
			Double_t SBscale;
			Double_t NSig = 6.;
			TF1 *fbg;/* = new TF1("fbg",gfpol3,Mean-3.*NSig*Sigma,Mean+3.*NSig*Sigma,3);
			fbg->SetParameters(1.,0.,0.);
			//fit only background excluding the signal area
			gFReject = kTRUE;
			//because ROOT is poorly designed we have to use globals now, careful (!)
			gFLeft = Mean-1.*NSig*Sigma; gFRight = Mean+1.*NSig*Sigma;
			TFitResultPtr r = hist->Fit(fbg,"0SQ");
			gFReject = kFALSE;
			//store 2 separate functions for visualization

			double a = Mean-NSig*Sigma;
			double b = Mean+NSig*Sigma;
			double c = Mean-2.*NSig*Sigma;
			double d = Mean+2.*NSig*Sigma;
			SBscale 	= fbg->Integral(a,b);
			Double_t SBscaleErr = fbg->IntegralError(a,b,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());
			SBscale /= (fbg->Integral(c,a)+fbg->Integral(b,d)) > 0 ? (fbg->Integral(c,a)+fbg->Integral(b,d)) : 1;
			
			SBscale = 1.;
			//cout << "sb scale is " << SBscale << endl;
			//cout << "sb scale is " << SBscale << endl;
			//if (hSidebandSF[spNumber]->GetBinContent(binNumber)==0) {
			hSidebandSF[iSp]->SetBinContent(iBin,SBscale);
			hSidebandSF[iSp]->SetBinError(iBin,0.5*SBscaleErr);
			//}
   			*/


			cFitsSB[iSp]->cd(iBin);
			RooPlot* plot1 = MassDT.frame(Title(" "));
			DT_set.plotOn(plot1,MarkerSize(0.4));
			if (empty!=0) {
				//fTotal.plotOn(plot1,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
				fTotal.plotOn(plot1,Components(fGaus),LineStyle(2),LineWidth(2),LineColor(kRed));
				fTotal.plotOn(plot1,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
			}
			plot1->SetMinimum(1e-05);
			plot1->SetMaximum(1.40*plot1->GetMaximum());
			plot1->GetXaxis()->SetTitleSize(0.05);
			plot1->GetYaxis()->SetTitleSize(0.05);
			plot1->Draw();
			//if (mHandler->GetFlagMC()) hist->Add(histrc,-1.);
			hist->GetXaxis()->SetRangeUser(-0.099,0.099);
			hist->SetStats(0);
			hist->Draw();


			//fbg->Draw("same");

			TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
			mHandler->MakeNiceLegend(leg1, 0.08, 1.);
			leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",XBINS[iBin-1],XBINS[iBin])," ");
			leg1->AddEntry((TObject*)0,Form("%s",SPECNAMES[iSp])," ");
			leg1->Draw();

		}

		//cFitsSB[iSp]->SaveAs(Form("tmp/%s.png",cFitsSB[iSp]->GetName()));
		cFitsSB[iSp]->Write();

		// INTERPOLATE PARAMETERS
		TString periodName;// = mHandler->filehist()->GetName();
		TF1 *fsigma;
		if (periodName.Contains("18n"))	{
			fsigma = new TF1(Form("fsigma_%i",iSp),"pol1",1e-3+XBINS[0],XBINS[NPTBINS]);
		} else {
			if (iSp==1) fsigma = new TF1(Form("fsigma_%i",iSp),"[0]+[1]*x+[2]/x",1e-3+XBINS[0],XBINS[NPTBINS]);
			else fsigma = new TF1(Form("fsigma_%i",iSp),"pol0",1e-3+XBINS[0],XBINS[NPTBINS]);
		} 
		//TF1 *fmean	= new TF1(Form("fsigma_%i",iSp),"[0]+[1]*x+[2]/x",XBINS[0],XBINS[NPTBINS]);
		TF1 *fsf;
		if (iSp == 1) {
			fsf = new TF1(Form("fsf_%i",iSp),"[0]+[1]/x+[2]*x^[3]",1e-3+XBINS[0],XBINS[NPTBINS]);
			fsf->SetParameters(-6.,0.003,7.,0.01); }
		else fsf = new TF1(Form("fsf_%i",iSp),"pol0",1e-3+XBINS[0],XBINS[NPTBINS]);

		hSidebandSigma[iSp]->Fit(fsigma,"W");

		//hSidebandSF[iSp]->Fit(fsf,"W");

	}
	
	{
		TCanvas* cSBmu = new TCanvas("cSBmu","",1800,900);
		cSBmu->Divide(2,1);
		TCanvas* cSBsig = new TCanvas("cSBsig","",900,900);
		for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
			cSBmu->cd(1);
			mHandler->MakeNiceHistogram(hSidebandMean[iSp],COLOURS[iSp]);
			mHandler->MakeNiceHistogram(hSidebandSigma[iSp],COLOURS[iSp]);
			hSidebandSigma[iSp]->GetFunction(Form("fsigma_%i",iSp))->SetLineColor(COLOURS[iSp]);
			cSBmu->cd(1);
			hSidebandMean[iSp]->Draw("same");
			cSBmu->cd(2);
			hSidebandSigma[iSp]->Draw("same");
		}
		TLegend* leg1 = new TLegend(0.19,0.13,0.9,0.25);//cFits[canCounter/NPTBINS]->BuildLegend();
		mHandler->MakeNiceLegend(leg1, 0.065, 3);
		for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
			leg1->AddEntry(hSidebandMean[iSp],Form("%s",SPECNAMES[iSp]),"pl");
		}
		leg1->Draw();
		//cSBmu->cd();
		cSBmu->Write();
		//cSBsig->cd();
		//leg1->Draw();
		cSBsig->Write();
	}
	mHandler->root()->SetBatch(kFALSE);

}

void MyAnalysisV0extract::TakeoverSidebands() {

	cout << "mpar " << mParMuK0s << endl; 

	mParSigK0s	= new TF1("funcSigK0s","[0]+[1]*x+[2]/x",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParSigK0s->SetParameters(cuts::K0S_PARSIG[0],cuts::K0S_PARSIG[1],cuts::K0S_PARSIG[2]);
	mParMuK0s	= new TF1("funcMuK0s","(x<=1.6)*([0]+[1]*x+[2]*x*x)+(x>1.6)*[3]",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParMuK0s->SetParameters(cuts::K0S_PARMU[0],cuts::K0S_PARMU[1],cuts::K0S_PARMU[2],cuts::K0S_PARMU[3]);
	
	mParSigL	= new TF1("funcSigL","[0]+[1]*x+[2]/x",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParSigL->SetParameters(cuts::L_PARSIG[0],cuts::L_PARSIG[1],cuts::L_PARSIG[2]);
	mParMuL	= new TF1("funcMuL","(x<=1.9)*([0]+[1]*x+[2]*x*x)+(x>1.9)*([3]+[4]*x)",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParMuL->SetParameters(cuts::L_PARMU[0],cuts::L_PARMU[1],cuts::L_PARMU[2],cuts::L_PARMU[3],cuts::L_PARMU[4]);
	
	mParSigK0s->Write();
	mParSigL->Write();
	mParMuK0s->Write();
	mParMuL->Write();
	
	//mParMuK0s	= (TF1*)mHandler->analysis(0)->dirFile()->Get("funcMuK0s");
	//mParSigK0s	= (TF1*)mHandler->analysis(0)->dirFile()->Get("funcSigK0s");
	//mParMuL		= (TF1*)mHandler->analysis(0)->dirFile()->Get("funcMuL");
	//mParSigL	= (TF1*)mHandler->analysis(0)->dirFile()->Get("funcSigL");

	cout << "mpar " << mParMuK0s << endl; 
}



void MyAnalysisV0extract::SetMCInputFile(const Char_t *name) {

	TString fileName = TString(name);
	if (fileName.Data() != "") {
		mFileMC = new TFile(fileName,"READ");
		printf("MC File %s loaded in. \n", fileName.Data()); }
	else {
		printf("No MC file loaded.");
	}
}


void MyAnalysisV0extract::ProducePtSpectraFromHists() {

	mHandler->root()->SetBatch(kTRUE);

	nBins 	= NPTBINS;
	xBins	= XBINS;
	Int_t binSize = 0;//-1 + TMath::Nint((Double_t)NPTBINS/nBins);	// this is buggy actually
	nPads = TMath::FloorNint(TMath::Sqrt(nBins));

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	iCan = 0;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{

		//if (iMu > 4 && (iSph < 3 && iSph)) continue;
		//if (iMu < 5 && iSph > 2) continue;
		if (iMu==0 && iSph) continue;

		printf("Extracting yield for pt spectrum iSp%i_iType%i_iMu%i_iSph%i \n",iSp,iType,iMu,iSph);

		cFits[iCan] = new TCanvas(Form("cFits_iSp%i_iType%i_iMu%i_iSph%i",iSp,iType,iMu,iSph),"",4200,3000);
		cFits[iCan]->Divide(nPads+2,nPads,0.00005,0.00005);
			
		Double_t* yield = 0;
		Int_t binCounter = 1;

		for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1+binSize)	{
		//for (int iBin = 10; iBin < 11; ++iBin)	{
			
			/*yield = ExtractYieldFit((TH1F*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, binCounter),
				iBin,iBin+binSize),		iType,	(iType==0&&iMu==3&&iSph==0) );		// should be perhaps changed to FindBin
			*/

			// FOR CLOSURE STUDY
			if (mHandler->GetFlagMC()) {
				
				yield = ExtractYieldSB((TH1F*)hV0IMvPtPrimary[iSp]->ProjectionY(
					Form("iSp%i_iBin%i", iSp, binCounter), iBin,iBin+binSize));
				hV0PtFitPrimary[iSp]->SetBinContent(binCounter,*(yield+0));
				hV0PtFitPrimary[iSp]->SetBinError(binCounter,*(yield+1));
				
				yield = ExtractYieldSB((TH1F*)hV0IMvPtPrimaryPDG[iSp]->ProjectionY(
					Form("iSp%i_iBin%i", iSp, binCounter), iBin,iBin+binSize));
				hV0PtFitPrimaryPDG[iSp]->SetBinContent(binCounter,*(yield+0));
				hV0PtFitPrimaryPDG[iSp]->SetBinError(binCounter,*(yield+1));
				
				yield = ExtractYieldSB((TH1F*)hV0IMvPtSecondary[iSp]->ProjectionY(
					Form("iSp%i_iBin%i", iSp, binCounter), iBin,iBin+binSize));
				hV0PtFitSecondary[iSp]->SetBinContent(binCounter,*(yield+0));
				hV0PtFitSecondary[iSp]->SetBinError(binCounter,*(yield+1));
				
				yield = ExtractYieldSB((TH1F*)hV0IMvPtSecondaryPDG[iSp]->ProjectionY(
					Form("iSp%i_iBin%i", iSp, binCounter), iBin,iBin+binSize));
				hV0PtFitSecondaryPDG[iSp]->SetBinContent(binCounter,*(yield+0));
				hV0PtFitSecondaryPDG[iSp]->SetBinError(binCounter,*(yield+1));
				
				yield = ExtractYieldSB((TH1F*)hV0IMvPtSecondaryXi[iSp]->ProjectionY(
					Form("iSp%i_iBin%i", iSp, binCounter), iBin,iBin+binSize));
				hV0PtFitSecondaryXi[iSp]->SetBinContent(binCounter,*(yield+0));
				hV0PtFitSecondaryXi[iSp]->SetBinError(binCounter,*(yield+1));
				
				yield = ExtractYieldSB((TH1F*)hV0IMvPtBackground[iSp]->ProjectionY(
					Form("iSp%i_iBin%i_Background", iSp, binCounter), iBin,iBin+binSize));
				hV0PtFitBackground[iSp]->SetBinContent(binCounter,*(yield+0));
				hV0PtFitBackground[iSp]->SetBinError(binCounter,*(yield+1));
				
			}

			// ACTUAL RAW YIELDS
			yield = ExtractYieldSB((TH1F*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, binCounter),
				iBin,iBin+binSize));

			hV0PtFit[iSp][iType][iMu][iSph]->SetBinContent(binCounter,*(yield+0));
			hV0PtFit[iSp][iType][iMu][iSph]->SetBinError(binCounter,*(yield+1));
			binCounter++;

		}

		//if (iType==0&&iMu==0&&iSph==0) cFits[iCan]->SaveAs(Form("tmp/%s.png",cFits[iCan]->GetName()));
		cFits[iCan]->Write();
		//cFits[iCan]->SaveAs(Form("plots/cFits_%i.png",iCan));
		iCan++;

		//return 0;

		if (mHandler->GetFlagMC()) {
			hV0PtFitPrimary[iSp]->Scale(1,"width");
			hV0PtFitPrimaryPDG[iSp]->Scale(1,"width");
			hV0PtFitSecondary[iSp]->Scale(1,"width");
			hV0PtFitSecondaryPDG[iSp]->Scale(1,"width");
			hV0PtFitSecondaryXi[iSp]->Scale(1,"width");
			hV0PtFitBackground[iSp]->Scale(1,"width");
		}

		hV0PtFit[iSp][iType][iMu][iSph]->Scale(1,"width");
		hV0PtFit[iSp][iType][iMu][iSph]->SetMarkerColor(1);
		hV0PtFit[iSp][iType][iMu][iSph]->SetMarkerStyle(20);
		hV0PtFit[iSp][iType][iMu][iSph]->SetLineColor(2);

	} } } }

	// Remove junk from the file memory
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
		}


	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0extract::ProducePtSpectraFromHistsRt() {

	mHandler->root()->SetBatch(kTRUE);

	nBins 	= NPTBINS;
	xBins	= XBINS;
	Int_t binSize = 0;//-1 + TMath::Nint((Double_t)NPTBINS/nBins);	// this is buggy actually
	nPads = TMath::FloorNint(TMath::Sqrt(nBins));
	canCounterRt = 0;

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	iCan = 0;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		for (int iNt = 1; iNt < hV0IMvPtNt[iSp][iType][iReg]->GetZaxis()->GetNbins()+1; iNt++)	{

			hV0IMvPtNt[iSp][iType][iReg]->GetZaxis()->SetRange(iNt,iNt);
			TH2F* hXY = (TH2F*)hV0IMvPtNt[iSp][iType][iReg]->Project3D("yx");			

			printf("Extracting yield for pt spectrum iSp%i_iType%i_iReg%i_iNtBin%i \n",iSp,iType,iReg,iNt);
			
			cFitsRt[iCan] = new TCanvas(Form("cFitsRt_iSp%i_iType%i_iReg%i_iNtBin%i",iSp,iType,iReg,iNt),"",4200,3000);
			cFitsRt[iCan]->Divide(nPads+2,nPads,0.00005,0.00005);
			
			Double_t* yield = 0;
			Int_t binCounter = 1;


			for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1)	{
			
				// ACTUAL RAW YIELDS
				yield = ExtractYieldSB((TH1F*)hXY->ProjectionY(
					Form("Rt_iSp%i_iType%i_iReg%i_iRtBin%i_iBin%i", iSp, iType, iReg, iNt, iBin),
					iBin,iBin));

				hV0PtNtFit[iSp][iType][iReg]->SetBinContent(iBin,iNt,*(yield+0));
				hV0PtNtFit[iSp][iType][iReg]->SetBinError(iBin,iNt,*(yield+1));
				binCounter++;

			}

			if (iCan%60==10) cFitsRt[iCan]->Write();
			//cFitsRt[iCan]->SaveAs(Form("plots/cFitsRt_%i.png",iCan));
			iCan++;

		} 

		// SCALE BY BIN WIDTH
		for (int iX = 1; iX < hV0PtNtFit[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNtFit[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float binwidth = hV0PtNtFit[iSp][iType][iReg]->GetXaxis()->GetBinWidth(iX);
				float binc = hV0PtNtFit[iSp][iType][iReg]->GetBinContent(iX,iY);
				float bine = hV0PtNtFit[iSp][iType][iReg]->GetBinError(iX,iY);
				if (binwidth>0) hV0PtNtFit[iSp][iType][iReg]->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) hV0PtNtFit[iSp][iType][iReg]->SetBinError(iX,iY,bine/binwidth);
			}
		}

	} } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{

		for (int iNt = 1; iNt < hV0IMvPtNtMin[iSp][iType][iReg]->GetZaxis()->GetNbins()+1; iNt++)	{

			hV0IMvPtNtMin[iSp][iType][iReg]->GetZaxis()->SetRange(iNt,iNt);
			TH2F* hXY = (TH2F*)hV0IMvPtNtMin[iSp][iType][iReg]->Project3D("yx");			

			printf("Extracting yield for pt spectrum RtMin iSp%i_iType%i_iReg%i_iNtBin%i \n",iSp,iType,iReg,iNt);
			
			Double_t* yield = 0;
			Int_t binCounter = 1;


			for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1)	{
			
				// ACTUAL RAW YIELDS
				yield = ExtractYieldSB((TH1F*)hXY->ProjectionY(
					Form("RtMin_iSp%i_iType%i_iReg%i_iRtBin%i_iBin%i", iSp, iType, iReg, iNt, iBin),
					iBin,iBin),false);

				hV0PtNtMinFit[iSp][iType][iReg]->SetBinContent(iBin,iNt,*(yield+0));
				hV0PtNtMinFit[iSp][iType][iReg]->SetBinError(iBin,iNt,*(yield+1));
				binCounter++;

			}

		} 

		// SCALE BY BIN WIDTH
		for (int iX = 1; iX < hV0PtNtMinFit[iSp][iType][iReg]->GetNbinsX()+1; iX++) {
			for (int iY = 1; iY < hV0PtNtMinFit[iSp][iType][iReg]->GetNbinsY()+1; iY++) {	
				float binwidth = hV0PtNtMinFit[iSp][iType][iReg]->GetXaxis()->GetBinWidth(iX);
				float binc = hV0PtNtMinFit[iSp][iType][iReg]->GetBinContent(iX,iY);
				float bine = hV0PtNtMinFit[iSp][iType][iReg]->GetBinError(iX,iY);
				if (binwidth>0) hV0PtNtMinFit[iSp][iType][iReg]->SetBinContent(iX,iY,binc/binwidth);
				if (binwidth>0) hV0PtNtMinFit[iSp][iType][iReg]->SetBinError(iX,iY,bine/binwidth);
			}
		}

	} } }

	// Remove junk from the file memory
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
		}


	mHandler->root()->SetBatch(kFALSE);
}


void MyAnalysisV0extract::MakeExclusiveS0Bins() { 

	// creates exclusive spherocity bins, i.e 0-1, 1-5, instead of 0-1, 0-5

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	enum { sphMB, Jetty20, Iso20, Jetty10, Iso10, Jetty5, Iso5, Jetty1, Iso1};


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 1; iMu < NMULTI; ++iMu)		{

		hV0PtFit[iSp][iType][iMu][sphMB]->Add(hV0PtFit[iSp][iType][iMu][Jetty20],-1.);
		hV0PtFit[iSp][iType][iMu][sphMB]->Add(hV0PtFit[iSp][iType][iMu][Iso20],-1.);
		hV0PtFit[iSp][iType][iMu][Jetty20]->Add(hV0PtFit[iSp][iType][iMu][Jetty10],-1.);
		hV0PtFit[iSp][iType][iMu][Jetty10]->Add(hV0PtFit[iSp][iType][iMu][Jetty5],-1.); 
		hV0PtFit[iSp][iType][iMu][Jetty5]->Add(hV0PtFit[iSp][iType][iMu][Jetty1],-1.); 
		hV0PtFit[iSp][iType][iMu][Iso20]->Add(hV0PtFit[iSp][iType][iMu][Iso10],-1.);
		hV0PtFit[iSp][iType][iMu][Iso10]->Add(hV0PtFit[iSp][iType][iMu][Iso5],-1.); 
		hV0PtFit[iSp][iType][iMu][Iso5]->Add(hV0PtFit[iSp][iType][iMu][Iso1],-1.); 

	}	}	}


}

void MyAnalysisV0extract::ProducePtSpectraFromTrees() {

	/*mHandler->root()->SetBatch(kTRUE);

	nPads = TMath::Nint(TMath::Sqrt(NPTBINS));
	Double_t rt_den = hNchTrans->GetMean();


	nBins 	= NPTBINS;
	xBins	= XBINS;

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	iCan = 0;
	canCounter = 0;
	
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{

		cFitsPtTree[iCan] = new TCanvas(Form("cFitsPtTree_iSp%i_iType%i_iReg%i_iRtBin%i",iSp,iType,iReg,iRtBin),"",2800,2000);
		cFitsPtTree[iCan]->Divide(nPads+2,nPads,0.00005,0.00005);
			
		Double_t* yield = 0;
		Float_t leftrtb =	RTBINS0[iRtBin]*rt_den;
		Float_t rightrtb =	RTBINS0[iRtBin+1]*rt_den;
		//		std::cout << "left " << leftrtb << " right " << rightrtb << std::endl;
		printf("Extracting yield for pt spectrum from trees iSp%i_iType%i_iReg%i_iRtBin%i \n",iSp,iType,iReg,iRtBin);
		TString treeName = tV0massRt[iSp][iType][iReg]->GetName();
		Int_t binCounter = 1;
		for (int iBin = 0; iBin < NPTBINS; ++iBin)	{

			tV0massRt[iSp][iType][iReg]->SetName(Form("treept_iSp%i_iBin%i",iSp,iBin+1));

			//yield = ExtractYieldFitPtTree((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
			//	Form("lPt>%f && lPt<%f && lNchTrans>=%f && lNchTrans<%f",XBINS2[iBin],XBINS2[iBin+1],
			//		leftrtb-1E-05,rightrtb)),     iType);
			TH1F* tmp = new TH1F("tmp",";V0 m (GeV/#it{c}^{2}); Entries", 1000, -0.1, 0.1);
			tV0massRt[iSp][iType][iReg]->Draw("MassDT>>tmp",Form("lPt>%f && lPt<%f && lNchTrans>=%f && lNchTrans<%f",XBINS2[iBin],XBINS2[iBin+1],
					leftrtb-1E-05,rightrtb),"goff");
			tmp->SetName(Form("treept_iSp%i_iBin%i",iSp,iBin+1));

			yield = ExtractYieldSB(tmp);

			tV0massRt[iSp][iType][iReg]->SetName(treeName);

			//std::cout << "yield " << *(yield+0) << std::endl;
			hV0PtRtFit[iSp][iType][iReg][iRtBin]->SetBinContent(iBin+1,*(yield+0));	//+1 for underflow bin
			hV0PtRtFit[iSp][iType][iReg][iRtBin]->SetBinError(iBin+1,*(yield+1));
			//hV0PtFit[iSp][iType][iMu][iSph]->SetBinError(binCounter,*(yield+1));
			binCounter++;
			delete tmp;
		}

		cFitsPtTree[iCan]->Write();
		iCan++;

		//return 0;
		hV0PtRtFit[iSp][iType][iReg][iRtBin]->Scale(1,"width");
		hV0PtRtFit[iSp][iType][iReg][iRtBin]->SetMarkerColor(1);
		hV0PtRtFit[iSp][iType][iReg][iRtBin]->SetMarkerStyle(20);
		hV0PtRtFit[iSp][iType][iReg][iRtBin]->SetLineColor(2);

	} } } }

	// Remove junk from the file memory
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("tV") ) {
			mDirFile->Remove(obj);	}
		}


	mHandler->root()->SetBatch(kFALSE);
	*/
}



void MyAnalysisV0extract::ProduceRtSpectraFromTrees() {

	/*mHandler->root()->SetBatch(kTRUE);

	Double_t rt_den = hNchTrans->GetMean();

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	iCan = 0;

	nchmax = hNchTrans->FindLastBinAbove();
	std::cout << "nchmax " << nchmax << std::endl;
	increm = 2;
	nPads = TMath::Nint(TMath::Sqrt((nchmax/increm)));

	canCounterRt = 0;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)			{
		
		// canvases w rt binning

		cFitsRt[iCan] = new TCanvas(Form("cFitsRt_iSp%i_iType%i_iReg%i_iPtBin%i",iSp,iType,iReg,iPtBin),"",2800,2000);
		cFitsRt[iCan]->Divide(nPads+1,nPads,0.00005,0.00005);


		//std::cout << " tV0massRt[iSp][0][0] " << tV0massRt[iSp][0][0] << std::endl;
		TString treeName = tV0massRt[iSp][iType][iReg]->GetName();

		Double_t* yield = 0;

		//for entire rt range
		tV0massRt[iSp][iType][iReg]->SetName(Form("tree_iSp%i_iType%i_iReg%i_iPtBin%i_iBin%i",iSp,iType,iReg,iPtBin,nchmax));//NRTBINS));
			
		yield = ExtractYieldFitRt((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
			Form("lPt>%f && lPt<%f && lNchTrans>%f && lNchTrans<%f",RT_PTRANGE[iPtBin][0],RT_PTRANGE[iPtBin][1],
					0-1E-4, (Double_t)nchmax)),			iType);
		hRtV0Yields[iType][iReg][iPtBin]->SetBinContent(1+iSp,*(yield+0));
		hRtV0Yields[iType][iReg][iPtBin]->SetBinError(1+iSp,*(yield+1));
		//

		tV0massRt[iSp][iType][iReg]->SetName(treeName);


		for (int iBin = 0; iBin < nchmax; iBin+=increm)	{
			
			tV0massRt[iSp][iType][iReg]->SetName(Form("tree_iSp%i_iType%i_iReg%i_iPtBin%i_iBin%i",iSp,iType,iReg,iPtBin,iBin));
			//std::cout << "ibin " << iBin << " jo " << (iBin < nchmax) << std::endl;
			yield = ExtractYieldFitRt((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
				Form("lPt>%f && lPt<%f && lNchTrans>%f && lNchTrans<%f",RT_PTRANGE[iPtBin][0],RT_PTRANGE[iPtBin][1],
					iBin-1E-4, iBin-1E-4+increm)),			iType);
					//RTBINS[iBin]*rt_den,RTBINS[iBin+1]*rt_den)),		0);
			
			hV0RtFit[iSp][iType][iReg][iPtBin]->SetBinContent(iBin+1,*(yield+0));
			hV0RtFit[iSp][iType][iReg][iPtBin]->SetBinError(iBin+1,*(yield+1));
		}
		hV0RtFit[iSp][iType][iReg][iPtBin]->Rebin(increm);

		//for entire rt range
		//tV0massRt[iSp][iType][iReg]->SetName(Form("tree_iSp%i_iType%i_iReg%i_iPtBin%i_iBin%i",iSp,iType,iReg,iPtBin,nchmax));//NRTBINS));
			
		//yield = ExtractYieldFitRt((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
		//	Form("lPt>%f && lPt<%f && lNchTrans>%f && lNchTrans<%f",RT_PTRANGE[iPtBin][0],RT_PTRANGE[iPtBin][1],
		//			0-1E-4, (Double_t)nchmax),			0) );
		//hRtV0Yields[iType][iReg][iPtBin]->SetBinContent(1+iSp,*(yield+0));
		//hRtV0Yields[iType][iReg][iPtBin]->SetBinError(1+iSp,*(yield+1));
		//

		//tV0massRt[iSp][iType][iReg]->SetName(treeName);

		

		cFitsRt[iCan]->Write();
		iCan++;

		hV0RtFit[iSp][iType][iReg][iPtBin]->SetMarkerColor(1);
		hV0RtFit[iSp][iType][iReg][iPtBin]->SetMarkerStyle(20);
		hV0RtFit[iSp][iType][iReg][iPtBin]->SetLineColor(2);

	}	}	}	}

	// Remove junk from the file memory
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("tV") ) {
			mDirFile->Remove(obj);	}
		}

	mHandler->root()->SetBatch(kFALSE);*/
}

Double_t* MyAnalysisV0extract::ExtractYieldSB(TH1F* hist, Bool_t willDraw) {



	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.1, fitMax = 0.1;

	for (int ib = 1; ib < hist->GetNbinsX()+1; ib++) hist->SetBinError(ib,TMath::Sqrt(hist->GetBinContent(ib)));

	TString histName(hist->GetName());
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();

	Bool_t isRt = histName.Contains("Rt");

	/*TString muName(histName(histName.Index("iMu")+3,1));
	Int_t muNumber = muName.Atoi();*/

	Int_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)));// == 0);

	Double_t Mean, Sigma, SF;
	
	if (mParMuK0s) {
		Mean	= (spNumber==1) ? mParMuK0s->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber])) : mParMuL->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber]));
		Sigma	= (spNumber==1) ? mParSigK0s->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber])) : mParSigL->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber]));
		SF = 1;
	} else {

		Mean = hSidebandMean[spNumber]->GetBinContent(hSidebandMean[spNumber]->FindBin(xBins[binNumber-1]));
		Sigma = hSidebandSigma[spNumber]->GetFunction(Form("fsigma_%i",spNumber))->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber]));
		SF = 1;//hSidebandSF[spNumber]->GetFunction(Form("fsf_%i",spNumber))->Eval(0.5*(0.6*xBins[binNumber-1]+0.4*xBins[binNumber]));

	}
	//if (spNumber==1) printf("bin %i mean %f sigma %f in bin %f and %f \n", binNumber, Mean, Sigma, xBins[binNumber-1], xBins[binNumber]);

	Double_t NSig = 6; //6;
	//hist->Rebin(8);  // IM bins reduced to 125 so no need to rebin anymore
	
	TF1 *fbg = new TF1("fbg",gfpol3,Mean-3.*NSig*Sigma,Mean+3.*NSig*Sigma,3);
	fbg->SetParameters(1.,0.,0.);
	//fit only the linear background excluding the signal area
	gFReject = kTRUE;
	//because ROOT is poorly designed we have to use globals now, careful (!)
	gFLeft = Mean-1.*NSig*Sigma; gFRight = Mean+1.*NSig*Sigma;
	//hist->Fit(fbg,"0");
	gFReject = kFALSE;
	//store 2 separate functions for visualization
   	/*TF1 *fleft = new TF1("fleft",fbg,Double_t(Mean-3.*NSig*Sigma),Double_t(Mean-1.*NSig*Sigma),3);
	fleft->SetParameters(fbg->GetParameters());
	hist->GetListOfFunctions()->Add(fleft);
	mHandler->root()->GetListOfFunctions()->Remove(fleft);
	TF1 *fright = new TF1("fright",fbg,Double_t(Mean+1.*NSig*Sigma),Double_t(Mean+3.*NSig*Sigma),3);
	fright->SetParameters(fbg->GetParameters());
	hist->GetListOfFunctions()->Add(fright);
	mHandler->root()->GetListOfFunctions()->Remove(fright);*/
	




	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",Mean-3*NSig*Sigma,Mean+3*NSig*Sigma);
	RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

	RooRealVar pPolBgA("pPolBgA","Pol. par. A",1.,-800.,800.);
	RooRealVar pPolBgB("pPolBgB","Pol. par. B",1.,-80,80);
	RooRealVar pPolBgC("pPolBgC","Pol. par. C",0,-80,80);
	RooPolynomial fPolBg = RooPolynomial("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA,pPolBgB,pPolBgC));
	//RooUniform fPolBg = RooUniform("fPolBg","fPolBg",MassDT);//,RooArgSet(pPolBgA));
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,9e08);
	MassDT.setRange("lowR",Mean-3*NSig*Sigma, Mean-NSig*Sigma);		//fitting B
	MassDT.setRange("highR",Mean+NSig*Sigma, Mean+3.*NSig*Sigma);
	MassDT.setRange("fullR",Mean-3*NSig*Sigma, Mean+3.*NSig*Sigma);
	MassDT.setRange("sigR",Mean-NSig*Sigma, Mean+NSig*Sigma);
	MassDT.setRange("lowBR",Mean-2*NSig*Sigma, Mean-NSig*Sigma);	// SB region
	MassDT.setRange("highBR",Mean+NSig*Sigma, Mean+2.*NSig*Sigma);

	RooAddPdf fTotal = RooAddPdf("fTotal","fTotal",RooArgList(fPolBg),RooArgList(nPolBg));

	RooFitResult* fR = 0; 
	//if (empty>3) fR = fTotal.chi2FitTo(DT_set,Range("lowR,highR"),Save(),PrintLevel(-1));
	//fTotal.getParameters(DT_set)->Print("V");
	//if (empty>3) fR = fTotal.fitTo(*blindedData,Save(),PrintLevel(-1));
	RooAbsReal* fInt = fTotal.createIntegral(MassDT,NormSet(MassDT),Range("sigR"));
	RooAbsReal* fIntTot = fTotal.createIntegral(MassDT,NormSet(MassDT),Range("lowBR,highBR"));
	RooAbsReal* fIntTot2 = fTotal.createIntegral(MassDT,NormSet(MassDT),Range("fullR"));
	RooAbsReal* fIntTot3 = fTotal.createIntegral(MassDT);


	Double_t N 	= hist->Integral(hist->FindBin(Mean-NSig*Sigma),hist->FindBin(Mean+NSig*Sigma));

	N = DT_set.sumEntries(0,"sigR");
	//cout << "N 2 " << N << endl;
	//Double_t Bg = hist->Integral(hist->FindBin(Mean-2*NSig*Sigma),hist->FindBin(Mean-NSig*Sigma));
	//		Bg += hist->Integral(hist->FindBin(Mean+NSig*Sigma),hist->FindBin(Mean+2*NSig*Sigma));

	Double_t Bg = DT_set.sumEntries(0,"lowBR");
	//cout << Bg << endl;
	Bg += DT_set.sumEntries(0,"highBR");
	//cout << "sb bg " << Bg << endl;

	double a = Mean-NSig*Sigma;
	double b = Mean+NSig*Sigma;
	double c = Mean-2.*NSig*Sigma;
	double d = Mean+2.*NSig*Sigma;
	Double_t SBscale = fbg->Integral(a,b);
	//cout << "sb scale is " << SBscale << endl;
	SBscale /= (fbg->Integral(c,a)+fbg->Integral(b,d)) > 0 ? (fbg->Integral(c,a)+fbg->Integral(b,d)) : 1;
	//cout << "sb scale is " << SBscale << endl;
	//SBscale = 1.;
	//if (hSidebandSF[spNumber]->GetBinContent(binNumber)==0) {
	//		hSidebandSF[spNumber]->SetBinContent(binNumber,SBscale);
//			hSidebandSF[spNumber]->SetBinError(binNumber,0.1*SBscale);
//	}
	//Bg *= SBscale;
	//Bg *= fInt->getVal();
	//cout << fInt->getVal() << endl;
	//Bg /= fIntTot->getVal();
	//cout << fIntTot->getVal() << endl;
	//cout << fIntTot2->getVal() << endl;
	//cout << fIntTot3->getVal() << endl;
	
	//cout << "bg val " << Bg << " compared to " << hist->Integral(hist->FindBin(Mean-2.*NSig*Sigma),hist->FindBin(Mean-NSig*Sigma)) + hist->Integral(hist->FindBin(Mean+NSig*Sigma),hist->FindBin(Mean+2.*NSig*Sigma)) << endl;
	//cout << "SF " << SF << endl;


	Bg *= SF;
	val[0] = N - Bg;
	val[1] = TMath::Sqrt(TMath::Abs(N + Bg));

	if (histName.Contains("Sys")) return val;


	/*std::cout << " ib " << binNumber << std::endl;
	std::cout << "h " << hist << " ib " << binNumber << " c " << canCounter/nBins << 
	" p " << canCounter%nBins << std::endl;*/

	//cout << NPTBINS << " " << nBins << endl;
	//cout << isRt << " " << cFitsRt[canCounterRt/nBins] << " " << canCounterRt/nBins << " " << canCounterRt%nBins << endl;
	//cout << "aa " << canCounterRt << " " << nBins << endl;


	if (!willDraw) return val;


	if (!isRt) cFits[canCounter/nBins]->cd(1+canCounter%nBins);
	else cFitsRt[canCounterRt/nBins]->cd(1+canCounterRt%nBins);



	/*RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (empty!=0) {
		fTotal.plotOn(plot1,Components(fPolBg),Range("fullR"),RooFit::NormRange("lowR,highR"), LineStyle(3),LineWidth(2),LineColor(kRed));
		fTotal.plotOn(plot1,Components(fPolBg),Range("lowR"),RooFit::NormRange("lowR,highR"), LineStyle(2),LineWidth(2),LineColor(kBlue));
		fTotal.plotOn(plot1,Components(fPolBg),Range("highR"),RooFit::NormRange("lowR,highR"), LineStyle(2),LineWidth(2),LineColor(kBlue));
		
	}
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();*/

	mHandler->MakeNiceHistogram(hist,kBlack);
	hist->SetMarkerSize(0.5); hist->GetXaxis()->SetLabelSize(0.055);
	//hist->Rebin(8);
	hist->GetXaxis()->SetRangeUser(fitMin,fitMax);
	if (spNumber>1) hist->GetXaxis()->SetRangeUser(-0.05,0.05);
	hist->GetYaxis()->SetRangeUser(1e-3, hist->GetMaximum()>0 ? 2.*hist->GetMaximum() : 2);
	hist->Draw();
	//fbg->Draw("same");

	if (!isRt) {
		cFits[canCounter/nBins]->Update();
		mHandler->DrawCut(Mean+2*NSig*Sigma,1,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins));
		mHandler->DrawCut(Mean+NSig*Sigma,2,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins));
		mHandler->DrawCut(Mean-NSig*Sigma,1,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins));
		mHandler->DrawCut(Mean-2*NSig*Sigma,2,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins)); }
	else {
		//cout << "UGH " << cFitsRt[canCounterRt/nBins] << " " << Mean << " " << Sigma << " " << cFitsRt[canCounterRt/nBins]->GetPad(1+canCounterRt%nBins) << "\n";
		cFitsRt[canCounterRt/nBins]->Update();
		mHandler->DrawCut(Mean+2*NSig*Sigma,1,cFitsRt[canCounterRt/nBins]->GetPad(1+canCounterRt%nBins));
		mHandler->DrawCut(Mean+NSig*Sigma,2,cFitsRt[canCounterRt/nBins]->GetPad(1+canCounterRt%nBins));
		mHandler->DrawCut(Mean-NSig*Sigma,1,cFitsRt[canCounterRt/nBins]->GetPad(1+canCounterRt%nBins));
		mHandler->DrawCut(Mean-2*NSig*Sigma,2,cFitsRt[canCounterRt/nBins]->GetPad(1+canCounterRt%nBins));
	}

	
	TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.08, 1.);
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",xBins[binNumber-1],xBins[binNumber])," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	if (!isRt && histName.Contains("_iType")) canCounter++;
	if (isRt && histName.Contains("_iType")) canCounterRt++;
	return val;
}

Double_t* MyAnalysisV0extract::ExtractYieldSBVarySigma(Double_t nsig, TH1F* hist) {


	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.1, fitMax = 0.1;

	for (int ib = 1; ib < hist->GetNbinsX()+1; ib++) hist->SetBinError(ib,TMath::Sqrt(hist->GetBinContent(ib)));

	TString histName(hist->GetName());
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();

	Bool_t isTree = histName.Contains("tree");

	/*TString muName(histName(histName.Index("iMu")+3,1));
	Int_t muNumber = muName.Atoi();*/

	Int_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)));// == 0);

	Double_t Mean, Sigma, SF;
	
	if (mParMuK0s) {
		Mean	= (spNumber==1) ? mParMuK0s->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber])) : mParMuL->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber]));
		Sigma	= (spNumber==1) ? mParSigK0s->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber])) : mParSigL->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber]));
		SF = 1;
	} else {

		Mean = hSidebandMean[spNumber]->GetBinContent(hSidebandMean[spNumber]->FindBin(xBins[binNumber-1]));
		Sigma = hSidebandSigma[spNumber]->GetFunction(Form("fsigma_%i",spNumber))->Eval(0.5*(xBins[binNumber-1]+xBins[binNumber]));
		SF = 1;//hSidebandSF[spNumber]->GetFunction(Form("fsf_%i",spNumber))->Eval(0.5*(0.6*xBins[binNumber-1]+0.4*xBins[binNumber]));

	}
	//if (spNumber==1) printf("bin %i mean %f sigma %f in bin %f and %f \n", binNumber, Mean, Sigma, xBins[binNumber-1], xBins[binNumber]);

	Double_t NSig = nsig;
	hist->Rebin(8);
	
	TF1 *fbg = new TF1("fbg",gfpol3,Mean-3.*NSig*Sigma,Mean+3.*NSig*Sigma,3);
	fbg->SetParameters(1.,0.,0.);
	//fit only the linear background excluding the signal area
	gFReject = kTRUE;
	//because ROOT is poorly designed we have to use globals now, careful (!)
	gFLeft = Mean-1.*NSig*Sigma; gFRight = Mean+1.*NSig*Sigma;
	//hist->Fit(fbg,"0");
	gFReject = kFALSE;
	//store 2 separate functions for visualization
   	/*TF1 *fleft = new TF1("fleft",fbg,Double_t(Mean-3.*NSig*Sigma),Double_t(Mean-1.*NSig*Sigma),3);
	fleft->SetParameters(fbg->GetParameters());
	hist->GetListOfFunctions()->Add(fleft);
	mHandler->root()->GetListOfFunctions()->Remove(fleft);
	TF1 *fright = new TF1("fright",fbg,Double_t(Mean+1.*NSig*Sigma),Double_t(Mean+3.*NSig*Sigma),3);
	fright->SetParameters(fbg->GetParameters());
	hist->GetListOfFunctions()->Add(fright);
	mHandler->root()->GetListOfFunctions()->Remove(fright);*/
	





	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",Mean-3*NSig*Sigma,Mean+3*NSig*Sigma);
	RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

	RooRealVar pPolBgA("pPolBgA","Pol. par. A",1.,-800.,800.);
	RooRealVar pPolBgB("pPolBgB","Pol. par. B",1.,-80,80);
	RooRealVar pPolBgC("pPolBgC","Pol. par. C",0,-80,80);
	RooPolynomial fPolBg = RooPolynomial("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA,pPolBgB,pPolBgC));
	//RooUniform fPolBg = RooUniform("fPolBg","fPolBg",MassDT);//,RooArgSet(pPolBgA));
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,9e08);
	MassDT.setRange("lowR",Mean-3*NSig*Sigma, Mean-NSig*Sigma);		//fitting B
	MassDT.setRange("highR",Mean+NSig*Sigma, Mean+3.*NSig*Sigma);
	MassDT.setRange("fullR",Mean-3*NSig*Sigma, Mean+3.*NSig*Sigma);
	MassDT.setRange("sigR",Mean-NSig*Sigma, Mean+NSig*Sigma);
	MassDT.setRange("lowBR",Mean-2*NSig*Sigma, Mean-NSig*Sigma);	// SB region
	MassDT.setRange("highBR",Mean+NSig*Sigma, Mean+2.*NSig*Sigma);

	RooAddPdf fTotal = RooAddPdf("fTotal","fTotal",RooArgList(fPolBg),RooArgList(nPolBg));

	RooFitResult* fR = 0; 
	//if (empty>3) fR = fTotal.chi2FitTo(DT_set,Range("lowR,highR"),Save(),PrintLevel(-1));
	//fTotal.getParameters(DT_set)->Print("V");
	//if (empty>3) fR = fTotal.fitTo(*blindedData,Save(),PrintLevel(-1));
	RooAbsReal* fInt = fTotal.createIntegral(MassDT,NormSet(MassDT),Range("sigR"));
	RooAbsReal* fIntTot = fTotal.createIntegral(MassDT,NormSet(MassDT),Range("lowBR,highBR"));
	RooAbsReal* fIntTot2 = fTotal.createIntegral(MassDT,NormSet(MassDT),Range("fullR"));
	RooAbsReal* fIntTot3 = fTotal.createIntegral(MassDT);


	Double_t N 	= hist->Integral(hist->FindBin(Mean-NSig*Sigma),hist->FindBin(Mean+NSig*Sigma));
	//cout << "N 1 " << N << endl;
	N = DT_set.sumEntries(0,"sigR");
	//cout << "N 2 " << N << endl;
	//Double_t Bg = hist->Integral(hist->FindBin(Mean-2*NSig*Sigma),hist->FindBin(Mean-NSig*Sigma));
	//		Bg += hist->Integral(hist->FindBin(Mean+NSig*Sigma),hist->FindBin(Mean+2*NSig*Sigma));

	Double_t Bg = DT_set.sumEntries(0,"lowBR");
	//cout << Bg << endl;
	Bg += DT_set.sumEntries(0,"highBR");
	//cout << "sb bg " << Bg << endl;

	double a = Mean-NSig*Sigma;
	double b = Mean+NSig*Sigma;
	double c = Mean-2.*NSig*Sigma;
	double d = Mean+2.*NSig*Sigma;
	Double_t SBscale = fbg->Integral(a,b);
	//cout << "sb scale is " << SBscale << endl;
	SBscale /= (fbg->Integral(c,a)+fbg->Integral(b,d)) > 0 ? (fbg->Integral(c,a)+fbg->Integral(b,d)) : 1;
	//cout << "sb scale is " << SBscale << endl;
	//SBscale = 1.;
	//if (hSidebandSF[spNumber]->GetBinContent(binNumber)==0) {
	//		hSidebandSF[spNumber]->SetBinContent(binNumber,SBscale);
//			hSidebandSF[spNumber]->SetBinError(binNumber,0.1*SBscale);
//	}
	//Bg *= SBscale;
	//Bg *= fInt->getVal();
	//cout << fInt->getVal() << endl;
	//Bg /= fIntTot->getVal();
	//cout << fIntTot->getVal() << endl;
	//cout << fIntTot2->getVal() << endl;
	//cout << fIntTot3->getVal() << endl;
	
	//cout << "bg val " << Bg << " compared to " << hist->Integral(hist->FindBin(Mean-2.*NSig*Sigma),hist->FindBin(Mean-NSig*Sigma)) + hist->Integral(hist->FindBin(Mean+NSig*Sigma),hist->FindBin(Mean+2.*NSig*Sigma)) << endl;
	//cout << "SF " << SF << endl;
	Bg *= SF;
	val[0] = N - Bg;
	val[1] = TMath::Sqrt(TMath::Abs(N + Bg));

	return val;
}

/*Double_t* MyAnalysisV0extract::ExtractYieldSBPtTree(TTree* hist) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.05, fitMax = 0.05;

	TString histName(hist->GetName());
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();
	TString muName(histName(histName.Index("iMu")+3,1));
	Int_t muNumber = muName.Atoi();

	Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

	Double_t Mean = hSidebandMean[spNumber]->GetBinContent(binNumber);
	Double_t Sigma = hSidebandSigma[spNumber]->GetFunction("pol1")->Eval(XBINS[binNumber-1]);

	Double_t NSig = 6;
	Double_t N 	= hist->Integral(hist->FindBin(Mean-NSig*Sigma),hist->FindBin(Mean+NSig*Sigma));
	Double_t Bg = hist->Integral(hist->FindBin(Mean-2*NSig*Sigma),hist->FindBin(Mean-NSig*Sigma));
			Bg += hist->Integral(hist->FindBin(Mean+NSig*Sigma),hist->FindBin(Mean+2*NSig*Sigma));

	val[0] = N - Bg;
	val[1] = TMath::Sqrt(N);


	cFits[canCounter/nBins]->cd(1+canCounter%nBins);
	mHandler->MakeNiceHistogram(hist,kBlack);
	hist->SetMarkerSize(0.5); hist->GetXaxis()->SetLabelSize(0.055);
	hist->Rebin(8);
	hist->Draw();
	hist->GetXaxis()->SetRangeUser(fitMin,fitMax);
	cFits[canCounter/nBins]->Update();
	mHandler->DrawCut(Mean+8.*Sigma,1,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins));
	mHandler->DrawCut(Mean+4.*Sigma,2,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins));
	mHandler->DrawCut(Mean-4.*Sigma,1,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins));
	mHandler->DrawCut(Mean-8.*Sigma,2,cFits[canCounter/nBins]->GetPad(1+canCounter%nBins));
	
	TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.10, 1.);
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",xBins[binNumber-1],xBins[binNumber])," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	canCounter++;
	return val;
}*/
/*Double_t* MyAnalysisV0extract::ExtractYieldSB(TH1F* hist) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	Double_t hmax = hist->GetMaximum();

	Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	hist->Rebin(8);
	RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

	TString histName(hist->GetName());
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();
	TString muName(histName(histName.Index("iMu")+3,1));
	Int_t muNumber = muName.Atoi();

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.005,0.005);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	//RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.03);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-2,2);
	RooChebychev fPolBg = RooChebychev("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA));
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);

	RooAddPdf fGaus("fGaus","fGaus",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
	RooRealVar nGaus("nGaus","N_{Gaus}",0.5*hmax,0,1e08);

	RooAddPdf fTotal = RooAddPdf("fTotal","fTotal",RooArgList(fGaus,fPolBg),RooArgList(nGaus,nPolBg));

	RooFitResult* fR = 0; 
	//if (!empty) fR = fTotal.chi2FitTo(DT_set,Save(),PrintLevel(-1));
	if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));

	// GENERATE RMS
	RooDataSet* histSigma = fGaus.generate(MassDT,100000);
	Double_t Sigma = histSigma->mean(MassDT);
	std::cout << "" << Sigma << std::endl;
		
	//	RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));
	//RooFormulaVar nGaus("nGaus","nGaus1",RooArgList(nGaus1));

	cFits[canCounter/nBins]->cd(1+canCounter%nBins);
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) {
		fTotal.plotOn(plot1,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotal.plotOn(plot1,Components(fGaus1),LineStyle(2),LineWidth(2),LineColor(kRed));
		//if (Type==0) fTotal.plotOn(plot1,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
		fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed)); }
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();

	Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);}

	TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.10, 1.);
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",xBins[binNumber-1],xBins[binNumber])," ");
	leg1->AddEntry((TObject*)0,Form("%s , #chi^{2}/ndf = %4.2f",SPECNAMES[spNumber],chi2ndf)," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	canCounter++;

	return val;
}*/

Double_t* MyAnalysisV0extract::ExtractYieldFit(TH1F* hist, Int_t Type, Int_t MB) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	Double_t hmax = hist->GetMaximum();

	Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	hist->Rebin(8);
	RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

	TString histName(hist->GetName());
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();
	TString muName(histName(histName.Index("iMu")+3,1));
	Int_t muNumber = muName.Atoi();

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.003,0.003);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.03);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-2,2);
	RooChebychev fPolBg
	 = (spNumber==2) ? RooChebychev("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA,pPolBgB))
	 : RooChebychev("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA));//,pPolBgB));
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);


	

	if (spNumber==1) nPolBg.setRange(0.,hmax);

	//constrain parameters
	/*Double_t fixv;
	fixv = fPol1_0[spNumber]->Eval(4.);
	//pGaus1A.setVal(fixv);
	//pGaus1A.setConstant(kTRUE);
	fixv = fPol1_1[spNumber]->Eval(XBINS[binNumber-1]);
	pGaus1B.setVal(fixv);
	if (muNumber>2) pGaus1B.setConstant(kTRUE);
	fixv = fPol1_2[spNumber]->Eval(XBINS[binNumber-1]);
	pGaus2B.setVal(fixv);
	if (muNumber>2) pGaus2B.setConstant(kTRUE);

	fixv = fPol1_3[spNumber]->Eval(XBINS[binNumber-1]);

	nGaus1.setVal(fixv); 
	if (muNumber>2) nGaus1.setConstant(kTRUE);
	nGaus2.setVal(1.-fixv); 
	if (muNumber>2) nGaus2.setConstant(kTRUE);*/

	RooAddPdf fGaus("fGaus","fGaus",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
	RooRealVar nGaus("nGaus","N_{Gaus}",0.5*hmax,0,1e08);

	RooAddPdf fTotal = (!Type) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus,fPolBg),RooArgList(nGaus,nPolBg))
					: RooAddPdf("fTotal","fTotal",RooArgList(fGaus),RooArgList(nGaus));


	RooFitResult* fR = 0; 
	//if (!empty) fR = fTotal.chi2FitTo(DT_set,Save(),PrintLevel(-1));
	if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
		
	//	RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));
	//RooFormulaVar nGaus("nGaus","nGaus1",RooArgList(nGaus1));

	cFits[canCounter/nBins]->cd(1+canCounter%nBins);
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) {
		fTotal.plotOn(plot1,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotal.plotOn(plot1,Components(fGaus1),LineStyle(2),LineWidth(2),LineColor(kRed));
		if (Type==0) fTotal.plotOn(plot1,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
		fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed)); }
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();

	Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);}

	if (MB==1 && !empty && 0) {
		Double_t fitparam1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getVal() : pGaus2B.getVal();
		//std::cout << "a: " << pGaus1B.getVal() << " b " << pGaus2B.getVal() << " ergo " << (pGaus1B.getVal() <= pGaus2B.getVal()) << " f " << fitparam1 << std::endl;
		Double_t fitparerr1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getError() : pGaus2B.getError(); 
		Double_t fitparam2= (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getVal() : pGaus1B.getVal();
		Double_t fitparerr2 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getError() : pGaus1B.getError(); 
		Double_t fitparam3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getVal()/nGaus.getVal() : nGaus2.getVal()/nGaus.getVal();
		Double_t fitparerr3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getError()/nGaus.getVal() : nGaus2.getError()/nGaus.getVal();

		hFitParam1[spNumber]->SetBinContent(binNumber,fitparam1);
		hFitParam1[spNumber]->SetBinError(binNumber,fitparerr1);
		hFitParam2[spNumber]->SetBinContent(binNumber,fitparam2);
		hFitParam2[spNumber]->SetBinError(binNumber,fitparerr2);
		hFitParam3[spNumber]->SetBinContent(binNumber,fitparam3);
		hFitParam3[spNumber]->SetBinError(binNumber,fitparerr3);
	}

	TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.10, 1.);
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",xBins[binNumber-1],xBins[binNumber])," ");
	leg1->AddEntry((TObject*)0,Form("%s , #chi^{2}/ndf = %4.2f",SPECNAMES[spNumber],chi2ndf)," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	canCounter++;


	//if (empty) return val;

	
	
	return val;
}

Double_t* MyAnalysisV0extract::ExtractYieldFitRt(TTree* tree, Int_t Type) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	Bool_t empty = (tree->GetEntries(Form("MassDT>%f&&MassDT<%f",fitMin,fitMax)) == 0);

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	RooDataSet DT_set("DT_set","DT_set",MassDT,Import(*tree)); 

	TString histName(tree->GetName());
	std::cout << histName.Data() << std::endl;
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();
	TString regName(histName(histName.Index("iReg")+4,1));
	Int_t Reg = regName.Atoi();
	TString ptbinName(histName(histName.Index("iPtBin")+6,1));
	Int_t PtBin = ptbinName.Atoi();
	Bool_t isTempl = (binNumber==nchmax);

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.003,0.003);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.02);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
	
	if (!isTempl) {
		if (!parRt0[spNumber][Type][Reg][PtBin]) {
			printf("Templates not yet created! Crashing.\n");
			return 0;
		}
		pGaus1A.setVal(parRt0[spNumber][Type][Reg][PtBin]);
		pGaus1A.setConstant(kTRUE);
		pGaus1B.setVal(parRt1[spNumber][Type][Reg][PtBin]);
		pGaus1B.setConstant(kTRUE);
		pGaus2B.setVal(parRt2[spNumber][Type][Reg][PtBin]);
		pGaus2B.setConstant(kTRUE);
		nGaus1.setVal(parRt3[spNumber][Type][Reg][PtBin]);
		nGaus1.setConstant(kTRUE);
		std::cout << " " << parRt3[spNumber][Type][Reg][PtBin] << " " << 1.-parRt3[spNumber][Type][Reg][PtBin] << std::endl;
		nGaus2.setVal(1.-parRt3[spNumber][Type][Reg][PtBin]);
		nGaus2.setConstant(kTRUE);
	}

	RooAddPdf fGaus("fGaus","fGaus",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
	RooRealVar nGaus("nGaus","N_{Gaus}",1,0,1e08);

	RooAddPdf fTotal = (!Type) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus,fPolBg),RooArgList(nGaus,nPolBg))
					: RooAddPdf("fTotal","fTotal",RooArgList(fGaus),RooArgList(nGaus));

	RooFitResult* fR = 0; 
	if (!empty) fR = fTotal.fitTo(DT_set,NumCPU(4),Save(),PrintLevel(-1));
		
	if (isTempl) {
		parRt0[spNumber][Type][Reg][PtBin] = pGaus1A.getVal();
		parRt1[spNumber][Type][Reg][PtBin] = pGaus1B.getVal();
		parRt2[spNumber][Type][Reg][PtBin] = pGaus2B.getVal();
		parRt3[spNumber][Type][Reg][PtBin] = nGaus1.getVal()/(nGaus1.getVal()+nGaus2.getVal());
		
	}
	//RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));

	std::cout << "cc " << canCounterRt << " den " << (nchmax-1)/increm+2 << " can " << cFitsRt[canCounterRt/((nchmax-1)/increm+2)] << std::endl; 
	cFitsRt[canCounterRt/((nchmax-1)/increm+2)]->cd(1+canCounterRt%((nchmax-1)/increm+2));
	
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) {
		fTotal.plotOn(plot1,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotal.plotOn(plot1,Components(fGaus1),LineStyle(2),LineWidth(2),LineColor(kRed));
		if (Type==0) fTotal.plotOn(plot1,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
		fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed)); }
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();

	Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);
	}

	TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.10, 1.);
	//Int_t rtleft	= (binNumber<NRTBINS) ? binNumber: 0;
	//Int_t rtright	= (binNumber<NRTBINS) ? binNumber+1 : NRTBINS;
	if (binNumber != nchmax) {
		leg1->AddEntry((TObject*)0,Form("%4.2f < N_ch < %4.2f",binNumber-1E-4, binNumber-1E-4+increm)," "); }
	else {
		leg1->AddEntry((TObject*)0,Form("%4.2f < N_ch < %4.2f",-1E-4, (Double_t)nchmax )," "); }
	leg1->AddEntry((TObject*)0,Form("%s , #chi^{2}/ndf = %4.2f",SPECNAMES[spNumber],chi2ndf)," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	canCounterRt++;


	//if (empty) return val;

	
	
	return val;
}

Double_t* MyAnalysisV0extract::ExtractYieldFitPtTree(TTree* tree, Int_t Type) {

	/*static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	
	Bool_t empty = (tree->GetEntries(Form("MassDT>%f&&MassDT<%f",fitMin,fitMax)) == 0);
	Double_t hmax = tree->GetEntries(Form("MassDT>%f&&MassDT<%f",0.5*fitMin,0.5*fitMax));

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	RooDataSet DT_set("DT_set","DT_set",MassDT,Import(*tree)); 

	TString histName(tree->GetName());
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.003,0.003);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.03);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-2,2);
	RooChebychev fPolBg
	 = (spNumber==2) ? RooChebychev("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA,pPolBgB))
	 : RooChebychev("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA));
	//RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
		
	if (spNumber==1) nPolBg.setRange(0.,hmax);

	Double_t fixv;
	fixv = fPol1_0[spNumber]->Eval(4.);
	//pGaus1A.setVal(fixv);
	//pGaus1A.setConstant(kTRUE);
	fixv = fPol1_1[spNumber]->Eval(XBINS[binNumber-1]);
	pGaus1B.setVal(fixv);
	pGaus1B.setConstant(kTRUE);
	fixv = fPol1_2[spNumber]->Eval(XBINS[binNumber-1]);
	pGaus2B.setVal(fixv);
	pGaus2B.setConstant(kTRUE);
	fixv = fPol1_3[spNumber]->Eval(XBINS[binNumber-1]);
	nGaus1.setVal(fixv); 
	nGaus1.setConstant(kTRUE);
	nGaus2.setVal(1.-fixv); 
	nGaus2.setConstant(kTRUE);

	RooAddPdf fGaus("fGaus","fGaus",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
	RooRealVar nGaus("nGaus","N_{Gaus}",0.5*hmax,0,1e08);

	RooAddPdf fTotal = (!Type) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus,fPolBg),RooArgList(nGaus,nPolBg))
					: RooAddPdf("fTotal","fTotal",RooArgList(fGaus),RooArgList(nGaus));


	//	RooAddPdf fTotal = (!Type) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg))
	//					: RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));

	RooFitResult* fR = 0; 
	if (!empty) fR = fTotal.fitTo(DT_set,NumCPU(4),Save(),PrintLevel(-1));
		
	//RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));

	//std::cout << "canCounter/(NPTBINS2) " << canCounter/(NPTBINS2) << " can " << cFitsPtTree[canCounter/(NPTBINS2)] << std::endl;
	cFitsPtTree[canCounter/(NPTBINS)]->cd(1+canCounter%(NPTBINS));
	
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) {
		fTotal.plotOn(plot1,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotal.plotOn(plot1,Components(fGaus1),LineStyle(2),LineWidth(2),LineColor(kRed));
		if (!Type) fTotal.plotOn(plot1,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
		fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed)); }
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();

	//TString histName(tree->GetName());
	//std::cout << histName.Data() << std::endl;
	//TString binName(histName(histName.Index("iBin")+4,2));
	//Int_t binNumber = binName.Atoi();
	//TString spName(histName(histName.Index("iSp")+3,1));
	//Int_t spNumber = spName.Atoi();
	//Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);
	}

	TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.10, 1.);
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",XBINS2[binNumber],XBINS2[binNumber+1])," ");
	leg1->AddEntry((TObject*)0,Form("%s , #chi^{2}/ndf = %4.2f",SPECNAMES[spNumber],chi2ndf)," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	canCounter++;


	//if (empty) return val;

	
	
	return val;*/
}

void MyAnalysisV0extract::DoClosureTest(Int_t opt) {

	// DIVIDING BLINDLY REC. MC DATA BY PDG-ID'D PT SPECTRA
	// NUM HAS FEEDDOWN, DEN DOESN'T -- really though ?


	mHandler->root()->SetBatch(kTRUE);
	for (Int_t iSp = 1; iSp < NSPECIES; iSp++)		{
		Int_t iMu = 0; Int_t iSph = 0;	
		hClosureTest[iSp]	= (TH1F*)hV0PtFit[iSp][0][iMu][iSph]->Clone(Form("hClosureTest_%s",SPECIES[iSp]));
		TH1D* hDen = (TH1D*)hV0Pt[iSp][1][iMu][iSph]->Clone(Form("hDen"));
		hDen->Scale(1.,"width");

		mHandler->MakeNiceHistogram(hClosureTest[iSp],kBlack);
		hClosureTest[iSp]->GetYaxis()->SetTitle("blind rec. / PID identified");
		hClosureTest[iSp]->GetYaxis()->SetRangeUser(0.7,1.3);
		hClosureTest[iSp]->Divide(hDen);
		TCanvas* cClosure = new TCanvas("cClosure","",900,900);
		hClosureTest[iSp]->Draw();
		if (iSp == 3) hClosureTest[3]->SetMarkerStyle(24);
		hClosureTest[iSp]->Draw();
		if (iSp == 3) hClosureTest[2]->Draw("same");
		TLegend* leg1 = new TLegend(0.171,0.67,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
		mHandler->MakeNiceLegend(leg1, 0.05, 1.);
		if (iSp == 3) leg1->AddEntry(hClosureTest[2],Form("%s",SPECNAMES[2]),"pl");
		leg1->AddEntry(hClosureTest[iSp],Form("%s",SPECNAMES[iSp]),"pl");
		leg1->Draw();
		cClosure->SaveAs(Form("plots/closure_%s.png",SPECIES[iSp]));
		
		delete hDen;


		
	}

	for (Int_t iSp = 1; iSp < NSPECIES; iSp++)		{
		Int_t iMu = 0; Int_t iSph = 0;	
		hClosureTest[iSp]	= (TH1F*)hV0PtFitSecondaryPDG[iSp]->Clone(Form("hClosureTestSigvPrSe_%s",SPECIES[iSp]));
		//hClosureTest[iSp]->Add(hV0PtFitSecondaryPDG[iSp],1.);
		TH1F* hDenPr = (TH1F*)hV0IMvPtSecondaryPDG[iSp]->ProjectionX(Form("hDenPr_%s",SPECIES[iSp]),0,-1);
		//TH1F* hDenSe = (TH1F*)hV0IMvPtSecondaryPDG[iSp]->ProjectionX(Form("hDenSe_%s",SPECIES[iSp]),0,-1);
		//hDenPr->Add(hDenSe);
		hDenPr->Scale(1,"width");

		mHandler->MakeNiceHistogram(hClosureTest[iSp],kBlack);
		//hClosureTest[iSp]->GetYaxis()->SetTitle("blind sig. ex. / PrimaryPDG + SecondaryPDG");
		hClosureTest[iSp]->GetYaxis()->SetTitle("SecondaryPDG signal ex. efficiency");
		hClosureTest[iSp]->GetYaxis()->SetRangeUser(0.7,1.3);
		hClosureTest[iSp]->Divide(hDenPr);
		TCanvas* cClosure = new TCanvas("cClosureSigvPrSe","",900,900);
		hClosureTest[iSp]->Draw();
		if (iSp == 3) hClosureTest[3]->SetMarkerStyle(24);
		hClosureTest[iSp]->Draw();
		if (iSp == 3) hClosureTest[2]->Draw("same");
		TLegend* leg1 = new TLegend(0.171,0.67,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
		mHandler->MakeNiceLegend(leg1, 0.05, 1.);
		if (iSp == 3) leg1->AddEntry(hClosureTest[2],Form("%s",SPECNAMES[2]),"pl");
		leg1->AddEntry(hClosureTest[iSp],Form("%s",SPECNAMES[iSp]),"pl");
		leg1->Draw();
		cClosure->SaveAs(Form("plots/closureSigvPrSe_%s.png",SPECIES[iSp]));
		
		//delete hDen;


		
	}

	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0extract::DrawPad(Int_t Sp, Int_t Type) {

	mHandler->root()->SetBatch(kTRUE);
	TCanvas* can = (TCanvas*)cFits[0];
	can->Draw();
	can->SaveAs("plots/fit2.png");
	mHandler->root()->SetBatch(kFALSE);
}

