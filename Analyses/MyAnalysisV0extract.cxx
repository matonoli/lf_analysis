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
#include <TTree.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include "MyAnalysisV0extract.h"
#include "../MyEvent.h"
#include "../MyTrack.h"
#include "../MyParticle.h"
#include "../MyV0.h"
#include "../MyHandler.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooChi2Var.h"

//#include <AliAnalysisPIDV0.h>
using namespace V0consts;
using namespace RooFit;

ClassImp(MyAnalysisV0extract)

MyAnalysisV0extract::MyAnalysisV0extract() {

}

Int_t MyAnalysisV0extract::Init() {

	TString dfName(this->GetName());
	dfName = Form("%s_%i",dfName.Data(),mHandler->nAnalysis());
	mDirFile = new TDirectoryFile(dfName,dfName,"",mHandler->file());
	mDirFile->cd();

	printf("Initialising analysis %s  \n", 
		this->GetName());

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

	hNchTrans	= (TH1D*)mHandler->analysis(0)->dirFile()->Get("hNchTrans");

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0IMvPt[iSp][iType][iMu][iSph] 
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));



	} } } }

	
	for (Int_t iSp = 1; iSp < NSPECIES; iSp++)		{
	for (Int_t iType = 0; iType < nType; iType++)	{
	for (Int_t iReg = 0; iReg < NREGIONS; iReg++)	{	
		
		//cout << iSp << " " << iType << " " << iReg << " " << Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]) << endl;
		
		tV0massRt[iSp][iType][iReg] = (TNtuple*)mHandler->analysis(0)->dirFile()->Get(Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		//cout << "tV0massRt[iSp][iType][iReg] " << tV0massRt[iSp][iType][iReg] << endl;
		//tV0massRt[iSp][iType][iReg] = (TNtuple*)mHandler->analysis(0)->dirFile()->Get(Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));

	}	}	}


	

}

Bool_t MyAnalysisV0extract::CreateHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{
		hRtV0Yields[iType][iReg][iPtBin]	= new TH1D(Form("hRtV0Yields_%s_%s_%i",TYPE[iType],REGIONS[iReg],iPtBin),";Species;R_{T} integrated yield",5,-0.5,4.5);

	}	}	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 
		hV0PtFit[iSp][iType][iMu][iSph]		= new TH1D(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);

		if (hV0PtFit[iSp][iType][iMu][iSph] && mHandler->IsRebinPt()) {
			TString tmpName = hV0PtFit[iSp][iType][iMu][iSph]->GetName();
			TH1D* htmp = (TH1D*)hV0PtFit[iSp][iType][iMu][iSph]->Rebin(NPTBINS2,"htmp",XBINS2);
			delete hV0PtFit[iSp][iType][iMu][iSph];
			hV0PtFit[iSp][iType][iMu][iSph] = (TH1D*)htmp->Clone(tmpName.Data());
			delete htmp;
		}

	} } } }

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iPtBin = 0; iPtBin < NRTPTBINS; ++iPtBin)	{

		printf("initing hV0RtFit_%s_%s_%s_%i \n",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iPtBin);
		hV0RtFit[iSp][iType][iReg][iPtBin] = new TH1D(Form("hV0RtFit_%s_%s_%s_%i",SPECIES[iSp],TYPE[iType],REGIONS[iReg],iPtBin),
			";R_{T}; Yield",								50, -0.5, 49.5);

	}	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{

		hFitParam0[iSp]	= new TH1D(Form("hFitParam0_%s", SPECIES[iSp]),"; p_{T} bin; gaus1 #mu;",NPTBINS,XBINS);
		hFitParam1[iSp]	= new TH1D(Form("hFitParam1_%s", SPECIES[iSp]),"; p_{T} bin; gaus1 #sigma;",NPTBINS,XBINS);
		hFitParam2[iSp]	= new TH1D(Form("hFitParam2_%s", SPECIES[iSp]),"; p_{T} bin; gaus2 #sigma;",NPTBINS,XBINS);
		hFitParam3[iSp]	= new TH1D(Form("hFitParam3_%s", SPECIES[iSp]),"; p_{T} bin; gaus1 fraction;",NPTBINS,XBINS);
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)		{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{
			
		hV0PtRtFit[iSp][iType][iReg][iRtBin]		= new TH1D(Form("hV0PtRtFit_%s_%s_%s_%1.1f-%1.1f",SPECIES[iSp],TYPE[iType],REGIONS[iReg],RTBINS0[iRtBin],RTBINS0[iRtBin+1]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS2,XBINS2);

	} } } }

}


Int_t MyAnalysisV0extract::Finish() {

	printf("Finishing analysis %s \n",this->GetName());

	mDirFile->cd();

	//StudyIMShapeRC();

	GetTemplates();

	ProducePtSpectraFromHists();
	//ProducePtSpectraFromTrees();
	//ProduceRtSpectraFromTrees();

	//DrawConstraints();
	
	return 0;	
}

void MyAnalysisV0extract::DrawConstraints() { 

	TCanvas* cParams = new TCanvas("cParams","",1000,900);
	mHandler->MakeNiceHistogram(hFitParam1[1],1);
	mHandler->MakeNiceHistogram(hFitParam2[1],2);
	hFitParam1[1]->GetYaxis()->SetRangeUser(0.0,0.03);
	hFitParam1[1]->Draw();
	hFitParam2[1]->Draw("same");
}


void MyAnalysisV0extract::GetTemplates() {


	TCanvas* cPars[NSPECIES];
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		int iType = 0; int iMu = 0; int iSph = 0;

		Float_t fitMin = -0.03, fitMax = 0.03;

		for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1)	{

			TH1D* hist = (TH1D*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, iBin),
				iBin,iBin);
			
			Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

			RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
			hist->Rebin(4);
			RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

			RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
			RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
			RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
			RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
				
			RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
			RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.03);
			RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
			RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
				
			RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
			RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-2,2);
			RooChebychev fPolBg("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA));//,pPolBgB));
			RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
				
			RooAddPdf fTotal = (!iType) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg))
							: RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
			

			RooFitResult* fR = 0; 
			//if (!empty) fR = fTotal.chi2FitTo(DT_set,Save(),PrintLevel(-1));
			if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
				
			RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));

			/*cFits[canCounter/nBins]->cd(1+canCounter%nBins);
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
			plot1->Draw();*/

			//Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

			/*if (!empty) {
				val[0] = nGaus.getVal();
				//printf("STATUS: int from fit is %f \n", val[0]);
				val[1] = nGaus.getPropagatedError(*fR);}*/

			if (!empty) {

				Double_t fitparam0 = pGaus1A.getVal();
				Double_t fitparerr0 = pGaus1A.getError();
				Double_t fitparam1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getVal() : pGaus2B.getVal();
				Double_t fitparerr1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getError() : pGaus2B.getError(); 
				Double_t fitparam2= (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getVal() : pGaus1B.getVal();
				Double_t fitparerr2 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getError() : pGaus1B.getError(); 
				Double_t fitparam3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getVal()/nGaus.getVal() : nGaus2.getVal()/nGaus.getVal();
				Double_t fitparerr3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getError()/nGaus.getVal() : nGaus2.getError()/nGaus.getVal();


				hFitParam0[iSp]->SetBinContent(iBin,fitparam0);
				hFitParam0[iSp]->SetBinError(iBin,fitparerr0);
				hFitParam1[iSp]->SetBinContent(iBin,fitparam1);
				hFitParam1[iSp]->SetBinError(iBin,fitparerr1);
				hFitParam2[iSp]->SetBinContent(iBin,fitparam2);
				hFitParam2[iSp]->SetBinError(iBin,fitparerr2);
				hFitParam3[iSp]->SetBinContent(iBin,fitparam3);
				hFitParam3[iSp]->SetBinError(iBin,fitparerr3);
			}
		}

		cPars[iSp]	= new TCanvas(Form("cPars_%s",SPECIES[iSp]),Form("cPars_%s",SPECIES[iSp]),1000,900);
		cPars[iSp]->Divide(2,2,1e-03,1e-03);

		fPol1_0[iSp] = new TF1(Form("fPol1_0_%s",SPECIES[iSp]),"[1]",XBINS[26],XBINS[51]);
		fPol1_1[iSp] = new TF1(Form("fPol1_1_%s",SPECIES[iSp]),"[0]*x+[1]",XBINS[21],XBINS[51]);
		fPol1_2[iSp] = new TF1(Form("fPol1_2_%s",SPECIES[iSp]),"[1]",XBINS[0],XBINS[NPTBINS]);
		fPol1_3[iSp] = new TF1(Form("fPol1_3_%s",SPECIES[iSp]),"[0]*x+[1]",XBINS[0],XBINS[NPTBINS]);
		
		Double_t boundsPar0[NPTBINS];
		for (int i=0; i<NPTBINS; i++) boundsPar0[i] = 0;
		Double_t boundsPar1[NPTBINS];
		for (int i=0; i<NPTBINS; i++) boundsPar1[i] = 0;
		Double_t boundsPar2[NPTBINS];
		for (int i=0; i<NPTBINS; i++) boundsPar2[i] = 0;
		Double_t boundsPar3[NPTBINS];
		for (int i=0; i<NPTBINS; i++) boundsPar3[i] = 0;
		


		cPars[iSp]->cd(1);
		TFitResultPtr fres = hFitParam0[iSp]->Fit(fPol1_0[iSp],"S R");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar0,0.95, false);
		mHandler->MakeNiceHistogram(hFitParam0[iSp],1);
		hFitParam0[iSp]->Draw("");
		fPol1_0[iSp]->Draw("same");

		cPars[iSp]->cd(2);
		fres = hFitParam1[iSp]->Fit(fPol1_1[iSp],"S");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar1,0.95, true);
		//cout << "b " << boundsPar1[5] << " " << boundsPar1[20] << " " << boundsPar1[35] << endl;
		//for (int i=0; i<NPTBINS; i++) boundsPar1[i] += fPol1_1->Eval(XBINS[i]);
		TGraph* tg = new TGraph(NPTBINS,XBINS,boundsPar1);
		tg->Write();

		mHandler->MakeNiceHistogram(hFitParam1[iSp],1);
		hFitParam1[iSp]->Draw("");
		fPol1_1[iSp]->Draw("same");
		//tg->Draw("same");


		cPars[iSp]->cd(3);
		fres = hFitParam2[iSp]->Fit(fPol1_2[iSp],"W S");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar2,0.95, true);
		
		mHandler->MakeNiceHistogram(hFitParam2[iSp],1);
		hFitParam2[iSp]->Draw("");
		fPol1_2[iSp]->Draw("same");
		//tg->Draw("same");


		cPars[iSp]->cd(4);
		fres = hFitParam3[iSp]->Fit(fPol1_3[iSp],"W S");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar3,0.95, true);
		
		mHandler->MakeNiceHistogram(hFitParam3[iSp],1);
		hFitParam3[iSp]->Draw("hist p");
		fPol1_3[iSp]->Draw("same");

		cPars[iSp]->Write();

		for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1)	{
		//Int_t iBin = 0; if (0) {
			TH1D* hist = (TH1D*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, iBin),
				iBin,iBin);
			
			Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

			RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
			hist->Rebin(8);
			RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

			RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
			RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
			RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
			RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
				
			RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
			RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.03);
			RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
			RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
				
			RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
			RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-2,2);
			RooChebychev fPolBg("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA,pPolBgB));//,pPolBgB));
			RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
				
			TString histName(hist->GetName());
			TString binName(histName(histName.Index("iBin")+4,2));
			Int_t binNumber = binName.Atoi();
			TString spName(histName(histName.Index("iSp")+3,1));
			Int_t spNumber = spName.Atoi();

			Double_t hmax = hist->GetMaximum();

			if (spNumber==1) nPolBg.setRange(0.,hmax);

			//constrain parameters
			Double_t fixv;
			fixv = fPol1_0[spNumber]->Eval(4.);
			//pGaus1A.setVal(fixv);
			//pGaus1A.setConstant(kTRUE);
			fixv = fPol1_1[spNumber]->Eval(XBINS[binNumber-1]);
			pGaus1B.setVal(fixv);
			pGaus1B.setConstant(kTRUE);
			fixv = fPol1_2[spNumber]->Eval(XBINS[binNumber-1]);
			pGaus2B.setVal(fixv);
			//pGaus2B.setConstant(kTRUE);

			fixv = fPol1_3[spNumber]->Eval(XBINS[binNumber-1]);


			RooAddPdf fTotal = (!iType) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg))
							: RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
			
			
			//fixv = fPol1_1->Eval(hist->GetBinCenter(iBin));
			//fixe = 1*(boundsPar1[iBin-1]);//+boundsPar1[iBin]);
			//cout << "fixe " << fixe << endl;
			//pGaus1B.setVal(fixv);
			//pGaus1B.setRange(fixv-fixe,fixv+fixe);
			//pGaus1A.setRange(fPol1_0->Eval(hist->GetBinCenter(iBin)) - (boundsPar0[30]),
			//	fPol1_0->Eval(hist->GetBinCenter(iBin)) + boundsPar0[30]);

			RooFitResult* fR = 0; 
			if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
				
			RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));


			if (!empty) {

				Double_t fitparam0 = pGaus1A.getVal();
				Double_t fitparerr0 = pGaus1A.getError();
				Double_t fitparam1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getVal() : pGaus2B.getVal();
				Double_t fitparerr1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getError() : pGaus2B.getError(); 
				Double_t fitparam2= (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getVal() : pGaus1B.getVal();
				Double_t fitparerr2 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getError() : pGaus1B.getError(); 
				Double_t fitparam3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getVal()/nGaus.getVal() : nGaus2.getVal()/nGaus.getVal();
				Double_t fitparerr3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getError()/nGaus.getVal() : nGaus2.getError()/nGaus.getVal();


				//hFitParam0[iSp]->SetBinContent(iBin,fitparam0);
				//hFitParam0[iSp]->SetBinError(iBin,fitparerr0);
				//hFitParam1[iSp]->SetBinContent(iBin,fitparam1);
				//hFitParam1[iSp]->SetBinError(iBin,fitparerr1);
				hFitParam2[iSp]->SetBinContent(iBin,fitparam2);
				hFitParam2[iSp]->SetBinError(iBin,fitparerr2);
				hFitParam3[iSp]->SetBinContent(iBin,fitparam3);
				hFitParam3[iSp]->SetBinError(iBin,fitparerr3);
			}
		}


		cPars[iSp]->cd(3);
		fres = hFitParam2[iSp]->Fit(fPol1_2[iSp],"W S");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar2,0.95, true);
		
		mHandler->MakeNiceHistogram(hFitParam2[iSp],1);
		hFitParam2[iSp]->Draw("");
		fPol1_2[iSp]->Draw("same");


		for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1)	{
		//Int_t iBin = 0; if (0) {
			TH1D* hist = (TH1D*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, iBin),
				iBin,iBin);
			
			Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

			RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
			hist->Rebin(8);
			RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

			RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
			RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
			RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
			RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
				
			RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
			RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.1);
			RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
			RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
				
			RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
			RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-2,2);
			RooChebychev fPolBg("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA,pPolBgB));//,pPolBgB));
			RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
				
			TString histName(hist->GetName());
			TString binName(histName(histName.Index("iBin")+4,2));
			Int_t binNumber = binName.Atoi();
			TString spName(histName(histName.Index("iSp")+3,1));
			Int_t spNumber = spName.Atoi();

			Double_t hmax = hist->GetMaximum();

			if (spNumber==1) nPolBg.setRange(0.,hmax);

			//constrain parameters
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


			RooAddPdf fTotal = (!iType) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg))
							: RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
			
			
			//fixv = fPol1_1->Eval(hist->GetBinCenter(iBin));
			//fixe = 1*(boundsPar1[iBin-1]);//+boundsPar1[iBin]);
			//cout << "fixe " << fixe << endl;
			//pGaus1B.setVal(fixv);
			//pGaus1B.setRange(fixv-fixe,fixv+fixe);
			//pGaus1A.setRange(fPol1_0->Eval(hist->GetBinCenter(iBin)) - (boundsPar0[30]),
			//	fPol1_0->Eval(hist->GetBinCenter(iBin)) + boundsPar0[30]);

			RooFitResult* fR = 0; 
			if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
				
			RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));


			if (!empty) {

				Double_t fitparam0 = pGaus1A.getVal();
				Double_t fitparerr0 = pGaus1A.getError();
				Double_t fitparam1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getVal() : pGaus2B.getVal();
				Double_t fitparerr1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getError() : pGaus2B.getError(); 
				Double_t fitparam2= (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getVal() : pGaus1B.getVal();
				Double_t fitparerr2 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getError() : pGaus1B.getError(); 
				Double_t fitparam3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getVal()/nGaus.getVal() : nGaus2.getVal()/nGaus.getVal();
				Double_t fitparerr3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getError()/nGaus.getVal() : nGaus2.getError()/nGaus.getVal();


				//hFitParam0[iSp]->SetBinContent(iBin,fitparam0);
				//hFitParam0[iSp]->SetBinError(iBin,fitparerr0);
				//hFitParam1[iSp]->SetBinContent(iBin,fitparam1);
				//hFitParam1[iSp]->SetBinError(iBin,fitparerr1);
				//hFitParam2[iSp]->SetBinContent(iBin,fitparam2);
				//hFitParam2[iSp]->SetBinError(iBin,fitparerr2);
				hFitParam3[iSp]->SetBinContent(iBin,fitparam3);
				hFitParam3[iSp]->SetBinError(iBin,fitparerr3);
			}
		}


		cPars[iSp]->cd(4);
		fres = hFitParam3[iSp]->Fit(fPol1_3[iSp],"W S");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar3,0.95, true);
		
		mHandler->MakeNiceHistogram(hFitParam3[iSp],1);
		hFitParam3[iSp]->GetYaxis()->SetRangeUser(-0.2,1.2);
		hFitParam3[iSp]->Draw("");
		fPol1_3[iSp]->Draw("same");

		/*
		cPars[iSp]->cd(2);
		fres = hFitParam1[iSp]->Fit(fPol1_1,"S W");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar1,0.95, true);
		//cout << "b " << boundsPar1[5] << " " << boundsPar1[20] << " " << boundsPar1[35] << endl;
		//for (int i=0; i<NPTBINS; i++) boundsPar1[i] += fPol1_1->Eval(XBINS[i]);
		

		mHandler->MakeNiceHistogram(hFitParam1[iSp],1);
		hFitParam[iSp]1[iSp]->Draw("");
		fPol1_1->Draw("same");
		//tg->Draw("same");


		cPars[iSp]->cd(3);
		fres = hFitParam2[iSp]->Fit(fPol1_2,"S");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar2,0.95, true);
		
		mHandler->MakeNiceHistogram(hFitParam2[iSp],1);
		hFitParam2[iSp]->Draw("");
		fPol1_2->Draw("same");
		//tg->Draw("same");


		cPars[iSp]->cd(4);
		fres = hFitParam3[iSp]->Fit(fPol1_3,"S");
		fres->GetConfidenceIntervals(NPTBINS,1,1,XBINS,boundsPar3,0.95, true);
		
		mHandler->MakeNiceHistogram(hFitParam3[iSp],1);
		hFitParam3[iSp]->Draw("hist p");
		fPol1_3->Draw("same");

		cPars[iSp]->Write();

		for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1)	{
		//Int_t iBin = 0; if (0) {
			TH1D* hist = (TH1D*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, iBin),
				iBin,iBin);
			
			Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

			RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
			hist->Rebin(8);
			RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

			RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
			RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0005,0.01);
			RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
			RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
				
			RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
			RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.1);
			RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
			RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
				
			RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
			RooRealVar pPolBgB("pPolBgB","Pol. par. B",0,-2,2);
			RooChebychev fPolBg("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA));//,pPolBgB));
			RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
				
			RooAddPdf fTotal = (!iType) ? RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg))
							: RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
			
			Double_t fixv = fPol1_0->Eval(hist->GetBinCenter(iBin));
			Double_t fixe = 3*boundsPar0[30];// 3e-4;// = 0.5*(boundsPar0[iBin-1]+boundsPar0[iBin];
			pGaus1A.setVal(fixv);
			pGaus1A.setRange(fixv-fixe,fixv+fixe);
			fixv = fPol1_1->Eval(hist->GetBinCenter(iBin));
			fixe = 1*(boundsPar1[iBin-1]);//+boundsPar1[iBin]);
			//cout << "fixe " << fixe << endl;
			pGaus1B.setVal(fixv);
			pGaus1B.setRange(fixv-fixe,fixv+fixe);
			//pGaus1A.setRange(fPol1_0->Eval(hist->GetBinCenter(iBin)) - (boundsPar0[30]),
			//	fPol1_0->Eval(hist->GetBinCenter(iBin)) + boundsPar0[30]);

			RooFitResult* fR = 0; 
			if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
				
			RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));


			if (!empty) {

				Double_t fitparam0 = pGaus1A.getVal();
				Double_t fitparerr0 = pGaus1A.getError();
				Double_t fitparam1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getVal() : pGaus2B.getVal();
				Double_t fitparerr1 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus1B.getError() : pGaus2B.getError(); 
				Double_t fitparam2= (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getVal() : pGaus1B.getVal();
				Double_t fitparerr2 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? pGaus2B.getError() : pGaus1B.getError(); 
				Double_t fitparam3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getVal()/nGaus.getVal() : nGaus2.getVal()/nGaus.getVal();
				Double_t fitparerr3 = (pGaus1B.getVal() <= pGaus2B.getVal()) ? nGaus1.getError()/nGaus.getVal() : nGaus2.getError()/nGaus.getVal();


				hFitParam0[iSp]->SetBinContent(iBin,fitparam0);
				hFitParam0[iSp]->SetBinError(iBin,fitparerr0);
				hFitParam1[iSp]->SetBinContent(iBin,fitparam1);
				hFitParam1[iSp]->SetBinError(iBin,fitparerr1);
				hFitParam2[iSp]->SetBinContent(iBin,fitparam2);
				hFitParam2[iSp]->SetBinError(iBin,fitparerr2);
				hFitParam3[iSp]->SetBinContent(iBin,fitparam3);
				hFitParam3[iSp]->SetBinError(iBin,fitparerr3);
			}
		}




	}*/



	}

	/*load mb data histograms 
	bin loop
		make template in the bin
		-get fixed widths and ratios 
		- perhaps for simplicity actually just pol1

	save widths and gaus ratios
	in an array or vector
	


	then for fitting -> load arrays and fix parameters to those values */
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


void MyAnalysisV0extract::StudyIMShapeRC() {

	if (!mFileMC) {
		printf("No MC file loaded! MC shape study not performed.\n");
		return;
	}

	TDirectoryFile* dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");
	
	for (int iSp = 1; iSp < NSPECIES; ++iSp)		{
		tV0massRCMB[iSp]	= (TNtuple*)dirFile1->Get(Form("tV0massRCMB_%s",SPECIES[iSp]));
	}

	Int_t iSp = 1;
	TCanvas* cShapeMC = new TCanvas("cShapeMC","",1000,900);
	TTree* tree = tV0massRCMB[iSp]->CopyTree("lPt>1.0 && lPt < 2.0");

	Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.04, fitMax = 0.04;

	Bool_t empty = (tree->GetEntries(Form("MassDT>%f&&MassDT<%f",fitMin,fitMax)) == 0);

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	MassDT.setRange("RangeK1",fitMin,-0.01);
	MassDT.setRange("RangeK2",0.01,fitMax);
	RooDataSet DT_set("DT_set","DT_set",MassDT,Import(*tree)); 

	RooRealVar pGaus1A("pGaus1A","Mean 1",0.,-0.01,0.01);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.0001,0.02);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",0.,0.,0.);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.0001,0.08);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);


	RooRealVar pGaus3A("pGaus3A","Mean 3",-0.,0.,0.);
	RooRealVar pGaus3B("pGaus3B","Sigma 3",0.0001,0.05);
	RooGaussian fGaus3("fGaus3","fGaus3",MassDT,pGaus1A,pGaus3B); 
	RooRealVar nGaus3("nGaus3","N_{Gaus3}",1,0,1e08);

	/*RooRealVar pVoigA("pVoigA","Mean",0.,-100.,100.);
	RooRealVar pVoigB("pVoigB","Sigma",0.001,-1000.,1000.);
	RooRealVar pVoigC("pVoigC","Width",0.001,-1000.,1000.);
	RooVoigtian fVoig("fVoig","fVoig",MassDT,pVoigA,pVoigB,pVoigC);
	RooRealVar nVoig("nVoig","N_{Voig}",1,0,1e08);*/  // voigtian actually sucks ass here a lot
		
	//RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	//RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);
	//RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
		
	//RooAddPdf fTotal = RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fGaus3),RooArgList(nGaus1,nGaus2,nGaus3));
	RooAddPdf fTotal = RooAddPdf("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
	//RooAddPdf fTotal = RooAddPdf("fTotal","fTotal",RooArgList(fVoig),RooArgList(nVoig));

	RooFitResult* fR = 0; 
	if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
		
	RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));
	//RooFormulaVar nGaus("nGaus","nGaus1+nGaus2+nGaus3",RooArgList(nGaus1,nGaus2,nGaus3));

	//cout << "canCounter/(NPTBINS2) " << canCounter/(NPTBINS2) << " can " << cFitsPtTree[canCounter/(NPTBINS2)] << endl;
	
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) {
		//fTotal.plotOn(plot1,Components(fGaus3),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotal.plotOn(plot1,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotal.plotOn(plot1,Components(fGaus1),LineStyle(2),LineWidth(2),LineColor(kRed));
		//if (!Type) fTotal.plotOn(plot1,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
		fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed));
		RooArgSet params = RooArgList(pGaus1A, pGaus1B, pGaus2A, pGaus2B, pGaus3A, pGaus3B);
		//RooArgSet params = RooArgList(pVoigA, pVoigB, pVoigC);
		fTotal.paramOn(plot1,Parameters(params),Layout(0.1,0.4,0.9), ShowConstants(kTRUE));

		 }
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

	TLegend* leg1 = new TLegend(0.13,0.25,0.17,0.50);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.035,1);
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",1.,2.)," ");
	leg1->AddEntry((TObject*)0,Form(" #chi^{2}/ndf = %4.2f",chi2ndf)," ");
	//leg1->AddEntry((TObject*)0,Form(" f1 pars: %2.2f ; %2.2f",pGaus1A.getVal(),pGaus1B.getVal())," ");
	//leg1->AddEntry((TObject*)0,Form(" f2 pars: %2.2f ; %2.2f",pGaus2A.getVal(),pGaus2B.getVal())," ");
	//leg1->AddEntry((TObject*)0,Form(" f3 pars: %2.2f ; %2.2f",pGaus3A.getVal(),pGaus3B.getVal())," ");
	leg1->AddEntry((TObject*)0,Form("%2.1f #pm %2.1f",val[0],val[1])," ");
	leg1->AddEntry((TObject*)0,Form("fractions: %2.2f ; %2.2f ; %2.2f",nGaus1.getVal()/val[0],nGaus2.getVal()/val[0],nGaus3.getVal()/val[0])," ");
	leg1->Draw();
	

	// DATA COMPARISON MB
	////////////////////
	TCanvas* cShapeD = new TCanvas("cShapeD","",1000,900);
	TH1D* hist = (TH1D*)hV0IMvPt[iSp][0][0][0]->ProjectionY(
				"hist",		hV0IMvPt[1][0][0][0]->GetXaxis()->FindBin(1.0),
				hV0IMvPt[1][0][0][0]->GetXaxis()->FindBin(2.0));
	

	empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);
	hist->Rebin(4);
	RooDataHist DT_hist("DT_set","DT_set",MassDT,Import(*hist));

	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0.,-100.,100.);//-2,2);
	RooRealVar pPolBgB("pPolBgB","Pol. par. B",0.,-100.,100.);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,RooArgSet(pPolBgA,pPolBgB));
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);

	if (!empty) fGaus2.fitTo(DT_hist,Save(),PrintLevel(-1),Range("RangeK1,RangeK2"));
	pPolBgA.setConstant(kTRUE); pPolBgB.setConstant(kTRUE);
	//pGaus1B.setConstant(kTRUE); pGaus2B.setConstant(kTRUE); pGaus3B.setConstant(kTRUE); 
	//nGaus1.setConstant(kTRUE); nGaus2.setConstant(kTRUE); nGaus3.setConstant(kTRUE); 
	RooAddPdf fTotalD = RooAddPdf("fTotalD","fTotalD",RooArgList(fTotal,fPolBg),RooArgList(nGaus,nPolBg));
	if (!empty) fR = fTotalD.fitTo(DT_hist,Save(),PrintLevel(-1));

	RooPlot* plot2 = MassDT.frame(Title(" "));
	DT_hist.plotOn(plot2,MarkerSize(0.4));
	if (!empty) {
		//fTotalD.plotOn(plot2,Components(fGaus3),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotalD.plotOn(plot2,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotalD.plotOn(plot2,Components(fGaus1),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotalD.plotOn(plot2,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
		fTotalD.plotOn(plot2,LineWidth(2),LineColor(kRed));
		RooArgSet params = RooArgList(pGaus1A, pGaus1B, pGaus2A, pGaus2B, pGaus3A, pGaus3B,pPolBgA,pPolBgB);
		//RooArgSet params = RooArgList(pVoigA, pVoigB, pVoigC);
		fTotalD.paramOn(plot2,Parameters(params),Layout(0.1,0.4,0.9), ShowConstants(kTRUE));

		 }
	plot2->SetMinimum(1e-05);
	plot2->SetMaximum(1.40*plot2->GetMaximum());
	plot2->GetXaxis()->SetTitleSize(0.05);
	plot2->GetYaxis()->SetTitleSize(0.05);
	plot2->Draw();

	chi2ndf = (!empty) ? plot2->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);
	}

	TLegend* leg2 = new TLegend(0.13,0.25,0.17,0.50);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg2, 0.035,1);
	leg2->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",1.,2.)," ");
	leg2->AddEntry((TObject*)0,Form(" #chi^{2}/ndf = %4.2f",chi2ndf)," ");
	//leg1->AddEntry((TObject*)0,Form(" f1 pars: %2.2f ; %2.2f",pGaus1A.getVal(),pGaus1B.getVal())," ");
	//leg1->AddEntry((TObject*)0,Form(" f2 pars: %2.2f ; %2.2f",pGaus2A.getVal(),pGaus2B.getVal())," ");
	//leg1->AddEntry((TObject*)0,Form(" f3 pars: %2.2f ; %2.2f",pGaus3A.getVal(),pGaus3B.getVal())," ");
	leg2->AddEntry((TObject*)0,Form("%2.1f #pm %2.1f",val[0],val[1])," ");
	leg2->AddEntry((TObject*)0,Form("fractions: %2.2f ; %2.2f ; %2.2f",nGaus1.getVal()/val[0],nGaus2.getVal()/val[0],nGaus3.getVal()/val[0])," ");
	leg2->Draw();

	Double_t bgFrac = nPolBg.getVal()/nGaus.getVal();


	//////////////////////////////////////
	// DATA COMPARISON  --- RT inc.
	////////////////////
	TCanvas* cShapeDRt = new TCanvas("cShapeDRt","",1000,900);
	hist = (TH1D*)hV0IMvPt[iSp][0][3][0]->ProjectionY(
				"hist",		hV0IMvPt[1][0][0][0]->GetXaxis()->FindBin(1.0),
				hV0IMvPt[1][0][0][0]->GetXaxis()->FindBin(2.0));
	

	empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);
	hist->Rebin(4);
	RooDataHist DT_histRt("DT_set","DT_set",MassDT,Import(*hist));


	pGaus1A.setConstant(kTRUE); pGaus1B.setConstant(kTRUE); pGaus2B.setConstant(kTRUE); pGaus3B.setConstant(kTRUE); 
	nGaus1.setConstant(kTRUE); nGaus2.setConstant(kTRUE); nGaus3.setConstant(kTRUE); 
	//pPolBgA.setConstant(kTRUE); pPolBgB.setConstant(kTRUE);
	//nPolBg.setConstant(kTRUE);
	//RooRealVar nPolBg2("nPolBg2","nPolBg2",bgFrac)
	RooAddPdf fTotalDRt = RooAddPdf("fTotalDRt","fTotalDRt",RooArgList(fTotal,fPolBg),RooArgList(nGaus,nPolBg));
	if (!empty) fR = fTotalDRt.fitTo(DT_histRt,Save(),PrintLevel(-1));

	RooPlot* plot3 = MassDT.frame(Title(" "));
	DT_histRt.plotOn(plot3,MarkerSize(0.4));
	if (!empty) {
		//fTotalDRt.plotOn(plot3,Components(fGaus3),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotalDRt.plotOn(plot3,Components(fGaus2),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotalDRt.plotOn(plot3,Components(fGaus1),LineStyle(2),LineWidth(2),LineColor(kRed));
		fTotalDRt.plotOn(plot3,Components(fPolBg),LineStyle(1),LineWidth(2),LineColor(kBlue));
		fTotalDRt.plotOn(plot3,LineWidth(2),LineColor(kRed));
		RooArgSet params = RooArgList(pGaus1A, pGaus1B, pGaus2A, pGaus2B, pGaus3A, pGaus3B,pPolBgA,pPolBgB);
		//RooArgSet params = RooArgList(pVoigA, pVoigB, pVoigC);
		fTotalDRt.paramOn(plot3,Parameters(params),Layout(0.1,0.4,0.9), ShowConstants(kTRUE));

		 }
	plot3->SetMinimum(1e-05);
	plot3->SetMaximum(1.40*plot3->GetMaximum());
	plot3->GetXaxis()->SetTitleSize(0.05);
	plot3->GetYaxis()->SetTitleSize(0.05);
	plot3->Draw();

	chi2ndf = (!empty) ? plot3->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);
	}

	TLegend* leg3 = new TLegend(0.13,0.25,0.17,0.50);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg3, 0.035,1);
	leg3->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",1.,2.)," ");
	leg3->AddEntry((TObject*)0,Form(" #chi^{2}/ndf = %4.2f",chi2ndf)," ");
	//leg1->AddEntry((TObject*)0,Form(" f1 pars: %2.2f ; %2.2f",pGaus1A.getVal(),pGaus1B.getVal())," ");
	//leg1->AddEntry((TObject*)0,Form(" f2 pars: %2.2f ; %2.2f",pGaus2A.getVal(),pGaus2B.getVal())," ");
	//leg1->AddEntry((TObject*)0,Form(" f3 pars: %2.2f ; %2.2f",pGaus3A.getVal(),pGaus3B.getVal())," ");
	leg3->AddEntry((TObject*)0,Form("%2.1f #pm %2.1f",val[0],val[1])," ");
	leg3->AddEntry((TObject*)0,Form("fractions: %2.2f ; %2.2f ; %2.2f",nGaus1.getVal()/val[0],nGaus2.getVal()/val[0],nGaus3.getVal()/val[0])," ");
	leg3->Draw();



}

void MyAnalysisV0extract::ProducePtSpectraFromHists() {

	mHandler->root()->SetBatch(kTRUE);

	nBins 	= (mHandler->IsRebinPt()) ? NPTBINS2 : NPTBINS;
	xBins	= (mHandler->IsRebinPt()) ? XBINS2 : XBINS;
	Int_t binSize = -1 + TMath::Nint((Double_t)NPTBINS/nBins);	// this is buggy actually
	nPads = TMath::FloorNint(TMath::Sqrt(nBins));

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	iCan = 0;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < nType; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{

		if (iMu > 2 && (iSph < 3 && iSph)) continue;
		if (iMu < 3 && iSph > 2) continue; 

		printf("Extracting yield for pt spectrum iSp%i_iType%i_iMu%i_iSph%i \n",iSp,iType,iMu,iSph);

		cFits[iCan] = new TCanvas(Form("cFits_iSp%i_iType%i_iMu%i_iSph%i",iSp,iType,iMu,iSph),"",2800,2000);
		cFits[iCan]->Divide(nPads+1,nPads,0.00005,0.00005);
			
		Double_t* yield = 0;
		Int_t binCounter = 1;
		for (int iBin = 1; iBin < NPTBINS+1; iBin=iBin+1+binSize)	{
		//for (int iBin = 10; iBin < 11; ++iBin)	{
			
			yield = ExtractYieldFit((TH1D*)hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, binCounter),
				iBin,iBin+binSize),		iType,	(iType==0&&iMu==3&&iSph==0) );		// should be perhaps changed to FindBin
			
			hV0PtFit[iSp][iType][iMu][iSph]->SetBinContent(binCounter,*(yield+0));
			hV0PtFit[iSp][iType][iMu][iSph]->SetBinError(binCounter,*(yield+1));
			binCounter++;
		}

		cFits[iCan]->Write();
		iCan++;

		//return 0;
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

void MyAnalysisV0extract::ProducePtSpectraFromTrees() {

	mHandler->root()->SetBatch(kTRUE);

	nPads = TMath::Nint(TMath::Sqrt(NPTBINS2));
	Double_t rt_den = hNchTrans->GetMean();

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	iCan = 0;
	canCounter = 0;
	for (int iSp = 1; iSp < NSPECIES; ++iSp)			{
	for (int iType = 0; iType < nType; ++iType)			{
	for (int iReg = 0; iReg < NREGIONS; ++iReg)			{
	for (int iRtBin = 0; iRtBin < NRTBINS0; ++iRtBin)	{

		cFitsPtTree[iCan] = new TCanvas(Form("cFitsPtTree_iSp%i_iType%i_iReg%i_iRtBin%i",iSp,iType,iReg,iRtBin),"",2800,2000);
		cFitsPtTree[iCan]->Divide(nPads+1,nPads,0.00005,0.00005);
			
		Double_t* yield = 0;
		Float_t leftrtb =	RTBINS0[iRtBin]*rt_den;
		Float_t rightrtb =	RTBINS0[iRtBin+1]*rt_den;
		//		cout << "left " << leftrtb << " right " << rightrtb << endl;
		printf("Extracting yield for pt spectrum from trees iSp%i_iType%i_iReg%i_iRtBin%i \n",iSp,iType,iReg,iRtBin);
		TString treeName = tV0massRt[iSp][iType][iReg]->GetName();
		for (int iBin = 0; iBin < NPTBINS2; ++iBin)	{

			
			tV0massRt[iSp][iType][iReg]->SetName(Form("treept_iSp%i_iBin%i",iSp,iBin));

			yield = ExtractYieldFitPtTree((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
				Form("lPt>%f && lPt<%f && lNchTrans>=%f && lNchTrans<%f",XBINS2[iBin],XBINS2[iBin+1],
					leftrtb-1E-05,rightrtb)),     iType);

			tV0massRt[iSp][iType][iReg]->SetName(treeName);

			//cout << "yield " << *(yield+0) << endl;
			hV0PtRtFit[iSp][iType][iReg][iRtBin]->SetBinContent(iBin+1,*(yield+0));	//+1 for underflow bin
			hV0PtRtFit[iSp][iType][iReg][iRtBin]->SetBinError(iBin+1,*(yield+1));
			//hV0PtFit[iSp][iType][iMu][iSph]->SetBinError(binCounter,*(yield+1));
			//binCounter++;
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
}

void MyAnalysisV0extract::ProduceRtSpectraFromTrees() {

	mHandler->root()->SetBatch(kTRUE);

	Double_t rt_den = hNchTrans->GetMean();

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
	iCan = 0;

	nchmax = hNchTrans->FindLastBinAbove();
	cout << "nchmax " << nchmax << endl;
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


		//cout << " tV0massRt[iSp][0][0] " << tV0massRt[iSp][0][0] << endl;
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
			//cout << "ibin " << iBin << " jo " << (iBin < nchmax) << endl;
			yield = ExtractYieldFitRt((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
				Form("lPt>%f && lPt<%f && lNchTrans>%f && lNchTrans<%f",RT_PTRANGE[iPtBin][0],RT_PTRANGE[iPtBin][1],
					iBin-1E-4, iBin-1E-4+increm)),			iType);
//					RTBINS[iBin]*rt_den,RTBINS[iBin+1]*rt_den)),		0);
			
			hV0RtFit[iSp][iType][iReg][iPtBin]->SetBinContent(iBin+1,*(yield+0));
			hV0RtFit[iSp][iType][iReg][iPtBin]->SetBinError(iBin+1,*(yield+1));
		}
		hV0RtFit[iSp][iType][iReg][iPtBin]->Rebin(increm);

		/*//for entire rt range
		tV0massRt[iSp][iType][iReg]->SetName(Form("tree_iSp%i_iType%i_iReg%i_iPtBin%i_iBin%i",iSp,iType,iReg,iPtBin,nchmax));//NRTBINS));
			
		yield = ExtractYieldFitRt((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
			Form("lPt>%f && lPt<%f && lNchTrans>%f && lNchTrans<%f",RT_PTRANGE[iPtBin][0],RT_PTRANGE[iPtBin][1],
					0-1E-4, (Double_t)nchmax),			0) );
		hRtV0Yields[iType][iReg][iPtBin]->SetBinContent(1+iSp,*(yield+0));
		hRtV0Yields[iType][iReg][iPtBin]->SetBinError(1+iSp,*(yield+1));*/
		//

		tV0massRt[iSp][iType][iReg]->SetName(treeName);

		

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

	mHandler->root()->SetBatch(kFALSE);
}

Double_t* MyAnalysisV0extract::ExtractYieldFit(TH1D* hist, Int_t Type, Int_t MB) {

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
	Double_t fixv;
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
	if (muNumber>2) nGaus2.setConstant(kTRUE);

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
		//cout << "a: " << pGaus1B.getVal() << " b " << pGaus2B.getVal() << " ergo " << (pGaus1B.getVal() <= pGaus2B.getVal()) << " f " << fitparam1 << endl;
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
	cout << histName.Data() << endl;
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
		cout << " " << parRt3[spNumber][Type][Reg][PtBin] << " " << 1.-parRt3[spNumber][Type][Reg][PtBin] << endl;
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

	cout << "cc " << canCounterRt << " den " << (nchmax-1)/increm+2 << " can " << cFitsRt[canCounterRt/((nchmax-1)/increm+2)] << endl; 
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

	static Double_t val[2];
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

	//cout << "canCounter/(NPTBINS2) " << canCounter/(NPTBINS2) << " can " << cFitsPtTree[canCounter/(NPTBINS2)] << endl;
	cFitsPtTree[canCounter/(NPTBINS2)]->cd(1+canCounter%(NPTBINS2));
	
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

	/*TString histName(tree->GetName());
	//cout << histName.Data() << endl;
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();*/
	Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

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

	
	
	return val;
}

