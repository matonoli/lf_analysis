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

	cout << "nreg " << NREGIONS << endl; 
	for (Int_t iSp = 1; iSp < NSPECIES; iSp++)		{
	for (Int_t iType = 0; iType < nType; iType++)	{
	for (Int_t iReg = 0; iReg < 3; iReg++)			{	
		
		//cout << iSp << " " << iType << " " << iReg << " " << Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]) << endl;
		
		tV0massRt[iSp][iType][iReg] = (TNtuple*)mHandler->analysis(0)->dirFile()->Get(Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));
		//cout << "tV0massRt[iSp][iType][iReg] " << tV0massRt[iSp][iType][iReg] << endl;
		//tV0massRt[iSp][iType][iReg] = (TNtuple*)mHandler->analysis(0)->dirFile()->Get(Form("tV0massRt_%s_%s_%s",SPECIES[iSp],TYPE[iType],REGIONS[iReg]));

	}	}	}

}

Bool_t MyAnalysisV0extract::CreateHistograms() {

	hRtV0Yields	= new TH1D("hRtV0Yields",";Species;R_{T} integrated yield",5,-0.5,4.5);

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;
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

		hV0RtFit[iSp][0][0] = new TH1D(Form("hV0RtFit_%s_%s_%s",SPECIES[iSp],TYPE[0],MULTI[3]),
			";R_{T}; Yield",								NRTBINS,RTBINS);
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

	ProducePtSpectraFromHists();
	ProducePtSpectraFromTrees();
	ProduceRtSpectraFromTrees();
	
	return 0;	
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
				iBin,iBin+binSize));		// should be perhaps changed to FindBin
			
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
		cout << "left " << leftrtb << " right " << rightrtb << endl;
		TString treeName = tV0massRt[iSp][iType][iReg]->GetName();
		for (int iBin = 0; iBin < NPTBINS2; ++iBin)	{

			
			tV0massRt[iSp][iType][iReg]->SetName(Form("treept_iSp%i_iBin%i",iSp,iBin));

			yield = ExtractYieldFitPtTree((TTree*)tV0massRt[iSp][iType][iReg]->CopyTree(
				Form("lPt>%f && lPt<%f && lNchTrans>=%f && lNchTrans<%f",XBINS2[iBin],XBINS2[iBin+1],
					leftrtb-1E-05,rightrtb)));

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
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
		}


	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0extract::ProduceRtSpectraFromTrees() {

	mHandler->root()->SetBatch(kTRUE);

	Double_t rt_den = hNchTrans->GetMean();

	iCan = 0;
	nPads = TMath::Nint(TMath::Sqrt(NRTBINS));
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		
		// canvases w rt binning

		cFitsRt[iCan] = new TCanvas(Form("cFitsRt_iSp%i_iType%i_iMu%i",iSp,0,3),"",2800,2000);
		cFitsRt[iCan]->Divide(nPads+1,nPads,0.00005,0.00005);

		cout << " tV0massRt[iSp][0][0] " << tV0massRt[iSp][0][0] << endl;
		TString treeName = tV0massRt[iSp][0][0]->GetName();

		Double_t* yield = 0;
		for (int iBin = 0; iBin < NRTBINS; ++iBin)	{
			
			tV0massRt[iSp][0][0]->SetName(Form("tree_iSp%i_iBin%i",iSp,iBin));
			
			yield = ExtractYieldFitRt((TTree*)tV0massRt[iSp][0][0]->CopyTree(
				Form("lPt>%f && lPt<%f && lNchTrans>%f && lNchTrans<%f",cuts::RT_PTRANGE[0],cuts::RT_PTRANGE[1],
					RTBINS[iBin]*rt_den,RTBINS[iBin+1]*rt_den)));
			
			hV0RtFit[iSp][0][0]->SetBinContent(iBin+1,*(yield+0));
			hV0RtFit[iSp][0][0]->SetBinError(iBin+1,*(yield+1));
		}

		//for entire rt range
		tV0massRt[iSp][0][0]->SetName(Form("tree_iSp%i_iBin%i",iSp,NRTBINS));
			
		yield = ExtractYieldFitRt((TTree*)tV0massRt[iSp][0][0]->CopyTree(
			Form("lPt>%f && lPt<%f && lNchTrans>%f && lNchTrans<%f",cuts::RT_PTRANGE[0],cuts::RT_PTRANGE[1],
				RTBINS[0]*rt_den,RTBINS[NRTBINS]*rt_den)));
		hRtV0Yields->SetBinContent(1+iSp,*(yield+0));
		hRtV0Yields->SetBinError(1+iSp,*(yield+1));
		//

		tV0massRt[iSp][0][0]->SetName(treeName);

		

		cFitsRt[iCan]->Write();
		iCan++;

		hV0RtFit[iSp][0][0]->SetMarkerColor(1);
		hV0RtFit[iSp][0][0]->SetMarkerStyle(20);
		hV0RtFit[iSp][0][0]->SetLineColor(2);

	}

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

Double_t* MyAnalysisV0extract::ExtractYieldFit(TH1D* hist) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	Bool_t empty = (hist->Integral(hist->FindBin(fitMin),hist->FindBin(fitMax)) == 0);

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	hist->Rebin(8);
	RooDataHist DT_set("DT_set","DT_set",MassDT,Import(*hist)); 

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.001,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.1);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
		
	RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
	RooFitResult* fR = 0; 
	if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
		
	RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));

	cFits[canCounter/nBins]->cd(1+canCounter%nBins);
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed));
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();

	TString histName(hist->GetName());
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();
	Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);
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

Double_t* MyAnalysisV0extract::ExtractYieldFitRt(TTree* tree) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	Bool_t empty = (tree->GetEntries(Form("MassDT>%f&&MassDT<%f",fitMin,fitMax)) == 0);

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	RooDataSet DT_set("DT_set","DT_set",MassDT,Import(*tree)); 

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.001,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.1);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
		
	RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
	RooFitResult* fR = 0; 
	if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
		
	RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));

	cFitsRt[canCounterRt/(NRTBINS+1)]->cd(1+canCounterRt%(NRTBINS+1));
	
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed));
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();

	TString histName(tree->GetName());
	cout << histName.Data() << endl;
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();
	Double_t chi2ndf = (!empty) ? plot1->chiSquare() : -1.;

	if (!empty) {
		val[0] = nGaus.getVal();
		//printf("STATUS: int from fit is %f \n", val[0]);
		val[1] = nGaus.getPropagatedError(*fR);
	}

	TLegend* leg1 = new TLegend(0.071,0.57,0.5,0.88);//cFits[canCounter/NPTBINS]->BuildLegend();
	mHandler->MakeNiceLegend(leg1, 0.10, 1.);
	Int_t rtleft	= (binNumber<NRTBINS) ? binNumber: 0;
	Int_t rtright	= (binNumber<NRTBINS) ? binNumber+1 : NRTBINS;
	leg1->AddEntry((TObject*)0,Form("%4.2f < R_{T} < %4.2f",RTBINS[rtleft],RTBINS[rtright])," ");
	leg1->AddEntry((TObject*)0,Form("%s , #chi^{2}/ndf = %4.2f",SPECNAMES[spNumber],chi2ndf)," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	canCounterRt++;


	//if (empty) return val;

	
	
	return val;
}

Double_t* MyAnalysisV0extract::ExtractYieldFitPtTree(TTree* tree) {

	static Double_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	Bool_t empty = (tree->GetEntries(Form("MassDT>%f&&MassDT<%f",fitMin,fitMax)) == 0);

	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	RooDataSet DT_set("DT_set","DT_set",MassDT,Import(*tree)); 

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.001,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e08);
		
	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0.1);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e08);
		
	RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-2,2);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e08);
		
	RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
	RooFitResult* fR = 0; 
	if (!empty) fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
		
	RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));

	//cout << "canCounter/(NPTBINS2) " << canCounter/(NPTBINS2) << " can " << cFitsPtTree[canCounter/(NPTBINS2)] << endl;
	cFitsPtTree[canCounter/(NPTBINS2)]->cd(1+canCounter%(NPTBINS2));
	
	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_set.plotOn(plot1,MarkerSize(0.4));
	if (!empty) fTotal.plotOn(plot1,LineWidth(2),LineColor(kRed));
	plot1->SetMinimum(1e-05);
	plot1->SetMaximum(1.40*plot1->GetMaximum());
	plot1->GetXaxis()->SetTitleSize(0.05);
	plot1->GetYaxis()->SetTitleSize(0.05);
	plot1->Draw();

	TString histName(tree->GetName());
	//cout << histName.Data() << endl;
	TString binName(histName(histName.Index("iBin")+4,2));
	Int_t binNumber = binName.Atoi();
	TString spName(histName(histName.Index("iSp")+3,1));
	Int_t spNumber = spName.Atoi();
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

