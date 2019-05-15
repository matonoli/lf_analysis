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
	BorrowHistograms();
	CreateHistograms();


	mList = (TList*)mHandler->directory()->GetList();

	return 0;
}

Int_t MyAnalysisV0extract::Make(Int_t iEv) {
	return 0;
}


Bool_t MyAnalysisV0extract::BorrowHistograms() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < 1; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		hV0IMvPt[iSp][iType][iMu][iSph] 
			= (TH2D*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPt_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]));

	} } } }

}

Bool_t MyAnalysisV0extract::CreateHistograms() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < 1; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
			
		hV0PtFit[iSp][iType][iMu][iSph]		= new TH1D(Form("hV0PtFit_%s_%s_%s_%s",SPECIES[iSp],TYPE[iType],MULTI[iMu],SPHERO[iSph]),
			";V0 Pt (GeV/#it{c}); Entries",								NPTBINS,XBINS);
		
	} } } }

}


Int_t MyAnalysisV0extract::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	mHandler->root()->SetBatch(kTRUE);

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iType = 0; iType < 1; ++iType)		{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{

		cFits[iCan] = new TCanvas(Form("cFits_iSp%i_iMu%i_iSph%i",iSp,iMu,iSph),"",2800,2000);
		cFits[iCan]->Divide(8,7,0.00005,0.00005);
			
		Double_t* yield = 0;
		for (int iBin = 1; iBin < NPTBINS+1; ++iBin)	{
		//for (int iBin = 10; iBin < 11; ++iBin)	{
			yield = ExtractYieldFit(hV0IMvPt[iSp][iType][iMu][iSph]->ProjectionY(
				Form("iSp%i_iType%i_iMu%i_iSph%i_iBin%i", iSp, iType, iMu, iSph, iBin),iBin,iBin));
			hV0PtFit[iSp][iType][iMu][iSph]->SetBinContent(iBin,*(yield+0));
			hV0PtFit[iSp][iType][iMu][iSph]->SetBinError(iBin,*(yield+1));
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
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") ) {
			mDirFile->Remove(obj);	}
		}


	mHandler->root()->SetBatch(kFALSE);
	return 0;	
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

	cFits[canCounter/NPTBINS]->cd(1+canCounter%NPTBINS);
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
	leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",XBINS[binNumber-1],XBINS[binNumber])," ");
	leg1->AddEntry((TObject*)0,Form("%s , #chi^{2}/ndf = %4.2f",SPECNAMES[spNumber],chi2ndf)," ");
	leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",val[0],val[1])," ");
	leg1->Draw();
	
	canCounter++;


	//if (empty) return val;

	
	
	return val;
}

