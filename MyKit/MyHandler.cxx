#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TH1D.h>
#include <TLegend.h>

#include "MyHandler.h"
#include "MyAnalysis.h"

#include <iostream>
#include <fstream>
#include <string>

#include <AliAnalysisPIDCascadeEvent.h>
#include <AliESDEvent.h>
#include <AliESDtrackCuts.h>

ClassImp(MyHandler)

MyHandler::MyHandler() : mOutName(""), mFile(0), mDir(0), mChain(0) {

	mFlagMC 	= false;
	mFlagHist	= false;
	mRebinPt	= false;
	
	for (int i = 0; i < 5; ++i)
	{
		mAnalysis[i] = 0;
	}
}

void MyHandler::AddAnalysis(MyAnalysis* ana) {

	if (!mChain) {
		printf("Chain not yet created!\n"); }
		//return; 	}
	
	this->mAnalysis[nAnalyses] = ana;
	mAnalysis[nAnalyses]->SetHandler(this);
	nAnalyses++;
}

Int_t MyHandler::Init() {

	if (mOutName.Length()) {
		printf("Creating mother output file %s \n", mOutName.Data());
		mFile = new TFile(mOutName.Data(),"RECREATE");
	}

	printf("nAnalyses is %i \n", nAnalyses);
	for (Int_t iAna = 0; iAna < nAnalyses; iAna++) {
		nAna = iAna;
		mAnalysis[iAna]->Init();
		if (mFile) mFile->cd();
	}
	return 0;
}

Int_t MyHandler::LoadInputTree(const Char_t *inputFile, const Char_t *chainName) {	// this should be only called when running on local trees
	
	TString inputFileStr(inputFile);
	mFlagMC = inputFileStr.Contains("MC");

	mChain = new TChain(chainName);
	std::string const dirFile = inputFileStr.Data();
	if (dirFile.find(".lis") != std::string::npos)	{			
		std::ifstream inputStream(dirFile.c_str());
		if (!inputStream)	{
			std::cout << "ERROR: Cannot open list file " << dirFile << std::endl;
			return 0;	}
		int nFile = 0;
		std::string file;

		while (getline(inputStream, file))	{
			if (file.find(".root") != std::string::npos)	{
				TFile* ftmp = TFile::Open(file.c_str());
				if (ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())	{
					std::cout << " Read in file " << file << std::endl;
					mChain->Add(file.c_str());
					++nFile;	}
				if (ftmp) ftmp->Close();
			}
		}

		std::cout << " Total " << nFile << " files have been read in. " << std::endl;
	}
	else if (dirFile.find(".root") != std::string::npos)	{
		mChain->Add(dirFile.c_str());	}
	else	{
		std::cout << " No good input file to read ... " << std::endl;
		return 1;	}

	mChain->SetBranchAddress("AnalysisEvent",&mEvent);
	mChain->SetBranchAddress("AnalysisTrack",&bTracks);
	mChain->SetBranchAddress("AnalysisV0Track",&bV0s);
	if (mFlagMC) mChain->SetBranchAddress("AnalysisParticle",&bParticles);
}

Int_t MyHandler::LoadInputHist(const Char_t *inputFile) {
	
	TString inputFileStr(inputFile);
	mFlagMC = inputFileStr.Contains("MC");
	mFlagHist = 1;

	mFileHist = new TFile(inputFileStr.Data(),"READ");	
	return 1;
}


Int_t MyHandler::Make(Int_t iEv) {

	Int_t iret = 0;
	//printf("Looping in handler %i \n", iEv);
	
	if (!mFlagHist && mChain) {
		if (!mChain->GetEntry(iEv)) iret++;
		//printf("iret is %i \n", iret); 
	}

	for (Int_t iAna = 0; iAna < nAnalyses; iAna++) {
		iret += mAnalysis[iAna]->Make(iEv); 
		//printf("iret is %i \n", iret);
	}

	return iret;	
}

Int_t MyHandler::Finish() {
	
	printf("Finishing handler \n");
	for (Int_t iAna = 0; iAna < nAnalyses; iAna++) {
		mAnalysis[iAna]->Finish();
	}

	if (mFile) {
		mFile->cd();
		mFile->Write();
		printf("File written \n");
	}

	return 0;	
}

void MyHandler::DrawCut(Double_t cut, Int_t direction, TVirtualPad* can) {

	Double_t x[2] = {cut, cut};
	Double_t y[2];
	y[0] = can->GetUymin(); y[1] = can->GetUymax();

	TGraph* gcut = new TGraph(2, x, y);
	if (direction==2) gcut->SetLineWidth(-202);
	if (direction==1) gcut->SetLineWidth(202);
	if (direction==0) gcut->SetLineWidth(2);
	gcut->SetLineColor(kRed);
	gcut->SetFillColor(kRed);
	gcut->SetFillStyle(3003);

	gcut->Draw("same");
}

void MyHandler::MakeNiceHistogram(TH1D* h, Int_t col) {

	h->SetLineColor(col);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(1.3);
	h->SetMarkerColor(col);
	h->SetStats(0);

	h->GetYaxis()->SetTitleOffset(1.1);
	h->GetYaxis()->SetLabelOffset(0.0025);
	h->GetYaxis()->SetLabelSize(0.03);

	//h->SetTopMargin(0.055);
}

void MyHandler::MakeNiceLegend(TLegend *leg, Float_t size, Int_t columns)	{
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetMargin(0.25);
	leg->SetTextSize(size);
	leg->SetEntrySeparation(0.5);
	leg->SetNColumns(columns);

}

void MyHandler::MakeRatioPlot(TH1D* hn, TH1D* hd, TCanvas* c, Double_t low, Double_t high, Double_t lowx, Double_t highx) {
	
	c->cd();

	// check for an already existent ratio plot
	Bool_t hasRatio = false;
	TObject* obj;
	TIter next(c->GetListOfPrimitives());
	while ( (obj = next()) ) {
		TString objName = obj->GetName();
		if (objName == Form("p2_%s",c->GetName())) {
			TVirtualPad* prat = (TVirtualPad*)obj;
			prat->cd();
			hasRatio = true;
		}
	}

	if (!hasRatio) {

		TCanvas* ctop = (TCanvas*)c->Clone("ctop");
		c->Clear();
		ctop->SetBottomMargin(0.005);
		c->cd();

		TPad* p1 = new TPad(Form("p1_%s",c->GetName()),"",0.,0.3,1.,1.);
		p1->SetBottomMargin(0.);
		p1->Draw();
		p1->cd();
		ctop->DrawClonePad();

		c->cd();
		TPad* p2 = new TPad(Form("p2_%s",c->GetName()),"",0.,0.00,1.,0.28);
		p2->SetTopMargin(0);
		p2->SetBottomMargin(0.32);
		p2->Draw();
		p2->cd();
	}

	TH1D* hr = (TH1D*)hn->Clone(Form("hr_%s",hn->GetName()));
	hr->SetMinimum(low);
	hr->SetMaximum(high);
	hr->GetXaxis()->SetRangeUser(lowx,highx);
	hr->Divide(hd);

	hr->GetYaxis()->SetTitle("ratio");
	hr->GetYaxis()->CenterTitle();
	hr->GetYaxis()->SetNdivisions(505);
	hr->GetYaxis()->SetTitleSize(25);
	//hr->GetYaxis()->SetTitleFont(43);
	hr->GetYaxis()->SetTitleOffset(1.55);
	hr->GetYaxis()->SetLabelFont(43); 
	hr->GetYaxis()->SetLabelSize(20);

	hr->GetXaxis()->SetTitleSize(25);
	hr->GetXaxis()->SetTitleFont(43);
	hr->GetXaxis()->SetTitleOffset(4.);
	hr->GetXaxis()->SetLabelFont(43); 
	hr->GetXaxis()->SetLabelSize(25);
	hr->GetXaxis()->SetTickLength(0.09);

	if (!hasRatio)	hr->Draw();
	else			hr->Draw("same");

	//c->SetCanvasSize()
	c->cd();

}


void MyHandler::MakeZoomPlot(TH1D* h, TCanvas* c, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax) {

	c->cd();
	TPad* p = new TPad(Form("p_%s",c->GetName()),"",0.13,0.5,0.5,0.89);
	p->SetLogy(0);
	p->Draw();
	p->cd();

	//Double_t prevXmin = h->GetXaxis()->GetR
	h->GetXaxis()->SetRangeUser(xmin,xmax);
	h->GetYaxis()->SetRangeUser(ymin,ymax);
	//

}

Double_t MyHandler::DeltaPhi(Double_t phi1, Double_t phi2) {
	
	Double_t dphi = phi2 - phi1;
	if ( dphi > TMath::Pi() )		dphi = dphi - 2*TMath::Pi();
	if ( dphi < -1.*TMath::Pi() )	dphi = dphi + 2*TMath::Pi();

	return dphi;
}
	

#if INPUTFORMAT == 2
void MyHandler::SetupTrackCuts() {
	
	mTrackCuts2010 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
	mTrackCuts2010->SetEtaRange(-0.8,0.8);

	mTrackCuts2011 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
	mTrackCuts2011->SetEtaRange(-0.8,0.8);

	mTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	mTrackCutsTPCOnly->SetEtaRange(-0.8,0.8);
	mTrackCutsTPCOnly->SetRequireTPCRefit(kTRUE);

	mTrackCuts2011sys =  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
	mTrackCuts2011sys->SetEtaRange(-0.8,0.8);
	mTrackCuts2011sys->SetMinNCrossedRowsTPC(60);
	mTrackCuts2011sys->SetMaxChi2PerClusterTPC(5);
	mTrackCuts2011sys->SetMaxDCAToVertexZ(3);

	mTrackCutsV0d = AliESDtrackCuts::GetStandardV0DaughterCuts();
	mTrackCutsV0d->SetEtaRange(-0.8,0.8);
}
#endif
