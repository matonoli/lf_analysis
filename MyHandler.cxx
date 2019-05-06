#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>

#include "MyHandler.h"
#include "MyAnalysis.h"

#include <iostream>
#include <fstream>

#include <AliAnalysisPIDEvent.h>

ClassImp(MyHandler)

void MyHandler::AddAnalysis(MyAnalysis* ana) {

	if (!mChain) {
		printf("Chain not yet created!\n"); }
		//return; 	}
	
	this->mAnalysis[nAnalyses] = ana;
	mAnalysis[nAnalyses]->SetHandler(this);
	nAnalyses++;
}

Int_t MyHandler::Init() {

	printf("Creating mother output file %s \n", mOutName.Data());
	mFile = new TFile(mOutName.Data(),"RECREATE");

	printf("nAnalyses is %i \n", nAnalyses);
	for (Int_t iAna = 0; iAna < nAnalyses; iAna++) {
		nAna = iAna;
		mAnalysis[iAna]->Init();
		mFile->cd();
	}
	return 0;
}

Int_t MyHandler::LoadInputTree(const Char_t *inputFile, const Char_t *chainName) {
	
	TString inputFileStr(inputFile);
	mFlagMC = inputFileStr.Contains("MC");

	mChain = new TChain(chainName);
	string const dirFile = inputFileStr.Data();
	if (dirFile.find(".lis") != string::npos)	{			
		ifstream inputStream(dirFile.c_str());
		if (!inputStream)	{
			cout << "ERROR: Cannot open list file " << dirFile << endl;
			return 0;	}
		int nFile = 0;
		string file;

		while (getline(inputStream, file))	{
			if (file.find(".root") != string::npos)	{
				TFile* ftmp = TFile::Open(file.c_str());
				if (ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())	{
					cout << " Read in file " << file << endl;
					mChain->Add(file.c_str());
					++nFile;	}
				if (ftmp) ftmp->Close();
			}
		}

		cout << " Total " << nFile << " files have been read in. " << endl;
	}
	else if (dirFile.find(".root") != string::npos)	{
		mChain->Add(dirFile.c_str());	}
	else	{
		cout << " No good input file to read ... " << endl;
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
	
	if (!mFlagHist) {
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

	mFile->cd();
	mFile->Write();
	printf("File written \n");

	return 0;	
}

void MyHandler::DrawCut(Double_t cut, Int_t direction, TCanvas* can) {

	Double_t x[2] = {cut, cut};
	Double_t y[2];
	y[0] = can->GetUymin(); y[1] = can->GetUymax();

	TGraph* gcut = new TGraph(2, x, y);
	if (direction) gcut->SetLineWidth(-402);
	else gcut->SetLineWidth(402);
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
