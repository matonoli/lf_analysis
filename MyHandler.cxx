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

ClassImp(MyHandler)

void MyHandler::AddAnalysis(MyAnalysis* ana) {

	if (!mChain) {
		printf("Chain not yet created!\n"); }
		//return; 	}
	
	this->mAnalysis = ana;
	mAnalysis->SetHandler(this);
}

Int_t MyHandler::Init() {

	mAnalysis->Init();

	return 0;
}

Int_t MyHandler::LoadInput(const Char_t *inputFile, const Char_t *chainName) {
	
	TString inputFileStr(inputFile);
	mFlagMC = inputFileStr.Contains("MC");
	mFlagHist = inputFileStr.BeginsWith("hist");

	if(!mFlagHist) {
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
		return 1;
		}

		else if (dirFile.find(".root") != string::npos)	{
			mChain->Add(dirFile.c_str());	
			return 1;	}

		else	{
			cout << " No good input file to read ... " << endl;
			return 1;	}

	} else {
		mFile = new TFile(inputFileStr.Data(),"READ");	}
	
	return 1;
}


Int_t MyHandler::Make(Int_t iEv) {

	Int_t iret = 0;
	//printf("Looping in handler %i \n", iEv);
	
	if (!mFlagHist) {
		if (!mChain->GetEntry(iEv)) iret++;
		iret += mAnalysis->Make(iEv);
	}

	return iret;	
}

Int_t MyHandler::Finish() {
	
	printf("Finishing handler \n");
	mAnalysis->Finish();
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
