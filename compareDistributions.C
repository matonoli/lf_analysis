
#include <iostream>

void compareDistributions(const Char_t *inputFile="test.list", const Char_t *MCinputFile="testMC.list") {

	gROOT->LoadMacro("MyKit/MyHandler.cxx+");

	MyHandler* handlerD 	= new MyHandler();
	handlerD->LoadInputTree(inputFile,"PIDTree");
	TChain* treeD 			= handlerD->chain();  

	MyHandler* handlerMC 	= new MyHandler();
	handlerMC->LoadInputTree(MCinputFile,"PIDTree");
	TChain* treeMC 			= handlerMC->chain();

	TH1::SetDefaultSumw2();
	const Int_t NDISTRO = 14;
	const Int_t NBINS = 200;
	const char* DISTRO[NDISTRO] = {"AnalysisV0Track.fRadius",
		"AnalysisV0Track.fV0CosinePA",
		"AnalysisV0Track.fDCAV0Daughters",
		"AnalysisV0Track.fDCAPV",

		"AnalysisV0Track.fEta",
		"AnalysisV0Track.fPt",
		"AnalysisV0Track.fPosAnalysisPIDTrack.fEta",
		"AnalysisV0Track.fPosAnalysisPIDTrack.fPt",

		"AnalysisV0Track.fPosAnalysisPIDTrack.fPhi",
		"AnalysisV0Track.fPosAnalysisPIDTrack.fImpactParameter",
		"AnalysisV0Track.fPosAnalysisPIDTrack.fTPCNcls",
		"AnalysisV0Track.fPosAnalysisPIDTrack.fTPCNclsF",
		
		"AnalysisV0Track.fPosAnalysisPIDTrack.fTPCNcr",
		"AnalysisV0Track.fPosAnalysisPIDTrack.nSigmaPionTPC"};

	const Double_t LOW[NDISTRO] = {-10.,	0.95,	-0.1,	-0.1,
									-1.1,	0.,		-1.1,	0.,	
									-0.1,	-0.1,	-1.,	-1.,
									-1.,	-8.};
	const Double_t HIGH[NDISTRO] = {110.,	1.01,	1.5,	1.5,
									1.1,	16.,	1.1,	16.,
									6.4,	5.,		200.,	200.,
									200.,	8.};

	TCanvas* c1 = new TCanvas("c1","",1800,900);
	c1->Divide(5,4,1e-4,1e-4);
	for (int iD = 0; iD < NDISTRO; ++iD)	{
		
		c1->cd(1+iD);
		c1->GetPad(1+iD)->SetLogy();
		c1->GetPad(1+iD)->SetGridx();
		TH1D* hD 	= new TH1D(Form("hD%i",iD),"",NBINS,LOW[iD],HIGH[iD]);
		TH1D* hMC 	= new TH1D(Form("hMC%i",iD),"",NBINS,LOW[iD],HIGH[iD]);
		TH1D* hDl 	= new TH1D(Form("hDl%i",iD),"",NBINS,LOW[iD],HIGH[iD]);
		TH1D* hMCl 	= new TH1D(Form("hMCl%i",iD),"",NBINS,LOW[iD],HIGH[iD]);
		treeD->Draw(Form("%s>>hD%i",DISTRO[iD],iD),"TMath::Abs(AnalysisV0Track.fInvMK0s)<0.002","");
		treeMC->Draw(Form("%s>>hMC%i",DISTRO[iD],iD),"TMath::Abs(AnalysisV0Track.fInvMK0s)<0.002","");
		treeD->Draw(Form("%s>>hDl%i",DISTRO[iD],iD),"TMath::Abs(AnalysisV0Track.fInvMK0s)<0.002 && AnalysisV0Track.fPt<0.5","");
		treeMC->Draw(Form("%s>>hMCl%i",DISTRO[iD],iD),"TMath::Abs(AnalysisV0Track.fInvMK0s)<0.002 && AnalysisV0Track.fPt<0.5","");
		handlerD->MakeNiceHistogram(hD,kBlack);
		handlerD->MakeNiceHistogram(hMC,kRed);
		handlerD->MakeNiceHistogram(hDl,kBlack);
		handlerD->MakeNiceHistogram(hMCl,kRed);
		hD->SetMarkerSize(0.3);
		hMC->SetMarkerSize(0.3);
		hDl->SetMarkerSize(0.3);
		hMCl->SetMarkerSize(0.3);

		hDl->SetMarkerStyle(27);
		hMCl->SetMarkerStyle(27);
		hD->Scale(1./hD->Integral());
		hMC->Scale(1./hMC->Integral());
		hDl->Scale(1./hDl->Integral());
		hMCl->Scale(1./hMCl->Integral());
		Double_t max = (hD->GetMaximum()>hMC->GetMaximum()) ? hD->GetMaximum() : hMC->GetMaximum();
		Double_t min = (hD->GetMinimum()<hMC->GetMinimum()) ? hD->GetMinimum() : hMC->GetMinimum();
		hD->GetYaxis()->SetRangeUser(-0.1*max,1.1*max);
		hD->GetYaxis()->SetRangeUser(0.0001,2.0);
		c1->SetGridx();
		c1->SetLogy();
		hD->Draw("hist p");
		hMC->Draw("hist p same");
		hDl->Draw("hist p same");
		hMCl->Draw("hist p same");

		TLegend *leg1 = new TLegend(0.25,0.60,0.75,0.85);
		handlerD->MakeNiceLegend(leg1,0.037,1);
		leg1->AddEntry((TObject*)0,Form("%s",DISTRO[iD]),"");
		leg1->AddEntry(hD,"Data","pl");
		leg1->AddEntry(hMC,"MC","pl");
		leg1->AddEntry(hDl,"<0.5 GeV","p");
		leg1->Draw();

		cout << "Finished drawing iD: " << iD << endl;
	}

	c1->SaveAs(Form("comparison/%s_%s.pdf",inputFile,MCinputFile));
	c1->SaveAs(Form("comparison/%s_%s.png",inputFile,MCinputFile));

}