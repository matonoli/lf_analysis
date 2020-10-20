{
	TFile* fileB = new TFile("testeta_b.root","READ");
	TFile* fileE = new TFile("testeta_esd.root","READ");

	TFile* fileB0 = new TFile("../rootOutputs/beast_200911_hist.root","READ");
	TFile* fileE0 = new TFile("../rootOutputs/pp16kP2_200921_AnalysisResults_hist.root","READ");

	TFile* fileBMC = new TFile("../rootOutputs/beastMC_200911_hist.root","READ");
	TFile* fileEMC = new TFile("../rootOutputs/pp16kP2MC_200921_AnalysisResults_hist.root","READ");


	// DATA FROM BEAST
	fileB0->cd();
	TDirectoryFile* MyAnalysisV0_0 = (TDirectoryFile*)fileB0->Get("MyAnalysisV0_0");
	MyAnalysisV0_0->cd();

	Double_t nEv_B = ((TH1D*)MyAnalysisV0_0->Get("hEventType"))->GetBinContent(24);

	// PT SPECTRA
	TH1D* hPtD_B	= (TH1D*)MyAnalysisV0_0->Get("hV0Pt_K0s_D_MB_MB");
	// ETA, Y
	TH1D* hEtaD_B 	= (TH1D*)MyAnalysisV0_0->Get("hV0Eta_K0s_D_MB_MB");
	TH1D* hYD_B 	= (TH1D*)MyAnalysisV0_0->Get("hV0Y_K0s_D_MB_MB");
	// SIGNAL SHAPE
	TH1D* hSigD_B	= (TH1D*)((TH2D*)MyAnalysisV0_0->Get("hV0IMvPt_K0s_D_MB_MB"))->ProjectionY("hSigD_B",3,3);

	// OTHER STUFF FROM BEAST
	fileB->cd();
	
	// YIELDS
	TDirectoryFile* MyAnalysisV0extract_1 = (TDirectoryFile*)fileB->Get("MyAnalysisV0extract_1");
	MyAnalysisV0extract_1->cd();
	TH1D* hRawD_B 	= (TH1D*)MyAnalysisV0extract_1->Get("hV0PtFit_K0s_D_MB_MB");

	// CORR. SPECTRA
	TDirectoryFile* MyAnalysisV0correct_2 = (TDirectoryFile*)fileB->Get("MyAnalysisV0correct_2");
	MyAnalysisV0correct_2->cd();
	TH1D* hCorrD_B 	= (TH1D*)MyAnalysisV0correct_2->Get("hV0PtFitCorr_K0s_D_MB_MB");	

	// EFFICIENCY
	fileBMC->cd();
	MyAnalysisV0_0 = (TDirectoryFile*)fileBMC->Get("MyAnalysisV0_0");
	MyAnalysisV0_0->cd();
	TH1D* hEffi_B 	= (TH1D*)MyAnalysisV0_0->Get("hV0Efficiency_K0s");

	// RC SPECTRA
	TH1D* hPtRC_B	= (TH1D*)MyAnalysisV0_0->Get("hV0Pt_K0s_RC_MB_MB");
	TH1D* hEtaRC_B	= (TH1D*)MyAnalysisV0_0->Get("hV0Eta_K0s_RC_MB_MB");
	Double_t nEv_BMC = ((TH1D*)MyAnalysisV0_0->Get("hEventType"))->GetBinContent(24);	

	// DATA FROM ESD
	fileE0->cd();
	MyAnalysisV0_0 = (TDirectoryFile*)fileE0->Get("MyAnalysisV0_0");
	MyAnalysisV0_0->cd();

	Double_t nEv_E = ((TH1D*)MyAnalysisV0_0->Get("hEventType"))->GetBinContent(24);

	// PT SPECTRA
	TH1D* hPtD_E	= (TH1D*)MyAnalysisV0_0->Get("hV0Pt_K0s_D_MB_MB");
	// ETA, Y
	TH1D* hEtaD_E 	= (TH1D*)MyAnalysisV0_0->Get("hV0Eta_K0s_D_MB_MB");
	TH1D* hYD_E 	= (TH1D*)MyAnalysisV0_0->Get("hV0Y_K0s_D_MB_MB");
	// SIGNAL SHAPE
	TH1D* hSigD_E	= (TH1D*)((TH2D*)MyAnalysisV0_0->Get("hV0IMvPt_K0s_D_MB_MB"))->ProjectionY("hSigD_E",3,3);

	// OTHER STUFF FROM BEAST
	fileE->cd();
	
	// YIELDS
	MyAnalysisV0extract_1 = (TDirectoryFile*)fileE->Get("MyAnalysisV0extract_1");
	MyAnalysisV0extract_1->cd();
	TH1D* hRawD_E 	= (TH1D*)MyAnalysisV0extract_1->Get("hV0PtFit_K0s_D_MB_MB");

	// CORR. SPECTRA
	MyAnalysisV0correct_2 = (TDirectoryFile*)fileE->Get("MyAnalysisV0correct_2");
	MyAnalysisV0correct_2->cd();
	TH1D* hCorrD_E 	= (TH1D*)MyAnalysisV0correct_2->Get("hV0PtFitCorr_K0s_D_MB_MB");	

	// EFFICIENCY
	fileEMC->cd();
	MyAnalysisV0_0 = (TDirectoryFile*)fileEMC->Get("MyAnalysisV0_0");
	MyAnalysisV0_0->cd();
	TH1D* hEffi_E 	= (TH1D*)MyAnalysisV0_0->Get("hV0Efficiency_K0s");

	// RC SPECTRA
	TH1D* hPtRC_E	= (TH1D*)MyAnalysisV0_0->Get("hV0Pt_K0s_RC_MB_MB");
	TH1D* hEtaRC_E	= (TH1D*)MyAnalysisV0_0->Get("hV0Eta_K0s_RC_MB_MB");
	Double_t nEv_EMC = ((TH1D*)MyAnalysisV0_0->Get("hEventType"))->GetBinContent(24);	


	////// RATIO PLOTS
	hPtD_B->Divide(hPtD_E); hPtD_B->Scale(nEv_E/nEv_B);
	hEtaD_B->Divide(hEtaD_E); hEtaD_B->Scale(nEv_E/nEv_B);
	hYD_B->Divide(hYD_E); hYD_B->Scale(nEv_E/nEv_B);
	hRawD_B->Divide(hRawD_E); hRawD_B->Scale(nEv_E/nEv_B);
	hCorrD_B->Divide(hCorrD_E);
	hEffi_B->Divide(hEffi_E);
	hPtRC_B->Divide(hPtRC_E); hPtRC_B->Scale(nEv_EMC/nEv_BMC);
	hEtaRC_B->Rebin(4);hEtaRC_E->Rebin(4);
	hEtaRC_B->Divide(hEtaRC_E); hEtaRC_B->Scale(nEv_EMC/nEv_BMC);
	hSigD_B->Rebin(20); hSigD_E->Rebin(20);
	//hSigD_B->Divide(hSigD_E); hSigD_B->Scale(nEv_E/nEv_B);

	

	TCanvas* c1 = new TCanvas("c1","",1200,900);
	c1->Divide(3,3,0.001,0.001);

	c1->cd(1);
	hPtD_B->SetMarkerStyle(20); hPtD_B->SetMarkerColor(kRed); 
	hPtD_B->Draw();
	c1->cd(2);
	hEtaD_B->SetMarkerStyle(20); hEtaD_B->SetMarkerColor(kRed); 
	hEtaD_B->Draw();
	c1->cd(3);
	hYD_B->SetMarkerStyle(20); hYD_B->SetMarkerColor(kRed); 
	hYD_B->Draw();
	c1->cd(4);
	hSigD_B->SetMarkerStyle(20); hSigD_B->SetMarkerColor(kRed); 
	hSigD_B->Draw();
	c1->cd(5);
	hRawD_B->SetMarkerStyle(20); hRawD_B->SetMarkerColor(kRed); 
	hRawD_B->Draw();
	c1->cd(6);
	hCorrD_B->SetMarkerStyle(20); hCorrD_B->SetMarkerColor(kRed); 
	hCorrD_B->Draw();
	c1->cd(7);
	hEffi_B->SetMarkerStyle(20); hEffi_B->SetMarkerColor(kRed);
	hEffi_B->GetYaxis()->SetRangeUser(0.75,1.25);
	hEffi_B->Draw();
	c1->cd(8);
	hPtRC_B->SetMarkerStyle(20); hPtRC_B->SetMarkerColor(kRed); 
	hPtRC_B->Draw();
	c1->cd(9);
	hEtaRC_B->SetMarkerStyle(20); hEtaRC_B->SetMarkerColor(kRed); 
	hEtaRC_B->Draw();
	




	// 

}