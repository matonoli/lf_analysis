void MakeNiceHistogram(TH1D* h, Int_t col) {

	h->SetLineColor(col);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(1.3);
	h->SetMarkerColor(col);
	h->SetStats(0);

	h->GetYaxis()->SetTitleOffset(1.1);
	h->GetYaxis()->SetLabelOffset(0.0025);
	h->GetYaxis()->SetLabelSize(0.03);

}

void MakeNiceLegend(TLegend *leg, Float_t size, Int_t columns)	{
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetMargin(0.25);
	leg->SetTextSize(size);
	leg->SetEntrySeparation(0.5);
	leg->SetNColumns(columns);

}

void MakeRatioPlot(TH1D* hn, TH1D* hd, TCanvas* c, Double_t low, Double_t high, Double_t lowx, Double_t highx) {
	
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

void compare18v16() {
	TFile* fD[12];
	TFile* fMC[12];
	TFile* fO[12];
	TDirectoryFile* dfD[12];
	TDirectoryFile* dfDx[12];
	TDirectoryFile* dfDc[12];
	TDirectoryFile* dfMC[12];
	TString name[] = {"b","d","e","f","g","l","m","n","o","p","all","16k"};
	Int_t col[] = {kAzure+2, kBlue+2, kViolet+2, kMagenta+2,
					kPink+2, kRed+2, kOrange+2, kCyan+2, kGreen+2,
					kTeal+2, kBlack, kRed};

	Double_t nEvD[12];
	Double_t nEvMC[12];
	TH1D* hRaw[12];
	TH1D* hEffi[12];
	TH1D* hMB[12];
	TH1D* hV0M[12];

	for (int i = 0; i < 11; ++i)	{
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		fD[i] 		= new TFile(Form("../rootOutputs/201020_2018x/beast_201020_2018_LHC18%s_hist.root",name[i].Data()),"READ");
		fMC[i] 		= new TFile(Form("../rootOutputs/201020_2018x/beastMC_201020_2018_LHC18%s_hist.root",name[i].Data()),"READ");
		fO[i] 		= new TFile(Form("18%s.root",name[i].Data()),"READ");
		dfD[i] 		= (TDirectoryFile*)fD[i]->Get("MyAnalysisV0_0");
		dfMC[i] 	= (TDirectoryFile*)fMC[i]->Get("MyAnalysisV0_0");
		dfDx[i] 	= (TDirectoryFile*)fO[i]->Get("MyAnalysisV0extract_1");
		dfDc[i] 	= (TDirectoryFile*)fO[i]->Get("MyAnalysisV0correct_2");

		//cout << "doing i " << i << " " << nEvD[i] << endl;
		nEvD[i] 	= ((TH1D*)dfD[i]->Get("hEventType"))->GetBinContent(24);
		nEvMC[i] 	= ((TH1D*)dfMC[i]->Get("hEventType"))->GetBinContent(24);

		hRaw[i]		= (TH1D*)dfDx[i]->Get("hV0PtFit_K0s_D_MB_MB");
		hEffi[i]	= (TH1D*)dfMC[i]->Get("hV0Efficiency_K0s");
		hMB[i]		= (TH1D*)dfDc[i]->Get("hV0PtFitCorr_K0s_D_MB_MB");
		hV0M[i]		= (TH1D*)dfDc[i]->Get("hV0PtFitCorr_K0s_D_V0M_MB");
	}

		//fD[11] 		= new TFile("../rootOutputs/beast_200911_hist.root","READ");
		//fMC[11] 	= new TFile("../rootOutputs/beastMC_200911_hist.root","READ");
		fD[11] 		= new TFile("../rootOutputs/beast_201017_2016_LHC16k_hist.root","READ");
		fMC[11] 	= new TFile("../rootOutputs/201020_2018x/beastMC_201020_2016_LHC16k_hist.root","READ");
		fO[11] 		= new TFile("16k.root","READ");
		dfD[11] 		= (TDirectoryFile*)fD[11]->Get("MyAnalysisV0_0");
		dfMC[11] 	= (TDirectoryFile*)fMC[11]->Get("MyAnalysisV0_0");
		dfDx[11] 	= (TDirectoryFile*)fO[11]->Get("MyAnalysisV0extract_1");
		dfDc[11] 	= (TDirectoryFile*)fO[11]->Get("MyAnalysisV0correct_2");

		nEvD[11] 	= ((TH1D*)dfD[11]->Get("hEventType"))->GetBinContent(24);
		cout << "doing i " << 11 << " " << nEvD[11] << endl;
		nEvMC[11] 	= ((TH1D*)dfMC[11]->Get("hEventType"))->GetBinContent(24);

		hRaw[11]	= (TH1D*)dfDx[11]->Get("hV0PtFit_K0s_D_MB_MB");
		hEffi[11]	= (TH1D*)dfMC[11]->Get("hV0Efficiency_K0s");
		hMB[11]		= (TH1D*)dfDc[11]->Get("hV0PtFitCorr_K0s_D_MB_MB");
		hV0M[11]	= (TH1D*)dfDc[11]->Get("hV0PtFitCorr_K0s_D_V0M_MB");

	////


	cout << "NUMBER OF EVENTS (data):" << endl;
	for (int i = 0; i < 12; ++i)	{
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		printf("18%s    : %3.1f\n",name[i].Data(),nEvD[i]);
	}
	cout << "NUMBER OF EVENTS (MC):" << endl;
	for (int i = 0; i < 12; ++i)	{
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		printf("18%s    : %3.1f\n",name[i].Data(),nEvMC[i]);
	}

	cout << "Drawing raw yields" << endl;
	TCanvas* cRaw = new TCanvas("cRaw","cRaw",1200,1200);
	TLegend* legRaw = new TLegend(0.48,0.64,0.81,0.85);
	MakeNiceLegend(legRaw,0.035,2);
	cRaw->SetLogy();
	MakeNiceHistogram(hRaw[10],col[10]); hRaw[10]->Scale(1./nEvD[10]);
	legRaw->AddEntry(hRaw[10],Form("LHC18%s",name[10].Data()),"pl");
	hRaw[10]->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeNiceHistogram(hRaw[i],col[i]); hRaw[i]->Scale(1./nEvD[i]);
		legRaw->AddEntry(hRaw[i],Form("LHC18%s",name[i].Data()),"pl");
		hRaw[i]->Draw("same");
	}
	legRaw->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeRatioPlot(hRaw[i],hRaw[10],cRaw,0.9,1.1,0.,14.);
	}
	cRaw->Draw();
	cRaw->SaveAs("cRaw.png");

	/////////////////

	cout << "Drawing efficiencies" << endl;
	TCanvas* cEffi = new TCanvas("cEffi","cEffi",1200,1200);
	TLegend* legEffi = new TLegend(0.48,0.64,0.81,0.85);
	MakeNiceLegend(legEffi,0.035,2);
	MakeNiceHistogram(hEffi[10],col[10]); hEffi[10]->Scale(1./4);
	legEffi->AddEntry(hEffi[10],Form("LHC18%s",name[10].Data()),"pl");
	hEffi[10]->GetYaxis()->SetRangeUser(0.0,1.0);
	hEffi[10]->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeNiceHistogram(hEffi[i],col[i]); 
		legEffi->AddEntry(hEffi[i],Form("LHC18%s",name[i].Data()),"pl");
		hEffi[i]->Draw("same");
	}
	legEffi->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeRatioPlot(hEffi[i],hEffi[10],cEffi,0.9,1.1,0.,14.);
	}
	cEffi->Draw();
	cEffi->SaveAs("cEffi.png");

	//////////////

	cout << "Drawing MB" << endl;
	TCanvas* cMB = new TCanvas("cMB","cMB",1200,1200);
	TLegend* legMB = new TLegend(0.48,0.64,0.81,0.85);
	MakeNiceLegend(legMB,0.035,2);
	cMB->SetLogy();
	MakeNiceHistogram(hMB[10],col[10]); hMB[10]->Scale(4.);
	legMB->AddEntry(hMB[10],Form("LHC18%s",name[10].Data()),"pl");
	hMB[10]->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeNiceHistogram(hMB[i],col[i]); 
		legMB->AddEntry(hMB[i],Form("LHC18%s",name[i].Data()),"pl");
		hMB[i]->Draw("same");
	}
	legMB->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeRatioPlot(hMB[i],hMB[10],cMB,0.9,1.1,0.,14.);
	}
	cMB->Draw();
	cMB->SaveAs("cMB.png");

	////////////

	cout << "Drawing V0M" << endl;
	TCanvas* cV0M = new TCanvas("cV0M","cV0M",1200,1200);
	TLegend* legV0M = new TLegend(0.48,0.64,0.81,0.85);
	MakeNiceLegend(legV0M,0.035,2);
	cV0M->SetLogy();
	MakeNiceHistogram(hV0M[10],col[10]); hV0M[10]->Scale(4.);
	legV0M->AddEntry(hV0M[10],Form("LHC18%s",name[10].Data()),"pl");
	hV0M[10]->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeNiceHistogram(hV0M[i],col[i]); 
		legV0M->AddEntry(hV0M[i],Form("LHC18%s",name[i].Data()),"pl");
		hV0M[i]->Draw("same");
	}
	legV0M->Draw();
	for (int i = 0; i < 12; ++i)	{
		if (i==10) continue;
		if (i==4||i==5||i==6||i==7||i==8||i==9) continue;
		MakeRatioPlot(hV0M[i],hV0M[10],cV0M,0.9,1.1,0.,14.);
	}
	cV0M->Draw();
	cV0M->SaveAs("cV0M.png");


}

	/*
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

}*/