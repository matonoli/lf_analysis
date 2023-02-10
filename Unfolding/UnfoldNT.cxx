#include "UnfoldNTclass.h"

TFile* fMC = nullptr;
TFile* fData = nullptr;
TList* lMC = nullptr;
TList* lData = nullptr;

void ReadFiles(const string&, const string&, const string&);
void FlipMatrix(TH2F* );
void FlipPtvNt(TH2F* );
void Undolf1D(const char* , bool );
void Undolf2DMC(const char* );
void UndolfV02D(const char* );
void Undolf2D(const char*,string,string,string);
void ComparisonPublished(TH1F* );
const bool eRM = 0;
const int NumberOfIters = 5;

void UnfoldNT(const char* fData, const char* fMC, const char* dOut, const bool isMC, string trackcuts, string detector, string charge)
{

	ReadFiles(fMC,fData,trackcuts);
	std::cout << "Files read in \n";

	if(strcmp(charge.c_str(),"Charged")==0) Undolf1D(dOut,isMC);
	std::cout << "Unfolding 1D performed \n";
	//UndolfV02D(dOut);
	//if(isMC) Undolf2DMC(dOut);
	//if(!isMC) Undolf2D(dOut,trackcuts,detector,charge);

}

void FlipMatrix(TH2F* h) {

	TH2F* htmp = (TH2F*)h->Clone("hRes");
	for (int i = 1; i < h->GetNbinsX()+1; i++) {
	for (int j = 1; j < h->GetNbinsY()+1; j++) {
			htmp->SetBinContent(i,j,h->GetBinContent(j,i));
			htmp->SetBinError(i,j,h->GetBinError(j,i));
	}	}
	htmp->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(h->GetXaxis()->GetTitle());
	delete h;
	h = htmp; h->SetName("hRes");
}

void FlipPtvNt(TH2F* h) {

	const Int_t NPTBINS = 16;
	const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 5.0, 
		6.5, 8.0 };

	TH2F* htmp = new TH2F(h->GetName(),";N_T;p_{T} (Gev/#it{c}); Entries", 50, -0.5, 49.5, NPTBINS, XBINS);
	for (int i = 1; i < h->GetNbinsX()+1; i++) {
	for (int j = 1; j < h->GetNbinsY()+1; j++) {
			htmp->SetBinContent(j,i,h->GetBinContent(i,j));
			htmp->SetBinError(j,i,h->GetBinError(i,j));
	}	}
	delete h;
	h = htmp;
}

void UndolfV02D(const char* dirOut)
{

	TFile* fOut = new TFile(Form("./%s/2D_newClass_mc.root",dirOut),"UPDATE");

	TList* lOut = new TList();
	lOut->SetOwner();

	UnfoldNTclass* dobj = new UnfoldNTclass();

	TH2F* hRM = (TH2F*)fMC->Get("MyAnalysisV0correct_2/hNchTransRCvMC");
	TH2F* htmp = (TH2F*)hRM->Clone("hRes");
	for (int i = 1; i < hRM->GetNbinsX()+1; i++) {
	for (int j = 1; j < hRM->GetNbinsY()+1; j++) {
			htmp->SetBinContent(i,j,hRM->GetBinContent(j,i));
			htmp->SetBinError(i,j,hRM->GetBinError(j,i));
	}	}
	htmp->GetXaxis()->SetTitle(hRM->GetYaxis()->GetTitle());
	htmp->GetYaxis()->SetTitle(hRM->GetXaxis()->GetTitle());
	delete hRM;
	hRM = htmp; hRM->SetName("hRes");

	if(eRM) dobj->ExtrapolateRM(hRM);

	TGraph* gYieldvsNT = nullptr;

	const char* Regions[3] = {"Toward","Away","Trans1D"};
	const char* V0RegNames[3] = {"Near","Away","Trans"};
	//const char* Regions[3] = {"Transverse","Toward","Away"};
	//const char* V0RegNames[3] = {"Trans","Near","Away"};

	for(int region = 0; region < 3; ++region){

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[region]);
		obj->SetnIter(NumberOfIters);
		obj->SetDetector("TPC");
		obj->SetPid("K0s");
		obj->SetPidIdx(4);
		obj->SetCharge("Charged");
		obj->SetMCAnalysis(kTRUE);
		obj->SaveSolutionNT(kFALSE);

		TH2F* h2DGe = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hV0PtNt_K0s_MC_%s",V0RegNames[region]));
		TH2F* h2DMC = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hV0PtNtFitCorr_K0s_RC_%s",V0RegNames[region]));

		TH2F* hLGen = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hV0PtNt_L_MC_%s",V0RegNames[region]));
		TH2F* hLRec = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hV0PtNtFitCorr_L_RC_%s",V0RegNames[region]));
		TH2F* hLbarGen = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hV0PtNt_Lbar_MC_%s",V0RegNames[region]));
		TH2F* hLbarRec = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hV0PtNtFitCorr_Lbar_RC_%s",V0RegNames[region]));
		hLGen->Add(hLbarGen); hLRec->Add(hLbarRec);

		TH2F* hPiGen = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hPiPtNtMC_%s",V0RegNames[region]));
		TH2F* hPiRec = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hPiPtNtRC_%s",V0RegNames[region]));
		TH2F* hKpmGen = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtMC_%s",V0RegNames[region]));
		TH2F* hKpmRec = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtRC_%s",V0RegNames[region]));

		/*if (region == 2) {
			hPiGen = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hPiPtNtMC_%s","Trans"));
			hPiRec = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hPiPtNtRC_%s","Trans"));
			hKpmGen = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtMC_%s","Trans"));
			hKpmRec = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtRC_%s","Trans"));
		}*/

		TH2F* hNchGen = (TH2F*)hPiGen->Clone("hNchGen"); hNchGen->Add(hKpmGen);
		TH2F* hNchRec = (TH2F*)hPiRec->Clone("hNchRec"); hNchRec->Add(hKpmRec);

		/*if (region == 2) {
			hKpmGen = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtMC_%s","TransMin"));
			hKpmRec = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtRC_%s","TransMin"));
		}*/

		//h2DGe = hNchGen; h2DMC = hNchRec;
		//h2DGe = hKpmGen; h2DMC = hKpmRec;
		//h2DGe = hPiGen; h2DMC = hPiRec;
		//h2DGe = hLGen; h2DMC = hLRec;

		//TH2F* h2DGe = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtMC_%s",V0RegNames[region]));
		//TH2F* h2DMC = (TH2F*)fMC->Get(Form("MyAnalysisV0correct_2/hKpmPtNtRC_%s",V0RegNames[region]));

		//TRUNCATE TO DIMENSIONS OF RM

		//FlipPtvNt(h2DGe);
		//FlipPtvNt(h2DMC);

		printf(" - Region : %s\n",obj->GetRegion());
		cout << "RM has " << hRM->GetNbinsX() << " x " << hRM->GetNbinsY() << endl;
		cout << "hMC has " << h2DGe->GetNbinsX() << endl;
		obj->UnfoldV02D(hRM,hNchRec,h2DMC);
		//obj->UnfoldV02D(hRM,h2DMC,h2DMC);
		//! Uncomment if MCclosure is desired
		//obj->GetMCclosureinNchBins(h2DGe);

		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		cout << "1 Finding " << h2DGe << " and " << (TH2F*)Arr->FindObject("_hPtvsNch") << "\n";
		obj->GetMCclosureinRTBins(h2DGe,(TH2F*)Arr->FindObject("_hPtvsNch"));	// first is Gen, second is Unf

		fOut->cd();


		TDirectory* dir = fOut->mkdir(Form("%s",Regions[region]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	// region

	delete dobj;

	//	lOut->Add(gYieldvsNT);
	lOut->Write();
	//fOut->Close();


}


void Undolf1D(const char* dOut, bool isMC)
{

	TList* lOut = new TList();
	lOut->SetOwner();

	std::cout << "Got here 1 \n";

	TH1F* hRec = (TH1F*)fData->Get("MyAnalysisV0correct_2/hNchTrans");
	TH1F* hMeas = (TH1F*)fMC->Get("MyAnalysisV0correct_2/hNchTransRC");
	TH1F* hGen = (TH1F*)fMC->Get("MyAnalysisV0correct_2/hNchTransMC");
	TH2F* hRes = (TH2F*)fMC->Get("MyAnalysisV0correct_2/hNchTransRCvMC");

	// X-AXIS NEEDS TO BE RC, Y-AXIS MC
	TH2F* hResFlip = (TH2F*)hRes->Clone("hResFlip");
	for (int i = 1; i < hRes->GetNbinsX()+1; i++) {
	for (int j = 1; j < hRes->GetNbinsY()+1; j++) {
			hResFlip->SetBinContent(i,j,hRes->GetBinContent(j,i));
			hResFlip->SetBinError(i,j,hRes->GetBinError(j,i));
	}	}
	hResFlip->GetXaxis()->SetTitle(hRes->GetYaxis()->GetTitle());
	hResFlip->GetYaxis()->SetTitle(hRes->GetXaxis()->GetTitle());
	delete hRes;
	hRes = hResFlip; hRes->SetName("hRes");

	// RES MATRIX IS OBTAINED FROM CORRELATION BY ROWWISE NORMALISATION
	/*Int_t nCols = hRes->GetNbinsX();
	Int_t nRows = hRes->GetNbinsY();
	for (int iR = 1; iR < nRows+1; ++iR)	{
		Double_t integral = hRes->Integral(1,nCols+1,iR,iR);
		for (int iC = 1; iC < nCols+1; ++iC)	{
			Double_t binContent = hRes->GetBinContent(iC,iR);
			if (binContent>0) hRes->SetBinContent(iC,iR,binContent/integral);
		}
	}*/ // ITS ALSO DONE IN THE UNFOLDING CLASS

	std::cout << "Histograms: " << hRec << " " << hMeas << " " << hGen << " " << hRes << "\n";

	std::cout << "Got here 2 \n";
	UnfoldNTclass* obj = new UnfoldNTclass();
	std::cout << "Got here 3 \n";

	if(eRM) obj->ExtrapolateRM(hRes);
	std::cout << "Got here 4 \n";

	obj->SetnIter(NumberOfIters);
	obj->SetError(hMeas);
	obj->SetError(hGen);
	obj->SaveSolutionNT(kTRUE);
	std::cout << "Got here 5 \n";
	if(isMC) obj->Setup(hMeas,hGen,hRes);
	else obj->Setup(hRec,hGen,hRes);
	std::cout << "Got here 6 \n";
	obj->Unfold();
	obj->V2H();
	std::cout << "Got here 7 \n";

	printf(" - Region : %s\n",obj->GetRegion());

	TH1F* hUnfolded = nullptr;
	TH1F* hClosure = nullptr;
	if(isMC){
		hUnfolded = (TH1F*)obj->GetUnfoldedDistH();
		hClosure = (TH1F*)obj->GetClosureH();
	}else{
		hUnfolded = (TH1F*)obj->GetUnfoldedDistH();
	}
	std::cout << "Got here 8 \n";
	TH1F* hNT = (TH1F*)hUnfolded->Clone("_hNT");
	TH1F* hRT = (TH1F*)obj->RebinNT2RT(hNT, kTRUE);
	hRT->Scale(1.0/hRT->Integral());
	std::cout << "Got here 9 \n";
	//! Relative statistical uncertainty
	//! NT distribution with final bins to plot
	TH1F* hRelStatUnc = (TH1F*)hNT->Clone("hRelStatUnc");
	hRelStatUnc->Reset();

	TH1F* hNch = (TH1F*)hNT->Clone("hNT");
	hNch->Reset();

	std::cout << "Got here 10 \n";
	for(int bin = 1; bin <= hRelStatUnc->GetNbinsX(); bin++){

		double yield = hNT->GetBinContent(bin);
		double error = hNT->GetBinError(bin);

		if( yield > 0. ) hRelStatUnc->SetBinContent(bin, error / yield);
//		if( bin > 32 ) continue;
		hNch->SetBinContent(bin, yield);
		hNch->SetBinError(bin, error);
	}

	std::cout << "Got here 11 \n";
	lOut->Add(hRes);
	if(!isMC) ComparisonPublished(hRT);
	if(isMC) lOut->Add(hClosure);
	std::cout << "Got here 12 \n";
	lOut->Add(hUnfolded);
	//lOut->Add(hNT);
	lOut->Add(hRT);
	lOut->Add(hNch);
	lOut->Add(hRelStatUnc);

	TGraph* gChi2 = nullptr;
	if(isMC) {
		obj->DrawNchClosure(hGen,hUnfolded,-0.5,30.0,"#it{N}_{t}","Unfolded/True",dOut);
		TVectorD iterX = (TVectorD)obj->GetChi2Vectors(kTRUE);
		TVectorD Chi2 = (TVectorD)obj->GetChi2Vectors(kFALSE);
		gChi2 = new TGraph(iterX,Chi2);
		gChi2->SetName("gChi2");
		gChi2->SetMarkerStyle(8);
		lOut->Add(gChi2);
	}
	std::cout << "Got here 13 \n";

	TFile* fOut = nullptr;
	if(isMC)
		fOut = new TFile(Form("./%s/1D_newClass_mc.root",dOut),"RECREATE");
	else
		fOut = new TFile(Form("./%s/1D_newClass_data.root",dOut),"RECREATE");

	std::cout << "Got here 14 \n";

	fOut->cd();
	lOut->Write();
	(obj->GetObjArray())->Write();
	fOut->Close();

	std::cout << "Got here 15 \n";

	delete obj;
	delete fOut;

	obj = nullptr;
	fOut = nullptr;

}

void ComparisonPublished(TH1F* hRT)
{

	TFile* fIn = new TFile("./results/HEPData-ins1762350-v1-Probability_distribution_as_function_of_R_T.root","READ");
	TDirectory* dIn = (TDirectory*)fIn->Get("Probability distribution as function of R_T");
	if(!dIn)printf("Probability distribution as function of R_T\n");

	TH1F* hstat = (TH1F*)dIn->Get("Hist1D_y1");
	TH1F* he = (TH1F*)dIn->Get("Hist1D_y1_e1");

	for(int bin = 1; bin <= hstat->GetNbinsX(); bin++)
		hstat->SetBinError(bin,he->GetBinError(bin));

	for(int bin = 1; bin <= hRT->GetNbinsX(); bin++)
		hRT->SetBinContent(bin,hRT->GetBinContent(bin)/hRT->GetBinWidth(bin));

	hRT->Divide(hRT,hstat,1.0,1.0,"");

}

void Undolf2DMC( const char* dirOut)
{

	TFile* fOut = new TFile(Form("./%s/2D_newClass_mc.root",dirOut),"UPDATE");

	TList* lOut = new TList();
	lOut->SetOwner();

	UnfoldNTclass* dobj = new UnfoldNTclass();

	TH2F* hRM = (TH2F*)lMC->FindObject("hNchResponse");
	if(eRM) dobj->ExtrapolateRM(hRM);

	TGraph* gYieldvsNT = nullptr;

	const char* Regions[3] = {"Toward","Away","Transverse"};
	for(int region = 0; region < 3; ++region){

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[region]);
		obj->SetnIter(NumberOfIters);
		obj->SetDetector("TPC");
		obj->SetPid("Charged");
		obj->SetCharge("Charged");
		obj->SetMCAnalysis(kTRUE);
		obj->SetPidIdx(0);
		obj->SaveSolutionNT(kFALSE);

		TH2F* h2DGe = (TH2F*)lMC->FindObject(Form("hNchGenVsPtGen_%s_Charged",Regions[region]));
		TH2F* h2DMC = (TH2F*)lMC->FindObject(Form("hNchVsPt_%s_Charged",Regions[region]));

		printf(" - Region : %s\n",obj->GetRegion());
		obj->Unfold2D(hRM,h2DMC);
		//! Uncomment if MCclosure is desired
		obj->GetMCclosureinNchBins(h2DGe);

		TObjArray* Arr = (TObjArray*)obj->GetObjArray();
		obj->GetMCclosureinRTBins(h2DGe,(TH2F*)Arr->FindObject("_hPtvsNch")); //_hPtvsNch contains the unfolded results

		TDirectory* dir = fOut->mkdir(Form("%s_Charged",Regions[region]));
		dir->cd();
		(obj->GetObjArray())->Write();

		delete obj;
		obj = nullptr;

	}	// region

	delete dobj;

	//	lOut->Add(gYieldvsNT);
	lOut->Write();
	fOut->Close();

}

void MCclosure(const char* Region)
{

	TFile* fGen = new TFile("~/Documents/DataAnalysis/pp_13TeV_spectra_RT/Data_MC/mc/one_hybrid/AnalysisResults.root","READ");
	TList* dGen = (TList*)fGen->Get("MyOutputContainer_PID_5GeV_Hybrid");

	TFile* fUn1D = new TFile("~/Documents/DataAnalysis/pp_13TeV_spectra_RT/analysis/results/1D_newClass_mc.root","READ");
	if(!fUn1D)printf("fUn1D\n");
	TH1F* hNTUnf = (TH1F*)fUn1D->Get("_hUnfDist");
	if(!hNTUnf) printf("hNTUnf\n");

	TFile* fUn = new TFile("~/Documents/DataAnalysis/pp_13TeV_spectra_RT/analysis/results/2D_newClass_mc.root","READ");
	TDirectory* dUn = (TDirectory*)fUn->Get(Form("%s_Charged",Region));
	TH2F* hUnf2D = (TH2F*)dUn->Get("_hPtvsNch");

	TH2F* hGen2D = (TH2F*)dGen->FindObject(Form("hNchGenVsPtGen_%s_Charged",Region));
	TH1F* hNTGen = (TH1F*)dGen->FindObject("hMultTSGen");
	TH1D* hUn = 0x0;
	TH1D* hpTgen = 0x0;
	TH1D* hpTunf = 0x0;
	TList* lOut = new TList();

	TH1D* hY = (TH1D*)hGen2D->ProjectionX("hY",1,1);
	hY->Reset();

	TGraphErrors* gr = new TGraphErrors();
	gr->SetName("IntegratedYield");
	gr->SetMarkerStyle(8);

	for(int bin = 1; bin <= hGen2D->GetNbinsX(); bin++){

		hpTunf = (TH1D*)hUnf2D->ProjectionX(Form("hPt_Corrected_Bin_nch_%d",bin),bin,bin);
		if(!hpTunf) continue;
		hpTgen = (TH1D*)hGen2D->ProjectionY(Form("hptgen_%d",bin),bin,bin);

		if(bin > 35) continue;

		double eyu = 0.0;
		double eyg = 0.0;

		double yu = 0.0;
		double yg = 0.0;
		double nevsge = hNTGen->GetBinContent(bin);
		double nevsun = hNTUnf->GetBinContent(bin);

		//yg = hpTgen->IntegralAndError(1,hpTgen->GetNbinsX(),eyg,"width");
		//yu = hpTunf->IntegralAndError(1,hpTunf->GetNbinsX(),eyu,"width");

		for(int j = 1; j <= hpTunf->GetNbinsX(); j++){
			yg += hpTgen->GetBinContent(j);
			eyg += TMath::Sqrt(hpTgen->GetBinContent(j)) * TMath::Sqrt(hpTgen->GetBinContent(j));
		}

		yg /= nevsge;

		for(int j = 1; j <= hpTunf->GetNbinsX(); j++){
			yu += hpTunf->GetBinContent(j);
			eyu += hpTunf->GetBinError(j) * hpTunf->GetBinError(j);
		}

		yu /= nevsun;
		eyu = TMath::Sqrt(eyu)/nevsun;
		eyg = TMath::Sqrt(eyg)/nevsge;

		//printf("Nch = %d Yu = %f yg = %f\n",bin-1, yu, yg);
		//printf("\t eYu = %f eYg = %f\n", eyu, eyg);
		gr->SetPoint(bin-1,yg,yu);
		gr->SetPointError(bin-1,0.0,eyg);

		if(bin>1){
			double ratio = yu/(bin - 1);
			double eratio = eyg/(bin - 1);
			//printf("ratio = %f @ bin = %d\n",ratio, bin );

			hY->SetBinContent(bin, ratio );
			hY->SetBinError(bin, eratio );
		}

	}


	lOut->Add(gr);
	lOut->Add(hY);

	TFile* fOut = new TFile("./results/extraMCclosure.root","UPDATE");
	fOut->cd();

	TDirectory* dOut = (TDirectory*)fOut->mkdir(Form("%s",Region));
	dOut->cd();
	lOut->Write();



	delete fOut;

}

void Undolf2D(const char* dirOut, string trackcuts, string detector, string charge)
{

	TFile* fOut = new TFile(Form("./%s/2D_newClass_data.root",dirOut),"UPDATE");

	TList* lOut = new TList();
	lOut->SetOwner();

	UnfoldNTclass* dobj = new UnfoldNTclass();

	TH2F* hRM = (TH2F*)lMC->FindObject("hNchResponse");
	if(eRM) dobj->ExtrapolateRM(hRM);
	const char* Regions[3] = {"Toward","Away","Transverse"};

	for(int region = 0; region < 3; ++region){

		UnfoldNTclass* obj = new UnfoldNTclass();

		obj->SetRegion(Regions[region]);
		obj->SetDetector("TPC");
		obj->SetPid("Charged");
		obj->SetCharge(charge.c_str());
		obj->SetnIter(NumberOfIters);
		obj->SetMCAnalysis(kFALSE);
		obj->SaveSolutionNT(kFALSE);
		obj->SetPidIdx(0);

		TH2F* hPtvsNchRec = nullptr;
		TH2F* hPtvsNchRecEta_02 = nullptr;
		TH2F* hPtvsNchRecEta_24 = nullptr;
		TH2F* hPtvsNchRecEta_46 = nullptr;
		TH2F* hPtvsNchRecEta_68 = nullptr;

		if(strcmp(charge.c_str(),"Charged")==0){

			TH2F* hPos_02 = (TH2F*)lData->FindObject(Form("hNchVsPtPosTPC_%s_02",Regions[region]));
			TH2F* hNeg_02 = (TH2F*)lData->FindObject(Form("hNchVsPtNegTPC_%s_02",Regions[region]));
			TH2F* hPos_24 = (TH2F*)lData->FindObject(Form("hNchVsPtPosTPC_%s_24",Regions[region]));
			TH2F* hNeg_24 = (TH2F*)lData->FindObject(Form("hNchVsPtNegTPC_%s_24",Regions[region]));
			TH2F* hPos_46 = (TH2F*)lData->FindObject(Form("hNchVsPtPosTPC_%s_46",Regions[region]));
			TH2F* hNeg_46 = (TH2F*)lData->FindObject(Form("hNchVsPtNegTPC_%s_46",Regions[region]));
			TH2F* hPos_68 = (TH2F*)lData->FindObject(Form("hNchVsPtPosTPC_%s_68",Regions[region]));
			TH2F* hNeg_68 = (TH2F*)lData->FindObject(Form("hNchVsPtNegTPC_%s_68",Regions[region]));

			hPtvsNchRecEta_02 = obj->SumCharges2D("02",hPos_02,hNeg_02);
			hPtvsNchRecEta_24 = obj->SumCharges2D("24",hPos_24,hNeg_24);
			hPtvsNchRecEta_46 = obj->SumCharges2D("46",hPos_46,hNeg_46);
			hPtvsNchRecEta_68 = obj->SumCharges2D("68",hPos_68,hNeg_68);

		}else{
			hPtvsNchRecEta_02 = (TH2F*)lData->FindObject(Form("hNchVsPt%sTPC_%s_02",charge.c_str(),Regions[region]));
			hPtvsNchRecEta_24 = (TH2F*)lData->FindObject(Form("hNchVsPt%sTPC_%s_24",charge.c_str(),Regions[region]));
			hPtvsNchRecEta_46 = (TH2F*)lData->FindObject(Form("hNchVsPt%sTPC_%s_46",charge.c_str(),Regions[region]));
			hPtvsNchRecEta_68 = (TH2F*)lData->FindObject(Form("hNchVsPt%sTPC_%s_68",charge.c_str(),Regions[region]));
		}

		hPtvsNchRec = (TH2F*)obj->SumDistPerEta2D(hPtvsNchRecEta_02,hPtvsNchRecEta_24,hPtvsNchRecEta_46,hPtvsNchRecEta_68);

		TDirectory* dOut = fOut->mkdir(Form("%s_PID_Charged_charge_%s",Regions[region],charge.c_str()));
		obj->Unfold2D(hRM,hPtvsNchRec);
		dOut->cd();
		(obj->GetObjArray())->Write();

		delete obj;

	}	// region

	delete dobj;

	lOut->Write();
	fOut->Close();
	delete fOut;

}

void ReadFiles(const string& pMC, const string& pData, const string& trackcuts)
{

	fMC = new TFile(Form("%s",pMC.c_str()),"READ");
	fData = new TFile(Form("%s",pData.c_str()),"READ");
	if(!fMC) printf("COULD NOT FIND: %s\n",pMC.c_str());
	if(!fData) printf("COULD NOT FIND: %s\n",pData.c_str());

	lMC = (TList*)fMC->Get(Form("MyOutputContainer_PID_5GeV_%s",trackcuts.c_str()));
	lData = (TList*)fData->Get(Form("MyOutputContainer_PID_5GeV_%s",trackcuts.data()));
	if(!lMC) printf("lMC NOT FOUND\n");
	if(!lData) printf("lData NOT FOUND\n");

}
