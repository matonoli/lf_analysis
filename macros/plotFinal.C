#include <iostream>

using namespace std;

void MakeRatioPlot(TH1D* hn, TH1D* hd, TCanvas* c, Double_t low, Double_t high, Double_t lowx, Double_t highx, const char* opt) {
  
  c->cd();
  TString strOpt(opt);

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

  // for ratios of S0 to HM, use the unc. uncertainties as systematics
  TString hrname(hr->GetName());
  if (hrname.Contains("Unc")) {
  for (int iBin = 1; iBin < hn->GetNbinsX()+1; iBin++) {
    double nerr = (double)hn->GetBinError(iBin)/hn->GetBinContent(iBin);
    hr->SetBinError(iBin,(double)hr->GetBinContent(iBin)*nerr);
  } }


  hr->GetYaxis()->SetTitle("ratio to ref.");
  hr->GetYaxis()->CenterTitle();
  hr->GetYaxis()->SetNdivisions(505);
  hr->GetYaxis()->SetTitleSize(20);
  hr->GetYaxis()->SetTitleFont(43);
  hr->GetYaxis()->SetLabelFont(43); 

  hr->GetYaxis()->SetTitleOffset(1.6);
  hr->GetYaxis()->SetLabelOffset(0.0025);
  hr->GetYaxis()->SetLabelSize(20);

  hr->GetXaxis()->SetTitleSize(25);
  hr->GetXaxis()->SetTitleFont(43);
  hr->GetXaxis()->SetTitleOffset(4.);
  hr->GetXaxis()->SetLabelFont(43); 
  hr->GetXaxis()->SetLabelSize(25);
  hr->GetXaxis()->SetTickLength(0.09);

  if (!hasRatio)  hr->Draw(strOpt.Data());
  else      hr->Draw(Form("same %s", strOpt.Data()));
  if (! hrname.Contains("hBlank")) hr->Write();

  //c->SetCanvasSize()
  c->cd();

}

TH1D* GetFrame(const char* name, double min, double max, bool isRatios)
{
    
    TH1D* hframe = new TH1D(Form("%s",name),Form("%s",name),(int)max-(int)min,min,max);
    
    hframe->GetYaxis()->SetTitleOffset(2);
    hframe->GetYaxis()->SetNdivisions(510,kTRUE);
    hframe->GetXaxis()->SetNdivisions(505,kTRUE);
    hframe->GetXaxis()->SetLabelFont(63);
    hframe->GetYaxis()->SetLabelFont(63);
    hframe->GetXaxis()->SetTitleFont(63);
    hframe->GetYaxis()->SetTitleFont(63);
    hframe->GetYaxis()->SetTitleSize( 30 );
    hframe->GetXaxis()->SetTitleSize( 30 );
    hframe->GetYaxis()->SetLabelSize( 25 );
    
    if(isRatios){
        hframe->GetYaxis()->SetNdivisions(505,kTRUE);
        hframe->GetYaxis()->SetTitleOffset(1.5);
        hframe->GetYaxis()->SetLabelSize( 20 );
        hframe->GetYaxis()->SetTitleSize( 30 );
    }
    
    if(strstr(name,"frame2")!=0){
        hframe->GetXaxis()->SetLabelSize( 25 );
        hframe->GetYaxis()->SetLabelSize( 25 );
        
        hframe->GetYaxis()->SetTitleOffset(2);
        hframe->GetXaxis()->SetTitleOffset(3);
        hframe->GetYaxis()->SetRangeUser(0.0,2.82);
        hframe->GetYaxis()->SetNdivisions(505,kTRUE);
        hframe->GetXaxis()->SetNdivisions(505,kTRUE);
        hframe->GetYaxis()->SetTitle("Ratio to #it{N}_{ch} #geq 10");
        hframe->GetXaxis()->SetTitle(titleX);
        
        if(isRatios){
            
            hframe->GetYaxis()->SetNdivisions(505,kTRUE);
            hframe->GetXaxis()->SetNdivisions(505,kTRUE);
            hframe->GetXaxis()->SetLabelSize( 20 );
            hframe->GetYaxis()->SetLabelSize( 20 );
            hframe->GetXaxis()->SetTitleSize( 24 );
            hframe->GetYaxis()->SetTitleSize( 25 );
            hframe->GetYaxis()->SetRangeUser(0.65,1.35);
            hframe->GetXaxis()->SetTitleOffset(3);
            hframe->GetYaxis()->SetTitleOffset(1.5);
            
        }
        
    }
    return hframe;
}
//______________________________________________________________
void myLegendSetUp(TLegend *currentLegend,float currentTextSize){
    currentLegend->SetTextFont(42);
    currentLegend->SetBorderSize(0);
    currentLegend->SetFillStyle(0);
    currentLegend->SetFillColor(0);
    currentLegend->SetMargin(0.25);
    currentLegend->SetTextSize(currentTextSize);
//    currentLegend->SetEntrySeparation(0.35);
    return;
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}

void MakeNiceHistogram(TH1D* h, Int_t col) {

  h->SetLineColor(col);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.0);
  h->SetMarkerColor(col);
  h->SetStats(0);
  h->SetFillStyle(0);

  h->GetYaxis()->SetTitleSize(30);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetLabelOffset(0.0025);
  h->GetYaxis()->SetLabelSize(20);
  h->GetYaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelFont(43);

  //h->SetTopMargin(0.055);
}


TString mult("V0M01"); int multInt = 0; int multOmar = 0; TString multPlot("V0M 0-1%");
//TString mult("NCharged01"); int multInt = 1; int multOmar = 0; TString multPlot("CL1 0-1%");
//TString mult("V0M"); int multInt = 2; int multOmar = 1; TString multPlot("V0M 0-10%");
//TString mult("NCharged"); int multInt = 3; int multOmar = 1; TString multPlot("CL1 0-10%");
//TString spher("20"); int sphInt = 3;
TString spher("10"); int sphInt = 2;
//TString spher("5"); int sphInt = 1;
//TString spher("1"); int sphInt = 0;

Bool_t IS_EXCLUSIVE = 0;


const Int_t NPTBINS = 13;
  const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
    1.0, 1.2, 
    1.4, 1.6, 1.8, 2.0, 2.2, 
    2.6, 3.0, 3.4, 4.0, 5.0, 
    6.5, 8.0 };


/*const Int_t NPTBINS = 16;
  const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
    0.4, 0.6, 0.8, 1.0, 1.2, 
    1.4, 1.6, 1.8, 2.0, 2.2, 
    2.6, 3.0, 3.4, 4.0, 5.0, 
    6.5, 8.0 };*/

const Int_t NPIONBINS = 51;   //Omar pi+- spectra
  const Double_t PIONBINS[NPIONBINS+1] = {
    0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 
    0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 
    0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 
    1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 
    2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.50, 5.00, 
    5.50, 6.00, 6.50, 7.00, 8.00, 10.0, 20.0 };

Int_t colors[3] = {kBlack, kRed, kBlue};


void plotFinal() {

  TH1::SetDefaultSumw2(1);

	TFile* fin = new TFile("output210826all.root", "READ");
  //fin = (IS_EXCLUSIVE) ? new TFile("output211114all_excl.root", "READ") : new TFile("output210826all.root", "READ");
  fin = (IS_EXCLUSIVE) ? new TFile("output220223_ex.root", "READ") : new TFile("output220223.root", "READ");
  TDirectoryFile* dirSyst = (TDirectoryFile*)fin->Get("MyAnalysisV0syst_3");
  TDirectoryFile* dirCorr = (TDirectoryFile*)fin->Get("MyAnalysisV0correct_2"); 

  TFile* fpion = new TFile("../official/pi_k_p_spectra.root","READ");
  const char* strM[4] = {"V0M", "Trks", "V0M", "Trks"};
  const char* strS[3] = {"Jetty", "Reference", "Isotropic"};
  TDirectoryFile* dirPionRef = (TDirectoryFile*)fpion->Get(Form("%s_Reference_Bin_mult_%i",strM[multInt],multOmar));
  TDirectoryFile* dirPionJet = (TDirectoryFile*)fpion->Get(Form("%s_Jetty_%i_Bin_mult_%i",strM[multInt],sphInt,multOmar));
  TDirectoryFile* dirPionIso = (TDirectoryFile*)fpion->Get(Form("%s_Isotropic_%i_Bin_mult_%i",strM[multInt],sphInt,multOmar));
  TDirectoryFile* dirPion[3] = { dirPionRef, dirPionJet, dirPionIso };

  TFile* fmc = new TFile("../official/Monash_Spectra_CL1_10_s0rt_01.root","READ");
  TH1D* hLLbarToPi_Ref_MC = (TH1D*)fmc->Get(Form("lambda_ref_ratio_CL110_S001"));
  TH1D* hLLbarToPi_Jet_MC = (TH1D*)fmc->Get(Form("lambda_jet_ratio_CL110_S001"));
  TH1D* hLLbarToPi_Iso_MC = (TH1D*)fmc->Get(Form("lambda_iso_ratio_CL110_S001"));
  TH1D* hK0s_Ref_MC = (TH1D*)fmc->Get("k0s_ref_CL110_S001");
  TH1D* hK0s_Jet_MC = (TH1D*)fmc->Get("k0s_jet_CL110_S001");
  TH1D* hK0s_Iso_MC = (TH1D*)fmc->Get("k0s_iso_CL110_S001");
  TH1D* hKpm_Ref_MC = (TH1D*)fmc->Get("pion_ref_CL110_S001");
  TH1D* hKpm_Jet_MC = (TH1D*)fmc->Get("pion_jet_CL110_S001");
  TH1D* hKpm_Iso_MC = (TH1D*)fmc->Get("pion_iso_CL110_S001");
  TH1D* hLLbar_Ref_MC = (TH1D*)fmc->Get("lambda_ref_CL110_S001");
  TH1D* hLLbar_Jet_MC = (TH1D*)fmc->Get("lambda_jet_CL110_S001");
  TH1D* hLLbar_Iso_MC = (TH1D*)fmc->Get("lambda_iso_CL110_S001");
  TH1D* hKMCs[6] = { hK0s_Ref_MC, hK0s_Jet_MC, hK0s_Iso_MC, hKpm_Ref_MC, hKpm_Jet_MC, hKpm_Iso_MC,};
  TH1D* hLLbarMCs[3] = { hLLbar_Ref_MC, hLLbar_Jet_MC, hLLbar_Iso_MC};
  for (int i=0; i<6; i++) {
    hKMCs[i]->Rebin(4);
    hKMCs[i]->Scale(1./4);
    MakeNiceHistogram(hKMCs[i],colors[i%3]);
    if (i>2) continue;
    hLLbarMCs[i]->Rebin(4);
    hLLbarMCs[i]->Scale(1./4);
    MakeNiceHistogram(hLLbarMCs[i],colors[i%3]);
  }
  //hKMCs[0]->Divide(hKMCs[3]);hKMCs[1]->Divide(hKMCs[4]);hKMCs[2]->Divide(hKMCs[5]);

  TFile* fout = (IS_EXCLUSIVE) ? new TFile(Form("k0s_l_spectra_%s_spher%s_excl.root",mult.Data(),spher.Data()),"RECREATE") : new TFile(Form("k0s_l_spectra_%s_spher%s.root",mult.Data(),spher.Data()),"RECREATE");

  // Load histograms
  TH1D* hK0s_Ref_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_K0s_D_%s_MB",mult.Data())); 
  TH1D* hK0s_Jet_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_K0s_D_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hK0s_Iso_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_K0s_D_%s_Iso%s",mult.Data(),spher.Data()));
  TH1D* hL_Ref_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_L_D_%s_MB",mult.Data()));
  TH1D* hL_Jet_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_L_D_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hL_Iso_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_L_D_%s_Iso%s",mult.Data(),spher.Data()));
  TH1D* hLbar_Ref_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_Lbar_D_%s_MB",mult.Data()));
  TH1D* hLbar_Jet_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_Lbar_D_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hLbar_Iso_Stat = (TH1D*)dirCorr->Get(Form("hV0PtFitCorr_Lbar_D_%s_Iso%s",mult.Data(),spher.Data()));

  TH1D* hK0s_Ref_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_K0s_%s_MB",mult.Data()));
  TH1D* hK0s_Jet_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_K0s_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hK0s_Iso_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_K0s_%s_Iso%s",mult.Data(),spher.Data()));
  TH1D* hL_Ref_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_L_%s_MB",mult.Data()));
  TH1D* hL_Jet_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_L_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hL_Iso_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_L_%s_Iso%s",mult.Data(),spher.Data()));
  TH1D* hLbar_Ref_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_Lbar_%s_MB",mult.Data()));
  TH1D* hLbar_Jet_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_Lbar_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hLbar_Iso_Syst = (TH1D*)dirSyst->Get(Form("hV0PtSysSum_Lbar_%s_Iso%s",mult.Data(),spher.Data()));

  TH1D* hK0s_Ref_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_K0s_%s_MB",mult.Data()));
  TH1D* hK0s_Jet_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_K0s_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hK0s_Iso_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_K0s_%s_Iso%s",mult.Data(),spher.Data()));
  TH1D* hL_Ref_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_L_%s_MB",mult.Data()));
  TH1D* hL_Jet_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_L_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hL_Iso_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_L_%s_Iso%s",mult.Data(),spher.Data()));
  TH1D* hLbar_Ref_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_Lbar_%s_MB",mult.Data()));
  TH1D* hLbar_Jet_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_Lbar_%s_Jetty%s",mult.Data(),spher.Data()));
  TH1D* hLbar_Iso_SystUnc = (TH1D*)dirSyst->Get(Form("hV0PtSysSumUnc_Lbar_%s_Iso%s",mult.Data(),spher.Data()));

  TH1D* hStats[9] = {
    hK0s_Ref_Stat, hK0s_Jet_Stat, hK0s_Iso_Stat,
    hL_Ref_Stat, hL_Jet_Stat, hL_Iso_Stat,
    hLbar_Ref_Stat, hLbar_Jet_Stat, hLbar_Iso_Stat };

  const char* statNames[9] = {
    "hK0s_Ref_Stat", "hK0s_Jet_Stat", "hK0s_Iso_Stat",
    "hL_Ref_Stat", "hL_Jet_Stat", "hL_Iso_Stat",
    "hLbar_Ref_Stat", "hLbar_Jet_Stat", "hLbar_Iso_Stat" };

  TH1D* hSysts[9] = {
    hK0s_Ref_Syst, hK0s_Jet_Syst, hK0s_Iso_Syst,
    hL_Ref_Syst, hL_Jet_Syst, hL_Iso_Syst,
    hLbar_Ref_Syst, hLbar_Jet_Syst, hLbar_Iso_Syst };

  const char* systNames[9] = {
    "hK0s_Ref_Syst", "hK0s_Jet_Syst", "hK0s_Iso_Syst",
    "hL_Ref_Syst", "hL_Jet_Syst", "hL_Iso_Syst",
    "hLbar_Ref_Syst", "hLbar_Jet_Syst", "hLbar_Iso_Syst" };

  TH1D* hSystsUnc[9] = {
    hK0s_Ref_SystUnc, hK0s_Jet_SystUnc, hK0s_Iso_SystUnc,
    hL_Ref_SystUnc, hL_Jet_SystUnc, hL_Iso_SystUnc,
    hLbar_Ref_SystUnc, hLbar_Jet_SystUnc, hLbar_Iso_SystUnc };

  const char* systUncNames[9] = {
    "hK0s_Ref_SystUnc", "hK0s_Jet_SystUnc", "hK0s_Iso_SystUnc",
    "hL_Ref_SystUnc", "hL_Jet_SystUnc", "hL_Iso_SystUnc",
    "hLbar_Ref_SystUnc", "hLbar_Jet_SystUnc", "hLbar_Iso_SystUnc" };


  for (int iH = 0; iH < 9; iH++) {

    for (int iB = 1; iB < hSysts[iH]->GetNbinsX()+1; iB++) {
      hSysts[iH]->SetBinError(iB, hSysts[iH]->GetBinContent(iB)*hStats[iH]->GetBinContent(iB));
      hSysts[iH]->SetBinContent(iB, hStats[iH]->GetBinContent(iB));

      hSystsUnc[iH]->SetBinError(iB, hSystsUnc[iH]->GetBinContent(iB)*hStats[iH]->GetBinContent(iB));
      hSystsUnc[iH]->SetBinContent(iB, hStats[iH]->GetBinContent(iB));
    }

    hStats[iH]->SetName(statNames[iH]); hSysts[iH]->SetName(systNames[iH]); hSystsUnc[iH]->SetName(systUncNames[iH]);
    hSysts[iH]->GetYaxis()->SetTitle("Entries"); 
    hSysts[iH]->GetYaxis()->UnZoom(); 
    //hStats[iH]->Write(); hSysts[iH]->Write();
  }    


  // Add L and Lbar together
  TH1D* hLLbar_Ref_Stat = (TH1D*)hL_Ref_Stat->Clone("hLLbar_Ref_Stat");
  TH1D* hLLbar_Jet_Stat = (TH1D*)hL_Jet_Stat->Clone("hLLbar_Jet_Stat");
  TH1D* hLLbar_Iso_Stat = (TH1D*)hL_Iso_Stat->Clone("hLLbar_Iso_Stat");
  TH1D* hLLbar_Ref_Syst = (TH1D*)hL_Ref_Syst->Clone("hLLbar_Ref_Syst");
  TH1D* hLLbar_Jet_Syst = (TH1D*)hL_Jet_Syst->Clone("hLLbar_Jet_Syst");
  TH1D* hLLbar_Iso_Syst = (TH1D*)hL_Iso_Syst->Clone("hLLbar_Iso_Syst");
  TH1D* hLLbar_Ref_SystUnc = (TH1D*)hL_Ref_SystUnc->Clone("hLLbar_Ref_SystUnc");
  TH1D* hLLbar_Jet_SystUnc = (TH1D*)hL_Jet_SystUnc->Clone("hLLbar_Jet_SystUnc");
  TH1D* hLLbar_Iso_SystUnc = (TH1D*)hL_Iso_SystUnc->Clone("hLLbar_Iso_SystUnc");
  hLLbar_Ref_Stat->Add(hLbar_Ref_Stat);
  hLLbar_Jet_Stat->Add(hLbar_Jet_Stat);
  hLLbar_Iso_Stat->Add(hLbar_Iso_Stat);
  hLLbar_Ref_Syst->Add(hLbar_Ref_Syst);
  hLLbar_Jet_Syst->Add(hLbar_Jet_Syst);
  hLLbar_Iso_Syst->Add(hLbar_Iso_Syst);
  hLLbar_Ref_SystUnc->Add(hLbar_Ref_SystUnc);
  hLLbar_Jet_SystUnc->Add(hLbar_Jet_SystUnc);
  hLLbar_Iso_SystUnc->Add(hLbar_Iso_SystUnc);

  TH1D* hLLbar_Stat[3] = {
    hLLbar_Ref_Stat, hLLbar_Jet_Stat, hLLbar_Iso_Stat };
  TH1D* hLLbar_Syst[3] = {
    hLLbar_Ref_Syst, hLLbar_Jet_Syst, hLLbar_Iso_Syst };
  TH1D* hLLbar_SystUnc[3] = {
    hLLbar_Ref_SystUnc, hLLbar_Jet_SystUnc, hLLbar_Iso_SystUnc };


  // Make L/K0s ratios
  TH1D* hLLbarToK0s_Ref_Stat = (TH1D*)hLLbar_Ref_Stat->Clone("hLLbarToK0s_Ref_Stat");
  TH1D* hLLbarToK0s_Jet_Stat = (TH1D*)hLLbar_Jet_Stat->Clone("hLLbarToK0s_Jet_Stat");
  TH1D* hLLbarToK0s_Iso_Stat = (TH1D*)hLLbar_Iso_Stat->Clone("hLLbarToK0s_Iso_Stat");
  TH1D* hLLbarToK0s_Ref_Syst = (TH1D*)hLLbar_Ref_Syst->Clone("hLLbarToK0s_Ref_Syst");
  TH1D* hLLbarToK0s_Jet_Syst = (TH1D*)hLLbar_Jet_Syst->Clone("hLLbarToK0s_Jet_Syst");
  TH1D* hLLbarToK0s_Iso_Syst = (TH1D*)hLLbar_Iso_Syst->Clone("hLLbarToK0s_Iso_Syst");
  TH1D* hLLbarToK0s_Ref_SystUnc = (TH1D*)hLLbar_Ref_SystUnc->Clone("hLLbarToK0s_Ref_SystUnc");
  TH1D* hLLbarToK0s_Jet_SystUnc = (TH1D*)hLLbar_Jet_SystUnc->Clone("hLLbarToK0s_Jet_SystUnc");
  TH1D* hLLbarToK0s_Iso_SystUnc = (TH1D*)hLLbar_Iso_SystUnc->Clone("hLLbarToK0s_Iso_SystUnc");

  TH1D* hKL_Stat[3] = {
    hLLbarToK0s_Ref_Stat, hLLbarToK0s_Jet_Stat, hLLbarToK0s_Iso_Stat };
  TH1D* hKL_Syst[3] = {
    hLLbarToK0s_Ref_Syst, hLLbarToK0s_Jet_Syst, hLLbarToK0s_Iso_Syst };
  TH1D* hKL_SystUnc[3] = {
    hLLbarToK0s_Ref_SystUnc, hLLbarToK0s_Jet_SystUnc, hLLbarToK0s_Iso_SystUnc };

  for (int iH = 0; iH < 3; ++iH) {
    hKL_Stat[iH]->Divide(hKL_Stat[iH],hStats[iH],.5);
    hKL_Syst[iH]->Divide(hKL_Syst[iH],hSysts[iH],.5);
    hKL_SystUnc[iH]->Divide(hKL_SystUnc[iH],hSystsUnc[iH],.5);
  }

  hLLbarMCs[0]->Divide(hKMCs[0]);hLLbarMCs[1]->Divide(hKMCs[1]);hLLbarMCs[2]->Divide(hKMCs[2]);
  hLLbarMCs[0]->Scale(0.5);hLLbarMCs[1]->Scale(0.5);hLLbarMCs[2]->Scale(0.5);


  // Make X/pi ratios
  TH1D* hK0sToPi_Ref_Stat = (TH1D*)hK0s_Ref_Stat->Clone("hK0sToPi_Ref_Stat");
  TH1D* hK0sToPi_Jet_Stat = (TH1D*)hK0s_Jet_Stat->Clone("hK0sToPi_Jet_Stat");
  TH1D* hK0sToPi_Iso_Stat = (TH1D*)hK0s_Iso_Stat->Clone("hK0sToPi_Iso_Stat");
  TH1D* hLLbarToPi_Ref_Stat = (TH1D*)hL_Ref_Stat->Clone("hLLbarToPi_Ref_Stat");
  TH1D* hLLbarToPi_Jet_Stat = (TH1D*)hL_Jet_Stat->Clone("hLLbarToPi_Jet_Stat");
  TH1D* hLLbarToPi_Iso_Stat = (TH1D*)hL_Iso_Stat->Clone("hLLbarToPi_Iso_Stat");

  TH1D* hXToPi_Stat[6] = {
    hK0sToPi_Ref_Stat, hK0sToPi_Jet_Stat, hK0sToPi_Iso_Stat,
    hLLbarToPi_Ref_Stat, hLLbarToPi_Jet_Stat, hLLbarToPi_Iso_Stat  };

  TH1D* hK0sToPi_Ref_Syst = (TH1D*)hK0s_Ref_Syst->Clone("hK0sToPi_Ref_Syst");
  TH1D* hK0sToPi_Jet_Syst = (TH1D*)hK0s_Jet_Syst->Clone("hK0sToPi_Jet_Syst");
  TH1D* hK0sToPi_Iso_Syst = (TH1D*)hK0s_Iso_Syst->Clone("hK0sToPi_Iso_Syst");
  TH1D* hLLbarToPi_Ref_Syst = (TH1D*)hL_Ref_Syst->Clone("hLLbarToPi_Ref_Syst");
  TH1D* hLLbarToPi_Jet_Syst = (TH1D*)hL_Jet_Syst->Clone("hLLbarToPi_Jet_Syst");
  TH1D* hLLbarToPi_Iso_Syst = (TH1D*)hL_Iso_Syst->Clone("hLLbarToPi_Iso_Syst");
  TH1D* hK0sToPi_Ref_SystUnc = (TH1D*)hK0s_Ref_SystUnc->Clone("hK0sToPi_Ref_SystUnc");
  TH1D* hK0sToPi_Jet_SystUnc = (TH1D*)hK0s_Jet_SystUnc->Clone("hK0sToPi_Jet_SystUnc");
  TH1D* hK0sToPi_Iso_SystUnc = (TH1D*)hK0s_Iso_SystUnc->Clone("hK0sToPi_Iso_SystUnc");
  TH1D* hLLbarToPi_Ref_SystUnc = (TH1D*)hL_Ref_SystUnc->Clone("hLLbarToPi_Ref_SystUnc");
  TH1D* hLLbarToPi_Jet_SystUnc = (TH1D*)hL_Jet_SystUnc->Clone("hLLbarToPi_Jet_SystUnc");
  TH1D* hLLbarToPi_Iso_SystUnc = (TH1D*)hL_Iso_SystUnc->Clone("hLLbarToPi_Iso_SystUnc");

  TH1D* hXToPi_Syst[6] = {
    hK0sToPi_Ref_Syst, hK0sToPi_Jet_Syst, hK0sToPi_Iso_Syst,
    hLLbarToPi_Ref_Syst, hLLbarToPi_Jet_Syst, hLLbarToPi_Iso_Syst  };  

  TH1D* hXToPi_SystUnc[6] = {
    hK0sToPi_Ref_SystUnc, hK0sToPi_Jet_SystUnc, hK0sToPi_Iso_SystUnc,
    hLLbarToPi_Ref_SystUnc, hLLbarToPi_Jet_SystUnc, hLLbarToPi_Iso_SystUnc  };  

  for (int iH = 0; iH < 6; ++iH) {
    
    TH1D* hPion_Stat = (TH1D*)dirPion[iH%3]->Get("hPion_stat");
    TH1D* hPion_Syst = (TH1D*)dirPion[iH%3]->Get("hPion_tot_sys");
    TH1D* hPion_SystUnc = (TH1D*)dirPion[iH%3]->Get("hPion_unc_sys_am");

    TH1D* hPionR_Stat = (TH1D*)hPion_Stat->Rebin(NPTBINS,hPion_Stat->GetName(),XBINS);
    TH1D* hPionR_Syst = (TH1D*)hPion_Syst->Rebin(NPTBINS,hPion_Syst->GetName(),XBINS);
    TH1D* hPionR_SystUnc = (TH1D*)hPion_SystUnc->Rebin(NPTBINS,hPion_SystUnc->GetName(),XBINS);
    for (int iB = 0; iB < NPTBINS +1; iB++) {
      int divF = 0;
      double errStat = 0;
      double errSyst = 0;
      double errSystUnc = 0;
      for (int iPB = 0; iPB < NPIONBINS+1; iPB++) {
        if (XBINS[iB] == PIONBINS[iPB]) {
          for (int iPB2 = iPB; iPB2 < NPIONBINS+1; iPB2++) {
            if (PIONBINS[iPB2] != XBINS[iB+1]) {
              divF++;
              errStat += hPion_Stat->GetBinError(1+iPB2)*hPion_Stat->GetBinError(1+iPB2);
              errSyst += hPion_Syst->GetBinError(1+iPB2);
              errSystUnc += hPion_SystUnc->GetBinError(1+iPB2);
            }
            else break;
          }
        }
      }
      
      hPionR_Stat->SetBinContent(1+iB,1./divF*hPionR_Stat->GetBinContent(1+iB));
      hPionR_Stat->SetBinError(1+iB,1./divF*TMath::Sqrt(errStat));
      hPionR_Syst->SetBinContent(1+iB,1./divF*hPionR_Syst->GetBinContent(1+iB));
      hPionR_Syst->SetBinError(1+iB,1./divF*errSyst);
      hPionR_SystUnc->SetBinContent(1+iB,1./divF*hPionR_SystUnc->GetBinContent(1+iB));
      hPionR_SystUnc->SetBinError(1+iB,1./divF*errSystUnc);
      
    }

    hXToPi_Stat[iH]->Divide(hPionR_Stat);
    hXToPi_Syst[iH]->Divide(hPionR_Syst);
    hXToPi_SystUnc[iH]->Divide(hPionR_SystUnc);
    //hXToPi_Stat[iH]->Write(); //hXToPi_Syst[iH]->Write();
  }



  // Produce plots -- spectra
  double rl = 0.; double rh = 9.;
  TH1D* hBlank = new TH1D("hBlank",";p_{T} (GeV/#it{c});K^{0}_{S}",10,rl,rh);
  MakeNiceHistogram(hBlank,kBlack);
  hK0s_Jet_Syst->SetTitle(";p_{T} (GeV/#it{c});");
  for (int i=1; i<11; i++) hBlank->SetBinContent(i,1000);

  for (int iH = 0; iH < 9; iH++) {
    hSysts[iH]->SetTitle("");
    hSysts[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hSysts[iH],colors[iH%3]);
    MakeNiceHistogram(hStats[iH],colors[iH%3]);
    hSysts[iH]->GetYaxis()->SetRangeUser(49e-5,100.*hSysts[iH]->GetMaximum());
    hSystsUnc[iH]->SetTitle("");
    hSystsUnc[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hSystsUnc[iH],colors[iH%3]);
    hStats[iH]->Write();
    hSysts[iH]->Write();
    hSystsUnc[iH]->Write();
  }

  TCanvas* cSpK0s = new TCanvas("cSpK0s","",800,1000);
  cSpK0s->SetLogy();

  hBlank->GetYaxis()->SetRangeUser(5e-4,5.01); hBlank->SetLineWidth(0);
  hBlank->Draw();
  hK0s_Ref_Syst->Draw("e2 same");
  hK0s_Jet_Syst->Draw("e2 same");
  hK0s_Iso_Syst->Draw("e2 same");
  hK0s_Ref_Stat->Draw("E X0 same");
  hK0s_Jet_Stat->Draw("E X0 same");
  hK0s_Iso_Stat->Draw("E X0 same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,0.8,Form("%s (N_{ch} #geq 10)",multPlot.Data()));
  latex.SetTextSize(0.030); latex.DrawLatex(0.5,0.76,"pp #sqrt{s} = 13 TeV, |y|<0.8");
  latex.DrawLatex(0.5,0.70,"Reference");
  latex.SetTextColor(colors[1]); latex.DrawLatex(0.5,0.67,Form("Jetty 0-%s%",spher.Data()));
  latex.SetTextColor(colors[2]); latex.DrawLatex(0.5,0.64,Form("Isotropic 0-%s%",spher.Data()));
  MakeRatioPlot(hBlank,hK0s_Ref_Syst,cSpK0s,0.1,2.1,rl,rh,"");
  MakeRatioPlot(hK0s_Jet_SystUnc,hK0s_Ref_SystUnc,cSpK0s,0.1,2.1,rl,rh,"e2");
  MakeRatioPlot(hK0s_Iso_SystUnc,hK0s_Ref_SystUnc,cSpK0s,0.1,2.1,rl,rh,"e2");
  MakeRatioPlot(hK0s_Jet_Stat,hK0s_Ref_Stat,cSpK0s,0.1,2.1,rl,rh,"e x0");
  MakeRatioPlot(hK0s_Iso_Stat,hK0s_Ref_Stat,cSpK0s,0.1,2.1,rl,rh,"e x0");


  for (int iH = 0; iH < 3; iH++) {
    hLLbar_Syst[iH]->SetTitle("");
    hLLbar_Syst[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hLLbar_Syst[iH],colors[iH%3]);
    MakeNiceHistogram(hLLbar_Stat[iH],colors[iH%3]);
    hLLbar_Syst[iH]->GetYaxis()->SetRangeUser(49e-5,100.*hLLbar_Syst[iH]->GetMaximum());
    hLLbar_SystUnc[iH]->SetTitle("");
    hLLbar_SystUnc[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hLLbar_SystUnc[iH],colors[iH%3]);
    hLLbar_Stat[iH]->Write();
    hLLbar_Syst[iH]->Write();
    hLLbar_SystUnc[iH]->Write();
  }

  TCanvas* cSpLLbar = new TCanvas("cSpLLbar","",800,1000);
  cSpLLbar->SetLogy();

  hBlank->GetYaxis()->SetTitle("#Lambda + #bar{#Lambda}");
  hBlank->GetYaxis()->SetRangeUser(3e-4,5.01);
  hBlank->Draw();
  hLLbar_Ref_Syst->Draw("e2 same");
  hLLbar_Jet_Syst->Draw("e2 same");
  hLLbar_Iso_Syst->Draw("e2 same");
  hLLbar_Ref_Stat->Draw("E X0 same");
  hLLbar_Jet_Stat->Draw("E X0 same");
  hLLbar_Iso_Stat->Draw("E X0 same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,0.8,Form("%s (N_{ch} #geq 10)",multPlot.Data()));
  latex.SetTextSize(0.030); latex.DrawLatex(0.5,0.76,"pp #sqrt{s} = 13 TeV, |y|<0.8");
  latex.DrawLatex(0.5,0.70,"Reference");
  latex.SetTextColor(colors[1]); latex.DrawLatex(0.5,0.67,Form("Jetty 0-%s%",spher.Data()));
  latex.SetTextColor(colors[2]); latex.DrawLatex(0.5,0.64,Form("Isotropic 0-%s%",spher.Data()));
  MakeRatioPlot(hBlank,hLLbar_Ref_Syst,cSpLLbar,0.1,2.1,rl,rh,"");
  MakeRatioPlot(hLLbar_Jet_SystUnc,hLLbar_Ref_SystUnc,cSpLLbar,0.1,2.1,rl,rh,"e2");
  MakeRatioPlot(hLLbar_Iso_SystUnc,hLLbar_Ref_SystUnc,cSpLLbar,0.1,2.1,rl,rh,"e2");
  MakeRatioPlot(hLLbar_Jet_Stat,hLLbar_Ref_Stat,cSpLLbar,0.1,2.1,rl,rh,"e x0");
  MakeRatioPlot(hLLbar_Iso_Stat,hLLbar_Ref_Stat,cSpLLbar,0.1,2.1,rl,rh,"e x0");


  // Produce plots -- L/K0s
  hLLbarToK0s_Jet_Syst->SetTitle(";p_{T} (GeV/#it{c});");
  
  for (int iH = 0; iH < 3; iH++) {
    hKL_Syst[iH]->SetTitle("");
    hKL_Syst[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hKL_Syst[iH],colors[iH%3]);
    MakeNiceHistogram(hKL_Stat[iH],colors[iH%3]);
    hKL_Syst[iH]->GetYaxis()->SetRangeUser(0,2.*hKL_Syst[iH]->GetMaximum());
    hKL_SystUnc[iH]->SetTitle("");
    hKL_SystUnc[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hKL_SystUnc[iH],colors[iH%3]);
    hKL_Stat[iH]->Write();
    hKL_Syst[iH]->Write();
  }

  TCanvas* cLLbarToK0s = new TCanvas("cLLbarToK0s","",800,1000);

  hBlank->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / 2K_{S}^{0}");
  hBlank->GetYaxis()->SetRangeUser(0,1.3); hBlank->SetLineWidth(0);
  hBlank->Draw();
  hLLbarToK0s_Ref_Syst->Draw("e2 same");
  hLLbarToK0s_Jet_Syst->Draw("e2 same");
  hLLbarToK0s_Iso_Syst->Draw("e2 same");
  hLLbarToK0s_Ref_Stat->Draw("E X0 same");
  hLLbarToK0s_Jet_Stat->Draw("E X0 same");
  hLLbarToK0s_Iso_Stat->Draw("E X0 same");
  if (multInt==3&&sphInt==0) hLLbar_Ref_MC->Draw("hist l same");hLLbar_Jet_MC->Draw("hist l same");hLLbar_Iso_MC->Draw("hist l same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,0.8,Form("%s (N_{ch} #geq 10)",multPlot.Data()));
  latex.SetTextSize(0.030); latex.DrawLatex(0.5,0.76,"pp #sqrt{s} = 13 TeV, |y|<0.8");
  latex.DrawLatex(0.5,0.70,"Reference");
  latex.SetTextColor(colors[1]); latex.DrawLatex(0.5,0.67,Form("Jetty 0-%s%",spher.Data()));
  latex.SetTextColor(colors[2]); latex.DrawLatex(0.5,0.64,Form("Isotropic 0-%s%",spher.Data()));
  MakeRatioPlot(hBlank,hLLbarToK0s_Ref_Syst,cLLbarToK0s,0.4,1.6,rl,rh,"");
  MakeRatioPlot(hLLbarToK0s_Jet_SystUnc,hLLbarToK0s_Ref_SystUnc,cLLbarToK0s,0.4,1.6,rl,rh,"e2");
  MakeRatioPlot(hLLbarToK0s_Iso_SystUnc,hLLbarToK0s_Ref_SystUnc,cLLbarToK0s,0.4,1.6,rl,rh,"e2");
  MakeRatioPlot(hLLbarToK0s_Jet_Stat,hLLbarToK0s_Ref_Stat,cLLbarToK0s,0.4,1.6,rl,rh,"e x0");
  MakeRatioPlot(hLLbarToK0s_Iso_Stat,hLLbarToK0s_Ref_Stat,cLLbarToK0s,0.4,1.6,rl,rh,"e x0");
  if (multInt==3&&sphInt==0) MakeRatioPlot(hLLbar_Jet_MC,hLLbar_Ref_MC,cLLbarToK0s,0.5,1.5,rl,rh,"hist l");MakeRatioPlot(hLLbar_Iso_MC,hLLbar_Ref_MC,cLLbarToK0s,0.5,1.5,rl,rh,"hist l");

  // Produce plots -- X/pi ratios
  hK0sToPi_Jet_Syst->SetTitle(";p_{T} (GeV/#it{c});");
  
  for (int iH = 0; iH < 6; iH++) {
    hXToPi_Syst[iH]->SetTitle("");
    hXToPi_Syst[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hXToPi_Syst[iH],colors[iH%3]);
    MakeNiceHistogram(hXToPi_Stat[iH],colors[iH%3]);
    hXToPi_Syst[iH]->GetYaxis()->SetRangeUser(0,2.*hXToPi_Syst[iH]->GetMaximum());
    hXToPi_SystUnc[iH]->SetTitle("");
    hXToPi_SystUnc[iH]->SetFillColor(kWhite);
    MakeNiceHistogram(hXToPi_SystUnc[iH],colors[iH%3]);
    hXToPi_Stat[iH]->Write();
    hXToPi_Syst[iH]->Write();
  }

  TCanvas* cK0sToPi = new TCanvas("cK0sToPi","",800,1000);

  hBlank->GetYaxis()->SetTitle("K_{S}^{0} / #pi");
  //hBlank->GetYaxis()->SetTitle("K_{S}^{0} / K^{#pm}");
  hBlank->GetYaxis()->SetTitleOffset(1.51);
  hBlank->GetYaxis()->SetRangeUser(0,0.5); hBlank->SetLineWidth(0);
  //hBlank->GetYaxis()->SetRangeUser(0.2,0.9); hBlank->SetLineWidth(0);
  hBlank->Draw();
  hK0sToPi_Ref_Syst->Draw("e2 same");
  hK0sToPi_Jet_Syst->Draw("e2 same");
  hK0sToPi_Iso_Syst->Draw("e2 same");
  hK0sToPi_Ref_Stat->Draw("E X0 same");
  hK0sToPi_Jet_Stat->Draw("E X0 same");
  hK0sToPi_Iso_Stat->Draw("E X0 same");
  //hK0s_Ref_MC->Draw("hist l same");hK0s_Jet_MC->Draw("hist l same");hK0s_Iso_MC->Draw("hist l same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,0.8,Form("%s (N_{ch} #geq 10)",multPlot.Data()));
  latex.SetTextSize(0.030); latex.DrawLatex(0.5,0.76,"pp #sqrt{s} = 13 TeV, |y|<0.8");
  latex.DrawLatex(0.5,0.70,"Reference");
  latex.SetTextColor(colors[1]); latex.DrawLatex(0.5,0.67,Form("Jetty 0-%s%",spher.Data()));
  latex.SetTextColor(colors[2]); latex.DrawLatex(0.5,0.64,Form("Isotropic 0-%s%",spher.Data()));
  MakeRatioPlot(hBlank,hK0sToPi_Ref_Syst,cK0sToPi,0.5,1.5,rl,rh,"");
  MakeRatioPlot(hK0sToPi_Jet_SystUnc,hK0sToPi_Ref_SystUnc,cK0sToPi,0.5,1.5,rl,rh,"e2");
  MakeRatioPlot(hK0sToPi_Iso_SystUnc,hK0sToPi_Ref_SystUnc,cK0sToPi,0.5,1.5,rl,rh,"e2");
  MakeRatioPlot(hK0sToPi_Jet_Stat,hK0sToPi_Ref_Stat,cK0sToPi,0.5,1.5,rl,rh,"e x0");
  MakeRatioPlot(hK0sToPi_Iso_Stat,hK0sToPi_Ref_Stat,cK0sToPi,0.5,1.5,rl,rh,"e x0");
  //MakeRatioPlot(hK0s_Jet_MC,hK0s_Ref_MC,cK0sToPi,0.5,1.5,rl,rh,"hist l");MakeRatioPlot(hK0s_Iso_MC,hK0s_Ref_MC,cK0sToPi,0.5,1.5,rl,rh,"hist l");

  

  hLLbarToPi_Jet_Syst->SetTitle(";p_{T} (GeV/#it{c});");
  TCanvas* cLLbarToPi = new TCanvas("cLLbarToPi","",800,1000);

  hBlank->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / #pi");
  //hBlank->GetYaxis()->SetTitle("(#Lambda + #bar{#Lambda}) / p");
  hBlank->GetYaxis()->SetTitleOffset(1.51);
  hBlank->GetYaxis()->SetRangeUser(0,0.301); hBlank->SetLineWidth(0);
  //hBlank->GetYaxis()->SetRangeUser(0.00,0.8); hBlank->SetLineWidth(0);
  hBlank->Draw();
  hLLbarToPi_Ref_Syst->Draw("e2 same");
  hLLbarToPi_Jet_Syst->Draw("e2 same");
  hLLbarToPi_Iso_Syst->Draw("e2 same");
  hLLbarToPi_Ref_Stat->Draw("E X0 same");
  hLLbarToPi_Jet_Stat->Draw("E X0 same");
  hLLbarToPi_Iso_Stat->Draw("E X0 same");

  hLLbarToPi_Ref_MC->SetLineWidth(2); hLLbarToPi_Jet_MC->SetLineWidth(2); hLLbarToPi_Iso_MC->SetLineWidth(2);
  hLLbarToPi_Ref_MC->Rebin(4); hLLbarToPi_Jet_MC->Rebin(4); hLLbarToPi_Iso_MC->Rebin(4);
  hLLbarToPi_Ref_MC->Scale(1./4); hLLbarToPi_Jet_MC->Scale(1./4); hLLbarToPi_Iso_MC->Scale(1./4);
  //hLLbarToPi_Ref_MC->Draw("hist l same"); hLLbarToPi_Jet_MC->Draw("hist l same"); hLLbarToPi_Iso_MC->Draw("hist l same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,0.8,Form("%s (N_{ch} #geq 10)",multPlot.Data()));
  latex.SetTextSize(0.030); latex.DrawLatex(0.5,0.76,"pp #sqrt{s} = 13 TeV, |y|<0.8");
  latex.DrawLatex(0.5,0.70,"Reference");
  latex.SetTextColor(colors[1]); latex.DrawLatex(0.5,0.67,Form("Jetty 0-%s%",spher.Data()));
  latex.SetTextColor(colors[2]); latex.DrawLatex(0.5,0.64,Form("Isotropic 0-%s%",spher.Data()));
  MakeRatioPlot(hBlank,hLLbarToPi_Ref_Syst,cLLbarToPi,0.5,1.5,rl,rh,"");
  MakeRatioPlot(hLLbarToPi_Jet_SystUnc,hLLbarToPi_Ref_SystUnc,cLLbarToPi,0.5,1.5,rl,rh,"e2");
  MakeRatioPlot(hLLbarToPi_Iso_SystUnc,hLLbarToPi_Ref_SystUnc,cLLbarToPi,0.5,1.5,rl,rh,"e2");
  MakeRatioPlot(hLLbarToPi_Jet_Stat,hLLbarToPi_Ref_Stat,cLLbarToPi,0.5,1.5,rl,rh,"e x0");
  MakeRatioPlot(hLLbarToPi_Iso_Stat,hLLbarToPi_Ref_Stat,cLLbarToPi,0.5,1.5,rl,rh,"e x0");

  //MakeRatioPlot(hLLbarToPi_Iso_MC,hLLbarToPi_Ref_MC,cLLbarToPi,0.5,1.5,rl,rh,"hist l"); MakeRatioPlot(hLLbarToPi_Jet_MC,hLLbarToPi_Ref_MC,cLLbarToPi,0.5,1.5,rl,rh,"hist l");


  cSpK0s->SaveAs(Form("sp_K0s_%s_spher%s.pdf",mult.Data(),spher.Data()));
  cSpLLbar->SaveAs(Form("sp_LLbar_%s_spher%s.pdf",mult.Data(),spher.Data()));
  
  cLLbarToK0s->SaveAs(Form("LLbarToK0s_%s_spher%s.pdf",mult.Data(),spher.Data()));
  
  cK0sToPi->SaveAs(Form("K0sToPi_%s_spher%s.pdf",mult.Data(),spher.Data()));
  cLLbarToPi->SaveAs(Form("LLbarToPi%s_spher%s.pdf",mult.Data(),spher.Data()));


//-------------------------
  TIter objIt(fout->GetList(),kIterForward);
  TObject* obj = 0;
  while ( (obj = objIt()) ) {
    TString objName(obj->GetName());
    if (objName.BeginsWith("hr_hBlank")) {      fout->Remove(obj);  }
      
  }


}

