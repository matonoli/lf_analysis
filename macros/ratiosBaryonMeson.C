#include <iostream>
#include "MakeCanvas.C"

using namespace std;

//// UNCOMMENT A LINE TO CHOOSE THE MULTIPLICITY SELECTION
//TString mult("V0M01"); int multInt = 0; int multOmar = 0; TString multPlot("V0M 0-1%");
TString mult("NCharged01"); int multInt = 1; int multOmar = 0; TString multPlot("CL1 0-1%");
//TString mult("V0M"); int multInt = 2; int multOmar = 1; TString multPlot("V0M 0-10%");
//TString mult("NCharged"); int multInt = 3; int multOmar = 1; TString multPlot("CL1 0-10%");

//// UNCOMMENT A LINE TO CHOOSE THE SPHEROCITY SELECTION
//TString spher("20"); int sphInt = 3;
TString spher("10"); int sphInt = 2;
//TString spher("5"); int sphInt = 1;
//TString spher("1"); int sphInt = 0;

// COMMON PATH FOR INPUTS
TString path("../files/paper");

// HELPER GLOBAL VARIABLES
enum {CL1, CL101, V0M, V0M01, NMULTI};
enum {Ref, Jetty, Iso, NSPHERO};
enum {Stat, Syst, SystUnc, Monash, Ropes, NHIST};
const char* strS[NSPHERO] = {"Ref", "Jetty", "Iso"};
const char* strH[NHIST] = {"Stat", "Syst", "SystUnc", "Monash", "Ropes"};
Int_t colors[3] = {kBlack, kRed, kBlue};

// PLOTTING FUNCTIONS

//______________________________________________________________
void myLegendSetUp(TLegend *currentLegend,float currentTextSize){
    
    currentLegend->SetBorderSize(0);
    currentLegend->SetFillStyle(0);
    currentLegend->SetFillColor(0);
    currentLegend->SetMargin(0.25);
    currentLegend->SetTextSize(currentTextSize);
//  currentLegend->SetEntrySeparation(0.35);
    return;
}

void TuneLatex(TLatex* latex, Float_t textSize){

  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextAngle(0);
  latex->SetTextFont(43);
  latex->SetTextSize(textSize);

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
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);
  h->SetMarkerColor(col);
  h->SetStats(0);
  h->SetFillStyle(0);

  //h->GetYaxis()->SetTitleSize(30);
  //h->GetYaxis()->SetTitleFont(43);
  //h->GetYaxis()->SetTitleOffset(1.0);
  //h->GetYaxis()->SetLabelOffset(0.0025);
  //h->GetYaxis()->SetLabelSize(20);
  //h->GetYaxis()->SetLabelFont(43);
  //h->GetXaxis()->SetLabelFont(43);

  //h->SetTopMargin(0.055);
}

TH1D* DivideSpline(TH1D* hn, TH1D* hd, Double_t scale) {

  // replaces the denominator with a spline and creates a clone with the numerator's binning
  // returns the ratio
  

  Int_t hd_nbins = hd->GetNbinsX();
  
  Double_t hd_x[400], hd_err[400]; 
  for (int iB = 0; iB < hd_nbins; ++iB)   {
    hd_x[iB]  = hd->GetBinCenter(iB+1);
    hd_err[iB]  = hd->GetBinError(iB+1);    }
  

  TSpline3* hd_spl    = new TSpline3(hd, 0, 1, 0);
  TSpline3* hd_splerr = new TSpline3(Form("splerr_%s",hd->GetName()), hd_x, hd_err, hd_nbins, 0, 1, 0);

  TH1D* hd2 = (TH1D*)hn->Clone(Form("hr_%s",hn->GetName()));
  for (int iB = 1; iB <= hd2->GetNbinsX(); ++iB)         {
    hd2->SetBinContent(iB, hd_spl->Eval(hd2->GetBinCenter(iB)));  
    hd2->SetBinError(iB, hd_splerr->Eval(hd2->GetBinCenter(iB))); }

  hn->Divide(hn,hd2,scale);
  delete hd_spl;
  delete hd_splerr;
  delete hd2;
  
  return hn;
}

void DrawHistograms(TH1D** h, Int_t nhist, Int_t col) {

  Int_t lineMonash = 1;
  Int_t lineRopes = 2;
  Int_t lineWidth = 2;
  Int_t rebinF = 10;
  h[Monash]->Rebin(rebinF); h[Monash]->Scale(1./rebinF);
  h[Ropes]->Rebin(rebinF); h[Ropes]->Scale(1./rebinF);

  for (int iH = 0; iH < nhist; ++iH) {
    
    MakeNiceHistogram(h[iH],col);
  }

  h[Monash]->SetLineStyle(lineMonash);
  h[Ropes]->SetLineStyle(lineRopes);
  h[Monash]->SetLineWidth(lineWidth);
  h[Ropes]->SetLineWidth(lineWidth);

  h[Monash]->GetXaxis()->SetRange(0.,1.+h[Stat]->GetBinCenter(h[Stat]->GetNbinsX()));
  h[Ropes]->GetXaxis()->SetRange(0.,1.+h[Stat]->GetBinCenter(h[Stat]->GetNbinsX()));

  h[Stat]->Draw("E X0 same");
  h[Syst]->Draw("e2 same");
  h[Monash]->Draw("hist l same");
  h[Ropes]->Draw("hist l same");

}

TH1D* MakeRatioHist(TH1D* hn, TH1D* hd) {

  TH1D* hr = (TH1D*)hn->Clone(Form("hr_%s",hn->GetName()));
  hr->Divide(hd);

  return hr;
}

void DrawHistogramsRatios(TH1D** hn, TH1D** hd, Int_t nhist, Int_t col) {

  TH1D* hr[NHIST];
  for (int iH = 0; iH < nhist; ++iH) {
    
    Int_t rebinF = hd[iH]->GetNbinsX() / hn[iH]->GetNbinsX();
    if (iH>SystUnc) {
      hd[iH]->Rebin(rebinF);
      hd[iH]->Scale(1./rebinF); }
    hr[iH] = MakeRatioHist(hn[iH],hd[iH]);
  }

  hr[Stat]->Draw("E X0 same");
  hr[SystUnc]->Draw("e2 same");
  hr[Monash]->Draw("hist l same");
  hr[Ropes]->Draw("hist l same");

}




void ratiosBaryonMeson() {

  TH1::SetDefaultSumw2(1);

  // LOAD FILES AND HISTOGRAMS
  TH1D* hK0s[NSPHERO][NHIST];
  TH1D* hL[NSPHERO][NHIST];
  TH1D* hP[NSPHERO][NHIST];
  TH1D* hPi[NSPHERO][NHIST];
  TH1D* hXi[NSPHERO][NHIST];
  TH1D* hPhi[NSPHERO][NHIST];

  TH1D* hPoverPi[NSPHERO][NHIST];
  TH1D* hLoverK0s[NSPHERO][NHIST];
  TH1D* hXioverPhi[NSPHERO][NHIST];

  // pions
  TFile* fPion = new TFile(Form("%s/pi_k_p_spectra.root",path.Data()),"READ");
  const char* strM[4] = {"V0M", "Trks", "V0M", "Trks"};
  hPi[Ref][Stat]      = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hPion_stat",strM[multInt],multOmar));
  hPi[Jetty][Stat]    = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hPion_stat",strM[multInt],sphInt,multOmar));
  hPi[Iso][Stat]      = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hPion_stat",strM[multInt],sphInt,multOmar));
  hPi[Ref][Syst]      = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hPion_tot_sys",strM[multInt],multOmar));
  hPi[Jetty][Syst]    = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hPion_tot_sys",strM[multInt],sphInt,multOmar));
  hPi[Iso][Syst]      = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hPion_tot_sys",strM[multInt],sphInt,multOmar));
  hPi[Ref][SystUnc]   = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hPion_unc_sys_am",strM[multInt],multOmar));
  hPi[Jetty][SystUnc] = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hPion_unc_sys_am",strM[multInt],sphInt,multOmar));
  hPi[Iso][SystUnc]   = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hPion_unc_sys_am",strM[multInt],sphInt,multOmar));

  // protons
  hP[Ref][Stat]      = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hProton_stat",strM[multInt],multOmar));
  hP[Jetty][Stat]    = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hProton_stat",strM[multInt],sphInt,multOmar));
  hP[Iso][Stat]      = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hProton_stat",strM[multInt],sphInt,multOmar));
  hP[Ref][Syst]      = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hProton_tot_sys",strM[multInt],multOmar));
  hP[Jetty][Syst]    = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hProton_tot_sys",strM[multInt],sphInt,multOmar));
  hP[Iso][Syst]      = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hProton_tot_sys",strM[multInt],sphInt,multOmar));
  hP[Ref][SystUnc]   = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hProton_unc_sys_am",strM[multInt],multOmar));
  hP[Jetty][SystUnc] = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hProton_unc_sys_am",strM[multInt],sphInt,multOmar));
  hP[Iso][SystUnc]   = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hProton_unc_sys_am",strM[multInt],sphInt,multOmar));

  // k0s
  TFile* fK0s = new TFile(Form("%s/k0s_l_spectra_%s_spher%s.root",path.Data(),mult.Data(),spher.Data()),"READ");
  hK0s[Ref][Stat]       = (TH1D*)fK0s->Get("hK0s_Ref_Stat");
  hK0s[Jetty][Stat]     = (TH1D*)fK0s->Get("hK0s_Jet_Stat");
  hK0s[Iso][Stat]       = (TH1D*)fK0s->Get("hK0s_Iso_Stat");
  hK0s[Ref][Syst]       = (TH1D*)fK0s->Get("hK0s_Ref_Syst");
  hK0s[Jetty][Syst]     = (TH1D*)fK0s->Get("hK0s_Jet_Syst");
  hK0s[Iso][Syst]       = (TH1D*)fK0s->Get("hK0s_Iso_Syst");
  hK0s[Ref][SystUnc]    = (TH1D*)fK0s->Get("hK0s_Ref_SystUnc");    // needs changing to unc.
  hK0s[Jetty][SystUnc]  = (TH1D*)fK0s->Get("hK0s_Jet_SystUnc");
  hK0s[Iso][SystUnc]    = (TH1D*)fK0s->Get("hK0s_Iso_SystUnc");

  // lambdas
  hL[Ref][Stat]       = (TH1D*)fK0s->Get("hLLbar_Ref_Stat");
  hL[Jetty][Stat]     = (TH1D*)fK0s->Get("hLLbar_Jet_Stat");
  hL[Iso][Stat]       = (TH1D*)fK0s->Get("hLLbar_Iso_Stat");
  hL[Ref][Syst]       = (TH1D*)fK0s->Get("hLLbar_Ref_Syst");
  hL[Jetty][Syst]     = (TH1D*)fK0s->Get("hLLbar_Jet_Syst");
  hL[Iso][Syst]       = (TH1D*)fK0s->Get("hLLbar_Iso_Syst");
  hL[Ref][SystUnc]    = (TH1D*)fK0s->Get("hLLbar_Ref_SystUnc");    // needs changing to unc.
  hL[Jetty][SystUnc]  = (TH1D*)fK0s->Get("hLLbar_Jet_SystUnc");
  hL[Iso][SystUnc]    = (TH1D*)fK0s->Get("hLLbar_Iso_SystUnc");

  // l to k0s we can get directly
  /*hLoverK0s[Ref][Stat]      = (TH1D*)fK0s->Get("hLLbarToK0s_Ref_Stat");
  hLoverK0s[Jetty][Stat]    = (TH1D*)fK0s->Get("hLLbarToK0s_Jet_Stat");
  hLoverK0s[Iso][Stat]      = (TH1D*)fK0s->Get("hLLbarToK0s_Iso_Stat");
  hLoverK0s[Ref][Syst]      = (TH1D*)fK0s->Get("hLLbarToK0s_Ref_Syst");
  hLoverK0s[Jetty][Syst]    = (TH1D*)fK0s->Get("hLLbarToK0s_Jet_Syst");
  hLoverK0s[Iso][Syst]      = (TH1D*)fK0s->Get("hLLbarToK0s_Iso_Syst");
  hLoverK0s[Ref][SystUnc]   = (TH1D*)fK0s->Get("hLLbarToK0s_Ref_Syst"); // needs changing to unc.
  hLoverK0s[Jetty][SystUnc] = (TH1D*)fK0s->Get("hLLbarToK0s_Jet_Syst");
  hLoverK0s[Iso][SystUnc]   = (TH1D*)fK0s->Get("hLLbarToK0s_Iso_Syst");*/
  
  // phi
  TFile* fPhi = new TFile(Form("%s/CorrectedCL1Spectra.root",path.Data()),"READ");
  hPhi[Ref][Stat]       = (TH1D*)fPhi->Get("hPhi_HM_Spectra_stat");
  hPhi[Jetty][Stat]     = (TH1D*)fPhi->Get("hPhi_Jetty_Spectra_stat");
  hPhi[Iso][Stat]       = (TH1D*)fPhi->Get("hPhi_Iso_Spectra_stat");
  hPhi[Ref][Syst]       = (TH1D*)fPhi->Get("hPhi_HM_Spectra_syst");
  hPhi[Jetty][Syst]     = (TH1D*)fPhi->Get("hPhi_Jetty_Spectra_syst");
  hPhi[Iso][Syst]       = (TH1D*)fPhi->Get("hPhi_Iso_Spectra_syst");
  hPhi[Ref][SystUnc]    = (TH1D*)fPhi->Get("hPhi_HM_Spectra_syst");
  hPhi[Jetty][SystUnc]  = (TH1D*)fPhi->Get("hPhi_Jetty_Spectra_syst");
  hPhi[Iso][SystUnc]    = (TH1D*)fPhi->Get("hPhi_Iso_Spectra_syst");

  //xi
  TFile* fXi = new TFile(Form("%s/xi_results_So_CL1_top1_aug_6.root",path.Data()),"READ");
  hXi[Ref][Stat]        = (TH1D*)fXi->Get("hXi_Ref_Spectra_stat");
  hXi[Jetty][Stat]      = (TH1D*)fXi->Get(Form("hXi_Jetty%s_Spectra_stat",spher.Data()));
  hXi[Iso][Stat]        = (TH1D*)fXi->Get(Form("hXi_Iso%s_Spectra_stat",spher.Data()));
  hXi[Ref][Syst]        = (TH1D*)fXi->Get("hXi_Ref_Spectra_syst");
  hXi[Jetty][Syst]      = (TH1D*)fXi->Get(Form("hXi_Jetty%s_Spectra_syst",spher.Data()));
  hXi[Iso][Syst]        = (TH1D*)fXi->Get(Form("hXi_Iso%s_Spectra_syst",spher.Data()));
  hXi[Ref][SystUnc]     = (TH1D*)fXi->Get("hXi_Ref_Spectra_syst");
  hXi[Jetty][SystUnc]   = (TH1D*)fXi->Get(Form("hXi_Jetty%s_Spectra_syst",spher.Data()));
  hXi[Iso][SystUnc]     = (TH1D*)fXi->Get(Form("hXi_Iso%s_Spectra_syst",spher.Data()));

  //monash
  TFile* fMonash = new TFile(Form("%s/Complete_Monash.root",path.Data()),"READ");
  const char* strM_MC[4] = { "V0M","CL1","V0M","CL1"};
  const char* strS_MC[4] = { "01","05","10","20"};
  hPi[Ref][Monash]    = (TH1D*)fMonash->Get(Form("fPions_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hPi[Jetty][Monash]  = (TH1D*)fMonash->Get(Form("fPions_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hPi[Iso][Monash]    = (TH1D*)fMonash->Get(Form("fPions_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hP[Ref][Monash]     = (TH1D*)fMonash->Get(Form("fProtons_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hP[Jetty][Monash]   = (TH1D*)fMonash->Get(Form("fProtons_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hP[Iso][Monash]     = (TH1D*)fMonash->Get(Form("fProtons_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Ref][Monash]   = (TH1D*)fMonash->Get(Form("fK0s_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hK0s[Jetty][Monash] = (TH1D*)fMonash->Get(Form("fK0s_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Iso][Monash]   = (TH1D*)fMonash->Get(Form("fK0s_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hL[Ref][Monash]     = (TH1D*)fMonash->Get(Form("fLambdas_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hL[Jetty][Monash]   = (TH1D*)fMonash->Get(Form("fLambdas_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hL[Iso][Monash]     = (TH1D*)fMonash->Get(Form("fLambdas_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hPhi[Ref][Monash]   = (TH1D*)fMonash->Get(Form("fPhis_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hPhi[Jetty][Monash] = (TH1D*)fMonash->Get(Form("fPhis_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hPhi[Iso][Monash]   = (TH1D*)fMonash->Get(Form("fPhis_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hXi[Ref][Monash]    = (TH1D*)fMonash->Get(Form("fXis_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hXi[Jetty][Monash]  = (TH1D*)fMonash->Get(Form("fXis_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hXi[Iso][Monash]    = (TH1D*)fMonash->Get(Form("fXis_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));

  //ropes
  TFile* fRopes = new TFile(Form("%s/Complete_Ropes.root",path.Data()),"READ");
  hPi[Ref][Ropes]    = (TH1D*)fRopes->Get(Form("fPions_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hPi[Jetty][Ropes]  = (TH1D*)fRopes->Get(Form("fPions_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hPi[Iso][Ropes]    = (TH1D*)fRopes->Get(Form("fPions_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hP[Ref][Ropes]     = (TH1D*)fRopes->Get(Form("fProtons_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hP[Jetty][Ropes]   = (TH1D*)fRopes->Get(Form("fProtons_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hP[Iso][Ropes]     = (TH1D*)fRopes->Get(Form("fProtons_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Ref][Ropes]   = (TH1D*)fRopes->Get(Form("fK0s_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hK0s[Jetty][Ropes] = (TH1D*)fRopes->Get(Form("fK0s_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Iso][Ropes]   = (TH1D*)fRopes->Get(Form("fK0s_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hL[Ref][Ropes]     = (TH1D*)fRopes->Get(Form("fLambdas_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hL[Jetty][Ropes]   = (TH1D*)fRopes->Get(Form("fLambdas_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hL[Iso][Ropes]     = (TH1D*)fRopes->Get(Form("fLambdas_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hPhi[Ref][Ropes]   = (TH1D*)fRopes->Get(Form("fPhis_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hPhi[Jetty][Ropes] = (TH1D*)fRopes->Get(Form("fPhis_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hPhi[Iso][Ropes]   = (TH1D*)fRopes->Get(Form("fPhis_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hXi[Ref][Ropes]    = (TH1D*)fRopes->Get(Form("fXis_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hXi[Jetty][Ropes]  = (TH1D*)fRopes->Get(Form("fXis_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hXi[Iso][Ropes]    = (TH1D*)fRopes->Get(Form("fXis_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  //fRopes->Close();

  TFile* fout = new TFile(Form("BMratios_%s_spher%s.root",mult.Data(),spher.Data()),"RECREATE");


  // CALCULATE RATIOS
  for (int iS = 0; iS < NSPHERO; ++iS)  {
  for (int iH = 0; iH < NHIST; ++iH)  {

    hPoverPi[iS][iH] = (TH1D*)hP[iS][iH]->Clone(Form("hPoverPi_%s_%s",strS[iS],strH[iH]));
    hLoverK0s[iS][iH] = (TH1D*)hL[iS][iH]->Clone(Form("hLoverK0s_%s_%s",strS[iS],strH[iH]));
    //hLoverK0s[iS][iH]->SetTitle(Form("hLoverK0s_%s_%s",strS[iS],strH[iH]));
    hXioverPhi[iS][iH] = (TH1D*)hXi[iS][iH]->Clone(Form("hXioverPhi_%s_%s",strS[iS],strH[iH]));
    
    hPoverPi[iS][iH]->Divide(hPoverPi[iS][iH],hPi[iS][iH],1.);
    hLoverK0s[iS][iH]->Divide(hLoverK0s[iS][iH],hK0s[iS][iH],0.33);
    hXioverPhi[iS][iH] = DivideSpline(hXioverPhi[iS][iH],hPhi[iS][iH],.5);
    hXioverPhi[iS][iH]->GetXaxis()->SetRange(1,6);

  } }

  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1024,640);
  C->SetFillStyle(4000);

  // Number of PADS
  const int Nx = 3;
  const int Ny = 2;

  // Margins
  float lMargin = 0.05;
  float rMargin = 0.05;
  float bMargin = 0.05;
  float tMargin = 0.05;

  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  int latexTextSize = 18; 
  float legendTextSize = 0.065;

  const double minYspectra = 0.001;
  const double maxYspectra = 0.501;

  const double minYratio = 0.67;
  const double maxYratio = 1.62;

  const double minXleft = 0.24;
  const double maxXleft = 29.;

  const double minXmiddle = 0.24;
  const double maxXmiddle = 29.;

  const double minXright = 0.24;
  const double maxXright = 29.;


  // Left frames
  TH1D* hframeSpectraLeft = new TH1D("hframeSpectraLeft","",100 , minXleft , maxXleft);
  hframeSpectraLeft->GetYaxis()->SetRangeUser(minYspectra,maxYspectra);

  TH1D* hframeRatioLeft = new TH1D("hframeRatioLeft","",100 , minXleft , maxXleft);
  hframeRatioLeft->GetYaxis()->SetRangeUser(minYratio,maxYratio);
  hframeRatioLeft->GetXaxis()->SetMoreLogLabels(kTRUE);
  //hframeRatioLeft->GetYaxis()->SetMoreLogLabels(kTRUE);
  hframeRatioLeft->GetXaxis()->SetNoExponent(kTRUE);

  // Middle frames
  TH1D* hframeSpectraMiddle = new TH1D("hframeSpectraMiddle","",100 , minXmiddle , maxXmiddle);
  hframeSpectraMiddle->GetYaxis()->SetRangeUser(minYspectra,maxYspectra);

  TH1D* hframeRatioMiddle = new TH1D("hframeRatioMiddle","",100 , minXmiddle , maxXmiddle);
  hframeRatioMiddle->GetYaxis()->SetRangeUser(minYratio,maxYratio);
  hframeRatioMiddle->GetXaxis()->SetMoreLogLabels(kTRUE);
  hframeRatioMiddle->GetXaxis()->SetNoExponent(kTRUE);

  // Rights frames
  TH1D* hframeSpectraRight = new TH1D("hframeSpectraRight","",100 , minXright , maxXright);
  hframeSpectraRight->GetYaxis()->SetRangeUser(minYspectra,maxYspectra);

  TH1D* hframeRatioRight = new TH1D("hframeRatioRight","",100 , minXright , maxXright);
  hframeRatioRight->GetYaxis()->SetRangeUser(minYratio,maxYratio);
  //hframeRatioRight->GetXaxis()->SetMoreLogLabels(kTRUE);
  hframeRatioRight->GetXaxis()->SetMoreLogLabels(kTRUE);
  hframeRatioRight->GetXaxis()->SetNoExponent(kTRUE);

  TLine* lineleft = new TLine(minXleft, 1., maxXleft, 1.);
  lineleft->SetLineStyle(2);

  TLine* linemiddle = new TLine(minXmiddle, 1., maxXmiddle, 1.);
  linemiddle->SetLineStyle(2);

  TLine* lineright = new TLine(minXright, 1., maxXright, 1.);
  lineright->SetLineStyle(6);

  TPad* padTitleX = (TPad*) gROOT->FindObject("pad_TitleX");
  TPad* padTitleY1 = (TPad*) gROOT->FindObject("pad_TitleY1");
  TPad* padTitleY2 = (TPad*) gROOT->FindObject("pad_TitleY2");

  TLatex* latexTitleX = new TLatex();
  TLatex* latexTitleY1 = new TLatex();
  TLatex* latexTitleY2 = new TLatex();
  TuneLatex(latexTitleX,1.5*latexTextSize);
  TuneLatex(latexTitleY1,1.5*latexTextSize);
  TuneLatex(latexTitleY2,1.5*latexTextSize);

  TLatex* latexSp = new TLatex();
  TuneLatex(latexSp,2.*latexTextSize);

  TLatex* latexSystem = new TLatex();
  TuneLatex(latexSystem,latexTextSize);

  TLegend* legendLeft = new TLegend(0.20,0.7,0.7,0.89);
  myLegendSetUp(legendLeft,legendTextSize);
  legendLeft->AddEntry(hPoverPi[Iso][Stat],"Isotropic","P");
  legendLeft->AddEntry(hPoverPi[Jetty][Stat],"Jetty","P");

  TLegend* legendMiddle = new TLegend(0.05,0.05,0.55,0.35);
  myLegendSetUp(legendMiddle,legendTextSize);

  TLegend* legendRight = new TLegend(0.57,0.7,0.93,0.88);
  myLegendSetUp(legendRight,legendTextSize);
  legendRight->AddEntry(hPoverPi[Ref][Monash],"Monash","L");
  legendRight->AddEntry(hPoverPi[Ref][Ropes],"Ropes","L");

  TPad *pad[Nx][Ny];

   TF1 *f2=new TF1("f2","exp(x)",0,2);
   TGaxis *A2 = new TGaxis(1,1,0,0,"f2");
   A2->SetTitle("exponential axis");
   A2->SetLabelSize(0.03);
   A2->SetTitleSize(0.03);
   A2->SetTitleOffset(1.2);
//   A2->Draw();

  for (Int_t i=0;i<Nx;i++) {
    for (Int_t j=0;j<Ny;j++) {
      C->cd(0);

      // Get the pads previously created.
      char pname[16];
      sprintf(pname,"pad_%i_%i",i,j);
      pad[i][j] = (TPad*) gROOT->FindObject(pname);
      pad[i][j]->Draw();
      pad[i][j]->SetFillStyle(1001);
      //if(!j)pad[i][j]->SetLogy(kTRUE);
      pad[i][j]->SetLogx(kTRUE);
      
    }
  }

  C->cd(0);
  pad[0][1]->Draw();
  pad[0][1]->cd();
  DrawHistograms(hPoverPi[Ref],NHIST,colors[Ref]);
  hframeSpectraLeft->Draw();
  //DrawHistograms(hPoverPi[Ref],NHIST,colors[Ref]);
  DrawHistograms(hPoverPi[Iso],NHIST,colors[Iso]);
  DrawHistograms(hPoverPi[Jetty],NHIST,colors[Jetty]);

  latexSp->DrawLatex(0.3,0.75,"p/#pi");
  

  C->cd(0);
  pad[1][1]->Draw();
  pad[1][1]->cd();
  hframeSpectraMiddle->Draw();
  DrawHistograms(hLoverK0s[Iso],NHIST,colors[Iso]);
  DrawHistograms(hLoverK0s[Jetty],NHIST,colors[Jetty]);


  latexSp->DrawLatex(0.16,0.75,"#Lambda/K^{0}_{S}");
  latexSystem->DrawLatex(0.08,0.75,"#frac{1}{3}#times");
  
  C->cd(0);
  pad[2][1]->Draw();
  pad[2][1]->cd();
  hframeSpectraRight->Draw();
  DrawHistograms(hXioverPhi[Iso],NHIST,colors[Iso]);
  DrawHistograms(hXioverPhi[Jetty],NHIST,colors[Jetty]);

  latexSp->DrawLatex(0.16,0.75,"#Xi/#phi");
  latexSystem->DrawLatex(0.08,0.75,"#frac{1}{2}#times");

  C->cd(0);
  pad[0][0]->Draw();
  pad[0][0]->cd();
  hframeRatioLeft->Draw();
  lineleft->Draw("same");
  DrawHistogramsRatios(hPoverPi[Iso],hPoverPi[Ref],NHIST,colors[Iso]);
  DrawHistogramsRatios(hPoverPi[Jetty],hPoverPi[Ref],NHIST,colors[Iso]);

  legendLeft->Draw();
  legendRight->Draw();


  C->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  hframeRatioMiddle->Draw();
  linemiddle->Draw("same");
  DrawHistogramsRatios(hLoverK0s[Iso],hLoverK0s[Ref],NHIST,colors[Iso]);
  DrawHistogramsRatios(hLoverK0s[Jetty],hLoverK0s[Ref],NHIST,colors[Iso]);

  latexSystem->DrawLatex(0.1,0.9,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.1,0.8,"#it{N}_{SPD} I (|#eta|<0.8)");
  latexSystem->DrawLatex(0.1,0.7,"#it{N}_{ch} #geq 10, #it{p}_{T} #geq 0.15, |#eta|<0.8");
  

  C->cd(0);
  pad[2][0]->Draw();
  pad[2][0]->cd();
  hframeRatioRight->Draw();
  lineright->Draw("same");
  DrawHistogramsRatios(hXioverPhi[Iso],hXioverPhi[Ref],NHIST,colors[Iso]);
  DrawHistogramsRatios(hXioverPhi[Jetty],hXioverPhi[Ref],NHIST,colors[Iso]);

  latexSystem->DrawLatex(0.1,0.85,"ALICE");
  
  C->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX = "#it{p}_{T} (GeV/#it{c})";
  latexTitleX->DrawLatex(0.45, 0.5, TitleX);

  C->cd(0);
  padTitleY1->Draw();
  padTitleY1->cd();
  const char* TitleY1 = "Ratio to S^{#it{p}_{T}=1}_{O}-int.";
  latexTitleY1->SetTextAngle(90);
  latexTitleY1->DrawLatex(0.5, 0.25, TitleY1);

  C->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY2 = "Particle ratio";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.5, 0.25, TitleY2);

fout->Write();
C->SaveAs("./ratiosBaryonMeson.pdf");
C->SaveAs("./ratiosBaryonMeson.png");
C->SaveAs("./ratiosBaryonMeson.root");
	
}

