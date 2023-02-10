#include <iostream>
#include "MakeCanvas.C"

using namespace std;


//// UNCOMMENT A LINE TO CHOOSE THE MULTIPLICITY SELECTION
TString mult("V0M01"); int multInt = 0; int multOmar = 0; TString multPlot("V0M I");
//TString mult("NCharged01"); int multInt = 1; int multOmar = 0; TString multPlot("#it{N}_{SPD} I (|#eta|<0.8)");
//TString mult("V0M"); int multInt = 2; int multOmar = 1; TString multPlot("V0M III");
//TString mult("NCharged"); int multInt = 3; int multOmar = 1; TString multPlot("#it{N}_{SPD} III (|#eta|<0.8)");

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



void plotFinalKtoK() {

  TH1::SetDefaultSumw2(1);

  TH1D* hK0s[NSPHERO][NHIST];
  TH1D* hK[NSPHERO][NHIST];
  TH1D* hK_r[NSPHERO][NHIST];
  TH1D* hK0soverK[NSPHERO][NHIST];

  // kpm
  TFile* fPion = new TFile(Form("%s/pi_k_p_spectra.root",path.Data()),"READ");
  const char* strM[4] = {"V0M", "Trks", "V0M", "Trks"};
  hK[Ref][Stat]      = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hKaon_stat",strM[multInt],multOmar));
  hK[Jetty][Stat]    = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hKaon_stat",strM[multInt],sphInt,multOmar));
  hK[Iso][Stat]      = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hKaon_stat",strM[multInt],sphInt,multOmar));
  hK[Ref][Syst]      = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hKaon_tot_sys",strM[multInt],multOmar));
  hK[Jetty][Syst]    = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hKaon_tot_sys",strM[multInt],sphInt,multOmar));
  hK[Iso][Syst]      = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hKaon_tot_sys",strM[multInt],sphInt,multOmar));
  hK[Ref][SystUnc]   = (TH1D*)fPion->Get(Form("%s_Reference_Bin_mult_%i/hKaon_unc_sys_am",strM[multInt],multOmar));
  hK[Jetty][SystUnc] = (TH1D*)fPion->Get(Form("%s_Jetty_%i_Bin_mult_%i/hKaon_unc_sys_am",strM[multInt],sphInt,multOmar));
  hK[Iso][SystUnc]   = (TH1D*)fPion->Get(Form("%s_Isotropic_%i_Bin_mult_%i/hKaon_unc_sys_am",strM[multInt],sphInt,multOmar));

  // k0s
  TFile* fK0s = new TFile(Form("%s/k0s_l_spectra_%s_spher%s.root",path.Data(),mult.Data(),spher.Data()),"READ");
  hK0s[Ref][Stat]       = (TH1D*)fK0s->Get("hK0s_Ref_Stat");
  hK0s[Jetty][Stat]     = (TH1D*)fK0s->Get("hK0s_Jet_Stat");
  hK0s[Iso][Stat]       = (TH1D*)fK0s->Get("hK0s_Iso_Stat");
  hK0s[Ref][Syst]       = (TH1D*)fK0s->Get("hK0s_Ref_Syst");
  hK0s[Jetty][Syst]     = (TH1D*)fK0s->Get("hK0s_Jet_Syst");
  hK0s[Iso][Syst]       = (TH1D*)fK0s->Get("hK0s_Iso_Syst");
  hK0s[Ref][SystUnc]    = (TH1D*)fK0s->Get("hK0s_Ref_SystUnc");   
  hK0s[Jetty][SystUnc]  = (TH1D*)fK0s->Get("hK0s_Jet_SystUnc");
  hK0s[Iso][SystUnc]    = (TH1D*)fK0s->Get("hK0s_Iso_SystUnc");

  //monash
  TFile* fMonash = new TFile(Form("%s/Complete_Monash.root",path.Data()),"READ");
  const char* strM_MC[4] = { "V0M","CL1","V0M","CL1"};
  const char* strS_MC[4] = { "01","05","10","20"};
  hK[Ref][Monash]    = (TH1D*)fMonash->Get(Form("fKaons_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hK[Jetty][Monash]  = (TH1D*)fMonash->Get(Form("fKaons_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK[Iso][Monash]    = (TH1D*)fMonash->Get(Form("fKaons_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Ref][Monash]   = (TH1D*)fMonash->Get(Form("fK0s_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hK0s[Jetty][Monash] = (TH1D*)fMonash->Get(Form("fK0s_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Iso][Monash]   = (TH1D*)fMonash->Get(Form("fK0s_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  
  //ropes
  TFile* fRopes = new TFile(Form("%s/Complete_Ropes.root",path.Data()),"READ");
  hK[Ref][Ropes]    = (TH1D*)fRopes->Get(Form("fKaons_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hK[Jetty][Ropes]  = (TH1D*)fRopes->Get(Form("fKaons_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK[Iso][Ropes]    = (TH1D*)fRopes->Get(Form("fKaons_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Ref][Ropes]   = (TH1D*)fRopes->Get(Form("fK0s_HM%s_%s",strM_MC[multInt],multInt>1?"10":"01"));
  hK0s[Jetty][Ropes] = (TH1D*)fRopes->Get(Form("fK0s_%s_%s_Jetty_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  hK0s[Iso][Ropes]   = (TH1D*)fRopes->Get(Form("fK0s_%s_%s_Iso_%s",strM_MC[multInt],multInt>1?"10":"01",strS_MC[sphInt]));
  //fRopes->Close();


  TFile* fout = new TFile(Form("./KtoK_%s_spher%s.root",mult.Data(),spher.Data()),"RECREATE");


  // REBIN KPM
  const Int_t NPTBINS = 13;
  const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
    1.0, 1.2, 
    1.4, 1.6, 1.8, 2.0, 2.2, 
    2.6, 3.0, 3.4, 4.0, 5.0, 
    6.5, 8.0 };

  const Int_t NPIONBINS = 51;   //Omar pi+- spectra
  const Double_t PIONBINS[NPIONBINS+1] = {
    0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 
    0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 
    0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 
    1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 
    2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.50, 5.00, 
    5.50, 6.00, 6.50, 7.00, 8.00, 10.0, 20.0 };

  for (int iS = 0; iS < NSPHERO; ++iS)  {
    
    hK_r[iS][Stat] = (TH1D*)hK[iS][Stat]->Rebin(NPTBINS,hK[iS][Stat]->GetName(),XBINS);
    hK_r[iS][Syst] = (TH1D*)hK[iS][Syst]->Rebin(NPTBINS,hK[iS][Syst]->GetName(),XBINS);
    hK_r[iS][SystUnc] = (TH1D*)hK[iS][SystUnc]->Rebin(NPTBINS,hK[iS][SystUnc]->GetName(),XBINS);

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
              errStat += hK[iS][Stat]->GetBinError(1+iPB2)*hK[iS][Stat]->GetBinError(1+iPB2);
              errSyst += hK[iS][Syst]->GetBinError(1+iPB2);
              errSystUnc += hK[iS][SystUnc]->GetBinError(1+iPB2);
            }
            else break;
          }
        }
      }
      
      hK_r[iS][Stat]->SetBinContent(1+iB,1./divF*hK_r[iS][Stat]->GetBinContent(1+iB));
      hK_r[iS][Stat]->SetBinError(1+iB,1./divF*TMath::Sqrt(errStat));
      hK_r[iS][Syst]->SetBinContent(1+iB,1./divF*hK_r[iS][Syst]->GetBinContent(1+iB));
      hK_r[iS][Syst]->SetBinError(1+iB,1./divF*errSyst);
      hK_r[iS][SystUnc]->SetBinContent(1+iB,1./divF*hK_r[iS][SystUnc]->GetBinContent(1+iB));
      hK_r[iS][SystUnc]->SetBinError(1+iB,1./divF*errSystUnc);
      
    }

    hK[iS][Stat] = hK_r[iS][Stat];
    hK[iS][Syst] = hK_r[iS][Syst];
    hK[iS][SystUnc] = hK_r[iS][SystUnc];
  }

  // CALCULATE RATIOS
  for (int iS = 0; iS < NSPHERO; ++iS)  {
  for (int iH = 0; iH < NHIST; ++iH)  {

    hK0soverK[iS][iH] = (TH1D*)hK0s[iS][iH]->Clone(Form("hK0soverK_%s_%s",strS[iS],strH[iH]));
    hK0soverK[iS][iH]->Divide(hK0soverK[iS][iH],hK[iS][iH],2.);
    

  } }


  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",700,600);
  C->SetFillStyle(4000);

  // Number of PADS
  const int Nx = 1;
  const int Ny = 1;

  // Margins
  float lMargin = 0.05;
  float rMargin = 0.05;
  float bMargin = 0.05;
  float tMargin = 0.05;

  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  int latexTextSize = 18; 
  float legendTextSize = 0.035;

  const double minYratio = 0.65;
  const double maxYratio = 1.55;

  const double minXleft = 0.49;
  const double maxXleft = 15.;


  // Left frames
  TH1D* hframeRatioLeft = new TH1D("hframeRatioLeft","",100 , minXleft , maxXleft);
  hframeRatioLeft->GetYaxis()->SetRangeUser(minYratio,maxYratio);
  hframeRatioLeft->GetXaxis()->SetMoreLogLabels(kTRUE);
  //hframeRatioLeft->GetYaxis()->SetMoreLogLabels(kTRUE);
  hframeRatioLeft->GetXaxis()->SetNoExponent(kTRUE);


  TLine* lineleft = new TLine(minXleft, 1., maxXleft, 1.);
  lineleft->SetLineStyle(2);


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

  TLegend* legendLeft = new TLegend(0.51,0.71,0.91,0.90);
  myLegendSetUp(legendLeft,legendTextSize);
  legendLeft->AddEntry(hK0soverK[Ref][Stat],"S_{0}^{p_{T}=1}-integrated","P");
  legendLeft->AddEntry(hK0soverK[Iso][Stat],"Isotropic","P");
  legendLeft->AddEntry(hK0soverK[Jetty][Stat],"Jetty","P");

  TLegend* legendRight = new TLegend(0.8,0.7,0.95,0.85);
  myLegendSetUp(legendRight,legendTextSize);
  legendRight->AddEntry(hK0soverK[Ref][Monash],"Monash","L");
  legendRight->AddEntry(hK0soverK[Ref][Ropes],"Ropes","L");

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
  pad[0][0]->Draw();
  pad[0][0]->cd();
  hframeRatioLeft->Draw();
  DrawHistograms(hK0soverK[Ref],NHIST,colors[Ref]);
  //DrawHistograms(hPoverPi[Ref],NHIST,colors[Ref]);
  DrawHistograms(hK0soverK[Iso],NHIST,colors[Iso]);
  DrawHistograms(hK0soverK[Jetty],NHIST,colors[Jetty]);

  latexSystem->DrawLatex(0.1,0.9,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.1,0.83,multPlot.Data());
  latexSystem->DrawLatex(0.1,0.76,"#it{N}_{ch} #geq 10, #it{p}_{T} #geq 0.15, |#eta|<0.8");
  latexSystem->DrawLatex(0.1,0.69,"ALICE");

  legendLeft->Draw();
  legendRight->Draw();

  C->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX = "#it{p}_{T} (GeV/#it{c})";
  latexTitleX->DrawLatex(0.45, 0.5, TitleX);

  C->cd(0);
  padTitleY1->Draw();
  padTitleY1->cd();
  const char* TitleY1 = "2x K_{S}^{0} /";
  latexTitleY1->SetTextAngle(90);
  latexTitleY1->DrawLatex(0.5, 0.70, TitleY1);
  C->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY2 = " K^{#pm}";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.5, 0.0, TitleY2);

  fout->Write();
  C->SaveAs(Form("./KtoK_%s_spher%s.pdf",mult.Data(),spher.Data()));
  C->SaveAs(Form("./KtoK_%s_spher%s.png",mult.Data(),spher.Data()));


}

