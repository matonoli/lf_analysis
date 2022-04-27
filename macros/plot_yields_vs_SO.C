#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TMath.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>

#include <iostream>
using namespace std;

/*
  Usage:
  .L plot_yields_vs_SO.C+

  draw_results(kFALSE, kTRUE)

 */

TFile* OpenFileFresh(const Char_t* fileName);
void SetHistErrorsToZero(TH1D* hist);
TH1D* GetErrorBand(TH1D* hist);
void ChangeBinHist(TH1D* oldBinning, TH1D* newBinning);
void ChangeBinGraph(TGraph* oldBinning, TH1D* newBinning);

void draw_results(Bool_t isV0M, Bool_t isVHM, Bool_t includeMC=kFALSE, Bool_t usePercentiles=kFALSE)
{
  gStyle->SetOptStat(0);
  //  gROOT->SetStyle("ATLAS");
  //  gStyle->SetErrorX(0.5);

  TFile* filePikp   = OpenFileFresh("./data_files/pi_k_p_spectra_ExclusiveBins.root");
  TFile* fileXI     = 0;
  TFile* fileL     = 0;
  TFile* file_Monash = 0;
  TFile* file_Ropes  = 0;

  const Char_t* figName = 0;

  // Pion references
  TH1D* hPi_stat = 0;
  TH1D* hPi_syst = 0;
  TH1D* hPi_ErrorBand = 0;

  // Double ratios
  TH1D* hPoverPi_stat = 0;
  TH1D* hPoverPi_syst = 0;

  TH1D* hXIoverPi_stat = 0;
  TH1D* hXIoverPi_syst = 0;

  TH1D* hLoverPi_stat = 0;
  TH1D* hLoverPi_syst = 0;

  
  {
    //
    // Get histograms
    //
    TH1D* hHelp    = 0; // Needed to make the piOverHm ratio
    
    if(isV0M) {
      
      if(isVHM) {

	figName = "_V0M_01";
	
	hPi_stat = (TH1D*)filePikp->Get("hYieldSO_Pion_V0M_binMult_0_sta");
	hPi_syst = (TH1D*)filePikp->Get("hYieldSO_Pion_V0M_binMult_0_sys");
	hHelp = (TH1D*)filePikp->Get("hYieldHM_Pion_V0M_binMult_0_sta");
	
	hPoverPi_stat = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_V0M_binMult_0_sta");
	hPoverPi_syst = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_V0M_binMult_0_sys");

        fileL = OpenFileFresh("./data_files/l_vs_s0_V0M01.root");

        fileXI   = OpenFileFresh("./data_files/xi_results_Vs_So_V0M_top1_nov_2.root");
        file_Monash = OpenFileFresh("./make_mc_graphs/Monash_V0M_VHM_oct_5.root");
        file_Ropes  = OpenFileFresh("./make_mc_graphs/Ropes_V0M_VHM_oct_5.root");
      } else {
	
	figName = "_V0M_10";

	hPi_stat = (TH1D*)filePikp->Get("hYieldSO_Pion_V0M_binMult_1_sta");
	hPi_syst = (TH1D*)filePikp->Get("hYieldSO_Pion_V0M_binMult_1_sys");
	hHelp = (TH1D*)filePikp->Get("hYieldHM_Pion_V0M_binMult_1_sta");
	
	hPoverPi_stat = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_V0M_binMult_1_sta");
	hPoverPi_syst = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_V0M_binMult_1_sys");

  fileL = OpenFileFresh("./data_files/l_vs_s0_V0M.root");

	fileXI   = OpenFileFresh("./data_files/xi_results_Vs_So_V0M_nov_2.root");
        file_Monash = OpenFileFresh("./make_mc_graphs/Monash_V0M_oct_5.root");
        file_Ropes  = OpenFileFresh("./make_mc_graphs/Ropes_V0M_oct_5.root");
      }
    } else {
      if(isVHM) {
	
	figName = "_CL1_01";

	hPi_stat = (TH1D*)filePikp->Get("hYieldSO_Pion_CL1_binMult_0_sta");
	hPi_syst = (TH1D*)filePikp->Get("hYieldSO_Pion_CL1_binMult_0_sys");
	hHelp = (TH1D*)filePikp->Get("hYieldHM_Pion_CL1_binMult_0_sta");
	
	hPoverPi_stat = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_Trks_binMult_0_sta");
	hPoverPi_syst = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_Trks_binMult_0_sys");

  fileL = OpenFileFresh("./data_files/l_vs_s0_NCharged01.root");

	fileXI   = OpenFileFresh("./data_files/xi_results_Vs_So_CL1_top1_nov_2.root");
        file_Monash = OpenFileFresh("./make_mc_graphs/Monash_CL1_VHM_oct_5.root");
        file_Ropes  = OpenFileFresh("./make_mc_graphs/Ropes_CL1_VHM_oct_5.root");
      } else {
	
	figName = "_CL1_10";
	
	hPi_stat = (TH1D*)filePikp->Get("hYieldSO_Pion_CL1_binMult_1_sta");
	hPi_syst = (TH1D*)filePikp->Get("hYieldSO_Pion_CL1_binMult_1_sys");
	hHelp = (TH1D*)filePikp->Get("hYieldHM_Pion_CL1_binMult_1_sta");
	
	hPoverPi_stat = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_Trks_binMult_1_sta");
	hPoverPi_syst = (TH1D*)filePikp->Get("hProton2Pion_DoubleRatio_Trks_binMult_1_sys");

  fileL = OpenFileFresh("./data_files/l_vs_s0_NCharged.root");

      	fileXI   = OpenFileFresh("./data_files/xi_results_Vs_So_CL1_nov_2.root");
        file_Monash = OpenFileFresh("./make_mc_graphs/Monash_CL1_oct_5.root");
        file_Ropes  = OpenFileFresh("./make_mc_graphs/Ropes_CL1_oct_5.root");
      }
    }

    // //
    // // To debug
    // //
    // {
    //   TH1D* hist_yield = (TH1D*)fileXI->Get("hist_yield");

    //   for(Int_t n=0; n < 9; n++) {

    // 	cout << "n: " << n << ", Npi: " << hPi_stat->GetBinContent(n+1)
    // 	     << ", Nxi: " << hist_yield->GetBinContent(n+1) << endl;
    //   }

    // }
    
    hPi_ErrorBand = GetErrorBand(hPi_syst);

    hPi_stat->Divide(hHelp);
    hPi_syst->Divide(hHelp);

    SetHistErrorsToZero(hPi_syst);

    hLoverPi_stat = (TH1D*)fileL->Get("hLLbarvS0_Stat");
    hLoverPi_syst = (TH1D*)fileL->Get("hLLbarvS0_Syst");
    
    hXIoverPi_stat = (TH1D*)fileXI->Get("hist_double_ratio_stat");
    hXIoverPi_syst = (TH1D*)fileXI->Get("hist_double_ratio_syst");
    
    if(usePercentiles) {

      ChangeBinHist(hXIoverPi_stat, hPi_stat);
      ChangeBinHist(hXIoverPi_syst, hPi_stat);
    } else {

      ChangeBinHist(hPi_ErrorBand, hXIoverPi_stat);
      ChangeBinHist(hPi_stat, hXIoverPi_stat);
      ChangeBinHist(hPi_syst, hXIoverPi_stat);

      ChangeBinHist(hLoverPi_stat, hXIoverPi_stat);
      ChangeBinHist(hLoverPi_syst, hXIoverPi_stat);

      ChangeBinHist(hPoverPi_stat, hXIoverPi_stat);
      ChangeBinHist(hPoverPi_syst, hXIoverPi_stat);
    }
    
    // This is how we do it:
    // (Xi_SO/pi_SO) / (XI_HM/pi_HM) =  (Xi_SO/XI_HM) / (pi_SO/pi_HM);
    hLoverPi_stat->Divide(hPi_stat);     
    hLoverPi_syst->Divide(hPi_syst);

    hXIoverPi_stat->Divide(hPi_stat);     
    hXIoverPi_syst->Divide(hPi_syst); 
  }

  
  hPi_ErrorBand->SetTitle("; Uncorrected #it{S}_{O}^{#it{p}_{T}=1}; Ratio to pions / (HM ratio)");
  if(usePercentiles)
    hPi_ErrorBand->GetXaxis()->SetTitle("#it{S}_{O}^{#it{p}_{T}=1} percentile");

  hPi_ErrorBand->GetYaxis()->SetRangeUser(0.75, 1.15);
  hPi_ErrorBand->GetYaxis()->CenterTitle();
  hPi_ErrorBand->GetYaxis()->SetTitleSize(0.06);
  hPi_ErrorBand->GetYaxis()->SetTitleOffset(1.1);
  hPi_ErrorBand->GetYaxis()->SetLabelSize(0.05);
  hPi_ErrorBand->GetXaxis()->SetTitleSize(0.06);
  hPi_ErrorBand->GetXaxis()->SetTitleOffset(0.95);
  hPi_ErrorBand->GetXaxis()->SetLabelSize(0.05);
  hPi_ErrorBand->GetYaxis()->SetNdivisions(505);
  hPi_ErrorBand->GetXaxis()->SetNdivisions(505);

  TCanvas * cRatio = new TCanvas("cRatio", "Ratio", 600, 600);
  cRatio->SetMargin(0.14, 0.02, 0.14, 0.02);
  //  cRatio->SetGridy();
  cRatio->cd();

  hPi_ErrorBand->SetFillColor(kGray+1);
  hPi_ErrorBand->SetFillStyle(1001);
  hPi_ErrorBand->Draw("E2");

  hPoverPi_syst->SetFillStyle(0);
  hPoverPi_syst->SetMarkerStyle(1);
  hPoverPi_syst->SetLineColor(1);
  hPoverPi_syst->SetLineWidth(1);
  hPoverPi_syst->Draw("E2 SAME");

  hLoverPi_syst->SetFillStyle(0);
  hLoverPi_syst->SetMarkerStyle(1);
  hLoverPi_syst->SetLineColor(kBlue+2);
  hLoverPi_syst->SetLineWidth(1);
  hLoverPi_syst->Draw("E2 SAME");

  hXIoverPi_syst->SetFillStyle(0);
  hXIoverPi_syst->SetMarkerStyle(1);
  hXIoverPi_syst->SetLineColor(kGreen+2);
  hXIoverPi_syst->SetLineWidth(1);
  hXIoverPi_syst->Draw("E2 SAME");

  hPoverPi_stat->SetFillStyle(0);
  hPoverPi_stat->SetLineColor(1);
  hPoverPi_stat->SetMarkerStyle(20);
  hPoverPi_stat->SetMarkerColor(1);
  hPoverPi_stat->SetMarkerSize(1.2);
  hPoverPi_stat->Draw("SAME");

  hLoverPi_stat->SetFillStyle(0);
  hLoverPi_stat->SetLineColor(1);
  hLoverPi_stat->SetLineWidth(1);
  hLoverPi_stat->SetMarkerStyle(20);
  hLoverPi_stat->SetMarkerColor(kBlue+2);
  hLoverPi_stat->SetMarkerSize(1.2);
  hLoverPi_stat->Draw("SAME");

  hXIoverPi_stat->SetFillStyle(0);
  hXIoverPi_stat->SetLineColor(1);
  hXIoverPi_stat->SetLineWidth(1);
  hXIoverPi_stat->SetMarkerStyle(20);
  hXIoverPi_stat->SetMarkerColor(kGreen+2);
  hXIoverPi_stat->SetMarkerSize(1.2);
  hXIoverPi_stat->Draw("SAME");

  if(includeMC) {

    TGraph* p_over_pi_Monash = (TGraph*)file_Monash->Get("g_p_over_pi");
    p_over_pi_Monash->SetLineStyle(2);

    TGraph* p_over_pi_Ropes = (TGraph*)file_Ropes->Get("g_p_over_pi");

    TGraph* la_over_pi_Monash = (TGraph*)file_Monash->Get("g_la_over_pi");
    la_over_pi_Monash->SetLineStyle(2);

    TGraph* la_over_pi_Ropes = (TGraph*)file_Ropes->Get("g_la_over_pi");

    TGraph* xi_over_pi_Monash = (TGraph*)file_Monash->Get("g_xi_over_pi");
    xi_over_pi_Monash->SetLineStyle(2);

    TGraph* xi_over_pi_Ropes = (TGraph*)file_Ropes->Get("g_xi_over_pi");

    if(!usePercentiles) {

      ChangeBinGraph(p_over_pi_Monash, hXIoverPi_stat);
      ChangeBinGraph(p_over_pi_Ropes, hXIoverPi_stat);
      
      ChangeBinGraph(la_over_pi_Monash, hXIoverPi_stat);
      ChangeBinGraph(la_over_pi_Ropes, hXIoverPi_stat);
      
      ChangeBinGraph(xi_over_pi_Monash, hXIoverPi_stat);
      ChangeBinGraph(xi_over_pi_Ropes, hXIoverPi_stat);
    }

    p_over_pi_Monash->Draw("L SAME");
    p_over_pi_Ropes->Draw("L SAME");

    la_over_pi_Monash->Draw("L SAME");
    la_over_pi_Ropes->Draw("L SAME");

    xi_over_pi_Monash->Draw("L SAME");
    xi_over_pi_Ropes->Draw("L SAME");
  }

  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextFont(42);

  latex.DrawLatex(0.17, 0.92, "#sqrt{#it{s}}=13 TeV, ");

  
  TString helpString(figName);

  if(helpString.Contains("CL1")) {
      
    if(helpString.Contains("10"))
      latex.DrawLatex(0.41, 0.92, "CL1 0-10%, |#eta|<0.8, N_{ch}#geq10");
    else
      latex.DrawLatex(0.41, 0.92, "CL1 0-1%, |#eta|<0.8, N_{ch}#geq10");
  } else {
    
    if(helpString.Contains("10"))
      latex.DrawLatex(0.41, 0.92, "V0M 0-10%, |#eta|<0.8, N_{ch}#geq10");
    else
      latex.DrawLatex(0.41, 0.92, "V0M 0-1%, |#eta|<0.8, N_{ch}#geq10");
  }
  
  latex.SetTextSize(0.03);
  latex.DrawLatex(0.65, 0.30, "N_{#pi}: 0.3 < #it{p}_{T} < 20 GeV/#it{c}");
  latex.DrawLatex(0.65, 0.26, "N_{p}: 0.45 < #it{p}_{T} < 20 GeV/#it{c}");
  latex.DrawLatex(0.65, 0.22, "N_{#Lambda}: 1.0 < #it{p}_{T} < 8 GeV/#it{c}");
  latex.DrawLatex(0.65, 0.18, "N_{#Xi}: 0.6 < #it{p}_{T} < 6.5 GeV/#it{c}");

  
  auto l = new TLegend(0.75,0.34,0.94,0.55);
  l->SetBorderSize(0.);
  l->SetTextFont(42);
  l->AddEntry(hPoverPi_stat, "N_{p} / N_{#pi}", "PL");
  l->AddEntry(hLoverPi_stat, "N_{#Lambda} / N_{#pi}", "PL");
  l->AddEntry(hXIoverPi_stat, "N_{#Xi} / N_{#pi}", "PL");
  l->Draw();

  
  cRatio->SaveAs(Form("xi_over_pi_vs_So%s.pdf", figName));
  cRatio->SaveAs(Form("xi_over_pi_vs_So%s.gif", figName));
}

//______________________________________________________________________
TFile* OpenFileFresh(const Char_t* fileName)
{
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    file->Close();
    delete file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

//______________________________________________________________________
void SetHistErrorsToZero(TH1D* hist)
{
  const Int_t nBinsX = hist->GetNbinsX();
  for(Int_t bin = 1; bin <= nBinsX; bin++) {

    hist->SetBinError(bin, 0);
  }
}

//______________________________________________________________________
TH1D* GetErrorBand(TH1D* hist)
{
  TH1D* hErrorBand = (TH1D*)hist->Clone("hErrorBand");
  const Int_t nBinsX = hist->GetNbinsX();
  for(Int_t bin = 1; bin <= nBinsX; bin++) {

    hErrorBand->SetBinContent(bin, 1);
    hErrorBand->SetBinError(bin, hist->GetBinError(bin)/hist->GetBinContent(bin));
  }

  return hErrorBand;
}

//______________________________________________________________________
void ChangeBinHist(TH1D* oldBinning, TH1D* newBinning)
{
  const Int_t nBinsX = newBinning->GetNbinsX();
  Double_t xNew[nBinsX+1];
  
  for(Int_t bin = 1; bin <= nBinsX; bin++) {
    
    xNew[bin-1] = newBinning->GetXaxis()->GetBinLowEdge(bin);
  }

  xNew[nBinsX] = newBinning->GetXaxis()->GetBinUpEdge(nBinsX);
  for(Int_t n = 0; n <= nBinsX; n++) {

    cout << n << ": " << xNew[n] << endl;
  }
  
  oldBinning->GetXaxis()->Set(nBinsX, xNew);
}

//______________________________________________________________________
void ChangeBinGraph(TGraph* oldBinning, TH1D* newBinning)
{
  const Int_t nBinsX = newBinning->GetNbinsX();
  Double_t xNew[nBinsX];
  
  for(Int_t bin = 1; bin <= nBinsX; bin++) {

    double x,y;
    oldBinning->GetPoint(bin-1,x,y);
    oldBinning->SetPoint(bin-1, newBinning->GetXaxis()->GetBinCenter(bin),y);
  }
}
