#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TLine.h>
#include <TPaveLabel.h>
#include <iostream>
#include <fstream>
#include <TPavesText.h>
#include <TPad.h>
#include <TLatex.h>

using namespace std;
/////////////////////////////////////////////////////////////////////////////
// USAGE:
// Compile: .L DrawBaryons.C+
// run DrawBaryon()
//Change any of the ShapeMy"X" and my"X" functions below to adjust general style changes.
//Default will plot all 3 variatns (HM, Jetty,Iso). Comment out the ones you dont want to plot.
//Will update into an array in the near future

void ShapeMyHist(TH1D **hist);
void ShapeMyHistMonash(TH1D **hist);
void ShapeMyHistRopes(TH1D **hist);
void myOptions();
void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.04);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);

void DrawMeson(){
//LOAD EVERYTHING
/////////////////////////////////////////////////////////////////////////////
  myOptions();
  gStyle->SetOptStat(0);
  TFile* Kstars = TFile::Open("Kstar_CL1_0_1_sp10.root");
  TH1D* hKstarHM = (TH1D*)Kstars->Get("h_kstarpi_minbias_stat");
  TH1D* hKstarJ = (TH1D*)Kstars->Get("h_kstarpi_jetty_stat");
  TH1D* hKstarI = (TH1D*)Kstars->Get("h_kstarpi_iso_stat");
  TFile* Kaons = TFile::Open("pi_k_p_spectra.root");
  TH1D* hKaonHM = (TH1D*)Kaons->Get("Trks_Reference_Bin_mult_1/hKaon_stat");
  TH1D* hKaonJ = (TH1D*)Kaons->Get("Trks_Jetty_2_Bin_mult_1/hKaon_stat");
  TH1D* hKaonI = (TH1D*)Kaons->Get("Trks_Isotropic_2_Bin_mult_1/hKaon_stat");
  TFile* Pis = TFile::Open("pi_k_p_spectra.root");
  TH1D* hPiHM = (TH1D*)Pis->Get("Trks_Reference_Bin_mult_1/hPion_stat");
  TH1D* hPiJ = (TH1D*)Pis->Get("Trks_Jetty_2_Bin_mult_1/hPion_stat");
  TH1D* hPiI = (TH1D*)Pis->Get("Trks_Isotropic_2_Bin_mult_1/hPion_stat");
  TFile* K0s = TFile::Open("k0s_l_spectra_NCharged01_spher10.root");
  TH1D* hK0HM = (TH1D*)K0s->Get("hK0sToPi_Ref_Stat");
  TH1D* hK0J = (TH1D*)K0s->Get("hK0sToPi_Jet_Stat");
  TH1D* hK0I = (TH1D*)K0s->Get("hK0sToPi_Iso_Stat");
  
  TFile* Monash = TFile::Open("Complete_Monash.root");
  TFile* Ropes = TFile::Open("Complete_Ropes.root");

  TH1D* hKstarHM_Monash = (TH1D*)Monash->Get("fKstars_HMCL1_01");
  TH1D* hKstarJ_Monash = (TH1D*)Monash->Get("fKstars_CL1_01_Jetty_10");
  TH1D* hKstarI_Monash = (TH1D*)Monash->Get("fKstars_CL1_01_Iso_10");
  TH1D* hPiHM_Monash = (TH1D*)Monash->Get("fPions_HMCL1_01");
  TH1D* hPiJ_Monash = (TH1D*)Monash->Get("fPions_CL1_01_Jetty_10");
  TH1D* hPiI_Monash = (TH1D*)Monash->Get("fPions_CL1_01_Iso_10");
  TH1D* hKstarHM_Ropes = (TH1D*)Ropes->Get("fKstars_HMCL1_01");
  TH1D* hKstarJ_Ropes = (TH1D*)Ropes->Get("fKstars_CL1_01_Jetty_10");
  TH1D* hKstarI_Ropes = (TH1D*)Ropes->Get("fKstars_CL1_01_Iso_10");
  TH1D* hPiHM_Ropes = (TH1D*)Ropes->Get("fPions_HMCL1_01");
  TH1D* hPiJ_Ropes = (TH1D*)Ropes->Get("fPions_CL1_01_Jetty_10");
  TH1D* hPiI_Ropes = (TH1D*)Ropes->Get("fPions_CL1_01_Iso_10");
  TH1D* hKaonHM_Ropes = (TH1D*)Ropes->Get("fKaons_HMCL1_01");
  TH1D* hKaonJ_Ropes = (TH1D*)Ropes->Get("fKaons_CL1_01_Jetty_10");
  TH1D* hKaonI_Ropes = (TH1D*)Ropes->Get("fKaons_CL1_01_Iso_10");
  TH1D* hKaonHM_Monash = (TH1D*)Monash->Get("fKaons_HMCL1_01");
  TH1D* hKaonJ_Monash = (TH1D*)Monash->Get("fKaons_CL1_01_Jetty_10");
  TH1D* hKaonI_Monash = (TH1D*)Monash->Get("fKaons_CL1_01_Iso_10");
  TH1D* hK0HM_Ropes = (TH1D*)Ropes->Get("fK0s_HMCL1_01");
  TH1D* hK0J_Ropes = (TH1D*)Ropes->Get("fK0s_CL1_01_Jetty_10");
  TH1D* hK0I_Ropes = (TH1D*)Ropes->Get("fK0s_CL1_01_Iso_10");
  TH1D* hK0HM_Monash = (TH1D*)Monash->Get("fK0s_HMCL1_01");
  TH1D* hK0J_Monash = (TH1D*)Monash->Get("fK0s_CL1_01_Jetty_10");
  TH1D* hK0I_Monash = (TH1D*)Monash->Get("fK0s_CL1_01_Iso_10");
//Done Loading
/////////////////////////////////////////////////////////////////////////////
//Set Line Color, and Enable the Style modifcations for each histo
  hKaonHM->SetLineColor(1);
  hKaonJ->SetLineColor(2);  
  hKaonHM->Divide(hPiHM);
  hKaonJ->Divide(hPiJ);
  hKaonI->Divide(hPiI);
  
  ShapeMyHistMonash(&hKstarHM_Monash);
  ShapeMyHistMonash(&hKstarJ_Monash);
  ShapeMyHistMonash(&hKstarI_Monash);
  hKstarHM_Monash->Divide(hPiHM_Monash);
  hKstarJ_Monash->Divide(hPiJ_Monash);
  hKstarI_Monash->Divide(hPiI_Monash);
  hKstarHM_Monash->SetLineColor(1);
  hKstarJ_Monash->SetLineColor(2);
  ShapeMyHistRopes(&hKstarHM_Ropes);
  ShapeMyHistRopes(&hKstarJ_Ropes);
  ShapeMyHistRopes(&hKstarI_Ropes);
  hKstarHM_Ropes->Divide(hPiHM_Ropes);
  hKstarJ_Ropes->Divide(hPiJ_Ropes);
  hKstarI_Ropes->Divide(hPiI_Ropes);
  hKstarHM_Ropes->SetLineColor(1);
  hKstarJ_Ropes->SetLineColor(2);


  ShapeMyHistMonash(&hKaonHM_Monash);
  ShapeMyHistMonash(&hKaonJ_Monash);
  ShapeMyHistMonash(&hKaonI_Monash);
  hKaonHM_Monash->Divide(hPiHM_Monash);
  hKaonJ_Monash->Divide(hPiJ_Monash);
  hKaonI_Monash->Divide(hPiI_Monash);
  hKaonHM_Monash->SetLineColor(1);
  hKaonJ_Monash->SetLineColor(2);
  ShapeMyHistRopes(&hKaonHM_Ropes);
  ShapeMyHistRopes(&hKaonJ_Ropes);
  ShapeMyHistRopes(&hKaonI_Ropes);
  hKaonHM_Ropes->Divide(hPiHM_Ropes);
  hKaonJ_Ropes->Divide(hPiJ_Ropes);
  hKaonI_Ropes->Divide(hPiI_Ropes);
  hKaonHM_Ropes->SetLineColor(1);
  hKaonJ_Ropes->SetLineColor(2);
 
  ShapeMyHistMonash(&hK0HM_Monash);
  ShapeMyHistMonash(&hK0J_Monash);
  ShapeMyHistMonash(&hK0I_Monash);
  hK0HM_Monash->Divide(hPiHM_Monash);
  hK0J_Monash->Divide(hPiJ_Monash);
  hK0I_Monash->Divide(hPiI_Monash);
  hK0HM_Monash->SetLineColor(1);
  hK0J_Monash->SetLineColor(2);
  ShapeMyHistRopes(&hK0HM_Ropes);
  ShapeMyHistRopes(&hK0J_Ropes);
  ShapeMyHistRopes(&hK0I_Ropes);
  hK0HM_Ropes->Divide(hPiHM_Ropes);
  hK0J_Ropes->Divide(hPiJ_Ropes);
  hK0I_Ropes->Divide(hPiI_Ropes);
  hK0HM_Ropes->SetLineColor(1);
  hK0J_Ropes->SetLineColor(2);

  ShapeMyHist(&hKstarHM);
  ShapeMyHist(&hKstarJ);
  ShapeMyHist(&hKstarI);  
  ShapeMyHist(&hKaonHM);
  ShapeMyHist(&hKaonJ);
  ShapeMyHist(&hKaonI);
  ShapeMyHist(&hK0HM);
  ShapeMyHist(&hK0J);
  ShapeMyHist(&hK0I);
  
//Histos prep'd, and ready to go
/////////////////////////////////////////////////////////////////////////////
//Rebin Monte-Carlo
  hKstarJ_Ropes->Rebin(8);
  hKstarJ_Ropes->Scale(0.125);
  hKstarI_Ropes->Rebin(8);
  hKstarI_Ropes->Scale(0.125);
  hKstarJ_Monash->Rebin(8);
  hKstarJ_Monash->Scale(0.125);
  hKstarI_Monash->Rebin(8);
  hKstarI_Monash->Scale(0.125);
  //
  hK0J_Ropes->Rebin(8);
  hK0J_Ropes->Scale(0.125);
  hK0I_Ropes->Rebin(8);
  hK0I_Ropes->Scale(0.125);
  hK0J_Monash->Rebin(8);
  hK0J_Monash->Scale(0.125);
  hK0I_Monash->Rebin(8);
  hK0I_Monash->Scale(0.125);
  //
  hKaonJ_Ropes->Rebin(8);
  hKaonJ_Ropes->Scale(0.125);
  hKaonI_Ropes->Rebin(8);
  hKaonI_Ropes->Scale(0.125);
  hKaonJ_Monash->Rebin(8);
  hKaonJ_Monash->Scale(0.125);
  hKaonI_Monash->Rebin(8);
  hKaonI_Monash->Scale(0.125);
/////////////////////////////////////////////////////////////////////////////
//Ready to draw on canvas
  TCanvas* c = new TCanvas("c1","c1",800,1080);
  auto l = new TLegend(0.354637,0.137827,0.997494,0.28169);
  myLegendSetUp(l);
   c->SetLogy();
   hKaonHM->Scale(100.0);
   hKaonI->Scale(100.0);
   hKaonJ->Scale(100.0);
   hKaonHM_Monash->Scale(100);
   hKaonI_Monash->Scale(100);
   hKaonJ_Monash->Scale(100);
   hKaonHM_Ropes->Scale(100);
   hKaonI_Ropes->Scale(100);
   hKaonJ_Ropes->Scale(100);

   hK0HM->Scale(15.0);
   hK0I->Scale(15.0);
   hK0J->Scale(15.0);
   hK0HM_Monash->Scale(15);
   hK0I_Monash->Scale(15);
   hK0J_Monash->Scale(15);
   hK0HM_Ropes->Scale(15);
   hK0I_Ropes->Scale(15);
   hK0J_Ropes->Scale(15);
  
  hKstarI->GetYaxis()->SetTitle("Ratio of yields to (#pi^{+} #pi^{-})");
  hKstarI->GetXaxis()->SetTitle("p_{T}");
  hKstarI->GetYaxis()->SetRangeUser(0.001,160.0);
  hKstarI->GetXaxis()->SetRangeUser(0.2,8.0);
  hKstarHM->GetYaxis()->SetTitle("Ratio of yields to (#pi^{+} #pi^{-})");
  hKstarHM->GetXaxis()->SetTitle("p_{T}");
  hKstarHM->GetYaxis()->SetRangeUser(0.001,160.0);
  hKstarHM->GetXaxis()->SetRangeUser(0.2,8.0);
  hKstarHM->SetTitle("");

  hKstarI->Draw();
   hKstarI_Monash->Draw("SAME HIST L");
   hKstarI_Ropes->Draw("SAME HIST L");

   hKstarHM->Draw("SAME");
   hKstarHM_Ropes->Draw("SAME HIST L");    
   hKstarHM_Monash->Draw("SAME HIST L");
  
   hKstarJ->Draw("SAME");
   hKstarJ_Monash->Draw("SAME HIST L");
   hKstarJ_Ropes->Draw("SAME HIST L");

   hKaonHM->Draw("SAME");
   hKaonHM_Monash->Draw("SAME HIST L");
   hKaonHM_Ropes->Draw("SAME HIST L");
   hKaonJ->Draw("SAME");
   hKaonJ_Monash->Draw("SAME HIST L");
   hKaonJ_Ropes->Draw("SAME HIST L");
   hKaonI->Draw("SAME");
   hKaonI_Monash->Draw("SAME HIST L");
   hKaonI_Ropes->Draw("SAME HIST L");

   hK0HM->Draw("SAME");
   hK0HM_Monash->Draw("SAME HIST L");
   hK0HM_Ropes->Draw("SAME HIST L");
   hK0J->Draw("SAME");
   hK0J_Monash->Draw("SAME HIST L");
   hK0J_Ropes->Draw("SAME HIST L");
   hK0I->Draw("SAME");
   hK0I_Monash->Draw("SAME HIST L");
   hK0I_Ropes->Draw("SAME HIST L");

  
   //Add all of the text 
   l->SetTextFont(62);
   l->SetHeader("ALICE","C");
   l->SetTextFont(42);
   l->AddEntry(hKstarHM, "S_{0}^{p_{T}=1} Integrated", "lep");
   l->AddEntry(hKstarI, "Isotropic (0-1%)", "lep");
   l->AddEntry(hKstarJ, "Jetty (0-1%)", "lep");
   l->AddEntry(hKstarHM_Monash, "PYTHIA 8.2 Monash", "lep");
   l->AddEntry(hKstarHM_Ropes, "PYTHIA 8.2 Ropes", "lep");
   l->Draw("SAME");

   TLatex TL;
   TL.DrawLatex(6.56892,0.644774, "K^{*}");
   TL.DrawLatex(6.56892,6.48791, "K^{0}_{s} (x15)");
   TL.DrawLatex(6.56892,100.893, "K (x100)");
 
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//END OF MAIN CODE
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void myLegendSetUp(TLegend *currentLegend,float currentTextSize){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}
void ShapeMyHist(TH1D **hist){
  TH1D *htemp = *hist;
  htemp->SetLineWidth(4);
  htemp->SetLabelSize(0.035,"xyz");
  htemp->SetLabelOffset(0.01,"y");
  htemp->SetLabelOffset(0.01,"x");
  htemp->SetLabelColor(kBlack,"xyz");
  htemp->SetTitleSize(0.039,"xyz");
  htemp->SetTitleOffset(1.2,"y");
  htemp->SetTitleOffset(1.2,"x");
  htemp->SetLabelFont(62,"xy");
  htemp->SetTitleFont(62,"xy");
  htemp->GetYaxis()->SetNdivisions(505);
}


void ShapeMyHistMonash(TH1D **hist){
  TH1D *htemp = *hist;
  htemp->SetLineWidth(2);
  htemp->SetLabelSize(0.035,"xyz");
  htemp->SetLabelOffset(0.01,"y");
  htemp->SetLabelOffset(0.01,"x");
  htemp->SetLabelColor(kBlack,"xyz");
  htemp->SetTitleSize(0.039,"xyz");
  htemp->SetTitleOffset(1.4,"y");
  htemp->SetTitleOffset(1.2,"x");
  htemp->SetLabelFont(62,"xy");
  htemp->SetTitleFont(62,"xy");
  htemp->GetYaxis()->SetNdivisions(505);

}
void ShapeMyHistRopes(TH1D **hist){
  TH1D *htemp = *hist;
  htemp->SetLineWidth(2);
  htemp->SetLineStyle(9);
  htemp->SetLabelSize(0.035,"xyz");
  htemp->SetLabelOffset(0.01,"y");
  htemp->SetLabelOffset(0.01,"x");
  htemp->SetLabelColor(kBlack,"xyz");
  htemp->SetTitleSize(0.039,"xyz");
  htemp->SetTitleOffset(1.4,"y");
  htemp->SetTitleOffset(1.2,"x");
  htemp->SetLabelFont(62,"xy");
  htemp->SetTitleFont(62,"xy");
  htemp->GetYaxis()->SetNdivisions(505);

}
void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
        currentPad->SetLeftMargin(currentLeft);
        currentPad->SetTopMargin(currentTop);
        currentPad->SetRightMargin(currentRight);
        currentPad->SetBottomMargin(currentBottom);
	currentPad->SetTickx(1);
	currentPad->SetTicky(1);
        return;
}
void myOptions(){
  gStyle->SetLineWidth(3);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(62);
  gStyle->SetStatFont(62);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(62,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(62,"xyz");
  gStyle->SetTitleOffset(1.0,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetMarkerSize(1);
  gStyle->SetPalette(1);
}
