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

void DrawBaryon(){
//LOAD EVERYTHING
/////////////////////////////////////////////////////////////////////////////
  myOptions();
  TFile* Phis = TFile::Open("CorrectedCL1Ratio.root");
  TH1D* hPhiHM = (TH1D*)Phis->Get("hPhiOverPi_HM_Ratio_stat");
  TH1D* hPhiJ = (TH1D*)Phis->Get("hPhiOverPi_Jetty_Ratio_stat");
  TH1D* hPhiI = (TH1D*)Phis->Get("hPhiOverPi_Iso_Ratio_stat");
  TFile* Xis = TFile::Open("xi_results_So_CL1_top1_aug_6.root");
  TH1D* hXiHM = (TH1D*)Xis->Get("hXiOverPi_Ref_Ratio_stat");
  TH1D* hXiJ = (TH1D*)Xis->Get("hXiOverPi_Jetty10_Ratio_stat");
  TH1D* hXiI = (TH1D*)Xis->Get("hXiOverPi_Iso10_Ratio_stat");
  TFile* Protons = TFile::Open("pi_k_p_spectra.root");
  TH1D* hProtonHM = (TH1D*)Protons->Get("Trks_Reference_Bin_mult_1/hProton_stat");
  TH1D* hProtonJ = (TH1D*)Protons->Get("Trks_Jetty_2_Bin_mult_1/hProton_stat");
  TH1D* hProtonI = (TH1D*)Protons->Get("Trks_Isotropic_2_Bin_mult_1/hProton_stat");
  TFile* Pis = TFile::Open("pi_k_p_spectra.root");
  TH1D* hPiHM = (TH1D*)Pis->Get("Trks_Reference_Bin_mult_1/hPion_stat");
  TH1D* hPiJ = (TH1D*)Pis->Get("Trks_Jetty_2_Bin_mult_1/hPion_stat");
  TH1D* hPiI = (TH1D*)Pis->Get("Trks_Isotropic_2_Bin_mult_1/hPion_stat");
  TFile* Lambdas = TFile::Open("k0s_l_spectra_NCharged01_spher10.root");
  TH1D* hLambdaHM = (TH1D*)Lambdas->Get("hLLbarToPi_Ref_Stat");
  TH1D* hLambdaJ = (TH1D*)Lambdas->Get("hLLbarToPi_Jet_Stat");
  TH1D* hLambdaI = (TH1D*)Lambdas->Get("hLLbarToPi_Iso_Stat");
  
  TFile* Monash = TFile::Open("Complete_Monash.root");
  TFile* Ropes = TFile::Open("Complete_Ropes.root");

  TH1D* hPhiHM_Monash = (TH1D*)Monash->Get("fPhis_HMCL1_01");
  TH1D* hPhiJ_Monash = (TH1D*)Monash->Get("fPhis_CL1_01_Jetty_10");
  TH1D* hPhiI_Monash = (TH1D*)Monash->Get("fPhis_CL1_01_Iso_10");
  TH1D* hPiHM_Monash = (TH1D*)Monash->Get("fPions_HMCL1_01");
  TH1D* hPiJ_Monash = (TH1D*)Monash->Get("fPions_CL1_01_Jetty_10");
  TH1D* hPiI_Monash = (TH1D*)Monash->Get("fPions_CL1_01_Iso_10");
  TH1D* hPhiHM_Ropes = (TH1D*)Ropes->Get("fPhis_HMCL1_01");
  TH1D* hPhiJ_Ropes = (TH1D*)Ropes->Get("fPhis_CL1_01_Jetty_10");
  TH1D* hPhiI_Ropes = (TH1D*)Ropes->Get("fPhis_CL1_01_Iso_10");
  TH1D* hPiHM_Ropes = (TH1D*)Ropes->Get("fPions_HMCL1_01");
  TH1D* hPiJ_Ropes = (TH1D*)Ropes->Get("fPions_CL1_01_Jetty_10");
  TH1D* hPiI_Ropes = (TH1D*)Ropes->Get("fPions_CL1_01_Iso_10");
  TH1D* hXiHM_Ropes = (TH1D*)Ropes->Get("fXis_HMCL1_01");
  TH1D* hXiJ_Ropes = (TH1D*)Ropes->Get("fXis_CL1_01_Jetty_10");
  TH1D* hXiI_Ropes = (TH1D*)Ropes->Get("fXis_CL1_01_Iso_10");
  TH1D* hXiHM_Monash = (TH1D*)Monash->Get("fXis_HMCL1_01");
  TH1D* hXiJ_Monash = (TH1D*)Monash->Get("fXis_CL1_01_Jetty_10");
  TH1D* hXiI_Monash = (TH1D*)Monash->Get("fXis_CL1_01_Iso_10");
  TH1D* hProtonHM_Ropes = (TH1D*)Ropes->Get("fProtons_HMCL1_01");
  TH1D* hProtonJ_Ropes = (TH1D*)Ropes->Get("fProtons_CL1_01_Jetty_10");
  TH1D* hProtonI_Ropes = (TH1D*)Ropes->Get("fProtons_CL1_01_Iso_10");
  TH1D* hProtonHM_Monash = (TH1D*)Monash->Get("fProtons_HMCL1_01");
  TH1D* hProtonJ_Monash = (TH1D*)Monash->Get("fProtons_CL1_01_Jetty_10");
  TH1D* hProtonI_Monash = (TH1D*)Monash->Get("fProtons_CL1_01_Iso_10");
  TH1D* hLambdaHM_Ropes = (TH1D*)Ropes->Get("fLambdas_HMCL1_01");
  TH1D* hLambdaJ_Ropes = (TH1D*)Ropes->Get("fLambdas_CL1_01_Jetty_10");
  TH1D* hLambdaI_Ropes = (TH1D*)Ropes->Get("fLambdas_CL1_01_Iso_10");
  TH1D* hLambdaHM_Monash = (TH1D*)Monash->Get("fLambdas_HMCL1_01");
  TH1D* hLambdaJ_Monash = (TH1D*)Monash->Get("fLambdas_CL1_01_Jetty_10");
  TH1D* hLambdaI_Monash = (TH1D*)Monash->Get("fLambdas_CL1_01_Iso_10");
//Done Loading
/////////////////////////////////////////////////////////////////////////////
//Set Line Color, and Enable the Style modifcations for each histo
  hProtonHM->SetLineColor(1);
  hProtonJ->SetLineColor(2);
  hXiHM->SetLineColor(1);
  hXiJ->SetLineColor(2);
  
  hProtonHM->Divide(hPiHM);
  hProtonJ->Divide(hPiJ);
  hProtonI->Divide(hPiI);
  
  ShapeMyHistMonash(&hPhiHM_Monash);
  ShapeMyHistMonash(&hPhiJ_Monash);
  ShapeMyHistMonash(&hPhiI_Monash);
  hPhiHM_Monash->Divide(hPiHM_Monash);
  hPhiJ_Monash->Divide(hPiJ_Monash);
  hPhiI_Monash->Divide(hPiI_Monash);
  hPhiHM_Monash->SetLineColor(1);
  hPhiJ_Monash->SetLineColor(2);
  ShapeMyHistRopes(&hPhiHM_Ropes);
  ShapeMyHistRopes(&hPhiJ_Ropes);
  ShapeMyHistRopes(&hPhiI_Ropes);
  hPhiHM_Ropes->Divide(hPiHM_Ropes);
  hPhiJ_Ropes->Divide(hPiJ_Ropes);
  hPhiI_Ropes->Divide(hPiI_Ropes);
  hPhiHM_Ropes->SetLineColor(1);
  hPhiJ_Ropes->SetLineColor(2);

  ShapeMyHistMonash(&hXiHM_Monash);
  ShapeMyHistMonash(&hXiJ_Monash);
  ShapeMyHistMonash(&hXiI_Monash);
  hXiHM_Monash->Divide(hPiHM_Monash);
  hXiJ_Monash->Divide(hPiJ_Monash);
  hXiI_Monash->Divide(hPiI_Monash);
  hXiHM_Monash->SetLineColor(1);
  hXiJ_Monash->SetLineColor(2);
  ShapeMyHistRopes(&hXiHM_Ropes);
  ShapeMyHistRopes(&hXiJ_Ropes);
  ShapeMyHistRopes(&hXiI_Ropes);
  hXiHM_Ropes->Divide(hPiHM_Ropes);
  hXiJ_Ropes->Divide(hPiJ_Ropes);
  hXiI_Ropes->Divide(hPiI_Ropes);
  hXiHM_Ropes->SetLineColor(1);
  hXiJ_Ropes->SetLineColor(2);

  ShapeMyHistMonash(&hProtonHM_Monash);
  ShapeMyHistMonash(&hProtonJ_Monash);
  ShapeMyHistMonash(&hProtonI_Monash);
  hProtonHM_Monash->Divide(hPiHM_Monash);
  hProtonJ_Monash->Divide(hPiJ_Monash);
  hProtonI_Monash->Divide(hPiI_Monash);
  hProtonHM_Monash->SetLineColor(1);
  hProtonJ_Monash->SetLineColor(2);
  ShapeMyHistRopes(&hProtonHM_Ropes);
  ShapeMyHistRopes(&hProtonJ_Ropes);
  ShapeMyHistRopes(&hProtonI_Ropes);
  hProtonHM_Ropes->Divide(hPiHM_Ropes);
  hProtonJ_Ropes->Divide(hPiJ_Ropes);
  hProtonI_Ropes->Divide(hPiI_Ropes);
  hProtonHM_Ropes->SetLineColor(1);
  hProtonJ_Ropes->SetLineColor(2);
 
  ShapeMyHistMonash(&hLambdaHM_Monash);
  ShapeMyHistMonash(&hLambdaJ_Monash);
  ShapeMyHistMonash(&hLambdaI_Monash);
  hLambdaHM_Monash->Divide(hPiHM_Monash);
  hLambdaJ_Monash->Divide(hPiJ_Monash);
  hLambdaI_Monash->Divide(hPiI_Monash);
  hLambdaHM_Monash->SetLineColor(1);
  hLambdaJ_Monash->SetLineColor(2);
  ShapeMyHistRopes(&hLambdaHM_Ropes);
  ShapeMyHistRopes(&hLambdaJ_Ropes);
  ShapeMyHistRopes(&hLambdaI_Ropes);
  hLambdaHM_Ropes->Divide(hPiHM_Ropes);
  hLambdaJ_Ropes->Divide(hPiJ_Ropes);
  hLambdaI_Ropes->Divide(hPiI_Ropes);
  hLambdaHM_Ropes->SetLineColor(1);
  hLambdaJ_Ropes->SetLineColor(2);

  ShapeMyHist(&hPhiHM);
  ShapeMyHist(&hPhiJ);
  ShapeMyHist(&hPhiI);  
  ShapeMyHist(&hProtonHM);
  ShapeMyHist(&hProtonJ);
  ShapeMyHist(&hProtonI);
  ShapeMyHist(&hXiHM);
  ShapeMyHist(&hXiJ);
  ShapeMyHist(&hXiI);
  ShapeMyHist(&hLambdaHM);
  ShapeMyHist(&hLambdaJ);
  ShapeMyHist(&hLambdaI);
  
//Histos prep'd, and ready to go
/////////////////////////////////////////////////////////////////////////////
//Rebin Monte-Carlo
  hPhiJ_Ropes->Rebin(8);
  hPhiJ_Ropes->Scale(0.125);
  hPhiI_Ropes->Rebin(8);
  hPhiI_Ropes->Scale(0.125);
  hPhiJ_Monash->Rebin(8);
  hPhiJ_Monash->Scale(0.125);
  hPhiI_Monash->Rebin(8);
  hPhiI_Monash->Scale(0.125);
  //
  hXiJ_Ropes->Rebin(8);
  hXiJ_Ropes->Scale(0.125);
  hXiI_Ropes->Rebin(8);
  hXiI_Ropes->Scale(0.125);
  hXiJ_Monash->Rebin(8);
  hXiJ_Monash->Scale(0.125);
  hXiI_Monash->Rebin(8);
  hXiI_Monash->Scale(0.125);
  //
  hLambdaJ_Ropes->Rebin(8);
  hLambdaJ_Ropes->Scale(0.125);
  hLambdaI_Ropes->Rebin(8);
  hLambdaI_Ropes->Scale(0.125);
  hLambdaJ_Monash->Rebin(8);
  hLambdaJ_Monash->Scale(0.125);
  hLambdaI_Monash->Rebin(8);
  hLambdaI_Monash->Scale(0.125);
  //
  hProtonJ_Ropes->Rebin(8);
  hProtonJ_Ropes->Scale(0.125);
  hProtonI_Ropes->Rebin(8);
  hProtonI_Ropes->Scale(0.125);
  hProtonJ_Monash->Rebin(8);
  hProtonJ_Monash->Scale(0.125);
  hProtonI_Monash->Rebin(8);
  hProtonI_Monash->Scale(0.125);
/////////////////////////////////////////////////////////////////////////////
//Ready to draw on canvas
  TCanvas* c = new TCanvas("c1","c1",800,1080);
  auto l = new TLegend(0.354637,0.137827,0.997494,0.28169);
  myLegendSetUp(l);
   c->SetLogy();
  hXiHM->Scale(30);
  hXiI->Scale(30);
  hXiJ->Scale(30);
  hXiHM_Monash->Scale(30);
  hXiI_Monash->Scale(30);
  hXiJ_Monash->Scale(30);
  hXiHM_Ropes->Scale(30);
  hXiI_Ropes->Scale(30);
  hXiJ_Ropes->Scale(30);
  
  hProtonHM->Scale(250.0);
  hProtonI->Scale(250.0);
  hProtonJ->Scale(250.0);
  hProtonHM_Monash->Scale(250);
  hProtonI_Monash->Scale(250);
  hProtonJ_Monash->Scale(250);
  hProtonHM_Ropes->Scale(250);
  hProtonI_Ropes->Scale(250);
  hProtonJ_Ropes->Scale(250);

  hLambdaHM->Scale(75.0);
  hLambdaI->Scale(75.0);
  hLambdaJ->Scale(75.0);
  hLambdaHM_Monash->Scale(75);
  hLambdaI_Monash->Scale(75);
  hLambdaJ_Monash->Scale(75);
  hLambdaHM_Ropes->Scale(75);
  hLambdaI_Ropes->Scale(75);
  hLambdaJ_Ropes->Scale(75);
  
  hPhiI->GetYaxis()->SetTitle("Ratio of yields to (#pi^{+} #pi^{-})");
  hPhiI->GetXaxis()->SetTitle("p_{T}");
  hPhiI->GetYaxis()->SetRangeUser(0.001,160.0);
  hPhiI->GetXaxis()->SetRangeUser(0.2,8.0);
  hPhiHM->GetYaxis()->SetTitle("Ratio of yields to (#pi^{+} #pi^{-})");
  hPhiHM->GetXaxis()->SetTitle("p_{T}");
  hPhiHM->GetYaxis()->SetRangeUser(0.001,160.0);
  hPhiHM->GetXaxis()->SetRangeUser(0.2,8.0);
  hPhiHM->SetTitle("");

  hPhiI->Draw();
  hPhiI_Monash->Draw("SAME HIST L");
  hPhiI_Ropes->Draw("SAME HIST L");

  hPhiHM->Draw("SAME");
  hPhiHM_Ropes->Draw("SAME HIST L");    
  hPhiHM_Monash->Draw("SAME HIST L");
  
  hPhiJ->Draw("SAME");
  hPhiJ_Monash->Draw("SAME HIST L");
  hPhiJ_Ropes->Draw("SAME HIST L");

  hProtonHM->Draw("SAME");
  hProtonHM_Monash->Draw("SAME HIST L");
  hProtonHM_Ropes->Draw("SAME HIST L");
  hProtonJ->Draw("SAME");
  hProtonJ_Monash->Draw("SAME HIST L");
  hProtonJ_Ropes->Draw("SAME HIST L");
  hProtonI->Draw("SAME");
  hProtonI_Monash->Draw("SAME HIST L");
  hProtonI_Ropes->Draw("SAME HIST L");

  hXiHM->SetMinimum(0.00000001);
  hXiJ->SetMinimum(0.00000001);
  hXiI->SetMinimum(0.00000001);
  hXiHM->Draw("SAME");
  hXiHM_Monash->Draw("SAME HIST L");
  hXiHM_Ropes->Draw("SAME HIST L");
  hXiJ->Draw("SAME");
  hXiJ_Monash->Draw("SAME HIST L");
  hXiJ_Ropes->Draw("SAME HIST L");
  hXiI->Draw("SAME");
  hXiI_Monash->Draw("SAME HIST L");
  hXiI_Ropes->Draw("SAME HIST L");

  hLambdaHM->Draw("SAME");
  hLambdaHM_Monash->Draw("SAME HIST L");
  hLambdaHM_Ropes->Draw("SAME HIST L");
  hLambdaJ->Draw("SAME");
  hLambdaJ_Monash->Draw("SAME HIST L");
  hLambdaJ_Ropes->Draw("SAME HIST L");
  hLambdaI->Draw("SAME");
  hLambdaI_Monash->Draw("SAME HIST L");
  hLambdaI_Ropes->Draw("SAME HIST L");

  
  //Add all of the text 
  l->SetTextFont(62);
  l->SetHeader("#it{ALICE}");
  l->SetTextFont(42);
  l->AddEntry(hPhiHM, "S_{0}^{p_{T}=1} Integrated", "lep");
  l->AddEntry(hPhiI, "Isotropic (0-1%)", "lep");
  l->AddEntry(hPhiJ, "Jetty (0-1%)", "lep");
  l->AddEntry(hPhiHM_Monash, "PYTHIA 8.2 Monash", "lep");
  l->AddEntry(hPhiHM_Ropes, "PYTHIA 8.2 Ropes", "lep");
  l->Draw("SAME");

  TLatex TL;
  TL.SetTextSize(0.0302415);
  TL.DrawLatex(6.56892,0.0886377, "#phi");
  TL.SetTextSize(0.0402415);
  TL.DrawLatex(6.56892,1.57086, "#Xi (x30)");
  TL.DrawLatex(6.56892,17.0459, "#Lambda (x75)");
  TL.DrawLatex(6.56892,100.893, "p (x250)");
 
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
