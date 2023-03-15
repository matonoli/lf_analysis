#include <iostream>
#include "MakeCanvas.C"

using namespace std;


//// UNCOMMENT A LINE TO CHOOSE THE MULTIPLICITY SELECTION


// COMMON PATH FOR INPUTS
TString path(".");

// HELPER GLOBAL VARIABLES
Int_t colors[5] = {kBlack, kRed, kBlack, kGreen+2, kBlue};
const int NRT = 4; 
enum {Stat, Syst, SystUnc, Monash, Ropes, NHIST};
const char* strH[NHIST] = {"Stat", "Syst", "SystUnc"};//, "Monash", "Ropes"};


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
  else gStyle->SetPalette(kRainBow);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
//  gStyle->SetHistLineWidth(1);
  //gStyle->SetHistLineColor(kRed);
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
  h->SetMarkerSize(0.9);
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


void DrawHistograms(TH1D** h, Int_t nhist) {

  
  Int_t lineWidth = 2;

  for (int iH = 0; iH < nhist; ++iH) {
    
    MakeNiceHistogram(h[iH],colors[iH]);
    h[iH]->Draw("E X0 hist same");
  }

}

TH1D* MakeRatioHist(TH1D* hn, TH1D* hd) {

  TH1D* hr = (TH1D*)hn->Clone(Form("hr_%s",hn->GetName()));
  hr->Divide(hd);

  return hr;
}

/*void DrawHistogramsRatios(TH1D** hn, TH1D** hd, Int_t nhist) {

  TH1D* hr[nhist];
  for (int iH = 1; iH < nhist; ++iH) {
    
    hr[iH] = MakeRatioHist(hn[iH][Stat],hd[Stat]);
    hr[iH]->Draw("same");
  }
}*/



void plotPtvRt_Pt(Int_t particle = 1) {

  TString part; TString partL;
  if (particle==1) {
    part = TString("K0s" );
    partL = TString("K^{0}_{S}");
  } else if (particle==2) {
    part = TString("L" );
    partL = TString("#Lambda");
  } else if (particle==3) {
    part = TString("Lbar" );
    partL = TString("#bar{#Lambda}");
  } else {
    printf("Wrong input \n");
    return;
  }
  TString dmc("D");

  TH1::SetDefaultSumw2(1);

  

  // DEFINE WORKING HISTOGRAMS HERE
  TH2D* hNtvPtT;
  TH1D* hNtT[17];
  TH1D* hPtT[65][NHIST];
  TH1D* hNtDT;
  TH1D* hPtDT[NHIST];
  TH2D* hNtvPtTMin;
  TH1D* hNtTMin[17];
  TH1D* hPtTMin[65][NHIST];
  TH1D* hNtDTMin;
  TH1D* hPtDTMin[NHIST];
  TH2D* hNtvPtTMax;
  TH1D* hNtTMax[17];
  TH1D* hPtTMax[65][NHIST];
  TH1D* hNtDTMax;
  TH1D* hPtDTMax[NHIST];
  TH2D* hNtvPtN;
  TH1D* hNtN[17];
  TH1D* hPtN[65][NHIST];
  TH1D* hNtDN;
  TH1D* hPtDN[NHIST];
  TH2D* hNtvPtA;
  TH1D* hNtA[17];
  TH1D* hPtA[65][NHIST];
  TH1D* hNtDA;
  TH1D* hPtDA[NHIST];

  

  // FETCH SOURCE HISTOGRAMS FROM FILES
  //TFile* fIn = new TFile(Form("%s/outMC_230123.root",path.Data()),"READ");
  TFile* fIn = new TFile(Form("%s/outD_230301_0xcsu_hist.root",path.Data()),"READ");

  TDirectoryFile* dirCorr = (TDirectoryFile*)fIn->Get("MyAnalysisV0unfold_4");
  
  for (int i = 0; i<NRT+1; i++) {
    hPtT[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Trans_%i",part.Data(),dmc.Data(),i));
    hPtTMin[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnf_%s_%s_TransMin_%i",part.Data(),dmc.Data(),i));
    hPtTMax[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnf_%s_%s_TransMax_%i",part.Data(),dmc.Data(),i));
    hPtN[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Near_%i",part.Data(),dmc.Data(),i));
    hPtA[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Away_%i",part.Data(),dmc.Data(),i));

    hPtT[i][Syst]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSyst_%s_%s_Trans_%i",part.Data(),dmc.Data(),i));
    hPtTMin[i][Syst] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnfSyst_%s_%s_TransMin_%i",part.Data(),dmc.Data(),i));
    hPtTMax[i][Syst] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnfSyst_%s_%s_TransMax_%i",part.Data(),dmc.Data(),i));
    hPtN[i][Syst]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSyst_%s_%s_Near_%i",part.Data(),dmc.Data(),i));
    hPtA[i][Syst]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSyst_%s_%s_Away_%i",part.Data(),dmc.Data(),i));

    hPtT[i][SystUnc]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSystUnc_%s_%s_Trans_%i",part.Data(),dmc.Data(),i));
    hPtTMin[i][SystUnc] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnfSystUnc_%s_%s_TransMin_%i",part.Data(),dmc.Data(),i));
    hPtTMax[i][SystUnc] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnfSystUnc_%s_%s_TransMax_%i",part.Data(),dmc.Data(),i));
    hPtN[i][SystUnc]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSystUnc_%s_%s_Near_%i",part.Data(),dmc.Data(),i));
    hPtA[i][SystUnc]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSystUnc_%s_%s_Away_%i",part.Data(),dmc.Data(),i));

}


  

  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1000,600);
  C->SetFillStyle(4000);

  TCanvas *C2 = (TCanvas*) gROOT->FindObject("C2");
  if (C2) delete C2;
  C2 = new TCanvas("C2","canvas",1000,600);
  C2->SetFillStyle(4000);

  // Number of PADS
  const int Nx = 3;
  const int Ny = 2;

  // Margins
  float lMargin = 0.03;
  float rMargin = 0.05;
  float bMargin = 0.06;
  float tMargin = 0.05;

  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
  CanvasPartition(C2,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  int latexTextSize = 18; 
  float legendTextSize = 0.06;

  const double minYspectra = 0.000007;
  const double maxYspectra = 11.;

  const double minYratio = 0.1;
  const double maxYratio = 5.2;

  const double minXleft = -0.2;
  const double maxXleft = 8.2;


  /// Left frames
  TH1D* hframeSpectraLeft = new TH1D("hframeSpectraLeft","",100 , minXleft , maxXleft);
  hframeSpectraLeft->GetYaxis()->SetRangeUser(minYspectra,maxYspectra);
  hframeSpectraLeft->GetXaxis()->SetMoreLogLabels(kTRUE);
  hframeSpectraLeft->GetXaxis()->SetNoExponent(kTRUE);

  TH1D* hframeRatioLeft = new TH1D("hframeRatioLeft","",100 , minXleft , maxXleft);
  hframeRatioLeft->GetYaxis()->SetRangeUser(minYratio,maxYratio);
  hframeRatioLeft->GetXaxis()->SetMoreLogLabels(kTRUE);
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
  TuneLatex(latexSp,1.4*latexTextSize);

  TLatex* latexR = new TLatex();
  TuneLatex(latexR,1.2*latexTextSize);

  TLatex* latexSystem = new TLatex();
  TuneLatex(latexSystem,0.9*latexTextSize);

  TLegend* legendLeft = new TLegend(0.04,0.7,0.55,0.93);
  myLegendSetUp(legendLeft,legendTextSize);
  legendLeft->AddEntry(hPtT[1][Stat],"0.0 < R_{T} < 0.8","PL");
  legendLeft->AddEntry(hPtT[2][Stat],"0.8 < R_{T} < 1.5","PL");
  legendLeft->AddEntry(hPtT[3][Stat],"1.5 < R_{T} < 2.5","PL");
  legendLeft->AddEntry(hPtT[4][Stat],"2.5 < R_{T} < 5.0","PL");

  TLegend* legendLeft2 = new TLegend(0.04,0.7,0.55,0.93);
  myLegendSetUp(legendLeft2,legendTextSize);
  legendLeft2->AddEntry(hPtT[1][Stat],"0.0 < R_{T,(-,min,max)} < 0.8","PL");
  legendLeft2->AddEntry(hPtT[2][Stat],"0.8 < R_{T,(-,min,max)} < 1.5","PL");
  legendLeft2->AddEntry(hPtT[3][Stat],"1.5 < R_{T,(-,min,max)} < 2.5","PL");
  legendLeft2->AddEntry(hPtT[4][Stat],"2.5 < R_{T,(-,min,max)} < 5.0","PL");

  TLegend* legendRight = new TLegend(0.8,0.7,0.95,0.85);
  myLegendSetUp(legendRight,legendTextSize);
  //legendRight->AddEntry(hK0soverK[Ref][Monash],"Monash","L");
  //legendRight->AddEntry(hK0soverK[Ref][Ropes],"Ropes","L");

  TPad *pad[Nx][Ny];

  for (Int_t i=0;i<Nx;i++) {
    for (Int_t j=0;j<Ny;j++) {
      C->cd(0);

      // Get the pads previously created.
      char pname[16];
      sprintf(pname,"pad_%i_%i",i,j);
      pad[i][j] = (TPad*) gROOT->FindObject(pname);
      pad[i][j]->Draw();
      pad[i][j]->SetFillStyle(1001);
      if(j) pad[i][j]->SetLogy(kTRUE);
      pad[i][j]->SetLogz(kTRUE);
      
    }
  }

  TPaletteAxis *palette;

  C->cd(0);
  pad[0][1]->Draw();
  pad[0][1]->cd();
  //pad[0][0]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  
  
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtN[iRt][Stat], colors[iRt]);
    MakeNiceHistogram(hPtN[iRt][Syst], colors[iRt]);
    MakeNiceHistogram(hPtN[iRt][SystUnc], colors[iRt]);
    hPtN[iRt][Stat]->Draw("E X0 same");
    hPtN[iRt][Syst]->Draw("e2 same");
  }

  latexSystem->DrawLatex(0.15,0.22,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.15,0.15,"|#eta|<0.8");
  latexSystem->DrawLatex(0.15,0.08,"ALICE Data");

  latexR->DrawLatex(0.20,0.75,"Toward");
  latexSp->DrawLatex(0.65,0.75,Form("#bf{%s}",partL.Data()));

  C->cd(0);
  pad[1][1]->Draw();
  pad[1][1]->cd();
  //pad[0][1]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtA[iRt][Stat], colors[iRt]);
    MakeNiceHistogram(hPtA[iRt][Syst], colors[iRt]);
    MakeNiceHistogram(hPtA[iRt][SystUnc], colors[iRt]);
    hPtA[iRt][Stat]->Draw("E X0 same");
    hPtA[iRt][Syst]->Draw("E2 same");
  }
  
  latexR->DrawLatex(0.2,0.75,"Away");

  C->cd(0);
  pad[2][1]->Draw();
  pad[2][1]->cd();
  //pad[0][2]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtT[iRt][Stat], colors[iRt]);
    MakeNiceHistogram(hPtT[iRt][Syst], colors[iRt]);
    MakeNiceHistogram(hPtT[iRt][SystUnc], colors[iRt]);
    hPtT[iRt][Stat]->Draw("E X0 same");
    hPtT[iRt][Syst]->Draw("E2 same");
  }
  latexR->DrawLatex(0.2,0.75,"Transverse");

  C->cd(0);
  pad[0][0]->Draw();
  pad[0][0]->cd();
  hframeRatioLeft->Draw();
  TH1D* hr[NRT];
  TH1D* hrs[NRT];
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtN[iRt][Stat],hPtN[0][Stat]);
    hr[iRt]->Draw("E X0 same");

    hrs[iRt] = MakeRatioHist(hPtN[iRt][SystUnc],hPtN[0][SystUnc]);
    hrs[iRt]->Draw("E2 same");
  }

  C->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtA[iRt][Stat],hPtA[0][Stat]);
    hr[iRt]->Draw("E X0 same");

    hrs[iRt] = MakeRatioHist(hPtA[iRt][SystUnc],hPtA[0][SystUnc]);
    hrs[iRt]->Draw("E2 same");
  }
  legendLeft->Draw();

  C->cd(0);
  pad[2][0]->Draw();
  pad[2][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtT[iRt][Stat],hPtT[0][Stat]);
    hr[iRt]->Draw("E X0 same");

    hrs[iRt] = MakeRatioHist(hPtT[iRt][SystUnc],hPtT[0][SystUnc]);
    hrs[iRt]->Draw("E2 same");
  }


  C->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX = "#it{p}_{T} (GeV/#it{c})";
  latexTitleX->DrawLatex(0.80, 0.65, TitleX);

   C->cd(0);
  padTitleY1->Draw();
  padTitleY1->cd();
  const char* TitleY1 = "Ratio to R_{T} > 0";
  latexTitleY1->SetTextAngle(90);
  latexTitleY1->DrawLatex(0.6, 0.27, TitleY1);

  C->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY2 = "1/#it{N}_{ev} d^{2}#it{N}/d#it{y}d#it{p}_{T}";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.6, 0.05, TitleY2);

  C->SaveAs(Form("./PtvRt_Pt_%s.pdf",part.Data()));
  C->SaveAs(Form("./PtvRt_Pt_%s.png",part.Data()));


  ///////////////////
  // SECOND CANVAS

  for (Int_t i=0;i<Nx;i++) {
    for (Int_t j=0;j<Ny;j++) {
      C2->cd(0);

      // Get the pads previously created.
      char pname[16];
      sprintf(pname,"pad_%i_%i",i,j);
      pad[i][j] = (TPad*) gROOT->FindObject(pname);
      pad[i][j]->Draw();
      pad[i][j]->SetFillStyle(1001);
      if(j) pad[i][j]->SetLogy(kTRUE);
      pad[i][j]->SetLogz(kTRUE);
      
    }
  }

  C2->cd(0);
  pad[0][1]->Draw();
  pad[0][1]->cd();
  //pad[0][0]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  
  
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtT[iRt][Stat], colors[iRt]);
    MakeNiceHistogram(hPtT[iRt][Syst], colors[iRt]);
    hPtT[iRt][Stat]->Draw("E X0 same");
    hPtT[iRt][Syst]->Draw("E2 same");
  }

  latexSystem->DrawLatex(0.15,0.22,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.15,0.15,"|#eta|<0.8");
  latexSystem->DrawLatex(0.15,0.08,"ALICE Data");

  latexR->DrawLatex(0.20,0.75,"Transverse");
  latexSp->DrawLatex(0.65,0.75,Form("#bf{%s}",partL.Data()));

  C2->cd(0);
  pad[1][1]->Draw();
  pad[1][1]->cd();
  //pad[0][1]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtTMin[iRt][Stat], colors[iRt]);
    MakeNiceHistogram(hPtTMin[iRt][Syst], colors[iRt]);
    MakeNiceHistogram(hPtTMin[iRt][SystUnc], colors[iRt]);
    hPtTMin[iRt][Stat]->Draw("E X0 same");
    hPtTMin[iRt][Syst]->Draw("E2 same");
  }
  
  latexR->DrawLatex(0.2,0.75,"Trans., min");

  C2->cd(0);
  pad[2][1]->Draw();
  pad[2][1]->cd();
  //pad[0][2]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtTMax[iRt][Stat], colors[iRt]);
    MakeNiceHistogram(hPtTMax[iRt][Syst], colors[iRt]);
    MakeNiceHistogram(hPtTMax[iRt][SystUnc], colors[iRt]);
    hPtTMax[iRt][Stat]->Draw("E X0 same");
    hPtTMax[iRt][Syst]->Draw("E2 same");
  }
  latexR->DrawLatex(0.2,0.75,"Trans., max");

  C2->cd(0);
  pad[0][0]->Draw();
  pad[0][0]->cd();
  hframeRatioLeft->Draw();
  
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtT[iRt][Stat],hPtT[0][Stat]);
    hr[iRt]->Draw("E X0 same");

    hrs[iRt] = MakeRatioHist(hPtT[iRt][SystUnc],hPtT[0][SystUnc]);
    hrs[iRt]->Draw("E2 same");
  }

  C2->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtTMin[iRt][Stat],hPtTMin[0][Stat]);
    hr[iRt]->Draw("E X0 same");

    hrs[iRt] = MakeRatioHist(hPtTMin[iRt][SystUnc],hPtTMin[0][SystUnc]);
    hrs[iRt]->Draw("E2 same");
  }
  legendLeft2->Draw();

  C2->cd(0);
  pad[2][0]->Draw();
  pad[2][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtTMax[iRt][Stat],hPtTMax[0][Stat]);
    hr[iRt]->Draw("E X0 same");

    hrs[iRt] = MakeRatioHist(hPtTMax[iRt][SystUnc],hPtTMax[0][SystUnc]);
    hrs[iRt]->Draw("E2 same");
  }


  C2->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX2 = "#it{p}_{T} (GeV/#it{c})";
  latexTitleX->DrawLatex(0.80, 0.65, TitleX2);

  C2->cd(0);
  padTitleY1->Draw();
  padTitleY1->cd();
  const char* TitleY12 = "Ratio to R_{T} > 0";
  latexTitleY1->SetTextAngle(90);
  latexTitleY1->DrawLatex(0.6, 0.27, TitleY12);

  C2->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY22 = "1/#it{N}_{ev} d^{2}#it{N}/d#it{y}d#it{p}_{T}";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.6, 0.05, TitleY22);

  C2->SaveAs(Form("./PtvRt_Pt2_%s.pdf",part.Data()));
  C2->SaveAs(Form("./PtvRt_Pt2_%s.png",part.Data()));


}

