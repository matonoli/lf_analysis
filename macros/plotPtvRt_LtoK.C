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
const char* strH[NHIST] = {"Stat", "Syst", "SystUnc", "Monash", "Ropes"};


// PLOTTING FUNCTIONS
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
  TString hrname(hr->GetName());
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
        //hframe->GetXaxis()->SetTitle(titleX);
        
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



void plotPtvRt_LtoK() {

  //TH1::SetDefaultSumw2(1);

  

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

  TH1D* hKPtT[65][NHIST];
  TH1D* hKPtTMin[65][NHIST];
  TH1D* hKPtTMax[65][NHIST];
  TH1D* hKPtN[65][NHIST];
  TH1D* hKPtA[65][NHIST];
  

  // FETCH SOURCE HISTOGRAMS FROM FILES
  //TFile* fIn = new TFile(Form("%s/outMC_230123.root",path.Data()),"READ");
  TFile* fIn = new TFile(Form("%s/outD_230220_0xcu_hist.root",path.Data()),"READ");
  //Bool_t isRC = true; Bool_t isRTMin = false; TString part("K0s" ); TString partL("K^{0}_{S}"); TString dmc("D");
  Bool_t isRC = true; Bool_t isRTMin = false; TString part("Lbar" ); TString partL("#Lambda");TString dmc("D");
  //Bool_t isRC = true; Bool_t isRTMin = false; TString part("Pion" ); TString partL("#pi^{#pm}");
  //Bool_t isRC = true; Bool_t isRTMin = false; TString part("Kaon" ); TString partL("K^{#pm}");

  TDirectoryFile* dirCorr = (TDirectoryFile*)fIn->Get("MyAnalysisV0unfold_3");

  for (int i = 0; i<NRT+1; i++) {
    hPtT[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Trans_%i",part.Data(),dmc.Data(),i));
    hPtTMin[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnf_%s_%s_TransMin_%i",part.Data(),dmc.Data(),i));
    hPtTMax[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnf_%s_%s_TransMax_%i",part.Data(),dmc.Data(),i));
    hPtN[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Near_%i",part.Data(),dmc.Data(),i));
    hPtA[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Away_%i",part.Data(),dmc.Data(),i));

    cout << "h is " << hPtT[i][Stat] << endl;
  }

  part = TString("K0s"); partL = TString("#Lambda / K^{0}_{S}");
  for (int i = 0; i<NRT+1; i++) {
    hKPtT[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Trans_%i",part.Data(),dmc.Data(),i));
    hKPtTMin[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnf_%s_%s_TransMin_%i",part.Data(),dmc.Data(),i));
    hKPtTMax[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnf_%s_%s_TransMax_%i",part.Data(),dmc.Data(),i));
    hKPtN[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Near_%i",part.Data(),dmc.Data(),i));
    hKPtA[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Away_%i",part.Data(),dmc.Data(),i));

    cout << "h is " << hPtT[i][Stat] << endl;
  }

  // PROCESS HISTOGRAMS
  for (int i = 0; i<NRT+1; i++) {
    hPtT[i][Stat]->Divide(hKPtT[i][Stat]);
    hPtTMin[i][Stat]->Divide(hKPtTMin[i][Stat]);
    hPtTMax[i][Stat]->Divide(hKPtTMax[i][Stat]);
    hPtN[i][Stat]->Divide(hKPtN[i][Stat]);
    hPtA[i][Stat]->Divide(hKPtA[i][Stat]);
    cout << "h is " << hPtT[i][Stat] << endl;
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

  const double minYspectra = 0.01;
  const double maxYspectra = 1.01;

  const double minYratio = 0.35;
  const double maxYratio = 2.05;

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
      //if(j) pad[i][j]->SetLogy(kTRUE);
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
    hPtN[iRt][Stat]->Draw("same");
  }

  
  latexR->DrawLatex(0.20,0.75,"Toward");
  latexSp->DrawLatex(0.65,0.75,Form("#bf{%s}",partL.Data()));

  C->cd(0);
  pad[1][1]->Draw();
  pad[1][1]->cd();
  //pad[0][1]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtA[iRt][Stat], colors[iRt]);
    hPtA[iRt][Stat]->Draw("same");
  }
  
  latexR->DrawLatex(0.2,0.75,"Away");

  C->cd(0);
  pad[2][1]->Draw();
  pad[2][1]->cd();
  //pad[0][2]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtT[iRt][Stat], colors[iRt]);
    hPtT[iRt][Stat]->Draw("same");
  }
  latexR->DrawLatex(0.2,0.75,"Transverse");

  C->cd(0);
  pad[0][0]->Draw();
  pad[0][0]->cd();
  hframeRatioLeft->Draw();
  TH1D* hr[NRT];
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtN[iRt][Stat],hPtN[0][Stat]);
    hr[iRt]->Draw("same");
  }
  latexSystem->DrawLatex(0.15,0.92,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.15,0.85,"|#eta|<0.8");
  latexSystem->DrawLatex(0.15,0.78,"ALICE Data");


  C->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtA[iRt][Stat],hPtA[0][Stat]);
    hr[iRt]->Draw("same");
  }
  legendLeft->Draw();

  C->cd(0);
  pad[2][0]->Draw();
  pad[2][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtT[iRt][Stat],hPtT[0][Stat]);
    hr[iRt]->Draw("same");
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
  const char* TitleY2 = "(#Lambda + #bar{#Lambda}) / 2K^{0}_{S}";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.6, 0.2, TitleY2);

  C->SaveAs(Form("./PtvRt_LtoK_%s.pdf",part.Data()));
  C->SaveAs(Form("./PtvRt_LtoK_%s.png",part.Data()));


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
      //if(j) pad[i][j]->SetLogy(kTRUE);
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
    hPtT[iRt][Stat]->Draw("same");
  }


  latexR->DrawLatex(0.20,0.75,"Transverse");
  latexSp->DrawLatex(0.65,0.75,Form("#bf{%s}",partL.Data()));

  C2->cd(0);
  pad[1][1]->Draw();
  pad[1][1]->cd();
  //pad[0][1]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtTMin[iRt][Stat], colors[iRt]);
    hPtTMin[iRt][Stat]->Draw("same");
  }
  
  latexR->DrawLatex(0.2,0.75,"Trans., min");

  C2->cd(0);
  pad[2][1]->Draw();
  pad[2][1]->cd();
  //pad[0][2]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {
    MakeNiceHistogram(hPtTMax[iRt][Stat], colors[iRt]);
    hPtTMax[iRt][Stat]->Draw("same");
  }
  latexR->DrawLatex(0.2,0.75,"Trans., max");

  C2->cd(0);
  pad[0][0]->Draw();
  pad[0][0]->cd();
  hframeRatioLeft->Draw();
  
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtT[iRt][Stat],hPtT[0][Stat]);
    hr[iRt]->Draw("same");
  }

  latexSystem->DrawLatex(0.15,0.92,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.15,0.85,"|#eta|<0.8");
  latexSystem->DrawLatex(0.15,0.78,"ALICE Data");

  C2->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtTMin[iRt][Stat],hPtTMin[0][Stat]);
    hr[iRt]->Draw("same");
  }
  legendLeft2->Draw();

  C2->cd(0);
  pad[2][0]->Draw();
  pad[2][0]->cd();
  hframeRatioLeft->Draw();
  for (int iRt = NRT; iRt > 0; --iRt)  {    
    hr[iRt] = MakeRatioHist(hPtTMax[iRt][Stat],hPtTMax[0][Stat]);
    hr[iRt]->Draw("same");
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
  const char* TitleY22 = "(#Lambda + #bar{#Lambda}) / 2K^{0}_{S}";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.6, 0.2, TitleY22);

  C2->SaveAs(Form("./PtvRt_LtoK2_%s.pdf",part.Data()));
  C2->SaveAs(Form("./PtvRt_LtoK2_%s.png",part.Data()));


}

