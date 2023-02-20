#include <iostream>
#include "MakeCanvas.C"

using namespace std;


//// UNCOMMENT A LINE TO CHOOSE THE MULTIPLICITY SELECTION


// COMMON PATH FOR INPUTS
TString path(".");

// HELPER GLOBAL VARIABLES

Int_t colors[3] = {kRed, kBlue, kBlack};



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

  //h->SetLineColor(col);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.9);
  //h->SetMarkerColor(col);
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

void DrawHistogramsRatios(TH1D** hn, TH1D** hd, Int_t nhist, Int_t col) {

  TH1D* hr[nhist];
  for (int iH = 0; iH < nhist; ++iH) {
    
    hr[iH] = MakeRatioHist(hn[iH],hd[iH]);
    hr[iH]->Draw("E X0 same");
  }
}



void plotPtvNt_Yield() {

  TH1::SetDefaultSumw2(1);

  // DEFINE WORKING HISTOGRAMS HERE
  TH2D* hNtvPtT;
  TH1D* hNtT[17];
  TH1D* hPtT[65];
  TH1D* hNtDT;
  TH1D* hPtDT;
  TH2D* hNtvPtN;
  TH1D* hNtN[17];
  TH1D* hPtN[65];
  TH1D* hNtDN;
  TH1D* hPtDN;
  TH2D* hNtvPtA;
  TH1D* hNtA[17];
  TH1D* hPtA[65];
  TH1D* hNtDA;
  TH1D* hPtDA;

  TH1D* hMeanPtT = new TH1D("hMeanPtT","",50,-0.5,49.5);
  TH1D* hMeanPtN = new TH1D("hMeanPtN","",50,-0.5,49.5);
  TH1D* hMeanPtA = new TH1D("hMeanPtA","",50,-0.5,49.5);
  

  // FETCH SOURCE HISTOGRAMS FROM FILES
  //TFile* fIn = new TFile(Form("%s/outMC_230123.root",path.Data()),"READ");
  TFile* fIn = new TFile(Form("%s/testout2.root",path.Data()),"READ");
  TDirectoryFile* dirCorr = (TDirectoryFile*)fIn->Get("MyAnalysisV0unfold_3");
  Bool_t isRC = true; Bool_t isRTMin = false; TString part("K0s" ); TString partL("K^{0}_{S}");
  //Bool_t isRC = true; Bool_t isRTMin = false; TString part("L" ); TString partL("#Lambda");
  //Bool_t isRC = true; Bool_t isRTMin = false; TString part("Pion" ); TString partL("#pi^{#pm}");
  //Bool_t isRC = true; Bool_t isRTMin = false; TString part("Kaon" ); TString partL("K^{#pm}");

  //hNtvPtT = (TH2D*)dirCorr->Get(Form("hKpmPtNtMinMC_TransMin",part.Data()));
  //hNtvPtN = (TH2D*)dirCorr->Get(Form("hKpmPtNtMinMC_Near",part.Data()));
  //hNtvPtA = (TH2D*)dirCorr->Get(Form("hKpmPtNtMinMC_Away",part.Data()));

  //hNtvPtT = (TH2D*)dirCorr->Get(Form("hPiPtNtMinMC_Trans",part.Data()));
  //hNtvPtN = (TH2D*)dirCorr->Get(Form("hPiPtNtMinMC_Near",part.Data()));
  //hNtvPtA = (TH2D*)dirCorr->Get(Form("hPiPtNtMinMC_Away",part.Data()));

  //hNtvPtT = (TH2D*)dirCorr->Get(Form("hV0PtNtMin_%s_MC_TransMin",part.Data()));
  //hNtvPtN = (TH2D*)dirCorr->Get(Form("hV0PtNtMin_%s_MC_Near",part.Data()));
  //hNtvPtA = (TH2D*)dirCorr->Get(Form("hV0PtNtMin_%s_MC_Away",part.Data()));

  hNtvPtT = (TH2D*)dirCorr->Get(Form("hV0PtNtFitCorrUnf_%s_RC_Trans",part.Data()));
  hNtvPtN = (TH2D*)dirCorr->Get(Form("hV0PtNtFitCorrUnf_%s_RC_Near",part.Data()));
  hNtvPtA = (TH2D*)dirCorr->Get(Form("hV0PtNtFitCorrUnf_%s_RC_Away",part.Data()));
  //hNtvPt = (TH2D*)dirCorr->Get("hPiPtNtRC_Trans");

  dirCorr = (TDirectoryFile*)fIn->Get("MyAnalysisV0correct_2");
  hPtDT = (TH1D*)((TH2D*)dirCorr->Get(Form("hV0PtNtFitCorr_%s_RC_Trans",part.Data())))->ProjectionX();
  hPtDN = (TH1D*)((TH2D*)dirCorr->Get(Form("hV0PtNtFitCorr_%s_RC_Near",part.Data())))->ProjectionX();
  hPtDA = (TH1D*)((TH2D*)dirCorr->Get(Form("hV0PtNtFitCorr_%s_RC_Away",part.Data())))->ProjectionX();
  

  // PROCESS HISTOGRAMS HERE
  Int_t nPt = hNtvPtT->GetNbinsX();
  Int_t nNt = hNtvPtT->GetNbinsY();
  hNtDT = (TH1D*)hNtvPtT->ProjectionY();
  Double_t meanPt = 0; Double_t meanPtErr = 0;
  Int_t blow = 1; Int_t bhigh = 17;

  for (int iNt = 1; iNt < nNt+1; ++iNt) {
    
    hPtT[iNt] = (TH1D*)hNtvPtT->ProjectionX(Form("hPtT_%i",iNt),iNt,iNt);
    meanPt = hPtT[iNt]->Integral(blow,bhigh);
    meanPtErr = 0.1*TMath::Sqrt(meanPt);
    meanPt = meanPt / hPtDT->Integral(blow,bhigh);
    hMeanPtT->SetBinContent(iNt,meanPt);
    hMeanPtT->SetBinError(iNt,meanPtErr);

    hPtN[iNt] = (TH1D*)hNtvPtN->ProjectionX(Form("hPtN_%i",iNt),iNt,iNt);
    meanPt = hPtN[iNt]->Integral(blow,bhigh);
    meanPtErr = 0.1*TMath::Sqrt(meanPt);
    meanPt = meanPt / hPtDN->Integral(blow,bhigh);
    hMeanPtN->SetBinContent(iNt,meanPt);
    hMeanPtN->SetBinError(iNt,meanPtErr);

    hPtA[iNt] = (TH1D*)hNtvPtA->ProjectionX(Form("hPtA_%i",iNt),iNt,iNt);
    meanPt = hPtA[iNt]->Integral(blow,bhigh);
    meanPtErr = 0.1*TMath::Sqrt(meanPt);
    meanPt = meanPt / hPtDA->Integral(blow,bhigh);
    hMeanPtA->SetBinContent(iNt,meanPt);
    hMeanPtA->SetBinError(iNt,meanPtErr);

  }

  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1000,400);
  C->SetFillStyle(4000);

  // Number of PADS
  const int Nx = 3;
  const int Ny = 1;

  // Margins
  float lMargin = 0.03;
  float rMargin = 0.05;
  float bMargin = 0.06;
  float tMargin = 0.05;

  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  int latexTextSize = 18; 
  float legendTextSize = 0.05;

  const double minYratio = 0.1;
  const double maxYratio = 8.2;

  const double minXleft = -0.5;
  const double maxXleft = 21.0;


  /// Left frames
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
  TuneLatex(latexSp,1.4*latexTextSize);

  TLatex* latexR = new TLatex();
  TuneLatex(latexR,1.2*latexTextSize);

  TLatex* latexSystem = new TLatex();
  TuneLatex(latexSystem,0.9*latexTextSize);

  TLegend* legendLeft = new TLegend(0.45,0.71,0.91,0.90);
  myLegendSetUp(legendLeft,legendTextSize);
  //legendLeft->AddEntry(hPhi[is2011],"ITSTPC2011","PL");

  TLegend* legendRight = new TLegend(0.8,0.7,0.95,0.85);
  myLegendSetUp(legendRight,legendTextSize);
  //legendRight->AddEntry(hK0soverK[Ref][Monash],"Monash","L");
  //legendRight->AddEntry(hK0soverK[Ref][Ropes],"Ropes","L");

  TPad *pad[Nx][Ny];

   //TF1 *f2=new TF1("f2","exp(x)",0,2);
   //TGaxis *A2 = new TGaxis(1,1,0,0,"f2");
   //A2->SetTitle("exponential axis");
   //A2->SetLabelSize(0.03);
   //A2->SetTitleSize(0.03);
   //A2->SetTitleOffset(1.2);
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
      pad[i][j]->SetLogz(kTRUE);
      
    }
  }

  TPaletteAxis *palette;

  //new TCanvas;
  //hframeRatioLeft->Draw("colz");
  //return;

  C->cd(0);
  pad[0][0]->Draw();
  pad[0][0]->cd();
  pad[0][0]->SetLogy(kFALSE);
  hframeRatioLeft->Draw();
  
  //MakeNiceHistogram(hNtD, colors[0]); 
  MakeNiceHistogram(hMeanPtN, colors[0]);
  hMeanPtN->Draw("same");
  

  latexSystem->DrawLatex(0.15,0.92,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.15,0.86,"|#eta|<0.8");
  latexSystem->DrawLatex(0.15,0.80,"ALICE MC");

  legendLeft->Draw();
  latexR->DrawLatex(0.15,0.7,"Toward");
  latexSp->DrawLatex(0.65,0.80,Form("#bf{%s}",partL.Data()));

  C->cd(0);
  pad[0][1]->Draw();
  pad[0][1]->cd();
  pad[0][1]->SetLogy(kFALSE);
  hframeRatioLeft->Draw();
  MakeNiceHistogram(hMeanPtA, colors[0]);
  hMeanPtA->Draw("same");

  //latexSystem->DrawLatex(0.04,0.92,"(blue) 0.15 GeV/#it{c} < #it{p}_{T} < 8 GeV/#it{c} (red)");
  latexR->DrawLatex(0.15,0.7,"Away");

  C->cd(0);
  pad[0][2]->Draw();
  pad[0][2]->cd();
  pad[0][2]->SetLogy(kFALSE);
  hframeRatioLeft->Draw();
  MakeNiceHistogram(hMeanPtT, colors[0]);
  hMeanPtT->Draw("same");

  latexR->DrawLatex(0.15,0.7,"Transverse");


  C->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX = "N_{T}";
  latexTitleX->DrawLatex(0.90, 0.45, TitleX);

   C->cd(0);
  padTitleY1->Draw();
  padTitleY1->cd();
  const char* TitleY1 = "Y(N_{T}) / ";
  latexTitleY1->SetTextAngle(90);
  latexTitleY1->DrawLatex(0.6, 0.5, TitleY1);

  C->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY2 = "#LT Y(N_{T}) #GT";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.6, 0.05, TitleY2);

  C->SaveAs(Form("./PtvNt_Yield_%s.pdf",part.Data()));
  C->SaveAs(Form("./PtvNt_Yield_%s.png",part.Data()));


}

