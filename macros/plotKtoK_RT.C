#include <iostream>
#include "MakeCanvas.C"

using namespace std;


//// UNCOMMENT A LINE TO CHOOSE THE MULTIPLICITY SELECTION


// COMMON PATH FOR INPUTS
TString path("../Unfolding");

// HELPER GLOBAL VARIABLES
enum {Transverse, Near, Away, NREGIONS};
enum {LowRT, HighRT, Int, NRT};

const char* strS[NREGIONS] = {"Transverse", "Near", "Away"};
const char* strH[NRT] = {"N_{T}<7", "N_{T}>14", "Integrated"};
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

void DrawHistograms(TH1D** h, Int_t nhist) {

  
  Int_t lineWidth = 2;

  for (int iH = 0; iH < nhist; ++iH) {
    
    MakeNiceHistogram(h[iH],colors[iH]);
  }

  h[Int]->Draw("E X0 same");
  h[LowRT]->Draw("E X0 same");
  h[HighRT]->Draw("E X0 same");

}
TH1D* MakeRatioHist(TH1D* hn, TH1D* hd) {

  TH1D* hr = (TH1D*)hn->Clone(Form("hr_%s",hn->GetName()));
  hr->Divide(hd);

  return hr;
}

void DrawHistogramsRatios(TH1D** hn, TH1D** hd, Int_t nhist, Int_t col) {

  TH1D* hr[NRT];
  for (int iH = 0; iH < nhist; ++iH) {
    
    hr[iH] = MakeRatioHist(hn[iH],hd[iH]);
  }

  hr[Int]->Draw("E X0 same");
  hr[LowRT]->Draw("E X0 same");
  hr[HighRT]->Draw("E X0 same");
}



void plotKtoK_RT() {

  TH1::SetDefaultSumw2(1);

  TH2F* hK0sNt[NREGIONS];
  TH2F* hKNt[NREGIONS];
  TH1D* hK0s[NREGIONS][NRT];
  TH1D* hK[NREGIONS][NRT];
  TH1D* hK_r[NREGIONS][NRT];
  TH1D* hK0soverK[NREGIONS][NRT];

  // kpm
  TFile* fIn = new TFile(Form("%s/outMC_221128.root",path.Data()),"READ");
  TDirectoryFile* dirCorr = (TDirectoryFile*)fIn->Get("MyAnalysisV0correct_2"); 
  hK0sNt[Transverse]  = (TH2F*)dirCorr->Get("hV0PtNtFitCorr_K0s_RC_Trans");
  hK0sNt[Near]        = (TH2F*)dirCorr->Get("hV0PtNtFitCorr_K0s_RC_Near");
  hK0sNt[Away]        = (TH2F*)dirCorr->Get("hV0PtNtFitCorr_K0s_RC_Away");
  
  hKNt[Transverse]  = (TH2F*)dirCorr->Get("hKpmPtNtRC_Trans");
  hKNt[Near]        = (TH2F*)dirCorr->Get("hKpmPtNtRC_Near");
  hKNt[Away]        = (TH2F*)dirCorr->Get("hKpmPtNtRC_Away");
  /*hK0sNt[Transverse]  = (TH2F*)dirCorr->Get("hV0PtNt_K0s_MC_Trans");
  hK0sNt[Near]        = (TH2F*)dirCorr->Get("hV0PtNt_K0s_MC_Near");
  hK0sNt[Away]        = (TH2F*)dirCorr->Get("hV0PtNt_K0s_MC_Away");
  
  hKNt[Transverse]  = (TH2F*)dirCorr->Get("hKpmPtNtMC_Trans");
  hKNt[Near]        = (TH2F*)dirCorr->Get("hKpmPtNtMC_Near");
  hKNt[Away]        = (TH2F*)dirCorr->Get("hKpmPtNtMC_Away");*/
  
  
  for (int iR = 0; iR < NREGIONS; ++iR)  {
    
    hK0s[iR][LowRT] = (TH1D*)hK0sNt[iR]->ProjectionX(Form("hK0s_%i_%i",iR,LowRT),1,6);
    hK0s[iR][HighRT] = (TH1D*)hK0sNt[iR]->ProjectionX(Form("hK0s_%i_%i",iR,HighRT),15,50);
    hK0s[iR][Int] = (TH1D*)hK0sNt[iR]->ProjectionX(Form("hK0s_%i_%i",iR,Int),0,50);

    hK[iR][LowRT] = (TH1D*)hKNt[iR]->ProjectionX(Form("hK_%i_%i",iR,LowRT),1,6);
    hK[iR][HighRT] = (TH1D*)hKNt[iR]->ProjectionX(Form("hK_%i_%i",iR,HighRT),15,50);
    hK[iR][Int] = (TH1D*)hKNt[iR]->ProjectionX(Form("hK_%i_%i",iR,Int),0,50);

  }



  // CALCULATE RATIOS
  for (int iR = 0; iR < NREGIONS; ++iR)  {
  for (int iN = 0; iN < NRT; ++iN)  {

    hK0soverK[iR][iN] = (TH1D*)hK0s[iR][iN]->Clone(Form("hK0soverK_%i_%i",iR,iN));
    hK0soverK[iR][iN]->Divide(hK0soverK[iR][iN],hK[iR][iN],2.);
    

  } }


  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1200,500);
  C->SetFillStyle(4000);

  // Number of PADS
  const int Nx = 3;
  const int Ny = 1;

  // Margins
  float lMargin = 0.03;
  float rMargin = 0.03;
  float bMargin = 0.05;
  float tMargin = 0.05;

  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  int latexTextSize = 18; 
  float legendTextSize = 0.05;

  const double minYratio = 0.7;
  const double maxYratio = 1.62;

  const double minXleft = 0.39;
  const double maxXleft = 5.21;


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
  TuneLatex(latexSp,1.6*latexTextSize);

  TLatex* latexSystem = new TLatex();
  TuneLatex(latexSystem,latexTextSize);

  TLegend* legendLeft = new TLegend(0.45,0.71,0.91,0.90);
  myLegendSetUp(legendLeft,legendTextSize);
  legendLeft->AddEntry(hK0soverK[Transverse][2],"N_{T}-integrated","P");
  legendLeft->AddEntry(hK0soverK[Transverse][0],"N_{T} < 7","P");
  legendLeft->AddEntry(hK0soverK[Transverse][1],"N_{T} > 14","P");

  TLegend* legendRight = new TLegend(0.8,0.7,0.95,0.85);
  myLegendSetUp(legendRight,legendTextSize);
  //legendRight->AddEntry(hK0soverK[Ref][Monash],"Monash","L");
  //legendRight->AddEntry(hK0soverK[Ref][Ropes],"Ropes","L");

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
      //pad[i][j]->SetLogx(kTRUE);
      
    }
  }

  C->cd(0);
  pad[0][0]->Draw();
  pad[0][0]->cd();
  hframeRatioLeft->Draw();
  lineleft->Draw("same");
  DrawHistograms(hK0soverK[Transverse],NRT);

  latexSystem->DrawLatex(0.15,0.92,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.15,0.86,"|#eta|<0.8");
  latexSystem->DrawLatex(0.15,0.80,"ALICE");

  latexSp->DrawLatex(0.15,0.7,"Transverse");

  C->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  hframeRatioLeft->Draw();
  lineleft->Draw("same");
  DrawHistograms(hK0soverK[Near],NRT);
  latexSp->DrawLatex(0.15,0.7,"Near");

  C->cd(0);
  pad[2][0]->Draw();
  pad[2][0]->cd();
  hframeRatioLeft->Draw();
  lineleft->Draw("same");
  DrawHistograms(hK0soverK[Away],NRT);
  latexSp->DrawLatex(0.15,0.7,"Away");
  legendLeft->Draw();

  C->cd(0);

  
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
  latexTitleY1->DrawLatex(0.5, 0.65, TitleY1);
  C->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY2 = " K^{#pm}";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.5, 0.0, TitleY2);

  C->SaveAs(Form("./KtoK.pdf"));


}

