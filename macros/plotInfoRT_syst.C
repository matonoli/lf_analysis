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



void plotInfoRT_syst(Int_t particle = 1, Bool_t isUnc = false) {

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
  TString unc("RatioToHM"); 

  enum sysSources {
    sysRadiusL, sysDCAdd, sysCPA, sysFastSignal,
    sysCompMass, sysLifetime, sysNSigmaTPC, sysDCAPVpos,
    sysDCAPVneg, sysNCluster, sysNClusterF, sysSizeof
  };
  const char* SYSTS[sysSizeof+2] = { 
    "sysRadiusL", "sysDCAdd", "sysCPA", "sysFastSignal",
    "sysCompMass", "sysLifetime", "sysNSigmaTPC", "sysDCAPVpos",
    "sysDCAPVneg", "sysNCluster", "sysNClusterF",
    "",""
  };
  const char* PLOTS_SYSTS[sysSizeof+2] = { 
    "radius", "DCA_{d-d}", "cos PA", "fast signals",
    "comp. mass #sigma's", "lifetime", "TPC #sigma's", "DCA_{PV-d+}",
    "DCA_{PV-d-}", "TPC crossed rows", "TPC find. ratio",
    "signal extr.", "total"
  };

  const Int_t NREGIONS = 3;//5;
  const char* REGIONS[] = {"Trans","Near","Away"};

  const int NRTBINS = 5;
  const double RTBINS[NRTBINS+1] = {5.0,0.0, 0.8, 1.5, 2.5, 5.0};


  // DEFINE WORKING HISTOGRAMS HERE
  TH1D* hSyst[NREGIONS][sysSizeof+2][NRTBINS];

  // FETCH SOURCE HISTOGRAMS FROM FILES
  //TFile* fIn = new TFile(Form("%s/outMC_230123.root",path.Data()),"READ");
  TFile* fIn = new TFile(Form("%s/outD_230301_0xcs_hist.root",path.Data()),"READ");

  //TString part("L" ); TString partL("#Lambda");
  //TString part("Lbar" ); TString partL("#bar{#Lambda}");

  TDirectoryFile* dirCorr = (TDirectoryFile*)fIn->Get("MyAnalysisV0syst_3");
  for (int iReg = 0; iReg < 3; ++iReg)  {
  for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)  {
    
    for (int iSo = 0; iSo < sysSizeof; ++iSo) {
      hSyst[iReg][iSo][iRtBin] = (TH1D*)dirCorr->Get(Form("hV0PtRtSysMaxD%s_%s_%s_%s_%i",isUnc?unc.Data():"",part.Data(),REGIONS[iReg],SYSTS[iSo],iRtBin));
    } 
    
    hSyst[iReg][sysSizeof][iRtBin] = (TH1D*)dirCorr->Get(Form("hV0PtRtSysMaxDSigEx%s_%s_%s_%i",isUnc?unc.Data():"",part.Data(),REGIONS[iReg],iRtBin));
    hSyst[iReg][sysSizeof+1][iRtBin] = (TH1D*)dirCorr->Get(Form("hV0PtRtSysSum%s_%s_%s_%i",isUnc?"Unc":"",part.Data(),REGIONS[iReg],iRtBin));

  } } 


  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1000,1100);
  C->SetFillStyle(4000);


  // Number of PADS
  const int Nx = 3;
  const int Ny = sysSizeof+3;

  // Margins
  float lMargin = 0.04;
  float rMargin = 0.04;
  float bMargin = 0.02;
  float tMargin = 0.022;

  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  int latexTextSize = 18; 
  float legendTextSize = 0.07;

  double minYratio = -0.01;
  double maxYratio = 0.09;
  if (part.Contains("L")) {
    cout << "Lambda range " << endl;
    minYratio = -0.02;
    maxYratio = 0.19;
  }

  const double minXleft = 0.35;
  const double maxXleft = 8.9;


  /// Left frames

  TH1D* hframeRatioLeft = new TH1D("hframeRatioLeft","",100 , minXleft , maxXleft);
  hframeRatioLeft->GetYaxis()->SetRangeUser(minYratio,maxYratio);
  hframeRatioLeft->GetXaxis()->SetMoreLogLabels(kTRUE);
  hframeRatioLeft->GetXaxis()->SetNoExponent(kTRUE);
  hframeRatioLeft->GetYaxis()->SetNdivisions(502, kTRUE);
  if (part.Contains("L")) hframeRatioLeft->GetYaxis()->SetNdivisions(504, kTRUE);
  hframeRatioLeft->GetXaxis()->SetTickSize(0.07);
  hframeRatioLeft->GetYaxis()->SetTickSize(0.03);
  hframeRatioLeft->SetLabelFont(43,"xyz");
  hframeRatioLeft->SetLabelSize(15,"xyz");
  //hframeRatioLeft->GetYaxis()->SetTitle("Relative deviation");
  

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

  TLegend* legendLeft = new TLegend(0.04,0.3,0.65,0.93);
  myLegendSetUp(legendLeft,legendTextSize);
  legendLeft->SetNColumns(3);
  legendLeft->SetTextFont(43);
  legendLeft->SetTextSize(20);
  legendLeft->AddEntry(hSyst[0][0][1],"0.0 < R_{T} < 0.8","PL");
  legendLeft->AddEntry(hSyst[0][0][2],"0.8 < R_{T} < 1.5","PL");
  if (!isUnc) legendLeft->AddEntry(hSyst[0][0][0],"0.0 < R_{T} < 5.0","PL");
  else legendLeft->AddEntry(hSyst[0][0][0],"","");
  legendLeft->AddEntry(hSyst[0][0][3],"1.5 < R_{T} < 2.5","PL");
  legendLeft->AddEntry(hSyst[0][0][4],"2.5 < R_{T} < 5.0","PL");

  TLegend* legendLeft2 = new TLegend(0.04,0.67,0.75,0.93);
  myLegendSetUp(legendLeft2,legendTextSize);

  TLegend* legendRight = new TLegend(0.8,0.7,0.95,0.85);
  myLegendSetUp(legendRight,legendTextSize);

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
      pad[i][j]->SetLogx(kTRUE);
      
    }
  }

  TPaletteAxis *palette;

  for (Int_t ir = 0; ir < Ny-1; ir++){

    //if (ir == Ny-1) hframeRatioLeft->GetYaxis()->SetTickLength(0.08);

    C->cd(0);
    pad[0][ir]->Draw();
    pad[0][ir]->cd();
    cout << pad[0][ir]->GetWw() << " X " << pad[0][ir]->GetWh() << endl;
    cout << pad[0][ir]->GetAbsWNDC() << " X " << pad[0][ir]->GetAbsHNDC() << endl;
    hframeRatioLeft->Draw("axis");

    MakeNiceHistogram(hSyst[0][ir][0], colors[0]);
    MakeNiceHistogram(hSyst[0][ir][1], colors[1]);
    MakeNiceHistogram(hSyst[0][ir][2], colors[2]);
    MakeNiceHistogram(hSyst[0][ir][3], colors[3]);
    MakeNiceHistogram(hSyst[0][ir][4], colors[4]);
    hSyst[0][ir][0]->SetMarkerStyle(24);
    
    hSyst[0][ir][1]->Draw("pl same");
    hSyst[0][ir][2]->Draw("pl same");
    hSyst[0][ir][3]->Draw("pl same");
    hSyst[0][ir][4]->Draw("pl same");
    if (!isUnc) hSyst[0][ir][0]->Draw("pl same");

    //latexR->DrawLatex(0.20,0.75,"Toward");
    //latexSp->DrawLatex(0.65,0.75,Form("#bf{%s}",partL.Data()));
    if (ir == sysSizeof+1) latexSystem->DrawLatex(0.2,isUnc?0.7:0.25,PLOTS_SYSTS[ir]);
    else if (ir == sysSizeof) latexSystem->DrawLatex(0.3,0.7,PLOTS_SYSTS[ir]);
    else if (ir == sysSizeof-1) latexSystem->DrawLatex(0.2,0.6,PLOTS_SYSTS[ir]);
    else latexSystem->DrawLatex(0.2,0.8-0.02*ir,PLOTS_SYSTS[ir]);

    C->cd(0);
    pad[1][ir]->Draw();
    pad[1][ir]->cd();
    hframeRatioLeft->Draw("axis");

    MakeNiceHistogram(hSyst[1][ir][0], colors[0]);
    MakeNiceHistogram(hSyst[1][ir][1], colors[1]);
    MakeNiceHistogram(hSyst[1][ir][2], colors[2]);
    MakeNiceHistogram(hSyst[1][ir][3], colors[3]);
    MakeNiceHistogram(hSyst[1][ir][4], colors[4]);
    hSyst[1][ir][0]->SetMarkerStyle(24);
    
    hSyst[1][ir][1]->Draw("pl same");
    hSyst[1][ir][2]->Draw("pl same");
    hSyst[1][ir][3]->Draw("pl same");
    hSyst[1][ir][4]->Draw("pl same");
    if (!isUnc) hSyst[1][ir][0]->Draw("pl same");
    

    C->cd(0);
    pad[2][ir]->Draw();
    pad[2][ir]->cd();
    hframeRatioLeft->Draw("axis");

    MakeNiceHistogram(hSyst[2][ir][0], colors[0]);
    MakeNiceHistogram(hSyst[2][ir][1], colors[1]);
    MakeNiceHistogram(hSyst[2][ir][2], colors[2]);
    MakeNiceHistogram(hSyst[2][ir][3], colors[3]);
    MakeNiceHistogram(hSyst[2][ir][4], colors[4]);
    hSyst[2][ir][0]->SetMarkerStyle(24);
    
    hSyst[2][ir][1]->Draw("pl same");
    hSyst[2][ir][2]->Draw("pl same");
    hSyst[2][ir][3]->Draw("pl same");
    hSyst[2][ir][4]->Draw("pl same");
    if (!isUnc) hSyst[2][ir][0]->Draw("pl same");

    
  }



  C->cd();
  pad[0][Ny-1]->cd();
  latexTitleX->DrawLatex(0.11, 0.25, "Transverse");
  latexTitleX->DrawLatex(0.85, 0.25, partL.Data());

  C->cd();
  pad[1][Ny-1]->cd();
  latexTitleX->DrawLatex(0.05, 0.25, "Toward");

  C->cd();
  pad[2][Ny-1]->cd();
  latexTitleX->DrawLatex(0.04, 0.25, "Away");


  C->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX = "#it{p}_{T} (GeV/#it{c})";
  latexTitleX->DrawLatex(0.80, 0.55, TitleX);
  legendLeft->Draw();

  if (!isUnc) {
    C->cd(0);
    padTitleY1->Draw();
    padTitleY1->cd();
    const char* TitleY1 = "Relative";
    latexTitleY1->SetTextAngle(90);
    latexTitleY1->DrawLatex(0.6, 0.78, TitleY1);

    C->cd(0);
    padTitleY2->Draw();
    padTitleY2->cd();
    const char* TitleY2 = "uncertainty";
    latexTitleY2->SetTextAngle(90);
    latexTitleY2->DrawLatex(0.62, 0.0, TitleY2);
  } else {
    C->cd(0);
    padTitleY1->Draw();
    padTitleY1->cd();
    const char* TitleY1 = "Uncorrelated unc.";
    latexTitleY1->SetTextAngle(90);
    latexTitleY1->DrawLatex(0.6, 0.57, TitleY1);

    C->cd(0);
    padTitleY2->Draw();
    padTitleY2->cd();
    const char* TitleY2 = "on ratios to R_{T} < 5.0";
    latexTitleY2->SetTextAngle(90);
    latexTitleY2->DrawLatex(0.62, 0.0, TitleY2);
  }

  C->SaveAs(Form("./InfoRT_syst%s_%s.pdf",isUnc?"Unc":"",part.Data()));
  C->SaveAs(Form("./InfoRT_syst%s_%s.png",isUnc?"Unc":"",part.Data()));


  
}

