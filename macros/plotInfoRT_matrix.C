#include <iostream>
#include "MakeCanvas.C"

using namespace std;


//// UNCOMMENT A LINE TO CHOOSE THE MULTIPLICITY SELECTION


// COMMON PATH FOR INPUTS
//TString path("results_unfolding"); TString RMname("hNtRM_flip");
//TString path("results_unfolding_min"); TString RMname("hNtMinRM_flip");
TString path("results_unfolding_max"); TString RMname("hNtMaxRM_flip");

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
  else gStyle->SetPalette(kBird);
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
  h->SetMarkerSize(0.6);
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

void DrawHistogramsRatios(TH1D** hn, TH1D** hd, Int_t nhist, Int_t col) {

  TH1D* hr[nhist];
  for (int iH = 0; iH < nhist; ++iH) {
    
    hr[iH] = MakeRatioHist(hn[iH],hd[iH]);
    hr[iH]->Draw("E X0 same");
  }
}



void plotInfoRT_matrix() {

  TH1::SetDefaultSumw2(1);

  // DEFINE WORKING HISTOGRAMS HERE
  TH2D* hRM;
  TH2D* hUnf;
  

  // FETCH SOURCE HISTOGRAMS FROM FILES
  TFile* fIn = new TFile(Form("%s/1D_newClass_mc.root",path.Data()),"READ");
  cout << "File found " << fIn << endl;

  hRM  = (TH2D*)fIn->Get(RMname.Data());
  hUnf  = (TH2D*)fIn->Get("_hSolutionNT");
  cout << "Histograms found " << hRM << " " << hUnf << endl;

  // PROCESS HISTOGRAMS HERE
  Int_t nCols = hRM->GetNbinsX();
  Int_t nRows = hRM->GetNbinsY();
  for (int iR = 1; iR < nRows+1; ++iR)  {
    Double_t integral = hRM->Integral(1,nCols+1,iR,iR);
    for (int iC = 1; iC < nCols+1; ++iC)  {
      Double_t binContent = hRM->GetBinContent(iC,iR);
      if (binContent>0) hRM->SetBinContent(iC,iR,binContent/integral);
    }
  }

  TH2D* hFlip = (TH2D*)hUnf->Clone("hFlip");
  for (int i = 1; i < hUnf->GetNbinsX()+1; i++) {
  for (int j = 1; j < hUnf->GetNbinsY()+1; j++) {
      hFlip->SetBinContent(i,j,hUnf->GetBinContent(j,i));
      hFlip->SetBinError(i,j,hUnf->GetBinError(j,i));
  } }
  hFlip->GetXaxis()->SetTitle(hUnf->GetYaxis()->GetTitle());
  hFlip->GetYaxis()->SetTitle(hUnf->GetXaxis()->GetTitle());
  delete hUnf;
  hUnf = hFlip; hUnf->SetName("hUnf");

  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",850,500);
  C->SetFillStyle(4000);

  // Number of PADS
  const int Nx = 2;
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

  const double minYratio = 0;
  const double maxYratio = 44;

  const double minXleft = 0;
  const double maxXleft = 44;


  // Left frames
  TH2D* hframeRatioLeft = new TH2D("hframeRatioLeft","",44,minYratio,maxYratio,44, minXleft , maxXleft);
  hframeRatioLeft->GetZaxis()->SetRangeUser(0.,1.);
  
  TLine* lineleft = new TLine(0., 0., hRM->GetXaxis()->GetBinLowEdge(hRM->GetNbinsX()), hRM->GetXaxis()->GetBinLowEdge(hRM->GetNbinsX()));
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
  TuneLatex(latexSp,1.2*latexTextSize);

  TLatex* latexSystem = new TLatex();
  TuneLatex(latexSystem,latexTextSize);

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
  hframeRatioLeft->Draw("colz");
  hRM->Draw("colz same");
  pad[0][0]->Modified(); pad[0][0]->Update();
  palette = (TPaletteAxis*)hRM->GetListOfFunctions()->FindObject("palette");
  cout << palette << endl;
  palette->SetX1NDC(0.85);
  palette->SetX2NDC(0.91);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.6);
  pad[0][0]->Modified();
  pad[0][0]->Update();
  lineleft->Draw("same");
  
  latexSp->DrawLatex(0.32,0.14,"Response matrix");

  C->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  hframeRatioLeft->Draw("colz");
  hUnf->Draw("colz same");
  pad[0][0]->Modified(); pad[0][0]->Update();
  palette = (TPaletteAxis*)hUnf->GetListOfFunctions()->FindObject("palette");
  cout << palette << endl;
  palette->SetX1NDC(0.75);
  palette->SetX2NDC(0.81);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.6);
  pad[0][0]->Modified();
  pad[0][0]->Update();
  lineleft->Draw("same");

  latexSystem->DrawLatex(0.07,0.92,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.07,0.86,"|#eta|<0.8, p_{T}>0.15 GeV/#it{c}");
  latexSystem->DrawLatex(0.07,0.80,"ALICE");
  latexSp->DrawLatex(0.28,0.14,"Unfolding matrix");
  
  //legendLeft->Draw();
  //latexSp->DrawLatex(0.15,0.7,"Transverse");

  C->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX = "N^{TS}_{acc,rec}";
  latexTitleX->DrawLatex(0.45, 0.5, TitleX);

  C->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY2 = "N^{TS}_{ch,true}";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.5, 0.0, TitleY2);
  

  C->SaveAs(Form("./InfoRT_matrix.pdf"));


}

