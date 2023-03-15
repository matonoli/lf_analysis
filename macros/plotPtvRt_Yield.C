#include <iostream>
#include "MakeCanvas.C"

using namespace std;


//// UNCOMMENT A LINE TO CHOOSE THE MULTIPLICITY SELECTION


// COMMON PATH FOR INPUTS
TString path(".");

// HELPER GLOBAL VARIABLES
Int_t colors[7] = {TColor::GetColor("#12B69B"),TColor::GetColor("#B53189"),kBlack, kRed, kBlack, kGreen+2, kBlue};
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

Double_t calcError(Double_t valA, Double_t valB, Double_t errA, Double_t errB) {

  Double_t res;
  res = errA*errA/(valA*valA) + errA*errB/(valB*valB);
  res = TMath::Sqrt(res);
  res = res*TMath::Abs(valA/valB);

  return res;
}

void CalculateSelfNormYield(Int_t i, TH1D* hy, TH1D* hn, TH1D* hd) {

  // integration range
  double low = 1.;
  double high = 99.;

  double d, dstat;
  double n, nstat;

  d = hd->IntegralAndError(low, high, dstat, "width");
  n = hn->IntegralAndError(low, high, nstat, "width");

  hy->SetBinContent(i,n/d);
  hy->SetBinError(i,calcError(n,d,nstat,dstat));

  cout << hy->GetName() << " in bin " << i << " gets " << n << " / " << d << " = " << n/d << endl;

}

void CalculateSelfNormYieldSyst(Int_t i, TH1D* hy, TH1D* hn, TH1D* hd) {

  // integration range
  double low = 1.;
  double high = 99.;

  double d, dstat, dsyst;
  double n, nstat, nsyst;

  d = hd->IntegralAndError(low, high, dstat, "width");
  n = hn->IntegralAndError(low, high, nstat, "width");

  TH1D* hd2 = (TH1D*)hd->Clone("");
  TH1D* hn2 = (TH1D*)hn->Clone("");

  for (int ib = 1; ib < hd->GetNbinsX()+1; ib++) {
    hd2->SetBinContent(ib, hd->GetBinContent(ib)+hd->GetBinError(ib));
  }

  for (int ib = 1; ib < hn->GetNbinsX()+1; ib++) {
    hn2->SetBinContent(ib, hn->GetBinContent(ib)+hn->GetBinError(ib));
  }
  nsyst = (hn2->Integral(low,high,"width") - n)/n;

  hy->SetBinContent(i,n/d);
  hy->SetBinError(i,nsyst*n/d);

  cout << hy->GetName() << " in bin " << i << " gets " << n << " / " << d << " = " << n/d << endl;
  delete hd2;
  delete hn2;
}


void plotPtvRt_Yield() {

  //TH1::SetDefaultSumw2(1);

  

  // DEFINE WORKING HISTOGRAMS HERE

  TH1D* hPtT[65][NHIST];
  TH1D* hPtTMin[65][NHIST];
  TH1D* hPtTMax[65][NHIST];
  TH1D* hPtN[65][NHIST];
  TH1D* hPtA[65][NHIST];

  TH1D* hLPtT[65][NHIST];
  TH1D* hLPtTMin[65][NHIST];
  TH1D* hLPtTMax[65][NHIST];
  TH1D* hLPtN[65][NHIST];
  TH1D* hLPtA[65][NHIST];

  TH1D* hYieldRT[NHIST];
  TH1D* hYieldRTMin[NHIST];
  TH1D* hYieldRTMax[NHIST];
  TH1D* hYieldRTN[NHIST];
  TH1D* hYieldRTA[NHIST];
  TH1D* hLYieldRT[NHIST];
  TH1D* hLYieldRTMin[NHIST];
  TH1D* hLYieldRTMax[NHIST];
  TH1D* hLYieldRTN[NHIST];
  TH1D* hLYieldRTA[NHIST];

  

  // FETCH SOURCE HISTOGRAMS FROM FILES
  TFile* fIn = new TFile(Form("%s/outD_230301_0xcsu_hist.root",path.Data()),"READ");
  Bool_t isRC = true; Bool_t isRTMin = false; TString part("K0s" ); TString partL("K^{0}_{S}"); TString dmc("D");
  
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

  part = TString("L");
  for (int i = 0; i<NRT+1; i++) {
    hLPtT[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Trans_%i","L",dmc.Data(),i));
    hLPtTMin[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnf_%s_%s_TransMin_%i","L",dmc.Data(),i));
    hLPtTMax[i][Stat] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnf_%s_%s_TransMax_%i","L",dmc.Data(),i));
    hLPtN[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Near_%i","L",dmc.Data(),i));
    hLPtA[i][Stat]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnf_%s_%s_Away_%i","L",dmc.Data(),i));

    hLPtT[i][Syst]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSyst_%s_%s_Trans_%i",part.Data(),dmc.Data(),i));
    hLPtTMin[i][Syst] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnfSyst_%s_%s_TransMin_%i",part.Data(),dmc.Data(),i));
    hLPtTMax[i][Syst] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnfSyst_%s_%s_TransMax_%i",part.Data(),dmc.Data(),i));
    hLPtN[i][Syst]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSyst_%s_%s_Near_%i",part.Data(),dmc.Data(),i));
    hLPtA[i][Syst]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSyst_%s_%s_Away_%i",part.Data(),dmc.Data(),i));

    hLPtT[i][SystUnc]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSystUnc_%s_%s_Trans_%i",part.Data(),dmc.Data(),i));
    hLPtTMin[i][SystUnc] = (TH1D*)dirCorr->Get(Form("hV0PtRtMinFitCorrUnfSystUnc_%s_%s_TransMin_%i",part.Data(),dmc.Data(),i));
    hLPtTMax[i][SystUnc] = (TH1D*)dirCorr->Get(Form("hV0PtRtMaxFitCorrUnfSystUnc_%s_%s_TransMax_%i",part.Data(),dmc.Data(),i));
    hLPtN[i][SystUnc]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSystUnc_%s_%s_Near_%i",part.Data(),dmc.Data(),i));
    hLPtA[i][SystUnc]    = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorrUnfSystUnc_%s_%s_Away_%i",part.Data(),dmc.Data(),i));

  }
  
  // PROCESS HISTOGRAMS HERE
  const int NRTBINS = 4;
  //const double RTBINS[NRTBINS+1] = {0.0, 0.5, 1.5, 2.5, 5.0};
  const double RTBINS[NRTBINS+1] = {0.0, 0.8, 1.5, 2.5, 5.0};
  hYieldRT[Stat] = new TH1D("hYieldRT[Stat]","",NRTBINS,RTBINS);
  hYieldRT[SystUnc] = new TH1D("hYieldRT[SystUnc]","",NRTBINS,RTBINS);
  hYieldRTMin[Stat] = new TH1D("hYieldRTMin[Stat]","",NRTBINS,RTBINS);
  hYieldRTMin[SystUnc] = new TH1D("hYieldRTMin[SystUnc]","",NRTBINS,RTBINS);
  hYieldRTMax[Stat] = new TH1D("hYieldRTMax[Stat]","",NRTBINS,RTBINS);
  hYieldRTMax[SystUnc] = new TH1D("hYieldRTMax[SystUnc]","",NRTBINS,RTBINS);
  hYieldRTN[Stat] = new TH1D("hYieldRTN[Stat]","",NRTBINS,RTBINS);
  hYieldRTN[SystUnc] = new TH1D("hYieldRTN[SystUnc]","",NRTBINS,RTBINS);
  hYieldRTA[Stat] = new TH1D("hYieldRTA[Stat]","",NRTBINS,RTBINS);
  hYieldRTA[SystUnc] = new TH1D("hYieldRTA[SystUnc]","",NRTBINS,RTBINS);

  hLYieldRT[Stat] = new TH1D("hLYieldRT[Stat]","",NRTBINS,RTBINS);
  hLYieldRT[SystUnc] = new TH1D("hLYieldRT[SystUnc]","",NRTBINS,RTBINS);
  hLYieldRTMin[Stat] = new TH1D("hLYieldRTMin[Stat]","",NRTBINS,RTBINS);
  hLYieldRTMin[SystUnc] = new TH1D("hLYieldRTMin[SystUnc]","",NRTBINS,RTBINS);
  hLYieldRTMax[Stat] = new TH1D("hLYieldRTMax[Stat]","",NRTBINS,RTBINS);
  hLYieldRTMax[SystUnc] = new TH1D("hLYieldRTMax[SystUnc]","",NRTBINS,RTBINS);
  hLYieldRTN[Stat] = new TH1D("hLYieldRTN[Stat]","",NRTBINS,RTBINS);
  hLYieldRTN[SystUnc] = new TH1D("hLYieldRTN[SystUnc]","",NRTBINS,RTBINS);
  hLYieldRTA[Stat] = new TH1D("hLYieldRTA[Stat]","",NRTBINS,RTBINS);
  hLYieldRTA[SystUnc] = new TH1D("hLYieldRTA[SystUnc]","",NRTBINS,RTBINS);
    

  for (int iRtBin = 1; iRtBin < NRT+1; iRtBin++) {

    cout << "aa " << hYieldRT[Stat] << " " << hPtT[iRtBin][Stat] << hPtT[0][Stat] << endl;
    cout << "aa " << hYieldRT[SystUnc] << " " << hPtT[iRtBin][SystUnc] << hPtT[0][SystUnc] << endl;
    CalculateSelfNormYield(iRtBin,hYieldRT[Stat],hPtT[iRtBin][Stat],hPtT[0][Stat]);
    CalculateSelfNormYield(iRtBin,hYieldRTMin[Stat],hPtTMin[iRtBin][Stat],hPtTMin[0][Stat]);
    CalculateSelfNormYield(iRtBin,hYieldRTMax[Stat],hPtTMax[iRtBin][Stat],hPtTMax[0][Stat]);
    CalculateSelfNormYield(iRtBin,hYieldRTN[Stat],hPtN[iRtBin][Stat],hPtN[0][Stat]);
    CalculateSelfNormYield(iRtBin,hYieldRTA[Stat],hPtA[iRtBin][Stat],hPtA[0][Stat]);

    CalculateSelfNormYieldSyst(iRtBin,hYieldRT[SystUnc],hPtT[iRtBin][SystUnc],hPtT[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hYieldRTMin[SystUnc],hPtTMin[iRtBin][SystUnc],hPtTMin[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hYieldRTMax[SystUnc],hPtTMax[iRtBin][SystUnc],hPtTMax[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hYieldRTN[SystUnc],hPtN[iRtBin][SystUnc],hPtN[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hYieldRTA[SystUnc],hPtA[iRtBin][SystUnc],hPtA[0][SystUnc]);

    CalculateSelfNormYield(iRtBin,hLYieldRT[Stat],hLPtT[iRtBin][Stat],hLPtT[0][Stat]);
    CalculateSelfNormYield(iRtBin,hLYieldRTMin[Stat],hLPtTMin[iRtBin][Stat],hLPtTMin[0][Stat]);
    CalculateSelfNormYield(iRtBin,hLYieldRTMax[Stat],hLPtTMax[iRtBin][Stat],hLPtTMax[0][Stat]);
    CalculateSelfNormYield(iRtBin,hLYieldRTN[Stat],hLPtN[iRtBin][Stat],hLPtN[0][Stat]);
    CalculateSelfNormYield(iRtBin,hLYieldRTA[Stat],hLPtA[iRtBin][Stat],hLPtA[0][Stat]);

    CalculateSelfNormYieldSyst(iRtBin,hLYieldRT[SystUnc],hLPtT[iRtBin][SystUnc],hLPtT[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hLYieldRTMin[SystUnc],hLPtTMin[iRtBin][SystUnc],hLPtTMin[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hLYieldRTMax[SystUnc],hLPtTMax[iRtBin][SystUnc],hLPtTMax[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hLYieldRTN[SystUnc],hLPtN[iRtBin][SystUnc],hLPtN[0][SystUnc]);
    CalculateSelfNormYieldSyst(iRtBin,hLYieldRTA[SystUnc],hLPtA[iRtBin][SystUnc],hLPtA[0][SystUnc]);

  }

  cout << hYieldRT[Stat] << endl;
  cout << hYieldRT[SystUnc] << endl;
  cout << hYieldRTMin[Stat] << endl;
  cout << hYieldRTMin[SystUnc] << endl;
  cout << hYieldRTMax[Stat] << endl;
  cout << hYieldRTMax[SystUnc] << endl;
  cout << hYieldRTN[Stat] << endl;
  cout << hYieldRTN[SystUnc] << endl;
  cout << hYieldRTA[Stat] << endl;
  cout << hYieldRTA[SystUnc] << endl;
  
  MakeNiceHistogram(hYieldRT[Stat],colors[0]);
  MakeNiceHistogram(hYieldRT[SystUnc],colors[0]);
  MakeNiceHistogram(hYieldRTMin[Stat],colors[0]);
  MakeNiceHistogram(hYieldRTMin[SystUnc],colors[0]);
  MakeNiceHistogram(hYieldRTMax[Stat],colors[0]);
  MakeNiceHistogram(hYieldRTMax[SystUnc],colors[0]);
  MakeNiceHistogram(hYieldRTN[Stat],colors[0]);
  MakeNiceHistogram(hYieldRTN[SystUnc],colors[0]);
  MakeNiceHistogram(hYieldRTA[Stat],colors[0]);
  MakeNiceHistogram(hYieldRTA[SystUnc],colors[0]);

  MakeNiceHistogram(hLYieldRT[Stat],colors[1]);
  MakeNiceHistogram(hLYieldRT[SystUnc],colors[1]);
  MakeNiceHistogram(hLYieldRTMin[Stat],colors[1]);
  MakeNiceHistogram(hLYieldRTMin[SystUnc],colors[1]);
  MakeNiceHistogram(hLYieldRTMax[Stat],colors[1]);
  MakeNiceHistogram(hLYieldRTMax[SystUnc],colors[1]);
  MakeNiceHistogram(hLYieldRTN[Stat],colors[1]);
  MakeNiceHistogram(hLYieldRTN[SystUnc],colors[1]);
  MakeNiceHistogram(hLYieldRTA[Stat],colors[1]);
  MakeNiceHistogram(hLYieldRTA[SystUnc],colors[1]);


  /*hLYieldRT[Stat]->GetYaxis()->SetRangeUser(0.2,3.1);
  hLYieldRT[Stat]->Draw();
  hLYieldRTMin[Stat]->SetMarkerStyle(21);
  hLYieldRTMax[Stat]->SetMarkerStyle(29);
  hLYieldRTMin[Stat]->Draw("same");
  hLYieldRTMax[Stat]->Draw("same");
  hLYieldRTN[Stat]->SetMarkerStyle(21);
  hLYieldRTA[Stat]->SetMarkerStyle(29);
  hLYieldRTN[Stat]->Draw("same");
  hLYieldRTA[Stat]->Draw("same");
  hYieldRTN[Stat]->Draw("same");
  hYieldRTA[Stat]->Draw("same");
  lineleft->Draw("same");*/
  //hYieldRTN[Stat]->Draw("same");
  //hYieldRTA[Stat]->Draw("same");


  // MAKE PLOTS
  SetStyle(kFALSE);

  TCanvas *C = (TCanvas*) gROOT->FindObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1000,500);
  C->SetFillStyle(4000);


  // Number of PADS
  const int Nx = 5;
  const int Ny = 1;

  // Margins
  float lMargin = 0.03;
  float rMargin = 0.05;
  float bMargin = 0.06;
  float tMargin = 0.05;

  // Canvas setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  int latexTextSize = 18; 
  float legendTextSize = 0.12;

  const double minYspectra = 0.41;
  const double maxYspectra = 2.79;

  const double minYratio = 0.41;
  const double maxYratio = 2.79;

  const double minXleft = -0.4;
  const double maxXleft = 5.4;


  /// Left frames
  TH1D* hframeSpectraLeft = new TH1D("hframeSpectraLeft","",100 , minXleft , maxXleft);
  hframeSpectraLeft->GetYaxis()->SetRangeUser(minYspectra,maxYspectra);
  hframeSpectraLeft->GetXaxis()->SetMoreLogLabels(kTRUE);
  hframeSpectraLeft->GetXaxis()->SetNoExponent(kTRUE);

  hframeSpectraLeft->GetXaxis()->SetNdivisions(510, kTRUE);
  hframeSpectraLeft->SetLabelFont(43,"xyz");
  hframeSpectraLeft->SetLabelSize(17,"xyz");

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
  TuneLatex(latexSp,0.8*latexTextSize);

  TLatex* latexR = new TLatex();
  TuneLatex(latexR,1.2*latexTextSize);

  TLatex* latexSystem = new TLatex();
  TuneLatex(latexSystem,0.9*latexTextSize);

  TLegend* legendLeft = new TLegend(0.075,0.68,0.94,0.81);
  myLegendSetUp(legendLeft,legendTextSize);
  legendLeft->AddEntry(hYieldRT[Stat]," K^{0}_{S}","P");
  legendLeft->AddEntry(hLYieldRT[Stat]," #Lambda+#bar{#Lambda}","P");


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
  pad[0][0]->Draw();
  pad[0][0]->cd();
  //pad[0][0]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  lineleft->Draw("same");
  
  hYieldRTN[Stat]->Draw("E X0 same");
  hLYieldRTN[Stat]->Draw("E X0 same");
  hYieldRTN[SystUnc]->Draw("E2 same");
  hLYieldRTN[SystUnc]->Draw("E2 same");
  
  latexSystem->DrawLatex(0.23,0.82,"pp, #sqrt{s} = 13 TeV");
  latexSystem->DrawLatex(0.23,0.77,"|#eta|<0.8");
  latexSystem->DrawLatex(0.23,0.72,"ALICE Data");

  latexR->DrawLatex(0.23,0.92,"Toward");

  
  C->cd(0);
  pad[1][0]->Draw();
  pad[1][0]->cd();
  //pad[0][0]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  lineleft->Draw("same");
  
  hYieldRTA[Stat]->Draw("E X0 same");
  hLYieldRTA[Stat]->Draw("E X0 same");
  hYieldRTA[SystUnc]->Draw("E2 same");
  hLYieldRTA[SystUnc]->Draw("E2 same");

  latexR->DrawLatex(0.18,0.92,"Away");
  latexSp->DrawLatex(0.12,0.84,"#it{N}_{h}: 0.4 < #it{p}_{T} < 8 GeV/#it{c}");
  legendLeft->Draw();


  C->cd(0);
  pad[2][0]->Draw();
  pad[2][0]->cd();
  //pad[0][0]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  lineleft->Draw("same");
  
  hYieldRT[Stat]->Draw("E X0 same");
  hLYieldRT[Stat]->Draw("E X0 same");
  hYieldRT[SystUnc]->Draw("E2 same");
  hLYieldRT[SystUnc]->Draw("E2 same");

  latexR->DrawLatex(0.16,0.92,"Transverse");


  C->cd(0);
  pad[3][0]->Draw();
  pad[3][0]->cd();
  //pad[0][0]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  lineleft->Draw("same");
  
  hYieldRTMin[Stat]->Draw("E X0 same");
  hLYieldRTMin[Stat]->Draw("E X0 same");
  hYieldRTMin[SystUnc]->Draw("E2 same");
  hLYieldRTMin[SystUnc]->Draw("E2 same");

  latexR->DrawLatex(0.13,0.92,"Trans., min");


  C->cd(0);
  pad[4][0]->Draw();
  pad[4][0]->cd();
  //pad[0][0]->SetLogy(kFALSE);
  hframeSpectraLeft->Draw();
  lineleft->Draw("same");
  
  hYieldRTMax[Stat]->Draw("E X0 same");
  hLYieldRTMax[Stat]->Draw("E X0 same");
  hYieldRTMax[SystUnc]->Draw("E2 same");
  hLYieldRTMax[SystUnc]->Draw("E2 same");

  latexR->DrawLatex(0.10,0.92,"Trans., max");

  

  C->cd(0);
  padTitleX->Draw();
  padTitleX->cd();
  const char* TitleX = "R_{T,(-,min,max)}";
  latexTitleX->DrawLatex(0.80, 0.55, TitleX);

   C->cd(0);
  padTitleY1->Draw();
  padTitleY1->cd();
  const char* TitleY1 = "#it{N}_{h} /";
  latexTitleY1->SetTextAngle(90);
  latexTitleY1->DrawLatex(0.6, 0.79, TitleY1);

  C->cd(0);
  padTitleY2->Draw();
  padTitleY2->cd();
  const char* TitleY2 = " #LT #it{N}_{h} #GT ";
  latexTitleY2->SetTextAngle(90);
  latexTitleY2->DrawLatex(0.6, 0.03, TitleY2);

  C->SaveAs(Form("./PtvRt_Yield_%s.pdf",part.Data()));
  C->SaveAs(Form("./PtvRt_Yield_%s.png",part.Data()));


}

