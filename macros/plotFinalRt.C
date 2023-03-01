#include <iostream>

using namespace std;

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

  // for ratios of S0 to HM, use the unc. uncertainties as systematics
  TString hrname(hr->GetName());
  if (hrname.Contains("Unc")) {
  for (int iBin = 1; iBin < hn->GetNbinsX()+1; iBin++) {
    double nerr = (double)hn->GetBinError(iBin)/hn->GetBinContent(iBin);
    hr->SetBinError(iBin,(double)hr->GetBinContent(iBin)*nerr);
  } }


  hr->GetYaxis()->SetTitle("ratio to R_{T} #geq 0");
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
    currentLegend->SetTextFont(42);
    currentLegend->SetBorderSize(0);
    currentLegend->SetFillStyle(0);
    currentLegend->SetFillColor(0);
    currentLegend->SetMargin(0.25);
    currentLegend->SetTextSize(currentTextSize);
//    currentLegend->SetEntrySeparation(0.35);
    return;
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
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.0);
  h->SetMarkerColor(col);
  h->SetStats(0);
  h->SetFillStyle(0);

  h->GetYaxis()->SetTitleSize(30);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetLabelOffset(0.0025);
  h->GetYaxis()->SetLabelSize(20);
  h->GetYaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelFont(43);

  //h->SetTopMargin(0.055);
}



Bool_t IS_EXCLUSIVE = 0;


const Int_t NPTBINS = 16;
  const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
    0.4, 0.6, 0.8, 1.0, 1.2, 
    1.4, 1.6, 1.8, 2.0, 2.2, 
    2.6, 3.0, 3.4, 4.0, 5.0, 
    6.5, 8.0 };


/*const Int_t NPTBINS = 16;
  const Double_t XBINS[NPTBINS+1] = { // divisible by omar's pions
    0.4, 0.6, 0.8, 1.0, 1.2, 
    1.4, 1.6, 1.8, 2.0, 2.2, 
    2.6, 3.0, 3.4, 4.0, 5.0, 
    6.5, 8.0 };*/

const Int_t NPIONBINS = 51;   //Omar pi+- spectra
  const Double_t PIONBINS[NPIONBINS+1] = {
    0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 
    0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 
    0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 
    1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 
    2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.50, 5.00, 
    5.50, 6.00, 6.50, 7.00, 8.00, 10.0, 20.0 };

Int_t colors[5] = {kBlack, kRed, kBlack, kGreen+2, kBlue};

TString reg("TransMin");


void plotFinalRt() {

  TH1::SetDefaultSumw2(1);

	TFile* fin = new TFile("output220503.root", "READ");
  TDirectoryFile* dirCorr = (TDirectoryFile*)fin->Get("MyAnalysisV0correct_2"); 
  
  fin = new TFile("../rootOutputs/beast_220503_201inc_hist.root","READ");
  TDirectoryFile* dirMain = (TDirectoryFile*)fin->Get("MyAnalysisV0_0");

  //hKMCs[0]->Divide(hKMCs[3]);hKMCs[1]->Divide(hKMCs[4]);hKMCs[2]->Divide(hKMCs[5]);

  TFile* fout = new TFile(Form("k0s_l_spectra_rt_%s.root",reg.Data()),"RECREATE");

  // Load histograms
  TH1D* hK0s_Ref_Stat = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorr_K0s_D_%s_0",reg.Data())); 
  TH1D* hK0s_0_Stat = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorr_K0s_D_%s_2",reg.Data())); 
  TH1D* hK0s_1_Stat = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorr_K0s_D_%s_3",reg.Data())); 
  TH1D* hK0s_2_Stat = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorr_K0s_D_%s_4",reg.Data())); 
  TH1D* hK0s_3_Stat = (TH1D*)dirCorr->Get(Form("hV0PtRtFitCorr_K0s_D_%s_5",reg.Data())); 
  


  TH1D* hStats[5] = {
    hK0s_Ref_Stat, hK0s_0_Stat, hK0s_1_Stat,
    hK0s_2_Stat, hK0s_3_Stat };

  const char* statNames[5] = {
    "hK0s_Ref_Stat", "hK0s_0_Stat", "hK0s_1_Stat",
    "hK0s_2_Stat", "hK0s_3_Stat"};

  
  for (int iH = 0; iH < 5; iH++) {

    hStats[iH]->SetName(statNames[iH]);
    //hStats[iH]->Write(); hSysts[iH]->Write();
  }    


  // Add L and Lbar together
/*  TH1D* hLLbar_Ref_Stat = (TH1D*)hL_Ref_Stat->Clone("hLLbar_Ref_Stat");
  TH1D* hLLbar_Jet_Stat = (TH1D*)hL_Jet_Stat->Clone("hLLbar_Jet_Stat");
  TH1D* hLLbar_Iso_Stat = (TH1D*)hL_Iso_Stat->Clone("hLLbar_Iso_Stat");
  TH1D* hLLbar_Ref_Syst = (TH1D*)hL_Ref_Syst->Clone("hLLbar_Ref_Syst");
  TH1D* hLLbar_Jet_Syst = (TH1D*)hL_Jet_Syst->Clone("hLLbar_Jet_Syst");
  TH1D* hLLbar_Iso_Syst = (TH1D*)hL_Iso_Syst->Clone("hLLbar_Iso_Syst");
  TH1D* hLLbar_Ref_SystUnc = (TH1D*)hL_Ref_SystUnc->Clone("hLLbar_Ref_SystUnc");
  TH1D* hLLbar_Jet_SystUnc = (TH1D*)hL_Jet_SystUnc->Clone("hLLbar_Jet_SystUnc");
  TH1D* hLLbar_Iso_SystUnc = (TH1D*)hL_Iso_SystUnc->Clone("hLLbar_Iso_SystUnc");
  hLLbar_Ref_Stat->Add(hLbar_Ref_Stat);
  hLLbar_Jet_Stat->Add(hLbar_Jet_Stat);
  hLLbar_Iso_Stat->Add(hLbar_Iso_Stat);
  hLLbar_Ref_Syst->Add(hLbar_Ref_Syst);
  hLLbar_Jet_Syst->Add(hLbar_Jet_Syst);
  hLLbar_Iso_Syst->Add(hLbar_Iso_Syst);
  hLLbar_Ref_SystUnc->Add(hLbar_Ref_SystUnc);
  hLLbar_Jet_SystUnc->Add(hLbar_Jet_SystUnc);
  hLLbar_Iso_SystUnc->Add(hLbar_Iso_SystUnc);

  TH1D* hLLbar_Stat[3] = {
    hLLbar_Ref_Stat, hLLbar_Jet_Stat, hLLbar_Iso_Stat };
  TH1D* hLLbar_Syst[3] = {
    hLLbar_Ref_Syst, hLLbar_Jet_Syst, hLLbar_Iso_Syst };
  TH1D* hLLbar_SystUnc[3] = {
    hLLbar_Ref_SystUnc, hLLbar_Jet_SystUnc, hLLbar_Iso_SystUnc };


  // Make L/K0s ratios
  TH1D* hLLbarToK0s_Ref_Stat = (TH1D*)hLLbar_Ref_Stat->Clone("hLLbarToK0s_Ref_Stat");
  TH1D* hLLbarToK0s_Jet_Stat = (TH1D*)hLLbar_Jet_Stat->Clone("hLLbarToK0s_Jet_Stat");
  TH1D* hLLbarToK0s_Iso_Stat = (TH1D*)hLLbar_Iso_Stat->Clone("hLLbarToK0s_Iso_Stat");
  TH1D* hLLbarToK0s_Ref_Syst = (TH1D*)hLLbar_Ref_Syst->Clone("hLLbarToK0s_Ref_Syst");
  TH1D* hLLbarToK0s_Jet_Syst = (TH1D*)hLLbar_Jet_Syst->Clone("hLLbarToK0s_Jet_Syst");
  TH1D* hLLbarToK0s_Iso_Syst = (TH1D*)hLLbar_Iso_Syst->Clone("hLLbarToK0s_Iso_Syst");
  TH1D* hLLbarToK0s_Ref_SystUnc = (TH1D*)hLLbar_Ref_SystUnc->Clone("hLLbarToK0s_Ref_SystUnc");
  TH1D* hLLbarToK0s_Jet_SystUnc = (TH1D*)hLLbar_Jet_SystUnc->Clone("hLLbarToK0s_Jet_SystUnc");
  TH1D* hLLbarToK0s_Iso_SystUnc = (TH1D*)hLLbar_Iso_SystUnc->Clone("hLLbarToK0s_Iso_SystUnc");

  TH1D* hKL_Stat[3] = {
    hLLbarToK0s_Ref_Stat, hLLbarToK0s_Jet_Stat, hLLbarToK0s_Iso_Stat };
  TH1D* hKL_Syst[3] = {
    hLLbarToK0s_Ref_Syst, hLLbarToK0s_Jet_Syst, hLLbarToK0s_Iso_Syst };
  TH1D* hKL_SystUnc[3] = {
    hLLbarToK0s_Ref_SystUnc, hLLbarToK0s_Jet_SystUnc, hLLbarToK0s_Iso_SystUnc };

  for (int iH = 0; iH < 3; ++iH) {
    hKL_Stat[iH]->Divide(hKL_Stat[iH],hStats[iH],.5);
    hKL_Syst[iH]->Divide(hKL_Syst[iH],hSysts[iH],.5);
    hKL_SystUnc[iH]->Divide(hKL_SystUnc[iH],hSystsUnc[iH],.5);
  }
*/


  // Produce plots -- spectra
  double rl = 0.; double rh = 9.;
  TH1D* hBlank = new TH1D("hBlank",";p_{T} (GeV/#it{c});K^{0}_{S}",10,rl,rh);
  MakeNiceHistogram(hBlank,kBlack);
  hK0s_0_Stat->SetTitle(";p_{T} (GeV/#it{c});");
  for (int i=1; i<11; i++) hBlank->SetBinContent(i,1000);

  for (int iH = 0; iH < 5; iH++) {
    MakeNiceHistogram(hStats[iH],colors[iH%5]);
    hStats[iH]->GetYaxis()->SetRangeUser(49e-5,100.*hStats[iH]->GetMaximum());
    hStats[iH]->Write();
      }

  TCanvas* cSpK0s = new TCanvas("cSpK0s","",800,1000);
  cSpK0s->SetLogy();

  hBlank->GetYaxis()->SetRangeUser(8e-6,5.01); hBlank->SetLineWidth(0);
  hBlank->Draw();
  //hK0s_Ref_Stat->Draw("E X0 same");
  hK0s_0_Stat->Draw("E X0 same");
  hK0s_1_Stat->Draw("E X0 same");
  hK0s_2_Stat->Draw("E X0 same");
  hK0s_3_Stat->Draw("E X0 same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,0.8,Form("%s",reg.Data()));
  latex.SetTextSize(0.030); latex.DrawLatex(0.5,0.76,"pp #sqrt{s} = 13 TeV, |y|<0.8");

  latex.SetTextColor(colors[1]); latex.DrawLatex(0.5,0.67,"0.0 #geq R_{T} < 0.5");
  latex.SetTextColor(colors[2]); latex.DrawLatex(0.5,0.64,"0.5 #geq R_{T} < 1.5");
  latex.SetTextColor(colors[3]); latex.DrawLatex(0.5,0.61,"1.5 #geq R_{T} < 2.5");
  latex.SetTextColor(colors[4]); latex.DrawLatex(0.5,0.58,"2.5 #geq R_{T} < 5.0");
  MakeRatioPlot(hBlank,hK0s_Ref_Stat,cSpK0s,0.1,5.6,rl,rh,"");
  MakeRatioPlot(hK0s_0_Stat,hK0s_Ref_Stat,cSpK0s,0.1,5.6,rl,rh,"e x0");
  MakeRatioPlot(hK0s_1_Stat,hK0s_Ref_Stat,cSpK0s,0.1,5.6,rl,rh,"e x0");
  MakeRatioPlot(hK0s_2_Stat,hK0s_Ref_Stat,cSpK0s,0.1,5.6,rl,rh,"e x0");
  MakeRatioPlot(hK0s_3_Stat,hK0s_Ref_Stat,cSpK0s,0.1,5.6,rl,rh,"e x0");
  

  cSpK0s->SaveAs(Form("sp_K0s_%s.pdf",reg.Data()));
  cSpK0s->SaveAs(Form("sp_K0s_%s.png",reg.Data()));
  

//-------------------------
  TIter objIt(fout->GetList(),kIterForward);
  TObject* obj = 0;
  while ( (obj = objIt()) ) {
    TString objName(obj->GetName());
    if (objName.BeginsWith("hr_hBlank")) {      fout->Remove(obj);  }
      
  }


}

