#include <iostream>

using namespace std;




TString mult("V0M01"); int multInt = 0; int multOmar = 0; TString multPlot("V0M 0-1%");
//TString mult("NCharged01"); int multInt = 1; int multOmar = 0; TString multPlot("CL1 0-1%");
//TString mult("V0M"); int multInt = 2; int multOmar = 1; TString multPlot("V0M 0-10%");
//TString mult("NCharged"); int multInt = 3; int multOmar = 1; TString multPlot("CL1 0-10%");

Double_t calcError(Double_t valA, Double_t valB, Double_t errA, Double_t errB) {

  Double_t res;
  res = errA*errA/(valA*valA) + errA*errB/(valB*valB);
  res = TMath::Sqrt(res);
  res = res*TMath::Abs(valA/valB);

  return res;
}

void makeLvS0() {

  TH1::SetDefaultSumw2(1);

  TH1D* hLLbarvS0_Stat = new TH1D("hLLbarvS0_Stat",";S0;Entries",9,0,9);
  TH1D* hLLbarvS0_Syst = new TH1D("hLLbarvS0_Syst",";S0;Entries",9,0,9);

  TH1D* hL_Stat;
  TH1D* hL_Syst;

  // integration range
  double low = 4;
  double high = 99.;

  // S0 unbiased
  double mb_Stat = 0;
  double mberr_Stat = 0;
  double mberr_Syst = 0;
  TFile* fin = new TFile(Form("k0s_l_spectra_%s_spher10.root",mult.Data()),"READ");
  hL_Stat = (TH1D*)fin->Get("hK0s_Ref_Stat");
  hL_Syst = (TH1D*)fin->Get("hK0s_Ref_SystUnc");

  mb_Stat = hL_Stat->IntegralAndError(low,high,mberr_Stat,"");
  for (int iB=1; iB<hL_Syst->GetNbinsX()+1; iB++) hL_Syst->SetBinContent(iB,hL_Syst->GetBinContent(iB)+hL_Syst->GetBinError(iB));
  mberr_Syst = hL_Syst->Integral(low,high) - mb_Stat;


  double s_Stat = 0;
  double serr_Stat = 0;
  double serr_Syst = 0;


  const char* fileN[9] = {Form("k0s_l_spectra_%s_spher1_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher5_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher10_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher20_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher20_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher20_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher10_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher5_excl.root",mult.Data()),
    Form("k0s_l_spectra_%s_spher1_excl.root",mult.Data()),
  };

  const char* sphN[9] = { "Jet","Jet","Jet","Jet","Ref","Iso","Iso","Iso","Iso" };
  TFile* fin;
  // jetty 1
  for (int iS = 1; iS < 10; iS++) {

    fin = new TFile(fileN[iS-1],"READ");
    hL_Stat = (TH1D*)fin->Get(Form("hK0s_%s_Stat",sphN[iS-1]));
    hL_Syst = (TH1D*)fin->Get(Form("hK0s_%s_SystUnc",sphN[iS-1]));
    s_Stat = hL_Stat->IntegralAndError(low,high,serr_Stat,"");
    for (int iB=1; iB<hL_Syst->GetNbinsX()+1; iB++) hL_Syst->SetBinContent(iB,hL_Syst->GetBinContent(iB)+hL_Syst->GetBinError(iB));
    serr_Syst = hL_Syst->Integral(low,high) - s_Stat;
    
    hLLbarvS0_Stat->SetBinContent(iS,s_Stat/mb_Stat);
    hLLbarvS0_Syst->SetBinContent(iS,s_Stat/mb_Stat);
    hLLbarvS0_Stat->SetBinError(iS,calcError(s_Stat,mb_Stat,serr_Stat,mberr_Stat));
    hLLbarvS0_Syst->SetBinError(iS,calcError(s_Stat,mb_Stat,serr_Syst,mberr_Syst));

    cout << "S0 bin " << iS << ": #Lambda = " << s_Stat/mb_Stat << " +- " << calcError(s_Stat,mb_Stat,serr_Stat,mberr_Stat) << endl;
  
  }

  TFile* fout = new TFile(Form("l_vs_s0_%s.root",mult.Data()),"RECREATE");
  hLLbarvS0_Stat->Write();
  hLLbarvS0_Syst->Write();

}