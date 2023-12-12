#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TList.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraphErrors.h>

#include "TVector.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace std;

#ifndef UNFOLNTCLASS_H
#define UNFOLNTCLASS_H

class UnfoldNTclass {

	public:

		UnfoldNTclass();	
		~UnfoldNTclass();
		void Setup(const TH1F* measured, const TH1F* truth, const TH2F* response);  // set up from already-filled histograms
		void SetnIter(int iter) { _niter = iter; }
		void UnfoldV02D(TH2F* , TH2F*, TH2F*);
		void Unfold2D(TH2F* , TH2F*);
		//void Unfold2D_Alternative_Approach(TH2F* , TH2F* );
		void CorrectPTdistributions();
		void CorrectPTdistributions(const TH2F* );
		void GetUnfoldingMatrix2D();
		void GetUnfoldingMatrix3D( int );
		void GetCovarianceMatrix2D();
		void FillCovarianceMatrix( int );
		void Unfold();
		void UnfoldTowardAway();
		void UnfoldPID(const TH1F* measured, const TH1F* truth, const TH2F* response); 
		void H2V(const TH1F*, TVectorD& );
		void H2M(const TH2F*, TMatrixD& );
		void GetUnfoldingMatrix();
		void GetUnfoldedDist();
		void FillMatrixFromNTSolution();
		void GetNewPrior();
		double GetChi2();
		TVectorD GetChi2Vectors(bool);
		void ExtrapolateRM(TH2F* );
		void V2H();
		TH1F* GetUnfoldedDistH();
		void GetMCclosureinNchBins(const TH2F* );
		void GetMCclosureinRTBins(const TH2F*, const TH2F*);
		void GetMCclosureinRTMinBins(const TH2F*, const TH2F*);
		void GetMCclosureinRTMaxBins(const TH2F*, const TH2F*);
		TH1F* RebinNT2RT(const TH1F* , bool);
		TH1F* RebinPtAxis(const TH1F* );
		void DrawNchClosure(const TH1F* ,const TH1F* ,const double& ,const double& ,const char* ,const char* ,const char*);
		void DrawClosureRatios(TH1F* ,TH1F* , const double& , const double& , const int& , const char* ,const char* , const double& , const double& , const int& , const bool , const bool, const char* );
		TPad* GetPad();
		TH1F* TuneFrame(const double& , const double& , const char* , const char* );
		void myLegendSetUp(TLegend*, float );
		void TuneLatex(TLatex* );
		void myOptions(int lStat);
		TH2F* SumCharges2D(const char*,const TH2F*,const TH2F*);
		TH2F* SumDistPerEta2D(const TH2F*,const TH2F*,const TH2F*,const TH2F*);

		//! Cloosure test
		void PerformClosure();
		TH1F* GetClosureH();

		//! Propagation of errors
		void Emeasured();   // Measured distribution errors as a TVectorD
		void H2VE();
		void GetMeasuredCov();
		void GetCovMatrixUnfolded();
		void SetError(TH2F* );	
		void SetError(TH1F* );	

		void SetMCAnalysis(bool ismc) { isMC = ismc; }
		void SetPidIdx(int idx) { _PidIdx = idx; }
		void SetRegion(const char*);
		void SetDetector(const char*);
		void SetCharge(const char*);
		void SetPid(const char*);
		void SaveSolutionNT(bool);
///		void TunePrior();
////		void LoadEventRecEfficiency();
		void LoadSolutionNT();
		void LoadSolutionNTMin();
		void LoadSolutionNTMax();

		void SetIsSystematics(bool issys) { isSys = issys; };
		void SetDirectory(const char* dir) { directory = dir; };
		void SetTrkID(string trkid) { TrkID = trkid; };

		const char* Pid() { return _Pid; }
		const char* GetRegion() { return _Region; } 
		const char* PidLatex(int p) {return ArrayPidLatex[p]; }
		const char* PidLatexPos(int p) {return ArrayPidLatexPos[p]; }
		const char* PidLatexNeg(int p) {return ArrayPidLatexNeg[p]; }
		const char* Eta(int e) { return EtaArray[e]; }

		int GetnNchBins() { return nNchBins; }
		int GetnPtBins() { return nPtBins; }
		double* GetPtBins() { return ptBins; }
		double* GetNchBins() { return NchBins; }   
		int GetnRTBins() { return nRTBins; }
		double* GetRTBins() { return RTBins; }

		void FillHistoVsRT(TH2F*, TH2F*);
		int GetRTBin(const int&,bool);
		int GetRTMinBin(const int&,bool);
		int GetRTMaxBin(const int&,bool);


		TObjArray* GetObjArray() { return ObjArray; }

	private:

		void NormRowWise(TMatrixD& );
		void NormRowWise(TVectorD& );
		void CheckNorm(TMatrixD& );

		TMatrixD _Res; // Response matrix
		TMatrixD _Mij; // unfolding matrix
		TMatrixD _Mij0; // Previous unfolding matrix
		TMatrixD _dnCidnEj; // measurement error propagation matrix
		TVectorD _UnfDist; // Unfolded dist.
		TVectorD _UnfDist0; // Previous unfolded dist.
		TVectorD _P0C; // Vector with prior
		TVectorD _nbarCi; // New estimate of true distribution
		TVectorD _Meas; // Vector with measured dist.
		TVectorD _True; // Vector with measured dist.
		TVectorD _eTrue; // Vector with true dist.

		TVectorD _eMes; // Vector with measured  uncertainties
		TMatrixD _covMes; // Measurement covariance matrix
		TMatrixD _covUnf; // Unfolded covariance matrix
		const TH2F* _res; // Response matrix (not owned)
		const TH1F* _meas;  // Measured distribution (not owned)
		const TH1F* _tru;  // True distribution (not owned)
		int _niter;
		int _nm;
		int _nt;
		double _N0C; // number of events in prior
		///double _nbartrue; // best estimate of number of true events
		TVectorD _vChi2x;
		TVectorD _vChi2y; // Vector with Chi2


		bool isMC;
		bool isSys;
		string TrkID;
		int _PidIdx = 0;
		const char* _Region;
		const char* _Detector;
		const char* _Pid;
		const char* _Charge;
		const char* directory;
		bool _SaveSolutionNT;
		TObjArray* ObjArray;

		/*const int nPtBins = 45;
		double ptBins[45+1] = {
			0.0, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
			0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
			0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
			1.80, 1.90, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
			3.60, 3.80, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0};*/

		// PT BINS
	const Int_t NPIPTBINS = 16;
	const Double_t PIXBINS[16+1] = { // divisible by omar's pions
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 5.0, 
		6.5, 8.0 };
		
	const Int_t NXIPTBINS = 7; 
  	const Double_t XIXBINS[7+1] = {
  		0.6, 1.2, 1.6, 2.2, 
  		2.8, 3.6, 5.0, 6.5 };

  	const Int_t NPHIPTBINS = 14; 
  	const Double_t PHIXBINS[14+1] = {
  		0.5, 0.7, 0.9, 1.2, 
  		1.4, 1.6, 1.8, 2.0,
  		2.2, 2.6, 3.0, 3.5,
  		4.0, 5.0, 8.0 };

		const int nPtBins = 16;
		double ptBins[16+1] = { // divisible by omar's pions
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 5.0, 
		6.5, 8.0 };


	const Int_t NSPECIES = 6;
  	const Int_t NPTBINS[7+6] = {
  		nPtBins, nPtBins, nPtBins, nPtBins, nPtBins, nPtBins, nPtBins, 
  		NPIPTBINS, NPIPTBINS, NXIPTBINS, NXIPTBINS, NPHIPTBINS, NPHIPTBINS
  	};

  	const Double_t* XBINS[7+6] = {
  		ptBins, ptBins, ptBins, ptBins, ptBins, ptBins, ptBins,
  		PIXBINS, PIXBINS, XIXBINS, XIXBINS, PHIXBINS, PHIXBINS
  	};

		

		/*const int nPtBins = 20;
		double ptBins[20+1] = { // test to study trigger bin edge effect
		0.4, 0.6, 0.8, 1.0, 1.2, 
		1.4, 1.6, 1.8, 2.0, 2.2, 
		2.6, 3.0, 3.4, 4.0, 4.5, 4.8, 5.0, 5.2, 5.5, 
		6.5, 8.0 };*/

		const int nNchBins = 50;
		double NchBins[50+1] = {
			-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
			9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,
			19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5,
			29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5,
			39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5};

		const int nRTBins = 4;
		//double RTBins[4+1] = {0.0, 0.5, 1.5, 2.5, 5.0};
		double RTBins[4+1] = {0.0, 0.85, 1.5, 2.5, 5.0};

		/*const int nRTBins = 1;
		double RTBins[1+1] = {0.0, 5.0};*/

		const char* const PidArray[13]
			= {"Charged","Pion","Kaon","Proton","K0s","L","Lbar","piKp","pr","XiInc","Xi","phiInc","phi"};
		const char* const ArrayPidLatex[13]
			= {"h^{+} + h^{-}","#pi^{+} + #pi^{-}","K^{+} + K^{-}","p + #bar{p}", "K^{0}_{S}","#Lambda","#bar{#Lambda}",
				"#pi/K/p","p + #bar{p}","#Xi inc.","#Xi w. DCA cut","#phi inc.","#phi #rightarrow K^{+}K^{-}"};
		const char* const ArrayPidLatexPos[7]
			= {"h^{+}","#pi^{+}","K^{+}","p", "K^{0}_{S}","#Lambda","#bar{#Lambda}"};
		const char* const ArrayPidLatexNeg[7]
			= {"h^{-}","#pi^{-}","K^{-}","#bar{p}", "K^{0}_{S}", "#Lambda","#bar{#Lambda}"};
		const char* const EtaArray[4]
			= {"02","24","46","68"};

	protected:

		TH1F* _hUnfDist;
		TH1F* _hClosure;
		TH2F* _hPtvsNacc;
		TH2F* _hPtvsNch;
		TH3F* _hUnfMatrix;
		TH3F* _hCovMatrix;
		TH2F* _hSolutionNT;
		TH2F* _hCovMatrixNT;
		TH1F* _hEventRecEff;
		TH2F* _hUnfMatrixNT;

};

inline void UnfoldNTclass::H2V(const TH1F* h, TVectorD& v )
{
	for (int i = 1; i <= _nm; i++)
		(v)(i-1)= h->GetBinContent(i);
}

inline void UnfoldNTclass::H2VE()
{

	//! Fill vector with statistical uncertainties
	//! from the measured distribution
	for (int i = 0; i < _nm; i++)
		_eMes[i] = TMath::Sqrt( TMath::Abs(_meas->GetBinContent(i+1) ));

	//! Fill vector with statistical uncertainties
	//! from the true distribution
	for (int i = 0; i < _nt; i++)
		_eTrue[i] = TMath::Sqrt( TMath::Abs(_tru->GetBinContent(i+1) ));

}

inline void UnfoldNTclass::H2M(const TH2F* h, TMatrixD& m )
{
	//! Returns TMatrixD of the Response matrix
	//! The structure is the following:
	//! The variable that runs along the x-axis is the RECONSTRUCTED Mult.
	//! The variable that runs along the y-axis is the GENERATED Mult.

	for (int row = 1; row <= _nt; row++) {
		for (int col = 1; col <= _nm; col++) {
			(m)(row-1,col-1) = h->GetBinContent(col,row);
		}
	}
}

inline void UnfoldNTclass::V2H()
{

	for(int i = 0; i < _nt; i++){
		_hUnfDist->SetBinContent( i+1, _UnfDist(i) );
		_hUnfDist->SetBinError( i+1, TMath::Sqrt(_covUnf(i,i)) );
	}

}

inline void UnfoldNTclass::SetError(TH1F* h)
{
	for(int bin = 1; bin <= h->GetNbinsX(); bin++)
		h->SetBinError(bin,TMath::Sqrt(h->GetBinContent(bin)));
}

inline void UnfoldNTclass::SetError(TH2F* h)
{
	for(int biny = 1; biny <= h->GetNbinsY(); biny++){
		for(int binx = 1; binx <= h->GetNbinsX(); binx++){
			h->SetBinError(biny,binx,TMath::Sqrt(h->GetBinContent(biny,binx)));
		}
	}
}

inline void UnfoldNTclass::NormRowWise(TMatrixD& m)
{

	for (int row = 0; row < _nm; row++){

		double norm = 0.0;
		for (int col = 0; col < _nm; col++)
			norm += (m)(row,col);

		for (int col = 0; col < _nm; col++)
			if(norm != 0.0) (m)(row,col) = (m)(row,col)/norm;

	}
}

inline void UnfoldNTclass::NormRowWise(TVectorD& v)
{

	double norm = v.Norm1();
	v *= 1.0/norm;
}

inline void UnfoldNTclass::CheckNorm(TMatrixD& m)
{

	for (int row = 0; row < _nm; row++){
		double sum = 0.0;

		for (int col = 0; col < _nm; col++)
			sum += (m)(col,row);

		printf("Sum equals %f at rwo %d\n",sum,row);
	}
}


inline TH1F* UnfoldNTclass::GetUnfoldedDistH()
{

	return _hUnfDist;

}

inline void UnfoldNTclass::Emeasured()
{

	H2VE();

}

inline void UnfoldNTclass::SetRegion(const char* region)
{

	_Region = region;

}

inline void UnfoldNTclass::SetDetector(const char* detector)
{

	_Detector = detector;

}

inline void UnfoldNTclass::SetPid(const char* pid)
{

	_Pid = pid;

}

inline void UnfoldNTclass::SetCharge(const char* charge)
{

	_Charge = charge;

}

inline void UnfoldNTclass::SaveSolutionNT(bool saveNT)
{

	_SaveSolutionNT = saveNT;

}

inline TVectorD UnfoldNTclass::GetChi2Vectors(bool isXvex)
{
	if(isXvex) return _vChi2x;
	else return _vChi2y;

}

inline void UnfoldNTclass::myOptions(int lStat=1)
{
	// Set gStyle
	int font = 42;
	// From plain
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetStatColor(0);
	gStyle->SetStatBorderSize(1);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetDrawBorder(0);
	gStyle->SetTextFont(font);
	gStyle->SetStatFont(font);
	gStyle->SetStatFontSize(0.05);
	gStyle->SetStatX(0.97);
	gStyle->SetStatY(0.98);
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.3);
	gStyle->SetTickLength(0.02,"y");
	gStyle->SetEndErrorSize(3);
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetLabelFont(font,"xyz");
	gStyle->SetLabelOffset(0.01,"xyz");
	gStyle->SetTitleFont(font,"xyz");
	gStyle->SetTitleOffset(1.1,"xyz");
	gStyle->SetTitleSize(0.06,"xyz");
	gStyle->SetNdivisions(505,"y");
	gStyle->SetMarkerSize(1);
	gStyle->SetPalette(1);
	if (lStat){
		gStyle->SetOptTitle(1);
		gStyle->SetOptStat(1111);
		gStyle->SetOptFit(1111);
	}else {
		gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
	}
}


#endif
