#include "UnfoldNTclass.h"

UnfoldNTclass::UnfoldNTclass():

	_res(0x0),
	_meas(0x0),
	_tru(0x0),
	_hUnfDist(0x0),
	_hClosure(0x0),
	_niter(3),
	_nm(0),
	_nt(0),
	_Region("NT"),
	_Detector("TPC"),
	_PidIdx(0),
	_Pid("Charged"),
	_Charge("Charged"),
	_SaveSolutionNT(kFALSE),
	_hSolutionNT(0x0),
	_hCovMatrixNT(0x0),
	_hEventRecEff(0x0),
	_hUnfMatrixNT(0x0),
	isSys(kFALSE),
	TrkID("_TrkID_9"),
	ObjArray { new TObjArray() }

{ }

UnfoldNTclass::~UnfoldNTclass() {

	delete ObjArray;

}

void UnfoldNTclass::Setup(const TH1F* measured, const TH1F* truth, const TH2F* response )
{

	_meas = (TH1F*) measured->Clone();
	_tru = (TH1F*) truth->Clone();
	_res= (TH2F*) response->Clone();
	_nm = _res->GetNbinsX();
	_nt = _res->GetNbinsY();

	

	_P0C.ResizeTo(_nm);
	_nbarCi.ResizeTo(_nm);
	_Meas.ResizeTo(_nm);
	_UnfDist.ResizeTo(_nt);
	_UnfDist0.ResizeTo(_nt);
	_True.ResizeTo(_nt);
	_Res.ResizeTo(_nm,_nt);
	_Mij.ResizeTo(_nt,_nm);
	_Mij0.ResizeTo(_nt,_nm);
	_dnCidnEj.ResizeTo(_nt,_nm);
	_covUnf.ResizeTo(_nt,_nm);

	_vChi2x.ResizeTo(_niter);
	_vChi2y.ResizeTo(_niter);

	//! Vector with uncertainties
	_eMes.ResizeTo(_nm);
	_eTrue.ResizeTo(_nt);

	//! Covariance matrix of measurements
	_covMes.ResizeTo(_nm,_nm);

	

	//! Get the statistical uncertanties from the
	//! measured and true distributions into the
	//! vectors _eMes and _eTrue.
	Emeasured();

	//! Fills the covariance matrix with the 
	//! measured variances 
	GetMeasuredCov();
	
	//! Fills vectors and matrix 
	H2V( _meas, _Meas );
	H2V( _tru, _True );
	H2M( _res, _Res );
	
	NormRowWise( _Res );

	//! Uses the measured distribution 
	//! as the initila prior
	H2V(_meas,_P0C); 

	//! Uses the true distribution 
	//! as the initila prior
	// H2V(_tru,_P0C);

	//! Initial distribution
	_N0C= _P0C.Sum();
	if( _N0C != 0.0 ) {
		_P0C *= 1.0/_N0C;
	}

	_hUnfDist = (TH1F*) measured->Clone("_hUnfDist");
	_hUnfDist->Reset();

}
//_________________________
void UnfoldNTclass::LoadSolutionNT()
{

	TFile* fIn = nullptr;

	if(isMC) fIn = new TFile("results_unfolding/1D_newClass_mc.root","READ");
	else fIn = new TFile("results_unfolding/1D_newClass_data.root","READ");
	if(!fIn) printf("Unfolding matrix of the 1D NT NOT FOUND\n");

	if(!isMC) printf("\t Loading the NT solution from Data !! \n");
	if(isMC) printf("\t Loading the NT solution from MC !! \n");

	_hUnfMatrixNT = (TH2F*)fIn->Get("_hSolutionNT");
	_hUnfMatrixNT->Clone("");
	_hUnfMatrixNT->SetDirectory(0);

	delete fIn;
	fIn = nullptr;

}
//_________________________
void UnfoldNTclass::LoadSolutionNTMin()
{

	TFile* fIn = nullptr;

	if(isMC) fIn = new TFile("results_unfolding_min/1D_newClass_mc.root","READ");
	else fIn = new TFile("results_unfolding_min/1D_newClass_data.root","READ");
	if(!fIn) printf("Unfolding matrix of the 1D NT NOT FOUND\n");

	if(!isMC) printf("\t Loading the NT solution from Data !! \n");
	if(isMC) printf("\t Loading the NT solution from MC !! \n");

	_hUnfMatrixNT = (TH2F*)fIn->Get("_hSolutionNT");
	_hUnfMatrixNT->Clone("");
	_hUnfMatrixNT->SetDirectory(0);

	delete fIn;
	fIn = nullptr;

}
//_________________________
void UnfoldNTclass::LoadSolutionNTMax()
{

	TFile* fIn = nullptr;

	if(isMC) fIn = new TFile("results_unfolding_max/1D_newClass_mc.root","READ");
	else fIn = new TFile("results_unfolding_max/1D_newClass_data.root","READ");
	if(!fIn) printf("Unfolding matrix of the 1D NT NOT FOUND\n");

	if(!isMC) printf("\t Loading the NT solution from Data !! \n");
	if(isMC) printf("\t Loading the NT solution from MC !! \n");

	_hUnfMatrixNT = (TH2F*)fIn->Get("_hSolutionNT");
	_hUnfMatrixNT->Clone("");
	_hUnfMatrixNT->SetDirectory(0);

	delete fIn;
	fIn = nullptr;

}
//_________________________
void UnfoldNTclass::UnfoldV02D(TH2F* hIniRM, TH2F* hPtRecvsNacc, TH2F* hPIDPtRecvsNacc)
{

	cout << "binning is " << NPTBINS[_PidIdx] << endl;

	_hPtvsNacc = new TH2F("_hPtvsNacc",";#it{p}_{T}^{rec} (GeV/#it{c});#it{N}_{m}",NPTBINS[_PidIdx],XBINS[_PidIdx],nNchBins,NchBins);
	TH2F* hPtPIDvsNacc = new TH2F("hPtPIDvsNacc",";#it{p}_{T}^{rec} (GeV/#it{c});#it{N}_{m}",NPTBINS[_PidIdx],XBINS[_PidIdx],nNchBins,NchBins);
	_hPtvsNch = new TH2F("_hPtvsNch",";#it{p}_{T} (GeV/#it{c});#it{N}_{t}",NPTBINS[_PidIdx],XBINS[_PidIdx],nNchBins,NchBins);
	cout << "nbins " << _hPtvsNch->GetNbinsX() << endl;
	_hUnfMatrix = new TH3F("_hUnfMatrix",";#it{N}_{ch};#it{N}_{acc};#it{p}_{T}",nNchBins,NchBins,nNchBins,NchBins,NPTBINS[_PidIdx],XBINS[_PidIdx]);
	_hCovMatrix = new TH3F("_hCovMatrix",";#it{N}_{ch};#it{N}_{acc};#it{p}_{T}",nNchBins,NchBins,nNchBins,NchBins,NPTBINS[_PidIdx],XBINS[_PidIdx]);

	ObjArray->AddLast(_hPtvsNacc);
	ObjArray->AddLast(hPtPIDvsNacc);
	ObjArray->AddLast(_hPtvsNch);

	//! Only for the Transverse region
	//! the Unfolding matrix in pT bins
	//! is stored
	if(strcmp(_Region,"Transverse")==0) {
		ObjArray->AddLast(_hUnfMatrix);
		ObjArray->AddLast(_hCovMatrix);
	} 

	TH2F* hRM = (TH2F*)hIniRM->Clone("hRM");
	TString rmname(hIniRM->GetName());

	if((strcmp(_Region,"Toward")==0) || (strcmp(_Region,"Away")==0) || (strcmp(_Region,"Trans1D")==0)
		|| (strcmp(_Region,"TransMin1D")==0) || (strcmp(_Region,"TransMax1D")==0)) {
		if (rmname.Contains("Min")) LoadSolutionNTMin();
		else if (rmname.Contains("Max")) LoadSolutionNTMax();
		else LoadSolutionNT();
	}
	//if((strcmp(_Region,"Transverse")==0)) LoadSolutionNT();

	if(strcmp(_Region,"Transverse")==0){
		for(int binx = 1; binx <= hIniRM->GetNbinsX(); binx++){
			for(int biny = 1; biny <= hIniRM->GetNbinsY(); biny++){
				//hRM->SetBinContent(binx,biny,hIniRM->GetBinContent(binx,biny)*1.);
				//hRM->SetBinContent(binx,biny,hIniRM->GetBinContent(binx,biny)*(binx-1));
				hRM->SetBinContent(binx,biny,hIniRM->GetBinContent(binx,biny)*( binx>1 ? binx-1 : 0.000)); // careful: here nt=0 bin alsogets emptied!
			}
		}
	}

	// NCHARGED
	for(int binNch = 1; binNch <= nNchBins; binNch++){

		TH1D* hPt = (TH1D*)hPtRecvsNacc->ProjectionX(Form("hPt_Corrected_Bin_nch_%d",binNch),binNch,binNch);

		for(int bin_pt = 1; bin_pt <= NPTBINS[_PidIdx]; bin_pt++){
			_hPtvsNacc->SetBinContent(bin_pt,binNch,hPt->GetBinContent(bin_pt));
			_hPtvsNacc->SetBinError(bin_pt,binNch,hPt->GetBinError(bin_pt));
		} 
	}

	// PID
	for(int binNch = 1; binNch <= nNchBins; binNch++){

		TH1D* hPt = (hPIDPtRecvsNacc) ? (TH1D*)hPIDPtRecvsNacc->ProjectionX(Form("hPtPID_Corrected_Bin_nch_%d",binNch),binNch,binNch) : 0x0;

		if (hPIDPtRecvsNacc) for(int bin_pt = 1; bin_pt <= NPTBINS[_PidIdx]; bin_pt++){
			hPtPIDvsNacc->SetBinContent(bin_pt,binNch,hPt->GetBinContent(bin_pt));
			hPtPIDvsNacc->SetBinError(bin_pt,binNch,hPt->GetBinError(bin_pt));
		} 
	}


	//! Unfold NT distributions in pT bins
	for(int bin_pt = 1; bin_pt <= NPTBINS[_PidIdx]; bin_pt++){

		double lowBound = _hPtvsNacc->GetXaxis()->GetBinLowEdge(bin_pt);
		double upBound  = _hPtvsNacc->GetXaxis()->GetBinUpEdge(bin_pt);

		TH1F* hProj  = (TH1F*)_hPtvsNacc->ProjectionY(Form("hNch_Unfolded_pt_%d",bin_pt),bin_pt,bin_pt);
		TH1F* hPIDProj  = (hPIDPtRecvsNacc) ? (TH1F*)hPtPIDvsNacc->ProjectionY(Form("hNch_PID_Unfolded_pt_%d",bin_pt),bin_pt,bin_pt) : 0x0;
		TH1F* hPrior = (TH1F*)_hPtvsNacc->ProjectionY(Form("hNch_Prior_pt_%d",bin_pt),bin_pt,bin_pt);
		//TH1F* hMultGen = (TH1F*)_hPtvsNacc->ProjectionX(Form("hGen_Prior_%d",bin_pt),bin_pt,bin_pt);
		TH1F* hMultGen = (TH1F*)_hPtvsNacc->ProjectionY(Form("hGen_Prior_%d",bin_pt),bin_pt,bin_pt);
		if(hProj->Integral()==0.0) continue;

		printf("Unfolding region %s ptbin %i \n", _Region, bin_pt);

		Setup(hProj, hMultGen, hRM );
		///Setup(hProj, hProj, hRM );
		if(strcmp(_Region,"Transverse")==0) { Unfold(); }		// Unfolding matrix is calculated from Nch
		if (hPIDProj) UnfoldPID(hPIDProj, hMultGen, hRM);	// Recalculating unfolded distribution from PID spectra instead of Nch
		
		if((strcmp(_Region,"Toward")==0) || (strcmp(_Region,"Away")==0) || (strcmp(_Region,"Trans1D")==0)
			|| (strcmp(_Region,"TransMin1D")==0) || (strcmp(_Region,"TransMax1D")==0) ) {
			//Setup(hProj, hMultGen, hRM );
			UnfoldTowardAway();
		}
		//UnfoldTowardAway();

		V2H();
		TH1F* hUnfolded = (TH1F*)GetUnfoldedDistH();
		hUnfolded->SetName(Form("hUnfolded_%s_%d",_Region,bin_pt));
		hUnfolded->SetTitle(Form("%.2f < #it{p}_{T} < %.2f",lowBound,upBound));
		ObjArray->AddLast(hUnfolded);

		//! Here I fill the Unfolding matrix in a 3D Histo.
		//! this Unfoldin matrix is only applicable
		//! in the Transverse region
		GetUnfoldingMatrix3D(bin_pt);

		//! Here I fill the Covariance matrices in a 3D Histo.
		//! this is applicable in the Transverse region.
		//! For the Toward and Away regions, I use the 
		//! Covariance matrix from the NT case
		FillCovarianceMatrix(bin_pt);

		for(int bin_nch = 1; bin_nch <= nNchBins; bin_nch++){
			if(hUnfolded->GetBinContent(bin_nch) >= 0.0){
				//double bin_width = _hPtvsNch->GetYaxis()->GetBinWidth(bin_pt);
				//_hPtvsNch->SetBinContent(bin_pt,bin_nch,hUnfolded->GetBinContent(bin_nch));
				//_hPtvsNch->SetBinError(bin_pt,bin_nch,hUnfolded->GetBinError(bin_nch));

				double bin_width = _hPtvsNch->GetYaxis()->GetBinWidth(bin_pt);
				_hPtvsNch->SetBinContent(bin_pt,bin_nch,hUnfolded->GetBinContent(bin_nch));
				_hPtvsNch->SetBinError(bin_pt,bin_nch,hUnfolded->GetBinError(bin_nch));


			}
		}

	}	

}

void UnfoldNTclass::Unfold2D(TH2F* hIniRM, TH2F* hPtRecvsNacc)
{

	_hPtvsNacc = new TH2F("_hPtvsNacc",";#it{p}_{T}^{rec} (GeV/#it{c});#it{N}_{m}",NPTBINS[_PidIdx],XBINS[_PidIdx],nNchBins,NchBins);
	_hPtvsNch = new TH2F("_hPtvsNch",";#it{p}_{T} (GeV/#it{c});#it{N}_{t}",NPTBINS[_PidIdx],XBINS[_PidIdx],nNchBins,NchBins);
	_hUnfMatrix = new TH3F("_hUnfMatrix",";#it{N}_{ch};#it{N}_{acc};#it{p}_{T}",nNchBins,NchBins,nNchBins,NchBins,NPTBINS[_PidIdx],XBINS[_PidIdx]);
	_hCovMatrix = new TH3F("_hCovMatrix",";#it{N}_{ch};#it{N}_{acc};#it{p}_{T}",nNchBins,NchBins,nNchBins,NchBins,NPTBINS[_PidIdx],XBINS[_PidIdx]);

	ObjArray->AddLast(_hPtvsNacc);
	ObjArray->AddLast(_hPtvsNch);

	//! Only for the Transverse region
	//! the Unfolding matrix in pT bins
	//! is stored
	if(strcmp(_Region,"Transverse")==0) {
		ObjArray->AddLast(_hUnfMatrix);
		ObjArray->AddLast(_hCovMatrix);
	} 

	TH2F* hRM = (TH2F*)hIniRM->Clone("hRM");

	if((strcmp(_Region,"Toward")==0) || (strcmp(_Region,"Away")==0) || (strcmp(_Region,"Trans1D")==0 )
		|| (strcmp(_Region,"TransMin1D")==0 ) || (strcmp(_Region,"TransMax1D")==0 ) ) LoadSolutionNT();

	if(strcmp(_Region,"Transverse")==0){

		// Multiplying each column N_T,m with the value of N_T,m
		for(int binx = 1; binx <= hIniRM->GetNbinsX(); binx++){
			for(int biny = 1; biny <= hIniRM->GetNbinsY(); biny++){
				hRM->SetBinContent(binx,biny,hIniRM->GetBinContent(binx,biny)*(binx-1));
			}
		}
	}

	//! Correct pT distributions by 
	//! tracking inefficiency and 
	//! secondary contamination
	CorrectPTdistributions(hPtRecvsNacc); // This fills corrected hPtRecvsNacc into _hPtvsNacc

	//! Unfold NT distributions in pT bins
	for(int bin_pt = 1; bin_pt <= NPTBINS[_PidIdx]; bin_pt++){

		double lowBound = _hPtvsNacc->GetXaxis()->GetBinLowEdge(bin_pt);
		double upBound  = _hPtvsNacc->GetXaxis()->GetBinUpEdge(bin_pt);

		TH1F* hProj  = (TH1F*)_hPtvsNacc->ProjectionY(Form("hNch_Unfolded_pt_%d",bin_pt),bin_pt,bin_pt); //nch
		TH1F* hPrior = (TH1F*)_hPtvsNacc->ProjectionY(Form("hNch_Prior_pt_%d",bin_pt),bin_pt,bin_pt); //nch
		TH1F* hMultGen = (TH1F*)_hPtvsNacc->ProjectionX(Form("hGen_Prior_%d",bin_pt),bin_pt,bin_pt); //pt
		if(hProj->Integral()==0.0) continue;

		Setup(hProj, hMultGen, hRM );
		if(strcmp(_Region,"Transverse")==0) Unfold();
		if((strcmp(_Region,"Toward")==0) || (strcmp(_Region,"Away")==0) || (strcmp(_Region,"Trans1D")==0)
			|| (strcmp(_Region,"TransMin1D")==0) || (strcmp(_Region,"TransMax1D")==0) ) UnfoldTowardAway();

		V2H();
		TH1F* hUnfolded = (TH1F*)GetUnfoldedDistH();
		hUnfolded->SetName(Form("hUnfolded_%s_%d",_Region,bin_pt));
		hUnfolded->SetTitle(Form("%.2f < #it{p}_{T} < %.2f",lowBound,upBound));
		//ObjArray->AddLast(hUnfolded);

		//! Here I fill the Unfolding matrix in a 3D Histo.
		//! this Unfoldin matrix is only applicable
		//! in the Transverse region
		GetUnfoldingMatrix3D(bin_pt);

		//! Here I fill the Covariance matrices in a 3D Histo.
		//! this is applicable in the Transverse region.
		//! For the Toward and Away regions, I use the 
		//! Covariance matrix from the NT case
		FillCovarianceMatrix(bin_pt);

		for(int bin_nch = 1; bin_nch <= nNchBins; bin_nch++){
			if(hUnfolded->GetBinContent(bin_nch) >= 0.0){
				double bin_width = _hPtvsNch->GetYaxis()->GetBinWidth(bin_pt);
				_hPtvsNch->SetBinContent(bin_pt,bin_nch,hUnfolded->GetBinContent(bin_nch));
				_hPtvsNch->SetBinError(bin_pt,bin_nch,hUnfolded->GetBinError(bin_nch));
			}
		}

	}	

}

//_________________________
void UnfoldNTclass::Unfold()
{

	for(int iter = 0; iter < _niter; ++iter){

		//cout << "Iteration : " << iter << endl;

		if(iter > 0){
			//! _P0C Update of the prior
			_P0C = _nbarCi;
			//! _N0C: total number of events
			//! in the unfolded distribution
			_N0C = _UnfDist.Norm1();
		}

		//! Computes the Unfolding matrix: _Mij
		GetUnfoldingMatrix();

		//! Gets the unfolded distribution: _UnfDist
		//! This si given as SUM_{ij} _M_{ij} * _Meas_{j}
		GetUnfoldedDist();

		//! Get Chi2
		_vChi2x(iter) = iter + 1;
		_vChi2y(iter) = GetChi2();

		//! Get the new estimate of true distribution
		GetNewPrior();

		//! Convariance matrix from the unfolded distribution 
		//! At the first iteration the covarianve matrix is equal 
		//! to the unfolding matrix from the first iteration 
		if(iter <= 0) {
			_dnCidnEj = _Mij;
		} else{

			TVectorD en(_nm), nr(_nm);
			for (Int_t i = 0 ; i < _nm ; i++) {
				if( _P0C[i] <= 0.0 ) continue;
				Double_t ni = 1.0 / ( _N0C * _P0C[i] );
				en[i]= -1.0 * ni;
				nr[i]=  ni * _UnfDist[i] ;
			}

			TMatrixD M1 = _dnCidnEj;
			M1.NormByColumn(nr,"M");

			TMatrixD dSum;
			dSum.ResizeTo(_nt,_nm);
			double sum = 0.0;

			for(int i = 0; i < _nt; i++){

				for(int j = 0; j < _nm; j++){

					for(int k = 0; k < _nm; k++){

						sum = 0.0;

						for(int l = 0; l < _nt; l++){

							if( _P0C[l] > 0.0 )
								sum += _Meas[k] * ( 1.0 / _N0C * _P0C[l] ) * _Mij(i,k) * _Mij(l,k) * _dnCidnEj(l,j); 

						}
					}
					dSum(i,j) = -1.0 * sum;
				}
			}

			_dnCidnEj = _Mij;
			_dnCidnEj += M1;
			_dnCidnEj += dSum;

		} //! Compute _dnCidnEj for iter > 0 

	} //! iter-loop

	GetCovMatrixUnfolded();

	//! Here I save the unfolding solution of NT
	//! which will be used to unfold the Toward/Away regions
	//! _SaveSolutionNT==kTRUE only when unfold the NT distribution
	if(_SaveSolutionNT){

		_hSolutionNT = new TH2F("_hSolutionNT",";#it{N}_{ch};#it{N}_{acc}",nNchBins,NchBins,nNchBins,NchBins);	
		_hCovMatrixNT = new TH2F("_hCovMatrixNT",";#it{N}_{ch};#it{N}_{acc}",nNchBins,NchBins,nNchBins,NchBins);	

		GetUnfoldingMatrix2D();
		ObjArray->AddLast(_hSolutionNT);

		GetCovarianceMatrix2D();
		ObjArray->AddLast(_hCovMatrixNT);

	} 

}
//_________________________
void UnfoldNTclass::UnfoldTowardAway()
{
	GetUnfoldedDist();

}
//_________________________
void UnfoldNTclass::GetUnfoldingMatrix()
{

	//! Response Matrix
	//! True NT runs downward along y-axis
	//! Rec  NT runs to the right along x-axis

	for(int j = 0; j < _nm; ++j){

		double weight = 0.0;
		for(int l = 0; l < _nt; l++)
			weight += (_Res)(l,j)*(_P0C)(l);

		for(int i = 0; i < _nt; ++i){

			double response = 0.0;
			double prior    = 0.0;

			response = (_Res)(i,j);
			prior    = (_P0C)(i);

			if(weight != 0.0) (_Mij)(j,i) = (response*prior)/weight;
		}
	}

}
//_________________________
void UnfoldNTclass::FillMatrixFromNTSolution()
{

	for(int i = 0; i < _nm; i++){

		for(int j = 0; j < _nm; j++)

			(_Mij)(j,i) = _hUnfMatrixNT->GetBinContent(i+1,j+1);
	}

}
//_________________________
void UnfoldNTclass::GetUnfoldedDist()
{

	if((strcmp(_Region,"Toward")==0) || (strcmp(_Region,"Away")==0) || (strcmp(_Region,"Trans1D")==0 )
	 || (strcmp(_Region,"TransMin1D")==0 )  || (strcmp(_Region,"TransMax1D")==0 ) )// || (strcmp(_Region,"Transverse")==0))
		FillMatrixFromNTSolution();

	//if(strcmp(_Region,"Toward")==0)
	//_Mij.Print();

	for(int i = 0; i < _nm; i++){

		double sum = 0.0;
		for(int j = 0; j < _nm; j++)
			sum += (_Mij)(j,i)*(_Meas)(j);

		(_UnfDist)(i) = sum;

	}

}
//_________________________
void UnfoldNTclass::UnfoldPID(const TH1F* measured, const TH1F* truth, const TH2F* response )
{

	_meas = (TH1F*) measured->Clone();
	_tru = (TH1F*) truth->Clone();
	_res= (TH2F*) response->Clone();
	_nm = _res->GetNbinsX();

	_Meas.ResizeTo(_nm);
	_True.ResizeTo(_nt);
	_Res.ResizeTo(_nm,_nt);

	//! Vector with uncertainties
	_eMes.ResizeTo(_nm);
	_eTrue.ResizeTo(_nt);

	//! Get the statistical uncertanties from the
	//! measured and true distributions into the
	//! vectors _eMes and _eTrue.
	Emeasured();
	
	//! Fills vectors and matrix 
	H2V( _meas, _Meas );
	H2V( _tru, _True );
	H2M( _res, _Res );
	
	//NormRowWise( _Res );


	for(int i = 0; i < _nm; i++){

		double sum = 0.0;
		for(int j = 0; j < _nm; j++) {
			//if (i==5&&j==5) printf("i %i j %i Mij %f Meas %f \n", i, j, (_Mij)(j,i), (_Meas)(j));
			sum += (_Mij)(j,i)*(_Meas)(j);
		}

		(_UnfDist)(i) = sum;

	}

}

//_________________________
double UnfoldNTclass::GetChi2()
{

	double nu = 0.;
	double ui2 = 0.0;
	double Sumui2 = 0.0;

	for(int i = 0; i < _nt; i++){

		double n = _True(i);
		double ni = _UnfDist(i);
		if( _eTrue[i] <= 0 ) continue;
		nu++;
		ui2 = TMath::Power( (ni - n), 2. ) / TMath::Power( _eTrue[i], 2. );

		Sumui2 += ui2;
	}

	return Sumui2 / nu;

} 
//_________________________
void UnfoldNTclass::GetNewPrior()
{

	double norm = _UnfDist.Norm1();

	_nbarCi = _UnfDist;

	//_nbartrue = norm;

	_nbarCi *= 1.0/norm;

}
//_________________________
void UnfoldNTclass::CorrectPTdistributions(const TH2F* hPtRecvsNacc)
{

	TFile* fEff = nullptr;

	if(!isSys) fEff = new TFile("~/Documents/DataAnalysis/pp_13TeV_spectra_RT/analysis/efficiencies/final_one_hybrid/MCEfficiencies_MB_PYTHIA.root","READ");
	else fEff = new TFile(Form("~/Documents/DataAnalysis/pp_13TeV_spectra_RT/analysis/efficiencies/final_one_hybrid/systematics/MCEfficiencies_MB_PYTHIA%s.root",TrkID.c_str()),"READ");
	if(!fEff)printf("no eff\n");

	TDirectory* dir  = (TDirectory*)fEff->Get("Tracking_Matching_Efficiencies");
	TDirectory* dNch = (TDirectory*)fEff->Get("Nch_Efficiencies");

	TH1D* hEff = (TH1D*)dir->Get(Form("%s_Eff_%s_%s_MB_Corrected",_Detector,_Pid,_Charge));
	TH1D* hSecCont = (TH1D*)dNch->Get("hChargedPrimaries");
	TH1D* hPiFD = (TH1D*)dNch->Get("hPionPrimaries");
	TH1D* hPFD = (TH1D*)dNch->Get("hProtonPrimaries");

	for(int binNch = 1; binNch <= nNchBins; binNch++){

		TH1D* hPt = (TH1D*)hPtRecvsNacc->ProjectionY(Form("hPt_Corrected_Bin_nch_%d",binNch),binNch,binNch);

		for(int bin_pt = 1; bin_pt <= NPTBINS[_PidIdx]; bin_pt++)
			hPt->SetBinError(bin_pt,TMath::Sqrt(hPt->GetBinContent(bin_pt)));

		//! Efficiency correction
		hPt->Divide(hEff);

		//! FD correction. Only charged particles
		if(isMC) hPt->Multiply(hSecCont);

		for(int bin_pt = 1; bin_pt <= NPTBINS[_PidIdx]; bin_pt++){
			_hPtvsNacc->SetBinContent(bin_pt,binNch,hPt->GetBinContent(bin_pt));
			_hPtvsNacc->SetBinError(bin_pt,binNch,hPt->GetBinError(bin_pt));
		} 
	}

}
//_________________________
void UnfoldNTclass::PerformClosure()
{

	_hClosure = (TH1F*)_hUnfDist->Clone("_hClosure");
	_hClosure->Divide(_tru);

}

void UnfoldNTclass::GetMCclosureinNchBins(const TH2F* h2D)
{

	for(int binnch = 1; binnch <= nNchBins; binnch++){

		//TH1F* hg = RebinPtAxis((TH1F*)h2D->ProjectionY(Form("hPtGen_Nch_%d",binnch),binnch,binnch));
		TH1F* hg = (TH1F*)h2D->ProjectionX(Form("hPtGen_Nch_%d",binnch),binnch,binnch);
		
		TH1F* hu = (TH1F*)_hPtvsNch->ProjectionX(Form("hPt_Corrected_Bin_nch_%d",binnch),binnch,binnch);
		cout << "nbins 2 " << hu->GetNbinsX() << endl;

		if((strcmp(_Region,"Toward")==0)||(strcmp(_Region,"Transverse")==0)){
			if(binnch==2||binnch==6||binnch==11||binnch==16||binnch==21)
				DrawClosureRatios(hg,hu,0.0,10.12,binnch,"#it{p}_{T} (GeV/#it{c})","Unfolded/True",h2D->GetXaxis()->GetBinLowEdge(binnch),h2D->GetXaxis()->GetBinUpEdge(binnch),binnch,kTRUE,kFALSE,"results_unfolding");
		}
	}
}

void UnfoldNTclass::GetMCclosureinRTBins(const TH2F* hGen, const TH2F* hUnf)
{

	for(int binrt = 0; binrt <= nRTBins; binrt++){

		int lowedge = GetRTBin(binrt,kTRUE);
		int upedge  = GetRTBin(binrt,kFALSE);

		TH1F* hg = (TH1F*)hGen->ProjectionX(Form("hPtGen_RT_%d",binrt),lowedge,upedge);
		
		TH1F* hu = (TH1F*)hUnf->ProjectionX(Form("hPtUnf_RT_%d",binrt),lowedge,upedge);

		cout << "2 Finding " << hg << " and " << hu << "\n";

		//	if((strcmp(_Region,"Toward")==0)||(strcmp(_Region,"Transverse")==0)){
		//DrawClosureRatios(hg,hu,0.0,10.12,binrt,"#it{p}_{T} (GeV/#it{c})","Unfolded/True",RTBins[binrt-1],RTBins[binrt],binrt,kFALSE,kTRUE,"results");
		DrawClosureRatios(hg,hu,0.0,8.12,binrt,"#it{p}_{T} (GeV/#it{c})","Unfolded/True",(binrt>0)?RTBins[binrt-1]:RTBins[0],(binrt>0)?RTBins[binrt]:RTBins[nRTBins],binrt,kFALSE,kTRUE,"results_unfolding");
		//}
	}
}

void UnfoldNTclass::GetMCclosureinRTMinBins(const TH2F* hGen, const TH2F* hUnf)
{

	for(int binrt = 0; binrt <= nRTBins; binrt++){

		int lowedge = GetRTMinBin(binrt,kTRUE);
		int upedge  = GetRTMinBin(binrt,kFALSE);

		TH1F* hg = (TH1F*)hGen->ProjectionX(Form("hPtGen_RT_%d",binrt),lowedge,upedge);
		
		TH1F* hu = (TH1F*)hUnf->ProjectionX(Form("hPtUnf_RT_%d",binrt),lowedge,upedge);

		cout << "2 Finding " << hg << " and " << hu << "\n";

		//	if((strcmp(_Region,"Toward")==0)||(strcmp(_Region,"Transverse")==0)){
		//DrawClosureRatios(hg,hu,0.0,10.12,binrt,"#it{p}_{T} (GeV/#it{c})","Unfolded/True",RTBins[binrt-1],RTBins[binrt],binrt,kFALSE,kTRUE,"results");
		DrawClosureRatios(hg,hu,0.0,8.12,binrt,"#it{p}_{T} (GeV/#it{c})","Unfolded/True",(binrt>0)?RTBins[binrt-1]:RTBins[0],(binrt>0)?RTBins[binrt]:RTBins[nRTBins],binrt,kFALSE,kTRUE,"results_unfolding_min");
		//}
	}
}

void UnfoldNTclass::GetMCclosureinRTMaxBins(const TH2F* hGen, const TH2F* hUnf)
{

	for(int binrt = 0; binrt <= nRTBins; binrt++){

		int lowedge = GetRTMaxBin(binrt,kTRUE);
		int upedge  = GetRTMaxBin(binrt,kFALSE);

		TH1F* hg = (TH1F*)hGen->ProjectionX(Form("hPtGen_RT_%d",binrt),lowedge,upedge);
		
		TH1F* hu = (TH1F*)hUnf->ProjectionX(Form("hPtUnf_RT_%d",binrt),lowedge,upedge);

		cout << "2 Finding " << hg << " and " << hu << "\n";

		//	if((strcmp(_Region,"Toward")==0)||(strcmp(_Region,"Transverse")==0)){
		//DrawClosureRatios(hg,hu,0.0,10.12,binrt,"#it{p}_{T} (GeV/#it{c})","Unfolded/True",RTBins[binrt-1],RTBins[binrt],binrt,kFALSE,kTRUE,"results");
		DrawClosureRatios(hg,hu,0.0,8.12,binrt,"#it{p}_{T} (GeV/#it{c})","Unfolded/True",(binrt>0)?RTBins[binrt-1]:RTBins[0],(binrt>0)?RTBins[binrt]:RTBins[nRTBins],binrt,kFALSE,kTRUE,"results_unfolding_max");
		//}
	}
}

TH1F* UnfoldNTclass::GetClosureH()
{

	PerformClosure();

	return _hClosure;

}

void UnfoldNTclass::GetMeasuredCov()
{

	//! Fills covariance matrix on measured distribution
	for (int i= 0 ; i <_nm; i++) {
		for (int j= 0 ; j <_nm; j++) {
			if( i == j ) (_covMes)(i,j) = _eMes[i] * _eMes[i];
			else (_covMes)(i,j) = 0.0;
		}
	}

}

void UnfoldNTclass::GetCovMatrixUnfolded()
{

	TMatrixD _dnCidnEjT(TMatrixD::kTransposed,_dnCidnEj);

	TMatrixD M1;
	M1.ResizeTo(_nm,_nm);
	M1.Mult(_dnCidnEj,_covMes);

	_covUnf.Mult(M1,_dnCidnEjT);

}

TH1F* UnfoldNTclass::RebinNT2RT(const TH1F* hIn, bool isSta)
{

	const int nRTBins_Re = 45;
	const double MeanNT = 7.356;

	double RTBins_Re[nRTBins_Re+1] = { 0.0 };
	for(int bin = 0; bin <= nRTBins_Re; bin++)
		RTBins_Re[bin] = (  (double)bin ) / MeanNT;

	TH1F* hReturn = new TH1F("_hRT","",nRTBins_Re,RTBins_Re);

	for(int i = 1; i <= hReturn->GetNbinsX(); i++){

		double lowedge = hReturn->GetBinLowEdge(i);
		double upedge  = hReturn->GetBinLowEdge(i) + hReturn->GetBinWidth(i);

		double y = 0.0;
		double ey = 0.0;
		int counter = 0;

		for(int bin = 1; bin <= hIn->GetNbinsX(); bin++){

			double NT = hIn->GetBinCenter(bin);
			double RT = NT / MeanNT;


			if( ( RT >= lowedge ) && ( RT < upedge ) ){

			if( i < 3) printf("low = %f RT = %f high = %f\n", lowedge, RT, upedge);

			y  += hIn->GetBinContent(bin);
			if(isSta) ey += hIn->GetBinError(bin) * hIn->GetBinError(bin);
			if(!isSta) ey += hIn->GetBinError(bin);

				}
		}

		hReturn->SetBinContent(i,y);
		if(isSta) hReturn->SetBinError(i,TMath::Sqrt(ey));
		if(!isSta) hReturn->SetBinError(i, ey / (double)counter );
	}

	return hReturn;

}

TH1F* UnfoldNTclass::RebinPtAxis(const TH1F* hIn)
{

	const int nPtBins_Re = 7;
	double ptBins_Re[nPtBins_Re+1] = { 0.0, 0.3, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0 };

	TH1F* hReturn = new TH1F("","",nPtBins_Re,ptBins_Re);
	hReturn->SetName(Form("%s_Re",hIn->GetName()));

	for(int i = 1; i <= hReturn->GetNbinsX(); i++){

		double low_edge = hReturn->GetBinLowEdge(i);
		double up_edge  = hReturn->GetBinLowEdge(i)+hReturn->GetBinWidth(i);

		double y = 0.0;
		double ey = 0.0;

		for(int bin = 1; bin <= hIn->GetNbinsX(); bin++){

			double pt = hIn->GetBinCenter(bin);

			if((pt<low_edge) || (pt>up_edge)) continue;

			y += hIn->GetBinContent(bin);
			ey += hIn->GetBinError(bin)*hIn->GetBinError(bin);

		}

		hReturn->SetBinContent(i,y);
		hReturn->SetBinError(i,TMath::Sqrt(ey));
	}

	return hReturn;

}

void UnfoldNTclass::DrawNchClosure(const TH1F* hGen,const TH1F* hUn, const double& min, const double& max, const char* titleX, const char* titleY, const char* dirOut)
{

	myOptions(0);
	TCanvas* cr = new TCanvas("Cr_Charged","",200,10,500,480);

	TPad* pad2 = GetPad();

	TH1F* hfrat  = TuneFrame(min,max,titleX,titleY);

	TLegend *legTmp = new TLegend(0.2,0.2,0.6,0.45);
	myLegendSetUp(legTmp,0.05);

	TLatex* latex = new TLatex();
	TuneLatex(latex);

	TH1D* hUnToGen  = (TH1D*)hUn->Clone("UnToGen");
	hUnToGen->Divide(hGen);
	hUnToGen->SetMarkerStyle(8);
	hUnToGen->SetMarkerColor(kBlack);
	hUnToGen->SetLineColor(kBlack);

	TLine* line0 = new TLine(min,1.0,max,1.0);
	TLine* line1 = new TLine(min,1.1,max,1.1);
	TLine* line2 = new TLine(min,0.9,max,0.9);
	line0->SetLineColor(1);
	line0->SetLineStyle(2);
	line1->SetLineColor(1);
	line1->SetLineStyle(2);
	line2->SetLineColor(1);
	line2->SetLineStyle(2);

	cr->cd();
	pad2->Draw();
	pad2->cd();
	pad2->SetFillColor(0);
	hfrat->Draw();
	hUnToGen->Draw("SAME P");
	line0->Draw("SAME");
	line1->Draw("SAME");
	line2->Draw("SAME");

	cr->SaveAs(Form("./%s/Nch_Closure.pdf",dirOut));
}

void UnfoldNTclass::DrawClosureRatios(TH1F* hGen,TH1F* hUn, const double& min, const double& max, const int& bin, const char* titleX, const char* titleY, const double& lowerBound, const double& upperBound, const int& binNch, const bool isNch, const bool isRT, const char* dirOut)
{

	myOptions(0);
	TCanvas* cr = 0x0;
	TPad* pad2 = GetPad();
	TH1F* hfrat  = TuneFrame(min,max,titleX,titleY);

	if(isNch){
		cr = new TCanvas(Form("Cr_%s_pid_%s_Nch_%d",_Region,_Pid,bin),"",200,10,500,480);
		hfrat->SetName(Form("Cr_%s_pid_%s_Nch_%d",_Region,_Pid,bin));                
	} else if(isRT){
		cr = new TCanvas(Form("Cr_%s_pid_%s_RT_%d",_Region,_Pid,bin),"",200,10,500,480);
		hfrat->SetName(Form("Cr_%s_pid_%s_RT_%d",_Region,_Pid,bin));                
	} else{
		cr = new TCanvas(Form("Cr_%s_pid_%s_pT_%d",_Region,_Pid,bin),"",200,10,500,480);
		hfrat->SetName(Form("Cr_%s_pid_%s_pT_%d",_Region,_Pid,bin));
	}

	TLegend *legTmp = new TLegend(0.2,0.2,0.6,0.45);
	myLegendSetUp(legTmp,0.05);

	TLatex* latex = new TLatex();
	TuneLatex(latex);

	TH1D* hUnToGen  = (TH1D*)hUn->Clone(Form("hUnToGen_%s_%s_%i",_Pid,_Region,bin));
	cout << "3.0 Finding " << hUnToGen << " " << hUnToGen->GetBinContent(8) << " and " << hUnToGen->GetNbinsX() << "\n";
	cout << "3.1 Finding " << hGen << " " << hGen->GetBinContent(8) << " and " << hGen->GetNbinsX() << "\n";
	hUnToGen->Divide(hGen);
	hUnToGen->SetMarkerStyle(8);
	hUnToGen->SetMarkerColor(kBlack);
	hUnToGen->SetLineColor(kBlack);

	cout << "3.15 Finding " << hfrat << "\n";

	cr->cd();
	pad2->Draw();
	pad2->cd();
	pad2->SetFillColor(0);
	hfrat->Draw();
	cout << "3.2 Finding " << hUnToGen << " " << hUnToGen->GetBinContent(8) << "\n";
	hUnToGen->Draw( strcmp(_Region,"Transverse")==0 ? "hist p same" : "same");

	TLine* line0 = new TLine(min,1.0,max,1.0);
	TLine* line1 = new TLine(min,1.1,max,1.1);
	TLine* line2 = new TLine(min,0.9,max,0.9);
	line0->SetLineColor(1);
	line0->SetLineStyle(2);
	line1->SetLineColor(1);
	line1->SetLineStyle(2);
	line2->SetLineColor(1);
	line2->SetLineStyle(2);
	line0->Draw("SAME");
	line1->Draw("SAME");
	line2->Draw("SAME");

	latex->SetTextAlign(12);
	latex->DrawLatex(0.2, 0.8,"pp, #sqrt{#it{s}} = 13 TeV, #it{This thesis}");
	latex->DrawLatex(0.2, 0.73, Form("%s, %s",_Region,PidLatex(_PidIdx)));
	gSystem->Exec(Form("mkdir -p ./%s/plots_differential",dirOut));
	latex->SetTextAlign(22);

	if(isNch){
		latex->DrawLatex(0.5, 0.25,Form("#it{N}_{ch}^{TS} = %d, |#eta|<0.8",int(binNch-1)));
		cr->SaveAs(Form("./%s/plots_differential/%s_PID_%s_nch_%d.pdf",dirOut,_Region,_Pid,int(binNch-1)));
	}
	else if(isRT){
		latex->DrawLatex(0.5, 0.25,Form("#bf{%.2f < #it{R}_{T} < %.2f}, |#eta|<0.8",lowerBound,upperBound));
		cr->SaveAs(Form("./%s/plots_differential/%s_PID_%s_RT_%d.pdf",dirOut,_Region,_Pid,bin));
		cr->SaveAs(Form("./%s/plots_differential/%s_PID_%s_RT_%d.png",dirOut,_Region,_Pid,bin));
	}
	else{
		cr->SaveAs(Form("./%s/plots_differential/%s_PID_%s_Bin_pt_%d.pdf",dirOut,_Region,_Pid,bin));
		latex->DrawLatex(0.4, 0.3,Form("%.2f < #it{p}_{T} < %.2f (GeV/#it{c})",lowerBound,upperBound));
		latex->DrawLatex(0.4, 0.22,Form("|#eta|<0.8"));
	}

	ObjArray->AddLast(hUnToGen);

}

TPad* UnfoldNTclass::GetPad()
{

	TPad *pad = new TPad();
	pad->SetLeftMargin(0.15);
	pad->SetBottomMargin(0.15);
	pad->SetTicky(1);
	pad->SetTickx(1);
	pad->SetLogy(kFALSE);
	pad->SetGridy(kFALSE);

	return pad;

}

TH1F* UnfoldNTclass::TuneFrame(const double& min, const double& max, const char* titleX, const char* titleY)
{

	TH1F* hframe = new TH1F("frame","",(int)max-(int)min,min,max);

	hframe->GetXaxis()->SetTitle(Form("%s",titleX));
	hframe->GetYaxis()->SetTitle(Form("%s",titleY));
	hframe->GetYaxis()->SetRangeUser(0.61,1.59);

	return hframe;

}

void UnfoldNTclass::myLegendSetUp(TLegend *legend, float size)
{

	legend->SetTextFont(42);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->SetFillColor(0);
	legend->SetMargin(0.25);
	legend->SetTextSize(size);
	legend->SetEntrySeparation(0.5);

	return;

}

void UnfoldNTclass::TuneLatex(TLatex* latex)
{

	latex->SetNDC();
	latex->SetTextAlign(22);
	latex->SetTextAngle(0);
	latex->SetTextFont(42);
	latex->SetTextSize(0.06);

}

TH2F* UnfoldNTclass::SumCharges2D(const char* Eta, const TH2F* hPos ,const TH2F* hNeg)
{

	TH2F* h2D = (TH2F*)hPos->Clone(Form("hSumPosNeg_Eta_%s",Eta));
	h2D->Reset();

	for(int biny = 1; biny <= h2D->GetNbinsY(); biny++){
		for(int binx = 1; binx <= h2D->GetNbinsX(); binx++){

			h2D->SetBinContent(binx,biny,hPos->GetBinContent(binx,biny)+hNeg->GetBinContent(binx,biny));

		}
	}

	return h2D;
}

TH2F* UnfoldNTclass::SumDistPerEta2D(const TH2F* h02, const TH2F* h24 ,const TH2F* h46 ,const TH2F* h68)
{

	TH2F* h2D = new TH2F("hFullAzimuth",";#it{N}_{ch};#it{p}_{T}^{rec} (GeV/#it{c})",nNchBins,NchBins,NPTBINS[_PidIdx],XBINS[_PidIdx]);

	for(int biny = 1; biny <= h2D->GetNbinsY(); biny++){
		for(int binx = 1; binx <= h2D->GetNbinsX(); binx++){

			double y02 = h02->GetBinContent(binx,biny);
			double y24 = h24->GetBinContent(binx,biny);
			double y46 = h46->GetBinContent(binx,biny);
			double y68 = h68->GetBinContent(binx,biny);

			h2D->SetBinContent(biny,binx,y02+y24+y46+y68);

		}
	}

	return h2D;
}
//_________________________________________
void UnfoldNTclass::GetUnfoldingMatrix2D()
{

	for(int nt = 0; nt < _nt; nt++){
		for(int nm = 0; nm < _nm; nm++){

			_hSolutionNT->SetBinContent(nm+1,nt+1,_Mij(nt,nm));

		}
	}

}
//_________________________________________
void UnfoldNTclass::GetCovarianceMatrix2D()
{

	for(int nt = 0; nt < _nt; nt++){
		for(int nm = 0; nm < _nm; nm++){

			_hCovMatrixNT->SetBinContent(nm+1,nt+1,_dnCidnEj(nt,nm));

		}
	}

}
//_________________________________________
void UnfoldNTclass::GetUnfoldingMatrix3D(int bin_pt)
{

	for(int nt = 0; nt < _nt; nt++){
		for(int nm = 0; nm < _nm; nm++){

			//_hUnfMatrix->SetBinContent(nt+1,nm+1,bin_pt,_Mij(nt,nm));
			_hUnfMatrix->SetBinContent(nm+1,nt+1,bin_pt,_Mij(nt,nm));

		}
	}
}
//_________________________________________
void UnfoldNTclass::FillCovarianceMatrix(int bin_pt)
{

	for(int nt = 0; nt < _nt; nt++){
		for(int nm = 0; nm < _nm; nm++){

			_hCovMatrix->SetBinContent(nm+1,nt+1,bin_pt,_dnCidnEj(nt,nm));

		}
	}

}
//_________________________________________
void UnfoldNTclass::FillHistoVsRT(TH2F* hVsNch, TH2F* hVsRT)
{

	for(int binrt = 1; binrt <= nRTBins; binrt++){

		int lowNchBin = GetRTBin(binrt,kTRUE);
		int upNchBin  = GetRTBin(binrt,kFALSE);

		TH1D* hProj = (TH1D*)hVsNch->ProjectionY(Form("hProj_from_Bin_%d_to_%d",lowNchBin,upNchBin),lowNchBin,upNchBin,"e");

		for(int binpt = 1; binpt <= NPTBINS[_PidIdx]; binpt++){

			hVsRT->SetBinContent(binrt,binpt,hProj->GetBinContent(binpt));
			hVsRT->SetBinError(binrt,binpt,hProj->GetBinError(binpt));

		}
	}

}

//_____________________________________________________________________________

int UnfoldNTclass::GetRTBin(const int& binRt, bool isLowEdge)
{

	int binNch = -1;

	//! <NT> = 7. 43
	
	if( binRt == 0 ){
		if(isLowEdge) binNch = 1;
		else binNch = 50;
	}
	else if( binRt == 1 ){ //! From NT = 0 to NT = 6
		if(isLowEdge) binNch = 1;
		else binNch = 7;//else binNch = 4;
	}
	else if( binRt == 2 ){ //! From NT = 6 to NT = 11
		if(isLowEdge) binNch = 8;
		else binNch = 12;
	}
	else if( binRt == 3 ){//! From NT = 12 to NT = 18
		if(isLowEdge) binNch = 13;
		else binNch = 19;
	}
	else if( binRt == 4 ){//! From NT = 19 to NT = 36
		if(isLowEdge) binNch = 20;
		else binNch = 38;
	}
	else{
		if(isLowEdge) binNch = 39;
		else binNch = 50;
	}

	return binNch;

}

//_____________________________________________________________________________

int UnfoldNTclass::GetRTMinBin(const int& binRt, bool isLowEdge)
{

	int binNch = -1;

	//! <NT> = 2.559
	
	if( binRt == 0 ){
		if(isLowEdge) binNch = 1;
		else binNch = 50;
	}
	else if( binRt == 1 ){ //! From NT = 0 to NT = 2
		if(isLowEdge) binNch = 1;
		else binNch = 3;//else binNch = 2;
	}
	else if( binRt == 2 ){ //! From NT = 3 to NT = 3
		if(isLowEdge) binNch = 4;
		else binNch = 4;
	}
	else if( binRt == 3 ){//! From NT = 4 to NT = 6
		if(isLowEdge) binNch = 5;
		else binNch = 7;
	}
	else if( binRt == 4 ){//! From NT = 7 to NT = 12
		if(isLowEdge) binNch = 8;
		else binNch = 13;
	}
	else{
		if(isLowEdge) binNch = 14;
		else binNch = 50;
	}

	return binNch;

}

int UnfoldNTclass::GetRTMaxBin(const int& binRt, bool isLowEdge)
{

	int binNch = -1;

	//! <NT> = 4.871
	
	if( binRt == 0 ){
		if(isLowEdge) binNch = 1;
		else binNch = 50;
	}
	else if( binRt == 1 ){ //! From NT = 0 to NT = 4
		if(isLowEdge) binNch = 1;
		else binNch = 5;//else binNch = 3;
	}
	else if( binRt == 2 ){ //! From NT = 4 to NT = 7
		if(isLowEdge) binNch = 6;
		else binNch = 8;
	}
	else if( binRt == 3 ){//! From NT = 8 to NT = 12
		if(isLowEdge) binNch = 9;
		else binNch = 13;
	}
	else if( binRt == 4 ){//! From NT = 13 to NT = 24
		if(isLowEdge) binNch = 14;
		else binNch = 25;
	}
	else{
		if(isLowEdge) binNch = 26;
		else binNch = 50;
	}

	return binNch;

}

//_____________________________________________________________________________



void UnfoldNTclass::ExtrapolateRM(TH2F* hMultRespMatrix)
{

	TRandom* fRandom = new TRandom(0);

	TGraphErrors* grMean = new TGraphErrors();
	grMean->SetName("grMean");
	grMean->SetTitle(";#it{N}_{ch}^{TS};#LT #it{N}_{acc}^{TS} #GT");
	grMean->SetMarkerStyle(21);
	TGraphErrors* grSigma = new TGraphErrors();
	grSigma->SetName("grSigma");
	grSigma->SetTitle(";#it{N}_{ch}^{TS};#sigma");
	grSigma->SetMarkerStyle(21);

	cout << "EXTR got here 1 \n";

	int point = 0;
	for(int binnch = 1; binnch <= hMultRespMatrix->GetNbinsY(); binnch++){

		if(binnch > 30) continue;
		cout << "EXTR got here 11 \n";
		TH1D* hproj = (TH1D*)hMultRespMatrix->ProjectionX(Form("hproj_%d",binnch),binnch,binnch);
		TF1* fgaus = new TF1(Form("fGaus_%d",binnch),"gaus",0.0,100.0);
		cout << hproj << endl;
		cout << fgaus << endl;
		fgaus->SetParameters(hproj->Integral(),hproj->GetMean(),hproj->GetStdDev());
		cout << hproj->GetMean() << " " << hproj->GetStdDev() << "\n";
		
		hproj->Fit(fgaus,"L");//,"","", hproj->GetMean() - 3.0 * hproj->GetStdDev(), hproj->GetMean() + 3.0 * hproj->GetStdDev());
		cout << "EXTR got here 12 \n";
		double nch = hMultRespMatrix->GetYaxis()->GetBinCenter(binnch);
/*
		grMean->SetPoint(point,nch,hproj->GetMean());
		grMean->SetPointError(point,0.0,hproj->GetMeanError());
		grSigma->SetPoint(point,nch,hproj->GetStdDev());
		grSigma->SetPointError(point,0.0,hproj->GetStdDevError());
*/
		cout << "EXTR got here 13 \n";
		grMean->SetPoint(point,nch,fgaus->GetParameter(1));
		grMean->SetPointError(point,0.0,fgaus->GetParError(1));
		grSigma->SetPoint(point,nch,fgaus->GetParameter(2));
		grSigma->SetPointError(point,0.0,fgaus->GetParError(2));
		cout << "EXTR got here 14 \n";
		point++;
		cout << "EXTR got here 15 \n";
		delete fgaus;

	}

	cout << "EXTR got here 2 \n";

	TF1* fMean = new TF1("fMean","pol1",0.0,50.0);
	TF1* fSigma = new TF1("fSigma","([0]+[1]*x+[2]*x*x)/(1+[3]*x)",0.0,50.0);

	fMean->SetLineWidth(5);
	fSigma->SetLineWidth(5);
	fSigma->SetLineColor(4);

	cout << "EXTR got here 3 \n";

	grMean->Fit(fMean,"","",0.0,25.0 );
	grSigma->Fit(fSigma,"","",0.0,25.0);

	cout << "EXTR got here 4 \n";

	//! Extrapolate the Response Matrix for Nch > 20
	for(int binnch = 1; binnch <= hMultRespMatrix->GetNbinsY(); binnch++){

		double nch = hMultRespMatrix->GetYaxis()->GetBinCenter(binnch);

		if( nch < 15 ){

			for(int binacc = 1; binacc <= hMultRespMatrix->GetNbinsX(); binacc++)
				hMultRespMatrix->SetBinContent(binacc,binnch,hMultRespMatrix->GetBinContent(binacc,binnch));

		}else{
			for(int r = 1; r <= 1E6; r++){

				double random = fRandom->Gaus(fMean->Eval(nch),fSigma->Eval(nch));	
				hMultRespMatrix->Fill(random,nch);
			}
		}
	}

	cout << "EXTR got here 5 \n";

	ObjArray->AddLast(grMean);
	ObjArray->AddLast(grSigma);

	delete fMean;
	delete fSigma;

}
//_________________________
/*
   void UnfoldNTclass::TunePrior()
   {
   for (int i = 0; i < _nm; i++){

   _P0C(i) = 1.0;
   }


   }
   */
//_________________________
/*void UnfoldNTclass::LoadEventRecEfficiency()
  {

  const char* Path = "/Users/omar/Documents/DataAnalysis/pp_13TeV_spectra_RT/task_MC/Phi_test/Add_Fakes/efficiencies/MCEfficiencies_MB_EPOS_HybridTracksSpectra_V0s_Phi.root";
  TFile* fIn = new TFile(Form("%s",Path),"READ");
  TDirectory* dIn = (TDirectory*)fIn->Get("Nch_Efficiencies");

  _hEventRecEff = (TH1F*)dIn->Get("EventReconstructionEfficiency");
//_hEventRecEff = (TH1F*)dIn->Get("hMissedEvents");
_hEventRecEff->Clone("");
_hEventRecEff->SetDirectory(0);

delete fIn;
fIn = nullptr;

}
*/
