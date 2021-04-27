// cuts based on f.fonda and p.kalinak
// except for proper lifetime (2d instead of 3d)
// not using upper bound of R
// using rapidity instead of eta
// pileup fast detector cut IS used
namespace cuts {
	// EVENT SELECTION

	// SPHEROCITY CUTS ARE pt=1; in order top/bottom 20%, 10%, 5%, 1%
	const Float_t EV_SPH_V0M_JETTY[4]			= {0.600, 0.529, 0.470, 0.357};
	const Float_t EV_SPH_V0M01_JETTY[4]			= {0.651, 0.589, 0.535, 0.433};
	const Float_t EV_SPH_NCharged_JETTY[4]		= {0.625, 0.561, 0.508, 0.408};
	const Float_t EV_SPH_NCharged01_JETTY[4]	= {0.680, 0.624, 0.577, 0.487};
	const Float_t EV_SPH_V0M_ISO[4]				= {0.823, 0.864, 0.891, 0.927};
	const Float_t EV_SPH_V0M01_ISO[4]			= {0.846, 0.882, 0.905, 0.936};
	const Float_t EV_SPH_NCharged_ISO[4]		= {0.833, 0.871, 0.896, 0.930};
	const Float_t EV_SPH_NCharged01_ISO[4]		= {0.859, 0.892, 0.913, 0.942};	

	const Float_t EV_MHM[2]			= {26., 999.};
	const Float_t EV_MHM01[2]		= {47., 999.};

	// V0 SELECTION
	const Float_t V0_ETA[2] 	= {-0.8, 0.8};
	const Float_t V0_Y[2]		= {-0.5, 0.5};
	const Float_t V0_PT[2] 		= {0.10, 15.};
	//const Float_t V0_DCADD		= 0.5;
	const Float_t V0_DCADD		= 1.0;
	const Float_t V0_CPA		= 0.97;
	const Float_t V0_R[2]		= {0.5, 9999.};
	const Int_t V0_FASTSIGNAL	= 1;

	// V0 DAUGHTER SELECTION
	const Float_t V0_D_ETA[2]	= {-0.8, 0.8};
	const Float_t V0_D_DCAPVXY 	= 0.06;
	const Bool_t V0_D_GOODTRACK	= 1;
	const Float_t V0_D_NCR		= 70.;
	const Float_t V0_D_NRATIO		= 0.8;

	// K0s, L, Lbar
	const Float_t K0S_AP 				= 0.2;
	const Float_t K0S_D_NSIGTPC[2] 		= {-5., 5.};
	const Float_t K0S_D_ETA[2]			= {-0.8, 0.8};
	//const Float_t K0S_D_DCAPVXY			= 0.05;
	const Float_t K0S_COMP_M			= 0.005;
	const Float_t K0S_TAU				= 20.;
	//const Float_t K0S_NSIGTOF[2] 		= {-3., 3.};
	const Float_t L_NSIGTPC[2] 			= {-5., 5.};
	//const Float_t L_NSIGTOF[2] 		= {-3., 3.}; 
	const Float_t LBAR_NSIGTPC[2] 		= {-5., 5.};
	//const Float_t LBAR_NSIGTOF[2] 	= {-3., 3.};
	const Float_t L_PT[2] 				= {0.4, 15.};
	const Float_t L_COMP_M				= 0.010;
	const Float_t L_TAU					= 30.;
	const Float_t L_CPA					= 0.995;

	// PRIMARY TRACK SELECTION
	const Float_t TR_PRIMARY_PAR[3]		= {0.0182, 0.0350, 1.01}; // dont use this for spherocity!


	const Float_t K0S_PARMU[4]	= {-2.302e-3, 4.568e-3, -1.689e-3, 8.383e-4 };
	const Float_t K0S_PARSIG[3]	= {3.274e-3, 3.539e-4, 1.193e-4};
	const Float_t L_PARMU[5]	= {-6.323e-4, 9.629e-4, -2.809e-4, 3.449e-4, -5.345e-05 };
	const Float_t L_PARSIG[3]	= {1.304e-3, 1.820e-4, 5.259e-4};
	const Float_t V0_COMP_NSIG = 4.;
}