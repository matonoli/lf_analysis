// cuts based on f.fonda and p.kalinak
// except for proper lifetime (2d instead of 3d)
// not using upper bound of R
// using rapidity instead of eta
// pileup fast detector cut IS used
namespace cuts {
	// EVENT SELECTION

	//const Float_t EV_SPH_JETTY 	= 0.446;		//mc results
	//const Float_t EV_SPH_ISO 	= 0.747;

	//const Float_t EV_SPH_JETTY[2] 	= {0.462, 0.484};	//data, [0] is v0m. [1] is ncharged 
	//const Float_t EV_SPH_ISO[2] 	= {0.755, 0.765};
	const Float_t EV_SPH_JETTY[4] 	= {0.600, 0.625, 0.651, 0.680};	//pt=1, V0M, CL1, V0M01, CL101
	const Float_t EV_SPH_ISO[4] 	= {0.823, 0.833, 0.846, 0.859};
	const Float_t EV_MHM[2]			= {26., 300.};
	const Float_t EV_MHM01[2]		= {47., 300.};

	// V0 SELECTION
	const Float_t V0_ETA[2] 	= {-0.8, 0.8};
	const Float_t V0_Y[2]		= {-0.5, 0.5};
	const Float_t V0_PT[2] 		= {0.10, 15.};
	//const Float_t V0_DCADD		= 0.5;
	const Float_t V0_DCADD		= 1.0;
	const Float_t V0_CPA		= 0.97;
	const Float_t V0_R[2]		= {0.5, 9999.};
	const Bool_t V0_FASTSIGNAL	= 1;

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
	const Float_t K0S_D_DCAPVXY			= 0.05;
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
	const Float_t L_CPA					= 0.997;

	// PRIMARY TRACK SELECTION
	const Float_t TR_PRIMARY_PAR[3]		= {0.0182, 0.0350, 1.01}; // dont use this for spherocity!

}