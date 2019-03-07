namespace cuts {
	// EVENT SELECTION

	// V0 SELECTION
	const Float_t V0_ETA[2] 	= {-0.8, 0.8};
	const Float_t V0_PT[2] 		= {1., 14.};
	const Float_t V0_DCADD		= 0.5;
	const Float_t V0_CPA		= 0.997;
	const Float_t V0_R[2]		= {5., 100.};

	// V0 DAUGHTER SELECTION
	const Float_t V0_D_ETA[2]	= {-0.8, 0.8};
	const Float_t V0_D_DCAPVXY 	= 0.05;

	// K0s, L, Lbar
	const Float_t K0S_AP 				= 0.2;
	const Float_t K0S_D_NSIGTPC[2] 		= {-3., 3.};
	const Float_t K0S_D_ETA[2]			= {-0.8, 0.8};
	const Float_t K0S_D_DCAPVXY			= 0.05;
	//const Float_t K0S_NSIGTOF[2] 		= {-3., 3.};
	const Float_t L_NSIGTPC[2] 			= {-3., 3.};
	//const Float_t L_NSIGTOF[2] 		= {-3., 3.}; 
	const Float_t LBAR_NSIGTPC[2] 		= {-3., 3.};
	//const Float_t LBAR_NSIGTOF[2] 	= {-3., 3.}; 
}