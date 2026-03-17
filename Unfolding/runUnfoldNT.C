{

	gROOT->LoadMacro("UnfoldNTclassModified.cxx++");
	gROOT->LoadMacro("UnfoldNT.cxx++");
	static string Charge[3] = {"Charged","Pos","Neg"};
	static string Detector[2] = {"TPC","TOF"};

	static const bool isMC = 1;
	static string Trackcuts = "Hybrid";

	if(!isMC){
		gSystem->Exec("rm -f ./results/2D_newClass_data.root");
		for(int ch = 0; ch < 1; ++ch){
			//gROOT->ProcessLine(Form("UnfoldNT\(\"../Data_MC/data/one_hybrid/AnalysisResults.root\",\"../Data_MC/mc/one_hybrid/AnalysisResults.root\",\"results\",%d,\"%s\",\"TPC\",\"%s\")",isMC,Trackcuts.data(),Charge[ch].c_str()));
			gROOT->ProcessLine(Form("UnfoldNT\(\"../HybridStudy/workdir2016kHybridStudy/beast_220510_2016_LHC16k00.root\",\"../rootOutputs/beastMC_220510_20167_hist.root\",\"results\",%d,\"%s\",\"TPC\",\"%s\")",isMC,Trackcuts.data(),Charge[ch].c_str()));
		}
	}else{

		gSystem->Exec("rm -f ./results/2D_newClass_mc.root");
		//gROOT->ProcessLine(Form("UnfoldNT\(\"../rootOutputs/beastMC_220510_20167_hist.root\",\"../rootOutputs/beastMC_220510_20167_hist.root\",\"results\",%d,\"%s\",\"TPC\",\"Charged\")",isMC,Trackcuts.data()));
		//gROOT->ProcessLine(Form("UnfoldNT\(\"../rootOutputs/beastMC_220719_1617_hist.root\",\"../rootOutputs/beastMC_220719_1617_hist.root\",\"results\",%d,\"%s\",\"TPC\",\"Charged\")",isMC,Trackcuts.data()));
		gROOT->ProcessLine(Form("UnfoldNT\(\"out_mc.root\",\"out_mc.root\",\"results\",%d,\"%s\",\"TPC\",\"Charged\")",isMC,Trackcuts.data()));
		//gROOT->ProcessLine(Form("UnfoldNT\(\"unfolding_mc.root\",\"unfolding_mc.root\",\"results\",%d,\"%s\",\"TPC\",\"Charged\")",isMC,Trackcuts.data()));
//		gSystem->Exec("rm -f ./results/extraMCclosure.root");
//		gROOT->ProcessLine("MCclosure\(\"Transverse\")");

	}
}
