{
	gSystem->Load("libPWGLFspectra");
	TInterpreter* root = gInterpreter;
	
	root->LoadMacro("TransverseSpherocity.cxx+");
	root->LoadMacro("MyAnalysis.cxx+");
	root->LoadMacro("MyHandler.cxx+");
	root->LoadMacro("MyEvent.cxx+");
	root->LoadMacro("MyTrack.cxx+");
	root->LoadMacro("MyParticle.cxx+");
	root->LoadMacro("MyV0.cxx++");
	root->LoadMacro("MyAnalysisV0.cxx++");
	root->LoadMacro("MyAnalysisV0extract.cxx+");
	root->LoadMacro("MyAnalysisV0correct.cxx+");
	root->LoadMacro("MyAnalysisV0syst.cxx+");
	root->LoadMacro("MyAnalysisV0plot.cxx+");

}