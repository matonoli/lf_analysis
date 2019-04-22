// V0 reconstruction macro using MyKit
// OliverM 2019 Lund

#include <iostream>
using namespace RooFit;

void doAnalysisV0(Int_t nEvents=10, const Char_t *inputFile="test.list", 
	const Char_t *outputFile="test.root", const Char_t *MCinputFile="") {

	// Loading ALICE libraries
	gROOT->LoadMacro("load_libraries.C");
	load_libraries();
	gROOT->LoadMacro("TransverseSpherocity/TransverseSpherocity.cxx+");

	// Loading and compiling custom MyKit libraries
	gROOT->LoadMacro("MyKit/MyEvent.cxx+");
	gROOT->LoadMacro("MyKit/MyTrack.cxx+");
	gROOT->LoadMacro("MyKit/MyParticle.cxx+");
	gROOT->LoadMacro("MyKit/MyV0.cxx+");
	gROOT->LoadMacro("MyKit/MyAnalysis.cxx+");
	gROOT->LoadMacro("MyKit/MyHandler.cxx++");
	gROOT->LoadMacro("MyKit/Analyses/MyAnalysisV0.cxx++");
	gROOT->LoadMacro("MyKit/Analyses/MyAnalysisV0extract.cxx++");
	gROOT->LoadMacro("MyKit/Analyses/MyAnalysisV0correct.cxx++");

	// Suppress RooFit spam
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	RooMsgService::instance().setSilentMode(true);

	// Set-up the analysis handler and load the input
	MyHandler* handler 	= new MyHandler();
	handler->SetOutputName(outputFile);
	handler->SetDirectory(gDirectory);
	handler->SetROOT(gROOT);

	TString inputFileStr(inputFile);
	if (inputFileStr.Contains("hist")) {
		if (!handler->LoadInputHist(inputFile)) {
			printf("ERROR: No files to analyse \n");
			return; }
		else {
			printf("Data (histograms) successfully loaded \n");	}
	} else {
		if (!handler->LoadInputTree(inputFile,"PIDTree")) {
			printf("ERROR: No files to analyse \n");
			return; }
		else {
			printf("Data (PIDTree) successfully loaded \n");	}
	}

	// Set-up analyses and link them to handler
	MyAnalysisV0* analysisV0 	= new MyAnalysisV0();
	analysisV0->SetOutputName(outputFile);
	MyAnalysisV0extract* analysisV0extract 	= new MyAnalysisV0extract();
	MyAnalysisV0correct* analysisV0correct 	= new MyAnalysisV0correct();	//load mc file too

	handler->AddAnalysis(analysisV0);
	handler->AddAnalysis(analysisV0extract);
	handler->AddAnalysis(analysisV0correct);

	// Initialise analyses
	handler->Init();
	analysisV0correct->SetMCInputFile(MCinputFile);

	// Start an event loop
	Int_t nEntries = (handler->GetFlagHist()) ? 0 : handler->chain()->GetEntries();
	printf("Total entries: %i \n", nEntries);
	nEvents = (nEvents < nEntries) ? nEvents : nEntries;

	// Count how often should echo happen
	Int_t evCounter = 0; Int_t nDigits = nEvents;
	do { nDigits /= 10; evCounter++; } while (nDigits != 0);
	evCounter = (evCounter < 3) ? 1 : 0.5*TMath::Power(10,evCounter-2);
	
	for (Int_t iEv = 0; iEv < nEvents; iEv++)	{
		
		if (iEv%evCounter==0) {
			printf("Processing event %i, out of total %i events...\n", iEv, nEvents); }
		
		Int_t iret = handler->Make(iEv);
		if (iret) { 
			printf("Bad return code %i at event %i . Exiting... \n", iret, iEv); 
			break;	}
	}

	// Finish up
	printf("Finishing up...\n");
	handler->Finish();

	printf("Analysis finished. Exiting... \n");
}




