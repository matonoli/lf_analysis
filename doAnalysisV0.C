// V0 reconstruction macro using MyKit
// OliverM 2019 Lund

#include <iostream>
using namespace RooFit;

void doAnalysisV0(Int_t nEvents=100, const Char_t *inputFile="test.list", 
	const Char_t *outputFile="test.root") {

	// Loading ALICE libraries
	gROOT->LoadMacro("load_libraries.C");
	load_libraries();
	gROOT->LoadMacro("TransverseSpherocity/TransverseSpherocity.cxx+");

	// Loading and compiling custom MyKit libraries
	gROOT->LoadMacro("MyEvent.cxx+");
	gROOT->LoadMacro("MyTrack.cxx+");
	gROOT->LoadMacro("MyParticle.cxx+");
	gROOT->LoadMacro("MyV0.cxx+");
	gROOT->LoadMacro("MyAnalysis.cxx+");
	gROOT->LoadMacro("MyHandler.cxx++");
	gROOT->LoadMacro("MyAnalysisV0.cxx++");


	// Suppress RooFit spam
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	RooMsgService::instance().setSilentMode(true);

	// Set-up the analysis handler and load the input
	MyHandler* handler 	= new MyHandler();
	if (!handler->LoadInput(inputFile,"PIDTree")) {
		printf("ERROR: No files to analyse \n");
		return; }
	else {
		printf("Data successfully loaded \n");	}

	// Set-up analyses and link them to handler
	MyAnalysisV0* analysisV0 	= new MyAnalysisV0();
	handler->AddAnalysis(analysisV0);
	analysisV0->SetDirectory(gDirectory);
	analysisV0->SetOutputName(outputFile);
	printf("Following analyses will be performed: %s \n", analysisV0->GetName());

	// Initialise analyses
	handler->Init();

	// Start an event loop
	Int_t nEntries = handler->Chain()->GetEntries();
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




