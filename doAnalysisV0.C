// V0 reconstruction macro using MyKit
// OliverM 2019 Lund

#include <iostream>
#include <fstream>
#include <string>
using namespace RooFit;

#if !defined (__CINT__) || defined (__CLING__)
#include "../load_libraries.C"

#include "AliAnalysisAlien.h"
#include "AliAnalysisGrid.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskMyTask.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TChain.h"
#include "MyKit/MyHandler.h"
#include "MyKit/Analyses/MyAnalysisV0.h"
#include "MyKit/Analyses/MyAnalysisV0extract.h"
#include "MyKit/Analyses/MyAnalysisV0correct.h"
#include "MyKit/Analyses/MyAnalysisV0plot.h"

	//using rootcl = TInterpreter;

#else
	//class TROOT;
	//using rootcl = TROOT;
#endif

// THIS MACRO SHOULD BE RUN FROM A WORKING DIRECTORY INSIDE AN ANALYSIS FOLDER

void doAnalysisV0(Long_t nEvents=10, const Char_t *flags = "0", const Char_t *inputFile="test.list", 
	const Char_t *outputFile="test.root", const Char_t *MCinputFile="") {

	// flags:	0...basic analysis
	//			x...signal extraction
	//			c...post-processing corrections
	//			p...plotting
	//			E...run on ESD files instead of local trees
	//			L...test analysis on an ESD file locally, inputfile must be AliESDs.root file
	//			G...test analysis on ESDs on the GRID, inputfile must be a run list file
	//			M...submit merging jobs on the GRID
	//			D...download merged file from the GRID

	#if !defined (__CINT__) || defined (__CLING__)
	std::cout << "---Launching this analysis using ROOT6 \n";
	TInterpreter* root = gInterpreter;
	#else
	std::cout << "---Launching this analysis using ROOT5 \n";
	TROOT* root = gROOT;
	#endif

	// CHECK INPUT FLAGS
	TString fl(flags);	enum { notset, trees, esds };	Int_t mode;
	if (fl.Contains("E")) mode = esds;
	else mode = trees;
	if (mode != esds && ( fl.Contains("L") || fl.Contains("G") || 
						fl.Contains("M") || fl.Contains("D") ))		{
		std::cout << "!!!Invalid flags, aborting. \n";
		return;	}

	// SPECIFY WHAT CUTS TO USE
	TString cutStr("cuts03.h");
	if (!gSystem->AccessPathName(Form("../%s",cutStr.Data())))
		std::cout << "---Using cuts specified in file " << cutStr.Data() << "\n";
	else {
		std::cout << "!!!Cannot load cuts from file " << cutStr.Data() << " , aborting. \n";
		return;	}

	// LOADING ALICE LIBRARIES
	#if defined (__CINT__) || !defined (__CLING__)
	gROOT->LoadMacro("../load_libraries.C");
	#endif
	load_libraries();
	gSystem->Load("libCore.so");
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libpythia6.so");
	gSystem->Load("libpythia6_4_21.so");
	gSystem->Load("libpythia6_4_25.so");
	gSystem->Load("libpythia6_4_28.so");
	gSystem->Load("libAliPythia6");			// order here matters

	// LOAD INCLUDES FOR HEADERS BEFORE COMPILING
	root->ProcessLine(".include $ROOTSYS/include");
	root->ProcessLine(".include $ALICE_ROOT/include");

	// COPY FILES
	gSystem->Exec("cp ../TransverseSpherocity/TransverseSpherocity.cxx .; cp ../TransverseSpherocity/TransverseSpherocity.h .");	// extra classes
	gSystem->Exec("cp ../MyKit/MyAnalysis.cxx .; cp ../MyKit/MyAnalysis.h ."); 		// MyKit
	gSystem->Exec("cp ../MyKit/MyHandler.cxx .; cp ../MyKit/MyHandler.h .");
	gSystem->Exec("cp ../MyKit/MyEvent.cxx .; cp ../MyKit/MyEvent.h .");
	gSystem->Exec("cp ../MyKit/MyTrack.cxx .; cp ../MyKit/MyTrack.h .");
	gSystem->Exec("cp ../MyKit/MyParticle.cxx .; cp ../MyKit/MyParticle.h .");
	gSystem->Exec("cp ../MyKit/MyV0.cxx .; cp ../MyKit/MyV0.h .");
	gSystem->Exec("cp ../MyKit/Analyses/MyAnalysisV0.cxx .; cp ../MyKit/Analyses/MyAnalysisV0.h .");	// MyKit analyses
	gSystem->Exec("cp ../MyKit/Analyses/MyAnalysisV0extract.cxx .; cp ../MyKit/Analyses/MyAnalysisV0extract.h .");
	gSystem->Exec("cp ../MyKit/Analyses/MyAnalysisV0correct.cxx .; cp ../MyKit/Analyses/MyAnalysisV0correct.h .");
	gSystem->Exec("cp ../MyKit/Analyses/MyAnalysisV0plot.cxx .; cp ../MyKit/Analyses/MyAnalysisV0plot.h .");
	gSystem->Exec("cp ../AliAnalysisTaskMyTask.cxx .; cp ../AliAnalysisTaskMyTask.h .");	// alice task
	gSystem->Exec("cp ../AddMyTask.C .");
	gSystem->Exec(Form("cp ../%s .", cutStr.Data()));	// cut file

	// DEFINE INSTRUCTIONS FOR COMPILING
	gSystem->Exec("echo '#ifndef COMPINSTRUCTIONS_H\n#define COMPINSTRUCTIONS_H' > compInstructions.h");
	gSystem->Exec(Form("echo '#define INPUTFORMAT %i' >> compInstructions.h", mode));
	gSystem->Exec(Form("echo '// %i : local aurora trees' >> compInstructions.h", trees));
	gSystem->Exec(Form("echo '// %i : esd files \n' >> compInstructions.h", esds));
	gSystem->Exec(Form("echo '#include \"%s\"' >> compInstructions.h", cutStr.Data()));
	gSystem->Exec("echo '#endif' >> compInstructions.h");

	// COMPILE CLASSES AND TASKS
	root->LoadMacro("TransverseSpherocity.cxx+");
	root->LoadMacro("MyAnalysis.cxx+");
	root->LoadMacro("MyHandler.cxx+");
	root->LoadMacro("MyEvent.cxx+");
	root->LoadMacro("MyTrack.cxx+");
	root->LoadMacro("MyParticle.cxx+");
	root->LoadMacro("MyV0.cxx++");
	if (fl.Contains("0")) root->LoadMacro("MyAnalysisV0.cxx++");
	if (fl.Contains("x")) root->LoadMacro("MyAnalysisV0extract.cxx+");
	if (fl.Contains("c")) root->LoadMacro("MyAnalysisV0correct.cxx+");
	if (fl.Contains("p")) root->LoadMacro("MyAnalysisV0plot.cxx+");
	if (mode == esds) {
		root->LoadMacro("AliAnalysisTaskMyTask.cxx+");
		root->LoadMacro("AddMyTask.C");
		root->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
		root->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
		root->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"); }

	// SILENCE ROOFIT SPAM
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	RooMsgService::instance().setSilentMode(true);

	// SETUP ANALYSIS ON TREES
	if (mode == trees) {
		// SETUP ANALYSIS HANDLER
		MyHandler* handler = new MyHandler();
		handler->SetOutputName(outputFile);
		handler->SetDirectory(gDirectory);
		handler->SetROOT(gROOT);
		handler->RebinPt(false);

		// CHECK INPUTFILE
		TString inputFileStr(inputFile);
		if (inputFileStr.Contains("hist")) {
			if (!handler->LoadInputHist(inputFile)) {					// this also determines the MC flag for analysis
				std::cout << "!!!ERROR: No hist file to analyse \n";
				return; }
			else {
				std::cout << "---Data (histograms) successfully loaded \n";	}
		} else {
			if (!handler->LoadInputTree(inputFile,"PIDTree")) {			// this also determines the MC flag for analysis
				std::cout << "!!!ERROR: No files to analyse \n";
				return; }
			else {
				std::cout << "---Data (PIDTree) successfully loaded \n";	}
		}

		// ADD SUB-ANALYSES
		if (fl.Contains("0")) {
			MyAnalysisV0* analysisV0					= new MyAnalysisV0();
			handler->AddAnalysis(analysisV0); }
		if (fl.Contains("x")) {
			MyAnalysisV0extract* analysisV0extract		= new MyAnalysisV0extract();
			analysisV0extract->SetMCInputFile(MCinputFile);
			handler->AddAnalysis(analysisV0extract);	}
		if (fl.Contains("c")) {
			MyAnalysisV0correct* analysisV0correct		= new MyAnalysisV0correct();
			analysisV0correct->SetMCInputFile(MCinputFile);
			handler->AddAnalysis(analysisV0correct); }
		if (fl.Contains("p")) {
			MyAnalysisV0plot* analysisV0plot			= new MyAnalysisV0plot();
			handler->AddAnalysis(analysisV0plot); }

		// INITIATE ANALYSES
		handler->Init();

		// PREPARING AN EVENT LOOP
		Long_t nEntries = (handler->GetFlagHist()) ? 0 : handler->chain()->GetEntries();
		printf("---Total entries: %i \n", nEntries);
		nEvents = (nEvents < nEntries) ? nEvents : nEntries;

		// COUNT HOW OFTEN SHOULD ECHO HAPPEN
		Int_t evCounter = 0; Int_t nDigits = nEvents;
		do { nDigits /= 10; evCounter++; } while (nDigits != 0);
		evCounter = (evCounter < 3) ? 1 : 0.5*TMath::Power(10,evCounter-2);
		
		// EVENT LOOP
		for (Long_t iEv = 0; iEv < nEvents; iEv++)	{
			
			if (iEv%evCounter==0) {
				printf("--Processing event %i, out of total %i events...\n", iEv, nEvents); }
			
			Int_t iret = handler->Make(iEv);
			if (iret) { 
				printf("!!!Bad return code %i at event %i . Exiting... \n", iret, iEv); 
				break;	}
		}

		// FINISH ANALYSES
		printf("---Finishing up...\n");
		handler->Finish();

		printf("---Analysis finished. Exiting... \n");
		root->ProcessLine("new TBrowser");
		

		return;	}
	else if (mode == esds) {

		// CREATE ANALYSIS MANAGER
		AliAnalysisManager* mgr		= new AliAnalysisManager("AnalysisTaskV0");
		AliESDInputHandler* esdH	= new AliESDInputHandler();
		mgr->SetInputEventHandler(esdH);

		// DETERMINE MC FLAG
		TString inputFileStr(inputFile);
		Bool_t mcflag = inputFileStr.Contains("MC");
		if (mcflag) {
			std::cout << "---Running a MC analysis on ESDs \n";
			AliMCEventHandler* mcH	= new AliMCEventHandler();
			mgr->SetMCtruthEventHandler(mcH);
		}

		// SETUP AND ADD TASKS
		Bool_t enablePileupCuts = true; // true for pp, false for PbPb
		#if !defined (__CINT__) || defined (__CLING__)
		// PHYSICS SELECTION
		AliAnalysisTask *physSelTask	= reinterpret_cast<AliAnalysisTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d,%d)",mcflag,enablePileupCuts)));
		mgr->AddTask(physSelTask);
		// MULTIPLICITY SELECTION
		AliAnalysisTask *multSelTask	= reinterpret_cast<AliAnalysisTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
		mgr->AddTask(multSelTask);
		// PID RESPONSE
		AliAnalysisTask *pidResTask		= reinterpret_cast<AliAnalysisTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%d)",mcflag)));
		mgr->AddTask(pidResTask);
		// ANALYSIS TASK
		AliAnalysisTaskMyTask *task		= reinterpret_cast<AliAnalysisTaskMyTask*>(gInterpreter->ExecuteMacro(Form("AddMyTask.C(\"myTaskV0\",%i)",mcflag)));
		#else
		// PHYSICS SELECTION
		AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(mcflag,enablePileupCuts);
		if (mcflag) physSelTask->GetPhysicsSelection()->SetAnalyzeMC();
		// MULTIPLICITY SELECTION
		AliMultSelectionTask* multSelTask = AddTaskMultSelection();
		// PID RESPONSE
		AliAnalysisTaskPIDResponse* pidResTask = AddTaskPIDResponse();
		// ANALYSIS TASK
		AliAnalysisTaskMyTask* task	= AddMyTask("myTaskV0", mcflag);
		#endif

		// INITIATE TASKS
		if (!mgr->InitAnalysis()) {
			std::cout << "!!!Error initiating tasks in the manager, aborting. \n";
			return;		} 
		mgr->SetDebugLevel(2);
		mgr->PrintStatus();
		mgr->SetUseProgressBar(1,25);

		if (fl.Contains("L")) {
			// RUNNING ANALYSIS LOCALLY
			TChain* chain = new TChain("esdTree");
			chain->Add(inputFile);
			std::cout << "---Running analysis locally on a file " << inputFileStr.Data() << "\n";
			mgr->StartAnalysis("local",chain);
		} else {
			// SETUP ANALYSIS ON THE GRID (TEST OR REAL)
			AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
			alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/PWGLF/SPECTRA/PiKaPr/TPCTOFfits/ -I$PWD/Scripts/");

			// INCLUDE ALL HEADERS AND SOURCE FILES
			TString strH("compInstructions.h libpythia6_4_28.so ");
			strH += cutStr.Data(); strH += " ";
			strH += "TransverseSpherocity.h ";
			strH += "MyAnalysis.h MyHandler.h MyEvent.h MyTrack.h MyParticle.h MyV0.h ";
			strH += "MyAnalysisV0.h AliAnalysisTaskMyTask.h ";
			strH += "TransverseSpherocity.cxx ";
			strH += "MyAnalysis.cxx MyHandler.cxx MyEvent.cxx MyTrack.cxx MyParticle.cxx MyV0.cxx ";
			strH += "MyAnalysisV0.cxx AliAnalysisTaskMyTask.cxx ";
			TString strS("");
			strS += "TransverseSpherocity.cxx ";
			strS += "MyAnalysis.cxx MyHandler.cxx MyEvent.cxx MyTrack.cxx MyParticle.cxx MyV0.cxx ";
			strS += "MyAnalysisV0.cxx AliAnalysisTaskMyTask.cxx ";
			alienHandler->SetAdditionalLibs(strH.Data());
			alienHandler->SetAnalysisSource(strS.Data());

			// SELECT ALIPHYSICS VERSION, OTHER PACKAGES ARE LOADED AUTOMATICALLY
			alienHandler->SetAliPhysicsVersion("vAN-20200304-1");
			// SET THE ALIEN API VERSION
			alienHandler->SetAPIVersion("V1.1x");

			// SPECIFY THE INPUT DATA
			if (!inputFileStr.Contains(".list"))	{
				std::cout << "!!!A wrong runlist selected for run on the grid, aborting. \n";
				return;	}
			ifstream inputStream(inputFileStr.Data());
			if (!inputStream)	{
				cout << "!!!ERROR: Cannot open list file " << inputFileStr.Data() << endl;
				return;	}
			std::vector<int> runNumbers;
			TString yearStr = inputFileStr(inputFileStr.First("20"),4);
			TString periStr = inputFileStr(inputFileStr.First("LHC"), (mcflag ? 7 : 6) );
			TString passStr = inputFileStr(inputFileStr.First("pa"),5);
			std::string number_as_string;
    		while (std::getline(inputStream, number_as_string, ', '))	{
    			TString tstr(number_as_string);
				runNumbers.push_back(tstr.Atoi());
    		}
    		std::cout << "---Specified " << runNumbers.size() << " runs in total belonging to period "
    			<<  periStr.Data() << " from " << yearStr.Data() << " with " << passStr.Data() << "\n";

			if (!mcflag) {
				TString gdd = "/alice/data/"; gdd += yearStr; gdd += "/"; gdd += periStr;
				alienHandler->SetGridDataDir(gdd.Data());
				TString dpstr = "*"; dpstr += passStr.Data(); dpstr += "/*ESDs.root"; 
				alienHandler->SetDataPattern(dpstr.Data());
				alienHandler->SetRunPrefix("000");
				for (int iR = 0; iR < runNumbers.size(); iR++) alienHandler->AddRunNumber(runNumbers[iR]);
			} else {
				TString gdd = "/alice/sim/"; gdd += yearStr; gdd += "/"; gdd += periStr;
				alienHandler->SetGridDataDir(gdd.Data());
				alienHandler->SetDataPattern("*ESDs.root");
				alienHandler->SetRunPrefix("");
				for (int iR = 0; iR < runNumbers.size(); iR++) alienHandler->AddRunNumber(runNumbers[iR]);
			}

			// SET NUMBER OF FILES PER SUBJOB
			alienHandler->SetSplitMaxInputFileNumber(120);
			// SET EXECUTABLE AND JDL NAME
			alienHandler->SetExecutable("myTask.sh");
			alienHandler->SetJDLName("myTask.jdl");
			// ESTIMATE HOW MANY SECONDS WILL THE JOB TAKE
			alienHandler->SetTTL(10000);
			// ADDITIONAL INFO
			alienHandler->SetCheckCopy(kTRUE);
			alienHandler->SetOutputToRunNo(kTRUE);
			alienHandler->SetKeepLogs(kTRUE);
			// SET MERGING INFO
			alienHandler->SetMaxMergeStages(1);
			if (fl.Contains("D")) alienHandler->SetMergeViaJDL(kFALSE);
			else alienHandler->SetMergeViaJDL(kTRUE);
			// SPECIFY WORKING AND OUTPUT DIRECTORIES
			alienHandler->SetGridWorkingDir("myWorkingDir");
			alienHandler->SetGridOutputDir("myOutputDir");
			// CONNECT THE ALIEN PLUGIN TO THE ANALYSIS MANAGER
			mgr->SetGridHandler(alienHandler);
			if (fl.Contains("G")) {
				// HOW MANY FILES FOR TEST RUN
				alienHandler->SetNtestFiles(2);
				// LAUNCH ANALYSIS
				alienHandler->SetRunMode("test");
				mgr->StartAnalysis("grid");
			} else {
				if (fl.Contains("M") || fl.Contains("D")) alienHandler->SetRunMode("terminate");
				else alienHandler->SetRunMode("full");
				mgr->StartAnalysis("grid");
			}


		}		

	}

	gSystem->Unlink("tmp.root");
}




