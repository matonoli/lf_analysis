/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include <iostream>

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TDirectoryFile.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
//#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h" //test
#include "AliMCEventHandler.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"


#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliAnalysisTaskMyTask.h"

#include "AliMultSelection.h"

#include "MyHandler.h"
#include "MyAnalysis.h"
#include "MyAnalysisV0.h"
#include "MyV0.h"
#include "TLorentzVector.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
//#include "MyKit/Analyses/MyAnalysisV0test.h"

class AliAnalysisTaskMyTask;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMyTask) // classimp: necessary for root

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask() : AliAnalysisTaskSE(), 
	fESD(0), fOutputList(0), fHistPt(0), fAnalysisMC(0), handler(0)
{
	// default constructor, don't allocate memory here!
	// this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) : AliAnalysisTaskSE(name),
	fESD(0), fOutputList(0), fHistPt(0), handler(0)
{
	// constructor
	DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
										// this chain is created by the analysis manager, so no need to worry about it, 
										// it does its work automatically
	DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
										// you can add more output objects by calling DefineOutput(2, classname::Class())
										// if you add more output objects, make sure to call PostData for all of them, and to
										// make changes to your AddTask macro!
	
	//DefineOutput(2, TDirectoryFile::Class());

	//DefineOutput(2, TList::Class());

	


	   

}

//_____________________________________________________________________________
void AliAnalysisTaskMyTask::DefineOutputs() {
	for (Int_t iAna = 0; iAna < handler->getNAnalyses(); iAna++) {
		DefineOutput(2+iAna, TList::Class());
	}
}

//_____________________________________________________________________________
AliAnalysisTaskMyTask::~AliAnalysisTaskMyTask()
{
	// destructor
	if(fOutputList) {
		delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserCreateOutputObjects()
{
	// create output objects
	//
	// this function is called ONCE at the start of your analysis (RUNTIME)
	// here you ceate the histograms that you want to use 
	//
	// the histograms are in this case added to a tlist, this list is in the end saved
	// to an output file
	//

	//OpenFile(1);

	/*AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if(man){
		fInputHandler = (AliESDInputHandler*)(man->GetInputEventHandler());
		
	}*/

	//crashes here so needs debugging and possibly move constructor
	// SETTING UP HANDLER FOR CUSTOM ANALYSIS KIT
	handler = new MyHandler();
	handler->SetOutputName("tmp.root");
	handler->SetFlagMC(fAnalysisMC);

	MyAnalysisV0* analysisV0                = new MyAnalysisV0();
	handler->AddAnalysis(analysisV0); 
	//task->DefineOutputs(); 

	//handler->SetDirectory(gDirectory);
	//handler->SetROOT(gROOT);    
	//handler->LoadInputHist("test.root");
	handler->SetupTrackCuts();
	handler->Init();

	fOutputList = new TList();          // this is a list which will contain all of your histograms
										// at the end of the analysis, the contents of this list are written
										// to the output file
	fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
										// if requested (dont worry about this now)

	// example of a histogram
	fHistPt = new TH1F("fHistPt", "fHistPt", 10000, -5000, 5000);       // create your histogra
	fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
										// your histogram in the output file, add it to the list!
	
	//PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
										// fOutputList object. the manager will in the end take care of writing your output to file
										// so it needs to know what's in the output

	cout << "bla \n";
	//handler->analysis(0)->dirFile()->GetList()->ls();

	for (Int_t iAna = 0; iAna < handler->getNAnalyses(); iAna++) {
		handler->analysis(iAna)->dirFile()->GetList()->SetOwner(kTRUE);
		PostData(1+iAna, handler->analysis(iAna)->dirFile()->GetList());
	}
	// FOR EACH ANALYSIS POSTDATA MDIRECTORY    
	//PostData(2, handler->analysis(0)->directory());
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
	// user exec
	// this function is called once for each event
	// the manager will take care of reading the events from file, and with the static function InputEvent() you 
	// have access to the current event. 
	// once you return from the UserExec function, the manager will retrieve the next event from the chain
	
	Long64_t iEv = Entry();
	fESD = dynamic_cast<AliESDEvent*>(InputEvent());    // get an event (called fAOD) from the input file
														// there's another event format (ESD) which works in a similar wya
														// but is more cpu/memory unfriendly. for now, we'll stick with aod's
	if(!fESD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event


	handler->SetEvent(fESD);
	if (fAnalysisMC) {
		AliMCEvent* mcevent = dynamic_cast<AliMCEvent*>(MCEvent());
		if (!mcevent) {
				Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}// else cout << "yay \n";

		AliStack* mcstack = mcevent->Stack();
		if (!mcstack) {
				Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}// else cout << "yay \n";

		handler->SetMCStack(mcstack);
	}

	//cout << "tracks " << handler->event()->GetNumberOfTracks() << endl;
	handler->Make(iEv); // HANDLER NEEDS NEW METHOD

	/*for (int iv = 0; iv < fESD->GetNumberOfV0s(); ++iv)
	{
		AliESDv0* v0 = static_cast<AliESDv0*>(fESD->GetV0(iv));
		if (!v0) continue;

		// pdg code
		//pos track
		//pos track label
		int posLabel = TMath::Abs(handler->track(TMath::Abs(v0->GetPindex()))->GetLabel());
		int negLabel = TMath::Abs(handler->track(TMath::Abs(v0->GetNindex()))->GetLabel());
		int mLabel = handler->mcstack()->Particle(TMath::Abs(posLabel))->GetFirstMother();
		int mLabelN = handler->mcstack()->Particle(TMath::Abs(negLabel))->GetFirstMother();

		if ((handler->mcstack()->Particle(mLabel)->GetPdgCode()==310 || handler->mcstack()->Particle(mLabelN)->GetPdgCode()==310) && !v0->GetOnFlyStatus()) {
			if (handler->track((TMath::Abs(v0->GetPindex())))->GetSign() < 0 ) printf(" SIGN IS WRONG \n");
			
			//hand calculate
			TLorentzVector a;
			a.SetPxPyPzE(handler->track((TMath::Abs(v0->GetPindex())))->Px(),
				handler->track((TMath::Abs(v0->GetPindex())))->Py(),
				handler->track((TMath::Abs(v0->GetPindex())))->Pz(),
				handler->track((TMath::Abs(v0->GetPindex())))->E());

			TLorentzVector b;
			b.SetPxPyPzE(handler->track((TMath::Abs(v0->GetNindex())))->Px(),
				handler->track((TMath::Abs(v0->GetNindex())))->Py(),
				handler->track((TMath::Abs(v0->GetNindex())))->Pz(),
				handler->track((TMath::Abs(v0->GetNindex())))->E());

			TLorentzVector m;
			m = a+b;

			AliKFVertex PrimaryVtxKF(*fESD->GetPrimaryVertex());
			AliKFParticle::SetField(fESD->GetMagneticField());


			AliKFParticle* negKF;
			AliKFParticle* posKF;
			negKF = new AliKFParticle(*(v0->GetParamN()), -211);
			posKF = new AliKFParticle(*(v0->GetParamP()), 211);		
			AliKFParticle V0KF;
			V0KF += (*posKF);
			V0KF += (*negKF);
			V0KF.SetProductionVertex(PrimaryVtxKF);

			float massorig = v0->GetEffMass(2,2);
			v0->ChangeMassHypothesis(310);
			//printf("ev vz %4.2f v0 r %4.2f v0 pt %4.2f kf pt %4.2f real pt %4.2f masses 1: %4.4f 2: %4.4f 3: %4.4f 4: %4.4f -- DP %i %i %i %i  DN %i %i %i %i K0S %i \n", fESD->GetPrimaryVertexTracks()->GetZ(), 
			//	TMath::Sqrt(v0->Xv()*v0->Xv()+v0->Yv()*v0->Yv()), v0->Pt(), V0KF.GetPt(), handler->mcstack()->Particle(mLabel)->Pt(),
			//	massorig-0.497614, v0->GetEffMass()-0.497614, m.M()-0.497614, V0KF.GetMass()-0.497614,
			//	handler->mcstack()->Particle(posLabel)->GetPdgCode(), handler->mcstack()->IsPhysicalPrimary(posLabel), handler->mcstack()->IsSecondaryFromWeakDecay(posLabel), handler->mcstack()->IsSecondaryFromMaterial(posLabel),
			//	handler->mcstack()->Particle(negLabel)->GetPdgCode(), handler->mcstack()->IsPhysicalPrimary(negLabel), handler->mcstack()->IsSecondaryFromWeakDecay(negLabel), handler->mcstack()->IsSecondaryFromMaterial(negLabel),
			//	handler->mcstack()->IsPhysicalPrimary(mLabel)	);
			printf("ev vz %4.2f kf pt %4.2f real pt P %4.2f real pt N %4.2f posML %i %i  negML %i %i DP %i %i %i %i  DN %i %i %i %i \n",
			 fESD->GetPrimaryVertexTracks()->GetZ(), V0KF.GetPt(), handler->mcstack()->Particle(mLabel)->Pt(),
			 handler->mcstack()->Particle(mLabelN)->Pt(), mLabel, handler->mcstack()->Particle(mLabel)->GetPdgCode(),
			 mLabelN, handler->mcstack()->Particle(mLabelN)->GetPdgCode(),
			 handler->mcstack()->Particle(posLabel)->GetPdgCode(), handler->mcstack()->IsPhysicalPrimary(posLabel), handler->mcstack()->IsSecondaryFromWeakDecay(posLabel), handler->mcstack()->IsSecondaryFromMaterial(posLabel),
				handler->mcstack()->Particle(negLabel)->GetPdgCode(), handler->mcstack()->IsPhysicalPrimary(negLabel), handler->mcstack()->IsSecondaryFromWeakDecay(negLabel), handler->mcstack()->IsSecondaryFromMaterial(negLabel) 
			 );
		}

	}*/

	/*Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i(0); i < iTracks; i++) {                 // loop ove rall these tracks
		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));         // get a track (type AliAODTrack) from the event
		//if(!track || !track->TestFilterBit(1)) continue;                            // if we failed, skip this track
		fHistPt->Fill(track->Pt());                     // plot the pt value of the track in a histogram
	}                                                   // continue until all the tracks are processed
	PostData(1, fOutputList);*/                           // stream the results the analysis of this event to
														// the output manager which will take care of writing
														// it to a file

	for (Int_t iAna = 0; iAna < handler->getNAnalyses(); iAna++) {
		PostData(1+iAna, handler->analysis(iAna)->dirFile()->GetList());
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::Terminate(Option_t *)
{
	// terminate
	// called at the END of the analysis (when all events are processed)

	//handler->Finish();
}
//_____________________________________________________________________________
