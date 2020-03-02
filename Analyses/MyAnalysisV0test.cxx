#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TList.h>
#include <TFile.h>
#include <TLegend.h>
#include <TNamed.h>
#include <THashList.h>
#include <TProfile.h>

#include "MyAnalysisV0test.h"
#include "../MyEvent.h"
#include "../MyTrack.h"
#include "../MyParticle.h"
#include "../MyV0.h"
#include "../MyHandler.h"

#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

//#include <AliAnalysisPIDV0.h>
using namespace V0consts;

ClassImp(MyAnalysisV0test)

MyAnalysisV0test::MyAnalysisV0test() {

}

Int_t MyAnalysisV0test::Init() {

	TString dfName(this->GetName());
	dfName = Form("%s_%i",dfName.Data(),mHandler->nAnalysis());
	mDirFile = new TDirectoryFile(dfName,dfName,"",mHandler->directory());
	mDirFile->cd();

	printf("Initialising analysis %s  \n", 
		this->GetName());

	TH1::SetDefaultSumw2(1);
	TH1::AddDirectory(kTRUE);
	CreateHistograms();
	

	return 0;
}

Int_t MyAnalysisV0test::Make(Int_t iEv) {

    if (!mHandler->event()) return 1;
    MyEvent event(mHandler->event());

    //cout << event.IsGoodAliEvent() << endl;
    if (!event.IsGoodAliEvent()) return 0;
    //hTest->Fill(event.GetZ());
    //cout << event.HasVertex() << " " << event.GetRefMult() << endl;
    //hTest->Fill(event.GetV0MCentrality())

    

    /*Int_t nTracks = mHandler->getNtracks();           // see how many tracks there are in the event
    for(Int_t iTr = 0; iTr < nTracks; iTr++) {                 // loop ove rall these tracks
        
        if (!mHandler->track(iTr)) continue;
        MyTrack t(mHandler->track(iTr));
        t.SetHandler(mHandler);
        

        //if (t.CalculateFlag()>0 && t.CalculateFlag()!=16) cout << t.GetEta() << " " <<  t.CalculateFlag() << endl;
        //hTest->Fill(t.GetDCApvXY());
        //if(!track || !track->TestFilterBit(1)) continue;                            // if we failed, skip this track
        //hTest->Fill(track->Pt());                     // plot the pt value of the track in a histogram
    }*/

    /*Int_t nV0s = mHandler->getNv0s();           // see how many tracks there are in the event
    for(Int_t iV0 = 0; iV0 < nV0s; iV0++) {

    	if (!mHandler->v0(iV0)) continue;
    	MyV0 v0(mHandler->v0(iV0));
    	v0.SetHandler(mHandler);

    	if (v0.HasFastSignal()) {cout << v0.GetPt() << endl;
    	cout << v0.HasFastSignal() << endl;
    	cout << "--" << endl;
    	hTest->Fill(v0.GetIMK0s());}
    }*/

	return 0;
}


Bool_t MyAnalysisV0test::CreateHistograms() {
	
	

	hTest		= new TH1D("hTest","",250,-5,5);
	

	cout << "Creating histograms \n";

}



Int_t MyAnalysisV0test::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();


	return 0;	
}
