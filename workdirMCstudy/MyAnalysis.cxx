#include <TChain.h>

#include "MyAnalysis.h"
#include "MyHandler.h"
#include <TDirectoryFile.h>
#include <TKey.h>
#include <TCollection.h>
#include <TClass.h>
#include <TFile.h>
#include <TROOT.h>

ClassImp(MyAnalysis)

Int_t MyAnalysis::SetHandler(MyHandler* h) {
	
	//printf("Setting up handler for the analysis...\n");
	this->mHandler = h;
	return 0;
}

Int_t MyAnalysis::SetDirectory() {
	
	printf("Setting up the directory for the analysis...\n");
	TString dfName(this->GetName());
	dfName = Form("%s_%i",dfName.Data(),mHandler->nAnalysis());
	mDirFile = new TDirectoryFile(dfName,dfName,"");//,mHandler->file());
	
	return 0;
}

Int_t MyAnalysis::TakeoverHistograms(const Char_t* sourceName) {
	
	TString sourceNameStr(sourceName);
	printf("Copying histograms for the analysis %s ...\n",sourceNameStr.Data());
	TDirectoryFile* source = (TDirectoryFile*)mHandler->filehist()->Get(sourceNameStr.Data());
	printf("File in %s and %s \n", mHandler->filehist()->GetName(),source->GetName());
	mDirFile->cd();
	printf("File out %s \n", mDirFile->GetName());
    //Loop on all entries of the directory
    TKey *key;
    TIter nextkey(source->GetListOfKeys());
    while ((key = (TKey*)nextkey())) {
    	const char *classname = key->GetClassName();
       	TClass *cl = mHandler->root()->GetClass(classname);
       	if (!cl) continue;
       	if (cl->InheritsFrom(TTree::Class())) {
       		TTree *T = (TTree*)source->Get(key->GetName());
       		mDirFile->cd();
         	TTree *newT = T->CloneTree(-1,"fast");
          	newT->Write();
       	} else {
        	source->cd();
        	TObject *obj = key->ReadObj();
        	mDirFile->cd();
        	obj->Write();
        	delete obj;
		}
	}
	mDirFile->SaveSelf(kTRUE);
	mHandler->file()->cd();

	return 0;
}
