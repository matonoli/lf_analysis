#include <TChain.h>

#include "MyAnalysis.h"
#include "MyHandler.h"

ClassImp(MyAnalysis)

Int_t MyAnalysis::SetHandler(MyHandler* h) {
	
	//printf("Setting up handler for the analysis...\n");
	this->mHandler = h;
	return 0;
}
