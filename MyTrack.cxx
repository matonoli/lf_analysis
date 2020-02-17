#include "MyTrack.h"
#include "MyHandler.h"
#include "AliESDtrackCuts.h"

ClassImp(MyTrack)

#if INPUTFORMAT == 2
Float_t MyTrack::GetDCApvXY() {
	Float_t dca_xy = 0;
	Float_t dca_z = 0;
	mAliTrack->GetImpactParameters(dca_xy,dca_z);
	return dca_xy;
}

Int_t MyTrack::CalculateFlag() {

	if (mFlagCal) return mFlag;
	mFlag = 0;
	if (mHandler->trackCuts2010()->AcceptTrack(mAliTrack)) mFlag += 1;
	if (mHandler->trackCuts2011()->AcceptTrack(mAliTrack)) mFlag += 2;
	if (mHandler->trackCutsTPCOnly()->AcceptTrack(mAliTrack)) mFlag += 4;
	if (mHandler->trackCuts2011sys()->AcceptTrack(mAliTrack)) mFlag += 8;
	if (mHandler->trackCutsV0d()->AcceptTrack(mAliTrack)) mFlag += 16;

	mFlagCal = true;
	return mFlag;
}

#endif