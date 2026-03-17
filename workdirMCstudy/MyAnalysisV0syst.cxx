#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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
#include <TNtuple.h>
#include <TCutG.h>
#include <TCut.h>
#include <TGaxis.h>

#include "MyAnalysisV0extract.h"
#include "MyAnalysisV0syst.h"
#include "MyEvent.h"
#include "MyTrack.h"
#include "MyParticle.h"
#include "MyV0.h"
#include "MyHandler.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooPlot.h"

//#include <AliAnalysisPIDV0.h>
using namespace V0consts;
using namespace RooFit;
using namespace std;

ClassImp(MyAnalysisV0syst)



MyAnalysisV0syst::MyAnalysisV0syst() {

}

Int_t MyAnalysisV0syst::Init() {

	mDirFile->cd();

	printf("Initialising analysis %s  \n", 
		this->GetName());

	TH1::SetDefaultSumw2(1);
	BorrowHistograms();
	CreateHistograms();
	


	mList = (TList*)mHandler->directory()->GetList();

	return 0;
}

Int_t MyAnalysisV0syst::Make(Int_t iEv) {
	return 0;
}


Bool_t MyAnalysisV0syst::BorrowHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	//for (int iMu = 0; iMu < NMULTI; ++iMu)		{
		hV0IMvRadiusL[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvRadiusL_%s",SPECIES[iSp]));
		hV0IMvDCAdd[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvDCAdd_%s",SPECIES[iSp]));
		hV0IMvCPA[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvCPA_%s",SPECIES[iSp]));
		hV0IMvFastSignal[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvFastSignal_%s",SPECIES[iSp]));
		hV0IMvCompMass[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvCompMass_%s",SPECIES[iSp]));
		if (iSp>1) hV0IMvLifetime[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvLifetime_%s",SPECIES[iSp]));
		hV0IMvNSigmaTPC[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvNSigmaTPC_%s",SPECIES[iSp]));
		hV0IMvDCAPVpos[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvDCAPVpos_%s",SPECIES[iSp]));
		hV0IMvDCAPVneg[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvDCAPVneg_%s",SPECIES[iSp]));
		hV0IMvNCluster[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvNCluster_%s",SPECIES[iSp]));
		hV0IMvNClusterF[iSp]			
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvNClusterF_%s",SPECIES[iSp]));
	}


	

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
		if (iMu == 0 && iSph > 0) continue;
	
		hV0IMvPtSys[iSp][iMu][iSph][sysRadiusL][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysRadiusL],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysRadiusL][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysRadiusL][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAdd][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAdd],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAdd][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysDCAdd][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysCPA][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysCPA],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysCPA][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysCPA][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysFastSignal][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysFastSignal],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysFastSignal][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysFastSignal][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysCompMass][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysCompMass],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysCompMass][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysCompMass][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysLifetime][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysLifetime],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysLifetime][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysLifetime][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysNSigmaTPC][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNSigmaTPC],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysNSigmaTPC][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysNSigmaTPC][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVpos][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAPVpos],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVpos][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVpos][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVneg][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAPVneg],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVneg][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysDCAPVneg][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysNCluster][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNCluster],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysNCluster][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysNCluster][iVar]);	
		hV0IMvPtSys[iSp][iMu][iSph][sysNClusterF][iVar]
			= (TH2F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNClusterF],SYSTVAR[iVar]));
		hV0IMvPtSys[iSp][iMu][iSph][sysNClusterF][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hV0IMvPtSys[iSp][iMu][iSph][sysNClusterF][iVar]);	
		
	}	}	}	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{

		hV0IMvPtNtSys[iSp][iReg][iSo][iVar] 
			= (TH3F*)mHandler->analysis(0)->dirFile()->Get(Form("hV0IMvPtNtSys_%s_%s_%s_%s",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo],SYSTVAR[iVar]));

		hV0IMvPtNtSys[iSp][iReg][iSo][iVar] = ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH3(hV0IMvPtNtSys[iSp][iReg][iSo][iVar]);	

	}	}	}	}

	// LOAD FROM MC
	if (!mFileMC) {
		printf("No MC file loaded! Comparisons of systematics with MC not performed.\n");
		return 0;
	}

	TDirectoryFile* dirFile1 = new TDirectoryFile("mcFile","mcFile","",mHandler->file());
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("TDirectoryFile")) {
		cout << "Doing systematics from a TDirectoryFile" << endl;
		dirFile1 = (TDirectoryFile*)mFileMC->Get("MyAnalysisV0_0");}
	if (mFileMC->Get("MyAnalysisV0_0")->ClassName() == string("THashList")) {
		cout << "Doing systematics from a THashList" << endl;
		THashList* hashList = (THashList*)mFileMC->Get("MyAnalysisV0_0");
		while (hashList->GetEntries()) {
			dirFile1->Append(hashList->First());
			hashList->RemoveFirst();
		}
	}
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	//for (int iMu = 0; iMu < NMULTI; ++iMu)		{
		hMCV0IMvRadiusL[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvRadiusL_%s",SPECIES[iSp]));
		hMCV0IMvDCAdd[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvDCAdd_%s",SPECIES[iSp]));
		hMCV0IMvCPA[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvCPA_%s",SPECIES[iSp]));
		hMCV0IMvFastSignal[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvFastSignal_%s",SPECIES[iSp]));
		hMCV0IMvCompMass[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvCompMass_%s",SPECIES[iSp]));
		if (iSp>1) hMCV0IMvLifetime[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvLifetime_%s",SPECIES[iSp]));
		hMCV0IMvNSigmaTPC[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvNSigmaTPC_%s",SPECIES[iSp]));
		hMCV0IMvDCAPVpos[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvDCAPVpos_%s",SPECIES[iSp]));
		hMCV0IMvDCAPVneg[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvDCAPVneg_%s",SPECIES[iSp]));
		hMCV0IMvNCluster[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvNCluster_%s",SPECIES[iSp]));
		hMCV0IMvNClusterF[iSp]			
			= (TH2F*)dirFile1->Get(Form("hV0IMvNClusterF_%s",SPECIES[iSp]));

	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < 1; ++iMu)	{
	for (int iSph = 0; iSph < 1; ++iSph)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iMu == 0 && iSph > 0) continue;

		hMCV0IMvPtSys[iSp][iMu][iSph][sysRadiusL][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysRadiusL],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysRadiusL][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysRadiusL][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAdd][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAdd],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAdd][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAdd][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysCPA][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysCPA],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysCPA][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysCPA][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysFastSignal][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysFastSignal],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysFastSignal][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysFastSignal][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysCompMass][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysCompMass],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysCompMass][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysCompMass][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysLifetime][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysLifetime],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysLifetime][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysLifetime][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysNSigmaTPC][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNSigmaTPC],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysNSigmaTPC][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysNSigmaTPC][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAPVpos][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAPVpos],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAPVpos][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAPVpos][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAPVneg][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysDCAPVneg],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAPVneg][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysDCAPVneg][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysNCluster][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNCluster],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysNCluster][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysNCluster][iVar]);	
		hMCV0IMvPtSys[iSp][iMu][iSph][sysNClusterF][iVar]
			= (TH2F*)dirFile1->Get(Form("hV0IMvPtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[sysNClusterF],SYSTVAR[iVar]));
		hMCV0IMvPtSys[iSp][iMu][iSph][sysNClusterF][iVar]	= ((MyAnalysisV0*)mHandler->analysis(0))->RebinTH2(hMCV0IMvPtSys[iSp][iMu][iSph][sysNClusterF][iVar]);	
		

	}	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hMCV0Pt[iSp] = (TH1F*)dirFile1->Get(Form("hV0Pt_%s_MC_MB_MB",SPECIES[iSp]));

		if (hMCV0Pt[iSp]->GetNbinsX() != NPTBINS) 
			hMCV0Pt[iSp] = (TH1F*)hMCV0Pt[iSp]->Rebin(NPTBINS,hMCV0Pt[iSp]->GetName(),XBINS);

	
	}

	return true;

}

Bool_t MyAnalysisV0syst::CreateHistograms() {

	// RAW YIELD FRACTION LOSS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hV0YieldvRadiusL[iSp] = (TH1F*)hV0IMvRadiusL[iSp]->ProjectionX(
			Form("hV0YieldvRadiusL_%s",SPECIES[iSp]),0,-1);
		hV0YieldvDCAdd[iSp] = (TH1F*)hV0IMvDCAdd[iSp]->ProjectionX(
			Form("hV0YieldvDCAdd_%s",SPECIES[iSp]),0,-1);
		hV0YieldvCPA[iSp] = (TH1F*)hV0IMvCPA[iSp]->ProjectionX(
			Form("hV0YieldvCPA_%s",SPECIES[iSp]),0,-1);
		hV0YieldvFastSignal[iSp] = (TH1F*)hV0IMvFastSignal[iSp]->ProjectionX(
			Form("hV0YieldvFastSignal_%s",SPECIES[iSp]),0,-1);
		hV0YieldvCompMass[iSp] = (TH1F*)hV0IMvCompMass[iSp]->ProjectionX(
			Form("hV0YieldvCompMass_%s",SPECIES[iSp]),0,-1);
		if (iSp>1) hV0YieldvLifetime[iSp] = (TH1F*)hV0IMvLifetime[iSp]->ProjectionX(
			Form("hV0YieldvLifetime_%s",SPECIES[iSp]),0,-1);
		hV0YieldvNSigmaTPC[iSp] = (TH1F*)hV0IMvNSigmaTPC[iSp]->ProjectionX(
			Form("hV0YieldvNSigmaTPC_%s",SPECIES[iSp]),0,-1);
		hV0YieldvDCAPVpos[iSp] = (TH1F*)hV0IMvDCAPVpos[iSp]->ProjectionX(
			Form("hV0YieldvDCAPVpos_%s",SPECIES[iSp]),0,-1);
		hV0YieldvDCAPVneg[iSp] = (TH1F*)hV0IMvDCAPVneg[iSp]->ProjectionX(
			Form("hV0YieldvDCAPVneg_%s",SPECIES[iSp]),0,-1);
		hV0YieldvNCluster[iSp] = (TH1F*)hV0IMvNCluster[iSp]->ProjectionX(
			Form("hV0YieldvNCluster_%s",SPECIES[iSp]),0,-1);
		hV0YieldvNClusterF[iSp] = (TH1F*)hV0IMvNClusterF[iSp]->ProjectionX(
			Form("hV0YieldvNClusterF_%s",SPECIES[iSp]),0,-1);
	}

	// CORRECTED VARIED YIELDS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iSp==1 && iSo==sysLifetime) continue;
	if (iMu == 0 && iSph > 0) continue;

		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);
		
		hV0PtSys[iSp][iMu][iSph][iSo][iVar] = (TH1F*)hpt->Clone(Form("hV0PtSys_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[iSo],SYSTVAR[iVar]));
		hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar] = (TH1F*)hpt->Clone(Form("hV0PtSysRatioToHM_%s_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[iSo],SYSTVAR[iVar]));

	}	}	}	}	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);
		hV0PtSysSigExLoose[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysSigExLoose_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysSigExTight[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysSigExTight_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysMaxDSigEx[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDSigEx_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysSigExLooseRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysSigExTightRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDSigExRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		hV0PtSysFeeddown[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysFeeddown_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysMaxDFeeddown[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDFeeddown_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysFeeddownXiErr[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysFeeddownXiErr_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysMaxDFeeddownXiErr[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDFeeddownXiErr_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDFeeddownTotal_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		hV0PtSysFeeddownRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysFeeddownRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysMaxDFeeddownRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDFeeddownRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysFeeddownXiErrRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysFeeddownXiErrRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysMaxDFeeddownXiErrRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDFeeddownXiErrRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hV0PtSysMaxDFeeddownTotalRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		hFracBudget[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracBudget_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracEffi[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracEffi_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracExpBias[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracExpBias_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracCuts[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracCuts_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracSigEx[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracSigEx_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracFD[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracFD_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		hFracEffiRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracEffiRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracExpBiasRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracExpBiasRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracCutsRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracCutsRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracSigExRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracSigExRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hFracFDRatioToHM[iSp][iMu][iSph]	= (TH1F*)hpt->Clone(Form("hFracFDRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

	}	}	}

	// FLAT SYSTEMATICS
	{
		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[1]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);
		hV0PtSysBudget = (TH1F*)hpt->Clone(Form("hV0PtSysBudget"));
		for (int i=1; i<hV0PtSysBudget->GetNbinsX()+1; i++) hV0PtSysBudget->SetBinContent(i,0.04);

		hV0PtSysEffi = (TH1F*)hpt->Clone(Form("hV0PtSysEffi"));
		for (int i=1; i<hV0PtSysEffi->GetNbinsX()+1; i++) hV0PtSysEffi->SetBinContent(i,0.02);

		hV0PtSysExpBiasJetty = (TH1F*)hpt->Clone(Form("hV0PtSysExpBiasJetty"));
		for (int i=1; i<hV0PtSysExpBiasJetty->GetNbinsX()+1; i++) hV0PtSysExpBiasJetty->SetBinContent(i, (XBINS[i-1]<1.)? 0.: 0.04);	

		hV0PtSysExpBiasIso = (TH1F*)hpt->Clone(Form("hV0PtSysExpBiasIso"));
		for (int i=1; i<hV0PtSysExpBiasIso->GetNbinsX()+1; i++) hV0PtSysExpBiasIso->SetBinContent(i, (XBINS[i-1]<1.)? 0.: 0.01);	
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	if (iSp==1 && iSo==sysLifetime) continue;
	if (iMu == 0 && iSph > 0) continue;

		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);
		
		hV0PtSysMaxD[iSp][iMu][iSph][iSo] = (TH1F*)hpt->Clone(Form("hV0PtSysMaxD_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[iSo]));
		hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo] = (TH1F*)hpt->Clone(Form("hV0PtSysMaxDRatioToHM_%s_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph],SYSTS[iSo]));

	}	}	}	}

	// ROGER BARLOW CHECKS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iVar!=loosest && iVar!=tightest) continue;

		hRBcheck[iSp][iSo][iVar] = new TH1F(Form("hRBcheck_%s_%s_%s",SPECIES[iSp],SYSTS[iSo],SYSTVAR[iVar]), 
			Form("%s;#Delta/#sigma_{cc};Entries",SYSTS[iSo]), 120,-6.,6.);

	}	}	}


	// MC

	if (!mFileMC) {
		printf("No MC file loaded! Comparisons of systematics with MC not performed.\n");
		return 0;
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		hMCV0YieldvRadiusL[iSp] = (TH1F*)hMCV0IMvRadiusL[iSp]->ProjectionX(
			Form("hMCV0YieldvRadiusL_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvDCAdd[iSp] = (TH1F*)hMCV0IMvDCAdd[iSp]->ProjectionX(
			Form("hMCV0YieldvDCAdd_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvCPA[iSp] = (TH1F*)hMCV0IMvCPA[iSp]->ProjectionX(
			Form("hMCV0YieldvCPA_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvFastSignal[iSp] = (TH1F*)hMCV0IMvFastSignal[iSp]->ProjectionX(
			Form("hMCV0YieldvFastSignal_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvCompMass[iSp] = (TH1F*)hMCV0IMvCompMass[iSp]->ProjectionX(
			Form("hMCV0YieldvCompMass_%s",SPECIES[iSp]),0,-1);
		if (iSp>1) hMCV0YieldvLifetime[iSp] = (TH1F*)hMCV0IMvLifetime[iSp]->ProjectionX(
			Form("hMCV0YieldvLifetime_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvNSigmaTPC[iSp] = (TH1F*)hMCV0IMvNSigmaTPC[iSp]->ProjectionX(
			Form("hMCV0YieldvNSigmaTPC_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvDCAPVpos[iSp] = (TH1F*)hMCV0IMvDCAPVpos[iSp]->ProjectionX(
			Form("hMCV0YieldvDCAPVpos_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvDCAPVneg[iSp] = (TH1F*)hMCV0IMvDCAPVneg[iSp]->ProjectionX(
			Form("hMCV0YieldvDCAPVneg_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvNCluster[iSp] = (TH1F*)hMCV0IMvNCluster[iSp]->ProjectionX(
			Form("hMCV0YieldvNCluster_%s",SPECIES[iSp]),0,-1);
		hMCV0YieldvNClusterF[iSp] = (TH1F*)hMCV0IMvNClusterF[iSp]->ProjectionX(
			Form("hMCV0YieldvNClusterF_%s",SPECIES[iSp]),0,-1);
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);
		hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hpt->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysSum[iSp][iMu][iSph]->GetYaxis()->SetTitle("Total relative syst. uncertainty");

		hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hpt->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtSysSumUnc[iSp][iMu][iSph]->GetYaxis()->SetTitle("Uncorrelated syst (w.r.t. spherocity)");

	}	}	}


	return true;
}


Bool_t MyAnalysisV0syst::CloneHistograms() {

	Int_t nType = (mHandler->GetFlagMC()) ? 2 : 1;


	// NT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
		
		TH2F* hptnt = (TH2F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtNtFit_%s_D_Trans",SPECIES[iSp]));
		hV0PtNtSys[iSp][iReg][iSo][iVar] = (TH2F*)hptnt->Clone(Form("hV0PtNtSys_%s_%s_%s_%s",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo],SYSTVAR[iVar]));

	}	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
		
		TH2F* hptnt = (TH2F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtNtFit_%s_D_Trans",SPECIES[iSp]));
		hV0PtNtSysMaxD[iSp][iReg][iSo] = (TH2F*)hptnt->Clone(Form("hV0PtNtSysMaxD_%s_%s_%s",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo]));

	}	}	}

	// RT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; iRtBin++) {
		
		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);

		//hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin] = (TH1F*)hpt->Clone(Form("hV0PtRtSys_%s_%s_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo],SYSTVAR[iVar],iRtBin));
		//hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->Reset();

	}	}	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; iRtBin++) {
		
		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);
		hpt->Reset();

		hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin] = (TH1F*)hpt->Clone(Form("hV0PtRtSysMaxD_%s_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo],iRtBin));
		hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin] = (TH1F*)hpt->Clone(Form("hV0PtRtSysMaxDRatioToHM_%s_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo],iRtBin));


	}	}	}	}


	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; iRtBin++) {

		TH1F* hpt = (TH1F*)mHandler->analysis(1)->dirFile()->Get(Form("hV0PtFit_%s_D_MB_MB",SPECIES[iSp]));
		if (hpt->GetNbinsX() != NPTBINS) hpt = (TH1F*)hpt->Rebin(NPTBINS,hpt->GetName(),XBINS);
		hpt->Reset();

		hV0PtRtSysSigExLoose[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSigExLoose_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSigExDeflt_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSigExTight[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSigExTight_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysMaxDSigEx[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysMaxDSigEx_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSigExLooseRatioToHM[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSigExLooseRatioToHM_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSigExDefltRatioToHM[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSigExDefltRatioToHM_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSigExTightRatioToHM[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSigExTightRatioToHM_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysMaxDSigExRatioToHM[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysMaxDSigExRatioToHM_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));

		hV0PtRtSysSum[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSum_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSumUnc[iSp][iReg][iRtBin]	= (TH1F*)hpt->Clone(Form("hV0PtRtSysSumUnc_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],iRtBin));


	}	}	}


	return true;


}

void MyAnalysisV0syst::SetMCInputFile(const Char_t *name) {

	TString fileName = TString(name);
	if (fileName.Data() != "") {
		mFileMC = new TFile(fileName,"READ");
		printf("MC File %s loaded in. \n", fileName.Data()); }
	else {
		printf("No MC file loaded.");
	}
}

Int_t MyAnalysisV0syst::Finish() {

	printf("Finishing analysis %s \n",this->GetName());
	mDirFile->cd();

	cout << "GGothere 1" << endl;
	if (!((MyAnalysisV0extract*)mHandler->analysis(2))->IsRun()) {
		cout << "Taking over sidebands" << endl;
		((MyAnalysisV0extract*)mHandler->analysis(2))->TakeoverSidebands();
		mDirFile->cd();	
	}

	CloneHistograms();

	
	//StudyRawYieldLoss();
	MakeEfficiencies();

	MakeCorrectedYields();
	BinHistogramsIntoRTrec();

	//MakeBarlowChecksPt();
	MakeDeviations();
	//CalculateAndPlotMaxDeviations();
	CalculateMaxDeviationsNt();

	//CalculateSignalExSys();
	CalculateSignalExSysNt();
	
	//CalculateFeeddownSys();

	//AddDeviations();
	AddDeviationsNt();


	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}

		if (objName.BeginsWith("hEffi_")) {
			mDirFile->Remove(obj);	}
	}
		
	return 0;	
}

void MyAnalysisV0syst::StudyRawYieldLoss() {

	enum { rising, sinking };
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		ProcessRawYieldLossHist(hV0IMvRadiusL[iSp],hV0YieldvRadiusL[iSp],iSp,0.3,sinking);
		ProcessRawYieldLossHist(hV0IMvDCAdd[iSp],hV0YieldvDCAdd[iSp],iSp,1.5,rising);
		ProcessRawYieldLossHist(hV0IMvCPA[iSp],hV0YieldvCPA[iSp],iSp,(iSp>1)?0.993:0.95,sinking);
		ProcessRawYieldLossHist(hV0IMvFastSignal[iSp],hV0YieldvFastSignal[iSp],iSp,1.,sinking);
		ProcessRawYieldLossHist(hV0IMvCompMass[iSp],hV0YieldvCompMass[iSp],iSp,2.5,sinking);
		if (iSp>1) ProcessRawYieldLossHist(hV0IMvLifetime[iSp],hV0YieldvLifetime[iSp],iSp,40.,rising);
		ProcessRawYieldLossHist(hV0IMvNSigmaTPC[iSp],hV0YieldvNSigmaTPC[iSp],iSp,6.5,rising);
		ProcessRawYieldLossHist(hV0IMvDCAPVpos[iSp],hV0YieldvDCAPVpos[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hV0IMvDCAPVneg[iSp],hV0YieldvDCAPVneg[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hV0IMvNCluster[iSp],hV0YieldvNCluster[iSp],iSp,70.,sinking);
		ProcessRawYieldLossHist(hV0IMvNClusterF[iSp],hV0YieldvNClusterF[iSp],iSp,0.8,sinking);
	}


	if (!mFileMC) {
		printf("No MC file loaded! Comparisons of systematics with MC not performed.\n");
		return;
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		ProcessRawYieldLossHist(hMCV0IMvRadiusL[iSp],hMCV0YieldvRadiusL[iSp],iSp,0.3,sinking);
		ProcessRawYieldLossHist(hMCV0IMvDCAdd[iSp],hMCV0YieldvDCAdd[iSp],iSp,1.5,rising);
		ProcessRawYieldLossHist(hMCV0IMvCPA[iSp],hMCV0YieldvCPA[iSp],iSp,(iSp>1)?0.993:0.95,sinking);
		ProcessRawYieldLossHist(hMCV0IMvFastSignal[iSp],hMCV0YieldvFastSignal[iSp],iSp,1.,sinking);
		ProcessRawYieldLossHist(hMCV0IMvCompMass[iSp],hMCV0YieldvCompMass[iSp],iSp,2.5,sinking);
		if (iSp>1) ProcessRawYieldLossHist(hMCV0IMvLifetime[iSp],hMCV0YieldvLifetime[iSp],iSp,40.,rising);
		ProcessRawYieldLossHist(hMCV0IMvNSigmaTPC[iSp],hMCV0YieldvNSigmaTPC[iSp],iSp,6.5,rising);
		ProcessRawYieldLossHist(hMCV0IMvDCAPVpos[iSp],hMCV0YieldvDCAPVpos[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hMCV0IMvDCAPVneg[iSp],hMCV0YieldvDCAPVneg[iSp],iSp,0.05,sinking);
		ProcessRawYieldLossHist(hMCV0IMvNCluster[iSp],hMCV0YieldvNCluster[iSp],iSp,70.,sinking);
		ProcessRawYieldLossHist(hMCV0IMvNClusterF[iSp],hMCV0YieldvNClusterF[iSp],iSp,0.8,sinking);
	}

	mHandler->root()->SetBatch(kTRUE);
	TCanvas* cRYL[NSPECIES]; 
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		cRYL[iSp] = new TCanvas(Form("cRYL_%s",SPECIES[iSp]), "", 2500, 1250);
		cRYL[iSp]->Divide(4,3,5e-5,5e-5);
		Int_t padC = 1; cRYL[iSp]->cd(padC++);

		DrawRawYieldLossHist(hV0YieldvRadiusL[iSp],hMCV0YieldvRadiusL[iSp],(iSp>1)?0.05:0.15); 
		//DrawVariation(0.3,kRed,cRYL[iSp]->GetPad(padC)); DrawVariation(0.4,kBlue,cRYL[iSp]->GetPad(padC)); DrawVariation(0.5,kBlack,cRYL[iSp]->GetPad(padC)); DrawVariation(0.6,kGreen+2,cRYL[iSp]->GetPad(padC));
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysRadiusL][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
		cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvDCAdd[iSp],hMCV0YieldvDCAdd[iSp],1.0); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysDCAdd][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvCPA[iSp],hMCV0YieldvCPA[iSp],0.2); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysCPA][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvFastSignal[iSp],hMCV0YieldvFastSignal[iSp],0.85); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysFastSignal][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvCompMass[iSp],hMCV0YieldvCompMass[iSp],(iSp>1)?0.20:0.15); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysCompMass][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		if (iSp>1) DrawRawYieldLossHist(hV0YieldvLifetime[iSp],hMCV0YieldvLifetime[iSp],(iSp>1)?1.:0.15); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysLifetime][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvNSigmaTPC[iSp],hMCV0YieldvNSigmaTPC[iSp],0.6); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysNSigmaTPC][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvDCAPVpos[iSp],hMCV0YieldvDCAPVpos[iSp],(iSp==2)?0.15:0.04); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysDCAPVpos][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvDCAPVneg[iSp],hMCV0YieldvDCAPVneg[iSp],(iSp==3)?0.15:0.04); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysDCAPVneg][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvNCluster[iSp],hMCV0YieldvNCluster[iSp],0.03); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysNCluster][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		DrawRawYieldLossHist(hV0YieldvNClusterF[iSp],hMCV0YieldvNClusterF[iSp],0.25); 
		for (int iVar = 0; iVar < sysVarsSizeof; iVar++) DrawVariation(sysVar[iSp-1][sysNClusterF][iVar],kBlue,2,cRYL[iSp]->GetPad(padC));
			cRYL[iSp]->cd(padC++);
		

		cRYL[iSp]->Write();
		cRYL[iSp]->SaveAs(Form("plots/cRYL_%s.png",SPECIES[iSp]));
	}

	mHandler->root()->SetBatch(kFALSE);
}



void MyAnalysisV0syst::ProcessRawYieldLossHist(TH2F* hist, TH1F* yieldhist, Int_t Sp, Double_t loose, Int_t opt) {

	Double_t* yieldLoose =
		((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSB((TH1F*)hist->ProjectionY(Form("Sys_iSp%i_iBin0",Sp),
			hist->GetXaxis()->FindBin(loose),hist->GetXaxis()->FindBin(loose)));
	//cout << *(yieldLoose) << " +- " << *(yieldLoose+1) << endl;
	Double_t cntL = *(yieldLoose);

	for (int iBin = 1; iBin < hist->GetNbinsX()+1; ++iBin)	{
		Double_t* yield =
			((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSB((TH1F*)hist->ProjectionY(Form("Sys_iSp%i_iBin0",Sp),
			iBin,iBin));

		Double_t cnt = *(yield);
		Double_t fractionY = ( cntL < cnt || cntL == 0) ? 0. :
			(1. - cnt/cntL);
		
		yieldhist->SetBinContent(iBin,fractionY);
		yieldhist->SetBinError(iBin,0);			 
	}

	// Remove junk from the file memory
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
		}
}

void MyAnalysisV0syst::DrawRawYieldLossHist(TH1F* da, TH1F* mc, Double_t ymax) {

	mHandler->MakeNiceHistogram(da,kBlack);
	mHandler->MakeNiceHistogram(mc,kRed);
	da->SetMarkerSize(0.9); mc->SetMarkerSize(0.9);
	da->GetYaxis()->SetTitle("Signal loss fraction w.r.t. loosest cut");

	da->GetYaxis()->SetRangeUser(0.,ymax);
	da->Draw("p");
	mc->Draw("p same");

}

void MyAnalysisV0syst::DrawVariation(Double_t cut, Int_t col, Int_t styl, TVirtualPad* can) {

	Double_t x[2] = {cut, cut};
	Double_t y[2];
	y[0] = can->GetUymin(); y[1] = can->GetUymax();

	TGraph* gcut = new TGraph(2, x, y);
	gcut->SetLineWidth(2);
	gcut->SetLineColor(col);
	gcut->SetLineStyle(styl);

	gcut->Draw("same");
}

void MyAnalysisV0syst::MakeEfficiencies() {

	mParSigK0s	= new TF1("funcSigK0s","[0]+[1]*x+[2]/x",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParSigK0s->SetParameters(cuts::K0S_PARSIG[0],cuts::K0S_PARSIG[1],cuts::K0S_PARSIG[2]);
	mParMuK0s	= new TF1("funcMuK0s","(x<=1.6)*([0]+[1]*x+[2]*x*x)+(x>1.6)*[3]",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParMuK0s->SetParameters(cuts::K0S_PARMU[0],cuts::K0S_PARMU[1],cuts::K0S_PARMU[2],cuts::K0S_PARMU[3]);
	
	mParSigL	= new TF1("funcSigL","[0]+[1]*x+[2]/x",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParSigL->SetParameters(cuts::L_PARSIG[0],cuts::L_PARSIG[1],cuts::L_PARSIG[2]);
	mParMuL	= new TF1("funcMuL","(x<=1.9)*([0]+[1]*x+[2]*x*x)+(x>1.9)*([3]+[4]*x)",1e-4+XBINS[0],XBINS[NPTBINS]);
	mParMuL->SetParameters(cuts::L_PARMU[0],cuts::L_PARMU[1],cuts::L_PARMU[2],cuts::L_PARMU[3],cuts::L_PARMU[4]);
	
	mParSigK0s->Write();
	mParSigL->Write();
	mParMuK0s->Write();
	mParMuL->Write();

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < 1; ++iMu)	{
	for (int iSph = 0; iSph < 1; ++iSph)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{


		if (iSp==1 && iSo==sysLifetime) continue;
		hV0EfficiencySys[iSp][iSo][iVar] = ProcessEfficiency(hMCV0IMvPtSys[iSp][iMu][iSph][iSo][iVar],hMCV0Pt[iSp],iSp);
		hV0EfficiencySys[iSp][iSo][iVar]->Write();
	}	}	}	}	}

}

TH1F* MyAnalysisV0syst::ProcessEfficiency(TH2F* hist, TH1F* hmc, Int_t Sp) {

	TString hname = "hEffi_";
	hname += hist->GetName();
	TH1F* hrc = (TH1F*)hmc->Clone(hname.Data());

	for (Int_t iBin = 1; iBin < hist->GetNbinsX()+1; iBin++) {
		Double_t pt = hist->GetXaxis()->GetBinCenter(iBin);
		Double_t hi = (Sp==1) ? mParMuK0s->Eval(pt) + cuts::V0_COMP_NSIG*mParSigK0s->Eval(pt)
			:	mParMuL->Eval(pt) + cuts::V0_COMP_NSIG*mParSigL->Eval(pt);
		Double_t lo = (Sp==1) ? mParMuK0s->Eval(pt) - cuts::V0_COMP_NSIG*mParSigK0s->Eval(pt)
			:	mParMuL->Eval(pt) - cuts::V0_COMP_NSIG*mParSigL->Eval(pt);

		Double_t yield =
			hist->Integral(iBin,iBin,hist->GetYaxis()->FindBin(lo),hist->GetYaxis()->FindBin(hi));
		hrc->SetBinContent(iBin,yield);
		hrc->SetBinError(iBin,TMath::Sqrt(yield));
	}

	hrc->Divide(hrc,hmc,1.,1.,"B");

	return hrc;
}

void MyAnalysisV0syst::MakeCorrectedYields() {

	// SPHEROCITY ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
		continue;
	if (iSp==1 && iSo==sysLifetime) continue;
	if (iMu == 0 && iSph > 0) continue;

		Double_t* yield = 0;
		cout << "Extracting histogram " << hV0IMvPtSys[iSp][iMu][iSph][iSo][iVar] << " " << hV0IMvPtSys[iSp][iMu][iSph][iSo][iVar]->GetName() << endl;
		for (int iBin = 1; iBin < NPTBINS+1; ++iBin)	{
			cout << "Extracting pt bin " << iBin << " of " << hV0IMvPtSys[iSp][iMu][iSph][iSo][iVar]->GetName() << endl;
			yield = ((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSB((TH1F*)hV0IMvPtSys[iSp][iMu][iSph][iSo][iVar]->ProjectionY(
					Form("Sys_iSp%i_iBin%i", iSp, iBin), iBin,iBin), false);
				hV0PtSys[iSp][iMu][iSph][iSo][iVar]->SetBinContent(iBin,*(yield+0));
				hV0PtSys[iSp][iMu][iSph][iSo][iVar]->SetBinError(iBin,*(yield+1));
		}

		hV0PtSys[iSp][iMu][iSph][iSo][iVar]->Scale(1,"width");
		hV0PtSys[iSp][iMu][iSph][iSo][iVar]->Divide(hV0EfficiencySys[iSp][iSo][iVar]);

	}	}	}	}	}


	// NT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iSp==1 && iSo==sysLifetime) continue;

		for (int iNt = 1; iNt < hV0IMvPtNtSys[iSp][iReg][iSo][iVar]->GetZaxis()->GetNbins(); ++iNt)	{

			hV0IMvPtNtSys[iSp][iReg][iSo][iVar]->GetZaxis()->SetRange(iNt,iNt);
			TH2F* hXY = (TH2F*)hV0IMvPtNtSys[iSp][iReg][iSo][iVar]->Project3D("yx");

			Double_t* yield = 0;
			cout << "Extracting histogram NtBin " << iNt << " " << hV0IMvPtNtSys[iSp][iReg][iSo][iVar] << " " << hV0IMvPtNtSys[iSp][iReg][iSo][iVar]->GetName() << endl;
			for (int iBin = 1; iBin < NPTBINS+1; ++iBin)	{
			
				yield = ((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSB((TH1F*)hXY->ProjectionY(
						Form("Sys_iNt%i_iSp%i_iBin%i", iNt, iSp, iBin), iBin,iBin),false);
					hV0PtNtSys[iSp][iReg][iSo][iVar]->SetBinContent(iBin,iNt,*(yield+0));
					hV0PtNtSys[iSp][iReg][iSo][iVar]->SetBinError(iBin,iNt,*(yield+1));
			}

			delete hXY;

		}

		//SCALE BY BIN WIDTH
		hV0PtNtSys[iSp][iReg][iSo][iVar] 
			= ((MyAnalysisV0*)mHandler->analysis(0))->ScaleWidthTH2(hV0PtNtSys[iSp][iReg][iSo][iVar]);

		// CORRECT W EFFICIENCY
		hV0PtNtSys[iSp][iReg][iSo][iVar] 
			= ((MyAnalysisV0*)mHandler->analysis(0))->DivideTH2ByTH1(hV0PtNtSys[iSp][iReg][iSo][iVar],hV0EfficiencySys[iSp][iSo][iVar]);

		// MEMORY OVERLOAD
		if (iSp == 2 && !iReg && !iSo && !iVar) {
			TIter objIt(mDirFile->GetList(),kIterForward);
			TObject* obj = 0;
			while ( (obj = objIt()) ) {
				TString objName(obj->GetName());
				if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
					mDirFile->Remove(obj);	}
			}
			
		}

	}	}	}	}

	// Remove junk from the file memory
	TIter objIt(mDirFile->GetList(),kIterForward);
	TObject* obj = 0;
	while ( (obj = objIt()) ) {
		TString objName(obj->GetName());
		if (!objName.BeginsWith("h") && !objName.BeginsWith("c") && !objName.BeginsWith("t") ) {
			mDirFile->Remove(obj);	}
		}



}

void MyAnalysisV0syst::BinHistogramsIntoRTrec() {

	mDirFile->cd();
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iSp==1 && iSo==sysLifetime) continue;
	
		for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

			// REBIN
			int lowedge = GetRTBinRec(iRtBin,kTRUE);
			int upedge  = GetRTBinRec(iRtBin,kFALSE);
			Double_t NormScaleNt = 1.;//hNtUnf->Integral(lowedge,upedge)/hNtUnf->Integral(1,51);

			if (hV0PtNtSys[iSp][iReg][iSo][iVar]) {
				hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]
					= (TH1F*)hV0PtNtSys[iSp][iReg][iSo][iVar]->ProjectionX(Form("hV0PtRtSys_%s_%s_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo],SYSTVAR[iVar],iRtBin),lowedge,upedge);
				hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->Scale(1./NormScaleNt);
			}


			hV0PtRtSysRatioToHM[iSp][iReg][iSo][iVar][iRtBin] 
				= (TH1F*)hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->Clone(Form("hV0PtRtSysRatioToHM_%s_%s_%s_%s_%i",SPECIES[iSp],REGIONS[iReg],SYSTS[iSo],SYSTVAR[iVar],iRtBin));

		}

	}	}	}	}


}

Int_t MyAnalysisV0syst::GetRTBinRec(const int& binRt, bool isLowEdge)
{

	int binNch = -1;

	//! <NT> = 5.742
	
	if( binRt == 0 ){
		if(isLowEdge) binNch = 1;
		else binNch = 50;
	}
	else if( binRt == 1 ){ //! From NT = 0 to NT = 4
		if(isLowEdge) binNch = 1;
		else binNch = 5;//else binNch = 4;
	}
	else if( binRt == 2 ){ //! From NT = 5 to NT = 8
		if(isLowEdge) binNch = 6;
		else binNch = 9;
	}
	else if( binRt == 3 ){//! From NT = 9 to NT = 14
		if(isLowEdge) binNch = 10;
		else binNch = 15;
	}
	else if( binRt == 4 ){//! From NT = 15 to NT = 28
		if(isLowEdge) binNch = 16;
		else binNch = 29;
	}
	else{
		if(isLowEdge) binNch = 30;
		else binNch = 50;
	}

	return binNch;

}

void MyAnalysisV0syst::MakeBarlowChecks() {
	
	mHandler->root()->SetBatch(kTRUE);

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{	
		if (iVar!=loosest && iVar!=tightest) continue;
		if (iSp==1 && iSo==sysLifetime) continue;
		Int_t iSph = 0;

		//cout << hV0PtSys[iSp][iMu][iSo][iVar] << " a " << hV0PtSys[iSp][iMu][iSo][deflt] << " a " << hRBcheck[iSp][iSo][iVar] << endl;

		for (int iBin = 1; iBin < hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetNbinsX()+1; iBin++) {
			Double_t errDiff = hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinError(iBin)*hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinError(iBin) - hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinError(iBin)*hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinError(iBin);
			Double_t sigma = TMath::Sqrt(TMath::Abs(errDiff));
			Double_t delta = hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinContent(iBin) - hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinContent(iBin);

			if (hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinContent(iBin)>0) hRBcheck[iSp][iSo][iVar]->Fill((sigma>0)?delta/sigma:0.);
		}

	}	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
		if (iVar!=loosest && iVar!=tightest) continue;
		
		TCanvas* cRB = new TCanvas(Form("cRB_%s_%s",SPECIES[iSp],SYSTVAR[iVar]),"",2800,2000);
		cRB->Divide(4,3,5e-5,5e-5);
		Int_t padC = 1; 

		for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
		cRB->cd(padC);
		if (iSp==1 && iSo==sysLifetime) continue;
			mHandler->MakeNiceHistogram(hRBcheck[iSp][iSo][iVar],kBlack);
			hRBcheck[iSp][iSo][iVar]->SetLineColor(kBlue+2);
			hRBcheck[iSp][iSo][iVar]->SetLineWidth(1);
			hRBcheck[iSp][iSo][iVar]->Draw("hist");

			cRB->Update();
			DrawVariation(0., kBlack, 1, cRB->GetPad(padC));
			DrawVariation(1., kBlack, 2, cRB->GetPad(padC));
			DrawVariation(-1., kBlack, 2, cRB->GetPad(padC));

			Double_t quantileValues[5] = {0.0, 0.158, 0.5, 0.842, 1.0};
			Double_t quantileCuts[5];
			hRBcheck[iSp][iSo][iVar]->GetQuantiles(5,quantileCuts,quantileValues);		

			DrawVariation(hRBcheck[iSp][iSo][iVar]->GetMean(), kRed, 1, cRB->GetPad(padC));
			DrawVariation(quantileCuts[3], kRed, 2, cRB->GetPad(padC));
			DrawVariation(quantileCuts[1], kRed, 2, cRB->GetPad(padC));
			padC++;
		}
		cRB->Write();
		cRB->SaveAs(Form("plots/cRB_%s_%s.png",SPECIES[iSp],SYSTVAR[iVar]));

	}	}

	mHandler->root()->SetBatch(kFALSE);

}

void MyAnalysisV0syst::MakeBarlowChecksPt() {
	
	mHandler->root()->SetBatch(kTRUE);

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{	
		if (iVar!=loosest && iVar!=tightest) continue;
		if (iSp==1 && iSo==sysLifetime) continue;
		Int_t iSph = 0;

		//cout << hV0PtSys[iSp][iMu][iSo][iVar] << " a " << hV0PtSys[iSp][iMu][iSo][deflt] << " a " << hRBcheck[iSp][iSo][iVar] << endl;

		hRBcheckPt[iSp][iMu][iSo][iVar] = DivideAndComputeRogerBarlow(hV0PtSys[iSp][iMu][iSph][iSo][iVar],hV0PtSys[iSp][iMu][iSph][iSo][deflt]);

		for (int iBin = 1; iBin < hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetNbinsX()+1; iBin++) {
			Double_t errDiff = hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinError(iBin)*hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinError(iBin) - hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinError(iBin)*hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinError(iBin);
			Double_t sigma = TMath::Sqrt(TMath::Abs(errDiff));
			Double_t delta = hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinContent(iBin) - hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinContent(iBin);

			if (hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinContent(iBin)>0) hRBcheck[iSp][iSo][iVar]->Fill((sigma>0)?delta/sigma:0.);
		}


	}	}	}	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
		if (iVar!=loosest && iVar!=tightest) continue;
		
		TCanvas* cRB = new TCanvas(Form("cRB_%s_%s",SPECIES[iSp],SYSTVAR[iVar]),"",2800,2000);
		cRB->Divide(4,3,5e-5,5e-5);
		Int_t padC = 1; 

		for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
		cRB->cd(padC);
		if (iSp==1 && iSo==sysLifetime) continue;
			mHandler->MakeNiceHistogram(hRBcheck[iSp][iSo][iVar],kBlack);
			hRBcheck[iSp][iSo][iVar]->SetLineColor(kBlue+2);
			hRBcheck[iSp][iSo][iVar]->SetLineWidth(1);
			//hRBcheck[iSp][iSo][iVar]->Draw("hist");

			//mHandler->MakeNiceHistogram(hRBcheckPt[iSp][3][iSo][iVar],kBlack);

			cRB->Update();
			//DrawVariation(0., kBlack, 1, cRB->GetPad(padC));
			//DrawVariation(1., kBlack, 2, cRB->GetPad(padC));
			//DrawVariation(-1., kBlack, 2, cRB->GetPad(padC));

			//Double_t quantileValues[5] = {0.0, 0.158, 0.5, 0.842, 1.0};
			//Double_t quantileCuts[5];
			//hRBcheck[iSp][iSo][iVar]->GetQuantiles(5,quantileCuts,quantileValues);		

			//DrawVariation(hRBcheck[iSp][iSo][iVar]->GetMean(), kRed, 1, cRB->GetPad(padC));
			//DrawVariation(quantileCuts[3], kRed, 2, cRB->GetPad(padC));
			//DrawVariation(quantileCuts[1], kRed, 2, cRB->GetPad(padC));
			padC++;
		}
		cRB->Write();
		cRB->SaveAs(Form("plots/cRB_%s_%s.png",SPECIES[iSp],SYSTVAR[iVar]));

	}	}

	mHandler->root()->SetBatch(kFALSE);

}

TH1F* MyAnalysisV0syst::DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){ 
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX(); 
  Double_t lh2NBins = h2->GetNbinsX(); 

  if( lh1NBins != lh2NBins ){ 
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return 0x0;
  }

  TH1F* hrb = (TH1F*)h1->Clone(Form("hRB_%s",h1->GetName()));

  Double_t lSigmaDelta[100]; 
  for( Int_t i=1; i<hrb->GetNbinsX()+1; i++){ 
    //Computation of roger barlow sigma_{delta} 
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(hrb->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot 
    //if ( h2->GetBinContent(i) > 1e-12 ){ 
    //  lSigmaDelta[i] /= h2->GetBinContent(i); 
    //}else{ 
    //  lSigmaDelta[i] = 0; 
    //}
  }
  //Regular Division 
  hrb->Add(h2, -1.); 
  //Replace Erorrs 
  for( Int_t i=1; i<hrb->GetNbinsX()+1; i++){ 
    hrb->SetBinError(i, lSigmaDelta[i]);
  }

  return hrb;  
}

/*TH1F* MyAnalysisV0syst::DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){ 
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX(); 
  Double_t lh2NBins = h2->GetNbinsX(); 

  if( lh1NBins != lh2NBins ){ 
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return 0x0;
  }

  TH1F* hrb = (TH1F*)h1->Clone(Form("hRB_%s",h1->GetName()));

  Double_t lSigmaDelta[100]; 
  for( Int_t i=1; i<hrb->GetNbinsX()+1; i++){ 
    //Computation of roger barlow sigma_{delta} 
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(hrb->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot 
    //if ( h2->GetBinContent(i) > 1e-12 ){ 
    //  lSigmaDelta[i] /= h2->GetBinContent(i); 
    //}else{ 
    //  lSigmaDelta[i] = 0; 
    //}
  }
  //Regular Division 
  hrb->Add(h2, -1.); 
  //Replace Erorrs 
  for( Int_t i=1; i<hrb->GetNbinsX()+1; i++){ 
    hrb->SetBinError(i, lSigmaDelta[i]);
  }

  return hrb;  
}*/

void MyAnalysisV0syst::MakeDeviations() {

	// SPHERO ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{	
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iVar==deflt) continue;
	if (iSp==1 && iSo==sysLifetime) continue;
	if (iMu == 0 && iSph > 0) continue;

		Double_t lSigmaDelta[100]; 
  		for( Int_t i=1; i<hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetNbinsX()+1; i++){ 
    		//Computation of roger barlow sigma_{delta} 
    		lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinError(i),2) - TMath::Power(hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinError(i),2) ) );
 		   //Computation of relationship to h2 for plotting in ratio plot 
    		if ( hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinContent(i) > 1e-12 )	{ 
    			lSigmaDelta[i] /= hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinContent(i);		}
    		else	{ 
    			lSigmaDelta[i] = 0; 	}
  		} 
		// dividing corrected varied yields by default
		hV0PtSys[iSp][iMu][iSph][iSo][iVar]->Divide(hV0PtSys[iSp][iMu][iSph][iSo][deflt]);
		// apply RB sigmas as errors
  		for( Int_t i=1; i<hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetNbinsX()+1; i++){ 
    		hV0PtSys[iSp][iMu][iSph][iSo][iVar]->SetBinError(i, lSigmaDelta[i]);
  		}

	}	}	}	}	}


	// Make double-ratios to study correlations of syst. between HM_S0 and HM_minbias
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{	
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iVar==deflt) continue;
	if (iSp==1 && iSo==sysLifetime) continue;
	if (iMu == 0 && iSph > 0) continue;


		hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->Divide(hV0PtSys[iSp][iMu][iSph][iSo][iVar],hV0PtSys[iSp][iMu][0][iSo][iVar],1,1,"");
		
		// and calculate uncertainties
		for( Int_t i=1; i<hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->GetNbinsX()+1; i++){ 
			Double_t errS0 = hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinError(i);
			Double_t errHM = hV0PtSys[iSp][iMu][0][iSo][iVar]->GetBinError(i);
			Double_t errR = TMath::Sqrt(TMath::Abs(errS0*errS0 - errHM*errHM));
			errR = (hV0PtSys[iSp][iMu][0][iSo][iVar]->GetBinContent(i) > 0) ? errR / hV0PtSys[iSp][iMu][0][iSo][iVar]->GetBinContent(i) : 0;
			hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->SetBinError(i,errR);
		}

	}	}	}	}	}


	// NT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iVar==deflt) continue;
	if (iSp==1 && iSo==sysLifetime) continue;

		Double_t lSigmaDelta[100][100]; 
		for (int iNt = 1; iNt < hV0PtNtSys[iSp][iReg][iSo][iVar]->GetYaxis()->GetNbins(); ++iNt)	{
  			for( Int_t iPt=1; iPt<hV0PtNtSys[iSp][iReg][iSo][iVar]->GetNbinsX()+1; iPt++){ 
    			//Computation of roger barlow sigma_{delta} 
    			lSigmaDelta[iPt][iNt] = TMath::Sqrt( TMath::Abs( TMath::Power(hV0PtNtSys[iSp][iReg][iSo][iVar]->GetBinError(iPt,iNt),2) - TMath::Power(hV0PtNtSys[iSp][iReg][iSo][deflt]->GetBinError(iPt,iNt),2) ) );
 		   		//Computation of relationship to h2 for plotting in ratio plot 
    			if ( hV0PtNtSys[iSp][iReg][iSo][deflt]->GetBinContent(iPt,iNt) > 1e-12 )	{ 
    				lSigmaDelta[iPt][iNt] /= hV0PtNtSys[iSp][iReg][iSo][deflt]->GetBinContent(iPt,iNt);		}
    			else	{ 
    				lSigmaDelta[iPt][iNt] = 0; 	}
  			} 
		}

		// dividing corrected varied yields by default
		hV0PtNtSys[iSp][iReg][iSo][iVar]->Divide(hV0PtNtSys[iSp][iReg][iSo][deflt]);
		// apply RB sigmas as errors
  		for (int iNt = 1; iNt < hV0PtNtSys[iSp][iReg][iSo][iVar]->GetYaxis()->GetNbins(); ++iNt)	{
  			for( Int_t iPt=1; iPt<hV0PtNtSys[iSp][iReg][iSo][iVar]->GetNbinsX()+1; iPt++){ 
    			hV0PtNtSys[iSp][iReg][iSo][iVar]->SetBinError(iPt, iNt, lSigmaDelta[iPt][iNt]);
  			}
  		}

	}	}	}	}

	// RT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{
	if (iVar==deflt) continue;
	if (iSp==1 && iSo==sysLifetime) continue;


		Double_t lSigmaDelta[100]; 
		for( Int_t iPt=1; iPt<hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->GetNbinsX()+1; iPt++){ 
    		//Computation of roger barlow sigma_{delta} 
    		lSigmaDelta[iPt] = TMath::Sqrt( TMath::Abs( TMath::Power(hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->GetBinError(iPt),2) - TMath::Power(hV0PtRtSys[iSp][iReg][iSo][deflt][iRtBin]->GetBinError(iPt),2) ) );
 			//Computation of relationship to h2 for plotting in ratio plot 
    		if ( hV0PtRtSys[iSp][iReg][iSo][deflt][iRtBin]->GetBinContent(iPt) > 1e-12 )	{ 
    			lSigmaDelta[iPt] /= hV0PtRtSys[iSp][iReg][iSo][deflt][iRtBin]->GetBinContent(iPt);		}
    		else	{ 
    			lSigmaDelta[iPt] = 0; 	}
  		}

		// dividing corrected varied yields by default
		hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->Divide(hV0PtRtSys[iSp][iReg][iSo][deflt][iRtBin]);
		// apply RB sigmas as errors
  		for( Int_t iPt=1; iPt<hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->GetNbinsX()+1; iPt++){ 
    		hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->SetBinError(iPt, lSigmaDelta[iPt]);
  		}

	}	}	}	}	}

	// MAKE RATIOS TO STUDY CORRELATIONS OF UNCERTAINTIES
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
	if (iVar==deflt) continue;
	if (iSp==1 && iSo==sysLifetime) continue;
	
		for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

			hV0PtRtSysRatioToHM[iSp][iReg][iSo][iVar][iRtBin]->Divide(hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin],hV0PtRtSys[iSp][iReg][iSo][iVar][0],1,1,"");
		
			// and calculate uncertainties
			for( Int_t i=1; i<hV0PtRtSysRatioToHM[iSp][iReg][iSo][iVar][iRtBin]->GetNbinsX()+1; i++){ 
				Double_t err1 = hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->GetBinError(i);
				Double_t errHM = hV0PtRtSys[iSp][iReg][iSo][iVar][0]->GetBinError(i);
				Double_t errR = TMath::Sqrt(TMath::Abs(err1*err1 - errHM*errHM));
				errR = (hV0PtRtSys[iSp][iReg][iSo][iVar][0]->GetBinContent(i) > 0) ? errR / hV0PtRtSys[iSp][iReg][iSo][iVar][0]->GetBinContent(i) : 0;
				hV0PtRtSysRatioToHM[iSp][iReg][iSo][iVar][iRtBin]->SetBinError(i,errR);
			}

		}

	}	}	}	}


	return;
}

void MyAnalysisV0syst::CalculateAndPlotMaxDeviations() {

	mHandler->root()->SetBatch(kTRUE);
	Int_t colors[] = { kBlue, kGreen+3, kBlack, kCyan+1, kMagenta};

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		//creating canvases for each species and hm and s0 configuration
		TCanvas* cDev = new TCanvas(Form("cDev_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),"",2500,1250);
		cDev->Divide(4,3,5e-5,5e-5); 
		Int_t padC = 1;

		TCanvas* cDevRatioToHM = new TCanvas(Form("cDevRatioToHM_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]),"",2500,1250);
		cDevRatioToHM->Divide(4,3,5e-5,5e-5);


		double lowrange = 0.849;

		//each pad shows different sources
		for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
			if (iSp==1 && iSo==sysLifetime) continue;
			//cDev->cd(padC);//->SetLogx(kTRUE);
				
			TLegend* leg1; TLegend* leg2;
			hV0PtSys[iSp][iMu][iSph][iSo][0]->GetYaxis()->SetRangeUser(lowrange,1.0701);
			hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][0]->GetYaxis()->SetRangeUser(lowrange,1.0701);
			leg1 = new TLegend(0.21,0.69,0.46,0.88);
			mHandler->MakeNiceLegend(leg1, 0.052, 1.);
			leg1->AddEntry((TObject*)0,Form("%s %s %s,  %s",SPECNAMES[iSp],PLOTS_MULTI[iMu],PLOTS_SPHERO[iSph],PLOTS_SYSTS[iSo])," ");

			// plot different deviations
			cDev->cd(padC);
			for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
			if (iVar==deflt) continue;

				mHandler->MakeNiceHistogram(hV0PtSys[iSp][iMu][iSph][iSo][iVar],colors[iVar]);
				hV0PtSys[iSp][iMu][iSph][iSo][iVar]->SetMarkerSize(1.0);
				hV0PtSys[iSp][iMu][iSph][iSo][iVar]->SetLineColor(colors[iVar]);
				hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetYaxis()->SetTitleSize(0.035);
				hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetYaxis()->SetTitle("Rel. deviation");
				hV0PtSys[iSp][iMu][iSph][iSo][iVar]->Draw((iVar)?"h same":"h ");		
			}
			leg1->Draw();
			if (iSo == 1) {
				leg2 = new TLegend(0.49,0.3,0.88,0.50);
				mHandler->MakeNiceLegend(leg2, 0.045, 2.);
				leg2->AddEntry(hV0PtSys[iSp][iMu][iSph][2][0],"loosest","pl");
				leg2->AddEntry(hV0PtSys[iSp][iMu][iSph][2][1],"loose","pl");
				leg2->AddEntry(hV0PtSys[iSp][iMu][iSph][2][3],"tight","pl");
				leg2->AddEntry(hV0PtSys[iSp][iMu][iSph][2][4],"tightest","pl");
				leg2->Draw();
			}

			// and the same for correlations to HM
			cDevRatioToHM->cd(padC);
			for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
			if (iVar==deflt) continue;

				mHandler->MakeNiceHistogram(hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar],colors[iVar]);
				hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->SetMarkerSize(1.0);
				hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->SetLineColor(colors[iVar]);
				hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->GetYaxis()->SetTitleSize(0.035);
				hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->GetYaxis()->SetTitle("Rel. deviation, ratio to HM");
				hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->Draw((iVar)?"h same":"h ");		
			}
			leg1->Draw();
			if (iSo == 1)	leg2->Draw();

			// calculate max deviations
			for (int iBin = 1; iBin<hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetNbinsX()+1;iBin++) {
				Double_t maxD = 0;
				for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
				if (iVar==deflt) continue;
					Double_t varD = TMath::Abs(hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinContent(iBin)-1.);
					if (varD>maxD && varD>nRBsigmas*hV0PtSys[iSp][iMu][iSph][iSo][iVar]->GetBinError(iBin)) maxD=varD;
					//applying only if dev larger than 1 RB sigma
				}
				if (hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinContent(iBin)>0) hV0PtSysMaxD[iSp][iMu][iSph][iSo]->SetBinContent(iBin,maxD);
			}

			// and for correlations
			for (int iBin = 1; iBin<hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][deflt]->GetNbinsX()+1;iBin++) {
				Double_t maxD = 0;
				for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
				if (iVar==deflt) continue;
					Double_t varD = TMath::Abs(hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->GetBinContent(iBin)-1.);
					if (varD>maxD && varD>nRBsigmas*hV0PtSysRatioToHM[iSp][iMu][iSph][iSo][iVar]->GetBinError(iBin)) maxD=varD;
					//applying only if dev larger than 1 RB sigma
				}
				if (hV0PtSys[iSp][iMu][iSph][iSo][deflt]->GetBinContent(iBin)>0) hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->SetBinContent(iBin,maxD);
			}

			
			//now superimpose max deviations in the same pads
			mHandler->MakeNiceHistogram(hV0PtSysMaxD[iSp][iMu][iSph][iSo],kRed);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetYaxis()->SetRangeUser(-0.0005,3.*hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetMaximum());
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->SetFillColor(kRed);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->SetLineColor(kRed);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->SetLineWidth(2);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->SetFillStyle(3005);
						
			//scale hint1 to the pad coordinates
			cDev->cd(padC)->Update();
			Float_t rightmax = 2.0*hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetMaximum();
			Float_t scale = (rightmax>0) ? (gPad->GetUymax()-gPad->GetUymin())/rightmax : 1.;
			TH1F* htmpMax = (TH1F*)hV0PtSysMaxD[iSp][iMu][iSph][iSo]->Clone("tmpMax");
			htmpMax->Scale(scale);
			for (int i=1; i < htmpMax->GetNbinsX()+1; i++) htmpMax->SetBinContent(i,lowrange+htmpMax->GetBinContent(i));	
			htmpMax->Draw("same");
			 
			//draw an axis on the right side
			TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
	    	gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
			axis->SetLineColor(kRed);
			axis->SetLabelColor(kRed);
			axis->SetVertical(1);
			axis->SetTitleSize(0.035);
			axis->SetTitleOffset(1.5);
			axis->SetTitleColor(kRed);
			axis->SetTitle("Max. deviation");
			axis->Draw();


			// and the same for correlations to HM
			mHandler->MakeNiceHistogram(hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo],kRed);
			hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->GetYaxis()->SetRangeUser(-0.0005,3.*hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->GetMaximum());
			hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->SetFillColor(kRed);
			hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->SetLineColor(kRed);
			hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->SetLineWidth(2);
			hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->SetFillStyle(3005);

						
			//scale hint1 to the pad coordinates
			cDevRatioToHM->cd(padC)->Update();
			rightmax = 2.0*hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->GetMaximum();
			scale = (rightmax>0) ? (gPad->GetUymax()-gPad->GetUymin())/rightmax : 1.;
			TH1F* htmpMax2 = (TH1F*)hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->Clone("tmpMax2");
			htmpMax2->Scale(scale);
			for (int i=1; i < htmpMax2->GetNbinsX()+1; i++) htmpMax2->SetBinContent(i,lowrange+htmpMax2->GetBinContent(i));	
			htmpMax2->Draw("same");
			 
			//draw an axis on the right side
			TGaxis *axis2 = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
	    	gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
			axis2->SetLineColor(kRed);
			axis2->SetLabelColor(kRed);
			axis2->SetVertical(1);
			axis2->SetTitleSize(0.035);
			axis2->SetTitleOffset(1.5);
			axis2->SetTitleColor(kRed);
			axis2->SetTitle("Max. deviation");
			axis2->Draw();


			//cDevMax->cd(padC)->SetLogx(kTRUE);
			//hV0PtSysMaxD[iSp][iMu][iSph][iSo]->Scale(1./scale);
			//hV0PtSysMaxD[iSp][iMu][iSph][iSo]->Draw("hist");
			//leg1->Draw();

			/*TVirtualPad* c = cDev->GetPad(padC);
			TPad* ctop = (TPad*)c->Clone("ctop");
			c->Clear();
			ctop->SetBottomMargin(0.005);
			c->cd();

			TPad* p1 = new TPad(Form("p1_%s",c->GetName()),"",0.,0.3,1.,1.);
			p1->SetBottomMargin(0.);
			p1->Draw();
			p1->cd();
			ctop->DrawClone();
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->Draw("hist");

			c->cd();
			TPad* p2 = new TPad(Form("p2_%s",c->GetName()),"",0.,0.00,1.,0.28);
			p2->SetTopMargin(0);
			p2->SetBottomMargin(0.32);
			p2->Draw();
			p2->cd();
			
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetYaxis()->CenterTitle();
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetYaxis()->SetNdivisions(505);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetYaxis()->SetTitleSize(25);
			//hr->GetYaxis()->SetTitleFont(43);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetYaxis()->SetTitleOffset(1.55);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetYaxis()->SetLabelFont(43); 
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetYaxis()->SetLabelSize(20);

			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetXaxis()->SetTitleSize(25);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetXaxis()->SetTitleFont(43);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetXaxis()->SetTitleOffset(4.);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetXaxis()->SetLabelFont(43); 
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetXaxis()->SetLabelSize(25);
			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetXaxis()->SetTickLength(0.09);

			hV0PtSysMaxD[iSp][iMu][iSph][iSo]->Draw();
			//c->SetCanvasSize()
			c->cd();*/


			padC++;
		}

		return;

		cDev->Write();
		cDev->SaveAs(Form("plotsys/cDev_%s_%s_%s.png",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		cDevRatioToHM->Write();
		cDevRatioToHM->SaveAs(Form("plotsys/cDevRatioToHM_%s_%s_%s.png",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		//cDevMax->Write();
		//cDevMax->SaveAs(Form("plotsys/cDevMax_%s_%s_%s.png",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

	}	}	}

	mHandler->root()->SetBatch(kFALSE);
}

void MyAnalysisV0syst::CalculateMaxDeviationsNt() {

	// NT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	if (iSp==1 && iSo==sysLifetime) continue;

		for (int iNt = 1; iNt < hV0PtNtSys[iSp][iReg][iSo][deflt]->GetYaxis()->GetNbins(); ++iNt)	{
  			for( Int_t iPt=1; iPt<hV0PtNtSys[iSp][iReg][iSo][deflt]->GetNbinsX()+1; iPt++){ 
    	
    			hV0PtNtSysMaxD[iSp][iReg][iSo]->SetBinContent(iPt,iNt,0);
    			hV0PtNtSysMaxD[iSp][iReg][iSo]->SetBinError(iPt,iNt,0);

    			Double_t maxD = 0;
				for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
				if (iVar==deflt) continue;

					Double_t varD = TMath::Abs(hV0PtNtSys[iSp][iReg][iSo][iVar]->GetBinContent(iPt,iNt)-1.);
					if (varD>maxD && varD>nRBsigmas*hV0PtNtSys[iSp][iReg][iSo][iVar]->GetBinError(iPt,iNt)) maxD=varD;
					//applying only if dev larger than 1 RB sigma
				}
				if (hV0PtNtSys[iSp][iReg][iSo][deflt]->GetBinContent(iPt,iNt)>0) hV0PtNtSysMaxD[iSp][iReg][iSo]->SetBinContent(iPt,iNt,maxD);
			}
		}

	}	}	}

	// RT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{
	if (iSp==1 && iSo==sysLifetime) continue;

		for( Int_t iPt=1; iPt<hV0PtRtSys[iSp][iReg][iSo][deflt][iRtBin]->GetNbinsX()+1; iPt++){ 
    	
    		hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin]->SetBinContent(iPt,0);
			hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin]->SetBinError(iPt,0);

   			Double_t maxD = 0;
			for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
			if (iVar==deflt) continue;

				Double_t varD = TMath::Abs(hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->GetBinContent(iPt)-1.);
				if (varD>maxD && varD>nRBsigmas*hV0PtRtSys[iSp][iReg][iSo][iVar][iRtBin]->GetBinError(iPt)) maxD=varD;
					//applying only if dev larger than 1 RB sigma
			}
			if (hV0PtRtSys[iSp][iReg][iSo][deflt][iRtBin]->GetBinContent(iPt)>0) hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin]->SetBinContent(iPt,maxD);
	
		}

	}	}	}	}

	// RT ANALYSIS -- RATIOS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{
	if (iSp==1 && iSo==sysLifetime) continue;

		for( Int_t iPt=1; iPt<hV0PtRtSysRatioToHM[iSp][iReg][iSo][deflt][iRtBin]->GetNbinsX()+1; iPt++){ 
    	
    		hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin]->SetBinContent(iPt,0);
			hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin]->SetBinError(iPt,0);

   			Double_t maxD = 0;
			for (int iVar = 0; iVar < sysVarsSizeof; ++iVar)	{
			if (iVar==deflt) continue;

				Double_t varD = TMath::Abs(hV0PtRtSysRatioToHM[iSp][iReg][iSo][iVar][iRtBin]->GetBinContent(iPt)-1.);
				if (varD>maxD && varD>nRBsigmas*hV0PtRtSysRatioToHM[iSp][iReg][iSo][iVar][iRtBin]->GetBinError(iPt)) maxD=varD;
					//applying only if dev larger than 1 RB sigma
			}
			if (hV0PtRtSys[iSp][iReg][iSo][deflt][iRtBin]->GetBinContent(iPt)>0) hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin]->SetBinContent(iPt,maxD);
	
		}

	}	}	}	}

}

void MyAnalysisV0syst::CalculateSignalExSys() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		// CALCULATE YIELDS WITH DIFFERENT SB SIGMAS
		Double_t* yield = 0;
		for (int iBin = 0; iBin < NPTBINS+1; ++iBin)	{
			cout << "Extracting pt bin " << iBin << " of " << hV0IMvPtSys[iSp][iMu][iSph][0][deflt]->GetName() << endl;
			yield = ((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSBVarySigma(7.,(TH1F*)hV0IMvPtSys[iSp][iMu][iSph][0][deflt]->ProjectionY(
					Form("iSp%i_iBin%i", iSp, iBin), iBin,iBin));
				hV0PtSysSigExLoose[iSp][iMu][iSph]->SetBinContent(iBin,*(yield+0));
				hV0PtSysSigExLoose[iSp][iMu][iSph]->SetBinError(iBin,*(yield+1));

			yield = ((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSBVarySigma(5.,(TH1F*)hV0IMvPtSys[iSp][iMu][iSph][0][deflt]->ProjectionY(
					Form("iSp%i_iBin%i", iSp, iBin), iBin,iBin));
				hV0PtSysSigExTight[iSp][iMu][iSph]->SetBinContent(iBin,*(yield+0));
				hV0PtSysSigExTight[iSp][iMu][iSph]->SetBinError(iBin,*(yield+1));
		}

		hV0PtSysSigExLoose[iSp][iMu][iSph]->Scale(1,"width");
		hV0PtSysSigExTight[iSp][iMu][iSph]->Scale(1,"width");
		hV0PtSysSigExLoose[iSp][iMu][iSph]->Divide(hV0EfficiencySys[iSp][0][deflt]);
		hV0PtSysSigExTight[iSp][iMu][iSph]->Divide(hV0EfficiencySys[iSp][0][deflt]);



		// CALCULATE RATIOS W RB SIGMAS
		Double_t lSigmaDelta[100]; 
  		for( Int_t i=1; i<hV0PtSysSigExLoose[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
    		//Computation of roger barlow sigma_{delta} 
    		lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(hV0PtSysSigExLoose[iSp][iMu][iSph]->GetBinError(i),2) - TMath::Power(hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinError(i),2) ) );
 		   //Computation of relationship to h2 for plotting in ratio plot 
    		if ( hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinContent(i) > 1e-12 )	{ 
    			lSigmaDelta[i] /= hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinContent(i);		}
    		else	{ 
    			lSigmaDelta[i] = 0; 	}
  		} 
		// dividing corrected varied yields by default
		hV0PtSysSigExLoose[iSp][iMu][iSph]->Divide(hV0PtSys[iSp][iMu][iSph][0][deflt]);
		// apply RB sigmas as errors
  		for( Int_t i=1; i<hV0PtSysSigExLoose[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
    		hV0PtSysSigExLoose[iSp][iMu][iSph]->SetBinError(i, lSigmaDelta[i]);
  		}
  		for( Int_t i=1; i<hV0PtSysSigExTight[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
    		//Computation of roger barlow sigma_{delta} 
    		lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(hV0PtSysSigExTight[iSp][iMu][iSph]->GetBinError(i),2) - TMath::Power(hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinError(i),2) ) );
 		   //Computation of relationship to h2 for plotting in ratio plot 
    		if ( hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinContent(i) > 1e-12 )	{ 
    			lSigmaDelta[i] /= hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinContent(i);		}
    		else	{ 
    			lSigmaDelta[i] = 0; 	}
  		} 
		// dividing corrected varied yields by default
		hV0PtSysSigExTight[iSp][iMu][iSph]->Divide(hV0PtSys[iSp][iMu][iSph][0][deflt]);
		// apply RB sigmas as errors
  		for( Int_t i=1; i<hV0PtSysSigExTight[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
    		hV0PtSysSigExTight[iSp][iMu][iSph]->SetBinError(i, lSigmaDelta[i]);
  		}


  		// CALCULATE MAX DEVIATIONS
  		for (int iBin = 1; iBin<hV0PtSysSigExLoose[iSp][iMu][iSph]->GetNbinsX()+1;iBin++) {
				Double_t maxD = 0;
				Double_t varD = TMath::Abs(hV0PtSysSigExLoose[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
				if (varD>maxD && varD>nRBsigmas*hV0PtSysSigExLoose[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
					//applying only if dev larger than 1 RB sigma

				varD = TMath::Abs(hV0PtSysSigExTight[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
				if (varD>maxD && varD>nRBsigmas*hV0PtSysSigExTight[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
				
				if (hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinContent(iBin)>0) hV0PtSysMaxDSigEx[iSp][iMu][iSph]->SetBinContent(iBin,maxD);
		}

		Int_t colors[] = { kBlue, kGreen+2, kBlack, kRed, kMagenta};
		mHandler->MakeNiceHistogram(hV0PtSysSigExLoose[iSp][iMu][iSph],colors[loose]);
		hV0PtSysSigExLoose[iSp][iMu][iSph]->SetMarkerSize(0.8);
		hV0PtSysSigExLoose[iSp][iMu][iSph]->SetLineColor(colors[loose]);
		mHandler->MakeNiceHistogram(hV0PtSysSigExTight[iSp][iMu][iSph],colors[tight]);
		hV0PtSysSigExTight[iSp][iMu][iSph]->SetMarkerSize(0.8);
		hV0PtSysSigExTight[iSp][iMu][iSph]->SetLineColor(colors[tight]);

		mHandler->MakeNiceHistogram(hV0PtSysMaxDSigEx[iSp][iMu][iSph],kBlack);
		hV0PtSysMaxDSigEx[iSp][iMu][iSph]->SetFillStyle(3002);
		hV0PtSysMaxDSigEx[iSp][iMu][iSph]->SetFillColor(kBlack);
		hV0PtSysMaxDSigEx[iSp][iMu][iSph]->GetYaxis()->SetRangeUser(-0.0005,0.2);
		hV0PtSysMaxDSigEx[iSp][iMu][iSph]->SetTitle(Form("%s %s %s; V0 p_{T} (GeV/#it{c}); Relative syst. uncertainty",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

	}	}	}

	// and study the effect on ratio S0/HM
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]->Divide(hV0PtSysSigExTight[iSp][iMu][iSph],hV0PtSysSigExTight[iSp][iMu][0],1,1,"");
		hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->Divide(hV0PtSysSigExLoose[iSp][iMu][iSph],hV0PtSysSigExLoose[iSp][iMu][0],1,1,"");
		
		// and calculate uncertainties
		for( Int_t i=1; i<hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
			Double_t errS0 = hV0PtSysSigExTight[iSp][iMu][iSph]->GetBinError(i);
			Double_t errHM = hV0PtSysSigExTight[iSp][iMu][0]->GetBinError(i);
			Double_t errR = TMath::Sqrt(TMath::Abs(errS0*errS0 - errHM*errHM));
			errR = (hV0PtSysSigExTight[iSp][iMu][0]->GetBinContent(i) > 0) ? errR / hV0PtSysSigExTight[iSp][iMu][0]->GetBinContent(i) : 0;
			hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]->SetBinError(i,errR);
		}
		for( Int_t i=1; i<hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
			Double_t errS0 = hV0PtSysSigExLoose[iSp][iMu][iSph]->GetBinError(i);
			Double_t errHM = hV0PtSysSigExLoose[iSp][iMu][0]->GetBinError(i);
			Double_t errR = TMath::Sqrt(TMath::Abs(errS0*errS0 - errHM*errHM));
			errR = (hV0PtSysSigExLoose[iSp][iMu][0]->GetBinContent(i) > 0) ? errR / hV0PtSysSigExLoose[iSp][iMu][0]->GetBinContent(i) : 0;
			hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->SetBinError(i,errR);
		}

		// calculate max deviations
  		for (int iBin = 1; iBin<hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->GetNbinsX()+1;iBin++) {
				Double_t maxD = 0;
				Double_t varD = TMath::Abs(hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
				if (varD>maxD && varD>nRBsigmas*hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
					//applying only if dev larger than 1 RB sigma

				varD = TMath::Abs(hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
				if (varD>maxD && varD>nRBsigmas*hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
				
				if (hV0PtSys[iSp][iMu][iSph][0][deflt]->GetBinContent(iBin)>0) hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,maxD);
		}

		Int_t colors[] = { kBlue, kGreen+2, kBlack, kRed, kMagenta};
		mHandler->MakeNiceHistogram(hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph],colors[loose]);
		hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->SetMarkerSize(0.8);
		hV0PtSysSigExLooseRatioToHM[iSp][iMu][iSph]->SetLineColor(colors[loose]);
		mHandler->MakeNiceHistogram(hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph],colors[tight]);
		hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]->SetMarkerSize(0.8);
		hV0PtSysSigExTightRatioToHM[iSp][iMu][iSph]->SetLineColor(colors[tight]);

		mHandler->MakeNiceHistogram(hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph],kBlack);
		hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->SetFillStyle(3002);
		hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->SetFillColor(kBlack);
		hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->GetYaxis()->SetRangeUser(-0.0005,0.2);
		hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->SetTitle(Form("%s %s %s; V0 p_{T} (GeV/#it{c}); Relative syst. uncertainty",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));


	}	}	}	

	mHandler->root()->SetBatch(kTRUE);
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	if (NMULTI<5) continue;
	if (NSPHERO<5) continue;
		TCanvas* cSig = new TCanvas(Form("cSig_%s",SPECIES[iSp]),"",1400,1000);
		cSig->Divide(3,2,5e-5,5e-5); 
		Int_t padC = 1;

		cout << hV0PtSysMaxDSigEx[iSp][3][0] << endl;
		cSig->cd(1)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][3][0]);
		cout << hV0PtSysMaxDSigEx[iSp][3][3] << endl;
		cSig->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][3][3]);
		cout << hV0PtSysMaxDSigEx[iSp][3][4] << endl;
		cSig->cd(3)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][3][4]);
		cout << hV0PtSysMaxDSigEx[iSp][4][0] << endl;
		cSig->cd(4)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][4][0]);
		cout << hV0PtSysMaxDSigEx[iSp][4][3] << endl;
		cSig->cd(5)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][4][3]);
		cout << hV0PtSysMaxDSigEx[iSp][4][4] << endl;
		cSig->cd(6)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][4][4]);
		cSig->Write();
		cSig->SaveAs(Form("plotsys/cSig_%s.png",SPECIES[iSp]));
	}
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	if (NMULTI<2) continue;
	if (NSPHERO<3) continue;
		TCanvas* cSigRatioToHM = new TCanvas(Form("cSigRatioToHM_%s",SPECIES[iSp]),"",1500,450);
		cSigRatioToHM->Divide(3,1,5e-5,5e-5); 
		Int_t padC = 1;

		cSigRatioToHM->cd(1)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigExRatioToHM[iSp][1][0]);	
		cSigRatioToHM->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigExRatioToHM[iSp][1][1]);
		cSigRatioToHM->cd(3)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigExRatioToHM[iSp][1][2]);
		cSigRatioToHM->Write();
		cSigRatioToHM->SaveAs(Form("plotsys/cSigRatioToHM_%s.png",SPECIES[iSp]));
	}
	mHandler->root()->SetBatch(kFALSE);

}

void MyAnalysisV0syst::CalculateSignalExSysNt() {

	// RT ANALYSIS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

		// RTREC BINS
		int lowedge = GetRTBinRec(iRtBin,kTRUE);
		int upedge  = GetRTBinRec(iRtBin,kFALSE);
		Double_t NormScaleNt = 1.;//hNtUnf->Integral(lowedge,upedge)/hNtUnf->Integral(1,51);

		hV0IMvPtNtSys[iSp][iReg][0][deflt]->GetZaxis()->SetRange(lowedge,upedge);
		TH2F* hXY = (TH2F*)hV0IMvPtNtSys[iSp][iReg][0][deflt]->Project3D("yx");

		// CALCULATE YIELDS WITH DIFFERENT SB SIGMAS
		Double_t* yield = 0;
		for (int iBin = 1; iBin < NPTBINS+1; ++iBin)	{

			//cout << "Extracting pt bin " << iBin << " of " << hV0IMvPtNtSys[iSp][iReg][0][deflt]->GetName() << endl;
			yield = ((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSBVarySigma(7.,
					(TH1F*)hXY->ProjectionY(
					Form("iSp%i_iBin%i", iSp, iBin), iBin,iBin));
			hV0PtRtSysSigExLoose[iSp][iReg][iRtBin]->SetBinContent(iBin,*(yield+0));
			hV0PtRtSysSigExLoose[iSp][iReg][iRtBin]->SetBinError(iBin,*(yield+1));

			yield = ((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSBVarySigma(6.,
					(TH1F*)hXY->ProjectionY(
					Form("iSp%i_iBin%i", iSp, iBin), iBin,iBin));
			hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]->SetBinContent(iBin,*(yield+0));
			hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]->SetBinError(iBin,*(yield+1));


			yield = ((MyAnalysisV0extract*)mHandler->analysis(1))->ExtractYieldSBVarySigma(5.,
					(TH1F*)hXY->ProjectionY(
					Form("iSp%i_iBin%i", iSp, iBin), iBin,iBin));
			hV0PtRtSysSigExTight[iSp][iReg][iRtBin]->SetBinContent(iBin,*(yield+0));
			hV0PtRtSysSigExTight[iSp][iReg][iRtBin]->SetBinError(iBin,*(yield+1));
		}

		delete hXY;

		/*hV0PtRtSysSigExLoose[iSp][iReg][iRtBin]->Scale(1,"width");
		hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]->Scale(1,"width");
		hV0PtRtSysSigExTight[iSp][iReg][iRtBin]->Scale(1,"width");
		hV0PtRtSysSigExLoose[iSp][iReg][iRtBin]->Divide(hV0EfficiencySys[iSp][0][deflt]);
		hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]->Divide(hV0EfficiencySys[iSp][0][deflt]);
		hV0PtRtSysSigExTight[iSp][iReg][iRtBin]->Divide(hV0EfficiencySys[iSp][0][deflt]);*/

		mDirFile->cd();
		// CALCULATE RATIOS W RB SIGMAS -- LOOSE
		hV0PtRtSysSigExLoose[iSp][iReg][iRtBin] 
			= CalculateRatiosRBSigmas(hV0PtRtSysSigExLoose[iSp][iReg][iRtBin],hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]);
		hV0PtRtSysSigExTight[iSp][iReg][iRtBin] 
			= CalculateRatiosRBSigmas(hV0PtRtSysSigExTight[iSp][iReg][iRtBin],hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]);


  		// CALCULATE MAX DEVIATIONS
		hV0PtRtSysMaxDSigEx[iSp][iReg][iRtBin]
			= CalculateSigExMaxD(hV0PtRtSysMaxDSigEx[iSp][iReg][iRtBin],
				hV0PtRtSysSigExLoose[iSp][iReg][iRtBin],
				hV0PtRtSysSigExTight[iSp][iReg][iRtBin],
				hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]);
  		
	}	}	}


	// RT ANALYSIS -- RATIOS
	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

		hV0PtRtSysSigExLooseRatioToHM[iSp][iReg][iRtBin]->Divide(hV0PtRtSysSigExLoose[iSp][iReg][iRtBin],hV0PtRtSysSigExLoose[iSp][iReg][0],1,1,"");
		hV0PtRtSysSigExTightRatioToHM[iSp][iReg][iRtBin]->Divide(hV0PtRtSysSigExTight[iSp][iReg][iRtBin],hV0PtRtSysSigExTight[iSp][iReg][0],1,1,"");
		hV0PtRtSysSigExDefltRatioToHM[iSp][iReg][iRtBin]->Divide(hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin],hV0PtRtSysSigExDeflt[iSp][iReg][0],1,1,"");
		
		// and calculate uncertainties
		for( Int_t i=1; i<hV0PtRtSysSigExDefltRatioToHM[iSp][iReg][iRtBin]->GetNbinsX()+1; i++){ 
			Double_t err1 = hV0PtRtSysSigExLoose[iSp][iReg][iRtBin]->GetBinError(i);
			Double_t errHM = hV0PtRtSysSigExLoose[iSp][iReg][0]->GetBinError(i);
			Double_t errR = TMath::Sqrt(TMath::Abs(err1*err1 - errHM*errHM));
			errR = (hV0PtRtSysSigExLoose[iSp][iReg][0]->GetBinContent(i) > 0) ? errR / hV0PtRtSysSigExLoose[iSp][iReg][0]->GetBinContent(i) : 0;
			hV0PtRtSysSigExLooseRatioToHM[iSp][iReg][iRtBin]->SetBinError(i,errR);

			err1 = hV0PtRtSysSigExTight[iSp][iReg][iRtBin]->GetBinError(i);
			errHM = hV0PtRtSysSigExTight[iSp][iReg][0]->GetBinError(i);
			errR = TMath::Sqrt(TMath::Abs(err1*err1 - errHM*errHM));
			errR = (hV0PtRtSysSigExTight[iSp][iReg][0]->GetBinContent(i) > 0) ? errR / hV0PtRtSysSigExTight[iSp][iReg][0]->GetBinContent(i) : 0;
			hV0PtRtSysSigExTightRatioToHM[iSp][iReg][iRtBin]->SetBinError(i,errR);

			err1 = hV0PtRtSysSigExDeflt[iSp][iReg][iRtBin]->GetBinError(i);
			errHM = hV0PtRtSysSigExDeflt[iSp][iReg][0]->GetBinError(i);
			errR = TMath::Sqrt(TMath::Abs(err1*err1 - errHM*errHM));
			errR = (hV0PtRtSysSigExDeflt[iSp][iReg][0]->GetBinContent(i) > 0) ? errR / hV0PtRtSysSigExDeflt[iSp][iReg][0]->GetBinContent(i) : 0;
			hV0PtRtSysSigExDefltRatioToHM[iSp][iReg][iRtBin]->SetBinError(i,errR);
		}

		// calculate max deviations
		hV0PtRtSysMaxDSigExRatioToHM[iSp][iReg][iRtBin]
			= CalculateSigExMaxD(hV0PtRtSysMaxDSigExRatioToHM[iSp][iReg][iRtBin],
				hV0PtRtSysSigExLooseRatioToHM[iSp][iReg][iRtBin],
				hV0PtRtSysSigExTightRatioToHM[iSp][iReg][iRtBin],
				hV0PtRtSysSigExDefltRatioToHM[iSp][iReg][iRtBin]);


	}	}	}

}

TH1F* MyAnalysisV0syst::CalculateRatiosRBSigmas(TH1F* ha, TH1F* hd) {

	Double_t lSigmaDelta[100]; 
  	for( Int_t i=1; i<ha->GetNbinsX()+1; i++){ 
    		//Computation of roger barlow sigma_{delta} 
    		lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(ha->GetBinError(i),2) - TMath::Power(hd->GetBinError(i),2) ) );
 		   //Computation of relationship to h2 for plotting in ratio plot 
    		if ( hd->GetBinContent(i) > 1e-12 )	{ 
    			lSigmaDelta[i] /= hd->GetBinContent(i);		}
    		else	{ 
    			lSigmaDelta[i] = 0; 	}
  		} 
		// dividing corrected varied yields by default
		ha->Divide(hd);
		// apply RB sigmas as errors
  		for( Int_t i=1; i<ha->GetNbinsX()+1; i++){ 
    		ha->SetBinError(i, lSigmaDelta[i]);
  		}

  	return ha;
}

TH1F* MyAnalysisV0syst::CalculateSigExMaxD(TH1F* hm, TH1F* h1, TH1F* h2, TH1F* hd) {

	for (int iBin = 1; iBin<h1->GetNbinsX()+1;iBin++) {
		Double_t maxD = 0;
		Double_t varD = TMath::Abs(h1->GetBinContent(iBin)-1.);
		//applying only if dev larger than 1 RB sigma
		if (varD>maxD && varD>nRBsigmas*h1->GetBinError(iBin)) maxD=varD;
		
		varD = TMath::Abs(h2->GetBinContent(iBin)-1.);
		if (varD>maxD && varD>nRBsigmas*h2->GetBinError(iBin)) maxD=varD;
				
		if (hd->GetBinContent(iBin)>0) hm->SetBinContent(iBin,maxD);
		}

  	return hm;
}


void MyAnalysisV0syst::AddDeviations() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		//cout << "Calculating relative systematic uncertainties from cut variations for " << hV0PtSysSum[iSp][iMu][iSph] << " " << hV0PtSysSum[iSp][iMu][iSph]->GetName() << endl;
		for (int iBin = 1; iBin < hV0PtSysSum[iSp][iMu][iSph]->GetNbinsX()+1; iBin++) {
			
			Double_t sum = 0;
			Double_t sumCuts = 0;
			Double_t sumUnc = 0;
			Double_t sumCutsUnc = 0;

			// sources from cut variations
			for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
			if (iSp==1 && iSo==sysLifetime) continue;
				//cout << hV0PtSysMaxD[iSp][iMu][iSph][iSo] << endl;
				//cout << hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetName() << endl;
					
				sum 	+= hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetBinContent(iBin)*hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetBinContent(iBin);
				sumCuts += hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetBinContent(iBin)*hV0PtSysMaxD[iSp][iMu][iSph][iSo]->GetBinContent(iBin);

				sumUnc 	+= hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->GetBinContent(iBin)*hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->GetBinContent(iBin);
				sumCutsUnc 	+= hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->GetBinContent(iBin)*hV0PtSysMaxDRatioToHM[iSp][iMu][iSph][iSo]->GetBinContent(iBin);
				// sum unc denote the uncorrelated uncertainty on ratios to HM, e.g unc of f = S0-biased/HM
			}

			// from material budget, fully correlated between s0
			sum += hV0PtSysBudget->GetBinContent(iBin)*hV0PtSysBudget->GetBinContent(iBin);
			
			// from using MB efficiency
			if (sum>0) {
				sum += hV0PtSysEffi->GetBinContent(iBin)*hV0PtSysEffi->GetBinContent(iBin);
				sumUnc += 1.*hV0PtSysEffi->GetBinContent(iBin)*hV0PtSysEffi->GetBinContent(iBin);
				// factor 2. comes from propagation on ratios
				// for f = A/B, (s_f/f)^2 = (s_A/A)^2 + (s_B/B)^2
				// but (s_A/A) == (s_B/B)
				//// decided to use more liberal 1.0 similarly to other ALICE analyses
			}

			// from experimental bias due to not unfolding s0 distribution
			TH1F* hexpbias = (iSph%2) ? hV0PtSysExpBiasJetty : hV0PtSysExpBiasIso;
			if (sum>0 && iSph>0) {
				sum += hexpbias->GetBinContent(iBin)*hexpbias->GetBinContent(iBin);
				sumUnc += 1.*hexpbias->GetBinContent(iBin)*hexpbias->GetBinContent(iBin);
			}	

			// from signal extraction
			if (sum>0) {
				sum += hV0PtSysMaxDSigEx[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDSigEx[iSp][iMu][iSph]->GetBinContent(iBin);
				sumUnc += hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin);
			}

			// from feeddown correction
			if (iSp>1 && sum>0) {
				sum += hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->GetBinContent(iBin);
				sumUnc += hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin);
			}

			


			hV0PtSysSum[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(sum));
			hV0PtSysSum[iSp][iMu][iSph]->SetFillStyle(3005);
			hV0PtSysSum[iSp][iMu][iSph]->SetFillColor(kBlack);
			hV0PtSysSum[iSp][iMu][iSph]->SetTitle(Form("%s %s %s; V0 p_{T} (GeV/#it{c}); Relative syst. uncertainty",SPECNAMES[iSp],PLOTS_MULTI[iMu],PLOTS_SPHERO[iSph]));
			hV0PtSysSum[iSp][iMu][iSph]->GetYaxis()->SetTitleOffset(1.52);
			hV0PtSysSum[iSp][iMu][iSph]->SetStats(0);

			hV0PtSysSumUnc[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(sumUnc));
			hV0PtSysSumUnc[iSp][iMu][iSph]->SetFillStyle(3005);
			hV0PtSysSumUnc[iSp][iMu][iSph]->SetFillColor(kBlack);
			hV0PtSysSumUnc[iSp][iMu][iSph]->SetTitle(Form("%s %s %s; V0 p_{T} (GeV/#it{c}); Unc. syst. uncertainty (w.r.t. S_{0})",SPECNAMES[iSp],PLOTS_MULTI[iMu],PLOTS_SPHERO[iSph]));
			hV0PtSysSumUnc[iSp][iMu][iSph]->GetYaxis()->SetTitleOffset(1.52);
			hV0PtSysSumUnc[iSp][iMu][iSph]->SetStats(0);

			Double_t part = 0;
			if (sum>0) {
				part += hV0PtSysBudget->GetBinContent(iBin)*hV0PtSysBudget->GetBinContent(iBin);
				hFracBudget[iSp][iMu][iSph]->SetBinContent(iBin,part/TMath::Sqrt(sum));
				
				part += hV0PtSysEffi->GetBinContent(iBin)*hV0PtSysEffi->GetBinContent(iBin);
				hFracEffi[iSp][iMu][iSph]->SetBinContent(iBin,part/TMath::Sqrt(sum));
				
				if (iSph>0) part += hexpbias->GetBinContent(iBin)*hexpbias->GetBinContent(iBin);
				hFracExpBias[iSp][iMu][iSph]->SetBinContent(iBin,part/TMath::Sqrt(sum));
				
				part += sumCuts;
				hFracCuts[iSp][iMu][iSph]->SetBinContent(iBin,part/TMath::Sqrt(sum));

				part += hV0PtSysMaxDSigEx[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDSigEx[iSp][iMu][iSph]->GetBinContent(iBin);
				hFracSigEx[iSp][iMu][iSph]->SetBinContent(iBin,part/TMath::Sqrt(sum));

				if (iSp>1) {
					part += hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->GetBinContent(iBin);
					hFracFD[iSp][iMu][iSph]->SetBinContent(iBin,part/TMath::Sqrt(sum));
				}
			}

			// and ratios
			Double_t partR = 0;
			if (sumUnc>0) {
				
				partR += 1.*hV0PtSysEffi->GetBinContent(iBin)*hV0PtSysEffi->GetBinContent(iBin);
				hFracEffiRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,partR/TMath::Sqrt(sumUnc));
				
				if (iSph>0) partR += hexpbias->GetBinContent(iBin)*hexpbias->GetBinContent(iBin);
				hFracExpBiasRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,partR/TMath::Sqrt(sumUnc));
				
				partR += sumCutsUnc;
				hFracCutsRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,partR/TMath::Sqrt(sumUnc));

				partR += hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDSigExRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin);
				hFracSigExRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,partR/TMath::Sqrt(sumUnc));

				if (iSp>1) {
					partR += hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin);
					hFracFDRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,partR/TMath::Sqrt(sumUnc));
				}
			}
			
		}

		hFracBudget[iSp][iMu][iSph]->SetLineColor(kMagenta); hFracBudget[iSp][iMu][iSph]->SetLineWidth(3);
		hFracEffi[iSp][iMu][iSph]->SetLineColor(kGreen+2); hFracEffi[iSp][iMu][iSph]->SetLineWidth(3);
		hFracExpBias[iSp][iMu][iSph]->SetLineColor(kCyan+1); hFracExpBias[iSp][iMu][iSph]->SetLineWidth(3);
		hFracCuts[iSp][iMu][iSph]->SetLineColor(kBlue); hFracCuts[iSp][iMu][iSph]->SetLineWidth(3);
		hFracSigEx[iSp][iMu][iSph]->SetLineColor(kRed); hFracSigEx[iSp][iMu][iSph]->SetLineWidth(3);
		hFracFD[iSp][iMu][iSph]->SetLineColor(kBlack); hFracFD[iSp][iMu][iSph]->SetLineWidth(3);

		hFracEffiRatioToHM[iSp][iMu][iSph]->SetLineColor(kGreen+2); hFracEffiRatioToHM[iSp][iMu][iSph]->SetLineWidth(3);
		hFracExpBiasRatioToHM[iSp][iMu][iSph]->SetLineColor(kCyan+1); hFracExpBiasRatioToHM[iSp][iMu][iSph]->SetLineWidth(3);
		hFracCutsRatioToHM[iSp][iMu][iSph]->SetLineColor(kBlue); hFracCutsRatioToHM[iSp][iMu][iSph]->SetLineWidth(3);
		hFracSigExRatioToHM[iSp][iMu][iSph]->SetLineColor(kRed); hFracSigExRatioToHM[iSp][iMu][iSph]->SetLineWidth(3);
		hFracFDRatioToHM[iSp][iMu][iSph]->SetLineColor(kBlack); hFracFDRatioToHM[iSp][iMu][iSph]->SetLineWidth(3);


	}	}	}

	return;

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		if (NMULTI<4) continue;
		if (NSPHERO<5) continue;
		TCanvas* cTot = new TCanvas(Form("cTot_%s",SPECIES[iSp]),"",1400,1000);
		cTot->Divide(3,2,5e-5,5e-5); 
		Int_t padC = 1;

		TLegend* leg = new TLegend(0.13,.68,0.88,0.88);
		mHandler->MakeNiceLegend(leg,0.027,2);	
		leg->AddEntry(hFracBudget[iSp][3][0],"mat. budget","l");
		leg->AddEntry(hFracEffi[iSp][3][0],"+ MB effi.","l");
		leg->AddEntry(hFracExpBias[iSp][3][0],"+ exp. bias","l");
		leg->AddEntry(hFracCuts[iSp][3][0],"+ cuts var.","l");
		leg->AddEntry(hFracSigEx[iSp][3][0],"+ signal extr.","l");
		if (iSp>1) leg->AddEntry(hFracFD[iSp][3][0],"+ feeddown corr.","l");


		cTot->cd(1)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSum[iSp][3][0]);
		DrawMirrored(hFracBudget[iSp][3][0],"same ][");DrawMirrored(hFracEffi[iSp][3][0],"same ][");DrawMirrored(hFracCuts[iSp][3][0],"same ][");DrawMirrored(hFracSigEx[iSp][3][0],"same ][");if(iSp>1)DrawMirrored(hFracFD[iSp][3][0],"same ][");
		
		leg->Draw();


		cTot->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSum[iSp][3][3]);
		DrawMirrored(hFracBudget[iSp][3][3],"same ][");DrawMirrored(hFracEffi[iSp][3][3],"same ][");DrawMirrored(hFracExpBias[iSp][3][3],"same ][");DrawMirrored(hFracCuts[iSp][3][3],"same ][");DrawMirrored(hFracSigEx[iSp][3][3],"same ][");if(iSp>1)DrawMirrored(hFracFD[iSp][3][3],"same ][");
		cTot->cd(3)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSum[iSp][3][4]);
		DrawMirrored(hFracBudget[iSp][3][4],"same ][");DrawMirrored(hFracEffi[iSp][3][4],"same ][");DrawMirrored(hFracExpBias[iSp][3][4],"same ][");DrawMirrored(hFracCuts[iSp][3][4],"same ][");DrawMirrored(hFracSigEx[iSp][3][4],"same ][");if(iSp>1)DrawMirrored(hFracFD[iSp][3][4],"same ][");
		cTot->cd(4)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSum[iSp][4][0]);
		DrawMirrored(hFracBudget[iSp][4][0],"same ][");DrawMirrored(hFracEffi[iSp][4][0],"same ][");DrawMirrored(hFracCuts[iSp][4][0],"same ][");DrawMirrored(hFracSigEx[iSp][4][0],"same ][");if(iSp>1)DrawMirrored(hFracFD[iSp][4][0],"same ][");
		cTot->cd(5)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSum[iSp][4][3]);
		DrawMirrored(hFracBudget[iSp][4][3],"same ][");DrawMirrored(hFracEffi[iSp][4][3],"same ][");DrawMirrored(hFracExpBias[iSp][4][3],"same ][");DrawMirrored(hFracCuts[iSp][4][3],"same ][");DrawMirrored(hFracSigEx[iSp][4][3],"same ][");if(iSp>1)DrawMirrored(hFracFD[iSp][4][3],"same ][");
		cTot->cd(6)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSum[iSp][4][4]);
		DrawMirrored(hFracBudget[iSp][4][4],"same ][");DrawMirrored(hFracEffi[iSp][4][4],"same ][");DrawMirrored(hFracExpBias[iSp][4][4],"same ][");DrawMirrored(hFracCuts[iSp][4][4],"same ][");DrawMirrored(hFracSigEx[iSp][4][4],"same ][");if(iSp>1)DrawMirrored(hFracFD[iSp][4][4],"same ][");
		cTot->Write();
		cTot->SaveAs(Form("plotsys/cTot_%s.png",SPECIES[iSp]));
	}

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		if (NMULTI<2) continue;
		if (NSPHERO<3) continue;
		TCanvas* cTotRatioToHM = new TCanvas(Form("cTotRatioToHM_%s",SPECIES[iSp]),"",1500,450);
		cTotRatioToHM->Divide(3,1,5e-5,5e-5); 
		Int_t padC = 1;

		TLegend* leg = new TLegend(0.13,.68,0.88,0.88);
		mHandler->MakeNiceLegend(leg,0.027,2);	
		
		leg->AddEntry(hFracEffiRatioToHM[iSp][1][0],"MB effi.","l");
		leg->AddEntry(hFracExpBiasRatioToHM[iSp][1][0],"+ exp. bias","l");
		leg->AddEntry(hFracCutsRatioToHM[iSp][1][0],"+ cuts var.","l");
		leg->AddEntry(hFracSigExRatioToHM[iSp][1][0],"+ signal extr.","l");
		if (iSp>1) leg->AddEntry(hFracFDRatioToHM[iSp][1][0],"+ feeddown corr.","l");


		cTotRatioToHM->cd(1)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSumUnc[iSp][1][0]);
		DrawMirrored(hFracEffiRatioToHM[iSp][1][0],"same ][");DrawMirrored(hFracExpBiasRatioToHM[iSp][1][0],"same ][");DrawMirrored(hFracCutsRatioToHM[iSp][1][0],"same ][");DrawMirrored(hFracSigExRatioToHM[iSp][1][0],"same ][");if(iSp>1)DrawMirrored(hFracFDRatioToHM[iSp][1][0],"same ][");
		
		leg->Draw();

		cTotRatioToHM->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSumUnc[iSp][1][1]);
		DrawMirrored(hFracEffiRatioToHM[iSp][1][1],"same ][");DrawMirrored(hFracExpBiasRatioToHM[iSp][1][1],"same ][");DrawMirrored(hFracCutsRatioToHM[iSp][1][1],"same ][");DrawMirrored(hFracSigExRatioToHM[iSp][1][1],"same ][");if(iSp>1)DrawMirrored(hFracFDRatioToHM[iSp][1][1],"same ][");

		cTotRatioToHM->cd(3)->SetLogx(kTRUE); DrawMirrored(hV0PtSysSumUnc[iSp][1][2]);
		DrawMirrored(hFracEffiRatioToHM[iSp][1][2],"same ][");DrawMirrored(hFracExpBiasRatioToHM[iSp][1][2],"same ][");DrawMirrored(hFracCutsRatioToHM[iSp][1][2],"same ][");DrawMirrored(hFracSigExRatioToHM[iSp][1][2],"same ][");if(iSp>1)DrawMirrored(hFracFDRatioToHM[iSp][1][2],"same ][");
		
		cTotRatioToHM->Write();
		cTotRatioToHM->SaveAs(Form("plotsys/cTotRatioToHM_%s.png",SPECIES[iSp]));
	}

	// IF NOT DOING FULL 4(HM)x8(S0) SYSTEMATIC STUDY, USE HM_01 S0_10 AS REFERENCES
	/*for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
		Int_t iMu = 0; Int_t iSph = 0;
		iMu=1;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][0]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=3;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=4;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iMu=3;	iSph=0; hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][0]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][3][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iMu=2;	iSph=0; hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][0]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=3;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=4;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iMu=4;	iSph=0; hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][0]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][3]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSum[iSp][iMu][iSph] = (TH1F*)hV0PtSysSum[iSp][4][4]->Clone(Form("hV0PtSysSum_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		iSph = 0;
		iMu=1;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][0]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=3;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=4;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iMu=3;	iSph=0; hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][0]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][3][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iMu=2;	iSph=0; hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][0]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=3;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=4;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iMu=4;	iSph=0; hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][0]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=1;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=2;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=5;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=6;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=7;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][3]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		iSph=8;	hV0PtSysSumUnc[iSp][iMu][iSph] = (TH1F*)hV0PtSysSumUnc[iSp][4][4]->Clone(Form("hV0PtSysSumUnc_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
	}*/
}


void MyAnalysisV0syst::AddDeviationsNt() {

	for (int iSp = 1; iSp < NSPECIES; ++iSp)	{
	for (int iReg = 0; iReg < 3; ++iReg)	{
	for (int iRtBin = 0; iRtBin < NRTBINS; ++iRtBin)	{

		//cout << "Calculating relative systematic uncertainties from cut variations for " << hV0PtSysSum[iSp][iMu][iSph] << " " << hV0PtSysSum[iSp][iMu][iSph]->GetName() << endl;
		for (int iBin = 1; iBin < hV0PtRtSysSum[iSp][iReg][iRtBin]->GetNbinsX()+1; iBin++) {
			
			Double_t sum = 0;
			Double_t sumCuts = 0;
			Double_t sumUnc = 0;
			Double_t sumCutsUnc = 0;

			// sources from cut variations
			for (int iSo = 0; iSo < sysSizeof; ++iSo)	{
			if (iSp==1 && iSo==sysLifetime) continue;
					
				sum 	+= hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin)*hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin);
				sumCuts += hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin)*hV0PtRtSysMaxD[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin);

				sumUnc 	+= hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin)*hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin);
				sumCutsUnc += hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin)*hV0PtRtSysMaxDRatioToHM[iSp][iReg][iSo][iRtBin]->GetBinContent(iBin);
				// sum unc denote the uncorrelated uncertainty on ratios to HM, e.g unc of f = RT-biased/RT-int.

			}

			// from material budget, fully correlated between s0
			sum += hV0PtSysBudget->GetBinContent(iBin)*hV0PtSysBudget->GetBinContent(iBin);
			
			// from using MB efficiency
			if (sum>0) {
				sum += hV0PtSysEffi->GetBinContent(iBin)*hV0PtSysEffi->GetBinContent(iBin);
				sumUnc += 1.*hV0PtSysEffi->GetBinContent(iBin)*hV0PtSysEffi->GetBinContent(iBin);
				// factor 2. comes from propagation on ratios
				// for f = A/B, (s_f/f)^2 = (s_A/A)^2 + (s_B/B)^2
				// but (s_A/A) == (s_B/B)
				//// decided to use more liberal 1.0 similarly to other ALICE analyses
			}

			// from experimental bias or unfolding
			
			
			// from signal extraction
			if (sum>0) {
				sum += hV0PtRtSysMaxDSigEx[iSp][iReg][iRtBin]->GetBinContent(iBin)*hV0PtRtSysMaxDSigEx[iSp][iReg][iRtBin]->GetBinContent(iBin);
				sumUnc += hV0PtRtSysMaxDSigExRatioToHM[iSp][iReg][iRtBin]->GetBinContent(iBin)*hV0PtRtSysMaxDSigExRatioToHM[iSp][iReg][iRtBin]->GetBinContent(iBin);
			}

			// from feeddown correction
			if (iSp>1 && sum>0) {
				//sum += hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->GetBinContent(iBin);
				//sumUnc += hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin);
			}

			


			hV0PtRtSysSum[iSp][iReg][iRtBin]->SetBinContent(iBin,TMath::Sqrt(sum));
			hV0PtRtSysSumUnc[iSp][iReg][iRtBin]->SetBinContent(iBin,TMath::Sqrt(sumUnc));
		}
		
		hV0PtRtSysSum[iSp][iReg][iRtBin]->SetFillStyle(3005);
		hV0PtRtSysSum[iSp][iReg][iRtBin]->SetFillColor(kBlack);
		hV0PtRtSysSum[iSp][iReg][iRtBin]->SetTitle(Form("%s %s R_{T} bin:%i; V0 p_{T} (GeV/#it{c}); Relative syst. uncertainty",SPECNAMES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSum[iSp][iReg][iRtBin]->GetYaxis()->SetTitleOffset(1.52);
		hV0PtRtSysSum[iSp][iReg][iRtBin]->SetStats(0);
		hV0PtRtSysSumUnc[iSp][iReg][iRtBin]->SetFillStyle(3005);
		hV0PtRtSysSumUnc[iSp][iReg][iRtBin]->SetFillColor(kBlack);
		hV0PtRtSysSumUnc[iSp][iReg][iRtBin]->SetTitle(Form("%s %s R_{T} bin:%i; V0 p_{T} (GeV/#it{c}); Relative syst. uncertainty",SPECNAMES[iSp],REGIONS[iReg],iRtBin));
		hV0PtRtSysSumUnc[iSp][iReg][iRtBin]->GetYaxis()->SetTitleOffset(1.52);
		hV0PtRtSysSumUnc[iSp][iReg][iRtBin]->SetStats(0);


	}	}	}

	return;

}

void MyAnalysisV0syst::DrawMirrored(TH1F* hist, const char* opt)	{

	TString hname(hist->GetName());
	double range = (hname.Contains("K0s") ? 0.15 : 0.25);
	if (hname.Contains("Unc") || hname.Contains("Ratio") ) range = 0.1; 
	//hist->GetYaxis()->SetRangeUser(-1.5*hist->GetMaximum(),1.5*hist->GetMaximum());
	hist->GetYaxis()->SetRangeUser(-1.*range,range);
	hist->Draw(opt);
	TH1F* htmpMir = (TH1F*)hist->Clone("dtmpMir");
	htmpMir->Scale(-1.);
	htmpMir->Draw(Form("same %s",opt));
}



void MyAnalysisV0syst::CalculateFeeddownSys() {

	//mHandler->analysis(2)->dirFile()->ls();
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)		{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		hV0PtFeeddown[iSp][iMu][iSph]			
			= (TH1F*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtFeeddown_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtFeeddownXi0[iSp][iMu][iSph]			
			= (TH1F*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtFeeddownXi0_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));
		hV0PtFeeddownXiErr[iSp][iMu][iSph]			
			= (TH1F*)mHandler->analysis(2)->dirFile()->Get(Form("hV0PtFeeddownXiErr_%s_%s_%s",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));

		for( Int_t i=1; i<hV0PtFeeddown[iSp][iMu][iSph]->GetNbinsX()+1; i++){	
  			hV0PtFeeddown[iSp][iMu][iSph]->SetBinContent(i,1.-hV0PtFeeddown[iSp][iMu][iSph]->GetBinContent(i));
  			hV0PtFeeddownXi0[iSp][iMu][iSph]->SetBinContent(i,1.-hV0PtFeeddownXi0[iSp][iMu][iSph]->GetBinContent(i));
  			hV0PtFeeddownXiErr[iSp][iMu][iSph]->SetBinContent(i,1.-hV0PtFeeddownXiErr[iSp][iMu][iSph]->GetBinContent(i));
  		}
	}	} }

	cout << "gothere" << endl;

	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

	cout << "Doing FD systematics for " << iMu << " " << iSph << endl;

		// FROM METHOD
		Double_t lSigmaDelta[100]; 
  		for( Int_t i=1; i<hV0PtFeeddown[iSp][iMu][iSph]->GetNbinsX()+1; i++){
    		//Computation of roger barlow sigma_{delta} 
    		lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(hV0PtFeeddown[iSp][iMu][iSph]->GetBinError(i),2) - TMath::Power(hV0PtFeeddownXi0[iSp][iMu][iSph]->GetBinError(i),2) ) );
 		   //Computation of relationship to h2 for plotting in ratio plot 
    		if ( hV0PtFeeddown[iSp][iMu][iSph]->GetBinContent(i) > 1e-12 )	{ 
    			lSigmaDelta[i] /= hV0PtFeeddown[iSp][iMu][iSph]->GetBinContent(i);		}
    		else	{ 
    			lSigmaDelta[i] = 0; 	}
  		} 
		// dividing corrected varied yields by default
		hV0PtSysFeeddown[iSp][iMu][iSph]->Divide(hV0PtFeeddown[iSp][iMu][iSph],hV0PtFeeddownXi0[iSp][iMu][iSph],1.,1.,"");
		// apply RB sigmas as errors
  		for( Int_t i=1; i<hV0PtFeeddown[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
    		hV0PtSysFeeddown[iSp][iMu][iSph]->SetBinError(i, lSigmaDelta[i]);
  		}

  		cout << "gothereee" << endl;

  		// FROM XI FIT VARIATION
  		for( Int_t i=1; i<hV0PtFeeddown[iSp][iMu][iSph]->GetNbinsX()+1; i++){
    		//Computation of roger barlow sigma_{delta} 
    		lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(hV0PtFeeddownXiErr[iSp][iMu][iSph]->GetBinError(i),2) - TMath::Power(hV0PtFeeddownXi0[iSp][iMu][iSph]->GetBinError(i),2) ) );
 		   //Computation of relationship to h2 for plotting in ratio plot 
    		if ( hV0PtFeeddownXi0[iSp][iMu][iSph]->GetBinContent(i) > 1e-12 )	{ 
    			lSigmaDelta[i] /= hV0PtFeeddownXi0[iSp][iMu][iSph]->GetBinContent(i);		}
    		else	{ 
    			lSigmaDelta[i] = 0; 	}
  		} 
		// dividing corrected varied yields by default
		hV0PtSysFeeddownXiErr[iSp][iMu][iSph]->Divide(hV0PtFeeddownXiErr[iSp][iMu][iSph],hV0PtFeeddownXi0[iSp][iMu][iSph],1.,1.,"");
		// apply RB sigmas as errors
  		for( Int_t i=1; i<hV0PtFeeddown[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
    		hV0PtSysFeeddownXiErr[iSp][iMu][iSph]->SetBinError(i, lSigmaDelta[i]);
  		}

  	
  		cout << "gothereeeeeeee" << endl;

  		// CALCULATE MAX DEVIATIONS
  		for (int iBin = 1; iBin<hV0PtSysFeeddown[iSp][iMu][iSph]->GetNbinsX()+1;iBin++) {
			Double_t maxD = 0;
			Double_t varD = TMath::Abs(hV0PtSysFeeddown[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
			if (varD>maxD && varD>nRBsigmas*hV0PtSysFeeddown[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
					
			if (hV0PtSysFeeddown[iSp][iMu][iSph]->GetBinContent(iBin)>0) hV0PtSysMaxDFeeddown[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(maxD*maxD));
			
			maxD = 0;
			varD = TMath::Abs(hV0PtSysFeeddownXiErr[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
			if (varD>maxD && varD>nRBsigmas*hV0PtSysFeeddownXiErr[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
					
			if (hV0PtSysFeeddownXiErr[iSp][iMu][iSph]->GetBinContent(iBin)>0) hV0PtSysMaxDFeeddownXiErr[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(maxD*maxD));

			hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(0.006*0.006 +
				hV0PtSysMaxDFeeddown[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddown[iSp][iMu][iSph]->GetBinContent(iBin) +
				hV0PtSysMaxDFeeddownXiErr[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownXiErr[iSp][iMu][iSph]->GetBinContent(iBin) ));			
			
		}

		cout << "gothereeeEEEEE" << endl;

		mHandler->MakeNiceHistogram(hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph],kBlack);
		hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->SetFillStyle(3002);
		hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->SetFillColor(kBlack);
		hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->GetYaxis()->SetRangeUser(-0.0005,0.2);
		hV0PtSysMaxDFeeddownTotal[iSp][iMu][iSph]->SetTitle(Form("%s %s %s; V0 p_{T} (GeV/#it{c}); Relative syst. uncertainty",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));	

	}	}	}


	// and now study the effect on s0/hm ratio
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	for (int iMu = 0; iMu < NMULTI; ++iMu)	{
	for (int iSph = 0; iSph < NSPHERO; ++iSph)	{
	if (iMu == 0 && iSph > 0) continue;

		// FROM METHOD
		hV0PtSysFeeddownRatioToHM[iSp][iMu][iSph]->Divide(hV0PtSysFeeddown[iSp][iMu][iSph],hV0PtSysFeeddown[iSp][iMu][0],1,1,"");
		
		// and calculate uncertainties
		for( Int_t i=1; i<hV0PtSysFeeddownRatioToHM[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
			Double_t errS0 = hV0PtSysFeeddown[iSp][iMu][iSph]->GetBinError(i);
			Double_t errHM = hV0PtSysFeeddown[iSp][iMu][0]->GetBinError(i);
			Double_t errR = TMath::Sqrt(TMath::Abs(errS0*errS0 - errHM*errHM));
			errR = (hV0PtSysFeeddown[iSp][iMu][0]->GetBinContent(i) > 0) ? errR / hV0PtSysFeeddown[iSp][iMu][0]->GetBinContent(i) : 0;
			hV0PtSysFeeddownRatioToHM[iSp][iMu][iSph]->SetBinError(i,errR);
		}

  		// FROM XI FIT VARIATION
  		hV0PtSysFeeddownXiErrRatioToHM[iSp][iMu][iSph]->Divide(hV0PtSysFeeddownXiErr[iSp][iMu][iSph],hV0PtSysFeeddownXiErr[iSp][iMu][0],1,1,"");
		for( Int_t i=1; i<hV0PtSysFeeddownXiErrRatioToHM[iSp][iMu][iSph]->GetNbinsX()+1; i++){ 
			Double_t errS0 = hV0PtSysFeeddownXiErr[iSp][iMu][iSph]->GetBinError(i);
			Double_t errHM = hV0PtSysFeeddownXiErr[iSp][iMu][0]->GetBinError(i);
			Double_t errR = TMath::Sqrt(TMath::Abs(errS0*errS0 - errHM*errHM));
			errR = (hV0PtSysFeeddownXiErr[iSp][iMu][0]->GetBinContent(i) > 0) ? errR / hV0PtSysFeeddownXiErr[iSp][iMu][0]->GetBinContent(i) : 0;
			hV0PtSysFeeddownXiErrRatioToHM[iSp][iMu][iSph]->SetBinError(i,errR);
		}  		
  	
  		// CALCULATE MAX DEVIATIONS
  		for (int iBin = 1; iBin<hV0PtSysFeeddownRatioToHM[iSp][iMu][iSph]->GetNbinsX()+1;iBin++) {
			Double_t maxD = 0;
			Double_t varD = TMath::Abs(hV0PtSysFeeddownRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
			if (varD>maxD && varD>nRBsigmas*hV0PtSysFeeddownRatioToHM[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
					
			if (hV0PtSysFeeddown[iSp][iMu][iSph]->GetBinContent(iBin)>0) hV0PtSysMaxDFeeddownRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(maxD*maxD));
			
			maxD = 0;
			varD = TMath::Abs(hV0PtSysFeeddownXiErrRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)-1.);
			if (varD>maxD && varD>nRBsigmas*hV0PtSysFeeddownXiErrRatioToHM[iSp][iMu][iSph]->GetBinError(iBin)) maxD=varD;
					
			if (hV0PtSysFeeddownXiErr[iSp][iMu][iSph]->GetBinContent(iBin)>0) hV0PtSysMaxDFeeddownXiErrRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(maxD*maxD));

			hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->SetBinContent(iBin,TMath::Sqrt(1.*0.006*0.006 +
				hV0PtSysMaxDFeeddownRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin) +
				hV0PtSysMaxDFeeddownXiErrRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin)*hV0PtSysMaxDFeeddownXiErrRatioToHM[iSp][iMu][iSph]->GetBinContent(iBin) ));			
			// factor 2.0 comes from error propagation on ratios
			// decided to use more liberal 1.0 similarly to other ALICE analyses
		}

		mHandler->MakeNiceHistogram(hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph],kBlack);
		hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->SetFillStyle(3002);
		hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->SetFillColor(kBlack);
		hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->GetYaxis()->SetRangeUser(-0.0005,0.2);
		hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][iMu][iSph]->SetTitle(Form("%s %s %s; V0 p_{T} (GeV/#it{c}); Relative syst. uncertainty",SPECIES[iSp],MULTI[iMu],SPHERO[iSph]));	

	}	}	}

	mHandler->root()->SetBatch(kTRUE);
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	if (NMULTI<5) continue;
	if (NSPHERO<5) continue;
		TCanvas* cFD = new TCanvas(Form("cFD_%s",SPECIES[iSp]),"",1600,700);
		cFD->Divide(2,1,5e-5,5e-5); 
		Int_t padC = 1;

		cout << "ggggGGGgothere" << endl;

		cFD->cd(1)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDFeeddownTotal[iSp][3][0]);
		cFD->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDFeeddownTotal[iSp][4][0]);
		//cSig->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][3][3]);
		cFD->Write();
		cFD->SaveAs(Form("plotsys/cFD_%s.png",SPECIES[iSp]));
	}
	mHandler->root()->SetBatch(kFALSE);


	mHandler->root()->SetBatch(kTRUE);
	for (int iSp = 2; iSp < NSPECIES; ++iSp)	{
	if (NMULTI<2) continue;
	if (NSPHERO<3) continue;
		TCanvas* cFDRatioToHM = new TCanvas(Form("cFDRatioToHM_%s",SPECIES[iSp]),"",1500,450);
		cFDRatioToHM->Divide(3,1,5e-5,5e-5); 
		Int_t padC = 1;

		cFDRatioToHM->cd(1)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][1][0]);
		cFDRatioToHM->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][1][1]);
		cFDRatioToHM->cd(3)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDFeeddownTotalRatioToHM[iSp][1][2]);
		//cSig->cd(2)->SetLogx(kTRUE); DrawMirrored(hV0PtSysMaxDSigEx[iSp][3][3]);
		cFDRatioToHM->Write();
		cFDRatioToHM->SaveAs(Form("plotsys/cFDRatioToHM_%s.png",SPECIES[iSp]));
	}
	mHandler->root()->SetBatch(kFALSE);

}