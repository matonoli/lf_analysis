// AliAnalysisPIDCascadeTrack wrapper class
// OliverM 2019 Lund

#include "compInstructions.h"

#ifndef MYTRACK_H
#define MYTRACK_H
#include <AliAnalysisPIDCascadeTrack.h>

#include <AliESDtrack.h>
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"

#include "TObject.h"
#include "TMath.h"

#if INPUTFORMAT == 1
	typedef AliAnalysisPIDCascadeTrack AnyTrack;
#elif INPUTFORMAT == 2
	typedef AliESDtrack AnyTrack;
#endif

class MyHandler;

class MyTrack: public TObject {

	public:
		MyTrack() : mAliTrack(0) { }
		MyTrack(AnyTrack* tr) : mAliTrack(tr) { }
		~MyTrack() { }
		void SetHandler(MyHandler* handler)		{ mHandler = handler;};
#if INPUTFORMAT == 1
		Float_t GetPt() 						const { return mAliTrack->GetPt();};
		Float_t GetPx() 						const { return mAliTrack->GetPt()*TMath::Cos(mAliTrack->GetPhi());};
		Float_t GetPy() 						const { return mAliTrack->GetPt()*TMath::Sin(mAliTrack->GetPhi());};
		Float_t GetEta() 						const { return mAliTrack->GetEta();};
		Float_t GetPhi() 						const { return mAliTrack->GetPhi();};
		Double_t GetSign()						const { return mAliTrack->GetSign();};
		Float_t GetNSigmaPionTPC() 				const { return mAliTrack->GetNSigmaPionTPC();};
		Float_t GetNSigmaProtonTPC() 			const { return mAliTrack->GetNSigmaProtonTPC();};
		Float_t GetDCApvXY() 					const { return mAliTrack->GetImpactParameter(0);};
		Bool_t IskITSrefit()					const { return (mAliTrack->GetStatus()&AliESDtrack::kITSrefit);};
		Bool_t IsITSTPC2011()					const { return (mAliTrack->GetTrackCutFlag()&2);};
		Bool_t IsTPCOnlyRefit()					const { return (mAliTrack->GetTrackCutFlag()&4);};
		Bool_t IsGoodV0daughter()				const { return (mAliTrack->GetTrackCutFlag()&16);};
		Bool_t IsITSTPC2011HybridOff()			const { return (mAliTrack->GetTrackCutFlag()&32);};
		Bool_t IsITSTPC2011HybridNone()			const { return (mAliTrack->GetTrackCutFlag()&64);};
		Bool_t HasTOFPID()						const { return mAliTrack->HasTOFPID();};
		Float_t GetTPCnc()						const { return mAliTrack->GetTPCnc();};
		UShort_t GetTPCNclsF()					const { return mAliTrack->GetTPCNclsF();};
		Bool_t IsMCPrimary()					const { return mAliTrack->IsMCPrimary();};
		Int_t GetMCPdgCode()					const { return mAliTrack->GetMCPdgCode();};
#elif INPUTFORMAT == 2		
		Float_t GetPt() 						const { return mAliTrack->Pt();};
		Float_t GetPx() 						const { return mAliTrack->Px();};
		Float_t GetPy() 						const { return mAliTrack->Py();};
		Float_t GetEta() 						const { return mAliTrack->Eta();};
		Float_t GetPhi() 						const { return mAliTrack->Phi();};
		Float_t GetNSigmaPionTPC() 				const { return ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse()->NumberOfSigmasTPC(mAliTrack, AliPID::kPion);}; // rather slow
		Float_t GetNSigmaProtonTPC() 			const { return ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse()->NumberOfSigmasTPC(mAliTrack, AliPID::kProton);}; // rather slow
		Float_t GetDCApvXY();
		Bool_t IskITSrefit()					const { return (mAliTrack->GetStatus()&AliESDtrack::kITSrefit);};
		Bool_t IsITSTPC2011()					{ return (!mFlagCal) ? CalculateFlag() &2 : mFlag&2;};
		Bool_t IsTPCOnlyRefit()					{ return (!mFlagCal) ? CalculateFlag() &4 : mFlag&4;};
		Bool_t IsGoodV0daughter()				{ return (!mFlagCal) ? CalculateFlag() &16 : mFlag&16;};
		Int_t CalculateFlag();
		Bool_t HasTOFPID();
		Float_t GetTPCnc()						const { return mAliTrack->GetTPCClusterInfo(2,1);};
		UShort_t GetTPCNclsF()					const { return mAliTrack->GetTPCNclsF();};
#endif
#if INPUTFORMAT == 1		



		// add safety measure to getters for if mAliTrack is an invalid pointer
#endif
		ClassDef(MyTrack,1);

	private:
		AnyTrack* mAliTrack;
		Bool_t mFlagCal = 0;
		Int_t mFlag = 0;

		MyHandler* mHandler = 0; 
};
#endif