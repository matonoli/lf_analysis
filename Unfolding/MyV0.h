// AliAnalysisPIDCascadeV0 wrapper class
// OliverM 2019 Lund

#include "compInstructions.h"

#ifndef MYV0_H
#define MYV0_H
#include <AliAnalysisPIDCascadeV0.h>
#include "TObject.h"
#include "TMath.h"

#include "MyHandler.h"

#include <AliESDtrack.h>
#include <AliESDv0.h>


#if INPUTFORMAT == 1
	typedef AliAnalysisPIDCascadeV0 AnyV0;
	typedef AliAnalysisPIDCascadeTrack AnyTrack;
#elif INPUTFORMAT == 2
	typedef AliESDv0 AnyV0;
	typedef AliESDtrack AnyTrack;
#endif

class MyV0: public TObject {

	public:
		MyV0() { }
		MyV0(AnyV0* v0);// : mAliV0(v0) { }	
		~MyV0() { }
		void SetHandler(MyHandler* handler)		{ mHandler = handler;};
		Double_t GetCut(Int_t Cut);

#if INPUTFORMAT == 1
		Double_t GetPt() 					const { return mAliV0->GetPt();};
		Double_t GetEta() 					const { return mAliV0->GetEta();};
		Double_t GetDCAdd() 				const { return mAliV0->GetDCAV0Daughters();};
		Double_t GetCPA() 					const { return mAliV0->GetV0CosinePA();};
		Double_t GetRadius() 				const { return mAliV0->GetRadius();};
		Double_t GetIMK0s() 				const { return mAliV0->GetIMK0s();};
		Double_t GetIML()					const { return mAliV0->GetIML();};
		Double_t GetIMLbar() 				const { return mAliV0->GetIMAL();};
		Double_t GetKFIMK0s() 				const { return -1;};	//not implemented yet
		Double_t GetKFIML()					const { return -1;};	//not implemented yet
		Double_t GetKFIMLbar() 				const { return -1;};	//not implemented yet
		AnyTrack* GetPosTrack() 	const { return mAliV0->GetPosAnalysisTrack();};
		AnyTrack* GetNegTrack() 	const { return mAliV0->GetNegAnalysisTrack();};
		Double_t* CalculateAP();
		Double_t GetPhi();
		Double_t CalculateY(Int_t Sp);
		Int_t HasFastSignal();
		Bool_t IsOffline()				const { return 1;};

		Int_t GetMCLabel() 				const { return mAliV0->GetPosAnalysisTrack()->GetMCMotherLabel();};
		Int_t IsMCPrimary() 				const { return mAliV0->GetPosAnalysisTrack()->GetMCMotherPrimary();};
		Int_t GetMCPrimaryPdgCode()			const { return mAliV0->GetPosAnalysisTrack()->GetMCPrimaryPdgCode();};
		Int_t GetMCPrimaryLabel()			const { return mAliV0->GetPosAnalysisTrack()->GetMCPrimaryLabel();};
		Int_t GetMCPrimaryLabelNeg()		const { return mAliV0->GetNegAnalysisTrack()->GetMCPrimaryLabel();};
		Int_t GetMCPdgCode()				const { return mAliV0->GetMCPdgCode();};
		Int_t GetPosTrackPdg()				{ return mAliV0->GetPosAnalysisTrack()->GetMCPdgCode();};
		Int_t GetNegTrackPdg()				{ return mAliV0->GetNegAnalysisTrack()->GetMCPdgCode();}; // not implemented yet

#elif INPUTFORMAT == 2
		Double_t GetPt() 					const { return mAliV0->Pt();};
		Double_t GetEta() 					const { return mAliV0->Eta();};
		Double_t GetDCAdd() 				const { return mAliV0->GetDcaV0Daughters();};
		Double_t GetCPA() 					const { return mAliV0->GetV0CosineOfPointingAngle();};
		Double_t GetRadius() 				const { return TMath::Sqrt(mAliV0->Xv()*mAliV0->Xv()+mAliV0->Yv()*mAliV0->Yv());};
		Double_t GetIMK0s();
		Double_t GetIML();
		Double_t GetIMLbar();
		Double_t GetKFIMK0s();
		Double_t GetKFIML();
		Double_t GetKFIMLbar();
		AnyTrack* GetPosTrack() 		const { return mHandler->track(TMath::Abs(mAliV0->GetPindex()));};
		AnyTrack* GetNegTrack() 		const { return mHandler->track(TMath::Abs(mAliV0->GetNindex()));};
		Double_t* CalculateAP();
		Double_t GetPhi()				const { return mAliV0->Phi();};
		Double_t CalculateY(Int_t Sp)	const { return ((Sp==1) ? mAliV0->RapK0Short() : ((Sp==2 || Sp==3) ? mAliV0->RapLambda() : -999));};
		Int_t HasFastSignal();
		Bool_t IsOffline()				const { return !mAliV0->GetOnFlyStatus();};

		Double_t CalculateY(Int_t Sp);

		Int_t GetMCLabel() const;
		Int_t GetMCPdgCode(); //				const { return mAliV0->GetPdgCode();}; // AliESDv0::GetPdgCode() always returns 310
		Int_t GetPosTrackPdg();
		Int_t GetNegTrackPdg();
		Int_t IsMCPrimary() 				const { return mHandler->mcstack()->IsPhysicalPrimary(this->GetMCLabel());};

		//Int_t GetMCPrimaryPdgCode()			const { return mAliV0->GetPosAnalysisTrack()->GetMCPrimaryPdgCode();};
		//Int_t GetMCPrimaryLabel()			const { return mAliV0->GetPosAnalysisTrack()->GetMCPrimaryLabel();};

#endif


		
		ClassDef(MyV0,1);

	private:
		AnyV0* mAliV0;
		MyHandler* mHandler = 0;
		Bool_t mAPcalculated;
		Double_t mAP[2];
		Bool_t mPhicalculated;
		Double_t mPhi;
		Bool_t mYcalculated[4];
		Double_t mY[4];

		enum sysCuts { NoCut,														//0
		DCAdd, CPA, RadiusL, RadiusH,											//4
		FastSignal, CompMassK0s, CompMassL,						//8
		CompMassLbar, LifetimeK0s,					//12
		LifetimeL, LifetimeLbar, 												//14
		cutsV0cuts,																//15
		NSigmaTPCposPiL, NSigmaTPCposPiH, NSigmaTPCnegPiL, NSigmaTPCnegPiH,		//19
		NSigmaTPCposPrL, NSigmaTPCposPrH, NSigmaTPCnegPrL, NSigmaTPCnegPrH,		//23
		DCAPVpos, NClusterpos, NClusterFpos, 									//26
		DCAPVneg, NClusterneg, NClusterFneg,  									//29
		cutsSizeof 																//30
	};
};
#endif