// AliAnalysisPIDCascadeV0 wrapper class
// OliverM 2019 Lund

#include "compInstructions.h"

#ifndef MYCASCADE_H
#define MYCASCADE_H
#include <AliAnalysisPIDCascade.h>
#include "TObject.h"
#include "TMath.h"

#include "MyHandler.h"



#if INPUTFORMAT == 1
	typedef AliAnalysisPIDCascade AnyCascade;
	typedef AliAnalysisPIDCascadeV0 AnyV0;
	typedef AliAnalysisPIDCascadeTrack AnyTrack;
#endif

class MyCascade: public TObject {

	public:
		MyCascade() { }
		MyCascade(AnyCascade* cascade);// : mAliCascade(v0) { }	
		~MyCascade() { }
		void SetHandler(MyHandler* handler)		{ mHandler = handler;};

#if INPUTFORMAT == 1
		Double_t GetPt() 					const { return mAliCascade->GetPtCasc();};
		Double_t GetEta() 					const { return mAliCascade->GetEtaCasc();};
		AnyTrack* GetBachTrack() 	const { return mAliCascade->GetBachAnalysisTrack();};
		AnyV0* GetV0()	 			const { return mAliCascade->GetV0();};
		Double_t GetPhi();
		
		Double_t GetBachDCApvXY()	const { return mAliCascade->GetBachAnalysisTrack()->GetImpactParameter(0);};
		Double_t GetV0posDCApvXY()	const { return mAliCascade->GetV0()->GetPosAnalysisTrack()->GetImpactParameter(0);};
		Double_t GetV0negDCApvXY()	const { return mAliCascade->GetV0()->GetNegAnalysisTrack()->GetImpactParameter(0);};

		Int_t GetMCLabel() 				const { return mAliCascade->GetBachAnalysisTrack()->GetMCMotherLabel();};
		Int_t IsMCPrimary() 				const { return mAliCascade->GetBachAnalysisTrack()->GetMCMotherPrimary();};
		Int_t GetMCPdgCode()				const { return mAliCascade->GetBachAnalysisTrack()->GetMCMotherPdgCode();};
		
#endif		

		ClassDef(MyCascade,1);

	private:
		AnyCascade* mAliCascade;
		MyHandler* mHandler = 0;
		Bool_t mPhicalculated;
		Double_t mPhi;

	
};
#endif