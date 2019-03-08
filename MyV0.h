// AliAnalysisPIDV0 wrapper class
// OliverM 2019 Lund

#ifndef __MyV0__
#define __MyV0__
#include <AliAnalysisPIDV0.h>
#include "TObject.h"

class AliAnalysisPIDTrack; // forward declaration

class MyV0: public TObject {

	public:
		MyV0() { }
		MyV0(AliAnalysisPIDV0* v0);// : mAliV0(v0) { }	
		~MyV0() { }
		Double_t GetPt() 					const { return mAliV0->GetPt();};
		Double_t GetEta() 					const { return mAliV0->GetEta();};
		Double_t GetDCAdd() 				const { return mAliV0->GetDCAV0Daughters();};
		Double_t GetCPA() 					const { return mAliV0->GetV0CosinePA();};
		Double_t GetRadius() 				const { return mAliV0->GetRadius();};
		Double_t GetIMK0s() 				const { return mAliV0->GetIMK0s();};
		Double_t GetIML()					const { return mAliV0->GetIML();};
		Double_t GetIMLbar() 				const { return mAliV0->GetIMAL();};
		AliAnalysisPIDTrack* GetPosTrack() 	const { return mAliV0->GetPosAnalysisTrack();};
		AliAnalysisPIDTrack* GetNegTrack() 	const { return mAliV0->GetNegAnalysisTrack();};

		Double_t GetMCLabel() 				const { return mAliV0->GetPosAnalysisTrack()->GetMCMotherLabel();};
		Double_t GetMCPdgCode()				const { return mAliV0->GetMCPdgCode();};

		Double_t* CalculateAP();
		
		ClassDef(MyV0,1);

	private:
		AliAnalysisPIDV0* mAliV0;
		Bool_t mAPcalculated;
		Double_t mAP[2];
};
#endif