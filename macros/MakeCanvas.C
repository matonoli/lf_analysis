#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"

/// \file
/// \ingroup tutorial_graphics
/// \notebook
/// Example of canvas partitioning.
/// Sometimes the Divide() method is not appropriate to divide a Canvas.
/// Because of the left and right margins, all the pads do not have the
/// same width and height. CanvasPartition does that properly. This
/// example also ensure that the axis labels and titles have the same
/// sizes and that the tick marks length is uniform.
///
/// \macro_image
/// \macro_code
///
/// \author Olivier Couet

void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
		Float_t lMargin = 0., Float_t rMargin = 0.,
		Float_t bMargin = 0.15, Float_t tMargin = 0.05);

TPad* GetTwoPadsTitleY(const Char_t* name = "pad_TitleX",
		Float_t hposlc = 0.1 , Float_t vposdc = 0.1);

void canvas2()
{

	gStyle->SetOptStat(0);

	TCanvas *C = (TCanvas*) gROOT->FindObject("C");
	if (C) delete C;
	C = new TCanvas("C","canvas",1024,640);
	C->SetFillStyle(4000);

	// Number of PADS
	const Int_t Nx = 3;
	const Int_t Ny = 2;

	// Margins
	Float_t lMargin = 0.08;
	Float_t rMargin = 0.05;
	Float_t bMargin = 0.08;
	Float_t tMargin = 0.05;

		// Canvas setup
	CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

	// Dummy histogram.
	TH1F *h = (TH1F*) gROOT->FindObject("histo");
	if (h) delete h;
	h = new TH1F("histo","",100,-5.0,5.0);
	h->FillRandom("gaus",10000);
	h->GetXaxis()->SetTitle("x axis");
	h->GetYaxis()->SetTitle("y axis");

	C->cd(0);

	TPad* pTitleX = (TPad*) gROOT->FindObject("pad_TitleX");
	if(!pTitleX) printf("pTitleX missing\n");
	pTitleX->SetFillStyle(1001);
	pTitleX->SetFillColor(8);
	pTitleX->Draw();
	pTitleX->cd();

	TPad* pTitleY1 = (TPad*) gROOT->FindObject("pad_TitleY1");
	if(!pTitleY1) printf("pTitleY1 missing\n");
	pTitleY1->SetFillStyle(1001);
	pTitleY1->SetFillColor(9);
	C->cd(0);
	pTitleY1->Draw();
	pTitleY1->cd();

	TPad* pTitleY2 = (TPad*) gROOT->FindObject("pad_TitleY2");
	if(!pTitleY2) printf("pTitleY2 missing\n");
	pTitleY2->SetFillStyle(1001);
	pTitleY2->SetFillColor(12);
	C->cd(0);
	pTitleY2->Draw();
	pTitleY2->cd();

	TPad *pad[Nx][Ny];

	for (Int_t i=0;i<Nx;i++) {
		for (Int_t j=0;j<Ny;j++) {
			C->cd(0);

			// Get the pads previously created.
			char pname[16];
			sprintf(pname,"pad_%i_%i",i,j);
			pad[i][j] = (TPad*) gROOT->FindObject(pname);
			pad[i][j]->Draw();
			pad[i][j]->SetFillStyle(1001);
			//		pad[i][j]->SetFrameFillStyle(0);
			pad[i][j]->SetFillColor(i+j+2);
			pad[i][j]->cd();

			// Size factors
			Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
			Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();

			char hname[16];
			sprintf(hname,"h_%i_%i",i,j);
			TH1F *hFrame = (TH1F*) h->Clone(hname);
			hFrame->Reset();
			hFrame->Draw();

			// y axis range
			hFrame->GetYaxis()->SetRangeUser(0.0001,1.2*h->GetMaximum());

			// Format for y axis
			hFrame->GetYaxis()->SetLabelFont(43);
			hFrame->GetYaxis()->SetLabelSize(16);
			hFrame->GetYaxis()->SetLabelOffset(0.02);
			hFrame->GetYaxis()->SetTitleFont(43);
			hFrame->GetYaxis()->SetTitleSize(16);
			hFrame->GetYaxis()->SetTitleOffset(5);

			hFrame->GetYaxis()->CenterTitle();
			hFrame->GetYaxis()->SetNdivisions(505);

			// TICKS Y Axis
			hFrame->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);

			// Format for x axis
			hFrame->GetXaxis()->SetLabelFont(43);
			hFrame->GetXaxis()->SetLabelSize(16);
			hFrame->GetXaxis()->SetLabelOffset(0.02);
			hFrame->GetXaxis()->SetTitleFont(43);
			hFrame->GetXaxis()->SetTitleSize(16);
			hFrame->GetXaxis()->SetTitleOffset(5);
			hFrame->GetXaxis()->CenterTitle();
			hFrame->GetXaxis()->SetNdivisions(505);

			// TICKS X Axis
			hFrame->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);

			h->Draw("same");
		}
	}
	C->cd();
}

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
		Float_t lMargin, Float_t rMargin,
		Float_t bMargin, Float_t tMargin)
{
	if (!C) return;

	printf("\t lMargin = %f \t rMargin = %f\n",lMargin,rMargin);

	const Float_t hposlc = 0.08;
	const Float_t vposdc = 0.08;

	// Setup Pad layout:
	Float_t vSpacing = 0.0;
	Float_t vStep  = (1.- vposdc - bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

	Float_t hSpacing = 0.0;
	Float_t hStep  = (1.- hposlc - lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
	printf("hStep = %f\n",hStep);

	Float_t vposd,vposu,vmard,vmaru,vfactor;
	Float_t hposl,hposr,hmarl,hmarr,hfactor;

	for (Int_t i=0;i<Nx;i++) {

		if (i==0) {
			hposl = hposlc;
			hposr = lMargin + hStep + hposlc;
			hfactor = hposr-hposl;
			hmarl = lMargin / hfactor;
			hmarr = 0.0;
		} else if (i == Nx-1) {
			hposl = hposr + hSpacing;
			hposr = hposl + hStep + rMargin;
			hfactor = hposr-hposl;
			hmarl = 0.0;
			hmarr = rMargin / (hposr-hposl);
		} else {
			hposl = hposr + hSpacing;
			hposr = hposl + hStep;
			hfactor = hposr-hposl;
			hmarl = 0.0;
			hmarr = 0.0;
		}

		printf("i = %d \t hposl = %f \t hposr = %f\n",i,hposl,hposr);
		printf("hmarl = %f \t hmarr = %f\n",hmarl,hmarr);

		for (Int_t j=0;j<Ny;j++) {

			if (j==0) {
				vposd = vposdc;
				vposu = bMargin + vStep + vposdc;
				vfactor = vposu-vposd;
				vmard = bMargin / vfactor;
				vmaru = 0.0;
			} else if (j == Ny-1) {
				vposd = vposu + vSpacing;
				vposu = vposd + vStep + tMargin;
				vfactor = vposu-vposd;
				vmard = 0.0;
				vmaru = tMargin / (vposu-vposd);
			} else {
				vposd = vposu + vSpacing;
				vposu = vposd + vStep;
				vfactor = vposu-vposd;
				vmard = 0.0;
				vmaru = 0.0;
			}

			C->cd(0);

			char name[16];
			sprintf(name,"pad_%i_%i",i,j);
			TPad *pad = (TPad*) gROOT->FindObject(name);
			if (pad) delete pad;
			pad = new TPad(name,"",hposl,vposd,hposr,vposu);
			pad->SetLeftMargin(hmarl);
			pad->SetRightMargin(hmarr);
			pad->SetBottomMargin(vmard);
			pad->SetTopMargin(vmaru);

			pad->SetFrameBorderMode(0);
			pad->SetBorderMode(0);
			pad->SetBorderSize(0);

			pad->Draw();
		}
	}

	TPad *padTitleX = GetTwoPadsTitleY("pad_TitleX",hposlc,vposdc);
	TPad *padTitleY1 = GetTwoPadsTitleY("pad_TitleY1",hposlc,vposdc);
	TPad *padTitleY2 = GetTwoPadsTitleY("pad_TitleY2",hposlc,vposdc);

	C->cd();
	padTitleX->Draw();
	padTitleY1->Draw();
	padTitleY2->Draw();

}

TPad* GetTwoPadsTitleY(const Char_t* name, Float_t hposlc, Float_t vposdc)
{

	TPad *pad = (TPad*) gROOT->FindObject(name);
	if (pad) delete pad;

	if(strcmp(name,"pad_TitleX")==0) pad = new TPad("pad_TitleX","",hposlc,0.,1.,vposdc);
	else if(strcmp(name,"pad_TitleY1")==0)pad = new TPad("pad_TitleY1","", 0., vposdc, hposlc, vposdc + (1.-vposdc)/2. );
	else pad = new TPad("pad_TitleY2","",0. ,vposdc + (1.-vposdc)/2., hposlc, 1.);

	pad->SetLeftMargin(0.);
	pad->SetRightMargin(0.);
	pad->SetBottomMargin(0.0);
	pad->SetTopMargin(0.);

	pad->SetFrameBorderMode(0);
	pad->SetBorderMode(0);
	pad->SetBorderSize(0);

	return pad;
}

