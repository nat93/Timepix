//ROOT
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TH2Poly.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TPaletteAxis.h"
//C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

int plot_taratin_extraction()
{
   	TString fileName1 = "output_function_1_taratin_extraction.root";
	cout<<"--> Infile: "<<fileName1<<endl;

	TFile *_file1 = TFile::Open(fileName1.Data());

	TH2D*	hh_1 = (TH2D*)_file1->Get("h_5");
	hh_1->SetName("hh_1");

	hh_1->Rebin2D(2);

	Double_t integral, integralErr;

	integral = hh_1->IntegralAndError(1,hh_1->GetNbinsX(),1,hh_1->GetNbinsY(),integralErr);
	scale2Dhisto(hh_1, integral,integralErr);

	TCanvas* c_0 = new TCanvas("c_0","c_0",1000,1000);
        c_0->cd();
        TPad *pad_0 = new TPad("pad_0","pad_0",0,0,1,1);
        pad_0->SetLeftMargin(0.15);
        pad_0->SetRightMargin(0.2);
        pad_0->Draw();
        pad_0->cd();
	gStyle->SetOptStat(0);
        gPad->SetGrid();
	hh_1->GetYaxis()->SetTitle("Vertical Axis [mm]");
	hh_1->GetXaxis()->SetTitle("Horizontal Axis [mm]");
        hh_1->GetZaxis()->SetTitle("Hits Density");
        hh_1->GetXaxis()->CenterTitle();
        hh_1->GetYaxis()->CenterTitle();
        hh_1->GetZaxis()->CenterTitle();
        hh_1->GetYaxis()->SetTitleOffset(1.2);
        hh_1->GetZaxis()->SetTitleOffset(2.1);
        hh_1->Draw("colz");

	TH1D* projhh_1 = hh_1->ProjectionX("projhh_1");
	
	projhh_1->SetTitle("");

	projhh_1->GetXaxis()->SetTitle("Horizontal Axis Projection [mm]");
        projhh_1->GetYaxis()->SetTitle("Hits Density");
        projhh_1->GetXaxis()->CenterTitle();
        projhh_1->GetYaxis()->CenterTitle();
        projhh_1->GetYaxis()->SetTitleOffset(1.8);

	projhh_1->SetLineColor(kBlack);

	projhh_1->SetLineStyle(9);

	Double_t par[2];

	TCanvas* c_1 = new TCanvas("c_1","c_1",1000,1000);
        c_1->cd();
        TPad *pad_1 = new TPad("pad_1","pad_1",0,0,1,1);
        pad_1->SetLeftMargin(0.15);
        pad_1->Draw();
        pad_1->cd();
	gStyle->SetOptStat(0);
        gPad->SetGrid();
        projhh_1->Draw("HIST");

	TF1* fit_1 = new TF1("fit_1","expo",8,14);
	TF1* draw_1 = new TF1("draw_1","expo",0,14);
	projhh_1->Fit(fit_1,"R0");
	fit_1->GetParameters(&par[0]);
	draw_1->SetParameters(par);
	draw_1->SetLineStyle(10);
	draw_1->Draw("Lsame");

	TH1D* projhh_1_clone = projhh_1->Clone("projhh_1_clone");

	for(Int_t i = 1; i <= projhh_1->GetNbinsX(); i++)
	{
		projhh_1_clone->SetBinContent(i,projhh_1->GetBinContent(i)-draw_1->Eval(projhh_1_clone->GetBinCenter(i)));
	}
	
	projhh_1_clone->SetLineStyle(1);
	projhh_1_clone->SetFillStyle(3345);
	projhh_1_clone->SetFillColor(40);
	projhh_1_clone->Draw("HISTsame");

	TF1* fit_2 = new TF1("fit_2","gaus",3,5);
	projhh_1_clone->Fit(fit_2,"R0");

	Double_t minChX = fit_2->GetParameter(1) - 3.0*fit_2->GetParameter(2);
	Double_t maxChX = fit_2->GetParameter(1) + 3.0*fit_2->GetParameter(2);

	if(minChX < 0) 		minChX = 0;
	if(maxChX > 256*0.055) 	maxChX = 256*0.055;

	TLine* lineLeft 	= new TLine(minChX,0,minChX,projhh_1->GetMaximum());
	TLine* lineRight 	= new TLine(maxChX,0,maxChX,projhh_1->GetMaximum());

	lineLeft->SetLineStyle(7);
	lineRight->SetLineStyle(7);
	lineLeft->SetLineColor(kBlue);
	lineRight->SetLineColor(kBlue);
	lineLeft->Draw("L SAME");
	lineRight->Draw("L SAME");

	Double_t intCH, intCHerr;
	Double_t intALL, intALLerr;
	
	intCH 	= projhh_1_clone->IntegralAndError(projhh_1_clone->GetXaxis()->FindBin(minChX),projhh_1_clone->GetXaxis()->FindBin(maxChX),intCHerr);
	intALL 	= projhh_1_clone->IntegralAndError(1,projhh_1_clone->GetNbinsX(),intALLerr);

	cout<<endl;
	cout<<"--> intCH = "<<intCH<<" +/- "<<intCHerr<<endl;
	cout<<"--> intALL = "<<intALL<<" +/- "<<intALLerr<<endl;
	cout<<"--> intCH/intALL = "<<intCH/intALL<<" +/- "<<TMath::Sqrt(TMath::Power(intCHerr/intALL,2) + TMath::Power(intCH*intALLerr/(intALL*intALL),2))<<endl;
	

/*
	TCanvas* c_126 = new TCanvas("c_126","c_126",1000,1000);
        c_126->cd();
        TPad *pad_126 = new TPad("pad_126","pad_6",0,0,1,1);
        pad_126->SetLeftMargin(0.15);
        pad_126->Draw();
        pad_126->cd();
	gStyle->SetOptStat(0);
        gPad->SetGrid();

	projhh_1_clone->SetMinimum(0);
	projhh_2_clone->SetMinimum(0);
	projhh_6_clone->SetMinimum(0);

	projhh_1_clone->SetFillStyle(3254);
	projhh_2_clone->SetFillStyle(3205);
	projhh_6_clone->SetFillStyle(3290);

	projhh_1_clone->SetFillColor(kBlack);
	projhh_2_clone->SetFillColor(kBlue);
	projhh_6_clone->SetFillColor(kRed);

	projhh_1_clone->SetLineWidth(2);
	projhh_2_clone->SetLineWidth(2);
	projhh_6_clone->SetLineWidth(2);

	projhh_2_clone->Draw("hist same");
	projhh_6_clone->Draw("hist same");
	projhh_1_clone->Draw("hist same");

	TLegend* leg = new TLegend(0.45,0.65,0.85,0.85);
   	leg->AddEntry(projhh_1_clone,"W/o Crystal2","f");
   	leg->AddEntry(projhh_2_clone,"Crystal2 in AM","f");
   	leg->AddEntry(projhh_6_clone,"Crystal2 in CH","f");
   	leg->Draw();

//======================================================================================================================================================//


	TGraphErrors* gr_12 = new TGraphErrors();
	TGraphErrors* gr_16 = new TGraphErrors();

	for(Int_t i = 1; i <= projhh_1_clone->GetNbinsX(); i++)
	{
		gr_12->SetPoint(gr_12->GetN(),projhh_1_clone->GetBinCenter(i),projhh_2_clone->GetBinContent(i) - projhh_1_clone->GetBinContent(i));
		gr_12->SetPointError(gr_12->GetN()-1,projhh_1_clone->GetBinWidth(i)/TMath::Sqrt(12.0),TMath::Sqrt( TMath::Power(projhh_2_clone->GetBinError(i),2) + TMath::Power(projhh_1_clone->GetBinError(i),2) ));

		gr_16->SetPoint(gr_16->GetN(),projhh_1_clone->GetBinCenter(i),projhh_6_clone->GetBinContent(i) - projhh_1_clone->GetBinContent(i));
		gr_16->SetPointError(gr_16->GetN()-1,projhh_1_clone->GetBinWidth(i)/TMath::Sqrt(12.0),TMath::Sqrt( TMath::Power(projhh_6_clone->GetBinError(i),2) + TMath::Power(projhh_1_clone->GetBinError(i),2) ));
	}

	TCanvas* c_126_2 = new TCanvas("c_126_2","c_126_2",1900,1200);
        c_126_2->cd();
	gStyle->SetOptStat(0);
        gPad->SetGrid();


	gr_12->SetLineColor(kBlue);
	gr_16->SetLineColor(kRed);
	gr_12->SetLineStyle(10);
	gr_16->SetLineStyle(1);
	gr_12->SetMarkerStyle(21);
	gr_16->SetMarkerStyle(20);
	gr_12->SetMarkerSize(0.5);
	gr_16->SetMarkerSize(0.5);

	gr_12->SetMinimum(-0.02);
	gr_12->SetMaximum(0.02);
	gr_12->GetXaxis()->SetLimits(0,256*0.055);
	gr_12->Draw("APC");
	gr_16->Draw("PCSAME");

	gr_12->GetXaxis()->SetTitle("Horizontal Axis Projection [mm]");
        gr_12->GetYaxis()->SetTitle("Hits Density");
        gr_12->GetXaxis()->CenterTitle();
        gr_12->GetYaxis()->CenterTitle();
        gr_12->GetYaxis()->SetTitleOffset(1.2);

	TLegend* leg2 = new TLegend(0.45,0.15,0.85,0.35);
   	leg2->AddEntry(gr_12,"Crystal2 in AM - W/o Crystal2","ple");
   	leg2->AddEntry(gr_16,"Crystal2 in CH - W/o Crystal2","ple");
   	leg2->Draw();

	TPaveText *textCH = new TPaveText(8,0.009,9,0.012);
        textCH->SetTextSize(0.04);
        textCH->AddText("CH");
        textCH->Draw();

	TPaveText *textDCH = new TPaveText(6,0.002,7,0.005);
        textDCH->SetTextSize(0.04);
        textDCH->AddText("DCH");
        textDCH->Draw();

	TPaveText *textINI = new TPaveText(1.5,0.015,5,0.018);
        textINI->SetTextSize(0.04);
        textINI->AddText("Nuclear Interaction");
        textINI->Draw();
	
	TPaveText *textCS = new TPaveText(1.5,-0.018,5,-0.015);
        textCS->SetTextSize(0.04);
        textCS->AddText("Crystal Shadow");
        textCS->Draw();

/**/
    return 0;
}

void scale2Dhisto(TH2D* histo, Double_t integral, Double_t interal_err)// to scale the 2D histo
{
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= histo->GetNbinsY(); j++)
        {
            Double_t val = histo->GetBinContent(i,j);
            Double_t val_err = histo->GetBinError(i,j);
            val_err = TMath::Sqrt(TMath::Power(val_err/integral,2) + TMath::Power(val*interal_err/(integral*integral),2));
            val = val/integral;

            histo->SetBinContent(i,j,val);
            histo->SetBinError(i,j,val_err);
        }
    }
}

void scale1Dhisto(TH1D* histo, Double_t integral, Double_t integral_err)// to scale the 1D histo
{
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        Double_t val = histo->GetBinContent(i);
        Double_t val_err = histo->GetBinError(i);
        val_err = TMath::Sqrt(TMath::Power(val_err/integral,2) + TMath::Power(val*integral_err/(integral*integral),2));
        val = val/integral;

        histo->SetBinContent(i,val);
        histo->SetBinError(i,val_err);
    }
}

Double_t fitfdouble(Double_t *x,Double_t *par) 
{
        Double_t par_erf1[4],par_erf2[4];
        par_erf1[0] = par[0];
        par_erf1[1] = par[1];
        par_erf1[2] = par[2];
        par_erf1[3] = par[3];
        par_erf2[0] = par[4];
        par_erf2[1] = par[5];
        par_erf2[2] = par[6];
        par_erf2[3] = par[7];
        Double_t fitval = fitf(x,par_erf1) + fitf(x,par_erf2);
        return fitval;
}

Double_t fitf(Double_t *x,Double_t *par) 
{
        Double_t arg = 0;
        if (par[2]!=0) arg = (x[0] - par[1])/(par[2]*TMath::Sqrt(2.0));
        Double_t fitval = par[0]*TMath::Erf(arg) + par[3];
        return fitval;
}

