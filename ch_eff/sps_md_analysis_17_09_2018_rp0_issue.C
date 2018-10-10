//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     Data analysis for SPS MD data
//--------------------------------------------------------------------------//

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

using namespace std;

const Int_t N_PIXELS                = 512;
const Double_t PIXEL_SIZE           = 0.055;

void removenoizypixelsY(TH1D* histo, Float_t factor_bad_pixels);
void removenoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels);
void scale1Dhisto(TH1D* histo, Double_t integral, Double_t interal_err);
void scale2Dhisto(TH2D* histo, Double_t integral, Double_t interal_err);
void definearea2D(TH2D* histo);
double fit_ch_dch(Double_t *x,Double_t *par);

int sps_md_analysis_17_09_2018_rp0_issue()
{
    cout<<"--> function_1() -- to plot rp0 profile during angular scan normalized by max"<<endl;
    return 0;
}

int function_1()
{
    const Int_t bin_ini = 1;
    const Int_t bin_fin = 4577;

    gStyle->SetOptStat(0);
    TFile *_file_rp0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/MD_2018_09_17_K6_HISTO_RP0I.root");
    TH2D* h_rp0 = (TH2D*)_file_rp0->Get("h_8");
    TH1D* h_y = (TH1D*)h_rp0->ProjectionY("h_y",bin_ini,bin_fin);
    removenoizypixelsY(h_y,35);
    for(Int_t j = 1; j <= h_y->GetNbinsX(); j++)
    {
        if(j == 208 || j == 246 || j == 209)
        {
            h_y->SetBinContent(j,0);
            h_y->SetBinError(j,0);
        }
    }

    Double_t integral_err = 0;
    Double_t integral = h_y->IntegralAndError(1,h_y->GetNbinsX(),integral_err);
    scale1Dhisto(h_y,integral,integral_err);

    TH1D* h_y_mm = new TH1D("h_y_mm","Projection on Y axis",h_y->GetNbinsX(),0.055,h_y->GetNbinsX()*0.055);
    for(Int_t i = 1; i <= h_y_mm->GetNbinsX(); i++)
    {
        h_y_mm->SetBinContent(i,h_y->GetBinContent(i));
        h_y_mm->SetBinError(i,h_y->GetBinError(i));
    }

    TCanvas* c_1 = new TCanvas("c_1","c_1",1200,600);
    c_1->cd();
    gPad->SetGrid();
    h_y_mm->GetXaxis()->SetRange(h_y_mm->FindBin(0),h_y_mm->FindBin(256*0.055));
    h_y_mm->SetLineWidth(2);
    h_y_mm->GetYaxis()->SetTitleOffset(2.0);
    h_y_mm->GetXaxis()-> CenterTitle();
    h_y_mm->Draw("hist");

    TF1 *f1 = new TF1("f1","x",61.98-0.055*256,61.98);
    TGaxis *axis1 = new TGaxis(0,-0.15e-2,0.055*256,-0.15e-2,"f1",510,"");
    axis1->SetLineColor(kRed);
    axis1->Draw();

    const Int_t nHist = bin_fin-bin_ini+1;
    TH1D* h_y_i[nHist];
    TH1D* h_y_i_mm[nHist];
    cout<<"--> nHist = "<<nHist<<endl;
/*
    TCanvas* c_2 = new TCanvas("c_2","c_2",1200,600);
    c_2->cd();
    gPad->SetGrid();

    for(Int_t i = bin_ini; i <= bin_fin; i+=2)
    {
        TString name = "h_y_i_";
        name += i;
        h_y_i[i] = (TH1D*)h_rp0->ProjectionY(name.Data(),i,i+1);
        removenoizypixelsY(h_y_i[i],1);
        for(Int_t j = 1; j <= h_y_i[i]->GetNbinsX(); j++)
        {
            if(j == 208 || j == 246 || j == 209)
            {
                h_y_i[i]->SetBinContent(j,0);
                h_y_i[i]->SetBinError(j,0);
            }
        }

        integral = h_y_i[i]->IntegralAndError(1,h_y_i[i]->GetNbinsX(),integral_err);
        scale1Dhisto(h_y_i[i],integral,integral_err);

        name = "h_y_i_mm_";
        name += i;
        h_y_i_mm[i] = new TH1D(name.Data(),"Projection on Y axis",h_y->GetNbinsX(),0.055,h_y->GetNbinsX()*0.055);
        for(Int_t j = 1; j <= h_y_i_mm[i]->GetNbinsX(); j++)
        {
            h_y_i_mm[i]->SetBinContent(j,h_y_i[i]->GetBinContent(j));
            h_y_i_mm[i]->SetBinError(j,h_y_i[i]->GetBinError(j));
        }

        h_y_i_mm[i]->GetXaxis()->SetRange(h_y_i_mm[i]->FindBin(0),h_y_i_mm[i]->FindBin(256*0.055));
        h_y_i_mm[i]->SetMaximum(6e-2);
        h_y_i_mm[i]->SetMinimum(-0.5e-2);
        h_y_i_mm[i]->SetLineWidth(2);
        h_y_i_mm[i]->GetYaxis()->SetTitleOffset(2.0);
        h_y_i_mm[i]->GetXaxis()-> CenterTitle();
//        h_y_i_mm[i]->SetLineColor(1+i-bin_ini);
        name = "Time: ";
        name += (Int_t)(h_rp0->GetXaxis()->GetBinCenter(i) - h_rp0->GetXaxis()->GetBinCenter(bin_ini));
        name += " [sec]";
        h_y_i_mm[i]->SetTitle(name.Data());
        h_y_i_mm[i]->Draw("hist & same");
    }
//    c_2->BuildLegend();
    TGaxis *axis2 = new TGaxis(0,-2.5e-3,0.055*256,-2.5e-3,"f1",510,"");
    axis2->SetLineColor(kRed);
    axis2->Draw();
*/
    return 0;
}

void removenoizypixelsY(TH1D* histo, Float_t factor_bad_pixels)// to remove noizy pixels for Y projection
{
    Int_t nMaskedPixels = 0;
    for(Int_t xi = 1; xi < histo->GetNbinsX(); xi++)
    {
        if(histo->GetBinContent(xi) - factor_bad_pixels*histo->GetBinError(xi) >
                histo->GetBinContent(xi+1) + factor_bad_pixels*histo->GetBinError(xi+1))
        {
            histo->SetBinContent(xi,0);
            histo->SetBinError(xi,0);
            nMaskedPixels++;
        }
    }
//    cout<<"--> nMaskedPixels = "<<nMaskedPixels<<endl;
}

void removenoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels)// to remove noizy pixels for XY image
{
    Int_t nMaskedPixels = 0;
    for(Int_t xi = 1; xi < histo->GetNbinsX(); xi++)
    {
        for(Int_t yi = 1; yi < histo->GetNbinsX(); yi++)
        {
            if(histo->GetBinContent(xi,yi) - factor_bad_pixels*histo->GetBinError(xi,yi) >
                    histo->GetBinContent(xi,yi+1) + factor_bad_pixels*histo->GetBinError(xi,yi+1))
            {
                histo->SetBinContent(xi,yi,0);
                histo->SetBinError(xi,yi,0);
                nMaskedPixels++;
            }
        }
    }
//    cout<<"--> nMaskedPixels = "<<nMaskedPixels<<endl;
}

void scale1Dhisto(TH1D* histo, Double_t integral, Double_t interal_err)// to scale the 1D histo
{
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        Double_t val = histo->GetBinContent(i);
        Double_t val_err = histo->GetBinError(i);
        val_err = TMath::Sqrt(TMath::Power(val_err/integral,2) + TMath::Power(val*interal_err/(integral*integral),2));
        val = val/integral;

        histo->SetBinContent(i,val);
        histo->SetBinError(i,val_err);
    }
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

void definearea2D(TH2D* histo)// to remove some area
{
    for(Int_t xi = 1; xi <= histo->GetNbinsX(); xi++)
    {
        for(Int_t yi = 1; yi <= histo->GetNbinsX(); yi++)
        {
            if(histo->GetYaxis()->GetBinCenter(yi) < 9.7-1 || histo->GetYaxis()->GetBinCenter(yi) > 9.7+1)
            {
                histo->SetBinContent(xi,yi,0);
                histo->SetBinError(xi,yi,0);
            }
        }
    }
}

double fit_ch_dch(Double_t *x,Double_t *par)
{
    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    Double_t fitval;
    if(x[0] < par[1]-0.5*par[2])
        fitval = par[0]*TMath::Exp(-0.5*arg*arg) + TMath::Exp(par[3] + par[4]*x[0]);
    else
        fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}
