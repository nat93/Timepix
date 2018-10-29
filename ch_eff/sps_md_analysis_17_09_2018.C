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
#include <algorithm>
#include <iostream>

using namespace std;

const Int_t N_PIXELS                = 512;
const Double_t PIXEL_SIZE           = 0.055;

void removenoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels);
void scale2Dhisto(TH2D* histo, Double_t integral, Double_t interal_err);
void definearea2D(TH2D* histo);
double fit_ch_dch(Double_t *x,Double_t *par);
void median_filter(TH2D* h_in, TH2D* h_out, Double_t factor_bad_pixels);
double findMedian(double a[], int n);

int sps_md_analysis_17_09_2018()
{
    cout<<"--> function_1() -- to plot rp1 profile normalized by rp0"<<endl;
    cout<<"--> function_2() -- to plot rp1 profile normalized by time"<<endl;
    cout<<"--> function_3(Int_t i) -- to fit the beam profile at rp1 (bkg fit expo)"<<endl;
    cout<<"--> function_4(Int_t i) -- to fit the beam profile at rp1 (bkg fit pol1)"<<endl;
    cout<<"--> function_5(Int_t i) -- to fit the beam profile at rp1 (no bkg fit)"<<endl;
    cout<<"--> function_6() -- to plot the beam profile with diff. rp1/0 position"<<endl;
    cout<<"--> function_7() -- to calculate the ratio CH/DCH [mm] "<<endl;
    cout<<"--> function_8(Int_t i) -- to fit the beam profile at rp1 (bkg fit expo)"<<endl;
    return 0;
}

int function_1()
{
    const Int_t nSets = 13;
    TString fileName_RP1[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP1I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP1I_RUN_6.root"
    };

    TString fileName_RP0[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP0I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP0I_RUN_6.root"
    };

    Double_t fit_lim_min[] = {7,9,9,9,9,9,9,0,0,7,9,9,9};
    Double_t fit_lim_max[] = {9,11,11,11,11,11,11,14,14,9,11,11,11};
    Float_t factor_bad_pixels[] = {200,200,200,200,200,30,4,100,30,200,200,100,50};

    TCanvas* c_1[nSets];
    TH1D* h_rp1_x[nSets];
    TH2D* h_rp0[nSets];
    TH2D* h_rp1[nSets];
    TString name;

    gStyle->SetOptStat(0);
    for(Int_t i = 0; i < nSets; i++)
    {
        cout<<endl<<"--> Set "<<i+1<<" <--"<<endl;
        TFile *_file_rp0 = TFile::Open(fileName_RP0[i].Data());
        TFile *_file_rp1 = TFile::Open(fileName_RP1[i].Data());

        h_rp0[i] = (TH2D*)_file_rp0->Get("h_6");
        h_rp1[i] = (TH2D*)_file_rp1->Get("h_6");

        removenoizypixelsXY(h_rp0[i],factor_bad_pixels[i]);
        removenoizypixelsXY(h_rp1[i],factor_bad_pixels[i]);

        Double_t integral_err_rp0, integral_err_rp1;
        Double_t integral_rp0 = h_rp0[i]->IntegralAndError(1,h_rp0[i]->GetNbinsX(),1,h_rp0[i]->GetNbinsY(),integral_err_rp0);
        Double_t integral_rp1 = h_rp1[i]->IntegralAndError(1,h_rp1[i]->GetNbinsX(),1,h_rp1[i]->GetNbinsY(),integral_err_rp1);
        cout<<"--> Integral RP1/RP0 = "<<integral_rp1/integral_rp0<<endl;
        scale2Dhisto(h_rp0[i], integral_rp0, integral_err_rp0);
        scale2Dhisto(h_rp1[i], integral_rp0, integral_err_rp0);

        name = "h_rp1_x";
        name += i+1;
        h_rp1_x[i] = h_rp1[i]->ProjectionX(name.Data());
        h_rp1_x[i]->SetLineWidth(2);
        TF1* fit_funct = new TF1("fit_funct","gaus",fit_lim_min[i],fit_lim_max[i]);

        name = "c_1_s";
        name += i+1;
        c_1[i] = new TCanvas(name.Data(),name.Data(),500,1000);
        c_1[i]->Divide(1,2);

        c_1[i]->cd(1);
        gPad->SetGrid();
        name = "RP1 (S";
        name += i+1;
        name += ") Normalized";
        h_rp1[i]->SetTitle(name.Data());
        h_rp1[i]->GetXaxis()->SetRange(1,h_rp1[i]->GetNbinsX()/2);
        h_rp1[i]->GetYaxis()->SetRange(1,h_rp1[i]->GetNbinsY()/2);
        h_rp1[i]->GetXaxis()-> CenterTitle();
        h_rp1[i]->GetYaxis()-> CenterTitle();
        h_rp1[i]->Draw("colz");

        c_1[i]->cd(2);
        gPad->SetGrid();
        name = "RP1 Projection on X axis (S";
        name += i+1;
        name += ")";
        h_rp1_x[i]->GetXaxis()->SetRange(1,h_rp1_x[i]->GetNbinsX()/2);
        h_rp1_x[i]->GetXaxis()-> CenterTitle();
        h_rp1_x[i]->GetYaxis()-> CenterTitle();
        h_rp1_x[i]->SetTitle(name.Data());
        h_rp1_x[i]->Draw("hist");

        if(i+1 != 8 && i+1 != 9)
        {
            h_rp1_x[i]->Fit(fit_funct,"R0Q+");
            fit_funct->Draw("same");
            Double_t minX = fit_funct->GetParameter(1)-3.0*fit_funct->GetParameter(2);
            Double_t meanX = fit_funct->GetParameter(1);
            Double_t maxX = fit_funct->GetParameter(1)+3.0*fit_funct->GetParameter(2);
            Double_t minY = h_rp1_x[i]->GetMinimum();
            Double_t maxY = h_rp1_x[i]->GetMaximum();
            TLine *l_min=new TLine(minX,minY,minX,maxY);
            l_min->SetLineColor(kRed);
            l_min->SetLineWidth(3);
            l_min->SetLineStyle(7);
            l_min->Draw();
            TLine *l_mean=new TLine(meanX,minY,meanX,maxY);
            l_mean->SetLineColor(kBlack);
            l_mean->SetLineWidth(3);
            l_mean->SetLineStyle(7);
            l_mean->Draw();
            TLine *l_max=new TLine(maxX,minY,maxX,maxY);
            l_max->SetLineColor(kGreen+2);
            l_max->SetLineWidth(3);
            l_max->SetLineStyle(7);
            l_max->Draw();
        }

        name = "c_rp1_rp0_s";
        name += i+1;
        name += ".png";
        c_1[i]->SaveAs(name.Data());
    }

    return 0;
}

int function_2()
{
    const Int_t nSets = 13;
    TString fileName_RP1[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP1I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP1I_RUN_6.root"
    };

    Double_t fit_lim_min[] = {7,9,9,9,9,9,9,0,0,7,9,9,9};
    Double_t fit_lim_max[] = {9,11,11,11,11,11,11,14,14,9,11,11,11};
//    Float_t factor_bad_pixels[] = {200,200,200,200,200,30,4,100,30,200,200,100,50};

    TCanvas* c_1[nSets];
    TH1D* h_rp1_x[nSets];
    TH2D* h_rp0[nSets];
    TH2D* h_rp1[nSets];
    TString name;

    gStyle->SetOptStat(0);
    for(Int_t i = 0; i < nSets; i++)
    {
        cout<<endl<<"--> Set "<<i+1<<" <--"<<endl;
        TFile *_file_rp1 = TFile::Open(fileName_RP1[i].Data());

        h_rp0[i] = (TH2D*)_file_rp1->Get("h_6");
        h_rp1[i] = (TH2D*)_file_rp1->Get("h_6");
        TH2D* h_rp1_t = (TH2D*)_file_rp1->Get("h_5");
        Double_t dtime = (h_rp1_t->GetXaxis()->GetBinCenter(h_rp1_t->GetNbinsX()) - h_rp1_t->GetXaxis()->GetBinCenter(1))*0.1;

        median_filter(h_rp0[i],h_rp1[i],1.2);
//        removenoizypixelsXY(h_rp1[i],factor_bad_pixels[i]);

        // To normalize by RP1
//        removenoizypixelsXY(h_rp0[i],factor_bad_pixels[i]);
//        Double_t integral_err_rp0, integral_err_rp1;
//        Double_t integral_rp0 = h_rp0[i]->IntegralAndError(1,h_rp0[i]->GetNbinsX(),1,h_rp0[i]->GetNbinsY(),integral_err_rp0);
//        Double_t integral_rp1 = h_rp1[i]->IntegralAndError(1,h_rp1[i]->GetNbinsX(),1,h_rp1[i]->GetNbinsY(),integral_err_rp1);
//        cout<<"--> Integral RP1/RP0 = "<<integral_rp1/integral_rp0<<endl;
//        scale2Dhisto(h_rp0[i], integral_rp0, integral_err_rp0);
//        scale2Dhisto(h_rp1[i], integral_rp0, integral_err_rp0);

        // To normalize by time
        cout<<"--> dTime = "<<dtime<<" [sec]"<<endl;
        scale2Dhisto(h_rp1[i], dtime, 0.0);

        name = "h_rp1_x";
        name += i+1;
        h_rp1_x[i] = h_rp1[i]->ProjectionX(name.Data());
        h_rp1_x[i]->SetLineWidth(2);
        TF1* fit_funct = new TF1("fit_funct","gaus",fit_lim_min[i],fit_lim_max[i]);

        name = "c_1_s";
        name += i+1;
        c_1[i] = new TCanvas(name.Data(),name.Data(),500,1000);
        c_1[i]->Divide(1,2);

        c_1[i]->cd(1);
        gPad->SetGrid();
        name = "RP1 (S";
        name += i+1;
        name += ") Normalized";
        h_rp1[i]->SetTitle(name.Data());
        h_rp1[i]->GetXaxis()->SetRange(1,h_rp1[i]->GetNbinsX()/2);
        h_rp1[i]->GetYaxis()->SetRange(1,h_rp1[i]->GetNbinsY()/2);
        h_rp1[i]->GetXaxis()-> CenterTitle();
        h_rp1[i]->GetYaxis()-> CenterTitle();
        h_rp1[i]->Draw("colz");

        c_1[i]->cd(2);
        gPad->SetGrid();
        name = "RP1 Projection on X axis (S";
        name += i+1;
        name += ")";
        h_rp1_x[i]->GetXaxis()->SetRange(1,h_rp1_x[i]->GetNbinsX()/2);
        h_rp1_x[i]->GetXaxis()-> CenterTitle();
        h_rp1_x[i]->GetYaxis()-> CenterTitle();
        h_rp1_x[i]->GetYaxis()->SetTimeOffset(1.2);
        h_rp1_x[i]->SetTitle(name.Data());        
        h_rp1_x[i]->Draw("hist");

        if(i+1 != 8 && i+1 != 9)
        {
            h_rp1_x[i]->Fit(fit_funct,"R0Q+");
            fit_funct->Draw("same");
            Double_t minX = fit_funct->GetParameter(1)-3.0*fit_funct->GetParameter(2);
            Double_t meanX = fit_funct->GetParameter(1);
            Double_t maxX = fit_funct->GetParameter(1)+3.0*fit_funct->GetParameter(2);
            Double_t minY = h_rp1_x[i]->GetMinimum();
            Double_t maxY = h_rp1_x[i]->GetMaximum();
            TLine *l_min=new TLine(minX,minY,minX,maxY);
            l_min->SetLineColor(kRed);
            l_min->SetLineWidth(3);
            l_min->SetLineStyle(7);
            l_min->Draw();
            TLine *l_mean=new TLine(meanX,minY,meanX,maxY);
            l_mean->SetLineColor(kBlack);
            l_mean->SetLineWidth(3);
            l_mean->SetLineStyle(7);
            l_mean->Draw();
            TLine *l_max=new TLine(maxX,minY,maxX,maxY);
            l_max->SetLineColor(kGreen+2);
            l_max->SetLineWidth(3);
            l_max->SetLineStyle(7);
            l_max->Draw();
        }

        name = "c_rp1_rp1_s";
        name += i+1;
        name += ".png";
        c_1[i]->SaveAs(name.Data());
    }

    return 0;
}

int function_3(Int_t i)
{
    const Int_t nSets = 13;

    Double_t mg_max[] = {/*0*/110.0,/*1*/110.0,/*2*/110.0,/*3*/110.0,/*4*/110.0,/*5*/110.0,/*6*/500.0,/*7*/500.0,/*8*/500.0,/*9*/70.00,/*10*/70.00,/*11*/70.00,/*12*/70.00};
    Double_t mg_min[] = {/*0*/20.00,/*1*/20.00,/*2*/20.00,/*3*/20.00,/*4*/20.00,/*5*/20.00,/*6*/0.000,/*7*/0.000,/*8*/0.000,/*9*/20.00,/*10*/20.00,/*11*/20.00,/*12*/20.00};

    TString fileName_RP0[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP0I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP0I_RUN_6.root"
    };

    TString fileName_RP1[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP1I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP1I_RUN_6.root"
    };

    Double_t fit_ch_lim_min[]   = {/*0*/7.00,/*1*/9.00,/*2*/9.00,/*3*/9.00,/*4*/9.00,/*5*/9.00,/*6*/9.00,/*7*/0.00,/*8*/0.00,/*9*/7.00,/*10*/9.00,/*11*/9.00,/*12*/9.00};
    Double_t fit_ch_lim_max[]   = {/*0*/9.00,/*1*/11.0,/*2*/11.0,/*3*/11.0,/*4*/11.0,/*5*/11.0,/*6*/11.0,/*7*/14.0,/*8*/14.0,/*9*/9.00,/*10*/11.0,/*11*/11.0,/*12*/11.0};
    Double_t fit_exp_lim_min[]  = {/*0*/12.0,/*1*/13.0,/*2*/13.0,/*3*/13.0,/*4*/13.0,/*5*/13.0,/*6*/13.0,/*7*/0.00,/*8*/0.00,/*9*/13.0,/*10*/13.0,/*11*/13.0,/*12*/13.0};
    Double_t fit_exp_lim_max[]  = {/*0*/13.8,/*1*/13.8,/*2*/13.8,/*3*/13.8,/*4*/13.8,/*5*/13.8,/*6*/13.8,/*7*/0.00,/*8*/0.00,/*9*/13.8,/*10*/13.8,/*11*/13.8,/*12*/13.8};
    Int_t expo_fit_bkg_status[] = {/*0*/1,   /*1*/1,   /*2*/0,   /*3*/1,   /*4*/1,   /*5*/1,   /*6*/1,   /*7*/0,   /*8*/0,   /*9*/1,   /*10*/1,   /*11*/1,   /*12*/1};       // 0 -- fit gauss+expo ; 1 -- fit expo ; -1 -- w/o fit
    Double_t fit_dch_lim_min[]  = {/*0*/2.50,/*1*/4.50,/*2*/4.50,/*3*/4.50,/*4*/4.50,/*5*/4.50,/*6*/0.50,/*7*/2.50,/*8*/2.50,/*9*/2.50,/*10*/4.50,/*11*/4.50,/*12*/4.50};
    Double_t fit_dch_lim_max[]  = {/*0*/4.00,/*1*/6.50,/*2*/6.50,/*3*/6.50,/*4*/6.50,/*5*/6.50,/*6*/6.50,/*7*/4.00,/*8*/4.00,/*9*/4.00,/*10*/6.50,/*11*/6.50,/*12*/6.50};
    Float_t factor_bad_pixels[] = {/*0*/200.,/*1*/200.,/*2*/200.,/*3*/200.,/*4*/200.,/*5*/30.0,/*6*/4.00,/*7*/100.,/*8*/30.0,/*9*/200.,/*10*/200.,/*11*/100.,/*12*/50.0};
    Double_t par_ch_bkg[nSets][5];

    // Gauss
    par_ch_bkg[0][0] =  0.60; par_ch_bkg[1][0] =  0.60; par_ch_bkg[2][0] =  0.06; par_ch_bkg[3][0] =  0.10;
    par_ch_bkg[4][0] =  0.60; par_ch_bkg[5][0] =  0.60; par_ch_bkg[6][0] =  0.60; par_ch_bkg[7][0] =  0.60;
    par_ch_bkg[8][0] =  0.60; par_ch_bkg[9][0] =  0.60; par_ch_bkg[10][0] = 0.60; par_ch_bkg[11][0] = 0.40;
    par_ch_bkg[12][0] = 0.60;

    par_ch_bkg[0][1] =  7.50; par_ch_bkg[1][1] =  7.50; par_ch_bkg[2][1] =  10.0; par_ch_bkg[3][1] =  7.50;
    par_ch_bkg[4][1] =  7.50; par_ch_bkg[5][1] =  7.50; par_ch_bkg[6][1] =  7.50; par_ch_bkg[7][1] =  10.00;
    par_ch_bkg[8][1] =  7.50; par_ch_bkg[9][1] =  7.50; par_ch_bkg[10][1] = 7.50; par_ch_bkg[11][1] = 10.00;
    par_ch_bkg[12][1] = 7.50;

    par_ch_bkg[0][2] =  0.80; par_ch_bkg[1][2] =  0.80; par_ch_bkg[2][2] =  0.70; par_ch_bkg[3][2] =  0.70;
    par_ch_bkg[4][2] =  0.80; par_ch_bkg[5][2] =  0.80; par_ch_bkg[6][2] =  0.80; par_ch_bkg[7][2] =  0.80;
    par_ch_bkg[8][2] =  0.80; par_ch_bkg[9][2] =  0.80; par_ch_bkg[10][2] = 0.80; par_ch_bkg[11][2] = 1.40;
    par_ch_bkg[12][2] = 0.80;

    // Expo
    par_ch_bkg[0][3] = -4.40; par_ch_bkg[1][3] = -3.70; par_ch_bkg[2][3] = -4.16; par_ch_bkg[3][3] = -3.80;
    par_ch_bkg[4][3] = -2.50; par_ch_bkg[5][3] = -6.00; par_ch_bkg[6][3] = -2.50; par_ch_bkg[7][3] = -2.50;
    par_ch_bkg[8][3] = -2.50; par_ch_bkg[9][3] = -2.50; par_ch_bkg[10][3] =-2.50; par_ch_bkg[11][3] =-1.60;
    par_ch_bkg[12][3] =-2.50;

    par_ch_bkg[0][4] = -0.20; par_ch_bkg[1][4] = -0.20; par_ch_bkg[2][4] = -0.20; par_ch_bkg[3][4] = -0.17;
    par_ch_bkg[4][4] = -0.06; par_ch_bkg[5][4] = -0.02; par_ch_bkg[6][4] = -0.06; par_ch_bkg[7][4] = -0.06;
    par_ch_bkg[8][4] = -0.06; par_ch_bkg[9][4] = -0.06; par_ch_bkg[10][4] =-0.06; par_ch_bkg[11][4] =-0.01;
    par_ch_bkg[12][4] =-0.06;


    TH2D* h_rp0[nSets];
    TH2D* h_rp1[nSets];

    gStyle->SetOptStat(0);

    cout<<endl<<"--> Set "<<i+1<<" <--"<<endl;
    TFile *_file_rp0 = TFile::Open(fileName_RP0[i].Data());
    TFile *_file_rp1 = TFile::Open(fileName_RP1[i].Data());

    h_rp0[i] = (TH2D*)_file_rp0->Get("h_6");
    h_rp1[i] = (TH2D*)_file_rp1->Get("h_6");

    removenoizypixelsXY(h_rp0[i],factor_bad_pixels[i]);
    removenoizypixelsXY(h_rp1[i],factor_bad_pixels[i]);

    Double_t integral_err_rp0;
    Double_t integral_rp0 = h_rp0[i]->IntegralAndError(1,h_rp0[i]->GetNbinsX(),1,h_rp0[i]->GetNbinsY(),integral_err_rp0);
    scale2Dhisto(h_rp0[i], integral_rp0, integral_err_rp0);
    scale2Dhisto(h_rp1[i], integral_rp0, integral_err_rp0);
    cout<<"--> Integral RP1: "<<h_rp1[i]->Integral()<<endl;

    h_rp1[i]->GetXaxis()->SetRange(1,h_rp1[i]->GetNbinsX()/2);
    h_rp1[i]->GetYaxis()->SetRange(1,h_rp1[i]->GetNbinsY()/2);
    h_rp1[i]->GetXaxis()-> CenterTitle();
    h_rp1[i]->GetYaxis()-> CenterTitle();

    TString name = "c_1_";
    name += i;
    TCanvas* c_1 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_1->cd();
    h_rp1[i]->Draw("colz");

    name = "h_rp1_x_";
    name += i+1;
    TH1D* h_rp1_x = h_rp1[i]->ProjectionX(name.Data());
    name += "_s";
    TH1D* h_rp1_x_s = h_rp1[i]->ProjectionX(name.Data());
    name += "_s";
    TH1D* h_rp1_x_s_s = h_rp1[i]->ProjectionX(name.Data());

    name = "c_2_";
    name += i;
    TCanvas* c_2 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_2->cd();
    h_rp1_x->SetMinimum(0);
    h_rp1_x->Draw("hist");
    TF1* fit_funct_ch_bkg;

    if(expo_fit_bkg_status[i] == 1)
    {
        fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","expo",fit_exp_lim_min[i],fit_exp_lim_max[i]);
        fit_funct_ch_bkg->SetParameter(0,par_ch_bkg[i][3]);
        fit_funct_ch_bkg->SetParameter(1,par_ch_bkg[i][4]);
        h_rp1_x->Fit(fit_funct_ch_bkg,"R+");
        par_ch_bkg[i][3] = fit_funct_ch_bkg->GetParameter(0);
        par_ch_bkg[i][4] = fit_funct_ch_bkg->GetParameter(1);
    }
    else
    {
        fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","gaus(0)+expo(3)",fit_ch_lim_min[i],fit_exp_lim_max[i]);
        fit_funct_ch_bkg->SetParameters(par_ch_bkg[i]);
        h_rp1_x->Fit(fit_funct_ch_bkg,"R+");
        fit_funct_ch_bkg->GetParameters(par_ch_bkg[i]);
    }

    name = "#chi^{2}/NDF = [";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare();
    name += "/";
    name += fit_funct_ch_bkg->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare()/fit_funct_ch_bkg->GetNDF();
    h_rp1_x->SetTitle(name.Data());
    fit_funct_ch_bkg->SetLineWidth(6);
    if(expo_fit_bkg_status[i] >= 0) fit_funct_ch_bkg->Draw("same");
    TF1* bkg_func = new TF1("bkg_func","expo",h_rp1_x->GetBinCenter(1),h_rp1_x->GetBinCenter(h_rp1_x->GetNbinsX()));
    bkg_func->SetParameter(0,par_ch_bkg[i][3]);
    bkg_func->SetParameter(1,par_ch_bkg[i][4]);
    bkg_func->SetLineColor(kGreen+2);
    if(expo_fit_bkg_status[i] >= 0) bkg_func->Draw("same");

    name = "c_3_";
    name += i;
    TCanvas* c_3 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_3->cd();
    h_rp1_x_s->SetMinimum(0);
    h_rp1_x_s->SetLineColor(kBlue);
    h_rp1_x_s->SetLineWidth(2);
    h_rp1_x_s->Draw("hist");
    for(Int_t j = 1; j <= h_rp1_x_s_s->GetNbinsX(); j++)
    {
        Double_t val;
        if(expo_fit_bkg_status[i] >= 0)
            val = h_rp1_x_s->GetBinContent(j) - bkg_func->Eval(h_rp1_x_s->GetBinCenter(j));
        else
            val = h_rp1_x_s->GetBinContent(j);



        if(val > 0)
            h_rp1_x_s_s->SetBinContent(j,val);
        else
            h_rp1_x_s_s->SetBinContent(j,0);
    }
    h_rp1_x_s_s->SetLineColor(kRed);
    h_rp1_x_s_s->SetLineWidth(2);
    h_rp1_x_s_s->Draw("same & hist");

    name = "c_4_";
    name += i;
    TCanvas* c_4 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_4->cd();
    h_rp1_x_s_s->SetLineWidth(8);
    h_rp1_x_s_s->Draw("hist");

    TF1* fit_funct_1 = new TF1("fit_funct_1","gaus",fit_ch_lim_min[i],fit_ch_lim_max[i]);
    h_rp1_x_s_s->Fit(fit_funct_1,"R0Q+");
    Double_t par_ch[3];
    fit_funct_1->GetParameters(par_ch);
    TF1* fit_funct_ch = new TF1("fit_funct_ch","gaus",h_rp1_x_s_s->GetBinCenter(1),h_rp1_x_s_s->GetBinCenter(h_rp1_x_s_s->GetNbinsX()));
    fit_funct_ch->SetParameters(par_ch);
    fit_funct_ch->SetLineColor(kBlue);
    fit_funct_ch->SetLineWidth(2);
//    fit_funct_ch->Draw("same");

    TF1* fit_funct_2 = new TF1("fit_funct_2","expo",fit_dch_lim_min[i],fit_dch_lim_max[i]);
    h_rp1_x_s_s->Fit(fit_funct_2,"R0Q+");
    Double_t par_dch[2];
    fit_funct_2->GetParameters(par_dch);
    TF1* fit_funct_dch = new TF1("fit_funct_dch","expo",h_rp1_x_s_s->GetBinCenter(1),h_rp1_x_s_s->GetBinCenter(h_rp1_x_s_s->GetNbinsX()));
    fit_funct_dch->SetParameters(par_dch);
    fit_funct_dch->SetLineColor(kBlack);
    fit_funct_dch->SetLineWidth(2);
//    fit_funct_dch->Draw("same");

    TLine *l_min_ch=new TLine(fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_min_ch->SetLineColor(kMagenta);
    l_min_ch->SetLineWidth(3);
    l_min_ch->SetLineStyle(7);
    l_min_ch->Draw();
    TLine *l_mean_ch=new TLine(fit_funct_ch->GetParameter(1),0,fit_funct_ch->GetParameter(1),h_rp1_x_s_s->GetMaximum());
    l_mean_ch->SetLineColor(kMagenta-2);
    l_mean_ch->SetLineWidth(3);
    l_mean_ch->SetLineStyle(7);
    l_mean_ch->Draw();
    TLine *l_max_ch=new TLine(fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_max_ch->SetLineColor(kMagenta);
    l_max_ch->SetLineWidth(3);
    l_max_ch->SetLineStyle(7);
    l_max_ch->Draw();

    TLine *l_min_dch=new TLine(fit_dch_lim_min[i],0,fit_dch_lim_min[i],h_rp1_x_s_s->GetMaximum());
    l_min_dch->SetLineColor(kGreen+1);
    l_min_dch->SetLineWidth(3);
    l_min_dch->SetLineStyle(7);
    l_min_dch->Draw();
    TLine *l_max_dch=new TLine(fit_dch_lim_max[i],0,fit_dch_lim_max[i],h_rp1_x_s_s->GetMaximum());
    l_max_dch->SetLineColor(kGreen+1);
    l_max_dch->SetLineWidth(3);
    l_max_dch->SetLineStyle(7);
    l_max_dch->Draw();

    TF1* fit_funct_ch_dch = new TF1("fit_funct_ch_dch",fit_ch_dch,fit_dch_lim_min[i],fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),5);
    fit_funct_ch_dch->SetParameter(0,par_ch[0]);
    fit_funct_ch_dch->SetParameter(1,par_ch[1]);
    fit_funct_ch_dch->SetParameter(2,par_ch[2]);
    fit_funct_ch_dch->SetParameter(3,par_dch[0]);
    fit_funct_ch_dch->SetParameter(4,par_dch[1]);
    fit_funct_ch_dch->SetLineColor(kBlack);
    fit_funct_ch_dch->SetLineWidth(4);
    fit_funct_ch_dch->Draw("same");

    name = "#chi^{2}/NDF (DCH) = [";
    name += (Int_t)fit_funct_2->GetChisquare();
    name += "/";
    name += fit_funct_2->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_2->GetChisquare()/fit_funct_2->GetNDF();
    name += " | #chi^{2}/NDF (CH) = [";
    name += (Int_t)fit_funct_1->GetChisquare();
    name += "/";
    name += fit_funct_1->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_1->GetChisquare()/fit_funct_1->GetNDF();
    h_rp1_x_s_s->SetTitle(name.Data());

    gPad->Modified();
    name = "profile_ch_dch_";
    name += i;
    name += ".png";
    c_4->SaveAs(name.Data());

    ofstream outputASCII_histo;
    name = "ASCII_DATA_PROFILE_METHOD3_S";
    name += i+1;
    name += ".dat";
    outputASCII_histo.open (name.Data());
    for(Int_t bini = 1; bini <= h_rp1_x_s_s->GetNbinsX(); bini++)
    {
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinCenter(bini)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinContent(bini)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinWidth(bini)/TMath::Sqrt(12.0)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinError(bini)<<"\n";
    }
    outputASCII_histo.close();


    //=====================================//
    // To calculate using histogram data
    //=====================================//
    const Int_t nSteps = 255;
    const Double_t min_x_pos = 0.0;// mm
    Double_t xi_l, xi_r, xi_err[nSteps] = {}, xi[nSteps] = {}, Ndch[nSteps] = {}, Ndch_err[nSteps] = {}, Nch, Nch_err, R[nSteps] = {}, R_err[nSteps] = {};
    // CH
    xi_l = fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2);
    xi_r = fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2);
    Nch = h_rp1_x_s_s->IntegralAndError(h_rp1_x_s_s->FindBin(xi_l),h_rp1_x_s_s->FindBin(xi_r),Nch_err);
    Int_t bin_start = h_rp1_x_s_s->FindBin(min_x_pos), dn_bins = 4;

    // DCH
    for(Int_t ii = 0; ii < nSteps; ii++)
    {
        xi[ii] = min_x_pos + ii*0.055;
        xi_err[ii] = 0.055/TMath::Sqrt(12.0);
        Ndch[ii] = h_rp1_x_s_s->IntegralAndError(bin_start+ii,bin_start+ii+dn_bins,Ndch_err[ii]);

        Ndch[ii] /= (dn_bins+1)*0.055;
        Ndch_err[ii] /= (dn_bins+1)*0.055;
        if(Ndch[ii] > 0)
        {
            R[ii] = Nch/Ndch[ii];
            R_err[ii] = TMath::Sqrt(TMath::Power(Nch_err/Ndch[ii],2) + TMath::Power(Nch*Ndch_err[ii]/(Ndch[ii]*Ndch[ii]),2));
        }
        else
        {
            R[ii] = 0.0;
            R_err[ii] = 0.0;
        }
    }
    //=====================================//
    // To calculate analytically
    //=====================================//
    Double_t Ndch_a[nSteps] = {}, Nch_a, R_a[nSteps] = {};
    // CH
    xi_l = fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2);
    xi_r = fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2);
    Nch_a = fit_funct_ch_dch->Integral(xi_l,xi_r);

    // DCH
    for(Int_t ii = 0; ii < nSteps; ii++)
    {
        xi_l = min_x_pos + ii*0.055;
        xi_r = min_x_pos + (ii+dn_bins)*0.055;
        Ndch_a[ii] = fit_funct_ch_dch->Integral(xi_l,xi_r);

        Ndch_a[ii] /= (dn_bins)*0.055;
        if(Ndch_a[ii] > 0)
        {
            R_a[ii] = Nch_a/Ndch_a[ii];
        }
        else
        {
            R_a[ii] = 0.0;
        }
    }
    //=====================================//
    // Ratio
    //=====================================//
    name = "c_5_";
    name += i;
    TCanvas* c_5 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_5->cd();
    gPad->SetGrid();
    TGraphErrors* gr_h = new TGraphErrors(nSteps,xi,R,xi_err,R_err);
    TGraphErrors* gr_a = new TGraphErrors(nSteps,xi,R_a,0,0);
    TMultiGraph* mg = new TMultiGraph();
    gr_h->SetLineWidth(2);
    gr_a->SetLineWidth(2);
    gr_a->SetLineColor(kBlack);
    gr_h->SetLineColor(kRed);
    gr_h->SetMarkerColor(kRed);
    mg->Add(gr_h,"AP");
    mg->Add(gr_a,"AL");
    mg->Draw("APL");
    gPad->Modified();
    mg->GetYaxis()->SetTitle("R^{CH/DCH} [per mm of septum width]");
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetXaxis()->SetTitle("X [mm]");
    mg->GetYaxis()->CenterTitle(1);
    mg->GetXaxis()->CenterTitle(1);
    mg->GetXaxis()->SetLimits(fit_dch_lim_min[i],fit_dch_lim_max[i]);
    mg->SetMaximum(mg_max[i]);
    mg->SetMinimum(mg_min[i]);
    gPad->Modified();
    name = "ratio_ch_dch_";
    name += i;
    name += ".png";
    c_5->SaveAs(name.Data());

    ofstream outputASCII_ratio;
    name = "ASCII_DATA_RATIO_METHOD3_S";
    name += i+1;
    name += ".dat";
    outputASCII_ratio.open (name.Data());
    for(Int_t pnti = 0; pnti < gr_h->GetN(); pnti++)
    {
        Double_t x, ex, y, ey;
        gr_h->GetPoint(pnti,x,y);
        ex = gr_h->GetErrorX(pnti);
        ey = gr_h->GetErrorY(pnti);

        if(x >= fit_dch_lim_min[i] && x <= fit_dch_lim_max[i])
        {
            outputASCII_ratio<<setw(10)<<x<<"\t";
            outputASCII_ratio<<setw(10)<<y<<"\t";
            outputASCII_ratio<<setw(10)<<ex<<"\t";
            outputASCII_ratio<<setw(10)<<ey<<"\n";
        }
    }
    outputASCII_ratio.close();

    return 0;
}

int function_4(Int_t i)
{
    const Int_t nSets = 13;

    Double_t mg_max[] = {/*0*/70.00,/*1*/70.00,/*2*/70.00,/*3*/70.00,/*4*/70.00,/*5*/70.00,/*6*/500.0,/*7*/500.0,/*8*/500.0,/*9*/70.00,/*10*/70.00,/*11*/70.00,/*12*/70.00};
    Double_t mg_min[] = {/*0*/20.00,/*1*/20.00,/*2*/20.00,/*3*/20.00,/*4*/20.00,/*5*/20.00,/*6*/0.000,/*7*/0.000,/*8*/0.000,/*9*/20.00,/*10*/20.00,/*11*/20.00,/*12*/20.00};

    TString fileName_RP0[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP0I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP0I_RUN_6.root"
    };

    TString fileName_RP1[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP1I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP1I_RUN_6.root"
    };

    Double_t fit_ch_lim_min[]   = {/*0*/7.00,/*1*/9.00,/*2*/9.00,/*3*/9.00,/*4*/9.00,/*5*/9.00,/*6*/9.00,/*7*/0.00,/*8*/0.00,/*9*/7.00,/*10*/9.00,/*11*/9.00,/*12*/9.00};
    Double_t fit_ch_lim_max[]   = {/*0*/9.00,/*1*/11.0,/*2*/11.0,/*3*/11.0,/*4*/11.0,/*5*/11.0,/*6*/11.0,/*7*/14.0,/*8*/14.0,/*9*/9.00,/*10*/11.0,/*11*/11.0,/*12*/11.0};
    Double_t fit_exp_lim_min[]  = {/*0*/12.0,/*1*/13.0,/*2*/13.0,/*3*/13.0,/*4*/13.0,/*5*/13.0,/*6*/13.0,/*7*/0.00,/*8*/0.00,/*9*/13.0,/*10*/13.0,/*11*/13.0,/*12*/13.0};
    Double_t fit_exp_lim_max[]  = {/*0*/13.8,/*1*/13.8,/*2*/13.8,/*3*/13.8,/*4*/13.8,/*5*/13.8,/*6*/13.8,/*7*/0.00,/*8*/0.00,/*9*/13.8,/*10*/13.8,/*11*/13.8,/*12*/13.8};
    Int_t expo_fit_bkg_status[] = {/*0*/1,   /*1*/1,   /*2*/1,   /*3*/1,   /*4*/1,   /*5*/1,   /*6*/1,   /*7*/1,   /*8*/1,   /*9*/1,   /*10*/1,   /*11*/1,   /*12*/1};       // 0 -- fit gauss+pol1 ; 1 -- fit pol1 ; -1 -- w/o fit
    Double_t fit_dch_lim_min[]  = {/*0*/2.50,/*1*/4.50,/*2*/4.50,/*3*/4.50,/*4*/4.50,/*5*/4.50,/*6*/0.50,/*7*/2.50,/*8*/2.50,/*9*/2.50,/*10*/4.50,/*11*/4.50,/*12*/4.50};
    Double_t fit_dch_lim_max[]  = {/*0*/4.00,/*1*/6.50,/*2*/6.50,/*3*/6.50,/*4*/6.50,/*5*/6.50,/*6*/6.50,/*7*/4.00,/*8*/4.00,/*9*/4.00,/*10*/6.50,/*11*/6.50,/*12*/6.50};
    Float_t factor_bad_pixels[] = {/*0*/200.,/*1*/200.,/*2*/200.,/*3*/200.,/*4*/200.,/*5*/30.0,/*6*/4.00,/*7*/100.,/*8*/30.0,/*9*/200.,/*10*/200.,/*11*/100.,/*12*/50.0};
    Double_t par_ch_bkg[nSets][5];

    // Gauss
    par_ch_bkg[0][0] =  0.60; par_ch_bkg[1][0] =  0.60; par_ch_bkg[2][0] =  0.06; par_ch_bkg[3][0] =  0.10;
    par_ch_bkg[4][0] =  0.60; par_ch_bkg[5][0] =  0.60; par_ch_bkg[6][0] =  0.60; par_ch_bkg[7][0] =  0.60;
    par_ch_bkg[8][0] =  0.60; par_ch_bkg[9][0] =  0.60; par_ch_bkg[10][0] = 0.60; par_ch_bkg[11][0] = 0.40;
    par_ch_bkg[12][0] = 0.60;

    par_ch_bkg[0][1] =  7.50; par_ch_bkg[1][1] =  7.50; par_ch_bkg[2][1] =  10.0; par_ch_bkg[3][1] =  7.50;
    par_ch_bkg[4][1] =  7.50; par_ch_bkg[5][1] =  7.50; par_ch_bkg[6][1] =  7.50; par_ch_bkg[7][1] =  10.00;
    par_ch_bkg[8][1] =  7.50; par_ch_bkg[9][1] =  7.50; par_ch_bkg[10][1] = 7.50; par_ch_bkg[11][1] = 10.00;
    par_ch_bkg[12][1] = 7.50;

    par_ch_bkg[0][2] =  0.80; par_ch_bkg[1][2] =  0.80; par_ch_bkg[2][2] =  0.70; par_ch_bkg[3][2] =  0.70;
    par_ch_bkg[4][2] =  0.80; par_ch_bkg[5][2] =  0.80; par_ch_bkg[6][2] =  0.80; par_ch_bkg[7][2] =  0.80;
    par_ch_bkg[8][2] =  0.80; par_ch_bkg[9][2] =  0.80; par_ch_bkg[10][2] = 0.80; par_ch_bkg[11][2] = 1.40;
    par_ch_bkg[12][2] = 0.80;

    // Pol1
    par_ch_bkg[0][3] = -4.40; par_ch_bkg[1][3] = -3.70; par_ch_bkg[2][3] = -4.16; par_ch_bkg[3][3] = -3.80;
    par_ch_bkg[4][3] = -2.50; par_ch_bkg[5][3] = -6.00; par_ch_bkg[6][3] = -2.50; par_ch_bkg[7][3] = -2.50;
    par_ch_bkg[8][3] = -2.50; par_ch_bkg[9][3] = -2.50; par_ch_bkg[10][3] =-2.50; par_ch_bkg[11][3] =-1.60;
    par_ch_bkg[12][3] =-2.50;

    par_ch_bkg[0][4] = -0.20; par_ch_bkg[1][4] = -0.20; par_ch_bkg[2][4] = -0.20; par_ch_bkg[3][4] = -0.17;
    par_ch_bkg[4][4] = -0.06; par_ch_bkg[5][4] = -0.02; par_ch_bkg[6][4] = -0.06; par_ch_bkg[7][4] = -0.06;
    par_ch_bkg[8][4] = -0.06; par_ch_bkg[9][4] = -0.06; par_ch_bkg[10][4] =-0.06; par_ch_bkg[11][4] =-0.01;
    par_ch_bkg[12][4] =-0.06;


    TH2D* h_rp0[nSets];
    TH2D* h_rp1[nSets];

    gStyle->SetOptStat(0);

    cout<<endl<<"--> Set "<<i+1<<" <--"<<endl;
    TFile *_file_rp0 = TFile::Open(fileName_RP0[i].Data());
    TFile *_file_rp1 = TFile::Open(fileName_RP1[i].Data());

    h_rp0[i] = (TH2D*)_file_rp0->Get("h_6");
    h_rp1[i] = (TH2D*)_file_rp1->Get("h_6");

    removenoizypixelsXY(h_rp0[i],factor_bad_pixels[i]);
    removenoizypixelsXY(h_rp1[i],factor_bad_pixels[i]);

    Double_t integral_err_rp0;
    Double_t integral_rp0 = h_rp0[i]->IntegralAndError(1,h_rp0[i]->GetNbinsX(),1,h_rp0[i]->GetNbinsY(),integral_err_rp0);
    scale2Dhisto(h_rp0[i], integral_rp0, integral_err_rp0);
    scale2Dhisto(h_rp1[i], integral_rp0, integral_err_rp0);
    cout<<"--> Integral RP1: "<<h_rp1[i]->Integral()<<endl;

    h_rp1[i]->GetXaxis()->SetRange(1,h_rp1[i]->GetNbinsX()/2);
    h_rp1[i]->GetYaxis()->SetRange(1,h_rp1[i]->GetNbinsY()/2);
    h_rp1[i]->GetXaxis()-> CenterTitle();
    h_rp1[i]->GetYaxis()-> CenterTitle();

    TString name = "c_1_";
    name += i;
    TCanvas* c_1 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_1->cd();
    h_rp1[i]->Draw("colz");

    name = "h_rp1_x_";
    name += i+1;
    TH1D* h_rp1_x = h_rp1[i]->ProjectionX(name.Data());
    name += "_s";
    TH1D* h_rp1_x_s = h_rp1[i]->ProjectionX(name.Data());
    name += "_s";
    TH1D* h_rp1_x_s_s = h_rp1[i]->ProjectionX(name.Data());

    name = "c_2_";
    name += i;
    TCanvas* c_2 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_2->cd();
    h_rp1_x->SetMinimum(0);
    h_rp1_x->Draw("hist");
    TF1* fit_funct_ch_bkg;

    if(expo_fit_bkg_status[i] == 1)
    {
        fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","pol1",fit_exp_lim_min[i],fit_exp_lim_max[i]);
        fit_funct_ch_bkg->SetParameter(0,par_ch_bkg[i][3]);
        fit_funct_ch_bkg->SetParameter(1,par_ch_bkg[i][4]);
        h_rp1_x->Fit(fit_funct_ch_bkg,"R+");
        par_ch_bkg[i][3] = fit_funct_ch_bkg->GetParameter(0);
        par_ch_bkg[i][4] = fit_funct_ch_bkg->GetParameter(1);
    }
    else
    {
        fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","gaus(0)+pol1(3)",fit_ch_lim_min[i],fit_exp_lim_max[i]);
        fit_funct_ch_bkg->SetParameters(par_ch_bkg[i]);
        h_rp1_x->Fit(fit_funct_ch_bkg,"R+");
        fit_funct_ch_bkg->GetParameters(par_ch_bkg[i]);
    }

    name = "#chi^{2}/NDF = [";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare();
    name += "/";
    name += fit_funct_ch_bkg->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare()/fit_funct_ch_bkg->GetNDF();
    h_rp1_x->SetTitle(name.Data());
    fit_funct_ch_bkg->SetLineWidth(6);
    if(expo_fit_bkg_status[i] >= 0) fit_funct_ch_bkg->Draw("same");
    TF1* bkg_func = new TF1("bkg_func","pol1",h_rp1_x->GetBinCenter(1),h_rp1_x->GetBinCenter(h_rp1_x->GetNbinsX()));
    bkg_func->SetParameter(0,par_ch_bkg[i][3]);
    bkg_func->SetParameter(1,par_ch_bkg[i][4]);
    bkg_func->SetLineColor(kGreen+2);
    if(expo_fit_bkg_status[i] >= 0) bkg_func->Draw("same");

    name = "c_3_";
    name += i;
    TCanvas* c_3 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_3->cd();
    h_rp1_x_s->SetMinimum(0);
    h_rp1_x_s->SetLineColor(kBlue);
    h_rp1_x_s->SetLineWidth(2);
    h_rp1_x_s->Draw("hist");
    for(Int_t j = 1; j <= h_rp1_x_s_s->GetNbinsX(); j++)
    {
        Double_t val;
        if(expo_fit_bkg_status[i] >= 0)
            val = h_rp1_x_s->GetBinContent(j) - bkg_func->Eval(h_rp1_x_s->GetBinCenter(j));
        else
            val = h_rp1_x_s->GetBinContent(j);



        if(val > 0)
            h_rp1_x_s_s->SetBinContent(j,val);
        else
            h_rp1_x_s_s->SetBinContent(j,0);
    }
    h_rp1_x_s_s->SetLineColor(kRed);
    h_rp1_x_s_s->SetLineWidth(2);
    h_rp1_x_s_s->Draw("same & hist");

    name = "c_4_";
    name += i;
    TCanvas* c_4 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_4->cd();
    h_rp1_x_s_s->SetLineWidth(8);
    h_rp1_x_s_s->Draw("hist");

    TF1* fit_funct_1 = new TF1("fit_funct_1","gaus",fit_ch_lim_min[i],fit_ch_lim_max[i]);
    h_rp1_x_s_s->Fit(fit_funct_1,"R0Q+");
    Double_t par_ch[3];
    fit_funct_1->GetParameters(par_ch);
    TF1* fit_funct_ch = new TF1("fit_funct_ch","gaus",h_rp1_x_s_s->GetBinCenter(1),h_rp1_x_s_s->GetBinCenter(h_rp1_x_s_s->GetNbinsX()));
    fit_funct_ch->SetParameters(par_ch);
    fit_funct_ch->SetLineColor(kBlue);
    fit_funct_ch->SetLineWidth(2);
//    fit_funct_ch->Draw("same");

    TF1* fit_funct_2 = new TF1("fit_funct_2","expo",fit_dch_lim_min[i],fit_dch_lim_max[i]);
    h_rp1_x_s_s->Fit(fit_funct_2,"R0Q+");
    Double_t par_dch[2];
    fit_funct_2->GetParameters(par_dch);
    TF1* fit_funct_dch = new TF1("fit_funct_dch","expo",h_rp1_x_s_s->GetBinCenter(1),h_rp1_x_s_s->GetBinCenter(h_rp1_x_s_s->GetNbinsX()));
    fit_funct_dch->SetParameters(par_dch);
    fit_funct_dch->SetLineColor(kBlack);
    fit_funct_dch->SetLineWidth(2);
//    fit_funct_dch->Draw("same");

    TLine *l_min_ch=new TLine(fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_min_ch->SetLineColor(kMagenta);
    l_min_ch->SetLineWidth(3);
    l_min_ch->SetLineStyle(7);
    l_min_ch->Draw();
    TLine *l_mean_ch=new TLine(fit_funct_ch->GetParameter(1),0,fit_funct_ch->GetParameter(1),h_rp1_x_s_s->GetMaximum());
    l_mean_ch->SetLineColor(kMagenta-2);
    l_mean_ch->SetLineWidth(3);
    l_mean_ch->SetLineStyle(7);
    l_mean_ch->Draw();
    TLine *l_max_ch=new TLine(fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_max_ch->SetLineColor(kMagenta);
    l_max_ch->SetLineWidth(3);
    l_max_ch->SetLineStyle(7);
    l_max_ch->Draw();

    TLine *l_min_dch=new TLine(fit_dch_lim_min[i],0,fit_dch_lim_min[i],h_rp1_x_s_s->GetMaximum());
    l_min_dch->SetLineColor(kGreen+1);
    l_min_dch->SetLineWidth(3);
    l_min_dch->SetLineStyle(7);
    l_min_dch->Draw();
    TLine *l_max_dch=new TLine(fit_dch_lim_max[i],0,fit_dch_lim_max[i],h_rp1_x_s_s->GetMaximum());
    l_max_dch->SetLineColor(kGreen+1);
    l_max_dch->SetLineWidth(3);
    l_max_dch->SetLineStyle(7);
    l_max_dch->Draw();

    TF1* fit_funct_ch_dch = new TF1("fit_funct_ch_dch",fit_ch_dch,fit_dch_lim_min[i],fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),5);
    fit_funct_ch_dch->SetParameter(0,par_ch[0]);
    fit_funct_ch_dch->SetParameter(1,par_ch[1]);
    fit_funct_ch_dch->SetParameter(2,par_ch[2]);
    fit_funct_ch_dch->SetParameter(3,par_dch[0]);
    fit_funct_ch_dch->SetParameter(4,par_dch[1]);
    fit_funct_ch_dch->SetLineColor(kBlack);
    fit_funct_ch_dch->SetLineWidth(4);
    fit_funct_ch_dch->Draw("same");

    name = "#chi^{2}/NDF (DCH) = [";
    name += (Int_t)fit_funct_2->GetChisquare();
    name += "/";
    name += fit_funct_2->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_2->GetChisquare()/fit_funct_2->GetNDF();
    name += " | #chi^{2}/NDF (CH) = [";
    name += (Int_t)fit_funct_1->GetChisquare();
    name += "/";
    name += fit_funct_1->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_1->GetChisquare()/fit_funct_1->GetNDF();
    h_rp1_x_s_s->SetTitle(name.Data());

    gPad->Modified();
    name = "profile_ch_dch_";
    name += i;
    name += ".png";
    c_4->SaveAs(name.Data());

    ofstream outputASCII_histo;
    name = "ASCII_DATA_PROFILE_METHOD4_S";
    name += i+1;
    name += ".dat";
    outputASCII_histo.open (name.Data());
    for(Int_t bini = 1; bini <= h_rp1_x_s_s->GetNbinsX(); bini++)
    {
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinCenter(bini)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinContent(bini)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinWidth(bini)/TMath::Sqrt(12.0)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s_s->GetBinError(bini)<<"\n";
    }
    outputASCII_histo.close();

    //=====================================//
    // To calculate using histogram data
    //=====================================//
    const Int_t nSteps = 255;
    const Double_t min_x_pos = 0.0;// mm
    Double_t xi_l, xi_r, xi_err[nSteps] = {}, xi[nSteps] = {}, Ndch[nSteps] = {}, Ndch_err[nSteps] = {}, Nch, Nch_err, R[nSteps] = {}, R_err[nSteps] = {};
    // CH
    xi_l = fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2);
    xi_r = fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2);
    Nch = h_rp1_x_s_s->IntegralAndError(h_rp1_x_s_s->FindBin(xi_l),h_rp1_x_s_s->FindBin(xi_r),Nch_err);
    Int_t bin_start = h_rp1_x_s_s->FindBin(min_x_pos), dn_bins = 4;

    // DCH
    for(Int_t ii = 0; ii < nSteps; ii++)
    {
        xi[ii] = min_x_pos + ii*0.055;
        xi_err[ii] = 0.055/TMath::Sqrt(12.0);
        Ndch[ii] = h_rp1_x_s_s->IntegralAndError(bin_start+ii,bin_start+ii+dn_bins,Ndch_err[ii]);

        Ndch[ii] /= (dn_bins+1)*0.055;
        Ndch_err[ii] /= (dn_bins+1)*0.055;
        if(Ndch[ii] > 0)
        {
            R[ii] = Nch/Ndch[ii];
            R_err[ii] = TMath::Sqrt(TMath::Power(Nch_err/Ndch[ii],2) + TMath::Power(Nch*Ndch_err[ii]/(Ndch[ii]*Ndch[ii]),2));
        }
        else
        {
            R[ii] = 0.0;
            R_err[ii] = 0.0;
        }
    }
    //=====================================//
    // To calculate analytically
    //=====================================//
    Double_t Ndch_a[nSteps] = {}, Nch_a, R_a[nSteps] = {};
    // CH
    xi_l = fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2);
    xi_r = fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2);
    Nch_a = fit_funct_ch_dch->Integral(xi_l,xi_r);

    // DCH
    for(Int_t ii = 0; ii < nSteps; ii++)
    {
        xi_l = min_x_pos + ii*0.055;
        xi_r = min_x_pos + (ii+dn_bins)*0.055;
        Ndch_a[ii] = fit_funct_ch_dch->Integral(xi_l,xi_r);

        Ndch_a[ii] /= (dn_bins)*0.055;
        if(Ndch_a[ii] > 0)
        {
            R_a[ii] = Nch_a/Ndch_a[ii];
        }
        else
        {
            R_a[ii] = 0.0;
        }
    }
    //=====================================//
    // Ratio
    //=====================================//
    name = "c_5_";
    name += i;
    TCanvas* c_5 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_5->cd();
    gPad->SetGrid();
    TGraphErrors* gr_h = new TGraphErrors(nSteps,xi,R,xi_err,R_err);
    TGraphErrors* gr_a = new TGraphErrors(nSteps,xi,R_a,0,0);
    TMultiGraph* mg = new TMultiGraph();
    gr_h->SetLineWidth(2);
    gr_a->SetLineWidth(2);
    gr_a->SetLineColor(kBlack);
    gr_h->SetLineColor(kRed);
    gr_h->SetMarkerColor(kRed);
    mg->Add(gr_h,"AP");
    mg->Add(gr_a,"AL");
    mg->Draw("APL");
    gPad->Modified();
    mg->GetYaxis()->SetTitle("R^{CH/DCH} [per mm of septum width]");
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetXaxis()->SetTitle("X [mm]");
    mg->GetYaxis()->CenterTitle(1);
    mg->GetXaxis()->CenterTitle(1);
    mg->GetXaxis()->SetLimits(fit_dch_lim_min[i],fit_dch_lim_max[i]);
    mg->SetMaximum(mg_max[i]);
    mg->SetMinimum(mg_min[i]);
    gPad->Modified();
    name = "ratio_ch_dch_";
    name += i;
    name += ".png";
    c_5->SaveAs(name.Data());

    ofstream outputASCII_ratio;
    name = "ASCII_DATA_RATIO_METHOD4_S";
    name += i+1;
    name += ".dat";
    outputASCII_ratio.open (name.Data());
    for(Int_t pnti = 0; pnti < gr_h->GetN(); pnti++)
    {
        Double_t x, ex, y, ey;
        gr_h->GetPoint(pnti,x,y);
        ex = gr_h->GetErrorX(pnti);
        ey = gr_h->GetErrorY(pnti);

        if(x >= fit_dch_lim_min[i] && x <= fit_dch_lim_max[i])
        {
            outputASCII_ratio<<setw(10)<<x<<"\t";
            outputASCII_ratio<<setw(10)<<y<<"\t";
            outputASCII_ratio<<setw(10)<<ex<<"\t";
            outputASCII_ratio<<setw(10)<<ey<<"\n";
        }
    }
    outputASCII_ratio.close();

    return 0;
}

int function_5(Int_t i)
{
    const Int_t nSets = 13;

    Double_t mg_max[] = {/*0*/40.00,/*1*/40.00,/*2*/40.00,/*3*/40.00,/*4*/40.00,/*5*/40.00,/*6*/500.0,/*7*/500.0,/*8*/500.0,/*9*/40.00,/*10*/40.00,/*11*/40.00,/*12*/40.00};
    Double_t mg_min[] = {/*0*/5.000,/*1*/5.000,/*2*/5.000,/*3*/5.000,/*4*/5.000,/*5*/5.000,/*6*/0.000,/*7*/0.000,/*8*/0.000,/*9*/5.000,/*10*/5.000,/*11*/5.000,/*12*/5.000};

    TString fileName_RP0[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP0I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP0I_RUN_6.root"
    };

    TString fileName_RP1[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP1I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP1I_RUN_6.root"
    };

    Double_t fit_ch_lim_min[]   = {/*0*/7.00,/*1*/9.00,/*2*/9.00,/*3*/9.00,/*4*/9.00,/*5*/9.00,/*6*/9.00,/*7*/0.00,/*8*/0.00,/*9*/7.00,/*10*/9.00,/*11*/9.00,/*12*/9.00};
    Double_t fit_ch_lim_max[]   = {/*0*/9.00,/*1*/11.0,/*2*/11.0,/*3*/11.0,/*4*/11.0,/*5*/11.0,/*6*/11.0,/*7*/14.0,/*8*/14.0,/*9*/9.00,/*10*/11.0,/*11*/11.0,/*12*/11.0};
    Double_t fit_dch_lim_min[]  = {/*0*/2.50,/*1*/4.50,/*2*/4.50,/*3*/4.50,/*4*/4.50,/*5*/4.50,/*6*/0.50,/*7*/2.50,/*8*/2.50,/*9*/2.50,/*10*/4.50,/*11*/4.50,/*12*/4.50};
    Double_t fit_dch_lim_max[]  = {/*0*/4.00,/*1*/6.50,/*2*/6.50,/*3*/6.50,/*4*/6.50,/*5*/6.50,/*6*/6.50,/*7*/4.00,/*8*/4.00,/*9*/4.00,/*10*/6.50,/*11*/6.50,/*12*/6.50};
    Float_t factor_bad_pixels[] = {/*0*/200.,/*1*/200.,/*2*/200.,/*3*/200.,/*4*/200.,/*5*/30.0,/*6*/4.00,/*7*/100.,/*8*/30.0,/*9*/200.,/*10*/200.,/*11*/100.,/*12*/50.0};
    Double_t par_ch_bkg[nSets][5];

    // Gauss
    par_ch_bkg[0][0] =  0.60; par_ch_bkg[1][0] =  0.60; par_ch_bkg[2][0] =  0.06; par_ch_bkg[3][0] =  0.10;
    par_ch_bkg[4][0] =  0.60; par_ch_bkg[5][0] =  0.60; par_ch_bkg[6][0] =  0.60; par_ch_bkg[7][0] =  0.60;
    par_ch_bkg[8][0] =  0.60; par_ch_bkg[9][0] =  0.60; par_ch_bkg[10][0] = 0.60; par_ch_bkg[11][0] = 0.40;
    par_ch_bkg[12][0] = 0.60;

    par_ch_bkg[0][1] =  7.50; par_ch_bkg[1][1] =  7.50; par_ch_bkg[2][1] =  10.0; par_ch_bkg[3][1] =  7.50;
    par_ch_bkg[4][1] =  7.50; par_ch_bkg[5][1] =  7.50; par_ch_bkg[6][1] =  7.50; par_ch_bkg[7][1] =  10.00;
    par_ch_bkg[8][1] =  7.50; par_ch_bkg[9][1] =  7.50; par_ch_bkg[10][1] = 7.50; par_ch_bkg[11][1] = 10.00;
    par_ch_bkg[12][1] = 7.50;

    par_ch_bkg[0][2] =  0.80; par_ch_bkg[1][2] =  0.80; par_ch_bkg[2][2] =  0.70; par_ch_bkg[3][2] =  0.70;
    par_ch_bkg[4][2] =  0.80; par_ch_bkg[5][2] =  0.80; par_ch_bkg[6][2] =  0.80; par_ch_bkg[7][2] =  0.80;
    par_ch_bkg[8][2] =  0.80; par_ch_bkg[9][2] =  0.80; par_ch_bkg[10][2] = 0.80; par_ch_bkg[11][2] = 1.40;
    par_ch_bkg[12][2] = 0.80;

    // Pol1
    par_ch_bkg[0][3] = -4.40; par_ch_bkg[1][3] = -3.70; par_ch_bkg[2][3] = -4.16; par_ch_bkg[3][3] = -3.80;
    par_ch_bkg[4][3] = -2.50; par_ch_bkg[5][3] = -6.00; par_ch_bkg[6][3] = -2.50; par_ch_bkg[7][3] = -2.50;
    par_ch_bkg[8][3] = -2.50; par_ch_bkg[9][3] = -2.50; par_ch_bkg[10][3] =-2.50; par_ch_bkg[11][3] =-1.60;
    par_ch_bkg[12][3] =-2.50;

    par_ch_bkg[0][4] = -0.20; par_ch_bkg[1][4] = -0.20; par_ch_bkg[2][4] = -0.20; par_ch_bkg[3][4] = -0.17;
    par_ch_bkg[4][4] = -0.06; par_ch_bkg[5][4] = -0.02; par_ch_bkg[6][4] = -0.06; par_ch_bkg[7][4] = -0.06;
    par_ch_bkg[8][4] = -0.06; par_ch_bkg[9][4] = -0.06; par_ch_bkg[10][4] =-0.06; par_ch_bkg[11][4] =-0.01;
    par_ch_bkg[12][4] =-0.06;


    TH2D* h_rp0[nSets];
    TH2D* h_rp1[nSets];

    gStyle->SetOptStat(0);

    cout<<endl<<"--> Set "<<i+1<<" <--"<<endl;
    TFile *_file_rp0 = TFile::Open(fileName_RP0[i].Data());
    TFile *_file_rp1 = TFile::Open(fileName_RP1[i].Data());

    h_rp0[i] = (TH2D*)_file_rp0->Get("h_6");
    h_rp1[i] = (TH2D*)_file_rp1->Get("h_6");

    removenoizypixelsXY(h_rp0[i],factor_bad_pixels[i]);
    removenoizypixelsXY(h_rp1[i],factor_bad_pixels[i]);

    Double_t integral_err_rp0;
    Double_t integral_rp0 = h_rp0[i]->IntegralAndError(1,h_rp0[i]->GetNbinsX(),1,h_rp0[i]->GetNbinsY(),integral_err_rp0);
    scale2Dhisto(h_rp0[i], integral_rp0, integral_err_rp0);
    scale2Dhisto(h_rp1[i], integral_rp0, integral_err_rp0);
    cout<<"--> Integral RP1: "<<h_rp1[i]->Integral()<<endl;

    h_rp1[i]->GetXaxis()->SetRange(1,h_rp1[i]->GetNbinsX()/2);
    h_rp1[i]->GetYaxis()->SetRange(1,h_rp1[i]->GetNbinsY()/2);
    h_rp1[i]->GetXaxis()-> CenterTitle();
    h_rp1[i]->GetYaxis()-> CenterTitle();

    TString name = "c_1_";
    name += i;
    TCanvas* c_1 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_1->cd();
    h_rp1[i]->Draw("colz");

    name = "h_rp1_x_";
    name += i+1;
    name += "_s";
    TH1D* h_rp1_x_s = h_rp1[i]->ProjectionX(name.Data());

    name = "c_3_";
    name += i;
    TCanvas* c_3 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_3->cd();
    h_rp1_x_s->SetMinimum(0);
    h_rp1_x_s->SetLineColor(kRed);
    h_rp1_x_s->SetLineWidth(2);
    h_rp1_x_s->Draw("same & hist");

    name = "c_4_";
    name += i;
    TCanvas* c_4 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_4->cd();
    h_rp1_x_s->SetLineWidth(8);
    h_rp1_x_s->SetTitle("");
    h_rp1_x_s->Draw("hist");

    TF1* fit_funct_1 = new TF1("fit_funct_1","gaus",fit_ch_lim_min[i],fit_ch_lim_max[i]);
    h_rp1_x_s->Fit(fit_funct_1,"R+");
    Double_t par_ch[3];
    fit_funct_1->GetParameters(par_ch);
    TF1* fit_funct_ch = new TF1("fit_funct_ch","gaus",h_rp1_x_s->GetBinCenter(1),h_rp1_x_s->GetBinCenter(h_rp1_x_s->GetNbinsX()));
    fit_funct_ch->SetParameters(par_ch);
    fit_funct_ch->SetLineColor(kBlue);
    fit_funct_ch->SetLineWidth(2);
    fit_funct_ch->Draw("same");

    if(fit_dch_lim_max[i] > fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2))
        fit_dch_lim_max[i] = fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2) - 0.055;
    if(fit_dch_lim_min[i] > fit_dch_lim_max[i])
        fit_dch_lim_min[i] = fit_dch_lim_max[i];

    TF1* fit_funct_2 = new TF1("fit_funct_2","expo",fit_dch_lim_min[i],fit_dch_lim_max[i]);
    h_rp1_x_s->Fit(fit_funct_2,"R0Q+");
    Double_t par_dch[2];
    fit_funct_2->GetParameters(par_dch);
    TF1* fit_funct_dch = new TF1("fit_funct_dch","expo",h_rp1_x_s->GetBinCenter(1),h_rp1_x_s->GetBinCenter(h_rp1_x_s->GetNbinsX()));
    fit_funct_dch->SetParameters(par_dch);
    fit_funct_dch->SetLineColor(kBlack);
    fit_funct_dch->SetLineWidth(2);
    fit_funct_dch->Draw("same");

    TLine *l_min_ch=new TLine(fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s->GetMaximum());
    l_min_ch->SetLineColor(kMagenta);
    l_min_ch->SetLineWidth(3);
    l_min_ch->SetLineStyle(7);
    l_min_ch->Draw();
    TLine *l_mean_ch=new TLine(fit_funct_ch->GetParameter(1),0,fit_funct_ch->GetParameter(1),h_rp1_x_s->GetMaximum());
    l_mean_ch->SetLineColor(kMagenta-2);
    l_mean_ch->SetLineWidth(3);
    l_mean_ch->SetLineStyle(7);
    l_mean_ch->Draw();
    TLine *l_max_ch=new TLine(fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s->GetMaximum());
    l_max_ch->SetLineColor(kMagenta);
    l_max_ch->SetLineWidth(3);
    l_max_ch->SetLineStyle(7);
    l_max_ch->Draw();

    TLine *l_min_dch=new TLine(fit_dch_lim_min[i],0,fit_dch_lim_min[i],h_rp1_x_s->GetMaximum());
    l_min_dch->SetLineColor(kGreen+1);
    l_min_dch->SetLineWidth(3);
    l_min_dch->SetLineStyle(7);
    l_min_dch->Draw();
    TLine *l_max_dch=new TLine(fit_dch_lim_max[i],0,fit_dch_lim_max[i],h_rp1_x_s->GetMaximum());
    l_max_dch->SetLineColor(kGreen+1);
    l_max_dch->SetLineWidth(3);
    l_max_dch->SetLineStyle(7);
    l_max_dch->Draw();

    name += ".png";
    c_4->SaveAs(name.Data());

    ofstream outputASCII_histo;
    name = "ASCII_DATA_PROFILE_METHOD5_S";
    name += i+1;
    name += ".dat";
    outputASCII_histo.open (name.Data());
    for(Int_t bini = 1; bini <= h_rp1_x_s->GetNbinsX(); bini++)
    {
        outputASCII_histo<<setw(10)<<h_rp1_x_s->GetBinCenter(bini)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s->GetBinContent(bini)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s->GetBinWidth(bini)/TMath::Sqrt(12.0)<<"\t";
        outputASCII_histo<<setw(10)<<h_rp1_x_s->GetBinError(bini)<<"\n";
    }
    outputASCII_histo.close();

    //=====================================//
    // To calculate using histogram data
    //=====================================//
    const Int_t nSteps = 255;
    const Double_t min_x_pos = 0.0;// mm
    Double_t xi_l, xi_r, xi_err[nSteps] = {}, xi[nSteps] = {}, Ndch[nSteps] = {}, Ndch_err[nSteps] = {}, Nch, Nch_err, R[nSteps] = {}, R_err[nSteps] = {};
    // CH
    xi_l = fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2);
    xi_r = fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2);
    Nch = h_rp1_x_s->IntegralAndError(h_rp1_x_s->FindBin(xi_l),h_rp1_x_s->FindBin(xi_r),Nch_err);
    Int_t bin_start = h_rp1_x_s->FindBin(min_x_pos), dn_bins = 4;

    // DCH
    for(Int_t ii = 0; ii < nSteps; ii++)
    {
        xi[ii] = min_x_pos + ii*0.055;
        xi_err[ii] = 0.055/TMath::Sqrt(12.0);
        Ndch[ii] = h_rp1_x_s->IntegralAndError(bin_start+ii,bin_start+ii+dn_bins,Ndch_err[ii]);

        Ndch[ii] /= (dn_bins+1)*0.055;
        Ndch_err[ii] /= (dn_bins+1)*0.055;
        if(Ndch[ii] > 0)
        {
            R[ii] = Nch/Ndch[ii];
            R_err[ii] = TMath::Sqrt(TMath::Power(Nch_err/Ndch[ii],2) + TMath::Power(Nch*Ndch_err[ii]/(Ndch[ii]*Ndch[ii]),2));
        }
        else
        {
            R[ii] = 0.0;
            R_err[ii] = 0.0;
        }
    }
    //=====================================//
    // Ratio
    //=====================================//
    name = "c_5_";
    name += i;
    TCanvas* c_5 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_5->cd();
    gPad->SetGrid();
    TGraphErrors* gr_h = new TGraphErrors(nSteps,xi,R,xi_err,R_err);
    TMultiGraph* mg = new TMultiGraph();
    gr_h->SetLineWidth(2);
    gr_h->SetLineColor(kRed);
    gr_h->SetMarkerColor(kRed);
    mg->Add(gr_h,"AP");
    mg->Draw("APL");
    gPad->Modified();
    mg->GetYaxis()->SetTitle("R^{CH/DCH} [per mm of septum width]");
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetXaxis()->SetTitle("X [mm]");
    mg->GetYaxis()->CenterTitle(1);
    mg->GetXaxis()->CenterTitle(1);
    mg->GetXaxis()->SetLimits(fit_dch_lim_min[i],fit_dch_lim_max[i]);
    mg->SetMaximum(mg_max[i]);
    mg->SetMinimum(mg_min[i]);
    gPad->Modified();
    name = "ratio_ch_dch_";
    name += i;
    name += ".png";
    c_5->SaveAs(name.Data());

    ofstream outputASCII_ratio;
    name = "ASCII_DATA_RATIO_METHOD5_S";
    name += i+1;
    name += ".dat";
    outputASCII_ratio.open (name.Data());
    for(Int_t pnti = 0; pnti < gr_h->GetN(); pnti++)
    {
        Double_t x, ex, y, ey;
        gr_h->GetPoint(pnti,x,y);
        ex = gr_h->GetErrorX(pnti);
        ey = gr_h->GetErrorY(pnti);

        if(x >= fit_dch_lim_min[i] && x <= fit_dch_lim_max[i])
        {
            outputASCII_ratio<<setw(10)<<x<<"\t";
            outputASCII_ratio<<setw(10)<<y<<"\t";
            outputASCII_ratio<<setw(10)<<ex<<"\t";
            outputASCII_ratio<<setw(10)<<ey<<"\n";
        }
    }
    outputASCII_ratio.close();

    return 0;
}

int function_6()
{
//    TString fileName_RP1 = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_06_18_HISTO_RP0I_M1_RUN_7.root";
    TString fileName_RP1 = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S17_HISTO_RP1I_RUN_8.root";
//    TString fileName_RP1 = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_M3_RUN_1.root";

    Float_t factor_bad_pixels = 500;

    gStyle->SetOptStat(0);

    TFile *_file_rp1 = TFile::Open(fileName_RP1.Data());

    TH2D* h_rp0 = (TH2D*)_file_rp1->Get("h_6");
    TH2D* h_rp1 = (TH2D*)_file_rp1->Get("h_6");

    removenoizypixelsXY(h_rp0,factor_bad_pixels);
    removenoizypixelsXY(h_rp1,factor_bad_pixels);

    Double_t integral_err_rp0;
    Double_t integral_rp0 = h_rp0->IntegralAndError(1,h_rp0->GetNbinsX(),1,h_rp0->GetNbinsY(),integral_err_rp0);
    scale2Dhisto(h_rp0, integral_rp0, integral_err_rp0);
    scale2Dhisto(h_rp1, integral_rp0, integral_err_rp0);

    h_rp1->GetXaxis()->SetRange(1,h_rp1->GetNbinsX()/2);
    h_rp1->GetYaxis()->SetRange(1,h_rp1->GetNbinsY()/2);
    h_rp1->GetXaxis()-> CenterTitle();
    h_rp1->GetYaxis()-> CenterTitle();

    TCanvas* c_1 = new TCanvas("c_1","c_1",1000,1000);
    c_1->cd();
    h_rp1->Draw("colz");

    TH1D* h_rp1_x_1 = h_rp1->ProjectionX("h_rp1_x_1");
    TH1D* h_rp1_x_2 = h_rp1->ProjectionX("h_rp1_x_2");
    TH1D* h_rp1_x_3 = h_rp1->ProjectionX("h_rp1_x_3");

    TCanvas* c_2 = new TCanvas("c_2","c_2",1000,500);
    c_2->cd();
    h_rp1_x_1->SetMinimum(0);
    h_rp1_x_1->SetLineColor(kBlue);
    h_rp1_x_1->Draw();


    return 0;
}

int function_7()
{
    TString fileName_RP1 = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP1I_RUN_6.root";
    TString fileName_RP0 = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP0I_RUN_6.root";

    TString fileName_RP1_AM = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP1I_RUN_6.root";
    TString fileName_RP0_AM = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP0I_RUN_6.root";

    Float_t factor_bad_pixels = 200;
    Float_t factor_bad_pixels_am = 30;

    gStyle->SetOptStat(0);

    TFile *_file_rp1 = TFile::Open(fileName_RP1.Data());
    TFile *_file_rp0 = TFile::Open(fileName_RP0.Data());

    TFile *_file_rp1_am = TFile::Open(fileName_RP1_AM.Data());
    TFile *_file_rp0_am = TFile::Open(fileName_RP0_AM.Data());

    TH2D* h_rp1 = (TH2D*)_file_rp1->Get("h_6");
    TH2D* h_rp0 = (TH2D*)_file_rp0->Get("h_6");

    TH2D* h_rp1_am = (TH2D*)_file_rp1_am->Get("h_6");
    TH2D* h_rp0_am = (TH2D*)_file_rp0_am->Get("h_6");

    removenoizypixelsXY(h_rp1,factor_bad_pixels);
    removenoizypixelsXY(h_rp0,factor_bad_pixels);

    removenoizypixelsXY(h_rp1_am,factor_bad_pixels_am);
    removenoizypixelsXY(h_rp0_am,factor_bad_pixels_am);

    Double_t integral_err_rp0;
    Double_t integral_err_rp0_am;
    Double_t integral_rp0 = h_rp0->IntegralAndError(1,h_rp0->GetNbinsX(),1,h_rp0->GetNbinsY(),integral_err_rp0);
    Double_t integral_rp0_am = h_rp0_am->IntegralAndError(1,h_rp0_am->GetNbinsX(),1,h_rp0_am->GetNbinsY(),integral_err_rp0_am);
//    Double_t integral_rp0 = h_rp1->IntegralAndError(1,h_rp0->GetNbinsX(),1,h_rp0->GetNbinsY(),integral_err_rp0);
//    Double_t integral_rp0_am = h_rp1_am->IntegralAndError(1,h_rp0_am->GetNbinsX(),1,h_rp0_am->GetNbinsY(),integral_err_rp0_am);

    cout<<"--> Ratio RP1/RP0 CRY2 in CH: "<<h_rp1->Integral(1,h_rp1->GetNbinsX(),1,h_rp1->GetNbinsY())/h_rp0->Integral(1,h_rp0->GetNbinsX(),1,h_rp0->GetNbinsY())<<endl;
    cout<<"--> Ratio RP1/RP0 CRY2 in AM: "<<h_rp1_am->Integral(1,h_rp1_am->GetNbinsX(),1,h_rp1_am->GetNbinsY())/h_rp0_am->Integral(1,h_rp0_am->GetNbinsX(),1,h_rp0_am->GetNbinsY())<<endl;

    scale2Dhisto(h_rp1, integral_rp0, integral_err_rp0);
    scale2Dhisto(h_rp0, integral_rp0, integral_err_rp0);

    scale2Dhisto(h_rp1_am, integral_rp0_am, integral_err_rp0_am);
    scale2Dhisto(h_rp0_am, integral_rp0_am, integral_err_rp0_am);

    h_rp1->GetXaxis()->SetRange(1,h_rp1->GetNbinsX()/2);
    h_rp1->GetYaxis()->SetRange(1,h_rp1->GetNbinsY()/2);
    h_rp1->GetXaxis()-> CenterTitle();
    h_rp1->GetYaxis()-> CenterTitle();

    h_rp0->GetXaxis()->SetRange(1,h_rp0->GetNbinsX()/2);
    h_rp0->GetYaxis()->SetRange(1,h_rp0->GetNbinsY()/2);
    h_rp0->GetXaxis()-> CenterTitle();
    h_rp0->GetYaxis()-> CenterTitle();

    h_rp1_am->GetXaxis()->SetRange(1,h_rp1_am->GetNbinsX()/2);
    h_rp1_am->GetYaxis()->SetRange(1,h_rp1_am->GetNbinsY()/2);
    h_rp1_am->GetXaxis()-> CenterTitle();
    h_rp1_am->GetYaxis()-> CenterTitle();

    h_rp0_am->GetXaxis()->SetRange(1,h_rp0_am->GetNbinsX()/2);
    h_rp0_am->GetYaxis()->SetRange(1,h_rp0_am->GetNbinsY()/2);
    h_rp0_am->GetXaxis()-> CenterTitle();
    h_rp0_am->GetYaxis()-> CenterTitle();

    TCanvas* c_1 = new TCanvas("c_1","c_1",1000,1000);
    c_1->Divide(2,2);
    c_1->cd(1);
    h_rp0->SetTitle("RP0 CRY2 in CH");
    h_rp0->SetMaximum(1.4e-3);
    h_rp0->Draw("colz");
    c_1->cd(2);
    h_rp1->SetTitle("RP1 CRY2 in CH");
    h_rp1->SetMaximum(1.4e-2);
    h_rp1->Draw("colz");
    c_1->cd(3);
    h_rp0_am->SetTitle("RP0 CRY2 in AM");
    h_rp0_am->SetMaximum(1.4e-3);
    h_rp0_am->Draw("colz");
    c_1->cd(4);
    h_rp1_am->SetTitle("RP1 CRY2 in AM");
    h_rp1_am->SetMaximum(1.4e-2);
    h_rp1_am->Draw("colz");

    TH1D* h_rp1_x_1 = h_rp1->ProjectionX("h_rp1_x_1");
    TH1D* h_rp1_x_2 = h_rp1_am->ProjectionX("h_rp1_x_2");
    TH1D* h_rp1_x_3 = h_rp1->ProjectionX("h_rp1_x_3");

    TCanvas* c_2 = new TCanvas("c_2","c_2",1000,500);
    c_2->cd();
    h_rp1_x_1->SetMinimum(0);
    h_rp1_x_1->SetLineColor(kBlue);
    h_rp1_x_1->SetLineWidth(2);
    h_rp1_x_1->Draw("hist");
    h_rp1_x_2->SetLineColor(kRed);
    h_rp1_x_2->SetLineWidth(2);
    h_rp1_x_2->Draw("same & hist");


    return 0;
}

int function_8(Int_t i)
{
    const Int_t nSets = 13;
//    Int_t i = 5; // setID-1: 0,3,4,5,9,10,11,12 -- good

    Double_t mg_max[] = {/*0*/110.0,/*1*/110.0,/*2*/110.0,/*3*/110.0,/*4*/110.0,/*5*/110.0,/*6*/500.0,/*7*/500.0,/*8*/500.0,/*9*/70.00,/*10*/70.00,/*11*/70.00,/*12*/70.00};
    Double_t mg_min[] = {/*0*/20.00,/*1*/20.00,/*2*/20.00,/*3*/20.00,/*4*/20.00,/*5*/20.00,/*6*/0.000,/*7*/0.000,/*8*/0.000,/*9*/20.00,/*10*/20.00,/*11*/20.00,/*12*/20.00};

    TString fileName_RP0[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP0I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP0I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP0I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP0I_RUN_6.root"
    };

    TString fileName_RP1[] = {
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S1_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S2_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S3_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S4_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S5_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S6_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S7_HISTO_RP1I_RUN_2.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S8_HISTO_RP1I_RUN_5.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S9_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S10_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S11_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S12_HISTO_RP1I_RUN_6.root",
        "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_S13_HISTO_RP1I_RUN_6.root"
    };

    Double_t fit_ch_lim_min[]   = {/*0*/7.00,/*1*/9.00,/*2*/9.00,/*3*/9.00,/*4*/9.00,/*5*/9.00,/*6*/9.00,/*7*/0.00,/*8*/0.00,/*9*/7.00,/*10*/9.00,/*11*/9.00,/*12*/9.00};
    Double_t fit_ch_lim_max[]   = {/*0*/9.00,/*1*/11.0,/*2*/11.0,/*3*/11.0,/*4*/11.0,/*5*/11.0,/*6*/11.0,/*7*/14.0,/*8*/14.0,/*9*/9.00,/*10*/11.0,/*11*/11.0,/*12*/11.0};
    Double_t fit_exp_lim_min[]  = {/*0*/12.0,/*1*/13.0,/*2*/13.0,/*3*/13.0,/*4*/13.0,/*5*/13.0,/*6*/13.0,/*7*/0.00,/*8*/0.00,/*9*/13.0,/*10*/13.0,/*11*/13.0,/*12*/13.0};
    Double_t fit_exp_lim_max[]  = {/*0*/13.8,/*1*/13.8,/*2*/13.8,/*3*/13.8,/*4*/13.8,/*5*/13.8,/*6*/13.8,/*7*/0.00,/*8*/0.00,/*9*/13.8,/*10*/13.8,/*11*/13.8,/*12*/13.8};
    Int_t expo_fit_bkg_status[] = {/*0*/1,   /*1*/1,   /*2*/0,   /*3*/1,   /*4*/1,   /*5*/1,   /*6*/1,   /*7*/0,   /*8*/0,   /*9*/1,   /*10*/1,   /*11*/1,   /*12*/1};       // 0 -- fit gauss+expo ; 1 -- fit expo ; -1 -- w/o fit
    Double_t fit_dch_lim_min[]  = {/*0*/2.50,/*1*/4.50,/*2*/4.50,/*3*/4.50,/*4*/4.50,/*5*/4.50,/*6*/0.50,/*7*/2.50,/*8*/2.50,/*9*/2.50,/*10*/4.50,/*11*/4.50,/*12*/4.50};
    Double_t fit_dch_lim_max[]  = {/*0*/4.00,/*1*/6.50,/*2*/6.50,/*3*/6.50,/*4*/6.50,/*5*/6.50,/*6*/6.50,/*7*/4.00,/*8*/4.00,/*9*/4.00,/*10*/6.50,/*11*/6.50,/*12*/6.50};
    Float_t factor_bad_pixels[] = {/*0*/200.,/*1*/200.,/*2*/200.,/*3*/200.,/*4*/200.,/*5*/30.0,/*6*/4.00,/*7*/100.,/*8*/30.0,/*9*/200.,/*10*/200.,/*11*/100.,/*12*/50.0};
    Double_t par_ch_bkg[nSets][5];

    // Gauss
    par_ch_bkg[0][0] =  0.60; par_ch_bkg[1][0] =  0.60; par_ch_bkg[2][0] =  0.06; par_ch_bkg[3][0] =  0.10;
    par_ch_bkg[4][0] =  0.60; par_ch_bkg[5][0] =  0.60; par_ch_bkg[6][0] =  0.60; par_ch_bkg[7][0] =  0.60;
    par_ch_bkg[8][0] =  0.60; par_ch_bkg[9][0] =  0.60; par_ch_bkg[10][0] = 0.60; par_ch_bkg[11][0] = 0.40;
    par_ch_bkg[12][0] = 0.60;

    par_ch_bkg[0][1] =  7.50; par_ch_bkg[1][1] =  7.50; par_ch_bkg[2][1] =  10.0; par_ch_bkg[3][1] =  7.50;
    par_ch_bkg[4][1] =  7.50; par_ch_bkg[5][1] =  7.50; par_ch_bkg[6][1] =  7.50; par_ch_bkg[7][1] =  10.00;
    par_ch_bkg[8][1] =  7.50; par_ch_bkg[9][1] =  7.50; par_ch_bkg[10][1] = 7.50; par_ch_bkg[11][1] = 10.00;
    par_ch_bkg[12][1] = 7.50;

    par_ch_bkg[0][2] =  0.80; par_ch_bkg[1][2] =  0.80; par_ch_bkg[2][2] =  0.70; par_ch_bkg[3][2] =  0.70;
    par_ch_bkg[4][2] =  0.80; par_ch_bkg[5][2] =  0.80; par_ch_bkg[6][2] =  0.80; par_ch_bkg[7][2] =  0.80;
    par_ch_bkg[8][2] =  0.80; par_ch_bkg[9][2] =  0.80; par_ch_bkg[10][2] = 0.80; par_ch_bkg[11][2] = 1.40;
    par_ch_bkg[12][2] = 0.80;

    // Expo
    par_ch_bkg[0][3] = -4.40; par_ch_bkg[1][3] = -3.70; par_ch_bkg[2][3] = -4.16; par_ch_bkg[3][3] = -3.80;
    par_ch_bkg[4][3] = -2.50; par_ch_bkg[5][3] = -6.00; par_ch_bkg[6][3] = -2.50; par_ch_bkg[7][3] = -2.50;
    par_ch_bkg[8][3] = -2.50; par_ch_bkg[9][3] = -2.50; par_ch_bkg[10][3] =-2.50; par_ch_bkg[11][3] =-1.60;
    par_ch_bkg[12][3] =-2.50;

    par_ch_bkg[0][4] = -0.20; par_ch_bkg[1][4] = -0.20; par_ch_bkg[2][4] = -0.20; par_ch_bkg[3][4] = -0.17;
    par_ch_bkg[4][4] = -0.06; par_ch_bkg[5][4] = -0.02; par_ch_bkg[6][4] = -0.06; par_ch_bkg[7][4] = -0.06;
    par_ch_bkg[8][4] = -0.06; par_ch_bkg[9][4] = -0.06; par_ch_bkg[10][4] =-0.06; par_ch_bkg[11][4] =-0.01;
    par_ch_bkg[12][4] =-0.06;


    TH2D* h_rp0[nSets];
    TH2D* h_rp1[nSets];

    gStyle->SetOptStat(0);

    cout<<endl<<"--> Set "<<i+1<<" <--"<<endl;
    TFile *_file_rp0 = TFile::Open(fileName_RP0[i].Data());
    TFile *_file_rp1 = TFile::Open(fileName_RP1[i].Data());

    h_rp0[i] = (TH2D*)_file_rp0->Get("h_6");
    h_rp1[i] = (TH2D*)_file_rp1->Get("h_6");

    removenoizypixelsXY(h_rp0[i],factor_bad_pixels[i]);
    removenoizypixelsXY(h_rp1[i],factor_bad_pixels[i]);

    Double_t integral_err_rp0;
    Double_t integral_rp0 = h_rp0[i]->IntegralAndError(1,h_rp0[i]->GetNbinsX(),1,h_rp0[i]->GetNbinsY(),integral_err_rp0);
    scale2Dhisto(h_rp0[i], integral_rp0, integral_err_rp0);
    scale2Dhisto(h_rp1[i], integral_rp0, integral_err_rp0);
    cout<<"--> Integral RP1: "<<h_rp1[i]->Integral()<<endl;

    h_rp1[i]->GetXaxis()->SetRange(1,h_rp1[i]->GetNbinsX()/2);
    h_rp1[i]->GetYaxis()->SetRange(1,h_rp1[i]->GetNbinsY()/2);
    h_rp1[i]->GetXaxis()-> CenterTitle();
    h_rp1[i]->GetYaxis()-> CenterTitle();

    TString name = "c_1_";
    name += i;
    TCanvas* c_1 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_1->cd();
    h_rp1[i]->Draw("colz");

    name = "h_rp1_x_";
    name += i+1;
    TH1D* h_rp1_x = h_rp1[i]->ProjectionX(name.Data());
    name += "_s";
    TH1D* h_rp1_x_s = h_rp1[i]->ProjectionX(name.Data());
    name += "_s";
    TH1D* h_rp1_x_s_s = h_rp1[i]->ProjectionX(name.Data());

    name = "c_2_";
    name += i;
    TCanvas* c_2 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_2->cd();
    h_rp1_x->SetMinimum(0);
    h_rp1_x->Draw("hist");
    TF1* fit_funct_ch_bkg;

    if(expo_fit_bkg_status[i] == 1)
    {
        fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","expo",fit_exp_lim_min[i],fit_exp_lim_max[i]);
        fit_funct_ch_bkg->SetParameter(0,par_ch_bkg[i][3]);
        fit_funct_ch_bkg->SetParameter(1,par_ch_bkg[i][4]);
        h_rp1_x->Fit(fit_funct_ch_bkg,"R+");
        par_ch_bkg[i][3] = fit_funct_ch_bkg->GetParameter(0);
        par_ch_bkg[i][4] = fit_funct_ch_bkg->GetParameter(1);
    }
    else
    {
        fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","gaus(0)+expo(3)",fit_ch_lim_min[i],fit_exp_lim_max[i]);
        fit_funct_ch_bkg->SetParameters(par_ch_bkg[i]);
        h_rp1_x->Fit(fit_funct_ch_bkg,"R+");
        fit_funct_ch_bkg->GetParameters(par_ch_bkg[i]);
    }

    name = "#chi^{2}/NDF = [";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare();
    name += "/";
    name += fit_funct_ch_bkg->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare()/fit_funct_ch_bkg->GetNDF();
    h_rp1_x->SetTitle(name.Data());
    fit_funct_ch_bkg->SetLineWidth(6);
    if(expo_fit_bkg_status[i] >= 0) fit_funct_ch_bkg->Draw("same");
    TF1* bkg_func = new TF1("bkg_func","expo",h_rp1_x->GetBinCenter(1),h_rp1_x->GetBinCenter(h_rp1_x->GetNbinsX()));
    bkg_func->SetParameter(0,par_ch_bkg[i][3]);
    bkg_func->SetParameter(1,par_ch_bkg[i][4]);
    bkg_func->SetLineColor(kGreen+2);
    if(expo_fit_bkg_status[i] >= 0) bkg_func->Draw("same");

    name = "c_3_";
    name += i;
    TCanvas* c_3 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_3->cd();
    h_rp1_x_s->SetMinimum(0);
    h_rp1_x_s->SetLineColor(kBlue);
    h_rp1_x_s->SetLineWidth(2);
    h_rp1_x_s->Draw("hist");
    for(Int_t j = 1; j <= h_rp1_x_s_s->GetNbinsX(); j++)
    {
        Double_t val;
        if(expo_fit_bkg_status[i] >= 0)
            val = h_rp1_x_s->GetBinContent(j) - bkg_func->Eval(h_rp1_x_s->GetBinCenter(j));
        else
            val = h_rp1_x_s->GetBinContent(j);



        if(val > 0)
            h_rp1_x_s_s->SetBinContent(j,val);
        else
            h_rp1_x_s_s->SetBinContent(j,0);
    }
    h_rp1_x_s_s->SetLineColor(kRed);
    h_rp1_x_s_s->SetLineWidth(2);
    h_rp1_x_s_s->Draw("same & hist");

    name = "c_4_";
    name += i;
    TCanvas* c_4 = new TCanvas(name.Data(),name.Data(),1000,500);
    c_4->cd();
    h_rp1_x_s_s->SetLineWidth(8);
    h_rp1_x_s_s->Draw("hist");

    TF1* fit_funct_1 = new TF1("fit_funct_1","gaus",fit_ch_lim_min[i],fit_ch_lim_max[i]);
    h_rp1_x_s_s->Fit(fit_funct_1,"R0Q+");
    Double_t par_ch[3];
    fit_funct_1->GetParameters(par_ch);
    TF1* fit_funct_ch = new TF1("fit_funct_ch","gaus",h_rp1_x_s_s->GetBinCenter(1),h_rp1_x_s_s->GetBinCenter(h_rp1_x_s_s->GetNbinsX()));
    fit_funct_ch->SetParameters(par_ch);
    fit_funct_ch->SetLineColor(kBlue);
    fit_funct_ch->SetLineWidth(2);
//    fit_funct_ch->Draw("same");

    TF1* fit_funct_2 = new TF1("fit_funct_2","expo",fit_dch_lim_min[i],fit_dch_lim_max[i]);
    h_rp1_x_s_s->Fit(fit_funct_2,"R0Q+");
    Double_t par_dch[2];
    fit_funct_2->GetParameters(par_dch);
    TF1* fit_funct_dch = new TF1("fit_funct_dch","expo",h_rp1_x_s_s->GetBinCenter(1),h_rp1_x_s_s->GetBinCenter(h_rp1_x_s_s->GetNbinsX()));
    fit_funct_dch->SetParameters(par_dch);
    fit_funct_dch->SetLineColor(kBlack);
    fit_funct_dch->SetLineWidth(2);
//    fit_funct_dch->Draw("same");

    TLine *l_min_ch=new TLine(fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)-3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_min_ch->SetLineColor(kMagenta);
    l_min_ch->SetLineWidth(3);
    l_min_ch->SetLineStyle(7);
    l_min_ch->Draw();
    TLine *l_mean_ch=new TLine(fit_funct_ch->GetParameter(1),0,fit_funct_ch->GetParameter(1),h_rp1_x_s_s->GetMaximum());
    l_mean_ch->SetLineColor(kMagenta-2);
    l_mean_ch->SetLineWidth(3);
    l_mean_ch->SetLineStyle(7);
    l_mean_ch->Draw();
    TLine *l_max_ch=new TLine(fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),0,fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_max_ch->SetLineColor(kMagenta);
    l_max_ch->SetLineWidth(3);
    l_max_ch->SetLineStyle(7);
    l_max_ch->Draw();

    TLine *l_min_dch=new TLine(fit_dch_lim_min[i],0,fit_dch_lim_min[i],h_rp1_x_s_s->GetMaximum());
    l_min_dch->SetLineColor(kGreen+1);
    l_min_dch->SetLineWidth(3);
    l_min_dch->SetLineStyle(7);
    l_min_dch->Draw();
    TLine *l_max_dch=new TLine(fit_dch_lim_max[i],0,fit_dch_lim_max[i],h_rp1_x_s_s->GetMaximum());
    l_max_dch->SetLineColor(kGreen+1);
    l_max_dch->SetLineWidth(3);
    l_max_dch->SetLineStyle(7);
    l_max_dch->Draw();

    TF1* fit_funct_ch_dch = new TF1("fit_funct_ch_dch",fit_ch_dch,fit_dch_lim_min[i],fit_funct_ch->GetParameter(1)+3.0*fit_funct_ch->GetParameter(2),5);
    fit_funct_ch_dch->SetParameter(0,par_ch[0]);
    fit_funct_ch_dch->SetParameter(1,par_ch[1]);
    fit_funct_ch_dch->SetParameter(2,par_ch[2]);
    fit_funct_ch_dch->SetParameter(3,par_dch[0]);
    fit_funct_ch_dch->SetParameter(4,par_dch[1]);
    fit_funct_ch_dch->SetLineColor(kBlack);
    fit_funct_ch_dch->SetLineWidth(4);
    fit_funct_ch_dch->Draw("same");

    name = "#chi^{2}/NDF (DCH) = [";
    name += (Int_t)fit_funct_2->GetChisquare();
    name += "/";
    name += fit_funct_2->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_2->GetChisquare()/fit_funct_2->GetNDF();
    name += " | #chi^{2}/NDF (CH) = [";
    name += (Int_t)fit_funct_1->GetChisquare();
    name += "/";
    name += fit_funct_1->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_1->GetChisquare()/fit_funct_1->GetNDF();
    h_rp1_x_s_s->SetTitle(name.Data());

    gPad->Modified();
    name = "profile_ch_dch_";
    name += i;
    name += ".png";
    c_4->SaveAs(name.Data());

    //=====================================//
    // To calculate using histogram data
    //=====================================//
    const Int_t nSteps = 125;
    const Double_t min_x_pos = 0.0;// mm
    Double_t xi_l, xi_r, xi_err[nSteps] = {}, xi[nSteps] = {}, Ndch[nSteps] = {}, Ndch_err[nSteps] = {}, Nch[nSteps] = {}, Nch_err[nSteps] = {},
            R[nSteps] = {}, R_err[nSteps] = {};

    Int_t bin_start = h_rp1_x_s_s->FindBin(min_x_pos), dn_bins = 4;

    for(Int_t ii = 0; ii < nSteps; ii++)
    {
        xi[ii] = min_x_pos + ii*0.055;
        xi_err[ii] = 0.055/TMath::Sqrt(12.0);
        Ndch[ii] = h_rp1_x_s_s->IntegralAndError(bin_start+ii,bin_start+ii+dn_bins,Ndch_err[ii]);
        Nch[ii] = h_rp1_x_s_s->IntegralAndError(bin_start+ii+dn_bins+1,h_rp1_x_s_s->FindBin(h_rp1_x_s_s->GetNbinsX()),Nch_err[ii]);

        Ndch[ii] /= (dn_bins+1)*0.055;
        Ndch_err[ii] /= (dn_bins+1)*0.055;
        if(Ndch[ii] > 0)
        {
            R[ii] = Nch[ii]/Ndch[ii];
            R_err[ii] = TMath::Sqrt(TMath::Power(Nch_err[ii]/Ndch[ii],2) + TMath::Power(Nch[ii]*Ndch_err[ii]/(Ndch[ii]*Ndch[ii]),2));
        }
        else
        {
            R[ii] = 0.0;
            R_err[ii] = 0.0;
        }
    }
    //=====================================//
    // To calculate analytically
    //=====================================//
    Double_t Ndch_a[nSteps] = {}, Nch_a[nSteps] = {}, R_a[nSteps] = {};

    for(Int_t ii = 0; ii < nSteps; ii++)
    {
        xi_l = min_x_pos + ii*0.055;
        xi_r = min_x_pos + (ii+dn_bins)*0.055;
        Ndch_a[ii] = fit_funct_ch_dch->Integral(xi_l,xi_r);
        Nch_a[ii] = fit_funct_ch_dch->Integral(xi_r,N_PIXELS*0.055/2);

        Ndch_a[ii] /= (dn_bins)*0.055;
        if(Ndch_a[ii] > 0)
        {
            R_a[ii] = Nch_a[ii]/Ndch_a[ii];
        }
        else
        {
            R_a[ii] = 0.0;
        }
    }
    //=====================================//
    // Ratio
    //=====================================//
    name = "c_5_";
    name += i;
    TCanvas* c_5 = new TCanvas(name.Data(),name.Data(),1000,1000);
    c_5->cd();
    gPad->SetGrid();
    TGraphErrors* gr_h = new TGraphErrors(nSteps,xi,R,xi_err,R_err);
    TGraphErrors* gr_a = new TGraphErrors(nSteps,xi,R_a,0,0);
    TMultiGraph* mg = new TMultiGraph();
    gr_h->SetLineWidth(2);
    gr_a->SetLineWidth(2);
    gr_a->SetLineColor(kBlack);
    gr_h->SetLineColor(kRed);
    gr_h->SetMarkerColor(kRed);
    mg->Add(gr_h,"AP");
    mg->Add(gr_a,"AL");
    mg->Draw("APL");
    gPad->Modified();
    mg->GetYaxis()->SetTitle("R^{CH/DCH} [per mm of septum width]");
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetXaxis()->SetTitle("X [mm]");
    mg->GetYaxis()->CenterTitle(1);
    mg->GetXaxis()->CenterTitle(1);
    mg->GetXaxis()->SetLimits(fit_dch_lim_min[i],fit_dch_lim_max[i]);
    mg->SetMaximum(mg_max[i]);
    mg->SetMinimum(mg_min[i]);
    gPad->Modified();
    name = "ratio_ch_dch_";
    name += i;
    name += ".png";
    c_5->SaveAs(name.Data());

    return 0;
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
    cout<<"--> nMaskedPixels = "<<nMaskedPixels<<endl;
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

void median_filter(TH2D* h_in, TH2D* h_out, Double_t factor_bad_pixels)
{
    const Int_t nBins = 3;
    if (nBins % 2 == 0)
    {
        cout<<"ERROR:: nBins should be odd!!!"<<endl;
        assert(0);
    }

    Double_t a[nBins*nBins];
    for(Int_t i = 1 + (nBins-1)/2; i <= h_in->GetNbinsX() - (nBins-1)/2; i++)
    {
        for(Int_t j = 1 + (nBins-1)/2; j <= h_in->GetNbinsY() - (nBins-1)/2; j++)
        {
            Int_t kk = 0;
            for(Int_t k = i - (nBins-1)/2; k<= i + (nBins-1)/2; k++)
            {
                for(Int_t l = j - (nBins-1)/2; l<= j + (nBins-1)/2; l++)
                {
                    a[kk] = h_in->GetBinContent(k,l);
                    kk++;
                }
            }
            Double_t val = findMedian(a, nBins*nBins);
            if(h_in->GetBinContent(i,j) < val*factor_bad_pixels)
                val = h_in->GetBinContent(i,j);

            h_out->SetBinContent(i,j,val);
        }
    }

    for(Int_t i = 1; i <= h_in->GetNbinsX()/2; i++)
    {
        for(Int_t j = 1; j <= h_in->GetNbinsY()/2; j++)
        {
            if(i == 1 || i == h_in->GetNbinsX()/2 || j == 1 || j == h_in->GetNbinsY()/2)
            {
                h_out->SetBinContent(i,j,0);
            }
        }
    }
}

double findMedian(double a[], int n)// Function for calculating median
{
    sort(a, a+n);
    return a[n/2];
}
