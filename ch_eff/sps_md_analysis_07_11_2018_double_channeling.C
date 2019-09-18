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
#include "TF2.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TFitResultPtr.h"
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

const Int_t N_PIXELS                = 512/2;
const Double_t PIXEL_SIZE           = 0.055;

void repairnoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels);
void repairnoizypixelsYT(TH2D* histo, Float_t factor_bad_pixels);
void scale1Dhisto(TH1D* histo, Double_t integral, Double_t integral_err);
void scale2Dhisto(TH2D* histo, Double_t integral, Double_t integral_err);
void scale2Dhisto_YT(TH2D* histo);
void rotateNscale(TH2D* h_in, TH2D* h_out);
void rewriteHisto(TH2D* h_in, TH2D* h_out, TH1D* h_fills);
double getAverage2D(TH2D* histo);
int median_filter(TH2D* h_in, TH2D* h_out, Double_t factor_bad_pixels);
double findMedian(double a[], int n);

int sps_md_analysis_07_11_2018_double_channeling()
{
    cout<<"--> function_1() -- to plot rp_1 2d beam image normalized by counts"<<endl;
    cout<<"--> function_2() -- to plot rp_1 projection"<<endl;
    return 0;
}

void function_1()
{
    //---------------------------------------------------------------------------------------//
    Int_t chipID = 1;// (1 - RP1I; 3 - RP0I)
    TChain* fChain;
    Long64_t nEntries;

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    UInt_t _AcquisType = -1, _TrigType = -1;
    Long64_t _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS*2][N_PIXELS*2];

    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_11_07_RUN_1.root";
    TString outFileName     = "output_chipid_1_MD_2018_11_07_RUN_1.root";

    TString tree_name;
    tree_name = "Tree_";
    tree_name += chipID;
    cout<<"--> Tree: "<<tree_name<<endl;
    fChain = new TChain(tree_name.Data());
    fChain->Add(inFileName1.Data());

    fChain->SetBranchAddress("UnixTime",    &_UnixTime);
    fChain->SetBranchAddress("TrigType",    &_TrigType);
    fChain->SetBranchAddress("AcquisType",  &_AcquisType);
    fChain->SetBranchAddress("DeltaTHR",    &_DeltaTHR);
    fChain->SetBranchAddress("Ikrum",       &_Ikrum);
    fChain->SetBranchAddress("Bias",        &_Bias);
    fChain->SetBranchAddress("Clock",       &_Clock);
    fChain->SetBranchAddress("Gate",        &_Gate);
    fChain->SetBranchAddress("Event",       &_event);
    fChain->SetBranchAddress("Time_ms",     &_Timems);
    fChain->SetBranchAddress("COUNTS",      _COUNTS);

    nEntries = (fChain->GetEntries());
    cout<<"--> ChipID: "<<chipID<<" (1 - RP1I; 3 - RP0I)"<<endl;
    cout<<"--> InputFileName: "<<inFileName1<<endl;
    cout<<"--> Number of Frames: "<<nEntries<<endl;
    //------------------------------------------------------------------------------//
    //------------------------------- HISTOGRAMS -----------------------------------//
    //------------------------------------------------------------------------------//

    fChain->GetEntry(nEntries-1);
    Double_t UT_max = _Timems;

    gStyle->SetOptStat(0);
    fChain->GetEntry(0);
    Double_t UT_min = _Timems;
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    TH2D* h_1 = new TH2D("h_1","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_2 = new TH2D("h_2","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_4   = new TH2D("h_4","RP1 Internal (normalized)",N_PIXELS,0,N_PIXELS*0.055,N_PIXELS,0,N_PIXELS*0.055);
    TH2D* h_8   = new TH2D("h_8","Y vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);

    h_4->GetXaxis()->SetTitle("Horizontal Axis [mm]");
    h_4->GetXaxis()->SetTimeOffset(1.2);
    h_4->GetXaxis()->CenterTitle();
    h_4->GetYaxis()->SetTitle("Vertical Axis [mm]");
    h_4->GetYaxis()->SetTimeOffset(1.2);
    h_4->GetYaxis()->CenterTitle();

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------- MAIN LOOP ------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//    
    for(Long64_t i = 0; i < nEntries; i++)
    {
        if(i%1 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntries);
            fflush(stdout);
        }

        fChain->GetEntry(i);

        if(_Timems < 1541594028000 || _Timems > 1541594036000) continue;

        for(Int_t xi = 1; xi < N_PIXELS - 1; xi++)
        {
            for(Int_t yi = 1; yi < N_PIXELS - 1; yi++)
            {
                if(_COUNTS[xi][yi] > 0 && _AcquisType == 0)
                {
                    h_1->Fill(xi,yi,_COUNTS[xi][yi]);
                    h_2->Fill(xi,yi,_COUNTS[xi][yi]);
                    h_8->Fill(_Timems/1000.0,yi);
                }
            }
        }
    }

    cout<<endl;
    cout<<"--> Filtered pixels = "<<median_filter(h_1,h_2,2.0)<<endl;;

    //----------------------------------------------------//
    // To check the algorithm
    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->Divide(2,1);

    c_1->cd(1);
    gPad->SetGrid();
    h_1->GetXaxis()->SetRange(1,256);
    h_1->GetYaxis()->SetRange(1,256);
    h_1->Draw("colz");

    c_1->cd(2);
    gPad->SetGrid();
    h_2->GetXaxis()->SetRange(1,256);
    h_2->GetYaxis()->SetRange(1,256);
    h_2->Draw("colz");
    //----------------------------------------------------//

    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File '"<<_file->GetName()<<"' opened successfully!"<<endl;

    Double_t integral, integral_err;
    integral = h_2->IntegralAndError(1,h_2->GetNbinsX(),1,h_2->GetNbinsY(),integral_err);
    scale2Dhisto(h_2,integral,integral_err);

    rotateNscale(h_2,h_4); // [pixels] --> [mm]

    TCanvas* c_2 = new TCanvas("c_2","c_2",900,900);
    c_2->cd();
    gPad->SetGrid();
    h_4->GetYaxis()->SetRange(1,256);
    h_4->GetXaxis()->SetRange(1,256);
    h_4->Draw("colz");

    TCanvas* c_3 = new TCanvas("c_3","c_3",900,900);
    c_3->cd();
    gPad->SetGrid();
    h_8->Draw("colz");

    c_1->Write();
    c_2->Write();
    h_4->Write();
    h_8->Write();

    _file->Close();
    //---------------------------------------------------------------------------------------//
}

void function_2()
{
    TFile *_file0 = TFile::Open("output_chipid_1_MD_2018_11_07_RUN_1.root");
    TH2D* hh = (TH2D*)_file0->Get("h_4");
    TH1D* hh_px = hh->ProjectionX("hh_px");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1000,1000);
    c_1->cd();
    hh->Draw("colz");

    TCanvas* c_2 = new TCanvas("c_2","c_2",1000,500);
    c_2->cd();
    hh_px->SetMarkerStyle(21);
    hh_px->SetMarkerSize(0.5);
    hh_px->SetTitle("");
    hh_px->Draw("e1");

    TF1* fit1 = new TF1("fit1","gaus",12.5,13.5);
    TF1* fit2 = new TF1("fit2","gaus",5.0,10.0);

    fit1->SetLineColor(kBlack);
    fit2->SetLineColor(kBlue);

    hh_px->Fit(fit1,"R+");
    hh_px->Fit(fit2,"R+");


}

void repairnoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels)// to repair noizy pixels for XY image
{
    Int_t nMaskedPixels = 0;
    for(Int_t xi = 1; xi <= histo->GetNbinsX()/2; xi++)
    {
        for(Int_t yi = 1; yi < histo->GetNbinsY()/2; yi++)
        {
            if(histo->GetBinContent(xi,yi) > factor_bad_pixels*histo->GetBinContent(xi,yi+1))
            {
                histo->SetBinContent(xi,yi,histo->GetBinContent(xi,yi+1));
                histo->SetBinError(xi,yi,histo->GetBinError(xi,yi+1));
                nMaskedPixels++;
            }
        }
    }
//    cout<<"--> nMaskedPixels = "<<nMaskedPixels<<endl;
}

void repairnoizypixelsYT(TH2D* histo, Float_t factor_bad_pixels)// to repair noizy pixels for YvsTime image
{
    Int_t nMaskedPixels = 0;
    for(Int_t xi = 1; xi <= histo->GetNbinsX(); xi++)
    {
        for(Int_t yi = 1; yi < histo->GetNbinsY()/2; yi++)
        {
            if(histo->GetBinContent(xi,yi) > factor_bad_pixels*histo->GetBinContent(xi,yi+1))
            {
                histo->SetBinContent(xi,yi,histo->GetBinContent(xi,yi+1));
                histo->SetBinError(xi,yi,histo->GetBinError(xi,yi+1));
                nMaskedPixels++;
            }
        }
    }
    cout<<"--> nMaskedPixels = "<<nMaskedPixels<<endl;
}

void scale2Dhisto(TH2D* histo, Double_t integral, Double_t integral_err)// to scale the 2D histo
{
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= histo->GetNbinsY(); j++)
        {
            Double_t val = histo->GetBinContent(i,j);
            Double_t val_err = histo->GetBinError(i,j);
            val_err = TMath::Sqrt(TMath::Power(val_err/integral,2) + TMath::Power(val*integral_err/(integral*integral),2));
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

void scale2Dhisto_YT(TH2D* histo)// to scale the 2D histo Y vs Time
{
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        Double_t integral, integral_err;
        integral = histo->IntegralAndError(i,i,1,histo->GetNbinsY(),integral_err);

        for(Int_t j = 1; j <= histo->GetNbinsY(); j++)
        {
            Double_t val = histo->GetBinContent(i,j);
            Double_t val_err = histo->GetBinError(i,j);
            val_err = TMath::Sqrt(TMath::Power(val_err/integral,2) + TMath::Power(val*integral_err/(integral*integral),2));
            val = val/integral;

            histo->SetBinContent(i,j,val);
            histo->SetBinError(i,j,val_err);
        }
    }
}

void rotateNscale(TH2D* h_in, TH2D* h_out)// to rotate and scale from pixels to mm for SPS RomanPot Internal
{
    // For 512x512
//    for(Int_t i = 1; i <= h_in->GetNbinsX(); i++)
//    {
//        for(Int_t j = 1; j <= h_in->GetNbinsY(); j++)
//        {
//            h_out->SetBinContent(h_in->GetNbinsX()/2-j+1,i,h_in->GetBinContent(i,j));
//            h_out->SetBinError(h_in->GetNbinsX()/2-j+1,i,h_in->GetBinError(i,j));
//        }
//    }
    // For 256x256
    for(Int_t i = 1; i <= h_in->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h_in->GetNbinsY(); j++)
        {
            h_out->SetBinContent(h_in->GetNbinsX()-j+1,i,h_in->GetBinContent(i,j));
            h_out->SetBinError(h_in->GetNbinsX()-j+1,i,h_in->GetBinError(i,j));
        }
    }
}

void rewriteHisto(TH2D* h_in, TH2D* h_out, TH1D *h_fills)// to rewrite a histogram
{
    Int_t nFills = 0;
    for(Int_t i = 1; i <= h_in->GetNbinsX(); i++)
    {
        nFills = h_fills->GetBinContent(i);
        if(nFills < 1) nFills = 1;

        for(Int_t j = 1; j <= h_in->GetNbinsY(); j++)
        {
            h_out->SetBinContent(i,j,h_in->GetBinContent(i,j)/nFills);
            h_out->SetBinError(i,j,h_in->GetBinError(i,j)/nFills);
        }
    }
}

double getAverage2D(TH2D* histo)// to get average value of the counts 2D
{
    Double_t average = 0;
    Int_t nBins = 0;
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= histo->GetNbinsY(); j++)
        {
            if(histo->GetBinContent(i,j) > 0)
            {
                average += histo->GetBinContent(i,j);
                nBins++;
            }
        }
    }
    return (Double_t)(average/nBins);
}

int median_filter(TH2D* h_in, TH2D* h_out, Double_t factor_bad_pixels)
{
    Int_t filtered_pixels = 0;
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
            filtered_pixels++;
            if(h_in->GetBinContent(i,j) < val*factor_bad_pixels && h_in->GetBinContent(i,j) > val/factor_bad_pixels)
            {
                val = h_in->GetBinContent(i,j);
                filtered_pixels--;
            }

            h_out->SetBinContent(i,j,val);
        }
    }
    return filtered_pixels;
}

double findMedian(double a[], int n)// Function for calculating median
{
    sort(a, a+n);
    return a[n/2];
}
