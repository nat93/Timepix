//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2019
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
#include "TProfile.h"
#include "TVirtualFitter.h"
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

const Int_t N_PIXELS                = 256;
const Double_t PIXEL_SIZE           = 0.055;

void function_1();
void function_2();
void function_3();
void function_4();
void function_5();
void function_6();
void function_7();
void function_8();

int median_filter(TH2D* h_in, TH2D* h_out, Double_t factor_bad_pixels);
double findMedian(double a[], int n);


int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout<<endl;
        cout<<">>> INPUT PARAMETERS <<<"<<endl;
        cout<<"[0] - script name"<<endl;
        cout<<"[1] - function number:"<<endl;
        cout<<"         1 --> function_1() - H8 Plot Frames"<<endl;
        cout<<endl;
    }
    else
    {
        switch (atoi(argv[1]))
        {
        case 1:
            function_1();
            break;
        default:
            cout<<endl<<"### NOTHING TO DO ;) ###"<<endl<<endl;
            break;
        }
    }
    return 0;
}

void function_1()
{
    cout<<endl<<">>> function_1() <<<"<<endl;

    // CRYSTAL3 Fine Linear Scan
    // 2018_09_17
    Long64_t minUnixTime_run    = 1537218019;
    Long64_t maxUnixTime_run    = 1537218120;
    Long64_t entryINI           = 12700;

    TString outputFileName      = "output_function_1_RP1.root";
    Int_t chipID = 0;

    //********************************************************//
    //********************************************************//

    Bool_t entryINIstatus       = kFALSE;
    //------------------------------------------------------------------------------//
    // READ TREE
    //------------------------------------------------------------------------------//
    TChain* fChain;
    Long64_t nEntries;

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    UInt_t _AcquisType = -1, _TrigType = -1;
    Long64_t _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS*2][N_PIXELS*2];

    TString _fileName   = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_RUN_8.root";

    TString tree_name;
    tree_name = "Tree_";
    tree_name += chipID;
    cout<<"--> Tree: "<<tree_name<<endl;
    fChain = new TChain(tree_name.Data());
    fChain->Add(_fileName.Data());

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
    cout<<"--> ChipID: "<<chipID<<" (0 B05-W0108-D)"<<endl;
    cout<<"--> InputFileName: "<<_fileName<<endl;
    cout<<"--> Number of Frames: "<<nEntries<<endl;

    //------------------------------------------------------------------------------//
    // CREATE HISTOGRAMS
    //------------------------------------------------------------------------------//

    fChain->GetEntry(0);
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    TH2D* h_1   = new TH2D("h_1","Hits XY [pixels]",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

    TCanvas* c_3 = new TCanvas("c_3","c_3",1800,900);
    c_3->Divide(2,1);
    //------------------------------------------------------------------------------//
    // MAIN LOOP
    //------------------------------------------------------------------------------//

    Int_t nFrames = 0;
    Double_t frameTime = 0;

    for(Long64_t iEntry = entryINI; iEntry < nEntries; iEntry++)
    {
        if(iEntry%1 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)iEntry/nEntries);
            fflush(stdout);
        }

        fChain->GetEntry(iEntry);

        frameTime = _Timems/1000.0; // [sec]

        if(frameTime < minUnixTime_run) continue;
        if(frameTime > maxUnixTime_run) break;

        if(!entryINIstatus)
        {
            entryINI = iEntry;
            cout<<endl<<"--> entryINI = "<<entryINI<<endl;
            entryINIstatus = kTRUE;
        }

        TH2D* h_1_temp = new TH2D("h_1_temp","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
        TH2D* h_2_temp = new TH2D("h_2_temp","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

        for(Int_t xi = 1; xi < N_PIXELS - 1; xi++)
        {
            for(Int_t yi = 1; yi < N_PIXELS - 1; yi++)
            {
                if(_COUNTS[xi][yi] > 0 && _AcquisType == 0)
                {
                    h_1_temp->Fill(xi,yi,_COUNTS[xi][yi]);
                    h_2_temp->Fill(xi,yi,_COUNTS[xi][yi]);
                }
            }
        }

        h_23->Fill(median_filter(h_1_temp,h_2_temp,filterFactor));

        //----------------------------------------------------//
        //----------- To check the algorythm -----------------//
        //----------------------------------------------------//
//        c_3->cd(1);
//        gPad->SetGrid();
//        h_1_temp->GetXaxis()->SetRange(1,256);
//        h_1_temp->GetYaxis()->SetRange(1,256);
//        h_1_temp->Draw("colz");

//        c_3->cd(2);
//        gPad->SetGrid();
//        h_2_temp->GetXaxis()->SetRange(1,256);
//        h_2_temp->GetYaxis()->SetRange(1,256);
//        h_2_temp->Draw("colz");
//        break;
        //----------------------------------------------------//
        //----------------------------------------------------//
        //----------------------------------------------------//

        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0)
                {
                    h_1->Fill(xi,yi,h_2_temp->GetBinContent(xi+1,yi+1));
                    h_2->Fill(frameTime,yi,h_2_temp->GetBinContent(xi+1,yi+1));
                    h_3->Fill(frameTime,xi,h_2_temp->GetBinContent(xi+1,yi+1));
                }
            }
        }

        h_22->Fill(frameTime,1);

        nFrames++;

        h_1_temp->Delete();
        h_2_temp->Delete();
    }
    for(Int_t binx = 1; binx <= N_PIXELS; binx++)
    {
        for(Int_t biny = 1; biny <= N_PIXELS; biny++)
        {
            h_5->SetBinContent(N_PIXELS-biny+1,binx,h_1->GetBinContent(binx,biny));
            h_5->SetBinError(N_PIXELS-biny+1,binx,h_1->GetBinError(binx,biny));
        }
    }

    cout<<endl<<"--> nFrames = "<<nFrames<<endl;
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    //------------------------------------------------------------------------------//
    // WRITE
    //------------------------------------------------------------------------------//

    TFile *_file = new TFile(outputFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File '"<<_file->GetName()<<"' opened successfully!"<<endl;

    c_3->Write();

    _file->Close();
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
            if(h_in->GetBinContent(i,j) <= (val*factor_bad_pixels + 1)  && h_in->GetBinContent(i,j) >= (val/factor_bad_pixels))
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
