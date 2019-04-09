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
        cout<<"         1 --> function_1() - Timepix frames integration (Medipix mode)"<<endl;
        cout<<"         2 --> function_2() - Convert crystal position/orientation or blm counts"<<endl;
        cout<<"         3 --> function_3() - Convert bctdc counts"<<endl;
        cout<<"         4 --> function_4() - BLM, BCTDC vs Crystal angle, Timepix vs Time"<<endl;
        cout<<endl;
    }
    else
    {
        switch (atoi(argv[1]))
        {
        case 1:
            function_1();
            break;
        case 2:
            function_2();
            break;
        case 3:
            function_3();
            break;
        case 4:
            function_4();
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

    //----------------------------------------------------------------------------------------//
    // ChipID: 0 RP0Ext G02-W0108, 1 RP1Int F04-W0108, 2 RP3Ext C08-W0255, 3 RP0Int I02-W0108 //
    //----------------------------------------------------------------------------------------//

    //********************************************************//
    //********************************************************//
    // CRYSTAL2 Angular Scan
    // 2018_09_17
    Long64_t minUnixTime_run    = 1537206118;
    Long64_t maxUnixTime_run    = 1537206480;
    Long64_t entryINI           = 0;

    //========= RomanPot 1 =========//
//    Double_t filterFactor       = 2.0;
//    TString outputFileName      = "output_function_1_RP1.root";
//    Int_t chipID = 1;

    //========= RomanPot 0 =========//
    Double_t filterFactor       = 2.0;
    TString outputFileName      = "output_function_1_RP0.root";
    Int_t chipID = 3;
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
    cout<<"--> ChipID: "<<chipID<<" (0 RP0Ext G02-W0108, 1 RP1Int F04-W0108, 2 RP3Ext C08-W0255, 3 RP0Int I02-W0108)"<<endl;
    cout<<"--> InputFileName: "<<_fileName<<endl;
    cout<<"--> Number of Frames: "<<nEntries<<endl;

    //------------------------------------------------------------------------------//
    // CREATE HISTOGRAMS
    //------------------------------------------------------------------------------//

    fChain->GetEntry(0);
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    TH2D* h_1   = new TH2D("h_1","RP Internal",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_2   = new TH2D("h_2","RP Internal Y vs Time",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run,N_PIXELS,0,N_PIXELS);
    TH2D* h_3   = new TH2D("h_3","RP Internal X vs Time",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run,N_PIXELS,0,N_PIXELS);
    TH2D* h_5   = new TH2D("h_5","RP Internal YX [mm]",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH1D* h_22  = new TH1D("h_22","RP Internal Time (Frames)",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run);
    TH1D* h_23  = new TH1D("h_23","Number of the filtered pixels per frame",10000,0,10000);

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
    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_5->Write();
    h_22->Write();

    _file->Close();
}

void function_2()
{
    TString Inputfilename   = "crystal2_angularscan_2018_09_17.dat";
    TString Outputfilename   = "crystal2_angularscan_2018_09_17.root";

    TFile *file = new TFile(Outputfilename.Data(),"recreate");
    TTree *tree = new TTree("Tree","A Root Tree");

    Double_t value;
    Double_t untime;

    TString value_ss = "value/D";
    TString untime_ss = "untime/D";

    tree->Branch("Value", &value, value_ss.Data());
    tree->Branch("UnixTime", &untime, untime_ss.Data());

    cout<<"--> input file: "<<Inputfilename<<endl;
    ifstream infile(Inputfilename.Data());

    while (infile>>untime)
    {
        infile>>value;

        tree->Fill();
    }

    file->Write();
    file->Close();

    cout<<"--> output file: "<<Outputfilename<<endl;
}

void function_3()
{
    TString Inputfilename   = "bctdc_crystal2_angularscan_2018_09_17.dat";
    TString Outputfilename   = "bctdc_crystal2_angularscan_2018_09_17.root";

    TFile *file = new TFile(Outputfilename.Data(),"recreate");
    TTree *tree = new TTree("Tree","A Root Tree");

    Double_t value;
    Double_t vstd;
    Double_t untime;

    TString value_ss = "value/D";
    TString vstd_ss = "vstd/D";
    TString untime_ss = "untime/D";

    tree->Branch("Value", &value, value_ss.Data());
    tree->Branch("Std", &vstd, vstd_ss.Data());
    tree->Branch("UnixTime", &untime, untime_ss.Data());

    cout<<"--> input file: "<<Inputfilename<<endl;
    ifstream infile(Inputfilename.Data());

    while (infile>>untime)
    {
        infile>>value;
        infile>>vstd;

        tree->Fill();
    }

    file->Write();
    file->Close();

    cout<<"--> output file: "<<Outputfilename<<endl;
}

void function_4()
{
    //=============================================================//
    // CRYSTAL2 Angular Scan
    // 2018_09_17
    TString output_file_name    = "output_function_4_RP1.root";
    TString input_file_tpx      = "output_function_1_RP1.root";
    TString input_file_motor    = "crystal2_angularscan_2018_09_17.root";
    TString input_file_blm      = "blm3_crystal2_angularscan_2018_09_17.root";
    TString input_file_bctdc    = "bctdc_crystal2_angularscan_2018_09_17.root";
    Long64_t minUnixTime_run    = 1537206118;
    Long64_t maxUnixTime_run    = 1537206480;
    Int_t angBins               = 700;
    Double_t angLimMin          = -6700;
    Double_t angLimMax          = -6000;

    //=============================================================//

    // Timepix
    TFile* _fileTpx = TFile::Open(input_file_tpx.Data());
    TH2D* h_1 = (TH2D*)_fileTpx->Get("h_1");
    TH2D* h_2 = (TH2D*)_fileTpx->Get("h_2");
    TH1D* h_22 = (TH1D*)_fileTpx->Get("h_22");

    // BLM
    Double_t blm;
    Double_t untime_blm;

    TChain *fChain2 = new TChain("Tree");
    fChain2->Add(input_file_blm);

    fChain2->SetBranchAddress("Value",  &blm);
    fChain2->SetBranchAddress("UnixTime",       &untime_blm);

    // Motor
    Double_t angle;
    Double_t untime_motor;

    TChain *fChain3 = new TChain("Tree");
    fChain3->Add(input_file_motor);

    fChain3->SetBranchAddress("Value",          &angle);
    fChain3->SetBranchAddress("UnixTime",       &untime_motor);

    // BCTDC
    Double_t bctdc;
    Double_t bctdc_std;
    Double_t untime_bctdc;

    TChain *fChain4 = new TChain("Tree");
    fChain4->Add(input_file_bctdc);

    fChain4->SetBranchAddress("Value",  &bctdc);
    fChain4->SetBranchAddress("Std",  &bctdc_std);
    fChain4->SetBranchAddress("UnixTime",       &untime_bctdc);


    //---------------//
    cout<<"--> Input file with tpx: "<<input_file_tpx<<endl;

    cout<<"--> Input file with blm: "<<input_file_blm<<endl;
    Double_t nEntries_2 = fChain2->GetEntries();
    cout<<"--> nEntries: "<<nEntries_2<<endl;

    cout<<"--> Input file with motor: "<<input_file_motor<<endl;
    Double_t nEntries_3 = fChain3->GetEntries();
    cout<<"--> nEntries: "<<nEntries_3<<endl;

    cout<<"--> Input file with bctdc: "<<input_file_bctdc<<endl;
    Double_t nEntries_4 = fChain4->GetEntries();
    cout<<"--> nEntries: "<<nEntries_4<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//

    TH2D* hh_1   = new TH2D("hh_1","RP1 Internal YX [mm]",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_2   = new TH2D("hh_2","RP1 Internal Y vs Time",maxUnixTime_run-minUnixTime_run,0,maxUnixTime_run-minUnixTime_run,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_10 = new TH2D("hh_10","CrystalAngle vs UnixTime",(maxUnixTime_run-minUnixTime_run),minUnixTime_run,maxUnixTime_run,angBins,angLimMin,angLimMax);

    TGraphErrors* gr_0 = new TGraphErrors();
    gr_0->SetName("gr_0");
    gr_0->SetTitle("Angle");
    Int_t gr_0_iter = 0;

    TGraphErrors* gr_1 = new TGraphErrors();
    gr_1->SetName("gr_1");
    gr_1->SetTitle("Bctdc");
    Int_t gr_1_iter = 0;

    TGraphErrors* gr_2 = new TGraphErrors();
    gr_2->SetName("gr_2");
    gr_2->SetTitle("gradI");
    Int_t gr_2_iter = 0;

    TGraphErrors* gr_3 = new TGraphErrors();
    gr_3->SetName("gr_3");
    gr_3->SetTitle("BlmVsAngle");
    Int_t gr_3_iter = 0;

    //--------------------------------------------------------------------------//

    // MOTOR
    for(Int_t eventID = 0; eventID < nEntries_3; eventID++)
    {
        fChain3->GetEntry(eventID);

        hh_10->Fill(untime_motor,angle);

        gr_0->SetPoint(gr_0_iter,untime_motor,angle);
        gr_0->SetPointError(gr_0_iter,0.001,0.1);
        gr_0_iter++;
    }
    cout<<endl;

    TCanvas* cFit = new TCanvas("cFit","cFit");
    cFit->cd();
    gr_0->Draw("AP");
    TF1* fitPol = new TF1("fitPol","pol1",minUnixTime_run,maxUnixTime_run);
    fitPol->SetParameters(-2.6e9,1.7);
    gr_0->Fit(fitPol,"R+");

    // BCTDC
    for(Int_t eventID = 0; eventID < nEntries_4; eventID++)
    {
        fChain4->GetEntry(eventID);

        gr_1->SetPoint(gr_1_iter,untime_bctdc,bctdc);
        gr_1->SetPointError(gr_1_iter,0.001,bctdc_std);
        gr_1_iter++;
    }
    cout<<endl;
    TCanvas* cc = new TCanvas("cc","cc");
    cc->cd();
    gr_1->Draw("APL");
    fChain4->GetEntry(0);
    Double_t timeINI = untime_bctdc;
    fChain4->GetEntry(nEntries_4-1);
    Double_t timeFIN = untime_bctdc;

    TF1* fit1 = new TF1("fit1","expo",timeINI,timeFIN);
    fit1->SetParameters(1e6,-6.4e-4);
    gr_1->Fit(fit1,"R+");

    TGraphErrors *grint = new TGraphErrors(gr_1->GetN());
    grint->SetTitle("Fitted line with .95 conf. band");
    for (Int_t i = 0; i < gr_1->GetN(); i++)
    {
        grint->SetPoint(i, gr_1->GetX()[i], 0);
    }
    /*Compute the confidence intervals at the x points of the created graph*/
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint);
    //Now the "grint" graph contains function values as its y-coordinates
    //and confidence intervals as the errors on these coordinates
    //Draw the graph, the function and the confidence intervals

    grint->SetLineColor(kRed);
    grint->Draw("Psame");
    cout<<endl;
    // BLM
    fChain2->GetEntry(0);
    Double_t timeOld = untime_blm;

    for(Int_t eventID = 1; eventID < nEntries_2; eventID++)
    {
        fChain2->GetEntry(eventID);

        Double_t angle_temp = fitPol->Eval(untime_blm);
        Double_t gradInt = TMath::Abs((fit1->Eval(timeOld) - fit1->Eval(untime_blm))/(untime_blm - timeOld));

        gr_3->SetPoint(gr_3_iter,angle_temp,blm/gradInt);
        gr_3->SetPointError(gr_3_iter,0.1,TMath::Sqrt(blm)/gradInt);
        gr_3_iter++;

        gr_2->SetPoint(gr_2_iter,untime_blm,gradInt);
        gr_2->SetPointError(gr_2_iter,0.001,0.0);
        gr_2_iter++;
    }
    Double_t x_min, y_min, x_max, y_max;
    gr_3->GetPoint(0,x_min,y_min);
    gr_3->GetPoint(gr_3->GetN()-1,x_max,y_max);
    cout<<endl<<x_min<<" <-> "<<x_max<<endl;
    TH1D* hgr_3 = new TH1D("hgr_3","BlmVsAngle",gr_3->GetN(),x_min,x_max);
    for(Int_t i = 0; i < gr_3->GetN(); i++)
    {
        Double_t x_temp, y_temp;
        gr_3->GetPoint(i,x_temp,y_temp);
        hgr_3->SetBinContent(i,y_temp);
        hgr_3->SetBinError(i,gr_3->GetErrorY(i));
    }
    cout<<endl;

    TH2D* hh_2_2   = new TH2D("hh_2_2","RP1 Internal Y vs Angle",maxUnixTime_run-minUnixTime_run,fitPol->Eval(minUnixTime_run),fitPol->Eval(maxUnixTime_run),N_PIXELS,0,N_PIXELS*PIXEL_SIZE);

    // TIMEPIX

    for(Int_t binx = 1; binx <= h_2->GetNbinsX(); binx++)
    {
        for(Int_t biny = 1; biny <= N_PIXELS; biny++)
        {
            hh_2->SetBinContent(binx,N_PIXELS-biny+1,h_2->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_2->SetBinError(binx,N_PIXELS-biny+1,TMath::Sqrt( TMath::Power(h_2->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
                                                                TMath::Power(h_2->GetBinContent(binx,biny)*h_22->GetBinError(binx)/(h_22->GetBinContent(binx)*h_22->GetBinContent(binx)),2) ));

            hh_2_2->SetBinContent(binx,N_PIXELS-biny+1,h_2->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_2_2->SetBinError(binx,N_PIXELS-biny+1,TMath::Sqrt( TMath::Power(h_2->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
                                                                TMath::Power(h_2->GetBinContent(binx,biny)*h_22->GetBinError(binx)/(h_22->GetBinContent(binx)*h_22->GetBinContent(binx)),2) ));
        }
    }
    for(Int_t binx = 1; binx <= N_PIXELS; binx++)
    {
        for(Int_t biny = 1; biny <= N_PIXELS; biny++)
        {
            hh_1->SetBinContent(N_PIXELS-biny+1,binx,h_1->GetBinContent(binx,biny));
            hh_1->SetBinError(N_PIXELS-biny+1,binx,h_1->GetBinError(binx,biny));
        }
    }
    cout<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name.Data(),"recreate");

    gr_0->Write();
    gr_1->Write();
    gr_2->Write();
    gr_3->Write();
    hgr_3->Write();

    cc->Write();

    hh_1->Write();
    hh_2->Write();
    hh_2_2->Write();
    hh_10->Write();

    cFit->Write();
    fitPol->Write();

    file->Write();
    //--------------------------------------------------------------------------//
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
