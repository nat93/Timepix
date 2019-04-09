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
void function_9();

Double_t fitexpo(Double_t *x, Double_t *par);

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout<<endl;
        cout<<">>> INPUT PARAMETERS <<<"<<endl;
        cout<<"[0] - script name"<<endl;
        cout<<"[1] - function number:"<<endl;
        cout<<"         1 --> function_1() - CRY2/4 angular scan Single-Channeling (Timepix RP1I ToA)"<<endl;
        cout<<"         2 --> function_2() - Convert crystal position/orientation or blm counts"<<endl;
        cout<<"         3 --> function_3() - CRY2 angular scan Single-Channeling (BLM vs Angle, Timepix)"<<endl;
        cout<<"         4 --> function_4() - CRY4 linear scan Double-Channeling (Timepix RP1I ToA)"<<endl;
        cout<<"         5 --> function_5() - CRY4 linear scan Double-Channeling (Timepix & Position)"<<endl;
        cout<<"         6 --> function_6() - CRY4 angular scan Double-Channeling (Timepix & Angle)"<<endl;
        cout<<"         7 --> function_7() - like function_3() without Timepix"<<endl;
        cout<<"         8 --> function_8() - Convert BCTDC counts"<<endl;
        cout<<"         9 --> function_9() - like function_7() with BCTDC"<<endl;
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
        case 5:
            function_5();
            break;
        case 6:
            function_6();
            break;
        case 7:
            function_7();
            break;
        case 8:
            function_8();
            break;
        case 9:
            function_9();
            break;
        default:
            break;
        }
    }
    return 0;
}

void function_1()
{
    cout<<endl<<">>> function_1() <<<"<<endl;
    //------------------------------------------------------------------------------//
    // CRYSTAL2 Angular Scan
    // 2017_10_17
    Long64_t minUnixTime_run    = 1508258897;
    Long64_t maxUnixTime_run    = 1508259288;
    Long64_t entryINI           = 271700;
    //------------------------------------------------------------------------------//
    // CRYSTAL4 Angular Scan
    // 2017_10_17
//    Long64_t minUnixTime_run    = 1508261894;
//    Long64_t maxUnixTime_run    = 1508262357;
//    Long64_t entryINI           = 295500;
    //------------------------------------------------------------------------------//
    // CRYSTAL2 in CH, step1
    // 2017_10_17
//    Long64_t minUnixTime_run    = 1508282585;
//    Long64_t maxUnixTime_run    = 1508282948;
//    Long64_t entryINI           = 438300;
    //------------------------------------------------------------------------------//
    // CRYSTAL2 in CH, step2
    // 2017_10_17
//    Long64_t minUnixTime_run    = 1508287980;
//    Long64_t maxUnixTime_run    = 1508288100;
//    Long64_t entryINI           = 488100;
    //------------------------------------------------------------------------------//
    // CRYSTAL2 in CH, step4
    // 2017_10_17
//    Long64_t minUnixTime_run    = 1508288322;
//    Long64_t maxUnixTime_run    = 1508288610;
//    Long64_t entryINI           = 491500;
    //------------------------------------------------------------------------------//
    // CRYSTAL2 in CH, step4
    // 2017_10_17
    // pos1
//    Long64_t minUnixTime_run    = 1508285169;
//    Long64_t maxUnixTime_run    = 1508285399;
//    Long64_t entryINI           = 462327;
    // pos2
//    Long64_t minUnixTime_run    = 1508286022;
//    Long64_t maxUnixTime_run    = 1508286262;
//    Long64_t entryINI           = 462327;
    // pos3
//    Long64_t minUnixTime_run    = 1508286830;
//    Long64_t maxUnixTime_run    = 1508287061;
//    Long64_t entryINI           = 462327;
    // pos4
//    Long64_t minUnixTime_run    = 1508287609;
//    Long64_t maxUnixTime_run    = 1508287772;
//    Long64_t entryINI           = 477794;
    // pos5
//    Long64_t minUnixTime_run    = 1508288322;
//    Long64_t maxUnixTime_run    = 1508288610;
//    Long64_t entryINI           = 462327;
    // pos6
//    Long64_t minUnixTime_run    = 1508289151;
//    Long64_t maxUnixTime_run    = 1508289339;
//    Long64_t entryINI           = 491537;
    // pos7
//    Long64_t minUnixTime_run    = 1508289965;
//    Long64_t maxUnixTime_run    = 1508290261;
//    Long64_t entryINI           = 498358;
    // pos8
//    Long64_t minUnixTime_run    = 1508291117;
//    Long64_t maxUnixTime_run    = 1508291321;
//    Long64_t entryINI           = 491537;
    // pos9
//    Long64_t minUnixTime_run    = 1508292393;
//    Long64_t maxUnixTime_run    = 1508292582;
//    Long64_t entryINI           = 505482;


    Bool_t entryINIstatus       = kFALSE;
    //------------------------------------------------------------------------------//
    // READ TREE
    //------------------------------------------------------------------------------//
    Int_t chipID = 2;// 2 - RP1I 2017
    TChain* fChain;
    Long64_t nEntries;

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    UInt_t _AcquisType = -1, _TrigType = -1;
    Long64_t _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS*2][N_PIXELS*2];

    TString _fileName   = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2017_10_17_RUN_1.root";

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
    cout<<"--> ChipID: "<<chipID<<" (2 - RP1I 2017)"<<endl;
    cout<<"--> InputFileName: "<<_fileName<<endl;
    cout<<"--> Number of Frames: "<<nEntries<<endl;

    //------------------------------------------------------------------------------//
    // CREATE HISTOGRAMS
    //------------------------------------------------------------------------------//

    fChain->GetEntry(0);
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    TH2D* h_1   = new TH2D("h_1","RP1 Internal",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_2   = new TH2D("h_2","RP1 Internal Y vs Time",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run,N_PIXELS,0,N_PIXELS);
    TH2D* h_3   = new TH2D("h_3","RP1 Internal X vs Time",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run,N_PIXELS,0,N_PIXELS);
    TH2D* h_4   = new TH2D("h_4","RP1 Internal X vs Time (cut)",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run,N_PIXELS,0,N_PIXELS);
    TH2D* h_5   = new TH2D("h_5","RP1 Internal YX [mm]",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* h_1_pl1   = new TH2D("h_1_pl1","RP1 Internal",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_1_pl2   = new TH2D("h_1_pl2","RP1 Internal",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_1_pl3   = new TH2D("h_1_pl3","RP1 Internal",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_1_pl4   = new TH2D("h_1_pl4","RP1 Internal",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

    TH1D* h_22  = new TH1D("h_22","RP1 Internal Time (Frames)",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run);

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

        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0)
                {
                    h_1->Fill(xi,yi,1);
                    h_2->Fill(frameTime,yi,1);
                    h_3->Fill(frameTime,xi,1);

                    if(yi < 160 && yi > 120)
                    {
                        h_4->Fill(frameTime,xi,1);
                    }

                    // For Angular Scan with Crystal2, close to the axial channeling
                    //plane 1
                    if(frameTime > 1508258960 && frameTime < 1508259000)
                    {
                        h_1_pl1->Fill(xi,yi,1.0/(1508259000 - 1508258960));
                        h_1_pl4->Fill(xi,yi,1.0/(1508259000 - 1508258960));
                    }
                    //plane 2
                    if(frameTime > 1508259060 && frameTime < 1508259140)
                    {
                        h_1_pl2->Fill(xi,yi,1.0/(1508259140 - 1508259060));
                        h_1_pl4->Fill(xi,yi,1.0/(1508259140 - 1508259060));
                    }
                    //plane 3
                    if(frameTime > 1508259250 && frameTime < 1508259280)
                    {
                        h_1_pl3->Fill(xi,yi,1.0/(1508259280 - 1508259250));
                        h_1_pl4->Fill(xi,yi,1.0/(1508259280 - 1508259250));
                    }
                }
            }
        }

        h_22->Fill(frameTime,1);

        nFrames++;
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

    TFile *_file = new TFile("output_function_1.root","RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File '"<<_file->GetName()<<"' opened successfully!"<<endl;

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_5->Write();
    h_1_pl1->Write();
    h_1_pl2->Write();
    h_1_pl3->Write();
    h_1_pl4->Write();
    h_22->Write();

    _file->Close();
}

void function_2()
{
    TString Inputfilename   = "collimator_linearscan_2017_10_17_pos72.7.dat";
    TString Outputfilename   = "collimator_linearscan_2017_10_17_pos72.7.root";

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
//        untime *= 1e6;

        tree->Fill();
    }

    file->Write();
    file->Close();

    cout<<"--> output file: "<<Outputfilename<<endl;
}

void function_3()
{
    //=============================================================//
    // CRYSTAL2 Angular Scan
    // 2017_10_17
    TString output_file_name    = "output_function_3.root";
    TString input_file_tpx      = "output_function_1.root";
    TString input_file_motor    = "crystal2_angularscan_2017_10_17.root";
    TString input_file_blm      = "blm_crystal2_angularscan_2017_10_17.root";
    Double_t blm_norm           = 1200;
    Double_t chorientation      = 2450.0; // [urad]
    Int_t angLimMax             = chorientation-2000;
    Int_t angLimMin             = chorientation-3000;
    Int_t angBins               = 1000;
    Long64_t minUnixTime_run    = 1508258897;
    Long64_t maxUnixTime_run    = 1508259288;
    //=============================================================//

    // Timepix
    TFile* _fileTpx = TFile::Open(input_file_tpx.Data());
    TH2D* h_1 = (TH2D*)_fileTpx->Get("h_1");
    TH2D* h_2 = (TH2D*)_fileTpx->Get("h_2");
    TH2D* h_3 = (TH2D*)_fileTpx->Get("h_3");
    TH2D* h_1_pl1 = (TH2D*)_fileTpx->Get("h_1_pl1");
    TH2D* h_1_pl2 = (TH2D*)_fileTpx->Get("h_1_pl2");
    TH2D* h_1_pl3 = (TH2D*)_fileTpx->Get("h_1_pl3");
    TH2D* h_1_pl4 = (TH2D*)_fileTpx->Get("h_1_pl4");

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

    fChain3->SetBranchAddress("Value",  &angle);
    fChain3->SetBranchAddress("UnixTime",       &untime_motor);

    cout<<"--> Input file with measurements: "<<input_file_tpx<<endl;

    cout<<"--> Input file with blm: "<<input_file_blm<<endl;
    Double_t nEntries_2 = fChain2->GetEntries();
    cout<<"--> nEntries: "<<nEntries_2<<endl;

    cout<<"--> Input file with motor: "<<input_file_motor<<endl;
    Double_t nEntries_3 = fChain3->GetEntries();
    cout<<"--> nEntries: "<<nEntries_3<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//

    TH2D* hh_1   = new TH2D("hh_1","RP1 Internal",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_2   = new TH2D("hh_2","RP1 Internal Y vs Time",maxUnixTime_run-minUnixTime_run+9,0,maxUnixTime_run-minUnixTime_run+9,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_3   = new TH2D("hh_3","RP1 Internal X vs Time",maxUnixTime_run-minUnixTime_run+9,0,maxUnixTime_run-minUnixTime_run+9,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_1_pl1   = new TH2D("hh_1_pl1","RP1 Internal",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_1_pl2   = new TH2D("hh_1_pl2","RP1 Internal",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_1_pl3   = new TH2D("hh_1_pl3","RP1 Internal",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_1_pl4   = new TH2D("hh_1_pl4","RP1 Internal",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);

    TH2D* hh_10 = new TH2D("hh_10","CrystalAngle vs UnixTime",(maxUnixTime_run-minUnixTime_run),minUnixTime_run,maxUnixTime_run,angBins,angLimMin,angLimMax);
    TH2D* hh_14 = new TH2D("hh_14","BLM vs Angle",angBins,angLimMin,angLimMax,3000,-1,2);
    TProfile* hprof_4  = new TProfile("hprof_4","Profile of BLM vs Angle",angBins/5,angLimMin,angLimMax,-1,2);

    TGraphErrors* gr_0 = new TGraphErrors();
    gr_0->SetName("gr_0");
    gr_0->SetTitle("Angle");
    Int_t gr_0_iter = 0;

    //--------------------------------------------------------------------------//

    // MOTOR    
    for(Int_t eventID = 0; eventID < nEntries_3; eventID++)
    {
        fChain3->GetEntry(eventID);

        angle = chorientation-TMath::Abs(angle);

        hh_10->Fill(untime_motor,angle);
        gr_0->SetPoint(gr_0_iter,untime_motor,angle);
        gr_0->SetPointError(gr_0_iter,0.001,0.1);
        gr_0_iter++;
    }
    cout<<endl;

    TCanvas* cFit = new TCanvas("cFit","cFit");
    cFit->cd();
    gr_0->Draw();
    fChain3->GetEntry(0);
    Double_t ttMin = untime_motor;
    fChain3->GetEntry(nEntries_3-1);
    Double_t ttMax = untime_motor;
    TF1* fitPol = new TF1("fitPol","pol1",ttMin,ttMax);
    gr_0->Fit(fitPol,"R+");

    // TIMEPIX
    for(Int_t xi = 1; xi <= N_PIXELS; xi++)
    {
        for(Int_t yi = 1; yi <= N_PIXELS; yi++)
        {
            hh_1->SetBinContent(N_PIXELS-yi+1,xi,h_1->GetBinContent(xi,yi));
            hh_1_pl1->SetBinContent(N_PIXELS-yi+1,xi,h_1_pl1->GetBinContent(xi,yi));
            hh_1_pl2->SetBinContent(N_PIXELS-yi+1,xi,h_1_pl2->GetBinContent(xi,yi));
            hh_1_pl3->SetBinContent(N_PIXELS-yi+1,xi,h_1_pl3->GetBinContent(xi,yi));
            hh_1_pl4->SetBinContent(N_PIXELS-yi+1,xi,h_1_pl4->GetBinContent(xi,yi));
        }
    }

    for(Int_t binx = 1; binx <= h_2->GetNbinsX(); binx++)
    {
        for(Int_t biny = 1; biny <= N_PIXELS; biny++)
        {
            hh_2->SetBinContent(binx,N_PIXELS-biny+1,h_2->GetBinContent(binx,biny));
            hh_3->SetBinContent(binx,biny,h_3->GetBinContent(binx,biny));
        }
    }
    cout<<endl;

    // BLM
    for(Int_t eventID = 0; eventID < nEntries_2; eventID++)
    {
        fChain2->GetEntry(eventID);

        Double_t angle_temp = fitPol->Eval(untime_blm);

        hh_14->Fill(angle_temp,blm/blm_norm);

        hprof_4->Fill(angle_temp,blm/blm_norm);
    }
    cout<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name.Data(),"recreate");

    hh_1->Write();
    hh_2->Write();
    hh_3->Write();
    hh_1_pl1->Write();
    hh_1_pl2->Write();
    hh_1_pl3->Write();
    hh_1_pl4->Write();
    hh_10->Write();
    hh_14->Write();
    hprof_4->Write();
    gr_0->Write();
    fitPol->Write();
    cFit->Write();

    file->Write();
    //--------------------------------------------------------------------------//
}

void function_4()
{
    cout<<endl<<">>> function_4() <<<"<<endl;
    //------------------------------------------------------------------------------//
    // CRYSTAL4 Linear Scan
    // 2017_10_17
    Long64_t minUnixTime_run    = 1508260568;
    Long64_t maxUnixTime_run    = 1508261376;
    Long64_t entryINI           = 271700;
    Bool_t entryINIstatus       = kFALSE;

    //------------------------------------------------------------------------------//
    // READ TREE
    //------------------------------------------------------------------------------//
    Int_t chipID = 2;// 2 - RP1I 2017
    TChain* fChain;
    Long64_t nEntries;

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    UInt_t _AcquisType = -1, _TrigType = -1;
    Long64_t _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS*2][N_PIXELS*2];

    TString _fileName   = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2017_10_17_RUN_1.root";

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
    cout<<"--> ChipID: "<<chipID<<" (2 - RP1I 2017)"<<endl;
    cout<<"--> InputFileName: "<<_fileName<<endl;
    cout<<"--> Number of Frames: "<<nEntries<<endl;

    //------------------------------------------------------------------------------//
    // CREATE HISTOGRAMS
    //------------------------------------------------------------------------------//

    fChain->GetEntry(0);
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    TH2D* h_1   = new TH2D("h_1","RP1 Internal",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_2   = new TH2D("h_2","RP1 Internal Y vs Time",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run,N_PIXELS,0,N_PIXELS);
    TH2D* h_3   = new TH2D("h_3","RP1 Internal X vs Time",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run,N_PIXELS,0,N_PIXELS);
    TH1D* h_22  = new TH1D("h_22","RP1 Internal Time (Frames)",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run);
    TH1D* h_23  = new TH1D("h_23","RP1 Internal Time (Pixels)",maxUnixTime_run-minUnixTime_run,minUnixTime_run,maxUnixTime_run);

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

        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0)
                {
                    h_1->Fill(xi,yi,1);
                    h_2->Fill(frameTime,yi,1);
                    h_3->Fill(frameTime,xi,1);

                    h_23->Fill(frameTime,1);
                }
            }
        }

        h_22->Fill(frameTime,1);

        nFrames++;
    }
    cout<<endl<<"--> nFrames = "<<nFrames<<endl;
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    //------------------------------------------------------------------------------//
    // WRITE
    //------------------------------------------------------------------------------//

    TFile *_file = new TFile("output_function_4.root","RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File '"<<_file->GetName()<<"' opened successfully!"<<endl;

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_22->Write();
    h_23->Write();

    _file->Close();
}

void function_5()
{
    //=============================================================//
    // CRYSTAL4 Linear Scan
    // 2017_10_17
    TString output_file_name    = "output_function_5.root";
    TString input_file_tpx      = "output_function_4.root";
    TString input_file_motor    = "crystal4_linearscan_2017_10_17.root";
    Long64_t minUnixTime_run    = 1508260568;
    Long64_t maxUnixTime_run    = 1508261376;
    Int_t posBins               = 2000;
    Double_t posLimMin          = 60;
    Double_t posLimMax          = 80;

    //=============================================================//

    // Timepix
    TFile* _fileTpx = TFile::Open(input_file_tpx.Data());
    TH2D* h_1 = (TH2D*)_fileTpx->Get("h_1");
    TH2D* h_2 = (TH2D*)_fileTpx->Get("h_2");
    TH1D* h_22 = (TH1D*)_fileTpx->Get("h_22");

    // Motor
    Double_t position;
    Double_t untime_motor;

    TChain *fChain3 = new TChain("Tree");
    fChain3->Add(input_file_motor);

    fChain3->SetBranchAddress("Value",          &position);
    fChain3->SetBranchAddress("UnixTime",       &untime_motor);

    cout<<"--> Input file with measurements: "<<input_file_tpx<<endl;

    cout<<"--> Input file with motor: "<<input_file_motor<<endl;
    Double_t nEntries_3 = fChain3->GetEntries();
    cout<<"--> nEntries: "<<nEntries_3<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//

    TH2D* hh_1   = new TH2D("hh_1","RP1 Internal YX [mm]",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_2   = new TH2D("hh_2","RP1 Internal Y vs Time",maxUnixTime_run-minUnixTime_run,0,maxUnixTime_run-minUnixTime_run,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_10 = new TH2D("hh_10","CrystalAngle vs UnixTime",(maxUnixTime_run-minUnixTime_run),minUnixTime_run,maxUnixTime_run,posBins,posLimMin,posLimMax);

    TGraphErrors* gr_0 = new TGraphErrors();
    gr_0->SetName("gr_0");
    gr_0->SetTitle("Position");
    Int_t gr_0_iter = 0;

    //--------------------------------------------------------------------------//

    // MOTOR
    for(Int_t eventID = 0; eventID < nEntries_3; eventID++)
    {
        fChain3->GetEntry(eventID);

        position = TMath::Abs(position);

        hh_10->Fill(untime_motor,position);
        gr_0->SetPoint(gr_0_iter,untime_motor,position);
        gr_0->SetPointError(gr_0_iter,0.001,0.001);
        gr_0_iter++;
    }
    cout<<endl;

    TH2D* hh_3   = new TH2D("hh_3","RP1 Internal Y vs Position",maxUnixTime_run-minUnixTime_run,gr_0->Eval(minUnixTime_run),gr_0->Eval(maxUnixTime_run),N_PIXELS,0,N_PIXELS*PIXEL_SIZE);

    // TIMEPIX

    for(Int_t binx = 1; binx <= h_2->GetNbinsX(); binx++)
    {
        for(Int_t biny = 1; biny <= N_PIXELS; biny++)
        {
            hh_2->SetBinContent(binx,N_PIXELS-biny+1,h_2->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_2->SetBinError(binx,N_PIXELS-biny+1,TMath::Sqrt( TMath::Power(h_2->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
                                                                TMath::Power(h_2->GetBinContent(binx,biny)*h_22->GetBinError(binx)/(h_22->GetBinContent(binx)*h_22->GetBinContent(binx)),2) ));

            hh_3->SetBinContent(binx,N_PIXELS-biny+1,h_2->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_3->SetBinError(binx,N_PIXELS-biny+1,TMath::Sqrt( TMath::Power(h_2->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
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

    hh_1->Write();
    hh_2->Write();
    hh_3->Write();
    hh_10->Write();
    gr_0->Write();

    file->Write();
    //--------------------------------------------------------------------------//
}

void function_6()
{
    //=============================================================//
    // CRYSTAL4 Angular Scan
    // 2017_10_17
    TString output_file_name    = "output_function_6.root";
    TString input_file_tpx      = "output_function_1.root";
    TString input_file_motor    = "crystal4_angularscan_2017_10_17.root";
    Long64_t minUnixTime_run    = 1508261894;
    Long64_t maxUnixTime_run    = 1508262357;
    Int_t angBins               = 2000;
    Double_t angLimMin          =  500;
    Double_t angLimMax          = 1500;

    //=============================================================//

    // Timepix
    TFile* _fileTpx = TFile::Open(input_file_tpx.Data());
    TH2D* h_1 = (TH2D*)_fileTpx->Get("h_1");
    TH2D* h_2 = (TH2D*)_fileTpx->Get("h_2");
    TH2D* h_4 = (TH2D*)_fileTpx->Get("h_4");
    TH1D* h_22 = (TH1D*)_fileTpx->Get("h_22");

    // Motor
    Double_t angle;
    Double_t untime_motor;

    TChain *fChain3 = new TChain("Tree");
    fChain3->Add(input_file_motor);

    fChain3->SetBranchAddress("Value",          &angle);
    fChain3->SetBranchAddress("UnixTime",       &untime_motor);

    cout<<"--> Input file with measurements: "<<input_file_tpx<<endl;

    cout<<"--> Input file with motor: "<<input_file_motor<<endl;
    Double_t nEntries_3 = fChain3->GetEntries();
    cout<<"--> nEntries: "<<nEntries_3<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//

    TH2D* hh_1   = new TH2D("hh_1","RP1 Internal YX [mm]",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_2   = new TH2D("hh_2","RP1 Internal Y vs Time",maxUnixTime_run-minUnixTime_run,0,maxUnixTime_run-minUnixTime_run,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_3   = new TH2D("hh_3","RP1 Internal X vs Time",maxUnixTime_run-minUnixTime_run,0,maxUnixTime_run-minUnixTime_run,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_10 = new TH2D("hh_10","CrystalAngle vs UnixTime",(maxUnixTime_run-minUnixTime_run),minUnixTime_run,maxUnixTime_run,angBins,angLimMin,angLimMax);

    TGraphErrors* gr_0 = new TGraphErrors();
    gr_0->SetName("gr_0");
    gr_0->SetTitle("Angle");
    Int_t gr_0_iter = 0;

    //--------------------------------------------------------------------------//

    // MOTOR
    for(Int_t eventID = 0; eventID < nEntries_3; eventID++)
    {
        fChain3->GetEntry(eventID);

        angle = TMath::Abs(angle);

        hh_10->Fill(untime_motor,angle);
        gr_0->SetPoint(gr_0_iter,untime_motor,angle);
        gr_0->SetPointError(gr_0_iter,0.001,0.1);
        gr_0_iter++;
    }
    cout<<endl;

    TH2D* hh_2_2   = new TH2D("hh_2_2","RP1 Internal Y vs Angle",maxUnixTime_run-minUnixTime_run,gr_0->Eval(minUnixTime_run),gr_0->Eval(maxUnixTime_run),N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* hh_3_2   = new TH2D("hh_3_2","RP1 Internal X vs Angle",maxUnixTime_run-minUnixTime_run,gr_0->Eval(minUnixTime_run),gr_0->Eval(maxUnixTime_run),N_PIXELS,0,N_PIXELS*PIXEL_SIZE);

    // TIMEPIX

    for(Int_t binx = 1; binx <= h_2->GetNbinsX(); binx++)
    {
        for(Int_t biny = 1; biny <= N_PIXELS; biny++)
        {
            hh_2->SetBinContent(binx,N_PIXELS-biny+1,h_2->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_2->SetBinError(binx,N_PIXELS-biny+1,TMath::Sqrt( TMath::Power(h_2->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
                                                                TMath::Power(h_2->GetBinContent(binx,biny)*h_22->GetBinError(binx)/(h_22->GetBinContent(binx)*h_22->GetBinContent(binx)),2) ));

            hh_3->SetBinContent(binx,biny,h_4->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_3->SetBinError(binx,biny,TMath::Sqrt( TMath::Power(h_4->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
                                                                TMath::Power(h_4->GetBinContent(binx,biny)*h_22->GetBinError(binx)/(h_22->GetBinContent(binx)*h_22->GetBinContent(binx)),2) ));

            hh_2_2->SetBinContent(binx,N_PIXELS-biny+1,h_2->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_2_2->SetBinError(binx,N_PIXELS-biny+1,TMath::Sqrt( TMath::Power(h_2->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
                                                                TMath::Power(h_2->GetBinContent(binx,biny)*h_22->GetBinError(binx)/(h_22->GetBinContent(binx)*h_22->GetBinContent(binx)),2) ));

            hh_3_2->SetBinContent(binx,biny,h_4->GetBinContent(binx,biny)/(h_22->GetBinContent(binx)));
            hh_3_2->SetBinError(binx,biny,TMath::Sqrt( TMath::Power(h_4->GetBinError(binx,biny)/(h_22->GetBinContent(binx)),2) +
                                                                TMath::Power(h_4->GetBinContent(binx,biny)*h_22->GetBinError(binx)/(h_22->GetBinContent(binx)*h_22->GetBinContent(binx)),2) ));
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

    hh_1->Write();
    hh_2->Write();
    hh_2_2->Write();
    hh_3->Write();
    hh_3_2->Write();
    hh_10->Write();
    gr_0->Write();

    file->Write();
    //--------------------------------------------------------------------------//
}

void function_7()
{
    //=============================================================//
    // CRYSTAL4 Angular Scan
    // 2017_10_17
    TString output_file_name    = "output_function_7.root";
    TString input_file_motor    = "crystal4_angularscan_2017_10_17_pos72.7.root";
    TString input_file_blm      = "blm_crystal4_angularscan_2017_10_17_pos72.7.root";
    Int_t angLimMax             = 2000;
    Int_t angLimMin             = 0;
    Int_t angBins               = 2000;
    //=============================================================//

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

    fChain3->SetBranchAddress("Value",  &angle);
    fChain3->SetBranchAddress("UnixTime",       &untime_motor);

    cout<<"--> Input file with blm: "<<input_file_blm<<endl;
    Double_t nEntries_2 = fChain2->GetEntries();
    cout<<"--> nEntries: "<<nEntries_2<<endl;

    cout<<"--> Input file with motor: "<<input_file_motor<<endl;
    Double_t nEntries_3 = fChain3->GetEntries();
    cout<<"--> nEntries: "<<nEntries_3<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- HISTOS ----------------------------------//
    //--------------------------------------------------------------------------//

    TH2D* hh_14 = new TH2D("hh_14","BLM vs Angle",angBins,angLimMin,angLimMax,1000,0,1000);
    TProfile* hprof_4  = new TProfile("hprof_4","Profile of BLM vs Angle",angBins/5,angLimMin,angLimMax,0,1000);

    TGraphErrors* gr_0 = new TGraphErrors();
    gr_0->SetName("gr_0");
    gr_0->SetTitle("Angle");
    Int_t gr_0_iter = 0;

    //--------------------------------------------------------------------------//

    // MOTOR
    for(Int_t eventID = 0; eventID < nEntries_3; eventID++)
    {
        fChain3->GetEntry(eventID);

        angle = TMath::Abs(angle);

        gr_0->SetPoint(gr_0_iter,untime_motor,angle);
        gr_0->SetPointError(gr_0_iter,0.001,0.1);
        gr_0_iter++;
    }
    cout<<endl;

    // BLM
    for(Int_t eventID = 0; eventID < nEntries_2; eventID++)
    {
        fChain2->GetEntry(eventID);

        Double_t angle_temp = gr_0->Eval(untime_blm);

        hh_14->Fill(angle_temp,blm);

        hprof_4->Fill(angle_temp,blm);
    }
    cout<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name.Data(),"recreate");

    hh_14->Write();
    hprof_4->Write();
    gr_0->Write();

    file->Write();
    //--------------------------------------------------------------------------//
}

void function_8()
{
    TString Inputfilename   = "bctdc_collimator_linearscan_2017_10_17_pos72.7.dat";
    TString Outputfilename   = "bctdc_collimator_linearscan_2017_10_17_pos72.7.root";

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
//        untime *= 1e6;

        tree->Fill();
    }

    file->Write();
    file->Close();

    cout<<"--> output file: "<<Outputfilename<<endl;
}

void function_9()
{
    //=============================================================//
    // 2017_10_17
    TString output_file_name    = "output_function_9.root";
    TString input_file_motor    = "collimator_linearscan_2017_10_17.root";
    TString input_file_blm      = "blm_collimator_linearscan_2017_10_17_pos72.7.root";
    TString input_file_bctdc    = "bctdc_collimator_linearscan_2017_10_17_pos72.7.root";
    Int_t posLimMax             = 0;
    Int_t posLimMin             = -20;
    Int_t posBins               = 2000;
    //=============================================================//

    // BLM
    Double_t blm;
    Double_t untime_blm;

    TChain *fChain2 = new TChain("Tree");
    fChain2->Add(input_file_blm);

    fChain2->SetBranchAddress("Value",  &blm);
    fChain2->SetBranchAddress("UnixTime",       &untime_blm);

    // Motor
    Double_t position;
    Double_t untime_motor;

    TChain *fChain3 = new TChain("Tree");
    fChain3->Add(input_file_motor);

    fChain3->SetBranchAddress("Value",  &position);
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

    TH2D* hh_14 = new TH2D("hh_14","BLM vs Angle",posBins,posLimMin,posLimMax,1000,0,1000);
    TProfile* hprof_4  = new TProfile("hprof_4","Profile of BLM vs Angle",posBins/5,posLimMin,posLimMax,0,1000);

    TGraphErrors* gr_0 = new TGraphErrors();
    gr_0->SetName("gr_0");
    gr_0->SetTitle("Position");
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
    gr_3->SetTitle("BlmVsPosition");
    Int_t gr_3_iter = 0;

    //--------------------------------------------------------------------------//

    // MOTOR
    for(Int_t eventID = 0; eventID < nEntries_3; eventID++)
    {
        fChain3->GetEntry(eventID);

        gr_0->SetPoint(gr_0_iter,untime_motor,position);
        gr_0->SetPointError(gr_0_iter,0.001,0.001);
        gr_0_iter++;
    }
    cout<<endl;

    // BCTDC
    for(Int_t eventID = 0; eventID < nEntries_4; eventID++)
    {
        fChain4->GetEntry(eventID);

        gr_1->SetPoint(gr_1_iter,untime_bctdc,bctdc);
        gr_1->SetPointError(gr_1_iter,0.001,bctdc_std);
        gr_1_iter++;
    }
    TCanvas* cc = new TCanvas("cc","cc");
    cc->cd();
    gr_1->Draw("APL");
    fChain4->GetEntry(0);
    Double_t timeINI = untime_bctdc;
    fChain4->GetEntry(nEntries_4-1);
    Double_t timeFIN = untime_bctdc;

    TF1* fit1 = new TF1("fit1","expo",timeINI,timeFIN);
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

        Double_t pos_temp = gr_0->Eval(untime_blm);
        Double_t gradInt = TMath::Abs((fit1->Eval(timeOld) - fit1->Eval(untime_blm))/(untime_blm - timeOld));

        gr_3->SetPoint(gr_3_iter,pos_temp,blm/gradInt);
        gr_3->SetPointError(gr_3_iter,0.001,TMath::Sqrt(blm)/gradInt);
        gr_3_iter++;

        gr_2->SetPoint(gr_2_iter,untime_blm,gradInt);
        gr_2->SetPointError(gr_2_iter,0.001,0.0);
        gr_2_iter++;
    }

    Double_t x_min, y_min, x_max, y_max;
    gr_3->GetPoint(0,x_min,y_min);
    gr_3->GetPoint(gr_3->GetN()-1,x_max,y_max);
    cout<<endl<<x_max<<endl;
    TH1D* hgr_3 = new TH1D("hgr_3","BlmVsPosition",gr_3->GetN(),x_min,x_max);
    for(Int_t i = 0; i < gr_3->GetN(); i++)
    {
        Double_t x_temp, y_temp;
        gr_3->GetPoint(i,x_temp,y_temp);
        hgr_3->SetBinContent(i,y_temp);
        hgr_3->SetBinError(i,gr_3->GetErrorY(i));
    }
    cout<<endl;

    //--------------------------------------------------------------------------//
    //-------------------------------- WRITE -----------------------------------//
    //--------------------------------------------------------------------------//
    cout<<"--> Output file: "<<output_file_name<<endl;
    TFile* file = new TFile(output_file_name.Data(),"recreate");

    hh_14->Write();
    hprof_4->Write();
    gr_0->Write();
    gr_1->Write();
    gr_2->Write();
    gr_3->Write();
    hgr_3->Write();
    cc->Write();

    file->Write();
    //--------------------------------------------------------------------------//
}

Double_t fitexpo(Double_t *x, Double_t *par)
{
    return par[0]*TMath::Exp(par[1]*x[0]);
}
/**/
