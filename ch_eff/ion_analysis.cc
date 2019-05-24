//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2019
// Description:     Analysis of the single chip. For H8 and SPS
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

const Int_t N_PIXELS                = 512;      // The maximum number of pixels per axis per chip
const Int_t N_MAX_CLUSTERS          = 10000;    // The maximum number of clusters per frame
const Int_t N_MAX_CLOCKS            = 11810;    // The maximum number of counts in the pixel
const Bool_t _cluster_analysis      = true;
const Double_t PIXEL_SIZE           = 0.055;    // Pixel size [mm]
const Double_t cluster_time_jitter  = 200.0;    // Time difference between fired pixels inside cluster [ns]

void epochtime2date(time_t in_time, string &out_time);
int clusteranalysis(Long64_t _matrix[][N_PIXELS], Int_t &cluster_num, Int_t *cluster_size, Double_t *cluster_pos_x, Double_t *cluster_pos_x_err, Double_t *cluster_pos_y, Double_t *cluster_pos_y_err, Double_t *cluster_clocks, Int_t AcqType, Double_t Clock);
void checkregionMPX(Long64_t _matrix[][N_PIXELS], Int_t **fired_matrix, Int_t xi, Int_t yi);
void checkregionTOA(Long64_t _matrix[][N_PIXELS], Int_t **fired_matrix, Int_t xi, Int_t yi, Double_t Clock);
bool checkcondition(Long64_t x, Long64_t ref_x, Double_t delta, Double_t Clock);

TH1D* h_com_1;
TH2D* h_com_2;
TH1D* h_com_3;
TH2D* h_com_4;

Bool_t oneCluster = kTRUE;

void function_1();
void function_2();

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout<<endl;
        cout<<">>> INPUT PARAMETERS <<<"<<endl;
        cout<<"[0] - script name"<<endl;
        cout<<"[1] - function number:"<<endl;
        cout<<"         1 --> function_1() - H8     (Initial time: Tue Nov 27 2018 13:21:57 | Final time: Wed Nov 28 2018 12:29:57)"<<endl;
        cout<<"         2 --> function_2() - SPS    (Initial time: Mon Nov 26 2018 17:47:37 | Final time: Tue Nov 27 2018 01:02:01)"<<endl;
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
        default:
            cout<<endl<<"### NOTHING TO DO ;) ###"<<endl<<endl;
            break;
        }
    }
    return 0;
}

void function_1()
{
    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_11_27_RUN_1.root";
    TString outFileName     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_11_27_RUN_1_outputfile_function_1.root";
    Int_t chipID            = 0;

    Long64_t minUnixTime_run    = 1543340483;
    Long64_t maxUnixTime_run    = 1543340484;
    Long64_t entryINI           = 35022;

    Bool_t entryINIstatus       = kFALSE;

    cout<<endl<<"--> function_1() <--"<<endl<<endl;

    time_t start_time, stop_time;
    start_time = time(NULL);


    // COMMON VARIABLES
    TChain* fChain;
    Long64_t nEntries;

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    UInt_t _AcquisType = -1, _TrigType = -1;
    Long64_t _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS][N_PIXELS];

    // For TIMEPIX

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
    cout<<"--> ChipID: "<<chipID<<endl;
    cout<<"--> InputFileName1: "<<inFileName1<<endl;
    cout<<"--> Number of nEntries: "<<nEntries<<endl;

    // For Output file
    cout<<"--> OutputFileName: "<<outFileName<<endl;
    TFile* outputfile = new TFile(outFileName.Data(),"RECREATE");
    if (outputfile->IsZombie())
    {
        cout<<"--> ERROR: The output file is zombie!"<<endl;
        assert(0);
    }

    //------------------------------------------------------------------------------//
    //------------------------ For output file with clusters infor -----------------//
    //------------------------------------------------------------------------------//
    Double_t unix_time_clinfo, pos_x_clinfo, pos_x_err_clinfo, pos_y_clinfo, pos_y_err_clinfo, clocks_clinfo;
    Int_t size_clinfo;
    Long64_t event_id_clinfo;

    TTree* tree = new TTree("Tree", "Cluster data");

    tree->Branch("UnixTime",        &unix_time_clinfo,      "unix_time_clinfo/D");
    tree->Branch("ClusterClocks",   &clocks_clinfo,         "clocks_clinfo/D");
    tree->Branch("ClusterSize",     &size_clinfo,           "size_clinfo/I");
    tree->Branch("ClusterPosX",     &pos_x_clinfo,          "pos_x_clinfo/D");
    tree->Branch("ClusterPosXerr",  &pos_x_err_clinfo,      "pos_x_err_clinfo/D");
    tree->Branch("ClusterPosY",     &pos_y_clinfo,          "pos_y_clinfo/D");
    tree->Branch("ClusterPosYerr",  &pos_y_err_clinfo,      "pos_y_err_clinfo/D");
    tree->Branch("EventID",         &event_id_clinfo,       "event_id_clinfo/L");
    tree->Branch("Clock",           &_Clock,                "_Clock/D");
    tree->Branch("Gate",            &_Gate,                 "_Gate/D");

    Int_t* cl_size          = new Int_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_x      = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_x_err  = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_y      = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_y_err  = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_clocks     = new Double_t[N_MAX_CLUSTERS];
    //------------------------------------------------------------------------------//
    //------------------------------- HISTOGRAMS -----------------------------------//
    //------------------------------------------------------------------------------//

    fChain->GetEntry(nEntries-1);
    Double_t UT_max = _Timems;

    fChain->GetEntry(0);
    Long64_t _event_ini = _event;
    Double_t UT_min = _Timems;
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    TH2D* h_1   = new TH2D("h_1","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH1D* h_2   = new TH1D("h_2","COUNTS vs X",N_PIXELS,0,N_PIXELS);
    TH1D* h_3   = new TH1D("h_3","COUNTS vs Y",N_PIXELS,0,N_PIXELS);
    TH2D* h_4   = new TH2D("h_4","X vs EventID",nEntries,0,nEntries,N_PIXELS,0,N_PIXELS);
    TH2D* h_5   = new TH2D("h_5","Y vs EventID",nEntries,0,nEntries,N_PIXELS,0,N_PIXELS);
    TH2D* h_6   = new TH2D("h_6","Y vs X vs C (in mm)",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* h_8   = new TH2D("h_8","Y vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH2D* h_9   = new TH2D("h_9","X vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH1D* h_10  = new TH1D("h_10","Hits vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0);
    TH1D* h_13  = new TH1D("h_13","Clusters number per frame",200,0,200);
    TH1D* h_14  = new TH1D("h_14","Hits vs EventID",nEntries,0,nEntries);
    TH1D* h_15  = new TH1D("h_15","Hits per frame",1e6,0,1e6);
    TH1D* h_16  = new TH1D("h_16","Cluster size",1e4,0,1e4);
    TH2D* h_20  = new TH2D("h_20","Y vs X vs C (inside the ch beam spot)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_21  = new TH2D("h_21","Y vs X vs C (inside the dch beam spot)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_22  = new TH2D("h_22","Y vs X vs C (outside the beam spot)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH1D* h_23  = new TH1D("h_23","#Delta Time between frames",10000,0,10000);
    TH1D* h_24  = new TH1D("h_24","_AcquisType vs EventID",nEntries,0,nEntries);
    TH1D* h_25  = new TH1D("h_25","_DeltaTHR vs EventID",nEntries,0,nEntries);
    TH1D* h_26  = new TH1D("h_26","_Gate vs EventID",nEntries,0,nEntries);
    TH1D* h_27  = new TH1D("h_27","_Bias vs EventID",nEntries,0,nEntries);
    TH2D* h_28  = new TH2D("h_28","Cluster Volume vs Cluster size",1e3,0,1e3,1e4,0,1e4);
    TH1D* h_29  = new TH1D("h_29","Number of Fired Pixels per frame",3e5,0,3e5);
    TH1D* h_31  = new TH1D("h_31","Cluster Volume",1e5,0,1e5);

    h_com_1     = new TH1D("h_com_1","#Delta Time inside cluster",2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_2     = new TH2D("h_com_2","#Delta Time inside cluster VS Cluster size",50,0,50,2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_3     = new TH1D("h_com_3","#Delta Time inside cluster (cluster size <= 2 pixels)",2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_4     = new TH2D("h_com_4","Y vs X vs ToT (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

    //------------------------------------------------------------------------------//
    Double_t zero_time = -999.999, event_time = -999.999, last_time = -999.999, previous_time = _Timems;
    Long64_t nFrames = 0;
    Long64_t frame_size = 0;

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------- MAIN LOOP ------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

//    nEntries /= 50;
    for(Long64_t i = entryINI; i < entryINI+1/*nEntries*/; i++)
    {
        if(i%1 == 0)
        {
            printf("\r--> Working: %3.3f %%",100*(Double_t)i/nEntries);
            fflush(stdout);
        }

//--------------------------------------------------------------//
//        To remove the first frame during the spill at H8
//        if(i > 0)
//        {
//            fChain->GetEntry(i);
//            Double_t t_now = _Timems;
//            fChain->GetEntry(i-1);
//            Double_t t_before = _Timems;
//            if(t_now - t_before > 5000) continue;
//        }
//--------------------------------------------------------------//

        fChain->GetEntry(i);

        //================================================//
        if(_Timems/1000.0 < minUnixTime_run) continue;
        if(_Timems/1000.0 > maxUnixTime_run) break;

        if(!entryINIstatus)
        {
            entryINI = i;
            cout<<endl<<"--> entryINI = "<<entryINI<<endl;
            entryINIstatus = kTRUE;
        }
        //================================================//

        if (_Timems > 0)
        {
            event_time  = _Timems;
        }
        else
        {
            continue;
        }

        if (zero_time < 0)
        {
            zero_time = _Timems;
        }
        last_time = _Timems;

        Int_t fired_pixels_num = 0;
        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0)
                {
                    fired_pixels_num++;
                }
            }
        }
        h_29->Fill(fired_pixels_num);

        h_24->Fill(_event-_event_ini,_AcquisType);
        h_25->Fill(_event-_event_ini,_DeltaTHR);
        h_26->Fill(_event-_event_ini,_Gate);
        h_27->Fill(_event-_event_ini,_Bias);

        frame_size = 0;
        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0)
                {
                    if(_AcquisType == 0 || _AcquisType == 1)
                    {
                        frame_size += _COUNTS[xi][yi];

                        h_1->Fill(xi,yi,_COUNTS[xi][yi]);
                        h_2->Fill(xi,_COUNTS[xi][yi]);
                        h_3->Fill(yi,_COUNTS[xi][yi]);
                        h_4->Fill(_event-_event_ini,xi,_COUNTS[xi][yi]);
                        h_5->Fill(_event-_event_ini,yi,_COUNTS[xi][yi]);
                        h_8->Fill(event_time/1000.0,yi,_COUNTS[xi][yi]);
                        h_9->Fill(event_time/1000.0,xi,_COUNTS[xi][yi]);
                    }
                    else
                    {
                        cout<<endl<<endl<<"--> ERROR:: _AcquisType is neither 0 nor 1"<<endl<<endl;
                        assert(0);
                    }
                }
            }
        }
        h_10->Fill(event_time/1000.0,frame_size);
        h_14->Fill(i,frame_size);
        h_15->Fill(frame_size);

        h_23->Fill(event_time - previous_time);
        previous_time = event_time;
        //---------------------------------------------------------------------------------------------------------------//
        //------------------------------------------ CLUSTER ANALYSIS ---------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        if(_cluster_analysis)
        {
            for(Int_t jk = 0; jk < N_MAX_CLUSTERS; jk++)
            {
                cl_size[jk]         = 0;
                cl_pos_x[jk]        = 0.0;
                cl_pos_x_err[jk]    = 0.0;
                cl_pos_y[jk]        = 0.0;
                cl_pos_y_err[jk]    = 0.0;
                cl_clocks[jk]       = 0.0;
            }

            Int_t cl_num = 0;

            if(clusteranalysis(_COUNTS,cl_num,cl_size,cl_pos_x,cl_pos_x_err,cl_pos_y,cl_pos_y_err,cl_clocks,_AcquisType,_Clock) != 0)
            {
                cout<<"--> ERROR: Cluster analysis did not check all pixels! ("<<clusteranalysis(_COUNTS,cl_num,cl_size,cl_pos_x,cl_pos_x_err,cl_pos_y,cl_pos_y_err,cl_clocks,_AcquisType,_Clock)<<")"<<endl;
                assert(0);
            }

            if(cl_num >= N_MAX_CLUSTERS)
            {
                cout<<endl<<"--> ERROR:: The maximum number of clusters per a frame was reached !!!"<<endl;
                assert(0);
            }

            for(Int_t yy = 0; yy < cl_num; yy++)
            {
                //------------------------------------------------------------------------------//
                //------------------------ For output file with clusters infor -----------------//
                //------------------------------------------------------------------------------//
                unix_time_clinfo    = event_time/1000.0;

                if(_AcquisType == 3) clocks_clinfo = (Double_t)N_MAX_CLOCKS - cl_clocks[yy];
                else clocks_clinfo = (Double_t)cl_clocks[yy];

                if(clocks_clinfo < 0) clocks_clinfo = 0;
                size_clinfo         = cl_size[yy];
                pos_x_clinfo        = cl_pos_x[yy];
                pos_x_err_clinfo    = cl_pos_x_err[yy];
                pos_y_clinfo        = cl_pos_y[yy];
                pos_y_err_clinfo    = cl_pos_y_err[yy];
                event_id_clinfo     = _event;
                if(_event != i) {event_id_clinfo = i;}

                tree->Fill();

                h_16->Fill(cl_size[yy]);
                h_28->Fill(cl_size[yy],cl_size[yy]*cl_clocks[yy]);
                h_31->Fill(cl_size[yy]*cl_clocks[yy]);

            }
            h_13->Fill(cl_num);
        }
        //---------------------------------------------------------------------------------------------------------------//
        nFrames++;
    }
    cout<<endl;

    delete [] cl_size;
    delete [] cl_pos_x;
    delete [] cl_pos_x_err;
    delete [] cl_pos_y;
    delete [] cl_pos_y_err;
    delete [] cl_clocks;

    cout<<endl;

    Double_t scfactor = _Gate*nFrames;
    cout<<"--> nFrames = "<<nFrames<<endl;
    cout<<"--> scaling factor = "<<scfactor<<" [sec]"<<endl;

    string zero_time_char;
    string last_time_char;
    epochtime2date(round(zero_time/1000),zero_time_char);
    epochtime2date(round(last_time/1000),last_time_char);

    cout<<endl;
    cout<<"--> Initial time: "<<zero_time_char<<" | Final time: "<<last_time_char<<endl;

    //------------------------------------------------------------------------------//
    //------------------------ For output file with clusters infor -----------------//
    //------------------------------------------------------------------------------//
    tree->Write();
    cout<<"--> The Tree with "<<tree->GetEntriesFast()<<" entries is written to the file: "<<outputfile->GetName()<<endl;
    //------------------------------------------------------------------------------//

    h_1->GetXaxis()->SetTitle("X [pixels]");
    h_1->GetYaxis()->SetTitle("Y [pixels]");
    h_2->GetXaxis()->SetTitle("X [pixels]");
    h_2->GetYaxis()->SetTitle("Counts");
    h_3->GetXaxis()->SetTitle("Y [pixels]");
    h_3->GetYaxis()->SetTitle("Counts");
    h_4->GetXaxis()->SetTitle("EventID");
    h_4->GetYaxis()->SetTitle("X [pixels]");
    h_5->GetXaxis()->SetTitle("EventID");
    h_5->GetYaxis()->SetTitle("Y [pixels]");
    h_6->GetXaxis()->SetTitle("X [mm]");
    h_6->GetYaxis()->SetTitle("Y [mm]");
    h_8->GetXaxis()->SetTitle("Time [sec]");
    h_8->GetYaxis()->SetTitle("Y [pixels]");
    h_9->GetXaxis()->SetTitle("Time [sec]");
    h_9->GetYaxis()->SetTitle("X [pixels]");
    h_10->GetXaxis()->SetTitle("Time [sec]");
    h_10->GetYaxis()->SetTitle("Hits");
    h_13->GetXaxis()->SetTitle("Clusters number per frame");
    h_13->GetYaxis()->SetTitle("Numer of frames");
    h_14->GetXaxis()->SetTitle("EventID");
    h_14->GetYaxis()->SetTitle("Hits");
    h_15->GetXaxis()->SetTitle("Hits per frame");
    h_15->GetYaxis()->SetTitle("Frames");
    h_16->GetXaxis()->SetTitle("Cluster size [pixels/cluster]");
    h_16->GetYaxis()->SetTitle("Number of clusters");
    h_20->GetXaxis()->SetTitle("X [pixels]");
    h_20->GetYaxis()->SetTitle("Y [pixels]");
    h_21->GetXaxis()->SetTitle("X [pixels]");
    h_21->GetYaxis()->SetTitle("Y [pixels]");
    h_22->GetXaxis()->SetTitle("X [pixels]");
    h_22->GetYaxis()->SetTitle("Y [pixels]");
    h_23->GetXaxis()->SetTitle("#Delta Time [ms]");
    h_23->GetYaxis()->SetTitle("Number of frames");
    h_24->GetXaxis()->SetTitle("EventID");
    h_24->GetYaxis()->SetTitle("_AcquisType");
    h_25->GetXaxis()->SetTitle("EventID");
    h_25->GetYaxis()->SetTitle("_DeltaTHR");
    h_26->GetXaxis()->SetTitle("EventID");
    h_26->GetYaxis()->SetTitle("_Gate");
    h_27->GetXaxis()->SetTitle("EventID");
    h_27->GetYaxis()->SetTitle("_Bias");
    h_28->GetXaxis()->SetTitle("Cluster size [pixels/cluster]");
    h_28->GetYaxis()->SetTitle("Cluster volume [counts/cluster]");
    h_31->GetYaxis()->SetTitle("Cluster volume [counts/cluster]");

    h_com_1->GetXaxis()->SetTitle("#Delta Time [#mus]");
    h_com_1->GetYaxis()->SetTitle("Pixels in a cluster");

    h_com_2->GetXaxis()->SetTitle("Cluster size [pixels]");
    h_com_2->GetYaxis()->SetTitle("#Delta Time [#mus]");

    h_com_3->GetXaxis()->SetTitle("#Delta Time [#mus]");
    h_com_3->GetYaxis()->SetTitle("Pixels in a cluster");

    for(Int_t i = 1; i <= N_PIXELS; i++)
    {
        for(Int_t j = 1; j <= N_PIXELS; j++)
        {
//            h_6->SetBinContent(N_PIXELS/2-j+1,i,h_1->GetBinContent(i,j));//SPS RomanPot Internal
//            h_6->SetBinContent(i,j,h_1->GetBinContent(i,j));
            h_6->SetBinContent(j,N_PIXELS-i+1,h_1->GetBinContent(i,j));//H8 september 2018
        }
    }

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_5->Write();
    h_6->Write();
    h_8->Write();
    h_9->Write();
    h_10->Write();
    h_13->Write();
    h_15->Write();
    h_16->Write();
    h_20->Write();
    h_21->Write();
    h_22->Write();
    h_23->Write();
    h_24->Write();
    h_25->Write();
    h_26->Write();
    h_27->Write();
    h_28->Write();
    h_29->Write();
    h_31->Write();

    h_com_1->Write();
    h_com_2->Write();
    h_com_3->Write();
    h_com_4->Write();

    outputfile->Close();
    cout<<"--> The output file '"<<outputfile->GetName()<<"' is closed."<<endl;
    fChain->Delete();

    stop_time = time(NULL);
    cout<<"--> Running time is : "<<stop_time - start_time<<" seconds"<<endl;
}

void function_2()
{
    // ChipID: 0 RP0Ext G02-W0108 | 1 RP1Int F04-W0108 | 2 RP3Ext C08-W0255 | 3 RP0Int I02-W0108

    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_RUN_1.root";
    TString outFileName     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_RUN_1_outputfile_function_2.root";
    Int_t chipID            = 1;

    // Stable Operation
//    Long64_t minUnixTime_run    = 1543248000;
//    Long64_t maxUnixTime_run    = 1543269600;
//    Long64_t entryINI           = 0;

    // Not Stable Operation
    Long64_t minUnixTime_run    = 1543271400;
    Long64_t maxUnixTime_run    = 1543277700;
    Long64_t entryINI           = 27700;

    Bool_t entryINIstatus       = kFALSE;

    cout<<endl<<"--> function_1() <--"<<endl<<endl;

    time_t start_time, stop_time;
    start_time = time(NULL);


    // COMMON VARIABLES
    TChain* fChain;
    Long64_t nEntries;

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    UInt_t _AcquisType = -1, _TrigType = -1;
    Long64_t _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS][N_PIXELS];

    // For TIMEPIX

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
    cout<<"--> ChipID: "<<chipID<<endl;
    cout<<"--> InputFileName1: "<<inFileName1<<endl;
    cout<<"--> Number of nEntries: "<<nEntries<<endl;

    // For Output file
    cout<<"--> OutputFileName: "<<outFileName<<endl;
    TFile* outputfile = new TFile(outFileName.Data(),"RECREATE");
    if (outputfile->IsZombie())
    {
        cout<<"--> ERROR: The output file is zombie!"<<endl;
        assert(0);
    }

    //------------------------------------------------------------------------------//
    //------------------------ For output file with clusters infor -----------------//
    //------------------------------------------------------------------------------//
    Double_t unix_time_clinfo, pos_x_clinfo, pos_x_err_clinfo, pos_y_clinfo, pos_y_err_clinfo, clocks_clinfo;
    Int_t size_clinfo;
    Long64_t event_id_clinfo;

    TTree* tree = new TTree("Tree", "Cluster data");

    tree->Branch("UnixTime",        &unix_time_clinfo,      "unix_time_clinfo/D");
    tree->Branch("ClusterClocks",   &clocks_clinfo,         "clocks_clinfo/D");
    tree->Branch("ClusterSize",     &size_clinfo,           "size_clinfo/I");
    tree->Branch("ClusterPosX",     &pos_x_clinfo,          "pos_x_clinfo/D");
    tree->Branch("ClusterPosXerr",  &pos_x_err_clinfo,      "pos_x_err_clinfo/D");
    tree->Branch("ClusterPosY",     &pos_y_clinfo,          "pos_y_clinfo/D");
    tree->Branch("ClusterPosYerr",  &pos_y_err_clinfo,      "pos_y_err_clinfo/D");
    tree->Branch("EventID",         &event_id_clinfo,       "event_id_clinfo/L");
    tree->Branch("Clock",           &_Clock,                "_Clock/D");
    tree->Branch("Gate",            &_Gate,                 "_Gate/D");

    Int_t* cl_size          = new Int_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_x      = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_x_err  = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_y      = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_pos_y_err  = new Double_t[N_MAX_CLUSTERS];
    Double_t* cl_clocks     = new Double_t[N_MAX_CLUSTERS];
    //------------------------------------------------------------------------------//
    //------------------------------- HISTOGRAMS -----------------------------------//
    //------------------------------------------------------------------------------//

    fChain->GetEntry(nEntries-1);
    Double_t UT_max = _Timems;

    fChain->GetEntry(0);
    Long64_t _event_ini = _event;
    Double_t UT_min = _Timems;
    cout<<"--> AcquisType = "<<_AcquisType<<" (0 - MPX, 1 - ToT, 3 - ToA)"<<endl;
    cout<<"--> Clock = "<<_Clock<<" [MHz]"<<endl;
    cout<<"--> Gate = "<<_Gate<<" [sec]"<<endl;

    TH2D* h_1   = new TH2D("h_1","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH1D* h_2   = new TH1D("h_2","COUNTS vs X",N_PIXELS,0,N_PIXELS);
    TH1D* h_3   = new TH1D("h_3","COUNTS vs Y",N_PIXELS,0,N_PIXELS);
    TH2D* h_4   = new TH2D("h_4","X vs EventID",nEntries,0,nEntries,N_PIXELS,0,N_PIXELS);
    TH2D* h_5   = new TH2D("h_5","Y vs EventID",nEntries,0,nEntries,N_PIXELS,0,N_PIXELS);
    TH2D* h_6   = new TH2D("h_6","Y vs X vs C (in mm)",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);
    TH2D* h_8   = new TH2D("h_8","Y vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH2D* h_9   = new TH2D("h_9","X vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH1D* h_10  = new TH1D("h_10","Hits vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0);
    TH1D* h_13  = new TH1D("h_13","Clusters number per frame",200,0,200);
    TH1D* h_14  = new TH1D("h_14","Hits vs EventID",nEntries,0,nEntries);
    TH1D* h_15  = new TH1D("h_15","Hits per frame",1e6,0,1e6);
    TH1D* h_16  = new TH1D("h_16","Cluster size",1e4,0,1e4);
    TH2D* h_20  = new TH2D("h_20","Y vs X vs C (inside the ch beam spot)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_21  = new TH2D("h_21","Y vs X vs C (inside the dch beam spot)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_22  = new TH2D("h_22","Y vs X vs C (outside the beam spot)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH1D* h_23  = new TH1D("h_23","#Delta Time between frames",10000,0,10000);
    TH1D* h_24  = new TH1D("h_24","_AcquisType vs EventID",nEntries,0,nEntries);
    TH1D* h_25  = new TH1D("h_25","_DeltaTHR vs EventID",nEntries,0,nEntries);
    TH1D* h_26  = new TH1D("h_26","_Gate vs EventID",nEntries,0,nEntries);
    TH1D* h_27  = new TH1D("h_27","_Bias vs EventID",nEntries,0,nEntries);
    TH2D* h_28  = new TH2D("h_28","Cluster Volume vs Cluster size",1e3,0,1e3,1e4,0,1e6);
    TH1D* h_29  = new TH1D("h_29","Number of Fired Pixels vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0);
    TH1D* h_31  = new TH1D("h_31","Cluster Volume",1e6,0,1e6);
    TH1D* h_32  = new TH1D("h_32","Number of Fired Pixels vs Time (cut)",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0);

    h_com_1     = new TH1D("h_com_1","#Delta Time inside cluster",2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_2     = new TH2D("h_com_2","#Delta Time inside cluster VS Cluster size",50,0,50,2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_3     = new TH1D("h_com_3","#Delta Time inside cluster (cluster size <= 2 pixels)",2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_4     = new TH2D("h_com_4","Y vs X vs ToT (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

    //------------------------------------------------------------------------------//
    Double_t zero_time = -999.999, event_time = -999.999, last_time = -999.999, previous_time = _Timems;
    Long64_t nFrames = 0;
    Long64_t frame_size = 0;

    string time_char;
    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------- MAIN LOOP ------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//

    for(Long64_t i = entryINI; i < nEntries; i++)
    {
        if(i%10 == 0)
        {
            epochtime2date(round(last_time/1000),time_char);
            printf("\r--> Working: %3.1f %% | Time: %s",100*(Double_t)i/nEntries,time_char.c_str());
            fflush(stdout);
        }

        fChain->GetEntry(i);

        //================================================//
        if(_Timems/1000.0 < minUnixTime_run) continue;
        if(_Timems/1000.0 > maxUnixTime_run) break;

        if(!entryINIstatus)
        {
            entryINI = i;
            cout<<endl<<"--> entryINI = "<<entryINI<<endl;
            entryINIstatus = kTRUE;
        }
        //================================================//

        if (_Timems > 0)
        {
            event_time  = _Timems;
        }
        else
        {
            continue;
        }

        if (zero_time < 0)
        {
            zero_time = _Timems;
        }
        last_time = _Timems;

        Int_t fired_pixels_num = 0;
        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0)
                {
                    fired_pixels_num++;
                }
            }
        }

        h_29->Fill(event_time/1000.0,fired_pixels_num);

//        if(fired_pixels_num > 1e4) continue;

        h_32->Fill(event_time/1000.0,fired_pixels_num);

        h_24->Fill(_event-_event_ini,_AcquisType);
        h_25->Fill(_event-_event_ini,_DeltaTHR);
        h_26->Fill(_event-_event_ini,_Gate);
        h_27->Fill(_event-_event_ini,_Bias);

        frame_size = 0;
        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0)
                {
                    if(_AcquisType == 0 || _AcquisType == 1)
                    {
                        frame_size += _COUNTS[xi][yi];

                        h_1->Fill(xi,yi,_COUNTS[xi][yi]);
                        h_2->Fill(xi,_COUNTS[xi][yi]);
                        h_3->Fill(yi,_COUNTS[xi][yi]);
                        h_4->Fill(_event-_event_ini,xi,_COUNTS[xi][yi]);
                        h_5->Fill(_event-_event_ini,yi,_COUNTS[xi][yi]);
                        h_8->Fill(event_time/1000.0,yi,_COUNTS[xi][yi]);
                        h_9->Fill(event_time/1000.0,xi,_COUNTS[xi][yi]);
                    }
                    else
                    {
                        cout<<endl<<endl<<"--> ERROR:: _AcquisType is neither 0 nor 1"<<endl<<endl;
                        assert(0);
                    }
                }
            }
        }
        h_10->Fill(event_time/1000.0,frame_size);
        h_14->Fill(i,frame_size);
        h_15->Fill(frame_size);

        h_23->Fill(event_time - previous_time);
        previous_time = event_time;
        //---------------------------------------------------------------------------------------------------------------//
        //------------------------------------------ CLUSTER ANALYSIS ---------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        if(_cluster_analysis)
        {
            for(Int_t jk = 0; jk < N_MAX_CLUSTERS; jk++)
            {
                cl_size[jk]         = 0;
                cl_pos_x[jk]        = 0.0;
                cl_pos_x_err[jk]    = 0.0;
                cl_pos_y[jk]        = 0.0;
                cl_pos_y_err[jk]    = 0.0;
                cl_clocks[jk]       = 0.0;
            }

            Int_t cl_num = 0;

            if(clusteranalysis(_COUNTS,cl_num,cl_size,cl_pos_x,cl_pos_x_err,cl_pos_y,cl_pos_y_err,cl_clocks,_AcquisType,_Clock) != 0)
            {
                cout<<"--> ERROR: Cluster analysis did not check all pixels! ("<<clusteranalysis(_COUNTS,cl_num,cl_size,cl_pos_x,cl_pos_x_err,cl_pos_y,cl_pos_y_err,cl_clocks,_AcquisType,_Clock)<<")"<<endl;
                assert(0);
            }

            if(cl_num >= N_MAX_CLUSTERS)
            {
                cout<<endl<<"--> ERROR:: The maximum number of clusters per a frame was reached !!!"<<endl;
                assert(0);
            }

            for(Int_t yy = 0; yy < cl_num; yy++)
            {
                //------------------------------------------------------------------------------//
                //------------------------ For output file with clusters infor -----------------//
                //------------------------------------------------------------------------------//
                unix_time_clinfo    = event_time/1000.0;

                if(_AcquisType == 3) clocks_clinfo = (Double_t)N_MAX_CLOCKS - cl_clocks[yy];
                else clocks_clinfo = (Double_t)cl_clocks[yy];

                if(clocks_clinfo < 0) clocks_clinfo = 0;
                size_clinfo         = cl_size[yy];
                pos_x_clinfo        = cl_pos_x[yy];
                pos_x_err_clinfo    = cl_pos_x_err[yy];
                pos_y_clinfo        = cl_pos_y[yy];
                pos_y_err_clinfo    = cl_pos_y_err[yy];
                event_id_clinfo     = _event;
                if(_event != i) {event_id_clinfo = i;}

                tree->Fill();

                h_16->Fill(cl_size[yy]);
                h_28->Fill(cl_size[yy],cl_size[yy]*cl_clocks[yy]);
                h_31->Fill(cl_size[yy]*cl_clocks[yy]);

            }
            h_13->Fill(cl_num);
        }
        //---------------------------------------------------------------------------------------------------------------//
        nFrames++;
    }
    cout<<endl;

    delete [] cl_size;
    delete [] cl_pos_x;
    delete [] cl_pos_x_err;
    delete [] cl_pos_y;
    delete [] cl_pos_y_err;
    delete [] cl_clocks;

    cout<<endl;

    Double_t scfactor = _Gate*nFrames;
    cout<<"--> nFrames = "<<nFrames<<endl;
    cout<<"--> scaling factor = "<<scfactor<<" [sec]"<<endl;

    string zero_time_char;
    string last_time_char;
    epochtime2date(round(zero_time/1000),zero_time_char);
    epochtime2date(round(last_time/1000),last_time_char);

    cout<<endl;
    cout<<"--> Initial time: "<<zero_time_char<<" | Final time: "<<last_time_char<<endl;

    //------------------------------------------------------------------------------//
    //------------------------ For output file with clusters infor -----------------//
    //------------------------------------------------------------------------------//
    tree->Write();
    cout<<"--> The Tree with "<<tree->GetEntriesFast()<<" entries is written to the file: "<<outputfile->GetName()<<endl;
    //------------------------------------------------------------------------------//

    h_1->GetXaxis()->SetTitle("X [pixels]");
    h_1->GetYaxis()->SetTitle("Y [pixels]");
    h_2->GetXaxis()->SetTitle("X [pixels]");
    h_2->GetYaxis()->SetTitle("Counts");
    h_3->GetXaxis()->SetTitle("Y [pixels]");
    h_3->GetYaxis()->SetTitle("Counts");
    h_4->GetXaxis()->SetTitle("EventID");
    h_4->GetYaxis()->SetTitle("X [pixels]");
    h_5->GetXaxis()->SetTitle("EventID");
    h_5->GetYaxis()->SetTitle("Y [pixels]");
    h_6->GetXaxis()->SetTitle("X [mm]");
    h_6->GetYaxis()->SetTitle("Y [mm]");
    h_8->GetXaxis()->SetTitle("Time [sec]");
    h_8->GetYaxis()->SetTitle("Y [pixels]");
    h_9->GetXaxis()->SetTitle("Time [sec]");
    h_9->GetYaxis()->SetTitle("X [pixels]");
    h_10->GetXaxis()->SetTitle("Time [sec]");
    h_10->GetYaxis()->SetTitle("Hits");
    h_13->GetXaxis()->SetTitle("Clusters number per frame");
    h_13->GetYaxis()->SetTitle("Numer of frames");
    h_14->GetXaxis()->SetTitle("EventID");
    h_14->GetYaxis()->SetTitle("Hits");
    h_15->GetXaxis()->SetTitle("Hits per frame");
    h_15->GetYaxis()->SetTitle("Frames");
    h_16->GetXaxis()->SetTitle("Cluster size [pixels/cluster]");
    h_16->GetYaxis()->SetTitle("Number of clusters");
    h_20->GetXaxis()->SetTitle("X [pixels]");
    h_20->GetYaxis()->SetTitle("Y [pixels]");
    h_21->GetXaxis()->SetTitle("X [pixels]");
    h_21->GetYaxis()->SetTitle("Y [pixels]");
    h_22->GetXaxis()->SetTitle("X [pixels]");
    h_22->GetYaxis()->SetTitle("Y [pixels]");
    h_23->GetXaxis()->SetTitle("#Delta Time [ms]");
    h_23->GetYaxis()->SetTitle("Number of frames");
    h_24->GetXaxis()->SetTitle("EventID");
    h_24->GetYaxis()->SetTitle("_AcquisType");
    h_25->GetXaxis()->SetTitle("EventID");
    h_25->GetYaxis()->SetTitle("_DeltaTHR");
    h_26->GetXaxis()->SetTitle("EventID");
    h_26->GetYaxis()->SetTitle("_Gate");
    h_27->GetXaxis()->SetTitle("EventID");
    h_27->GetYaxis()->SetTitle("_Bias");
    h_28->GetXaxis()->SetTitle("Cluster size [pixels/cluster]");
    h_28->GetYaxis()->SetTitle("Cluster volume [counts/cluster]");
    h_31->GetYaxis()->SetTitle("Cluster volume [counts/cluster]");

    h_com_1->GetXaxis()->SetTitle("#Delta Time [#mus]");
    h_com_1->GetYaxis()->SetTitle("Pixels in a cluster");

    h_com_2->GetXaxis()->SetTitle("Cluster size [pixels]");
    h_com_2->GetYaxis()->SetTitle("#Delta Time [#mus]");

    h_com_3->GetXaxis()->SetTitle("#Delta Time [#mus]");
    h_com_3->GetYaxis()->SetTitle("Pixels in a cluster");

    for(Int_t i = 1; i <= N_PIXELS; i++)
    {
        for(Int_t j = 1; j <= N_PIXELS; j++)
        {
//            h_6->SetBinContent(N_PIXELS/2-j+1,i,h_1->GetBinContent(i,j));//SPS RomanPot Internal
//            h_6->SetBinContent(i,j,h_1->GetBinContent(i,j));
            h_6->SetBinContent(j,N_PIXELS-i+1,h_1->GetBinContent(i,j));//H8 september 2018
        }
    }

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_5->Write();
    h_6->Write();
    h_8->Write();
    h_9->Write();
    h_10->Write();
    h_13->Write();
    h_15->Write();
    h_16->Write();
    h_20->Write();
    h_21->Write();
    h_22->Write();
    h_23->Write();
    h_24->Write();
    h_25->Write();
    h_26->Write();
    h_27->Write();
    h_28->Write();
    h_29->Write();
    h_31->Write();
    h_32->Write();

    h_com_1->Write();
    h_com_2->Write();
    h_com_3->Write();
    h_com_4->Write();

    outputfile->Close();
    cout<<"--> The output file '"<<outputfile->GetName()<<"' is closed."<<endl;
    fChain->Delete();

    stop_time = time(NULL);
    cout<<"--> Running time is : "<<stop_time - start_time<<" seconds"<<endl;
}

void epochtime2date(time_t in_time, string &out_time)
{
    const char default_format[] = "%a %b %d %Y %H:%M:%S";
    time_t t = in_time;
    const char *format = default_format;

    struct tm lt;
    char res[32];

    (void) localtime_r(&t, &lt);

    if (strftime(res, sizeof(res), format, &lt) == 0)
    {
        (void) fprintf(stderr,  "strftime(3): cannot format supplied "
                                "date/time into buffer of size %lu "
                                "using: '%s'\n",
                       sizeof(res), format);
    }

    out_time = res;
}

int clusteranalysis(Long64_t _matrix[][N_PIXELS], Int_t &cluster_num, Int_t *cluster_size, Double_t *cluster_pos_x, Double_t *cluster_pos_x_err, Double_t *cluster_pos_y, Double_t *cluster_pos_y_err, Double_t *cluster_clocks, Int_t AcqType, Double_t Clock)
{
    Int_t** fired_matrix = new Int_t*[N_PIXELS];// matrix of the fired pixels for each cluster
    Int_t** fired_matrix_full = new Int_t*[N_PIXELS];// matrix of the fired pixels for frame
    for(Int_t i = 0; i < N_PIXELS; i++)
    {
        fired_matrix[i] = new Int_t[N_PIXELS];
        fired_matrix_full[i] = new Int_t[N_PIXELS];
        for(Int_t j = 0; j < N_PIXELS; j++)
        {
            fired_matrix[i][j] = 0;
            fired_matrix_full[i][j] = 0;
        }
    }

    cluster_num = 0;
    for(Int_t yi = 0; yi < N_PIXELS; yi++)
    {
        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            if(_matrix[xi][yi] > 0 && fired_matrix_full[xi][yi] == 0)
            {
                switch (AcqType)
                {
                case 0:
                    checkregionMPX(_matrix,fired_matrix,xi,yi);
                    break;
                case 1:
                    checkregionMPX(_matrix,fired_matrix,xi,yi);
                    break;
                case 3:
                    checkregionTOA(_matrix,fired_matrix,xi,yi,Clock);
                    break;
                default:
                    cout<<endl<<endl<<"--> ERROR:: _AcquisType < 0"<<endl<<endl;
                    assert(0);
                }

                // Projections of the cluster
                TH1D* h_x = new TH1D("h_x","h_x",N_PIXELS,0,N_PIXELS);
                TH1D* h_y = new TH1D("h_y","h_y",N_PIXELS,0,N_PIXELS);

                Long64_t ref_clock = -999;
                Int_t xj_ref_clock = -1, yj_ref_clock = -1;
                for(Int_t yj = 0; yj < N_PIXELS; yj++)
                {
                    for(Int_t xj = 0; xj < N_PIXELS; xj++)
                    {
                        if(fired_matrix[xj][yj] > 0)
                        {
                            ref_clock = _matrix[xj][yj];
                            xj_ref_clock = xj;
                            yj_ref_clock = yj;
                        }
                        if(ref_clock != -999) {break;}
                    }
                    if(ref_clock != -999) {break;}
                }

                for(Int_t yj = 0; yj < N_PIXELS; yj++)
                {
                    for(Int_t xj = 0; xj < N_PIXELS; xj++)
                    {
                        if(fired_matrix[xj][yj] > 0)
                        {
                            cluster_size[cluster_num]++;
                        }
                    }
                }

                for(Int_t yj = 0; yj < N_PIXELS; yj++)
                {
                    for(Int_t xj = 0; xj < N_PIXELS; xj++)
                    {
                        if(fired_matrix[xj][yj] > 0)
                        {
                            h_x->Fill(xj,_matrix[xj][yj]);
                            h_y->Fill(yj,_matrix[xj][yj]);

                            fired_matrix_full[xj][yj] = 1;
                            cluster_clocks[cluster_num] += _matrix[xj][yj];
                            if(AcqType == 3)
                            {
                                if(xj == xj_ref_clock && yj == yj_ref_clock)
                                {;}
                                else
                                {
                                    h_com_1->Fill((ref_clock-_matrix[xj][yj])/Clock);
                                    h_com_2->Fill(cluster_size[cluster_num],(ref_clock-_matrix[xj][yj])/Clock);
                                    if(cluster_size[cluster_num] <= 2) h_com_3->Fill((ref_clock-_matrix[xj][yj])/Clock);
                                }
                            }
                        }
                    }
                }

                if(cluster_clocks[cluster_num] > 100e3 && oneCluster)
                {
                    for(Int_t yj = 0; yj < N_PIXELS; yj++)
                    {
                        for(Int_t xj = 0; xj < N_PIXELS; xj++)
                        {
                            if(fired_matrix[xj][yj] > 0)
                            {
                                h_com_4->Fill(xj,yj,_matrix[xj][yj]);
                            }
                        }
                    }

                    oneCluster = kFALSE;
                }

                for(Int_t yj = 0; yj < N_PIXELS; yj++)
                {
                    for(Int_t xj = 0; xj < N_PIXELS; xj++)
                    {
                        fired_matrix[xj][yj] = 0;
                    }
                }

                cluster_clocks[cluster_num]     = (Double_t)cluster_clocks[cluster_num]/cluster_size[cluster_num];
                cluster_pos_x[cluster_num]      = h_x->GetMean();
                cluster_pos_x_err[cluster_num]  = h_x->GetMeanError();
                cluster_pos_y[cluster_num]      = h_y->GetMean();
                cluster_pos_y_err[cluster_num]  = h_y->GetMeanError();

                h_x->Delete();
                h_y->Delete();

                cluster_num++;
            }
        }
    }

    Int_t check_value = 0;
    // Check if cluster analysis took all pixels
    for(Int_t yi = 0; yi < N_PIXELS; yi++)
    {
        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            if(_matrix[xi][yi] != 0 && fired_matrix_full[xi][yi] == 0)
                check_value++;
            if(_matrix[xi][yi] == 0 && fired_matrix_full[xi][yi] != 0)
                check_value++;
        }
    }

    for(Int_t i = 0; i < N_PIXELS; i++)
    {
        delete [] fired_matrix[i];
        delete [] fired_matrix_full[i];
    }
    delete [] fired_matrix;
    delete [] fired_matrix_full;

    return check_value;
}

void checkregionMPX(Long64_t _matrix[][N_PIXELS], Int_t **fired_matrix, Int_t xi, Int_t yi)
{
    fired_matrix[xi][yi] = 1;
    // Check neighbor pixels
    if(xi-1 >= 0 && yi-1 >=0)
        if(_matrix[xi-1][yi-1]  > 0 && fired_matrix[xi-1][yi-1] == 0)   checkregionMPX(_matrix,fired_matrix,xi-1,yi-1);
    if(yi-1 >= 0)
        if(_matrix[xi][yi-1]    > 0 && fired_matrix[xi][yi-1] == 0)     checkregionMPX(_matrix,fired_matrix,xi,yi-1);
    if(xi+1 < N_PIXELS  && yi-1 >=0)
        if(_matrix[xi+1][yi-1]  > 0 && fired_matrix[xi+1][yi-1] == 0)   checkregionMPX(_matrix,fired_matrix,xi+1,yi-1);
    if(xi-1 >= 0)
        if(_matrix[xi-1][yi]    > 0 && fired_matrix[xi-1][yi] == 0)     checkregionMPX(_matrix,fired_matrix,xi-1,yi);
    if(xi+1 < N_PIXELS)
        if(_matrix[xi+1][yi]    > 0 && fired_matrix[xi+1][yi] == 0)     checkregionMPX(_matrix,fired_matrix,xi+1,yi);
    if(xi-1 >= 0 && yi+1 < N_PIXELS)
        if(_matrix[xi-1][yi+1]  > 0 && fired_matrix[xi-1][yi+1] == 0)   checkregionMPX(_matrix,fired_matrix,xi-1,yi+1);
    if(yi+1 < N_PIXELS)
        if(_matrix[xi][yi+1]    > 0 && fired_matrix[xi][yi+1] == 0)     checkregionMPX(_matrix,fired_matrix,xi,yi+1);
    if(xi+1 < N_PIXELS  && yi+1 < N_PIXELS)
        if(_matrix[xi+1][yi+1]  > 0 && fired_matrix[xi+1][yi+1] == 0)   checkregionMPX(_matrix,fired_matrix,xi+1,yi+1);
}

void checkregionTOA(Long64_t _matrix[][N_PIXELS], Int_t **fired_matrix, Int_t xi, Int_t yi, Double_t Clock)
{
    fired_matrix[xi][yi] = 1;
    // Check neighbor pixels
    if(xi-1 >= 0 && yi-1 >=0)
        if(checkcondition(_matrix[xi-1][yi-1],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi-1][yi-1] == 0)    checkregionTOA(_matrix,fired_matrix,xi-1,yi-1,Clock);
    if(yi-1 >= 0)
        if(checkcondition(_matrix[xi][yi-1],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi][yi-1] == 0)        checkregionTOA(_matrix,fired_matrix,xi,yi-1,Clock);
    if(xi+1 < N_PIXELS  && yi-1 >=0)
        if(checkcondition(_matrix[xi+1][yi-1],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi+1][yi-1] == 0)    checkregionTOA(_matrix,fired_matrix,xi+1,yi-1,Clock);
    if(xi-1 >= 0)
        if(checkcondition(_matrix[xi-1][yi],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi-1][yi] == 0)        checkregionTOA(_matrix,fired_matrix,xi-1,yi,Clock);
    if(xi+1 < N_PIXELS)
        if(checkcondition(_matrix[xi+1][yi],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi+1][yi] == 0)        checkregionTOA(_matrix,fired_matrix,xi+1,yi,Clock);
    if(xi-1 >= 0 && yi+1 < N_PIXELS)
        if(checkcondition(_matrix[xi-1][yi+1],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi-1][yi+1] == 0)    checkregionTOA(_matrix,fired_matrix,xi-1,yi+1,Clock);
    if(yi+1 < N_PIXELS)
        if(checkcondition(_matrix[xi][yi+1],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi][yi+1] == 0)        checkregionTOA(_matrix,fired_matrix,xi,yi+1,Clock);
    if(xi+1 < N_PIXELS  && yi+1 < N_PIXELS)
        if(checkcondition(_matrix[xi+1][yi+1],(_matrix[xi][yi]),cluster_time_jitter,Clock) && fired_matrix[xi+1][yi+1] == 0)    checkregionTOA(_matrix,fired_matrix,xi+1,yi+1,Clock);
}

bool checkcondition(Long64_t x, Long64_t ref_x, Double_t delta, Double_t Clock)
{
    Double_t x_time = x*1e3/Clock;
    Double_t x_ref_time = ref_x*1e3/Clock;
    if(TMath::Abs(x_time - x_ref_time) < delta && x > 0) return true;
//    if(x > 0) return true;
    return false;
}
