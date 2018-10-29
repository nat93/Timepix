//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
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

int main(int argc, char *argv[])
{
    cout<<endl<<"--> ANALYSIS COMMON"<<endl<<endl;
    if(argc != 4)
    {
        cout<<"--> ERROR: Wrong number of input parameters("<<argc<<")!"<<endl<<endl;

        cout<<"--> [0] ./script_name"<<endl;
        cout<<"--> [1] Input file name"<<endl;
        cout<<"--> [2] Output file name"<<endl;
        cout<<"--> [3] ChipID"<<endl;

        return -1;
    }

    time_t start_time, stop_time;
    start_time = time(NULL);

    Int_t chipID = atoi(argv[3]);

    // COMMON VARIABLES
    TChain* fChain;
    Long64_t nEntries;

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    UInt_t _AcquisType = -1, _TrigType = -1;
    Long64_t _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS][N_PIXELS];

    // For TIMEPIX
    TString inFileName1     = argv[1];

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
    TString outFileName     = argv[2];
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
    TH1D* h_7   = new TH1D("h_7","Int.FrameSize vs ToA",N_MAX_CLOCKS,0,_Gate*1000);
    TH2D* h_8   = new TH2D("h_8","Y vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH2D* h_9   = new TH2D("h_9","X vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH1D* h_10  = new TH1D("h_10","Hits vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0);
    TH1D* h_12  = new TH1D("h_12","TOA (only one frame)",N_MAX_CLOCKS,0,_Gate*1e6);
    TH1D* h_13  = new TH1D("h_13","Clusters number per frame",200,0,200);
    TH1D* h_14  = new TH1D("h_14","Hits vs EventID",nEntries,0,nEntries);
    TH1D* h_15  = new TH1D("h_15","Hits per frame",1e6,0,1e6);
    TH1D* h_16  = new TH1D("h_16","Cluster size",1e4,0,1e4);
    TH1D* h_17  = new TH1D("h_17","Int.FrameSize vs ToA (inside the ch beam spot)",N_MAX_CLOCKS,0,_Gate*1000);
    TH1D* h_18  = new TH1D("h_18","Int.FrameSize vs ToA (inside the dch beam spot)",N_MAX_CLOCKS,0,_Gate*1000);
    TH1D* h_19  = new TH1D("h_19","Int.FrameSize vs ToA (outside the beam spot)",N_MAX_CLOCKS,0,_Gate*1000);
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
    TH1D* h_30  = new TH1D("h_30","Number of Fired Pixels per frame (cut)",3e5,0,3e5);
    TH1D* h_31  = new TH1D("h_31","Cluster Volume",1e5,0,1e5);

    h_com_1     = new TH1D("h_com_1","#Delta Time inside cluster",2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_2     = new TH2D("h_com_2","#Delta Time inside cluster VS Cluster size",50,0,50,2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);
    h_com_3     = new TH1D("h_com_3","#Delta Time inside cluster (cluster size <= 2 pixels)",2*N_MAX_CLOCKS-1,-_Gate*1e6,_Gate*1e6);

    //------------------------------------------------------------------------------//
    Double_t zero_time = -999.999, event_time = -999.999, last_time = -999.999, previous_time = _Timems;
    Long64_t nFrames = 0;
    Long64_t frame_size = 0;

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------- MAIN LOOP ------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
//    nEntries /= 50;
    for(Long64_t i = 0; i < nEntries; i++)
    {
        if(i%1 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntries);
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

        if(fired_pixels_num < 1e4 || _AcquisType == 0)
        {
            h_30->Fill(fired_pixels_num);
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
                        else if(_AcquisType == 3)
                        {
    //                        _Clock = 4.8;
                            Double_t TOA = (N_MAX_CLOCKS-_COUNTS[xi][yi])*1e-6/_Clock;  // [sec]
                            frame_size++;

    //                        if(TOA <= 0.000246)
                            if(1)
                            {
                                h_1->Fill(xi,yi,1);
                                h_2->Fill(xi,1);
                                h_3->Fill(yi,1);
                                h_4->Fill(_event-_event_ini,xi,1);
                                h_5->Fill(_event-_event_ini,yi,1);

                                for(Int_t fff = 1; fff <= h_7->GetNbinsX(); fff++)
                                {
                                    if(TOA*1000 < h_7->GetBinCenter(fff))
                                    {
                                        h_7->SetBinContent(fff,h_7->GetBinContent(fff)+1);

                                        //------------------------------------------------------//
                                        // For SPS RP0I Timepix
                                        if(xi >= 100 && yi >= 20 && yi <= 150) // inside CH
                                        {
                                            h_17->SetBinContent(fff,h_17->GetBinContent(fff)+1);
                                            h_20->Fill(xi,yi,1);
                                        }
                                        else if(xi >= 150 && yi > 150)   // inside DCH
                                        {
                                            h_18->SetBinContent(fff,h_18->GetBinContent(fff)+1);
                                            h_21->Fill(xi,yi,1);
                                        }
                                        else    // outside the beam spot
                                        {
                                            h_19->SetBinContent(fff,h_19->GetBinContent(fff)+1);
                                            h_22->Fill(xi,yi,1);
                                        }
                                        //------------------------------------------------------//
                                    }
                                }

                                h_8->Fill(event_time/1000.0,yi,1);
                                h_9->Fill(event_time/1000.0,xi,1);
                                if(i == 190) h_12->Fill(TOA*1e6);
                            }
                        }
                        else
                        {
                            cout<<endl<<endl<<"--> ERROR:: _AcquisType is not 0 or 3"<<endl<<endl;
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
    h_7->GetXaxis()->SetTitle("ToA [ms]");
    h_7->GetYaxis()->SetTitle("Integrated Frame Size [counts]");
    h_8->GetXaxis()->SetTitle("Time [sec]");
    h_8->GetYaxis()->SetTitle("Y [pixels]");
    h_9->GetXaxis()->SetTitle("Time [sec]");
    h_9->GetYaxis()->SetTitle("X [pixels]");
    h_10->GetXaxis()->SetTitle("Time [sec]");
    h_10->GetYaxis()->SetTitle("Hits");
    h_12->GetXaxis()->SetTitle("ToA [#mus]");
    h_12->GetYaxis()->SetTitle("Number of fired pixels");
    h_13->GetXaxis()->SetTitle("Clusters number per frame");
    h_13->GetYaxis()->SetTitle("Numer of frames");
    h_14->GetXaxis()->SetTitle("EventID");
    h_14->GetYaxis()->SetTitle("Hits");
    h_15->GetXaxis()->SetTitle("Hits per frame");
    h_15->GetYaxis()->SetTitle("Frames");
    h_16->GetXaxis()->SetTitle("Cluster size [pixels/cluster]");
    h_16->GetYaxis()->SetTitle("Number of clusters");
    h_17->GetXaxis()->SetTitle("ToA [ms]");
    h_17->GetYaxis()->SetTitle("Integrated Frame Size [counts]");
    h_18->GetXaxis()->SetTitle("ToA [ms]");
    h_18->GetYaxis()->SetTitle("Integrated Frame Size [counts]");
    h_19->GetXaxis()->SetTitle("ToA [ms]");
    h_19->GetYaxis()->SetTitle("Integrated Frame Size [counts]");
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
            h_6->SetBinContent(i,j,h_1->GetBinContent(i,j));
//            h_6->SetBinContent(j,N_PIXELS-i+1,h_1->GetBinContent(i,j));//H8 september 2018
        }
    }

    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_5->Write();
    h_6->Write();
    h_7->Write();
    h_8->Write();
    h_9->Write();
    h_10->Write();
    h_12->Write();
    h_13->Write();
    h_14->Write();
    h_15->Write();
    h_16->Write();
    h_17->Write();
    h_18->Write();
    h_19->Write();
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
    h_30->Write();
    h_31->Write();

    h_com_1->Write();
    h_com_2->Write();
    h_com_3->Write();

    outputfile->Close();
    cout<<"--> The output file '"<<outputfile->GetName()<<"' is closed."<<endl;
    fChain->Delete();

    stop_time = time(NULL);
    cout<<"--> Running time is : "<<stop_time - start_time<<" seconds"<<endl;
    return 0;
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
