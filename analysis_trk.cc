//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     Cluster Analysis for the Tracking at H8 for 4 planes.
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

const Int_t N_PIXELS            = 512;          // The maximum number of pixels per axis per chip
const Int_t N_MAX_CHIP          = 4;            // The maximum number of chips
const Int_t N_MAX_CLUSTERS_1    = 10000;        // The maximum number of clusters per frame
const Int_t N_MAX_CLUSTERS_2    = 30;           // The maximum number of clusters per clock
const Int_t N_MAX_CLOCKS        = 11810;        // The maximum number of counts in the pixel
const Bool_t _cluster_analysis  = true;
const Double_t cluster_time_jitter  = 200.0;    // Time difference between fired pixels inside cluster [ns]

void epochtime2date(time_t in_time, string &out_time);
int clusteranalysis(Long64_t _matrix[][N_PIXELS], Int_t &cluster_num, Int_t *cluster_size, Double_t *cluster_pos_x, Double_t *cluster_pos_x_err, Double_t *cluster_pos_y, Double_t *cluster_pos_y_err, Double_t *cluster_clocks, Int_t AcqType, Double_t Clock);
void checkregionMPX(Long64_t _matrix[][N_PIXELS], Int_t **fired_matrix, Int_t xi, Int_t yi);
void checkregionTOA(Long64_t _matrix[][N_PIXELS], Int_t **fired_matrix, Int_t xi, Int_t yi, Double_t Clock);
bool checkcondition(Long64_t x, Long64_t ref_x, Double_t delta, Double_t Clock);

int main(int argc, char *argv[])
{
    cout<<endl<<"--> ANALYSIS TRACKING"<<endl<<endl;
    if(argc != 3)
    {
        cout<<"--> ERROR: Wrong number of input parameters("<<argc<<")!"<<endl<<endl;

        cout<<"--> [0] ./script_name"<<endl;
        cout<<"--> [1] Input file name"<<endl;
        cout<<"--> [2] Output file name"<<endl;

        return -1;
    }

    time_t start_time, stop_time;
    start_time = time(NULL);

    // COMMON VARIABLES
    TChain* fChain[N_MAX_CHIP];
    Long64_t nEntries[N_MAX_CHIP];

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate;
    Int_t _TrigType, _AcquisType, _event;
    Double_t _Timems;
    Long64_t _COUNTS[N_PIXELS][N_PIXELS];

    // For TIMEPIX
    TString inFileName1     = argv[1];

    TString tree_name;
    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        tree_name = "Tree_";
        tree_name += i;
        fChain[i] = new TChain(tree_name.Data());
        fChain[i]->Add(inFileName1.Data());

        fChain[i]->SetBranchAddress("UnixTime",    &_UnixTime);
        fChain[i]->SetBranchAddress("TrigType",    &_TrigType);
        fChain[i]->SetBranchAddress("AcquisType",  &_AcquisType);
        fChain[i]->SetBranchAddress("DeltaTHR",    &_DeltaTHR);
        fChain[i]->SetBranchAddress("Ikrum",       &_Ikrum);
        fChain[i]->SetBranchAddress("Bias",        &_Bias);
        fChain[i]->SetBranchAddress("Clock",       &_Clock);
        fChain[i]->SetBranchAddress("Gate",        &_Gate);
        fChain[i]->SetBranchAddress("Event",       &_event);
        fChain[i]->SetBranchAddress("Time_ms",     &_Timems);
        fChain[i]->SetBranchAddress("COUNTS",      _COUNTS);

        nEntries[i] = (fChain[i]->GetEntries());
        cout<<"--> ChipID: "<<i<<endl;
        cout<<"--> InputFileName1: "<<inFileName1<<endl;
        cout<<"--> Number of nEntries: "<<nEntries[i]<<endl<<endl;
    }

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
    Int_t size_clinfo, event_id_clinfo;

    TTree* tree[N_MAX_CHIP];
    UInt_t cluster_id[N_MAX_CHIP] = {};

    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        tree_name = "Tree_"; tree_name += i;
        tree[i] = new TTree(tree_name.Data(), "Cluster data");

        tree[i]->Branch("UnixTime",        &unix_time_clinfo,      "unix_time_clinfo/D");
        tree[i]->Branch("ClusterClocks",   &clocks_clinfo,         "clocks_clinfo/D");
        tree[i]->Branch("ClusterSize",     &size_clinfo,           "size_clinfo/I");
        tree[i]->Branch("ClusterPosX",     &pos_x_clinfo,          "pos_x_clinfo/D");
        tree[i]->Branch("ClusterPosXerr",  &pos_x_err_clinfo,      "pos_x_err_clinfo/D");
        tree[i]->Branch("ClusterPosY",     &pos_y_clinfo,          "pos_y_clinfo/D");
        tree[i]->Branch("ClusterPosYerr",  &pos_y_err_clinfo,      "pos_y_err_clinfo/D");
        tree[i]->Branch("EventID",         &event_id_clinfo,       "event_id_clinfo/L");
        tree[i]->Branch("EventID",         &event_id_clinfo,       "event_id_clinfo/L");
        tree[i]->Branch("Clock",           &_Clock,                "_Clock/D");
        tree[i]->Branch("Gate",            &_Gate,                 "_Gate/D");
    }

    UInt_t clock_info[N_MAX_CHIP][N_MAX_CLOCKS+1];
    UInt_t clock_info_0[N_MAX_CLOCKS+1][N_MAX_CLUSTERS_2];
    UInt_t clock_info_1[N_MAX_CLOCKS+1][N_MAX_CLUSTERS_2];
    UInt_t clock_info_2[N_MAX_CLOCKS+1][N_MAX_CLUSTERS_2];
    UInt_t clock_info_3[N_MAX_CLOCKS+1][N_MAX_CLUSTERS_2];

    Int_t* cl_size          = new Int_t[N_MAX_CLUSTERS_1];
    Double_t* cl_pos_x      = new Double_t[N_MAX_CLUSTERS_1];
    Double_t* cl_pos_x_err  = new Double_t[N_MAX_CLUSTERS_1];
    Double_t* cl_pos_y      = new Double_t[N_MAX_CLUSTERS_1];
    Double_t* cl_pos_y_err  = new Double_t[N_MAX_CLUSTERS_1];
    Double_t* cl_clocks     = new Double_t[N_MAX_CLUSTERS_1];

    TString clock_info_cs     = "clock_info["; clock_info_cs   += N_MAX_CHIP;     clock_info_cs   += "]["; clock_info_cs   += N_MAX_CLOCKS+1; clock_info_cs     += "]/i";
    TString clock_info_0_cs   = "clock_info_0["; clock_info_0_cs += N_MAX_CLOCKS+1; clock_info_0_cs += "]["; clock_info_0_cs += N_MAX_CLUSTERS_2;  clock_info_0_cs   += "]/i";
    TString clock_info_1_cs   = "clock_info_1["; clock_info_1_cs += N_MAX_CLOCKS+1; clock_info_1_cs += "]["; clock_info_1_cs += N_MAX_CLUSTERS_2;  clock_info_1_cs   += "]/i";
    TString clock_info_2_cs   = "clock_info_2["; clock_info_2_cs += N_MAX_CLOCKS+1; clock_info_2_cs += "]["; clock_info_2_cs += N_MAX_CLUSTERS_2;  clock_info_2_cs   += "]/i";
    TString clock_info_3_cs   = "clock_info_3["; clock_info_3_cs += N_MAX_CLOCKS+1; clock_info_3_cs += "]["; clock_info_3_cs += N_MAX_CLUSTERS_2;  clock_info_3_cs   += "]/i";

    TTree* tree_clock = new TTree("Tree_Clock", "Cluster data for Clocks");
    Int_t clock_tic = 0;
    tree_clock->Branch("clock_info",clock_info,clock_info_cs.Data());
    tree_clock->Branch("clock_info_0",clock_info_0,clock_info_0_cs.Data());
    tree_clock->Branch("clock_info_1",clock_info_1,clock_info_1_cs.Data());
    tree_clock->Branch("clock_info_2",clock_info_2,clock_info_2_cs.Data());
    tree_clock->Branch("clock_info_3",clock_info_3,clock_info_3_cs.Data());
    //------------------------------------------------------------------------------//
    //------------------------------- HISTOGRAMS -----------------------------------//
    //------------------------------------------------------------------------------//

    TString hist_name;

    TH2D* h_1[N_MAX_CHIP];

    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        hist_name = "h_1_"; hist_name += i;
        h_1[i] = new TH2D(hist_name.Data(),hist_name.Data(),N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    }
    //------------------------------------------------------------------------------//
    Double_t zero_time = -999.999, event_time = -999.999, last_time = -999.999;
    Long64_t nFrames[N_MAX_CHIP] = {};

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------- MAIN LOOP ------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    Long64_t nEntriesCommon = 0;
    for(Int_t i = 1; i < N_MAX_CHIP; i++)
    {
        if(nEntries[0] != nEntries[i])
        {
            cout<<"--> ERROR:: Different number of entries for the chips!!!"<<endl;
            assert(0);
        }
    }

    UInt_t** yyy = new UInt_t*[N_MAX_CHIP];
    for(Int_t yyy_id1 = 0; yyy_id1 < N_MAX_CHIP; yyy_id1++)
    {
        yyy[yyy_id1] = new UInt_t[N_MAX_CLOCKS+1];
    }

    nEntriesCommon = nEntries[0];
    cout<<endl<<"--> nEntriesCommon: "<<nEntriesCommon<<endl;
    UInt_t nWarnings = 0;

    for(Int_t i = 0; i < nEntriesCommon; i++)
    {
        if(i%1 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntriesCommon);
            fflush(stdout);
        }

        for(Int_t chipID = 0; chipID < N_MAX_CHIP; chipID++)
        {
            for(clock_tic = 0; clock_tic < N_MAX_CLOCKS+1; clock_tic++)
            {
                clock_info[chipID][clock_tic] = 0;
            }

            fChain[chipID]->GetEntry(i);

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

            Long64_t frame_size = 0;
            for(Int_t xi = 0; xi < N_PIXELS; xi++)
            {
                for(Int_t yi = 0; yi < N_PIXELS; yi++)
                {
                    if(_COUNTS[xi][yi] > 0)
                    {
                        if(_AcquisType == 0)
                        {
                            frame_size += _COUNTS[xi][yi];

                            h_1[chipID]->Fill(xi,yi,_COUNTS[xi][yi]);
                        }
                        else if(_AcquisType == 3)
                        {
                            frame_size++;

                            h_1[chipID]->Fill(xi,yi,1);
                        }
                        else
                        {
                            cout<<endl<<endl<<"--> ERROR:: _AcquisType is not 0 or 3"<<endl<<endl;
                            assert(0);
                        }
                    }
                }
            }
            //---------------------------------------------------------------------------------------------------------------//
            //------------------------------------------ CLUSTER ANALYSIS ---------------------------------------------------//
            //---------------------------------------------------------------------------------------------------------------//
            if(_cluster_analysis)
            {
                for(Int_t jk = 0; jk < N_MAX_CLUSTERS_1; jk++)
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

                if(cl_num >= N_MAX_CLUSTERS_1)
                {
                    cout<<endl<<"--> ERROR:: The maximum number of clusters per a frame was reached !!!"<<endl;
                    assert(0);
                }

                for(Int_t yyy_id1 = 0; yyy_id1 < N_MAX_CHIP; yyy_id1++)
                {
                    for(Int_t yyy_id2 = 0; yyy_id2 < N_MAX_CLOCKS+1; yyy_id2++)
                    {
                        yyy[yyy_id1][yyy_id2] = 0;
                    }
                }

                for(Int_t yy = 0; yy < cl_num; yy++)
                {
                    //------------------------------------------------------------------------------//
                    //------------------------ For output file with clusters infor -----------------//
                    //------------------------------------------------------------------------------//
                    unix_time_clinfo    = event_time/1000.0;
                    clocks_clinfo       = (Double_t)N_MAX_CLOCKS - cl_clocks[yy];
                    if(clocks_clinfo < 0) clocks_clinfo = 0;
                    size_clinfo         = cl_size[yy];
                    pos_x_clinfo        = cl_pos_x[yy];
                    pos_x_err_clinfo    = cl_pos_x_err[yy];
                    pos_y_clinfo        = cl_pos_y[yy];
                    pos_y_err_clinfo    = cl_pos_y_err[yy];
                    event_id_clinfo     = _event;
                    if(_event != i) {event_id_clinfo = i;}

                    tree[chipID]->Fill();

                    clock_tic = round(clocks_clinfo);
                    clock_info[chipID][clock_tic]++;

                    if(clock_tic > N_MAX_CLOCKS) clock_tic = N_MAX_CLOCKS;

                    if(yyy[chipID][clock_tic] >= N_MAX_CLUSTERS_2)
                    {
                        cout<<endl;
                        cout<<"--> EventID: "<<i<<endl;
                        cout<<"--> chipID: "<<chipID<<" | clock_tic: "<<clock_tic<<endl;
                        cout<<endl<<"--> WARNING:: The maximum number of clusters per clock tic was reached !!!"<<endl;
                        nWarnings++;
                    }
                    else
                    {
                        switch (chipID)
                        {
                        case 0:
                            clock_info_0[clock_tic][yyy[chipID][clock_tic]] = cluster_id[chipID];
                            break;
                        case 1:
                            clock_info_1[clock_tic][yyy[chipID][clock_tic]] = cluster_id[chipID];
                            break;
                        case 2:
                            clock_info_2[clock_tic][yyy[chipID][clock_tic]] = cluster_id[chipID];
                            break;
                        case 3:
                            clock_info_3[clock_tic][yyy[chipID][clock_tic]] = cluster_id[chipID];
                            break;
                        default:
                            cout<<endl<<"--> ERROR:: Wrong chipID !!!"<<endl;
                            assert(0);
                            break;
                        }
                        yyy[chipID][clock_tic]++;
                    }

                    cluster_id[chipID]++;
                    //------------------------------------------------------------------------------//
                }
            }            
            //---------------------------------------------------------------------------------------------------------------//
            nFrames[chipID]++;
        }
        tree_clock->Fill();
    }
    cout<<endl;
    cout<<endl<<"--> Number of Warnings: "<<nWarnings<<endl<<endl;

    for(Int_t yyy_id1 = 0; yyy_id1 < N_MAX_CHIP; yyy_id1++)
    {
        delete[] yyy[yyy_id1];
    }
    delete[] yyy;

    delete [] cl_size;
    delete [] cl_pos_x;
    delete [] cl_pos_x_err;
    delete [] cl_pos_y;
    delete [] cl_pos_y_err;
    delete [] cl_clocks;

    cout<<endl;

    Double_t scfactor;
    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        scfactor = _Gate*nFrames[i];
        cout<<"--> nFrames["<<i<<"] = "<<nFrames[i]<<endl;
        cout<<"--> scaling factor = "<<scfactor<<" [sec]"<<endl;
    }

    string zero_time_char;
    string last_time_char;
    epochtime2date(round(zero_time/1000),zero_time_char);
    epochtime2date(round(last_time/1000),last_time_char);

    cout<<endl;
    cout<<"--> Initial time: "<<zero_time_char<<" | Final time: "<<last_time_char<<endl;

    //------------------------------------------------------------------------------//
    //------------------------ For output file with clusters infor -----------------//
    //------------------------------------------------------------------------------//
    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        tree[i]->Write();
        cout<<"--> The tree["<<i<<"] with "<<tree[i]->GetEntriesFast()<<" entries is written to the file: "<<outputfile->GetName()<<endl;
    }
    tree_clock->Write();
    cout<<"--> The tree_clock with "<<tree_clock->GetEntriesFast()<<" entries is written to the file: "<<outputfile->GetName()<<endl;
    //------------------------------------------------------------------------------//

    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        h_1[i]->GetXaxis()->SetTitle("X [pixels]");
        h_1[i]->GetYaxis()->SetTitle("Y [pixels]");


        h_1[i]->Write();
    }

    outputfile->Close();    
    cout<<"--> The output file '"<<outputfile->GetName()<<"' is closed."<<endl;
    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        fChain[i]->Delete();
    }

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
    Double_t mean_x, mean_y;

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
                case 3:
                    checkregionTOA(_matrix,fired_matrix,xi,yi,Clock);
                    break;
                default:
                    cout<<endl<<endl<<"--> ERROR:: _AcquisType < 0"<<endl<<endl;
                    assert(0);
                }

                // mean value
                mean_x = 0;
                mean_y = 0;

                for(Int_t yj = 0; yj < N_PIXELS; yj++)
                {
                    for(Int_t xj = 0; xj < N_PIXELS; xj++)
                    {
                        if(fired_matrix[xj][yj] > 0)
                        {
                            cluster_size[cluster_num]++;

                            mean_x += xj;
                            mean_y += yj;

                            fired_matrix_full[xj][yj] = 1;
                            cluster_clocks[cluster_num] += _matrix[xj][yj];
                        }
                        fired_matrix[xj][yj] = 0;
                    }
                }

                cluster_clocks[cluster_num]     = (Double_t)cluster_clocks[cluster_num]/cluster_size[cluster_num];
                cluster_pos_x[cluster_num]      = (Double_t)mean_x/cluster_size[cluster_num];
                cluster_pos_x_err[cluster_num]  = (Double_t)(1.0/TMath::Sqrt(12.0))/cluster_size[cluster_num];
                cluster_pos_y[cluster_num]      = (Double_t)mean_y/cluster_size[cluster_num];
                cluster_pos_y_err[cluster_num]  = (Double_t)(1.0/TMath::Sqrt(12.0))/cluster_size[cluster_num];

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
