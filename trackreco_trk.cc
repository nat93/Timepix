//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     Track reconstruction at H8 for 4 planes.
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
#include "TProfile.h"
//C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/*
 *  devCHIP_2=plane_0		#I04
 *  devCHIP_0=plane_1		#J05/L07
 *  devCHIP_1=plane_2		#B04
 *  devCHIP_3=plane_3		#B05
 */

//TOA= (N_MAX_CLOCKS-_COUNTS)*1e-6/_Clock //[sec]

const Int_t N_PIXELS        = 512;      // The maximum number of pixels per axis per chip
const Int_t N_MAX_CHIP      = 4;        // The maximum number of chips
const Int_t N_MAX_CLOCKS    = 11810;    // The maximum number of counts in the pixel
const Int_t N_MAX_CLUSTERS  = 30;       // The maximum number of clusters per clock
const Double_t PIXEL_SIZE   = 0.055;    // Pixel size [mm]

void epochtime2date(time_t in_time, string &out_time);
int GetPlaneID(Int_t chipID);
void GetRotaion(Int_t planeID, Double_t &x_pos, Double_t &y_pos, Double_t &x_pos_err, Double_t &y_pos_err);

int main(int argc, char *argv[])
{
    //------------------------------------------------------------------------------//
    //---------------------- Parameters for track reconstructions ------------------//
    //------------------------------------------------------------------------------//

    Int_t   nIterations                     = 10;
    Double_t offset_x_ChipID[N_MAX_CHIP], offset_errx_ChipID[N_MAX_CHIP], offset_y_ChipID[N_MAX_CHIP], offset_erry_ChipID[N_MAX_CHIP];

    Double_t offset_xres_Plane[N_MAX_CHIP]      = {};
    Double_t offset_errxres_Plane[N_MAX_CHIP]   = {};
    Double_t offset_yres_Plane[N_MAX_CHIP]      = {};
    Double_t offset_erryres_Plane[N_MAX_CHIP]   = {};

    //---------------------------------------------//
    //---------- Residual corrections -------------//
    //---------------------------------------------//
    // After 10 iteration
    offset_xres_Plane[0]    = -0.213174;
    offset_errxres_Plane[0] =  0.001565;
    offset_yres_Plane[0]    = -0.149726;
    offset_erryres_Plane[0] =  0.001886;

    offset_xres_Plane[1]    =  0.290287;
    offset_errxres_Plane[1] =  0.001539;
    offset_yres_Plane[1]    =  0.213782;
    offset_erryres_Plane[1] =  0.001902;

    offset_xres_Plane[2]    = -0.060535;
    offset_errxres_Plane[2] =  0.001830;
    offset_yres_Plane[2]    =  0.012324;
    offset_erryres_Plane[2] =  0.001840;

    offset_xres_Plane[3]    = -0.072363;
    offset_errxres_Plane[3] =  0.001367;
    offset_yres_Plane[3]    = -0.100682;
    offset_erryres_Plane[3] =  0.001518;
    //---------------------------------------------//

    Short_t clock_jitter                    = 2;                                            // [clocks]
    Double_t plane_position[]               = {0.0,305.0,305.0+600.0,305.0+600.0+270.0};    // [mm]
    Double_t plane_position_err[]           = {10.0,10.0,10.0,10.0};                        // [mm]
    Double_t timecut                        = 10000.0;                                      // [ns]
    Double_t projection_plane_position      = 310.0;                                        // [mm]
    Double_t projection_plane_position_err  = 10.0;                                         // [mm]
    Double_t chi2_x_max                     = 1000;                                         // [chi2/NDF]
    Double_t chi2_y_max                     = 1000;                                         // [chi2/NDF]


    //------------------------------------------------------------------------------//
    //--------------------------- Input parameters/variables -----------------------//
    //------------------------------------------------------------------------------//
    cout<<endl<<"--> TRACK RECONSTRUCTION"<<endl<<endl;
    if(argc != 4)
    {
        cout<<"--> ERROR: Wrong number of input parameters("<<argc<<")!"<<endl<<endl;

        cout<<"--> [0] ./script_name"<<endl;
        cout<<"--> [1] Input file name"<<endl;
        cout<<"--> [2] Output file name"<<endl;
        cout<<"--> [3] Run mode"<<endl;

        return -1;
    }

    TString inFileName  = argv[1];
    TString outFileName = argv[2];
    Int_t Run_Mode   = atoi(argv[3]);

    time_t start_time, stop_time;
    start_time = time(NULL);

    //------------------------------------------------------------------------------//
    //------------------------ For input file with clusters infor ------------------//
    //------------------------------------------------------------------------------//
    Double_t unix_clocks_clinfo, pos_x_clinfo, pos_x_err_clinfo, pos_y_clinfo, pos_y_err_clinfo, _Clock, clocks_clinfo, _Gate;
    Int_t size_clinfo, event_id_clinfo;    

    TString tree_name;
    TChain* fChain[N_MAX_CHIP];
    Long64_t nEntries[N_MAX_CHIP];

    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        tree_name = "Tree_";
        tree_name += i;
        fChain[i] = new TChain(tree_name.Data());
        fChain[i]->Add(inFileName.Data());

        fChain[i]->SetBranchAddress("UnixTime",        &unix_clocks_clinfo);
        fChain[i]->SetBranchAddress("ClusterClocks",   &clocks_clinfo);
        fChain[i]->SetBranchAddress("ClusterSize",     &size_clinfo);
        fChain[i]->SetBranchAddress("ClusterPosX",     &pos_x_clinfo);
        fChain[i]->SetBranchAddress("ClusterPosXerr",  &pos_x_err_clinfo);
        fChain[i]->SetBranchAddress("ClusterPosY",     &pos_y_clinfo);
        fChain[i]->SetBranchAddress("ClusterPosYerr",  &pos_y_err_clinfo);
        fChain[i]->SetBranchAddress("EventID",         &event_id_clinfo);
        fChain[i]->SetBranchAddress("Clock",           &_Clock);
        fChain[i]->SetBranchAddress("Gate",            &_Gate);
    }
    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        fChain[i]->GetEntry(0);
        nEntries[i] = (fChain[i]->GetEntries());
        cout<<"--> ChipID: "<<i<<endl;
        cout<<"--> InputFileName: "<<inFileName<<endl;
        cout<<"--> Number of nEntries: "<<nEntries[i]<<endl<<endl;
    }

    UInt_t clock_info[N_MAX_CHIP][N_MAX_CLOCKS+1];
    UInt_t clock_info_0[N_MAX_CLOCKS+1][N_MAX_CLUSTERS];
    UInt_t clock_info_1[N_MAX_CLOCKS+1][N_MAX_CLUSTERS];
    UInt_t clock_info_2[N_MAX_CLOCKS+1][N_MAX_CLUSTERS];
    UInt_t clock_info_3[N_MAX_CLOCKS+1][N_MAX_CLUSTERS];

    Long64_t nEntriesClock;
    TChain* fChain_clock = new TChain("Tree_Clock");
    fChain_clock->Add(inFileName.Data());

    fChain_clock->SetBranchAddress("clock_info",clock_info);
    fChain_clock->SetBranchAddress("clock_info_0",clock_info_0);
    fChain_clock->SetBranchAddress("clock_info_1",clock_info_1);
    fChain_clock->SetBranchAddress("clock_info_2",clock_info_2);
    fChain_clock->SetBranchAddress("clock_info_3",clock_info_3);

    nEntriesClock = fChain_clock->GetEntries();
    cout<<"--> A Tree with ClockInfo"<<endl;
    cout<<"--> Number of nEntries: "<<nEntriesClock<<endl<<endl;
    //------------------------------------------------------------------------------//
    //------------------------- For output file ------------------------------------//
    //------------------------------------------------------------------------------//
    cout<<"--> OutputFileName: "<<outFileName<<endl<<endl;
    TFile* outputfile = new TFile(outFileName.Data(),"RECREATE");
    if (outputfile->IsZombie())
    {
        cout<<"--> ERROR: The output file is zombie!"<<endl;
        assert(0);
    }
    //------------------------------------------------------------------------------//
    //------------------------------- HISTOGRAMS -----------------------------------//
    //------------------------------------------------------------------------------//

    //---------//
    // Reco 0
    //---------//

    TH1D* h1 = new TH1D("h1","#DeltaTime Chip0-Chip2",2*N_MAX_CLOCKS-1,-_Gate*1e9,_Gate*1e9);
    TH1D* h2 = new TH1D("h2","#DeltaTime Chip0-Chip1",2*N_MAX_CLOCKS-1,-_Gate*1e9,_Gate*1e9);
    TH1D* h3 = new TH1D("h3","#DeltaTime Chip0-Chip3",2*N_MAX_CLOCKS-1,-_Gate*1e9,_Gate*1e9);

    TH1D* h4 = new TH1D("h4","X position without appl. offset Chip_0",2000,-1000,1000);
    TH1D* h5 = new TH1D("h5","X position with appl. offset Chip_0",2000,-1000,1000);
    TH1D* h6 = new TH1D("h6","X position without appl. offset Chip_1",2000,-1000,1000);
    TH1D* h7 = new TH1D("h7","X position with appl. offset Chip_1",2000,-1000,1000);
    TH1D* h8 = new TH1D("h8","X position without appl. offset Chip_2",2000,-1000,1000);
    TH1D* h9 = new TH1D("h9","X position with appl. offset Chip_2",2000,-1000,1000);
    TH1D* h10 = new TH1D("h10","X position without appl. offset Chip_3",2000,-1000,1000);
    TH1D* h11 = new TH1D("h11","X position with appl. offset Chip_3",2000,-1000,1000);

    TH1D* h12 = new TH1D("h12","Y position without appl. offset Chip_0",2000,-1000,1000);
    TH1D* h13 = new TH1D("h13","Y position with appl. offset Chip_0",2000,-1000,1000);
    TH1D* h14 = new TH1D("h14","Y position without appl. offset Chip_1",2000,-1000,1000);
    TH1D* h15 = new TH1D("h15","Y position with appl. offset Chip_1",2000,-1000,1000);
    TH1D* h16 = new TH1D("h16","Y position without appl. offset Chip_2",2000,-1000,1000);
    TH1D* h17 = new TH1D("h17","Y position with appl. offset Chip_2",2000,-1000,1000);
    TH1D* h18 = new TH1D("h18","Y position without appl. offset Chip_3",2000,-1000,1000);
    TH1D* h19 = new TH1D("h19","Y position with appl. offset Chip_3",2000,-1000,1000);

    TH2D* h21 = new TH2D("h21","X position vs Z position",160,-100,1500,2000,-1000*PIXEL_SIZE,1000*PIXEL_SIZE);
    TH2D* h22 = new TH2D("h22","Y position vs Z position",160,-100,1500,2000,-1000*PIXEL_SIZE,1000*PIXEL_SIZE);

    //---------//
    // Reco 1
    //---------//

    TH2D* h25 = new TH2D("h25","#Chi^{2} in X direction vs Iteration",nIterations+1,0,nIterations+1,100000,0,10000);
    TH2D* h26 = new TH2D("h26","#Chi^{2}/NDF in X direction vs Iteration",nIterations+1,0,nIterations+1,100000,0,10000);
    TH2D* h27 = new TH2D("h27","#Chi^{2} in Y direction vs Iteration",nIterations+1,0,nIterations+1,100000,0,10000);
    TH2D* h28 = new TH2D("h28","#Chi^{2}/NDF in Y direction vs Iteration",nIterations+1,0,nIterations+1,100000,0,10000);

    TH1D* h29 = new TH1D("h29","Residual in X for Plane 0",10000,-10,10);
    TH1D* h30 = new TH1D("h30","Residual in X for Plane 1",10000,-10,10);
    TH1D* h31 = new TH1D("h31","Residual in X for Plane 2",10000,-10,10);
    TH1D* h32 = new TH1D("h32","Residual in X for Plane 3",10000,-10,10);

    TH1D* h33 = new TH1D("h33","Residual in Y for Plane 0",10000,-10,10);
    TH1D* h34 = new TH1D("h34","Residual in Y for Plane 1",10000,-10,10);
    TH1D* h35 = new TH1D("h35","Residual in Y for Plane 2",10000,-10,10);
    TH1D* h36 = new TH1D("h36","Residual in Y for Plane 3",10000,-10,10);

    TH1D* h37 = new TH1D("h37","#Chi^{2} in X direction (with Residual corrections)",100000,0,10000);
    TH1D* h38 = new TH1D("h38","#Chi^{2}/NDF in X direction (with Residual corrections)",100000,0,10000);
    TH1D* h39 = new TH1D("h39","#Chi^{2} in Y direction (with Residual corrections)",100000,0,10000);
    TH1D* h40 = new TH1D("h40","#Chi^{2}/NDF in Y direction (with Residual corrections)",100000,0,10000);
    TH1D* h41 = new TH1D("h41","ToA of clusters for Chip_0",11811,0,11811);
    TH1D* h42 = new TH1D("h42","ToA of clusters for Chip_1",11811,0,11811);
    TH1D* h43 = new TH1D("h43","ToA of clusters for Chip_2",11811,0,11811);
    TH1D* h44 = new TH1D("h44","ToA of clusters for Chip_3",11811,0,11811);
    TH1D* h45 = new TH1D("h45","Size of clusters for Chip_0",1000,0,100);
    TH1D* h46 = new TH1D("h46","Size of clusters for Chip_1",1000,0,100);
    TH1D* h47 = new TH1D("h47","Size of clusters for Chip_2",1000,0,100);
    TH1D* h48 = new TH1D("h48","Size of clusters for Chip_3",1000,0,100);

    TH2D* h49 = new TH2D("h49","Residual in X vs Iteraction for Plane 0",nIterations+1,0,nIterations+1,10000,-10,10);
    TH2D* h50 = new TH2D("h50","Residual in X vs Iteraction for Plane 1",nIterations+1,0,nIterations+1,10000,-10,10);
    TH2D* h51 = new TH2D("h51","Residual in X vs Iteraction for Plane 2",nIterations+1,0,nIterations+1,10000,-10,10);
    TH2D* h52 = new TH2D("h52","Residual in X vs Iteraction for Plane 3",nIterations+1,0,nIterations+1,10000,-10,10);

    TH2D* h120 = new TH2D("h120","Residual in Y vs Iteraction for Plane 0",nIterations+1,0,nIterations+1,10000,-10,10);
    TH2D* h121 = new TH2D("h121","Residual in Y vs Iteraction for Plane 1",nIterations+1,0,nIterations+1,10000,-10,10);
    TH2D* h122 = new TH2D("h122","Residual in Y vs Iteraction for Plane 2",nIterations+1,0,nIterations+1,10000,-10,10);
    TH2D* h123 = new TH2D("h123","Residual in Y vs Iteraction for Plane 3",nIterations+1,0,nIterations+1,10000,-10,10);

    TH2D* h124 = new TH2D("h124","Residual in X vs Interpolated Y",10000,-50,50,2000,-10,10);
    TH2D* h125 = new TH2D("h125","Residual in X vs Interpolated Y",10000,-50,50,2000,-10,10);
    TH2D* h126 = new TH2D("h126","Residual in X vs Interpolated Y",10000,-50,50,2000,-10,10);
    TH2D* h127 = new TH2D("h127","Residual in X vs Interpolated Y",10000,-50,50,2000,-10,10);

    TH2D* h128 = new TH2D("h128","Residual in Y vs Interpolated X",10000,-50,50,2000,-10,10);
    TH2D* h129 = new TH2D("h129","Residual in Y vs Interpolated X",10000,-50,50,2000,-10,10);
    TH2D* h130 = new TH2D("h130","Residual in Y vs Interpolated X",10000,-50,50,2000,-10,10);
    TH2D* h131 = new TH2D("h131","Residual in Y vs Interpolated X",10000,-50,50,2000,-10,10);

    //---------//
    // Reco 2
    //---------//

    TH1D* h53 = new TH1D("h53","Track Angle in X Planes 0-1 (II)",100000,-1,1);
    TH1D* h54 = new TH1D("h54","Track Angle in Y Planes 0-1 (II)",100000,-1,1);
    TH1D* h55 = new TH1D("h55","Track Angle in X Planes 2-3 (II)",100000,-1,1);
    TH1D* h56 = new TH1D("h56","Track Angle in Y Planes 2-3 (II)",100000,-1,1);
    TH1D* h57 = new TH1D("h57","dX",10000,-100,100);
    TH1D* h58 = new TH1D("h58","dY",10000,-100,100);
    TH1D* h59 = new TH1D("h59","X Arm 1",10000,-100,100);
    TH1D* h60 = new TH1D("h60","X Arm 2",10000,-100,100);
    TH1D* h61 = new TH1D("h61","Y Arm 1",10000,-100,100);
    TH1D* h62 = new TH1D("h62","Y Arm 2",10000,-100,100);
    TH1D* h63 = new TH1D("h63","#Chi^{2}_{X}",1100,-100,1000);
    TH1D* h64 = new TH1D("h64","#Chi^{2}_{Y}",1100,-100,1000);
    TH2D* h65 = new TH2D("h65","#Chi^{2}_{X} vs dX",1000,-100,100,1100,-100,1000);
    TH2D* h66 = new TH2D("h66","#Chi^{2}_{Y} vs dY",1000,-100,100,1100,-100,1000);
    TH2D* h67 = new TH2D("h67","#Delta_{X} Plane3-0 vs X position on Plane 0",400,-20,20,400,-20,20);
    TH2D* h68 = new TH2D("h68","#Delta_{X} Plane3-0 vs Y position on Plane 0",400,-20,20,400,-20,20);
    TH2D* h69 = new TH2D("h69","#Delta_{Y} Plane3-0 vs X position on Plane 0",400,-20,20,400,-20,20);
    TH2D* h70 = new TH2D("h70","#Delta_{Y} Plane3-0 vs Y position on Plane 0",400,-20,20,400,-20,20);
    TH2D* h71 = new TH2D("h71","#Delta_{X} Plane3-1 vs X position on Plane 1",400,-20,20,400,-20,20);
    TH2D* h72 = new TH2D("h72","#Delta_{X} Plane3-1 vs Y position on Plane 1",400,-20,20,400,-20,20);
    TH2D* h73 = new TH2D("h73","#Delta_{Y} Plane3-1 vs X position on Plane 1",400,-20,20,400,-20,20);
    TH2D* h74 = new TH2D("h74","#Delta_{Y} Plane3-1 vs Y position on Plane 1",400,-20,20,400,-20,20);
    TH2D* h75 = new TH2D("h75","#Delta_{X} Plane3-2 vs X position on Plane 2",400,-20,20,400,-20,20);
    TH2D* h76 = new TH2D("h76","#Delta_{X} Plane3-2 vs Y position on Plane 2",400,-20,20,400,-20,20);
    TH2D* h77 = new TH2D("h77","#Delta_{Y} Plane3-2 vs X position on Plane 2",400,-20,20,400,-20,20);
    TH2D* h78 = new TH2D("h78","#Delta_{Y} Plane3-2 vs Y position on Plane 2",400,-20,20,400,-20,20);
    TH2D* h79 = new TH2D("h79","Y position on Plane 3 vs Y position on Plane 0",400,-20,20,400,-20,20);
    TH2D* h80 = new TH2D("h80","Y position on Plane 3 vs Y position on Plane 1",400,-20,20,400,-20,20);
    TH2D* h81 = new TH2D("h81","Y position on Plane 3 vs Y position on Plane 2",400,-20,20,400,-20,20);

    TH1D* h83 = new TH1D("h83","#Theta_{X}^{Arm 1} - #Theta_{X}^{Arm 2}",2000,-1000,1000);
    TH1D* h84 = new TH1D("h84","#Theta_{Y}^{Arm 1} - #Theta_{Y}^{Arm 2}",2000,-1000,1000);

    TH2D* h85 = new TH2D("h85","#Delta#Theta_{X} vs X position",200,-10,10,100000,-0.1,0.1);
    TProfile *h86 = new TProfile("h86","Profile: #Delta_{X} Plane3-0 vs Y position on Plane 0",400,-20,20,-20,20);
    TProfile *h87 = new TProfile("h87","Profile: #Delta_{Y} Plane3-0 vs X position on Plane 0",400,-20,20,-20,20);
    TProfile *h88 = new TProfile("h88","Profile: #Delta_{X} Plane3-1 vs Y position on Plane 1",400,-20,20,-20,20);
    TProfile *h89 = new TProfile("h89","Profile: #Delta_{Y} Plane3-1 vs X position on Plane 1",400,-20,20,-20,20);
    TProfile *h90 = new TProfile("h90","Profile: #Delta_{X} Plane3-2 vs Y position on Plane 2",400,-20,20,-20,20);
    TProfile *h91 = new TProfile("h91","Profile: #Delta_{Y} Plane3-2 vs X position on Plane 2",400,-20,20,-20,20);

    TH1D* h92 = new TH1D("h92","#Phi on Plane 0",2000,-2.0,2.0);
    TH1D* h93 = new TH1D("h93","#Phi on Plane 1",2000,-2.0,2.0);
    TH1D* h94 = new TH1D("h94","#Phi on Plane 2",2000,-2.0,2.0);
    TH1D* h95 = new TH1D("h95","#Phi on Plane 3",2000,-2.0,2.0);

    TH1D* h96 = new TH1D("h96","#Delta#Phi Plane 0 - Plane 3",4000,-4.0,4.0);
    TH1D* h97 = new TH1D("h97","#Delta#Phi Plane 1 - Plane 3",4000,-4.0,4.0);
    TH1D* h98 = new TH1D("h98","#Delta#Phi Plane 2 - Plane 3",4000,-4.0,4.0);
    TH1D* h99 = new TH1D("h99","X coord Plane 0",400,-20.0,20.0);
    TH1D* h100 = new TH1D("h100","X coord Plane 1",400,-20.0,20.0);
    TH1D* h101 = new TH1D("h101","X coord Plane 2",400,-20.0,20.0);
    TH1D* h102 = new TH1D("h102","X coord Plane 3",400,-20.0,20.0);
    TH1D* h103 = new TH1D("h103","Y coord Plane 0",400,-20.0,20.0);
    TH1D* h104 = new TH1D("h104","Y coord Plane 1",400,-20.0,20.0);
    TH1D* h105 = new TH1D("h105","Y coord Plane 2",400,-20.0,20.0);
    TH1D* h106 = new TH1D("h106","Y coord Plane 3",400,-20.0,20.0);
    TH2D* h107 = new TH2D("h107","X Plane 0 vs X Plane 3",400,-20.0,20.0,400,-20.0,20.0);
    TH2D* h108 = new TH2D("h108","X Plane 1 vs X Plane 3",400,-20.0,20.0,400,-20.0,20.0);
    TH2D* h109 = new TH2D("h109","X Plane 2 vs X Plane 3",400,-20.0,20.0,400,-20.0,20.0);
    TH2D* h110 = new TH2D("h110","Y vs X Plane 0",400,-20,20,400,-20,20);
    TH2D* h111 = new TH2D("h111","Y vs X Plane 1",400,-20,20,400,-20,20);
    TH2D* h112 = new TH2D("h112","Y vs X Plane 2",400,-20,20,400,-20,20);
    TH2D* h113 = new TH2D("h113","Y vs X Plane 3",400,-20,20,400,-20,20);
    TH2D* h114 = new TH2D("h114","Y vs X Plane 0 (cut on #Delta#Theta > 0.4 mrad))",400,-20,20,400,-20,20);
    TH2D* h115 = new TH2D("h115","Y vs X Plane 0 (cut on #Delta#Theta < 0.4 mrad)",400,-20,20,400,-20,20);
    TH2D* h116 = new TH2D("h116","#Delta#Theta_{Y} vs X position",200,-10,10,100000,-0.1,0.1);
    TH2D* h117 = new TH2D("h117","#Delta#Theta_{Y} vs Y position",200,-10,10,100000,-0.1,0.1);
    TH2D* h118 = new TH2D("h118","Y vs X Plane 3 (cut on #Delta#Theta > 0.4 mrad))",400,-20,20,400,-20,20);
    TH2D* h119 = new TH2D("h119","Y vs X Plane 3 (cut on #Delta#Theta < 0.4 mrad)",400,-20,20,400,-20,20);
    //------------------------------------------------------------------------------------------------//
    //------------------------------- DELTA TIME DISTRIBUTION ----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    Int_t i_initial_0 = 0, i_initial_1 = 0, i_initial_2 = 0, i_initial_3 = 0;
    if(Run_Mode == 0)
    {
        cout<<endl<<"--> Delta time distribution <--"<<endl;

        for(Int_t i = i_initial_0; i < nEntries[0]; i++)
        {
            if(i%1 == 0)
            {
                printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntries[0]);
                fflush(stdout);
            }

            // Chip_0
            fChain[0]->GetEntry(i);
            if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

            Int_t event_ref = event_id_clinfo;
            Double_t clocks_clinfo_ref = clocks_clinfo;

            // Chip_2
            for(Int_t j = i_initial_2; j < nEntries[2]; j++)
            {
                fChain[2]->GetEntry(j);
                if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

                if(event_id_clinfo == event_ref)
                {
                    h1->Fill((clocks_clinfo - clocks_clinfo_ref)*1.0e3/_Clock);
                    i_initial_2 = j;
                }
                else if(event_id_clinfo > event_ref) break;
            }
            // Chip_1
            for(Int_t j = i_initial_1; j < nEntries[1]; j++)
            {
                fChain[1]->GetEntry(j);
                if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

                if(event_id_clinfo == event_ref)
                {
                    h2->Fill((clocks_clinfo - clocks_clinfo_ref)*1.0e3/_Clock);
                    i_initial_1 = j;
                }
                else if(event_id_clinfo > event_ref) break;
            }
            // Chip_3
            for(Int_t j = i_initial_3; j < nEntries[3]; j++)
            {
                fChain[3]->GetEntry(j);
                if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

                if(event_id_clinfo == event_ref)
                {
                    h3->Fill((clocks_clinfo - clocks_clinfo_ref)*1.0e3/_Clock);
                    i_initial_3 = j;
                }
                else if(event_id_clinfo > event_ref) break;
            }
        }
    }

    //------------------------------------------------------------------------------------------------//
    //------------------------ INTEGRATED IMAGE OF THE BEAM ------------------------------------------//
    //------------------------------------------------------------------------------------------------//

    cout<<endl<<"--> Integrating image of the beam <--"<<endl<<endl;
    for(Int_t i = 0; i < nEntries[0]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(1.1): %3.1f %%",100*(Double_t)i/nEntries[0]);
            fflush(stdout);
        }

        fChain[0]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = pos_x_clinfo;
        Double_t y_pos = pos_y_clinfo;
        //-------------------//

        h4->Fill(x_pos);
        h12->Fill(y_pos);

        h41->Fill(clocks_clinfo);
        h45->Fill(size_clinfo);
    }
    cout<<endl;
    for(Int_t i = 0; i < nEntries[0]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(1.2): %3.1f %%",100*(Double_t)i/nEntries[0]);
            fflush(stdout);
        }

        fChain[0]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = pos_x_clinfo;
        Double_t y_pos = pos_y_clinfo;
        //-------------------//

        offset_x_ChipID[0]      = h4->GetMean();
        offset_errx_ChipID[0]   = h4->GetMeanError();
        offset_y_ChipID[0]      = h12->GetMean();
        offset_erry_ChipID[0]   = h12->GetMeanError();

        h5->Fill(x_pos-offset_x_ChipID[0]);
        h13->Fill(y_pos-offset_y_ChipID[0]);

        h21->Fill(plane_position[GetPlaneID(0)],(x_pos-offset_x_ChipID[0])*PIXEL_SIZE);
        h22->Fill(plane_position[GetPlaneID(0)],(y_pos-offset_y_ChipID[0])*PIXEL_SIZE);

        h41->Fill(clocks_clinfo);
        h45->Fill(size_clinfo);
    }
    cout<<endl;
    for(Int_t i = 0; i < nEntries[1]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(2.1): %3.1f %%",100*(Double_t)i/nEntries[1]);
            fflush(stdout);
        }

        fChain[1]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = (-1.0)*pos_y_clinfo;
        Double_t y_pos = pos_x_clinfo;
        //-------------------//

        h6->Fill(x_pos);
        h14->Fill(y_pos);

        h42->Fill(clocks_clinfo);
        h46->Fill(size_clinfo);
    }
    cout<<endl;
    for(Int_t i = 0; i < nEntries[1]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(2.2): %3.1f %%",100*(Double_t)i/nEntries[1]);
            fflush(stdout);
        }

        fChain[1]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = (-1.0)*pos_y_clinfo;
        Double_t y_pos = pos_x_clinfo;
        //-------------------//

        offset_x_ChipID[1]      = h6->GetMean();
        offset_errx_ChipID[1]   = h6->GetMeanError();
        offset_y_ChipID[1]      = h14->GetMean();
        offset_erry_ChipID[1]   = h14->GetMeanError();

        h7->Fill(x_pos-offset_x_ChipID[1]);
        h15->Fill(y_pos-offset_y_ChipID[1]);

        h21->Fill(plane_position[GetPlaneID(1)],(x_pos-offset_x_ChipID[1])*PIXEL_SIZE);
        h22->Fill(plane_position[GetPlaneID(1)],(y_pos-offset_y_ChipID[1])*PIXEL_SIZE);

        h42->Fill(clocks_clinfo);
        h46->Fill(size_clinfo);
    }
    cout<<endl;
    for(Int_t i = 0; i < nEntries[2]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(3.1): %3.1f %%",100*(Double_t)i/nEntries[2]);
            fflush(stdout);
        }

        fChain[2]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = pos_x_clinfo;
        Double_t y_pos = pos_y_clinfo;
        //-------------------//

        h8->Fill(x_pos);
        h16->Fill(y_pos);

        h43->Fill(clocks_clinfo);
        h47->Fill(size_clinfo);
    }
    cout<<endl;
    for(Int_t i = 0; i < nEntries[2]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(3.2): %3.1f %%",100*(Double_t)i/nEntries[2]);
            fflush(stdout);
        }

        fChain[2]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = pos_x_clinfo;
        Double_t y_pos = pos_y_clinfo;
        //-------------------//

        offset_x_ChipID[2]      = h8->GetMean();
        offset_errx_ChipID[2]   = h8->GetMeanError();
        offset_y_ChipID[2]      = h16->GetMean();
        offset_erry_ChipID[2]   = h16->GetMeanError();

        h9->Fill(x_pos-offset_x_ChipID[2]);
        h17->Fill(y_pos-offset_y_ChipID[2]);

        h21->Fill(plane_position[GetPlaneID(2)],(x_pos-offset_x_ChipID[2])*PIXEL_SIZE);
        h22->Fill(plane_position[GetPlaneID(2)],(y_pos-offset_y_ChipID[2])*PIXEL_SIZE);

        h43->Fill(clocks_clinfo);
        h47->Fill(size_clinfo);
    }
    cout<<endl;
    for(Int_t i = 0; i < nEntries[3]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(4.1): %3.1f %%",100*(Double_t)i/nEntries[3]);
            fflush(stdout);
        }

        fChain[3]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = (-1.0)*pos_y_clinfo;
        Double_t y_pos = pos_x_clinfo;
        //-------------------//

        h10->Fill(x_pos);
        h18->Fill(y_pos);

        h44->Fill(clocks_clinfo);
        h48->Fill(size_clinfo);
    }
    cout<<endl;
    for(Int_t i = 0; i < nEntries[3]; i++)
    {
        if(i%100 == 0)
        {
            printf("\r--> Working(4.2): %3.1f %%",100*(Double_t)i/nEntries[3]);
            fflush(stdout);
        }

        fChain[3]->GetEntry(i);
        if(clocks_clinfo*1.0e3/_Clock < timecut) continue;

        //--- Orientation ---//
        Double_t x_pos = (-1.0)*pos_y_clinfo;
        Double_t y_pos = pos_x_clinfo;
        //-------------------//

        offset_x_ChipID[3]      = h10->GetMean();
        offset_errx_ChipID[3]   = h10->GetMeanError();
        offset_y_ChipID[3]      = h18->GetMean();
        offset_erry_ChipID[3]   = h18->GetMeanError();

        h11->Fill(x_pos-offset_x_ChipID[3]);
        h19->Fill(y_pos-offset_y_ChipID[3]);

        h21->Fill(plane_position[GetPlaneID(3)],(x_pos-offset_x_ChipID[3])*PIXEL_SIZE);
        h22->Fill(plane_position[GetPlaneID(3)],(y_pos-offset_y_ChipID[3])*PIXEL_SIZE);

        h44->Fill(clocks_clinfo);
        h48->Fill(size_clinfo);
    }
    cout<<endl;

    //------------------------------------------------------------------------------------------------//
    //--------------------------------- TRACKING SECTION ---------------------------------------------//
    //------------------------------------------------------------------------------------------------//

    TF1* fit_function = new TF1("fit_function","pol1");
    Double_t slope_arm1_x, slope_arm2_x, slope_arm1_y, slope_arm2_y;
    Double_t slope_arm1_x_err, slope_arm2_x_err, slope_arm1_y_err, slope_arm2_y_err/*, offset_arm1_x_err, offset_arm2_x_err, offset_arm1_y_err, offset_arm2_y_err*/;
    Double_t x_arm1, x_arm2, y_arm1, y_arm2, x_arm1_err, x_arm2_err, y_arm1_err, y_arm2_err;

    Double_t chi2_x, chi2_y;
    Double_t fit_slope_x  = 1.0e-3;
    Double_t fit_offset_x = 0.0;
    Double_t fit_slope_y  = 1.0e-3;
    Double_t fit_offset_y = 0.0;

    TGraphErrors* gr_x = new TGraphErrors();
    TGraphErrors* gr_y = new TGraphErrors();

    if(Run_Mode == 1)
    {
        cout<<endl<<"--> Alignment <--"<<endl;

        //-----------//
        // Iteration //
        //-----------//

        for(Int_t iteratorID = 1; iteratorID <= nIterations; iteratorID++)
        {
            //----------------//
            // Loop by frames //
            //----------------//
            cout<<endl;

            TH1D* h_res_plane0_x = new TH1D("h_res_plane0_x","Residual in X for Plane 0",10000,-10,10);
            TH1D* h_res_plane1_x = new TH1D("h_res_plane1_x","Residual in X for Plane 1",10000,-10,10);
            TH1D* h_res_plane2_x = new TH1D("h_res_plane2_x","Residual in X for Plane 2",10000,-10,10);
            TH1D* h_res_plane3_x = new TH1D("h_res_plane3_x","Residual in X for Plane 3",10000,-10,10);

            TH1D* h_res_plane0_y = new TH1D("h_res_plane0_y","Residual in Y for Plane 0",10000,-10,10);
            TH1D* h_res_plane1_y = new TH1D("h_res_plane1_y","Residual in Y for Plane 1",10000,-10,10);
            TH1D* h_res_plane2_y = new TH1D("h_res_plane2_y","Residual in Y for Plane 2",10000,-10,10);
            TH1D* h_res_plane3_y = new TH1D("h_res_plane3_y","Residual in Y for Plane 3",10000,-10,10);


            TH1D* h_fit_offset_x = new TH1D("h_fit_offset_x","Fit Offset in X",100000,-100,100);
            TH1D* h_fit_slope_x = new TH1D("h_fit_slope_x","Fit Slope in X",10000000,-10,10);
            TH1D* h_fit_offset_y = new TH1D("h_fit_offset_y","Fit Offset in Y",100000,-100,100);
            TH1D* h_fit_slope_y = new TH1D("h_fit_slope_y","Fit Slope in Y",10000000,-10,10);

            for(Int_t i = 0; i < nEntriesClock; i++)
            {
                if(i%10 == 0)
                {
                    printf("\r--> Iteration: %5d/%5d | Working: %3.1f %%",iteratorID,nIterations,100*(Double_t)i/nEntriesClock);
                    fflush(stdout);
                }

                fChain_clock->GetEntry(i);

                //----------------//
                // Loop by clocks //
                //----------------//
                for(Int_t clockID = 0; clockID < N_MAX_CLOCKS+1; clockID++)
                {
                    gr_x->Set(0);
                    gr_y->Set(0);

                    UInt_t gr_point = 0;
                    Int_t chip_ready[N_MAX_CHIP] = {};

                    if(clockID*1.0e3/_Clock < timecut) continue;

                    //----------------//
                    // Loop by chips  //
                    //----------------//
                    for(Int_t chipID_temp = 0; chipID_temp < N_MAX_CHIP; chipID_temp++)
                    {
                        if(clock_info[chipID_temp][clockID] > 0)
                        {
                            for(Int_t chipID = 0; chipID < N_MAX_CHIP; chipID++)
                            {
                                //-----------------------------//
                                // Loop by clocks inside range //
                                //-----------------------------//
                                for(Int_t kkl = clockID-clock_jitter; kkl <= clockID+clock_jitter; kkl++)
                                {
                                    if(kkl < 0 || kkl > N_MAX_CLOCKS) continue;

                                    //-------------------------------//
                                    // Loop by clusters inside clock //
                                    //-------------------------------//
                                    for(UInt_t kkll = 0; kkll < clock_info[chipID][kkl]; kkll++)
                                    {
                                        chip_ready[chipID]++;
                                        UInt_t clusterID = -1;
                                        switch (chipID)
                                        {
                                        case 0:
                                            clusterID = clock_info_0[kkl][kkll];
                                            break;
                                        case 1:
                                            clusterID = clock_info_1[kkl][kkll];
                                            break;
                                        case 2:
                                            clusterID = clock_info_2[kkl][kkll];
                                            break;
                                        case 3:
                                            clusterID = clock_info_3[kkl][kkll];
                                            break;
                                        }
                                        fChain[chipID]->GetEntry(clusterID);

                                        //--- Orientation ---//
                                        Double_t x_pos      = pos_x_clinfo;
                                        Double_t x_pos_err  = pos_x_err_clinfo;
                                        Double_t y_pos      = pos_y_clinfo;
                                        Double_t y_pos_err  = pos_y_err_clinfo;

                                        if(chipID == 1 || chipID == 3)
                                        {
                                            x_pos       = (-1.0)*pos_y_clinfo;
                                            x_pos_err   = pos_y_err_clinfo;
                                            y_pos       = pos_x_clinfo;
                                            y_pos_err   = pos_x_err_clinfo;
                                        }

                                        GetRotaion(GetPlaneID(chipID),x_pos,y_pos,x_pos_err,y_pos_err);
                                        //-------------------//

                                        //-----------//
                                        // In X dir. //
                                        //-----------//
                                        gr_x->SetPoint(GetPlaneID(chipID),plane_position[GetPlaneID(chipID)],(x_pos-offset_x_ChipID[chipID])*PIXEL_SIZE-offset_xres_Plane[GetPlaneID(chipID)]);
                                        gr_x->SetPointError(GetPlaneID(chipID),plane_position_err[GetPlaneID(chipID)],
                                                TMath::Sqrt(TMath::Power(x_pos_err*PIXEL_SIZE,2) + TMath::Power(offset_errx_ChipID[chipID]*PIXEL_SIZE,2) +
                                                            TMath::Power(offset_errxres_Plane[GetPlaneID(chipID)],2)));

                                        //-----------//
                                        // In Y dir. //
                                        //-----------//
                                        gr_y->SetPoint(GetPlaneID(chipID),plane_position[GetPlaneID(chipID)],(y_pos-offset_y_ChipID[chipID])*PIXEL_SIZE-offset_yres_Plane[GetPlaneID(chipID)]);
                                        gr_y->SetPointError(GetPlaneID(chipID),plane_position_err[GetPlaneID(chipID)],
                                                TMath::Sqrt(TMath::Power(y_pos_err*PIXEL_SIZE,2) + TMath::Power(offset_erry_ChipID[chipID]*PIXEL_SIZE,2) +
                                                            TMath::Power(offset_erryres_Plane[GetPlaneID(chipID)],2)));

                                        gr_point++;
                                    }
                                }
                            }
                            clockID += clock_jitter + 1;
                        }
                    }

                    Bool_t track_is_ready = true;
                    for(Int_t chipID = 0; chipID < N_MAX_CHIP; chipID++)
                    {
                        if(chip_ready[chipID] == 0)
                            track_is_ready = false;
                    }

                    if(track_is_ready && gr_point == N_MAX_CHIP)
                    {
                        Double_t x_res[N_MAX_CHIP] = {};
                        Double_t x_mes[N_MAX_CHIP] = {};
                        Double_t y_mes[N_MAX_CHIP] = {};
                        Double_t y_res[N_MAX_CHIP] = {};
                        Int_t plane_res_ready[N_MAX_CHIP] = {};
                        //-----------//
                        // In X dir. //
                        //-----------//
                        Double_t fit_offset_x_temp  = fit_offset_x;
                        Double_t fit_slope_x_temp   = fit_slope_x;

                        fit_function->SetParameter(0,fit_offset_x_temp);
                        fit_function->SetParameter(1,fit_slope_x_temp);

                        gr_x->Fit(fit_function,"Q");
                        fit_offset_x_temp = fit_function->GetParameter(0);
                        fit_slope_x_temp = fit_function->GetParameter(1);
                        chi2_x = fit_function->GetChisquare()/fit_function->GetNDF();

                        if(chi2_x < chi2_x_max)
                        {
                            h_fit_offset_x->Fill(fit_offset_x_temp);
                            h_fit_slope_x->Fill(fit_slope_x_temp);

                            h25->Fill(iteratorID,chi2_x*fit_function->GetNDF());
                            h26->Fill(iteratorID,chi2_x);

                            for(Int_t gr_point_id = 0; gr_point_id < gr_x->GetN(); gr_point_id++)
                            {
                                Double_t xxx, zzz;
                                gr_x->GetPoint(gr_point_id,zzz,xxx);

                                if(round(zzz) == round(plane_position[0]))
                                {
                                    h_res_plane0_x->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h29->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                        x_res[0] = xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        x_mes[0] = (fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        plane_res_ready[0]++;
                                    }
                                    h49->Fill(iteratorID,xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                }
                                if(round(zzz) == round(plane_position[1]))
                                {
                                    h_res_plane1_x->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h30->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                        x_res[1] = xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        x_mes[1] = (fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        plane_res_ready[1]++;
                                    }
                                    h50->Fill(iteratorID,xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                }
                                if(round(zzz) == round(plane_position[2]))
                                {
                                    h_res_plane2_x->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h31->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                        x_res[2] = xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        x_mes[2] = (fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        plane_res_ready[2]++;
                                    }
                                    h51->Fill(iteratorID,xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                }
                                if(round(zzz) == round(plane_position[3]))
                                {
                                    h_res_plane3_x->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h32->Fill(xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                        x_res[3] = xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        x_mes[3] = (fit_offset_x_temp + fit_slope_x_temp*zzz);
                                        plane_res_ready[3]++;
                                    }
                                    h52->Fill(iteratorID,xxx-(fit_offset_x_temp + fit_slope_x_temp*zzz));
                                }
                            }
                        }
                        //-----------//
                        // In Y dir. //
                        //-----------//
                        Double_t fit_offset_y_temp = fit_offset_y;
                        Double_t fit_slope_y_temp = fit_slope_y;

                        fit_function->SetParameter(0,fit_offset_y_temp);
                        fit_function->SetParameter(1,fit_slope_y_temp);

                        gr_y->Fit(fit_function,"Q");
                        fit_offset_y_temp = fit_function->GetParameter(0);
                        fit_slope_y_temp = fit_function->GetParameter(1);
                        chi2_y = fit_function->GetChisquare()/fit_function->GetNDF();

                        if(chi2_y < chi2_y_max)
                        {
                            h_fit_offset_y->Fill(fit_offset_y_temp);
                            h_fit_slope_y->Fill(fit_slope_y_temp);

                            h27->Fill(iteratorID,chi2_y*fit_function->GetNDF());
                            h28->Fill(iteratorID,chi2_y);

                            for(Int_t gr_point_id = 0; gr_point_id < gr_y->GetN(); gr_point_id++)
                            {
                                Double_t yyy, zzz;
                                gr_y->GetPoint(gr_point_id,zzz,yyy);

                                if(round(zzz) == round(plane_position[0]))
                                {
                                    h_res_plane0_y->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h33->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                        y_res[0] = yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        y_mes[0] = (fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        plane_res_ready[0]++;
                                    }
                                    h120->Fill(iteratorID,yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                }
                                if(round(zzz) == round(plane_position[1]))
                                {
                                    h_res_plane1_y->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h34->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                        y_res[1] = yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        y_mes[1] = (fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        plane_res_ready[1]++;
                                    }
                                    h121->Fill(iteratorID,yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                }
                                if(round(zzz) == round(plane_position[2]))
                                {
                                    h_res_plane2_y->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h35->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                        y_res[2] = yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        y_mes[2] = (fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        plane_res_ready[2]++;
                                    }
                                    h122->Fill(iteratorID,yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                }
                                if(round(zzz) == round(plane_position[3]))
                                {
                                    h_res_plane3_y->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                    if(iteratorID == nIterations)
                                    {
                                        h36->Fill(yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                        y_res[3] = yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        y_mes[3] = (fit_offset_y_temp + fit_slope_y_temp*zzz);
                                        plane_res_ready[3]++;
                                    }
                                    h123->Fill(iteratorID,yyy-(fit_offset_y_temp + fit_slope_y_temp*zzz));
                                }
                            }
                        }

                        //-----------//
                        // Rotation. //
                        //-----------//

                        if(chi2_x < chi2_x_max && chi2_y < chi2_y_max)
                        {
                            if(iteratorID == nIterations)
                            {
                                if(plane_res_ready[0] == 2)
                                {
                                    h124->Fill(y_mes[0],x_res[0]);
                                    h128->Fill(x_mes[0],y_res[0]);
                                }
                                if(plane_res_ready[1] == 2)
                                {
                                    h125->Fill(y_mes[1],x_res[1]);
                                    h129->Fill(x_mes[1],y_res[1]);
                                }
                                if(plane_res_ready[2] == 2)
                                {
                                    h126->Fill(y_mes[2],x_res[2]);
                                    h130->Fill(x_mes[2],y_res[2]);
                                }
                                if(plane_res_ready[3] == 2)
                                {
                                    h127->Fill(y_mes[3],x_res[3]);
                                    h131->Fill(x_mes[3],y_res[3]);
                                }
                            }
                        }
                    }
                }
            }

            //-----------//
            // In X dir. //
            //-----------//
            offset_xres_Plane[0]    += h_res_plane0_x->GetMean();
            offset_errxres_Plane[0] = TMath::Sqrt(TMath::Power(offset_errxres_Plane[0],2) + TMath::Power(h_res_plane0_x->GetMeanError(),2));
            offset_xres_Plane[1]    += h_res_plane1_x->GetMean();
            offset_errxres_Plane[1] = TMath::Sqrt(TMath::Power(offset_errxres_Plane[1],2) + TMath::Power(h_res_plane1_x->GetMeanError(),2));
            offset_xres_Plane[2]    += h_res_plane2_x->GetMean();
            offset_errxres_Plane[2] = TMath::Sqrt(TMath::Power(offset_errxres_Plane[2],2) + TMath::Power(h_res_plane2_x->GetMeanError(),2));
            offset_xres_Plane[3]    += h_res_plane3_x->GetMean();
            offset_errxres_Plane[3] = TMath::Sqrt(TMath::Power(offset_errxres_Plane[3],2) + TMath::Power(h_res_plane3_x->GetMeanError(),2));
            fit_offset_x            = h_fit_offset_x->GetMean();
            fit_slope_x             = h_fit_slope_x->GetMean();

            //-----------//
            // In Y dir. //
            //-----------//
            offset_yres_Plane[0]    += h_res_plane0_y->GetMean();
            offset_erryres_Plane[0] = TMath::Sqrt(TMath::Power(offset_erryres_Plane[0],2) + TMath::Power(h_res_plane0_y->GetMeanError(),2));
            offset_yres_Plane[1]    += h_res_plane1_y->GetMean();
            offset_erryres_Plane[1] = TMath::Sqrt(TMath::Power(offset_erryres_Plane[1],2) + TMath::Power(h_res_plane1_y->GetMeanError(),2));
            offset_yres_Plane[2]    += h_res_plane2_y->GetMean();
            offset_erryres_Plane[2] = TMath::Sqrt(TMath::Power(offset_erryres_Plane[2],2) + TMath::Power(h_res_plane2_y->GetMeanError(),2));
            offset_yres_Plane[3]    += h_res_plane3_y->GetMean();
            offset_erryres_Plane[3] = TMath::Sqrt(TMath::Power(offset_erryres_Plane[3],2) + TMath::Power(h_res_plane3_y->GetMeanError(),2));
            fit_offset_y            = h_fit_offset_y->GetMean();
            fit_slope_y             = h_fit_slope_y->GetMean();

            // Clean histos
            h_res_plane0_x->Delete();
            h_res_plane1_x->Delete();
            h_res_plane2_x->Delete();
            h_res_plane3_x->Delete();
            h_res_plane0_y->Delete();
            h_res_plane1_y->Delete();
            h_res_plane2_y->Delete();
            h_res_plane3_y->Delete();
            h_fit_offset_x->Delete();
            h_fit_slope_x->Delete();
            h_fit_offset_y->Delete();
            h_fit_slope_y->Delete();
        }

        cout<<endl<<"--> Residual corrections: "<<endl;
        for(Int_t i = 0; i < N_MAX_CHIP; i++)
        {
            cout<<">> Plane["<<i<<"][X] = "<<offset_xres_Plane[i]<<" +/- "<<offset_errxres_Plane[i]<<" [mm]"<<endl;
            cout<<">> Plane["<<i<<"][Y] = "<<offset_yres_Plane[i]<<" +/- "<<offset_erryres_Plane[i]<<" [mm]"<<endl;
        }
    }

    gr_x->Delete();
    gr_y->Delete();

    if(Run_Mode == 2)
    {
        cout<<endl<<"--> Track reconstruction 2 <--"<<endl;

        //----------------//
        // Loop by frames //
        //----------------//
        for(Int_t i = 0; i < nEntriesClock; i++)
        {
            if(i%10 == 0)
            {
                printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntriesClock);
                fflush(stdout);
            }

            fChain_clock->GetEntry(i);

            //----------------//
            // Loop by clocks //
            //----------------//
            for(Int_t clockID = 0; clockID < N_MAX_CLOCKS+1; clockID++)
            {
                UInt_t gr_point = 0;
                Int_t chip_ready[N_MAX_CHIP] = {};
                Double_t x_coord[N_MAX_CHIP], y_coord[N_MAX_CHIP], errx_coord[N_MAX_CHIP], erry_coord[N_MAX_CHIP];

                if(clockID*1.0e3/_Clock < timecut) continue;

                //----------------//
                // Loop by chips  //
                //----------------//
                for(Int_t chipID_temp = 0; chipID_temp < N_MAX_CHIP; chipID_temp++)
                {
                    if(clock_info[chipID_temp][clockID] > 0)
                    {
                        for(Int_t chipID = 0; chipID < N_MAX_CHIP; chipID++)
                        {
                            //-----------------------------//
                            // Loop by clocks inside range //
                            //-----------------------------//
                            for(Int_t kkl = clockID-clock_jitter; kkl <= clockID+clock_jitter; kkl++)
                            {
                                if(kkl < 0 || kkl > N_MAX_CLOCKS) continue;

                                //-------------------------------//
                                // Loop by clusters inside clock //
                                //-------------------------------//
                                for(UInt_t kkll = 0; kkll < clock_info[chipID][kkl]; kkll++)
                                {
                                    chip_ready[chipID]++;
                                    UInt_t clusterID = -1;
                                    switch (chipID)
                                    {
                                    case 0:
                                        clusterID = clock_info_0[kkl][kkll];
                                        break;
                                    case 1:
                                        clusterID = clock_info_1[kkl][kkll];
                                        break;
                                    case 2:
                                        clusterID = clock_info_2[kkl][kkll];
                                        break;
                                    case 3:
                                        clusterID = clock_info_3[kkl][kkll];
                                        break;
                                    }
                                    fChain[chipID]->GetEntry(clusterID);

                                    //--- Orientation ---//
                                    Double_t x_pos      = pos_x_clinfo;
                                    Double_t x_pos_err  = pos_x_err_clinfo;
                                    Double_t y_pos      = pos_y_clinfo;
                                    Double_t y_pos_err  = pos_y_err_clinfo;

                                    if(chipID == 1 || chipID == 3)
                                    {
                                        x_pos       = (-1.0)*pos_y_clinfo;
                                        x_pos_err   = pos_y_err_clinfo;
                                        y_pos       = pos_x_clinfo;
                                        y_pos_err   = pos_x_err_clinfo;
                                    }

                                    GetRotaion(GetPlaneID(chipID),x_pos,y_pos,x_pos_err,y_pos_err);
                                    //-------------------//

                                    //-----------//
                                    // In X dir. //
                                    //-----------//
                                    x_coord[GetPlaneID(chipID)] = (x_pos-offset_x_ChipID[chipID])*PIXEL_SIZE-offset_xres_Plane[GetPlaneID(chipID)];
                                    errx_coord[GetPlaneID(chipID)] = TMath::Sqrt(TMath::Power(x_pos_err*PIXEL_SIZE,2) +
                                                                                 TMath::Power(offset_errx_ChipID[chipID]*PIXEL_SIZE,2) +
                                                                                 TMath::Power(offset_errxres_Plane[GetPlaneID(chipID)],2));

                                    //-----------//
                                    // In Y dir. //
                                    //-----------//
                                    y_coord[GetPlaneID(chipID)] = (y_pos-offset_y_ChipID[chipID])*PIXEL_SIZE-offset_yres_Plane[GetPlaneID(chipID)];
                                    erry_coord[GetPlaneID(chipID)] = TMath::Sqrt(TMath::Power(y_pos_err*PIXEL_SIZE,2) +
                                                                                 TMath::Power(offset_erry_ChipID[chipID]*PIXEL_SIZE,2) +
                                                                                 TMath::Power(offset_erryres_Plane[GetPlaneID(chipID)],2));

                                    gr_point++;
                                }
                            }
                        }
                        clockID += clock_jitter + 1;
                    }
                }

                Bool_t track_is_ready = true;
                for(Int_t chipID = 0; chipID < N_MAX_CHIP; chipID++)
                {
                    if(chip_ready[chipID] == 0)
                        track_is_ready = false;
                }

                if(track_is_ready && gr_point == N_MAX_CHIP)
                {
                    //-----------//
                    // In X dir. //
                    //-----------//
                    slope_arm1_x = TMath::ATan((x_coord[1] - x_coord[0])/(plane_position[1] - plane_position[0]));
                    slope_arm2_x = TMath::ATan((x_coord[3] - x_coord[2])/(plane_position[3] - plane_position[2]));
                    slope_arm1_x_err = TMath::Sqrt(TMath::Power(errx_coord[1]/(plane_position[1] - plane_position[0]),2) +
                            TMath::Power(errx_coord[0]/(plane_position[1] - plane_position[0]),2) +
                            TMath::Power(plane_position_err[1]/TMath::Power(plane_position[1] - plane_position[0],2),2) +
                            TMath::Power(plane_position_err[0]/TMath::Power(plane_position[1] - plane_position[0],2),2))/(1.0 + TMath::Power((x_coord[1] - x_coord[0])/(plane_position[1] - plane_position[0]),2));
                    slope_arm2_x_err = TMath::Sqrt(TMath::Power(errx_coord[3]/(plane_position[3] - plane_position[2]),2) +
                            TMath::Power(errx_coord[2]/(plane_position[3] - plane_position[2]),2) +
                            TMath::Power(plane_position_err[3]/TMath::Power(plane_position[3] - plane_position[2],2),2) +
                            TMath::Power(plane_position_err[2]/TMath::Power(plane_position[3] - plane_position[2],2),2))/(1.0 + TMath::Power((x_coord[3] - x_coord[2])/(plane_position[3] - plane_position[2]),2));
                    /*offset_arm1_x = x_coord[1] - TMath::Tan(slope_arm1_x)*plane_position[1];
                    offset_arm2_x = x_coord[3] - TMath::Tan(slope_arm2_x)*plane_position[3];
                    offset_arm1_x_err = TMath::Sqrt(TMath::Power(errx_coord[1],2) +
                            TMath::Power(slope_arm1_x_err*plane_position[1]/(TMath::Cos(slope_arm1_x)*TMath::Cos(slope_arm1_x)),2) +
                            TMath::Power(slope_arm1_x*plane_position_err[1],2));
                    offset_arm2_x_err = TMath::Sqrt(TMath::Power(errx_coord[3],2) +
                            TMath::Power(slope_arm2_x_err*plane_position[3]/(TMath::Cos(slope_arm2_x)*TMath::Cos(slope_arm2_x)),2) +
                            TMath::Power(slope_arm2_x*plane_position_err[3],2));
                    x_arm1 = TMath::Tan(slope_arm1_x)*projection_plane_position + offset_arm1_x;
                    x_arm2 = TMath::Tan(slope_arm2_x)*projection_plane_position + offset_arm2_x;*/
                    x_arm1 = TMath::Tan(slope_arm1_x)*projection_plane_position;
                    x_arm2 = TMath::Tan(slope_arm2_x)*projection_plane_position;
                    x_arm1_err = TMath::Sqrt(TMath::Power(slope_arm1_x_err*projection_plane_position/(TMath::Cos(slope_arm1_x)*TMath::Cos(slope_arm1_x)),2) +
                                             TMath::Power(slope_arm1_x*projection_plane_position_err,2) /*+
                                             TMath::Power(offset_arm1_x_err,2)*/);
                    x_arm2_err = TMath::Sqrt(TMath::Power(slope_arm2_x_err*projection_plane_position/(TMath::Cos(slope_arm2_x)*TMath::Cos(slope_arm2_x)),2) +
                                             TMath::Power(slope_arm2_x*projection_plane_position_err,2) /*+
                                             TMath::Power(offset_arm2_x_err,2)*/);
                    chi2_x = TMath::Power(x_arm1-x_arm2,2)/(x_arm1_err*x_arm1_err + x_arm2_err*x_arm2_err);

                    h53->Fill(slope_arm1_x);
                    h55->Fill(slope_arm2_x);
                    h57->Fill(x_arm1-x_arm2);
                    h59->Fill(x_arm1);
                    h60->Fill(x_arm2);
                    h63->Fill(chi2_x);
                    h65->Fill(x_arm1-x_arm2,chi2_x);

                    //-----------//
                    // In Y dir. //
                    //-----------//
                    slope_arm1_y = TMath::ATan((y_coord[1] - y_coord[0])/(plane_position[1] - plane_position[0]));
                    slope_arm2_y = TMath::ATan((y_coord[3] - y_coord[2])/(plane_position[3] - plane_position[2]));
                    slope_arm1_y_err = TMath::Sqrt(TMath::Power(erry_coord[1]/(plane_position[1] - plane_position[0]),2) +
                            TMath::Power(erry_coord[0]/(plane_position[1] - plane_position[0]),2) +
                            TMath::Power(plane_position_err[1]/TMath::Power(plane_position[1] - plane_position[0],2),2) +
                            TMath::Power(plane_position_err[0]/TMath::Power(plane_position[1] - plane_position[0],2),2))/(1.0 + TMath::Power((y_coord[1] - y_coord[0])/(plane_position[1] - plane_position[0]),2));
                    slope_arm2_y_err = TMath::Sqrt(TMath::Power(erry_coord[3]/(plane_position[3] - plane_position[2]),2) +
                            TMath::Power(erry_coord[2]/(plane_position[3] - plane_position[2]),2) +
                            TMath::Power(plane_position_err[3]/TMath::Power(plane_position[3] - plane_position[2],2),2) +
                            TMath::Power(plane_position_err[2]/TMath::Power(plane_position[3] - plane_position[2],2),2))/(1.0 + TMath::Power((y_coord[3] - y_coord[2])/(plane_position[3] - plane_position[2]),2));
                    /*offset_arm1_y = y_coord[1] - TMath::Tan(slope_arm1_y)*plane_position[1];
                    offset_arm2_y = y_coord[3] - TMath::Tan(slope_arm2_y)*plane_position[3];
                    offset_arm1_y_err = TMath::Sqrt(TMath::Power(erry_coord[1],2) +
                            TMath::Power(slope_arm1_y_err*plane_position[1]/(TMath::Cos(slope_arm1_y)*TMath::Cos(slope_arm1_y)),2) +
                            TMath::Power(slope_arm1_y*plane_position_err[1],2));
                    offset_arm2_y_err = TMath::Sqrt(TMath::Power(erry_coord[3],2) +
                            TMath::Power(slope_arm2_y_err*plane_position[3]/(TMath::Cos(slope_arm2_y)*TMath::Cos(slope_arm2_y)),2) +
                            TMath::Power(slope_arm2_y*plane_position_err[3],2));
                    y_arm1 = TMath::Tan(slope_arm1_y)*projection_plane_position + offset_arm1_y;
                    y_arm2 = TMath::Tan(slope_arm2_y)*projection_plane_position + offset_arm2_y;*/
                    y_arm1 = TMath::Tan(slope_arm1_y)*projection_plane_position;
                    y_arm2 = TMath::Tan(slope_arm2_y)*projection_plane_position;
                    y_arm1_err = TMath::Sqrt(TMath::Power(slope_arm1_y_err*projection_plane_position/(TMath::Cos(slope_arm1_y)*TMath::Cos(slope_arm1_y)),2) +
                                             TMath::Power(slope_arm1_y*projection_plane_position_err,2) /*+
                                             TMath::Power(offset_arm1_y_err,2)*/);
                    y_arm2_err = TMath::Sqrt(TMath::Power(slope_arm2_y_err*projection_plane_position/(TMath::Cos(slope_arm2_y)*TMath::Cos(slope_arm2_y)),2) +
                                             TMath::Power(slope_arm2_y*projection_plane_position_err,2) /*+
                                             TMath::Power(offset_arm2_y_err,2)*/);
                    chi2_y = TMath::Power(y_arm1-y_arm2,2)/(y_arm1_err*y_arm1_err + y_arm2_err*y_arm2_err);

                    h54->Fill(slope_arm1_y);
                    h56->Fill(slope_arm2_y);
                    h58->Fill(y_arm1-y_arm2);
                    h61->Fill(y_arm1);
                    h62->Fill(y_arm2);
                    h64->Fill(chi2_y);
                    h66->Fill(y_arm1-y_arm2,chi2_y);


                    h67->Fill(x_coord[0],x_coord[3] - x_coord[0]);
                    h68->Fill(y_coord[0],x_coord[3] - x_coord[0]);
                    h69->Fill(x_coord[0],y_coord[3] - y_coord[0]);
                    h70->Fill(y_coord[0],y_coord[3] - y_coord[0]);
                    h71->Fill(x_coord[1],x_coord[3] - x_coord[1]);
                    h72->Fill(y_coord[1],x_coord[3] - x_coord[1]);
                    h73->Fill(x_coord[1],y_coord[3] - y_coord[1]);
                    h74->Fill(y_coord[1],y_coord[3] - y_coord[1]);
                    h75->Fill(x_coord[2],x_coord[3] - x_coord[2]);
                    h76->Fill(y_coord[2],x_coord[3] - x_coord[2]);
                    h77->Fill(x_coord[2],y_coord[3] - y_coord[2]);
                    h78->Fill(y_coord[2],y_coord[3] - y_coord[2]);
                    h79->Fill(y_coord[0],y_coord[3]);
                    h80->Fill(y_coord[1],y_coord[3]);
                    h81->Fill(y_coord[2],y_coord[3]);

                    h86->Fill(y_coord[0],x_coord[3] - x_coord[0]);
                    h87->Fill(x_coord[0],y_coord[3] - y_coord[0]);
                    h88->Fill(y_coord[1],x_coord[3] - x_coord[1]);
                    h89->Fill(x_coord[1],y_coord[3] - y_coord[1]);
                    h90->Fill(y_coord[2],x_coord[3] - x_coord[2]);
                    h91->Fill(x_coord[2],y_coord[3] - y_coord[2]);


                    h83->Fill((slope_arm1_x - slope_arm2_x)*1e6);
                    h84->Fill((slope_arm1_y - slope_arm2_y)*1e6);


                    h85->Fill(x_arm1,slope_arm2_x - slope_arm1_x);
                    h116->Fill(x_arm1,slope_arm2_y - slope_arm1_y);
                    h117->Fill(y_arm1,slope_arm2_y - slope_arm1_y);

                    h92->Fill(TMath::ATan(x_coord[0]/y_coord[0]));
                    h93->Fill(TMath::ATan(x_coord[1]/y_coord[1]));
                    h94->Fill(TMath::ATan(x_coord[2]/y_coord[2]));
                    h95->Fill(TMath::ATan(x_coord[3]/y_coord[3]));
                    h96->Fill(TMath::ATan(x_coord[0]/y_coord[0]) - TMath::ATan(x_coord[3]/y_coord[3]));
                    h97->Fill(TMath::ATan(x_coord[1]/y_coord[1]) - TMath::ATan(x_coord[3]/y_coord[3]));
                    h98->Fill(TMath::ATan(x_coord[2]/y_coord[2]) - TMath::ATan(x_coord[3]/y_coord[3]));

                    h99->Fill(x_coord[0]);
                    h100->Fill(x_coord[1]);
                    h101->Fill(x_coord[2]);
                    h102->Fill(x_coord[3]);
                    h103->Fill(y_coord[0]);
                    h104->Fill(y_coord[1]);
                    h105->Fill(y_coord[2]);
                    h106->Fill(y_coord[3]);
                    h107->Fill(x_coord[3],x_coord[0]);
                    h108->Fill(x_coord[3],x_coord[1]);
                    h109->Fill(x_coord[3],x_coord[2]);
                    h110->Fill(x_coord[0],y_coord[0]);
                    h111->Fill(x_coord[1],y_coord[1]);
                    h112->Fill(x_coord[2],y_coord[2]);
                    h113->Fill(x_coord[3],y_coord[3]);

                    if(slope_arm2_x - slope_arm1_x > 0.0004)
                    {
                        h114->Fill(x_coord[0],y_coord[0]);
                        h118->Fill(x_coord[3],y_coord[3]);
                    }
                    else
                    {
                        h115->Fill(x_coord[0],y_coord[0]);
                        h119->Fill(x_coord[3],y_coord[3]);
                    }
                }
            }
        }
    }
    //------------------------------------------------------------------------------------------------//
    cout<<endl;
    //------------------------------------------------------------------------------------------------//
    h1->GetXaxis()->SetTitle("#Delta Time [ns]");
    h2->GetXaxis()->SetTitle("#Delta Time [ns]");
    h3->GetXaxis()->SetTitle("#Delta Time [ns]");
    h4->GetXaxis()->SetTitle("[pixels]");
    h5->GetXaxis()->SetTitle("[pixels]");
    h6->GetXaxis()->SetTitle("[pixels]");
    h7->GetXaxis()->SetTitle("[pixels]");
    h8->GetXaxis()->SetTitle("[pixels]");
    h9->GetXaxis()->SetTitle("[pixels]");
    h10->GetXaxis()->SetTitle("[pixels]");
    h11->GetXaxis()->SetTitle("[pixels]");
    h12->GetXaxis()->SetTitle("[pixels]");
    h13->GetXaxis()->SetTitle("[pixels]");
    h14->GetXaxis()->SetTitle("[pixels]");
    h15->GetXaxis()->SetTitle("[pixels]");
    h16->GetXaxis()->SetTitle("[pixels]");
    h17->GetXaxis()->SetTitle("[pixels]");
    h18->GetXaxis()->SetTitle("[pixels]");
    h19->GetXaxis()->SetTitle("[pixels]");
    h21->GetYaxis()->SetTitle("X position [mm]");
    h21->GetXaxis()->SetTitle("Z position [mm]");
    h22->GetYaxis()->SetTitle("Y position [mm]");
    h22->GetXaxis()->SetTitle("Z position [mm]");

    h29->GetXaxis()->SetTitle("Residual X [mm]");
    h30->GetXaxis()->SetTitle("Residual X [mm]");
    h31->GetXaxis()->SetTitle("Residual X [mm]");
    h32->GetXaxis()->SetTitle("Residual X [mm]");
    h33->GetXaxis()->SetTitle("Residual Y [mm]");
    h34->GetXaxis()->SetTitle("Residual Y [mm]");
    h35->GetXaxis()->SetTitle("Residual Y [mm]");
    h36->GetXaxis()->SetTitle("Residual Y [mm]");

    h41->GetXaxis()->SetTitle("[clocks]");
    h42->GetXaxis()->SetTitle("[clocks]");
    h43->GetXaxis()->SetTitle("[clocks]");
    h44->GetXaxis()->SetTitle("[clocks]");
    h45->GetXaxis()->SetTitle("[pixels]");
    h46->GetXaxis()->SetTitle("[pixels]");
    h47->GetXaxis()->SetTitle("[pixels]");
    h48->GetXaxis()->SetTitle("[pixels]");

    h49->GetXaxis()->SetTitle("IteratorID");
    h49->GetYaxis()->SetTitle("Residual X [mm]");
    h50->GetXaxis()->SetTitle("IteratorID");
    h50->GetYaxis()->SetTitle("Residual X [mm]");
    h51->GetXaxis()->SetTitle("IteratorID");
    h51->GetYaxis()->SetTitle("Residual X [mm]");
    h52->GetXaxis()->SetTitle("IteratorID");
    h52->GetYaxis()->SetTitle("Residual X [mm]");

    h53->GetXaxis()->SetTitle("[rad]");
    h54->GetXaxis()->SetTitle("[rad]");
    h55->GetXaxis()->SetTitle("[rad]");
    h56->GetXaxis()->SetTitle("[rad]");
    h57->GetXaxis()->SetTitle("[mm]");
    h58->GetXaxis()->SetTitle("[mm]");
    h59->GetXaxis()->SetTitle("[mm]");
    h60->GetXaxis()->SetTitle("[mm]");
    h61->GetXaxis()->SetTitle("[mm]");
    h62->GetXaxis()->SetTitle("[mm]");
    h65->GetXaxis()->SetTitle("dX [mm]");
    h65->GetYaxis()->SetTitle("#Chi^{2}_{X}");
    h66->GetXaxis()->SetTitle("dY [mm]");
    h66->GetYaxis()->SetTitle("#Chi^{2}_{Y}");
    h67->GetXaxis()->SetTitle("X [mm]");
    h67->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h68->GetXaxis()->SetTitle("Y [mm]");
    h68->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h69->GetXaxis()->SetTitle("X [mm]");
    h69->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h70->GetXaxis()->SetTitle("Y [mm]");
    h70->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h71->GetXaxis()->SetTitle("X [mm]");
    h71->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h72->GetXaxis()->SetTitle("Y [mm]");
    h72->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h73->GetXaxis()->SetTitle("X [mm]");
    h73->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h74->GetXaxis()->SetTitle("Y [mm]");
    h74->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h75->GetXaxis()->SetTitle("X [mm]");
    h75->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h76->GetXaxis()->SetTitle("Y [mm]");
    h76->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h77->GetXaxis()->SetTitle("X [mm]");
    h77->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h78->GetXaxis()->SetTitle("Y [mm]");
    h78->GetYaxis()->SetTitle("#Delta_{Y} [mm]");

    h83->GetXaxis()->SetTitle("[#murad]");
    h84->GetXaxis()->SetTitle("[#murad]");


    h85->GetXaxis()->SetTitle("X [mm]");
    h85->GetYaxis()->SetTitle("#Delta#Theta_{X} [rad]");

    h86->GetXaxis()->SetTitle("Y [mm]");
    h86->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h87->GetXaxis()->SetTitle("X [mm]");
    h87->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h88->GetXaxis()->SetTitle("Y [mm]");
    h88->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h89->GetXaxis()->SetTitle("X [mm]");
    h89->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h90->GetXaxis()->SetTitle("Y [mm]");
    h90->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    h91->GetXaxis()->SetTitle("X [mm]");
    h91->GetYaxis()->SetTitle("#Delta_{Y} [mm]");
    h92->GetXaxis()->SetTitle("[rad]");
    h93->GetXaxis()->SetTitle("[rad]");
    h94->GetXaxis()->SetTitle("[rad]");
    h95->GetXaxis()->SetTitle("[rad]");
    h96->GetXaxis()->SetTitle("[rad]");
    h97->GetXaxis()->SetTitle("[rad]");
    h98->GetXaxis()->SetTitle("[rad]");
    h99->GetXaxis()->SetTitle("[mm]");
    h100->GetXaxis()->SetTitle("[mm]");
    h101->GetXaxis()->SetTitle("[mm]");
    h102->GetXaxis()->SetTitle("[mm]");
    h103->GetXaxis()->SetTitle("[mm]");
    h104->GetXaxis()->SetTitle("[mm]");
    h105->GetXaxis()->SetTitle("[mm]");
    h106->GetXaxis()->SetTitle("[mm]");
    h107->GetXaxis()->SetTitle("[mm]");
    h107->GetYaxis()->SetTitle("[mm]");
    h108->GetXaxis()->SetTitle("[mm]");
    h108->GetYaxis()->SetTitle("[mm]");
    h109->GetXaxis()->SetTitle("[mm]");
    h109->GetYaxis()->SetTitle("[mm]");
    h110->GetXaxis()->SetTitle("[mm]");
    h110->GetYaxis()->SetTitle("[mm]");
    h111->GetXaxis()->SetTitle("[mm]");
    h111->GetYaxis()->SetTitle("[mm]");
    h112->GetXaxis()->SetTitle("[mm]");
    h112->GetYaxis()->SetTitle("[mm]");
    h113->GetXaxis()->SetTitle("[mm]");
    h113->GetYaxis()->SetTitle("[mm]");
    h114->GetXaxis()->SetTitle("[mm]");
    h114->GetYaxis()->SetTitle("[mm]");
    h115->GetXaxis()->SetTitle("[mm]");
    h115->GetYaxis()->SetTitle("[mm]");
    h116->GetXaxis()->SetTitle("X [mm]");
    h116->GetYaxis()->SetTitle("#Delta#Theta_{Y} [rad]");
    h117->GetXaxis()->SetTitle("Y [mm]");
    h117->GetYaxis()->SetTitle("#Delta#Theta_{Y} [rad]");
    h118->GetXaxis()->SetTitle("[mm]");
    h118->GetYaxis()->SetTitle("[mm]");
    h119->GetXaxis()->SetTitle("[mm]");
    h119->GetYaxis()->SetTitle("[mm]");
    h120->GetXaxis()->SetTitle("IteratorID");
    h120->GetYaxis()->SetTitle("Residual Y [mm]");
    h121->GetXaxis()->SetTitle("IteratorID");
    h121->GetYaxis()->SetTitle("Residual Y [mm]");
    h122->GetXaxis()->SetTitle("IteratorID");
    h122->GetYaxis()->SetTitle("Residual Y [mm]");
    h123->GetXaxis()->SetTitle("IteratorID");
    h123->GetYaxis()->SetTitle("Residual Y [mm]");
    h124->GetXaxis()->SetTitle("Y Impact Position [mm]");
    h124->GetYaxis()->SetTitle("Residual in X [mm]");
    h125->GetXaxis()->SetTitle("Y Impact Position [mm]");
    h125->GetYaxis()->SetTitle("Residual in X [mm]");
    h126->GetXaxis()->SetTitle("Y Impact Position [mm]");
    h126->GetYaxis()->SetTitle("Residual in X [mm]");
    h127->GetXaxis()->SetTitle("Y Impact Position [mm]");
    h127->GetYaxis()->SetTitle("Residual in X [mm]");
    h128->GetXaxis()->SetTitle("X Impact Position [mm]");
    h128->GetYaxis()->SetTitle("Residual in Y [mm]");
    h129->GetXaxis()->SetTitle("X Impact Position [mm]");
    h129->GetYaxis()->SetTitle("Residual in Y [mm]");
    h130->GetXaxis()->SetTitle("X Impact Position [mm]");
    h130->GetYaxis()->SetTitle("Residual in Y [mm]");
    h131->GetXaxis()->SetTitle("X Impact Position [mm]");
    h131->GetYaxis()->SetTitle("Residual in Y [mm]");

    cout<<endl<<"--> Writing a histograms to the output file ..."<<endl;
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h5->Write();
    h6->Write();
    h7->Write();
    h8->Write();
    h9->Write();
    h10->Write();
    h11->Write();
    h12->Write();
    h13->Write();
    h14->Write();
    h15->Write();
    h16->Write();
    h17->Write();
    h18->Write();
    h19->Write();
    h21->Write();
    h22->Write();
    h25->Write();
    h26->Write();
    h27->Write();
    h28->Write();
    h29->Write();
    h30->Write();
    h31->Write();
    h32->Write();
    h33->Write();
    h34->Write();
    h35->Write();
    h36->Write();
    h37->Write();
    h38->Write();
    h39->Write();
    h40->Write();
    h41->Write();
    h42->Write();
    h43->Write();
    h44->Write();
    h45->Write();
    h46->Write();
    h47->Write();
    h48->Write();
    h49->Write();
    h50->Write();
    h51->Write();
    h52->Write();
    h53->Write();
    h54->Write();
    h55->Write();
    h56->Write();
    h57->Write();
    h58->Write();
    h59->Write();
    h60->Write();
    h61->Write();
    h62->Write();
    h63->Write();
    h64->Write();
    h65->Write();
    h66->Write();
    h67->Write();
    h68->Write();
    h69->Write();
    h70->Write();
    h71->Write();
    h72->Write();
    h73->Write();
    h74->Write();
    h75->Write();
    h76->Write();
    h77->Write();
    h78->Write();
    h79->Write();
    h80->Write();
    h81->Write();


    h83->Write();
    h84->Write();


    h85->Write();
    h86->Write();
    h87->Write();
    h88->Write();
    h89->Write();
    h90->Write();
    h91->Write();
    h92->Write();
    h93->Write();
    h94->Write();
    h95->Write();
    h96->Write();
    h97->Write();
    h98->Write();
    h99->Write();
    h100->Write();
    h101->Write();
    h102->Write();
    h103->Write();
    h104->Write();
    h105->Write();
    h106->Write();
    h107->Write();
    h108->Write();
    h109->Write();
    h110->Write();
    h111->Write();
    h112->Write();
    h113->Write();
    h114->Write();
    h115->Write();
    h116->Write();
    h117->Write();
    h118->Write();
    h119->Write();
    h120->Write();
    h121->Write();
    h122->Write();
    h123->Write();
    h124->Write();
    h125->Write();
    h126->Write();
    h127->Write();
    h128->Write();
    h129->Write();
    h130->Write();
    h131->Write();

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

int GetPlaneID(Int_t chipID)
{
    if(chipID < 0 || chipID >= N_MAX_CHIP)
    {
        cout<<endl<<"--> ERROR:: Wrong chipID!!!"<<endl<<endl;
        assert(0);
    }
    switch (chipID)
    {
    case 0:
        return 1;
    case 1:
        return 2;
    case 2:
        return 0;
    case 3:
        return 3;
    }
    return -1;
}

void GetRotaion(Int_t planeID, Double_t &x_pos, Double_t &y_pos, Double_t &x_pos_err, Double_t &y_pos_err)
{
    Double_t p1, p1_err;
    switch (planeID)
    {
    case 0:
        p1      =   -0.018214;
        p1_err  =    0.000138;
        break;
    case 1:
        p1      =    0.0228305;
        p1_err  =    0.0001347;
        break;
    case 2:
        p1      =   -0.0020796;
        p1_err  =    0.0001481;
        break;
    default:
        p1      =   -0.0065317;
        p1_err  =    0.0001104;
        break;
    }

    Double_t x_pos_old      = x_pos;
    Double_t y_pos_old      = y_pos;
    Double_t x_pos_err_old  = x_pos_err;
    Double_t y_pos_err_old  = y_pos_err;

    x_pos = x_pos_old - p1*y_pos_old;
    y_pos = y_pos_old + p1*x_pos_old;

    x_pos_err = TMath::Sqrt(TMath::Power(x_pos_err_old,2) + TMath::Power(p1_err*y_pos_old,2) + TMath::Power(p1*y_pos_err_old,2));
    y_pos_err = TMath::Sqrt(TMath::Power(y_pos_err_old,2) + TMath::Power(p1_err*x_pos_old,2) + TMath::Power(p1*x_pos_err_old,2));

//    Double_t p0x, p1x, p0y, p1y;
//    switch (planeID)
//    {
//    case 0:
//        p0x =  0.1694;
//        p1x =  0.0006226;
//        p0y =  0.09406;
//        p1y = -0.001413;
//        break;
//    case 1:
//        p0x =  0.1719;
//        p1x = -0.03652;
//        p0y =  0.07538;
//        p1y =  0.03714;
//        break;
//    case 2:
//        p0x =  0.06158;
//        p1x = -0.005974;
//        p0y =  0.01573;
//        p1y =  0.006431;
//        break;
//    default:
//        p0x = 0.0;
//        p1x = 0.0;
//        p0y = 0.0;
//        p1y = 0.0;
//        break;
//    }

//    Double_t x_pos_old      = x_pos;
//    Double_t y_pos_old      = y_pos;

//    x_pos = x_pos_old + (p0x + p1x*y_pos_old);
//    y_pos = y_pos_old + (p0y + p1y*x_pos_old);
}
