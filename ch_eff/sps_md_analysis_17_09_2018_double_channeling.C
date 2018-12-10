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
void rewriteHisto(TH2D* h_in, TH2D* h_out);
double getAverage2D(TH2D* histo);
int median_filter(TH2D* h_in, TH2D* h_out, Double_t factor_bad_pixels);
double findMedian(double a[], int n);

int sps_md_analysis_17_09_2018_double_channeling()
{
    cout<<"--> function_1() -- [L1-L19] to plot rp_i 2d beam image normalized by rp_i counts"<<endl;
    cout<<"--> function_2() -- [L1] to plot Timepix data vs CRY3 angle"<<endl;
    cout<<"--> function_3() -- [L3] to plot Timepix data vs CRY2 angle"<<endl;
    cout<<"--> function_4() -- [L2,L7-L10] to plot Timepix projection"<<endl;
    cout<<"--> function_5() -- [L2] double CH stability"<<endl;
    cout<<"--> function_6() -- [L1] VR projection"<<endl;
    cout<<"--> function_7() -- [L5] to plot Timepix data vs CRY3 position"<<endl;
    cout<<"--> function_8() -- [L6] to plot Timepix data vs CRY3 position"<<endl;
    cout<<"--> function_9() -- [L5,L6] spatial distribution of the first CH beam"<<endl;
    cout<<"--> function_10() -- [L1,L8] to calibrate the RP1I w.r.t. the RP0"<<endl;
    cout<<"--> function_11() -- [L1] to plot Timepix data vs CRY3 angle"<<endl;
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

    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L2_RUN_8.root";
    TString outFileName     = "output_L2_chipid_1.root";

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

    TH2D* h_3   = new TH2D("h_3","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
    TH2D* h_8   = new TH2D("h_8","Y vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH2D* h_9   = new TH2D("h_9","X vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);
    TH2D* h_4   = new TH2D("h_4","RP1 Internal (normalized)",N_PIXELS,0,N_PIXELS*0.055,N_PIXELS,0,N_PIXELS*0.055);
    TH1D* h_5   = new TH1D("h_5","Number of the filtered pixels per frame",10000,0,10000);

    h_4->GetXaxis()->SetTitle("Horizontal Axis [mm]");
    h_4->GetXaxis()->SetTimeOffset(1.2);
    h_4->GetXaxis()->CenterTitle();
    h_4->GetYaxis()->SetTitle("Vertical Axis [mm]");
    h_4->GetYaxis()->SetTimeOffset(1.2);
    h_4->GetYaxis()->CenterTitle();

    //------------------------------------------------------------------------------//
    Double_t event_time = -999.999;

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------- MAIN LOOP ------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    for(Long64_t i = 0; i < nEntries; i++)
    {
        TH2D* h_1 = new TH2D("h_1","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
        TH2D* h_2 = new TH2D("h_2","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

        if(i%1 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntries);
            fflush(stdout);
        }


        fChain->GetEntry(i);
        if (_Timems > 0)
        {
            event_time  = _Timems;
        }
        else
        {
            continue;
        }

        for(Int_t xi = 1; xi < N_PIXELS - 1; xi++)
        {
            for(Int_t yi = 1; yi < N_PIXELS - 1; yi++)
            {
                if(_COUNTS[xi][yi] > 0 && _AcquisType == 0)
                {
                    h_1->Fill(xi,yi,_COUNTS[xi][yi]);
                    h_2->Fill(xi,yi,_COUNTS[xi][yi]);
                }
            }
        }

        h_5->Fill(median_filter(h_1,h_2,2.0));

        //----------------------------------------------------//
        // To check the algorythm

//        TCanvas* c_3 = new TCanvas("c_3","c_3",1800,900);
//        c_3->Divide(2,1);

//        c_3->cd(1);
//        gPad->SetGrid();
//        h_1->GetXaxis()->SetRange(1,256);
//        h_1->GetYaxis()->SetRange(1,256);
//        h_1->Draw("colz");

//        c_3->cd(2);
//        gPad->SetGrid();
//        h_2->GetXaxis()->SetRange(1,256);
//        h_2->GetYaxis()->SetRange(1,256);
//        h_2->Draw("colz");
//        break;
        //----------------------------------------------------//

        Double_t integral, integral_err;
        integral= h_2->IntegralAndError(1,h_2->GetNbinsX(),1,h_2->GetNbinsY(),integral_err);

        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                h_3->Fill(xi,yi,h_2->GetBinContent(xi+1,yi+1));
                h_8->Fill(event_time/1000.0,yi,h_2->GetBinContent(xi+1,yi+1));

                if(yi < 100)
                    h_9->Fill(event_time/1000.0,xi,h_2->GetBinContent(xi+1,yi+1));
            }
        }

        for(Int_t yi = 1; yi <= h_8->GetNbinsY(); yi++)
        {
            h_8->SetBinContent(h_8->GetXaxis()->FindBin(event_time/1000.0),yi,h_8->GetBinContent(h_8->GetXaxis()->FindBin(event_time/1000.0),yi)/integral);
            h_8->SetBinError(h_8->GetXaxis()->FindBin(event_time/1000.0),yi,
                             TMath::Sqrt(TMath::Power(h_8->GetBinError(h_8->GetXaxis()->FindBin(event_time/1000.0),yi)/integral,2) +
                                         TMath::Power(h_8->GetBinContent(h_8->GetXaxis()->FindBin(event_time/1000.0),yi)*integral_err/(integral*integral),2)));
        }

        for(Int_t yi = 1; yi <= h_9->GetNbinsY(); yi++)
        {
            h_9->SetBinContent(h_9->GetXaxis()->FindBin(event_time/1000.0),yi,h_9->GetBinContent(h_9->GetXaxis()->FindBin(event_time/1000.0),yi)/integral);
            h_9->SetBinError(h_9->GetXaxis()->FindBin(event_time/1000.0),yi,
                             TMath::Sqrt(TMath::Power(h_9->GetBinError(h_9->GetXaxis()->FindBin(event_time/1000.0),yi)/integral,2) +
                                         TMath::Power(h_9->GetBinContent(h_9->GetXaxis()->FindBin(event_time/1000.0),yi)*integral_err/(integral*integral),2)));
        }

        // To save frames in the folder

        /*Int_t dtime = ((event_time - UT_min)/1000 + 0.5);
        TString name = "./frames_l6_chip1/rp1int_";
        name += i;
        TString title = "dTime: ";
        title += dtime;
        title += " [sec]";
        h_2->SetTitle(title.Data());
//        h_2->SetMaximum(1000);
        TCanvas* cc = new TCanvas(name.Data(),"cc",1000,1000);
        cc->cd();
        gPad->SetGrid();
        h_2->GetXaxis()->SetRange(1,256);
        h_2->GetYaxis()->SetRange(1,256);
        h_2->Draw("colz");
        name += ".png";
        cc->SaveAs(name.Data());*/
        h_1->Delete();
        h_2->Delete();
    }
    //---------------------------------------------------------------------------------------//

    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File '"<<_file->GetName()<<"' opened successfully!"<<endl;

    TH2D* h_rp1_6 = new TH2D("h_rp1_6","RP1 Internal (normalized)",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS*0.055);
    h_rp1_6->GetXaxis()->SetTitle("Time");
    h_rp1_6->GetXaxis()->SetTimeOffset(1.2);
    h_rp1_6->GetXaxis()->CenterTitle();
    h_rp1_6->GetYaxis()->SetTitle("Projection on Horizontal Axis [mm]");
    h_rp1_6->GetYaxis()->SetTimeOffset(1.2);
    h_rp1_6->GetYaxis()->CenterTitle();

    TH2D* h_rp1_7 = new TH2D("h_rp1_7","RP1 Internal (normalized)",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS*0.055);
    h_rp1_7->GetXaxis()->SetTitle("Time");
    h_rp1_7->GetXaxis()->SetTimeOffset(1.2);
    h_rp1_7->GetXaxis()->CenterTitle();
    h_rp1_7->GetYaxis()->SetTitle("Projection on Vertical Axis [mm]");
    h_rp1_7->GetYaxis()->SetTimeOffset(1.2);
    h_rp1_7->GetYaxis()->CenterTitle();

    Double_t integral, integral_err;
    integral = h_3->IntegralAndError(1,h_3->GetNbinsX(),1,h_3->GetNbinsY(),integral_err);
    scale2Dhisto(h_3,integral,integral_err);

    rewriteHisto(h_8,h_rp1_6);
    rewriteHisto(h_9,h_rp1_7);

    rotateNscale(h_3,h_4);

    TCanvas* c_0 = new TCanvas("c_0","c_0",900,900);
    c_0->cd();
    gPad->SetGrid();
    h_3->GetYaxis()->SetRange(1,256);
    h_3->GetXaxis()->SetRange(1,256);
    h_3->Draw("colz");

    TCanvas* c_00 = new TCanvas("c_00","c_00",900,900);
    c_00->cd();
    gPad->SetGrid();
    h_4->GetYaxis()->SetRange(1,256);
    h_4->GetXaxis()->SetRange(1,256);
    h_4->Draw("colz");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    h_8->GetYaxis()->SetRange(1,256);
    h_8->GetXaxis()->SetTimeDisplay(1);
    h_8->Draw("colz");

    TCanvas* c_11 = new TCanvas("c_11","c_11",1800,900);
    c_11->cd();
    gPad->SetGrid();
    h_9->GetYaxis()->SetRange(1,256);
    h_9->GetXaxis()->SetTimeDisplay(1);
    h_9->Draw("colz");

    TCanvas* c_2 = new TCanvas("c_2","c_2",1800,900);
    c_2->cd();
    gPad->SetGrid();
    h_rp1_6->GetYaxis()->SetRange(1,256);
    h_rp1_6->GetXaxis()->SetTimeDisplay(1);
    h_rp1_6->Draw("colz");

    TCanvas* c_22 = new TCanvas("c_22","c_22",1800,900);
    c_22->cd();
    gPad->SetGrid();
    h_rp1_7->GetYaxis()->SetRange(1,256);
    h_rp1_7->GetXaxis()->SetTimeDisplay(1);
    h_rp1_7->Draw("colz");

    c_0->Write();
    c_00->Write();
    c_1->Write();
    c_2->Write();
    h_3->Write();
    h_4->Write();
    h_8->Write();
    h_9->Write();
    h_5->Write();
    h_rp1_6->Write();
    h_rp1_7->Write();

    _file->Close();
    //---------------------------------------------------------------------------------------//
}

void function_2()
{
    gStyle->SetOptStat(0);
    TFile *_file_chipid_chip1 = TFile::Open("output_L1_chipid_1.root");
    TGraph *_cry3_angularscan = new TGraph("crystal3_angularscan.dat","%lg %lg");
    _cry3_angularscan->SetName("_cry3_angularscan");
    _cry3_angularscan->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    _cry3_angularscan->GetXaxis()->CenterTitle();
    _cry3_angularscan->GetXaxis()->SetTitle("Time");

    Double_t time_start = 1537238370;   // Tuesday, 18 September 2018, 04:39:30
    Double_t time_stop  = 1537238880;   // Tuesday, 18 September 2018, 04:48:00

//    to fit with splines
//    TSpline3 *_cry3_angularscan_spline = new TSpline3("_cry3_angularscan_spline",_cry3_angularscan);

//    to fit with pol1
    TF1 *_cry3_angularscan_spline = new TF1("_cry3_angularscan_spline","pol1",time_start,time_stop);
    _cry3_angularscan->Fit(_cry3_angularscan_spline,"R+");

    TH2D* h_chip1_1 = (TH2D*)_file_chipid_chip1->Get("h_rp1_6");
    TH2D* h_chip1_2 = new TH2D("h_chip1_HvsA","h_chip1_HvsA",600,1400,2000,N_PIXELS,0,N_PIXELS*0.055);

    for(Int_t i = 1; i <= h_chip1_1->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h_chip1_1->GetNbinsY(); j++)
        {
            Double_t _time      = h_chip1_1->GetXaxis()->GetBinCenter(i);
            Double_t _angle     = (-1)*_cry3_angularscan_spline->Eval(_time);
            Double_t _position  = h_chip1_1->GetYaxis()->GetBinCenter(j);
            Double_t _value     = h_chip1_1->GetBinContent(i,j);
            Double_t _value_err = h_chip1_1->GetBinError(i,j);

            if(_time > time_start && _time < time_stop)
            {
                h_chip1_2->SetBinContent(h_chip1_2->GetXaxis()->FindBin(_angle),h_chip1_2->GetYaxis()->FindBin(_position),_value);
                h_chip1_2->SetBinError(h_chip1_2->GetXaxis()->FindBin(_angle),h_chip1_2->GetYaxis()->FindBin(_position),_value_err);

//                h_chip1_2->Fill(_angle,_position,_value);
            }
        }
    }

    TH1D* h_chip1_2_projection = (TH1D*)h_chip1_2->ProjectionY("h_chip1_2_projection",h_chip1_2->GetXaxis()->FindBin(1760),h_chip1_2->GetXaxis()->FindBin(1880));
    TH1D* h_chip1_3_projection = (TH1D*)h_chip1_2->ProjectionX("h_chip1_3_projection",h_chip1_2->GetYaxis()->FindBin(1.5),h_chip1_2->GetYaxis()->FindBin(5.5));

    TCanvas* c_0 = new TCanvas("c_0","c_0",1000,1000);
    c_0->cd();
    gPad->SetGrid();
    _cry3_angularscan->SetMarkerStyle(21);
    _cry3_angularscan->Draw("AP");
    _cry3_angularscan_spline->SetLineColor(kRed);
    _cry3_angularscan_spline->SetLineWidth(2);
    _cry3_angularscan_spline->Draw("same & L");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    h_chip1_2->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(1760),h_chip1_2->GetXaxis()->FindBin(1880));
    h_chip1_2->SetTitle("");
    h_chip1_2->GetXaxis()->CenterTitle();
    h_chip1_2->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    h_chip1_2->GetYaxis()->CenterTitle();
    h_chip1_2->GetYaxis()->SetTitle("Projection on Horizontal Axis [mm]");
    h_chip1_2->Draw("colz");

    TCanvas* c_2 = new TCanvas("c_2","c_2",900,900);
    c_2->cd();
    gPad->SetGrid();
    h_chip1_2_projection->Draw("hist");

    TCanvas* c_3 = new TCanvas("c_3","c_3",1800,900);
    c_3->Divide(2,1);
    c_3->cd(1);
    gPad->SetGrid();
    h_chip1_3_projection->SetTitle("");
    h_chip1_3_projection->SetLineWidth(2);
    h_chip1_3_projection->SetFillColor(kCyan);
    h_chip1_3_projection->Draw("bar1 & e");
    h_chip1_3_projection->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(1760),h_chip1_2->GetXaxis()->FindBin(1880));

    TF1* fit_1 = new TF1("fit_1","gaus",1780,1860);
    h_chip1_3_projection->Fit(fit_1,"RQ0+");
    fit_1->SetParameter(1,1825);
    fit_1->SetLineStyle(7);
    fit_1->SetLineWidth(4);
    fit_1->SetLineColor(kBlack);
//    fit_1->Draw("same");

    TH1D* hh = new TH1D("hh","hh",h_chip1_3_projection->GetNbinsX(),h_chip1_3_projection->GetBinLowEdge(1),h_chip1_3_projection->GetBinLowEdge(h_chip1_3_projection->GetNbinsX())+h_chip1_3_projection->GetBinWidth(1));
    for(Int_t i = 1; i <= h_chip1_3_projection->GetNbinsX(); i++)
    {
        hh->SetBinContent(i,h_chip1_3_projection->GetBinContent(i));
        hh->SetBinError(i,h_chip1_3_projection->GetBinError(i));
        if(h_chip1_3_projection->GetBinContent(i) - 10.0*h_chip1_3_projection->GetBinError(i) > fit_1->Eval(h_chip1_3_projection->GetBinCenter(i)))
        {
            hh->SetBinContent(i,fit_1->Eval(h_chip1_3_projection->GetBinCenter(i)));
        }
        if(h_chip1_3_projection->GetBinContent(i) + 10.0*h_chip1_3_projection->GetBinError(i) < fit_1->Eval(h_chip1_3_projection->GetBinCenter(i)))
        {
            hh->SetBinContent(i,fit_1->Eval(h_chip1_3_projection->GetBinCenter(i)));
        }
    }

    hh->SetTitle("");
    hh->SetLineWidth(2);
    hh->SetLineColor(kGreen+2);
    hh->Draw("hist & same & e");
    h_chip1_3_projection->GetXaxis()->CenterTitle();
    h_chip1_3_projection->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");

    c_3->cd(2);
    hh->Draw("hist & e");
    hh->GetXaxis()->SetRange(hh->FindBin(1760),hh->FindBin(1880));
    TF1* fit_2 = new TF1("fit_2","gaus",1760,1880);
    hh->Fit(fit_2,"RQ0+");
    fit_2->SetLineStyle(7);
    fit_2->SetLineWidth(4);
    fit_2->SetLineColor(kBlack);
    fit_2->Draw("same");
    hh->GetXaxis()->CenterTitle();
    hh->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");

    cout<<"--> chi2/ndf = "<<fit_2->GetChisquare()<<"/"<<fit_2->GetNDF()<<endl;
    cout<<"--> Constant = "<<fit_2->GetParameter(0)<<" +/- "<<fit_2->GetParError(0)<<endl;
    cout<<"--> Mean = "<<fit_2->GetParameter(1)<<" +/- "<<fit_2->GetParError(1)<<endl;
    cout<<"--> Sigma = "<<fit_2->GetParameter(2)<<" +/- "<<fit_2->GetParError(2)<<endl;
}

void function_3()
{
    gStyle->SetOptStat(0);
    TFile *_file_chipid_chip1 = TFile::Open("output_L3_chipid_1.root");
    TGraph *_cry2_angularscan = new TGraph("crystal2_angularscan.dat","%lg %lg");
    _cry2_angularscan->SetName("_cry2_angularscan");
    _cry2_angularscan->SetTitle("CRYSTAL2 LVDT Angle [urad]");

    TSpline3 *_cry2_angularscan_spline = new TSpline3("_cry2_angularscan_spline",_cry2_angularscan);

    TH2D* h_chip1_1 = (TH2D*)_file_chipid_chip1->Get("h_rp1_6");
    TH2D* h_chip1_2 = new TH2D("h_chip1_HvsA","h_chip1_HvsA",7000,-6700,-6000,N_PIXELS,0,N_PIXELS*0.055);

    for(Int_t i = 1; i <= h_chip1_1->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h_chip1_1->GetNbinsY(); j++)
        {
            Double_t _time      = h_chip1_1->GetXaxis()->GetBinCenter(i);
            Double_t _angle     = _cry2_angularscan_spline->Eval(_time);
            Double_t _position  = h_chip1_1->GetYaxis()->GetBinCenter(j);
            Double_t _value     = h_chip1_1->GetBinContent(i,j);

            if(_angle > -6660 && _angle < -6060)
            {
                h_chip1_2->Fill(_angle,_position,_value);
            }
        }
    }

    TString outFileName     = "output_function_2.root";
    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File opened successfully!"<<endl;

    TCanvas* c_0 = new TCanvas("c_0","c_0",1000,1000);
    c_0->cd();
    gPad->SetGrid();
    _cry2_angularscan->SetMarkerStyle(21);
    _cry2_angularscan->Draw("AP");
    _cry2_angularscan_spline->SetLineColor(kRed);
    _cry2_angularscan_spline->SetLineWidth(2);
    _cry2_angularscan_spline->Draw("same & L");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    h_chip1_2->Draw("colz");

    c_0->Write();
    _cry2_angularscan->Write();
    _file->Close();
}

void function_4()
{
    gStyle->SetOptStat(0);
    TString fileName_RP1 = "output_L2_chipid_1.root";

    Double_t fit_bkg_lim_min = 12.7;
    Double_t fit_bkg_lim_max = 13.2;
    Double_t fit_doublech_lim_min = 10.5;
    Double_t fit_doublech_lim_max = 12.2;
    Double_t fit_singlech_lim_min =  2.0;
    Double_t fit_singlech_lim_max =  4.5;
    Double_t par_bkg[2];

    TFile *_file_rp1 = TFile::Open(fileName_RP1.Data());
    TH2D* h_rp1 = (TH2D*)_file_rp1->Get("h_4");

    cout<<"--> Integral RP1: "<<h_rp1->Integral()<<endl;

    TCanvas* c_1 = new TCanvas("c_1","c_1",1000,1000);
    c_1->cd();
    h_rp1->Draw("colz");

    TH1D* h_rp1_x       = h_rp1->ProjectionX("h_rp1_x");
    TH1D* h_rp1_x_s     = h_rp1->ProjectionX("h_rp1_x_s");
    TH1D* h_rp1_x_s_s   = h_rp1->ProjectionX("h_rp1_x_s_s");

    TCanvas* c_2 = new TCanvas("c_2","c_2",1000,500);
    c_2->cd();
    h_rp1_x->SetMinimum(0);
    h_rp1_x->Draw("hist");
    TF1* fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","pol0",fit_bkg_lim_min,fit_bkg_lim_max);
    fit_funct_ch_bkg->SetParameter(0, 2e-3);
    h_rp1_x->Fit(fit_funct_ch_bkg,"R+");
    fit_funct_ch_bkg->GetParameters(&par_bkg[0]);

    TString name = "#chi^{2}/NDF = [";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare();
    name += "/";
    name += fit_funct_ch_bkg->GetNDF();
    name += "] = ";
    name += (Int_t)fit_funct_ch_bkg->GetChisquare()/fit_funct_ch_bkg->GetNDF();
    h_rp1_x->SetTitle(name.Data());
    fit_funct_ch_bkg->SetLineWidth(6);
    fit_funct_ch_bkg->Draw("same");

    TF1* bkg_func = new TF1("bkg_func","pol0",h_rp1_x->GetBinCenter(1),h_rp1_x->GetBinCenter(h_rp1_x->GetNbinsX()));
    bkg_func->SetParameters(par_bkg);
    bkg_func->SetLineColor(kGreen+2);
    bkg_func->Draw("same");

    TCanvas* c_3 = new TCanvas("c_3","c_3",1000,500);
    c_3->cd();
    h_rp1_x_s->SetMinimum(0);
    h_rp1_x_s->SetLineColor(kBlue);
    h_rp1_x_s->SetLineWidth(2);
    h_rp1_x_s->Draw("hist");
    for(Int_t j = 1; j <= h_rp1_x_s_s->GetNbinsX(); j++)
    {
        Double_t val = h_rp1_x_s->GetBinContent(j) - bkg_func->Eval(h_rp1_x_s->GetBinCenter(j));

        if(val > 0)
            h_rp1_x_s_s->SetBinContent(j,val);
        else
            h_rp1_x_s_s->SetBinContent(j,0);
    }
    h_rp1_x_s_s->SetLineColor(kRed);
    h_rp1_x_s_s->SetLineWidth(2);
    h_rp1_x_s_s->Draw("same & hist");

    TCanvas* c_4 = new TCanvas("c_4","c_4",1000,500);
    c_4->cd();
    gPad->SetGrid();
    h_rp1_x_s_s->SetLineWidth(2);
    h_rp1_x_s_s->SetTitle("");
    h_rp1_x_s_s->GetYaxis()->SetTitle("Normalized Counts");
    h_rp1_x_s_s->GetXaxis()-> CenterTitle();
    h_rp1_x_s_s->GetYaxis()-> CenterTitle();
    h_rp1_x_s_s->Draw("hist");

    TF1* fit_funct_1 = new TF1("fit_funct_1","gaus",fit_singlech_lim_min,fit_singlech_lim_max);
    h_rp1_x_s_s->Fit(fit_funct_1,"R+");
    fit_funct_1->SetLineColor(kBlue);
    fit_funct_1->SetLineWidth(4);
    fit_funct_1->Draw("same");

    TF1* fit_funct_2 = new TF1("fit_funct_2","gaus",fit_doublech_lim_min,fit_doublech_lim_max);
    h_rp1_x_s_s->Fit(fit_funct_2,"R+");
    fit_funct_2->SetLineColor(kBlack);
    fit_funct_2->SetLineWidth(4);
    fit_funct_2->Draw("same");

    TLine *l_min_ch=new TLine(fit_funct_2->GetParameter(1)-3.0*fit_funct_2->GetParameter(2),0,fit_funct_2->GetParameter(1)-3.0*fit_funct_2->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_min_ch->SetLineColor(kMagenta);
    l_min_ch->SetLineWidth(3);
    l_min_ch->SetLineStyle(7);
    l_min_ch->Draw();
    TLine *l_mean_ch=new TLine(fit_funct_2->GetParameter(1),0,fit_funct_2->GetParameter(1),h_rp1_x_s_s->GetMaximum());
    l_mean_ch->SetLineColor(kMagenta-2);
    l_mean_ch->SetLineWidth(3);
    l_mean_ch->SetLineStyle(7);
    l_mean_ch->Draw();
    TLine *l_max_ch=new TLine(fit_funct_2->GetParameter(1)+3.0*fit_funct_2->GetParameter(2),0,fit_funct_2->GetParameter(1)+3.0*fit_funct_2->GetParameter(2),h_rp1_x_s_s->GetMaximum());
    l_max_ch->SetLineColor(kMagenta);
    l_max_ch->SetLineWidth(3);
    l_max_ch->SetLineStyle(7);
    l_max_ch->Draw();

    Double_t integral_total, integral_total_err;
    Double_t integral_doublech, integral_doublech_err;

    integral_total = h_rp1_x_s_s->IntegralAndError(1,h_rp1_x_s_s->GetNbinsX(),integral_total_err);
    integral_doublech = h_rp1_x_s_s->IntegralAndError(h_rp1_x_s_s->FindBin(fit_funct_2->GetParameter(1)-3.0*fit_funct_2->GetParameter(2)),h_rp1_x_s_s->FindBin(fit_funct_2->GetParameter(1)+3.0*fit_funct_2->GetParameter(2)),integral_doublech_err);

    cout<<"--> integral_total = "<<integral_total<<" +/- "<<integral_total_err<<endl;
    cout<<"--> integral_doublech = "<<integral_doublech<<" +/- "<<integral_doublech_err<<endl;
    cout<<"--> integral_doublech/integral_total = "<<(Double_t)integral_doublech/integral_total<<" +/- "<<
          TMath::Sqrt(TMath::Power(integral_doublech_err/integral_total,2) + TMath::Power(integral_doublech*integral_total_err/(integral_total*integral_total),2))<<endl;

}

void function_5()
{
//    gStyle->SetOptStat(0);
    TString fileName_RP1 = "output_L2_chipid_1.root";

    Double_t fit_bkg_lim_max = 256*0.055-12.7;
    Double_t fit_bkg_lim_min = 256*0.055-13.2;
    Double_t fit_doublech_lim_max = 256*0.055-(1.13274e+01-3.0*3.16233e-01);
    Double_t fit_doublech_lim_min = 256*0.055-(1.13274e+01+3.0*3.16233e-01);
    Double_t par_bkg[1];

    TFile *_file_rp1 = TFile::Open(fileName_RP1.Data());
    TH2D* h_rp1 = (TH2D*)_file_rp1->Get("h_4");
    TH2D* h_rp1_x = (TH2D*)_file_rp1->Get("h_rp1_6");

    TCanvas* c_00 = new TCanvas("c_00","c_00",1800,900);
    c_00->cd();
    h_rp1_x->Draw("colz");

    TCanvas* c_22 = new TCanvas("c_22","c_22",1000,1000);
    c_22->cd();
    h_rp1->Draw("colz");

    TF1* fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","pol0",fit_bkg_lim_min,fit_bkg_lim_max);
    TF1* bkg_func = new TF1("bkg_func","pol0",h_rp1_x->GetYaxis()->GetBinCenter(1),h_rp1_x->GetYaxis()->GetBinCenter(h_rp1_x->GetNbinsY()));

    Double_t integral_total, integral_total_err;
    Double_t integral_doublech, integral_doublech_err;
    vector<Double_t> ratio;
    vector<Double_t> ratio_err;
    vector<Double_t> unix_time;
    vector<Double_t> unix_time_err;
    vector<Double_t> chi2ndf;

    TH1D* ratio_histo = new TH1D("ratio_histo","ratio_histo",16,0.04,0.12);
    TH1D* chi2ndf_histo = new TH1D("chi2ndf_histo","chi2ndf_histo",200,-5,15);

    Int_t notcutted_slices = 0, total_slices = 0;

    for(Int_t i = 1; i <= h_rp1_x->GetNbinsX(); i++)
    {
        total_slices++;
        TString name = "h_x_bin";
        name += i;
        TH1D* h_x   = h_rp1_x->ProjectionY(name.Data(),i,i);

        fit_funct_ch_bkg->SetParameter(0, 2e-3);
        h_x->Fit(fit_funct_ch_bkg,"R0Q+");
        par_bkg[0] = fit_funct_ch_bkg->GetParameter(0);
        bkg_func->SetParameter(0,par_bkg[0]);

        for(Int_t j = 1; j <= h_x->GetNbinsX(); j++)
        {
            Double_t val = h_x->GetBinContent(j) - bkg_func->Eval(h_x->GetBinCenter(j));

            if(val > 0)
                h_x->SetBinContent(j,val);
            else
                h_x->SetBinContent(j,0);
        }

        integral_total = h_x->IntegralAndError(1,h_x->GetNbinsX(),integral_total_err);
        integral_doublech = h_x->IntegralAndError(h_x->FindBin(fit_doublech_lim_min),h_x->FindBin(fit_doublech_lim_max),integral_doublech_err);

        if(integral_doublech/integral_total < 0.12)
        {
            notcutted_slices++;
            unix_time.push_back(h_rp1_x->GetXaxis()->GetBinCenter(i));
            unix_time_err.push_back(h_rp1_x->GetXaxis()->GetBinWidth(i)/TMath::Sqrt(12.0));
            chi2ndf.push_back(fit_funct_ch_bkg->GetChisquare()/fit_funct_ch_bkg->GetNDF());
            chi2ndf_histo->Fill(fit_funct_ch_bkg->GetChisquare()/fit_funct_ch_bkg->GetNDF());
            ratio.push_back(integral_doublech/integral_total);
            ratio_histo->Fill(integral_doublech/integral_total);
            ratio_err.push_back(TMath::Sqrt(TMath::Power(integral_doublech_err/integral_total,2) + TMath::Power(integral_doublech*integral_total_err/(integral_total*integral_total),2)));
        }
    }

    cout<<endl<<"--> Cutted slices = "<<total_slices-notcutted_slices<<" ("<<(Int_t)((total_slices-notcutted_slices)*100.0/total_slices + 0.5)<<"%)"<<endl;

    Int_t dim = unix_time.size();

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,450);
    c_1->cd();
    gPad->SetGrid();
    TGraphErrors *gr_1 = new TGraphErrors(dim, &unix_time[0], &chi2ndf[0], &unix_time_err[0], 0);
    gr_1->SetMarkerStyle(21);
    gr_1->SetMarkerColor(kBlack);
    gr_1->SetLineWidth(2);
    gr_1->SetLineColor(kBlue);
    gr_1->SetMaximum(15);
    gr_1->SetMinimum(-5);
    gr_1->SetTitle("");
    gr_1->GetXaxis()->SetTitle("Time");
    gr_1->GetYaxis()->SetTitle("Background #chi^{2}/NDF");
    gr_1->GetXaxis()->SetTimeDisplay(1);
    gr_1->GetYaxis()->CenterTitle();
    gr_1->GetXaxis()->CenterTitle();
    gr_1->Draw("AP");

    TCanvas* c_11 = new TCanvas("c_11","c_11",450,450);
    c_11->cd();
    gPad->SetGrid();
    chi2ndf_histo->SetLineWidth(2);
    chi2ndf_histo->Draw("hist");

    TCanvas* c_0 = new TCanvas("c_0","c_0",1800,450);
    c_0->cd();
    gPad->SetGrid();
    TGraphErrors *gr_0 = new TGraphErrors(dim, &unix_time[0], &ratio[0], &unix_time_err[0], &ratio_err[0]);
    gr_0->SetMarkerStyle(21);
    gr_0->SetMarkerColor(kBlack);
    gr_0->SetLineWidth(2);
    gr_0->SetLineColor(kBlue);
    gr_0->SetMaximum(0.12);
    gr_0->SetMinimum(0.04);
    gr_0->SetTitle("");
    gr_0->GetXaxis()->SetTitle("Time");
    gr_0->GetYaxis()->SetTitle("I2/I1");
    gr_0->GetXaxis()->SetTimeDisplay(1);
    gr_0->GetYaxis()->CenterTitle();
    gr_0->GetXaxis()->CenterTitle();
    gr_0->Draw("AP");

    TCanvas* c_01 = new TCanvas("c_01","c_01",900,900);
    c_01->cd();
    gPad->SetGrid();
    ratio_histo->SetLineWidth(2);
    ratio_histo->Draw("hist");
}

void function_6()
{
    TString fileName_RP1 = "output_L1_chipid_1.root";
    TFile *_file_rp1 = TFile::Open(fileName_RP1.Data());
    TH2D* h_rp1_6 = (TH2D*)_file_rp1->Get("h_rp1_6");
    Double_t start_time = 1537238620; // Tuesday, 18 September 2018, 04:43:40
    Double_t stop_time  = 1537238700; // Tuesday, 18 September 2018, 04:45:00

    TH1D* h_x = (TH1D*)h_rp1_6->ProjectionY("h_x",h_rp1_6->GetXaxis()->FindBin(start_time),h_rp1_6->GetXaxis()->FindBin(stop_time));
    TH1D* h_x_x = (TH1D*)h_rp1_6->ProjectionY("h_x_x",h_rp1_6->GetXaxis()->FindBin(start_time),h_rp1_6->GetXaxis()->FindBin(stop_time));

    TCanvas* c_0 = new TCanvas("c_0","c_0",1800,900);
    c_0->cd();
    gPad->SetGrid();
    h_rp1_6->Draw("colz");
    h_rp1_6->GetXaxis()->SetRange(h_rp1_6->GetXaxis()->FindBin(start_time),h_rp1_6->GetXaxis()->FindBin(stop_time));

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    h_x->Draw();
    h_x->SetLineColor(kBlue);
    h_x->SetLineWidth(2);
    TF1* g_BK = new TF1("g_BK","pol0",1.00,4.00);
    h_x->Fit(g_BK,"R+");

    for(Int_t i = 1; i <= h_x->GetNbinsX(); i++)
    {
        Double_t val = h_x->GetBinContent(i) - g_BK->Eval(h_x->GetBinCenter(i));
        if(val > 0)
            h_x_x->SetBinContent(i,val);
        else
            h_x_x->SetBinContent(i,0);
    }

    h_x_x->SetLineColor(kRed);
    h_x_x->SetLineWidth(2);
    h_x->SetMinimum(0.0);
    h_x_x->Draw("same");

    TCanvas* c_2 = new TCanvas("c_2","c_2",1800,900);
    c_2->cd();
    gPad->SetGrid();
    h_x_x->Draw();

    TF1* g_AM = new TF1("g_AM","gaus",9.90,11.44);
    TF1* g_VR = new TF1("g_VR","gaus",11.88,13.64);
    g_AM->SetLineColor(kBlue);
    g_VR->SetLineColor(kBlue);

    h_x_x->Fit(g_AM,"RQ0+");
    h_x_x->Fit(g_VR,"RQ0+");

    Double_t g_AM_VR_par[6];
    Double_t g_AM_VR_par_err[6];
    g_AM->GetParameters(&g_AM_VR_par[0]);
    g_VR->GetParameters(&g_AM_VR_par[3]);
    TF1* g_AM_VR = new TF1("g_AM_VR_BK","gaus(0)+gaus(3)",10.0,13.5);
    g_AM_VR->SetParameters(g_AM_VR_par);
    g_AM_VR->SetLineColor(kBlack);
    h_x_x->Fit(g_AM_VR,"R+");
    g_AM_VR->GetParameters(g_AM_VR_par);
    g_AM_VR_par_err[0] = g_AM_VR->GetParError(0);
    g_AM_VR_par_err[1] = g_AM_VR->GetParError(1);
    g_AM_VR_par_err[2] = g_AM_VR->GetParError(2);
    g_AM_VR_par_err[3] = g_AM_VR->GetParError(3);
    g_AM_VR_par_err[4] = g_AM_VR->GetParError(4);
    g_AM_VR_par_err[5] = g_AM_VR->GetParError(5);
    h_x_x->SetMinimum(0.0);

    TF1* g_AM_f = new TF1("g_AM_f","gaus",0.0,0.5*N_PIXELS*0.055);
    TF1* g_VR_f = new TF1("g_VR_f","gaus",0.0,0.5*N_PIXELS*0.055);
    g_AM_f->SetParameters(&g_AM_VR_par[0]);
    g_AM_f->SetParErrors(&g_AM_VR_par_err[0]);
    g_VR_f->SetParameters(&g_AM_VR_par[3]);
    g_VR_f->SetParErrors(&g_AM_VR_par_err[3]);
    g_AM_f->SetLineColor(kMagenta);
    g_VR_f->SetLineColor(kBlue);
    g_AM_f->SetLineWidth(2);
    g_VR_f->SetLineWidth(2);
    g_AM_f->Draw("same");
    g_VR_f->Draw("same");

    Double_t integral_total, integral_total_err;
    Double_t integral_vr = 0;
    integral_total = h_x_x->IntegralAndError(1,h_x_x->GetNbinsX(),integral_total_err)*h_x_x->GetBinWidth(1);
    integral_total_err = TMath::Sqrt(TMath::Power(integral_total_err,2) + TMath::Power(h_x_x->GetBinWidth(1)/TMath::Sqrt(12.0),2));
    TGraph* gr = new TGraph();
    for(Int_t i = 0; i < 20000; i++)
    {
        integral_vr += g_VR_f->Eval(i*(0.5*N_PIXELS*0.055)/10000.0)*(0.5*N_PIXELS*0.055)/10000.0;
        gr->SetPoint(i,i*(0.5*N_PIXELS*0.055)/10000.0,g_VR_f->Eval(i*(0.5*N_PIXELS*0.055)/10000.0));
    }
    gr->SetLineStyle(7);
    gr->SetLineWidth(4);
    gr->Draw("same");
    cout<<"--> Chi2/NDF = "<<g_AM_VR->GetChisquare()<<"/"<<g_AM_VR->GetNDF()<<endl;
    cout<<"--> I1 = "<<integral_total<<" +/- "<<integral_total_err<<endl;
    cout<<"--> I3 = "<<integral_vr<<endl;
    cout<<"--> I3/I1 = "<<integral_vr/integral_total<<" +/- "<<integral_vr*integral_total_err/(integral_total*integral_total)<<endl;
}

void function_7()
{
//    gStyle->SetOptStat(0);
    TFile *_file_chipid_chip1 = TFile::Open("output_L5_chipid_1.root");
    TGraph *_cry3_linearscan = new TGraph("crystal3_linearscan_fill2.dat","%lg %lg");
    _cry3_linearscan->SetName("_cry3_linearscan");
    _cry3_linearscan->SetTitle("CRYSTAL3 LVDT Position [mm]");
    _cry3_linearscan->GetXaxis()->CenterTitle();
    _cry3_linearscan->GetXaxis()->SetTitle("Time");

    TCanvas* ccc_000 = new TCanvas();
    ccc_000->cd();
    _cry3_linearscan->Draw("APL");


    Double_t time_start = 1537218040;   // Tuesday, 17 September 2018, 23:00:40
    Double_t time_stop  = 1537218100;   // Tuesday, 17 September 2018, 23:01:40

    Double_t x_ini,y_ini,x_fin,y_fin;
    Int_t nPoints = _cry3_linearscan->GetN();
    _cry3_linearscan->GetPoint(0, x_ini, y_ini);
    _cry3_linearscan->GetPoint(nPoints-1, x_fin, y_fin);
    TH2* h = new TH2D("h","h",nPoints,x_ini,x_fin,nPoints,_cry3_linearscan->GetMinimum(),_cry3_linearscan->GetMaximum());
    for(Int_t i = 0; i < nPoints; i++)
    {
        Double_t x,y;
        _cry3_linearscan->GetPoint(i,x,y);
        h->Fill(x,y);
    }
    TProfile* pfy = h->ProfileY("_pfy",h->GetXaxis()->FindBin(time_start),h->GetXaxis()->FindBin(time_stop));


//    to fit with splines
//    TSpline3 *_cry3_linearscan_spline = new TSpline3("_cry3_linearscan_spline",_cry3_linearscan);

//    to fit with pol1
    TF1 *_cry3_linearscan_spline = new TF1("_cry3_linearscan_spline","pol1",60,80);
    pfy->Fit(_cry3_linearscan_spline,"RQ0+");

    TH2D* h_chip1_1 = (TH2D*)_file_chipid_chip1->Get("h_rp1_6");
    TH2D* h_chip1_2 = new TH2D("h_chip1_HvsA","h_chip1_HvsA",200,63,83,N_PIXELS,0,N_PIXELS*0.055);

    Double_t slope = _cry3_linearscan_spline->GetParameter(1);
    Double_t offset = _cry3_linearscan_spline->GetParameter(0);

    for(Int_t i = 1; i <= h_chip1_1->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h_chip1_1->GetNbinsY(); j++)
        {
            Double_t _time      = h_chip1_1->GetXaxis()->GetBinCenter(i);
            Double_t _linpos    = (_time - offset)/slope;
            Double_t _position  = h_chip1_1->GetYaxis()->GetBinCenter(j);
            Double_t _value     = h_chip1_1->GetBinContent(i,j);
            Double_t _value_err = h_chip1_1->GetBinError(i,j);

            if(_time > time_start && _time < time_stop)
            {
                h_chip1_2->SetBinContent(h_chip1_2->GetXaxis()->FindBin(_linpos),h_chip1_2->GetYaxis()->FindBin(_position),_value);
                h_chip1_2->SetBinError(h_chip1_2->GetXaxis()->FindBin(_linpos),h_chip1_2->GetYaxis()->FindBin(_position),_value_err);

//                h_chip1_2->Fill(_linpos,_position,_value);
            }
        }
    }

    TH1D* h_chip1_2_projection = (TH1D*)h_chip1_2->ProjectionY("h_chip1_2_projection",h_chip1_2->GetXaxis()->FindBin(63),h_chip1_2->GetXaxis()->FindBin(83));
    TH1D* h_chip1_3_projection = (TH1D*)h_chip1_2->ProjectionX("h_chip1_3_projection",h_chip1_2->GetYaxis()->FindBin(1.5),h_chip1_2->GetYaxis()->FindBin(5.5));
    TH1D* h_pr_y_wo_cry = (TH1D*)h_chip1_2->ProjectionY("h_pr_y_wo_cry",h_chip1_2->GetXaxis()->FindBin(76),h_chip1_2->GetXaxis()->FindBin(78));
    TH1D* h_pr_y_wo_cry_n_bkg = (TH1D*)h_chip1_2->ProjectionY("h_pr_y_wo_cry_n_bkg",h_chip1_2->GetXaxis()->FindBin(76),h_chip1_2->GetXaxis()->FindBin(78));

    TString outFileName     = "output_function_2.root";
    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File opened successfully!"<<endl;

    TCanvas* c_0 = new TCanvas("c_0","c_0",1000,1000);
    c_0->cd();
    gPad->SetGrid();
    pfy->SetMarkerStyle(21);
    pfy->SetMinimum(time_start-100);
    pfy->SetMaximum(time_stop+100);
    pfy->Draw();
    _cry3_linearscan_spline->SetLineColor(kRed);
    _cry3_linearscan_spline->SetLineWidth(2);
    _cry3_linearscan_spline->Draw("same & L");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    h_chip1_2->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(63),h_chip1_2->GetXaxis()->FindBin(83));
    h_chip1_2->SetTitle("");
    h_chip1_2->GetXaxis()->CenterTitle();
    h_chip1_2->GetXaxis()->SetTitle("CRYSTAL3 LVDT Position [mm]");
    h_chip1_2->GetYaxis()->CenterTitle();
    h_chip1_2->GetYaxis()->SetTitle("Projection on Horizontal Axis [mm]");
    h_chip1_2->Draw("colz");

    TCanvas* c_2 = new TCanvas("c_2","c_2",900,900);
    c_2->cd();
    gPad->SetGrid();
    h_chip1_2_projection->Draw("hist");

    TCanvas* c_3 = new TCanvas("c_3","c_3",900,900);
    c_3->cd();
    gPad->SetGrid();
    h_chip1_3_projection->Draw("hist");

    const Int_t nPnts = 5;
    Double_t mpx_pos[] = {7.3,  7.5,    11.0,   12.0,   12.5};
    Double_t mpx_err[] = {0.1,  0.1,    0.1,    0.1,    0.1};
    Double_t cry_pos[] = {64.8, 65.0,   69.7,   71.1,   71.8};
    Double_t cry_err[] = {0.1,  0.1,    0.1,    0.1,    0.1};

    TCanvas* c_4 = new TCanvas("c_4","c_4",1000,1000);
    c_4->cd();
    gPad->SetGrid();
    TGraphErrors* gr_1 = new TGraphErrors(nPnts,mpx_pos,cry_pos,mpx_err,cry_err);
    gr_1->SetMarkerStyle(21);
    gr_1->SetLineWidth(2);
    gr_1->SetLineColor(kBlue);
    gr_1->GetXaxis()->SetTitle("Position on Timepix [mm]");
    gr_1->GetYaxis()->SetTitle("CRYSTAL3 LVDT Position [mm]");
    gr_1->GetXaxis()->CenterTitle();
    gr_1->GetYaxis()->CenterTitle();
    gr_1->Draw("AP");

    TF1* fit_1 = new TF1("fit_1","pol1",mpx_pos[0]-1,mpx_pos[nPnts-1]+1);
    gr_1->Fit(fit_1,"R+");
    fit_1->SetLineWidth(3);

    TCanvas* c_5 = new TCanvas("c_5","c_5",1800,900);
    c_5->Divide(1,2);
    c_5->cd(1);
    h_pr_y_wo_cry->Draw();
    Double_t par[5];
    TF1* fit_pol = new TF1("fit_pol","pol1",0,4);
    TF1* fit_gaus = new TF1("fit_gaus","gaus",6.5,8.5);
    TF1* fit_pol_gaus = new TF1("fit_pol_gaus","pol1(0) + gaus(2)",0.0,8.5);
    h_pr_y_wo_cry->Fit(fit_pol,"RQ0+");
    h_pr_y_wo_cry->Fit(fit_gaus,"RQ0+");
    fit_pol->GetParameters(&par[0]);
    fit_gaus->GetParameters(&par[2]);
    fit_pol_gaus->SetParameters(par);
    h_pr_y_wo_cry->Fit(fit_pol_gaus,"R+");
    fit_pol_gaus->GetParameters(&par[0]);
    TF1* func_bkg = new TF1("func_bkg","pol1",0,h_pr_y_wo_cry->GetBinCenter(h_pr_y_wo_cry->GetNbinsX()));
    func_bkg->SetParameters(&par[0]);
    func_bkg->Draw("same");

    for(Int_t i = 1; i <= h_pr_y_wo_cry_n_bkg->GetNbinsX(); i++)
    {
        Double_t val = h_pr_y_wo_cry_n_bkg->GetBinContent(i) - func_bkg->Eval(h_pr_y_wo_cry_n_bkg->GetBinCenter(i));
        if(val < 0) val = 0;
        h_pr_y_wo_cry_n_bkg->SetBinContent(i,val);
    }
    c_5->cd(2);
    h_pr_y_wo_cry_n_bkg->Draw();

    TGraphErrors* gr_projection = new TGraphErrors();
    for(Int_t i = 1; i <= h_pr_y_wo_cry_n_bkg->GetNbinsX(); i++)
    {
        Double_t x = fit_1->Eval(h_pr_y_wo_cry_n_bkg->GetBinCenter(i));
        Double_t y = h_pr_y_wo_cry_n_bkg->GetBinContent(i);
        Double_t ey = h_pr_y_wo_cry_n_bkg->GetBinError(i);
        Double_t ex = 0.055/TMath::Sqrt(12.0);
        gr_projection->SetPoint(i-1,x,y);
        gr_projection->SetPointError(i-1,ex,ey);
    }

    TCanvas* c_6 = new TCanvas("c_6","c_6",1800,900);
    c_6->cd();
    gPad->SetGrid();
    gr_projection->SetLineColor(kBlue);
    gr_projection->SetLineWidth(2);
    gr_projection->Draw("AP");


    c_0->Write();
    c_1->Write();
    c_2->Write();
    c_3->Write();
    c_4->Write();
    c_5->Write();
    c_6->Write();
    h_pr_y_wo_cry->Write();
    gr_1->Write();
    gr_projection->Write();
    _cry3_linearscan->Write();
    _file->Close();
}

void function_8()
{
//    gStyle->SetOptStat(0);
    TFile *_file_chipid_chip1 = TFile::Open("output_L6_chipid_1.root");
    TGraph *_cry3_linearscan = new TGraph("crystal3_linearscan_fill3.dat","%lg %lg");
    _cry3_linearscan->SetName("_cry3_linearscan");
    _cry3_linearscan->SetTitle("CRYSTAL3 LVDT Position [mm]");
    _cry3_linearscan->GetXaxis()->CenterTitle();
    _cry3_linearscan->GetXaxis()->SetTitle("Time");


    Double_t time_start = 1537231066;   // Tuesday, 18 September 2018, 02:37:46
    Double_t time_stop  = 1537231105;   // Tuesday, 18 September 2018, 02:38:25

    Double_t x_ini,y_ini,x_fin,y_fin;
    Int_t nPoints = _cry3_linearscan->GetN();
    _cry3_linearscan->GetPoint(0, x_ini, y_ini);
    _cry3_linearscan->GetPoint(nPoints-1, x_fin, y_fin);
    TH2* h = new TH2D("h","h",nPoints,x_ini,x_fin,nPoints,_cry3_linearscan->GetMinimum(),_cry3_linearscan->GetMaximum());
    for(Int_t i = 0; i < nPoints; i++)
    {
        Double_t x,y;
        _cry3_linearscan->GetPoint(i,x,y);
        h->Fill(x,y);
    }
    TProfile* pfy = h->ProfileY("_pfy",h->GetXaxis()->FindBin(time_start),h->GetXaxis()->FindBin(time_stop));

//    to fit with splines
//    TSpline3 *_cry3_linearscan_spline = new TSpline3("_cry3_linearscan_spline",_cry3_linearscan);

//    to fit with pol1
    TF1 *_cry3_linearscan_spline = new TF1("_cry3_linearscan_spline","pol1",60,70);
    pfy->Fit(_cry3_linearscan_spline,"R+");

    TH2D* h_chip1_1 = (TH2D*)_file_chipid_chip1->Get("h_rp1_6");
    TH2D* h_chip1_2 = new TH2D("h_chip1_HvsA","h_chip1_HvsA",200,60,70,N_PIXELS,0,N_PIXELS*0.055);

    Double_t slope = _cry3_linearscan_spline->GetParameter(1);
    Double_t offset = _cry3_linearscan_spline->GetParameter(0);

    for(Int_t i = 1; i <= h_chip1_1->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h_chip1_1->GetNbinsY(); j++)
        {
            Double_t _time      = h_chip1_1->GetXaxis()->GetBinCenter(i);
            Double_t _linpos    = (_time - offset)/slope;
            Double_t _position  = h_chip1_1->GetYaxis()->GetBinCenter(j);
            Double_t _value     = h_chip1_1->GetBinContent(i,j);
            Double_t _value_err = h_chip1_1->GetBinError(i,j);

            if(_time > time_start && _time < time_stop)
            {
                h_chip1_2->SetBinContent(h_chip1_2->GetXaxis()->FindBin(_linpos),h_chip1_2->GetYaxis()->FindBin(_position),_value);
                h_chip1_2->SetBinError(h_chip1_2->GetXaxis()->FindBin(_linpos),h_chip1_2->GetYaxis()->FindBin(_position),_value_err);

//                h_chip1_2->Fill(_linpos,_position,_value);
            }
        }
    }

    TH1D* h_chip1_2_projection = (TH1D*)h_chip1_2->ProjectionY("h_chip1_2_projection",h_chip1_2->GetXaxis()->FindBin(63),h_chip1_2->GetXaxis()->FindBin(83));
    TH1D* h_chip1_3_projection = (TH1D*)h_chip1_2->ProjectionX("h_chip1_3_projection",h_chip1_2->GetYaxis()->FindBin(1.5),h_chip1_2->GetYaxis()->FindBin(5.5));

    TString outFileName     = "output_function_2.root";
    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File opened successfully!"<<endl;

    TCanvas* c_0 = new TCanvas("c_0","c_0",1000,1000);
    c_0->cd();
    gPad->SetGrid();
    pfy->SetMarkerStyle(21);
    pfy->SetMinimum(time_start-100);
    pfy->SetMaximum(time_stop+100);
    pfy->Draw();
    _cry3_linearscan_spline->SetLineColor(kRed);
    _cry3_linearscan_spline->SetLineWidth(2);
    _cry3_linearscan_spline->Draw("same & L");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    h_chip1_2->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(63),h_chip1_2->GetXaxis()->FindBin(83));
    h_chip1_2->SetTitle("");
    h_chip1_2->GetXaxis()->CenterTitle();
    h_chip1_2->GetXaxis()->SetTitle("CRYSTAL3 LVDT Position [mm]");
    h_chip1_2->GetYaxis()->CenterTitle();
    h_chip1_2->GetYaxis()->SetTitle("Projection on Horizontal Axis [mm]");
    h_chip1_2->Draw("colz");

    TCanvas* c_2 = new TCanvas("c_2","c_2",900,900);
    c_2->cd();
    gPad->SetGrid();
    h_chip1_2_projection->Draw("hist");

    TCanvas* c_3 = new TCanvas("c_3","c_3",900,900);
    c_3->cd();
    gPad->SetGrid();
    h_chip1_3_projection->Draw("hist");
//    h_chip1_3_projection->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(1760),h_chip1_2->GetXaxis()->FindBin(1880));

    c_0->Write();
    _cry3_linearscan->Write();
    _file->Close();
}

void function_9()
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

    //--- SPS FILL 2 ---//
//    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L5_RUN_8.root";
//    Double_t start_time = 1537218030, stop_time = 1537218040; // w/o CRYSTAL3
//    Double_t start_time = 1537218052, stop_time = 1537218059; // with CRYSTAL3
    //------------------//


    //--- SPS FILL 3 ---//
    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L6_RUN_8.root";
//    Double_t start_time = 1537231040, stop_time = 1537231060; // w/o CRYSTAL3
    Double_t start_time = 1537231094, stop_time = 1537231100; // with CRYSTAL3
    //------------------//

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

    TH1D* h_py_pixels = new TH1D("h_py_pixels","Projection on Horizontal Axis (pixels)",N_PIXELS,0,N_PIXELS);
    h_py_pixels->GetXaxis()->SetTitle("Horizontal Axis [pixels]");
    h_py_pixels->GetXaxis()->SetTimeOffset(1.2);
    h_py_pixels->GetXaxis()->CenterTitle();
    h_py_pixels->GetYaxis()->SetTitle("Normalized Counts");
    h_py_pixels->GetYaxis()->SetTimeOffset(1.2);
    h_py_pixels->GetYaxis()->CenterTitle();

    TH2D* h_8   = new TH2D("h_8","Y vs Time",(UT_max-UT_min)/1000.0,UT_min/1000.0,UT_max/1000.0,N_PIXELS,0,N_PIXELS);

    //------------------------------------------------------------------------------//
    Double_t event_time = -999.999;

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
        if (_Timems > 0)
        {
            event_time  = _Timems;
        }
        else
        {
            continue;
        }

        if(event_time/1000.0 < start_time) continue;
        if(event_time/1000.0 > stop_time) break;


        TH2D* h_1 = new TH2D("h_1","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
        TH2D* h_2 = new TH2D("h_2","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                if(_COUNTS[xi][yi] > 0 && _AcquisType == 0)
                {
                    h_1->Fill(xi,yi,_COUNTS[xi][yi]);
                    h_2->Fill(xi,yi,_COUNTS[xi][yi]);
                }
            }
        }

        median_filter(h_1,h_2,2.0);

        Double_t integral, integral_err;
        integral= h_2->IntegralAndError(1,h_2->GetNbinsX(),1,h_2->GetNbinsY(),integral_err);

        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                h_py_pixels->Fill(yi,h_2->GetBinContent(xi+1,yi+1));
                h_8->Fill(event_time/1000.0,yi,h_2->GetBinContent(xi+1,yi+1));
            }
        }

        for(Int_t yi = 1; yi <= h_8->GetNbinsY(); yi++)
        {
            h_8->SetBinContent(h_8->GetXaxis()->FindBin(event_time/1000.0),yi,h_8->GetBinContent(h_8->GetXaxis()->FindBin(event_time/1000.0),yi)/integral);
            h_8->SetBinError(h_8->GetXaxis()->FindBin(event_time/1000.0),yi,
                             TMath::Sqrt(TMath::Power(h_8->GetBinError(h_8->GetXaxis()->FindBin(event_time/1000.0),yi)/integral,2) +
                                         TMath::Power(h_8->GetBinContent(h_8->GetXaxis()->FindBin(event_time/1000.0),yi)*integral_err/(integral*integral),2)));
        }

        h_1->Delete();
        h_2->Delete();
    }
    //---------------------------------------------------------------------------------------//

    TH1D* h_py_mm = new TH1D("h_py_mm","Projection on Horizontal Axis (mm)",N_PIXELS,0,N_PIXELS*0.055);
    h_py_mm->GetXaxis()->SetTitle("Horizontal Axis [mm]");
    h_py_mm->GetXaxis()->SetTimeOffset(1.2);
    h_py_mm->GetXaxis()->CenterTitle();
    h_py_mm->GetYaxis()->SetTitle("Normalized Counts");
    h_py_mm->GetYaxis()->SetTimeOffset(1.2);
    h_py_mm->GetYaxis()->CenterTitle();

    Double_t integral, integral_err;
    integral = h_py_pixels->IntegralAndError(1,h_py_pixels->GetNbinsX(),integral_err);
    scale1Dhisto(h_py_pixels,integral,integral_err);

    for(Int_t j = 1; j <= h_py_pixels->GetNbinsX(); j++)
    {
        h_py_mm->SetBinContent(j,h_py_pixels->GetBinContent(j));
        h_py_mm->SetBinError(j,h_py_pixels->GetBinError(j));
    }

    TCanvas* c_2 = new TCanvas("c_2","c_2",1800,900);
    c_2->cd();
    gPad->SetGrid();
    h_8->Draw("colz");

    TCanvas* c_0 = new TCanvas("c_0","c_0",900,900);
    c_0->cd();
    gPad->SetGrid();
    h_py_pixels->GetXaxis()->SetRange(1,N_PIXELS);
    h_py_pixels->Draw();

    TCanvas* c_1 = new TCanvas("c_1","c_1",900,900);
    c_1->cd();
    gPad->SetGrid();
    h_py_mm->GetXaxis()->SetRange(1,N_PIXELS);
    h_py_mm->Draw("hist");

    Double_t par[5];
    TF1* f_1 = new TF1("f_1","expo",1.0,7.0);
    TF1* f_2 = new TF1("f_2","gaus",9.5,13);
    h_py_mm->Fit(f_1,"R0Q+");
    h_py_mm->Fit(f_2,"R0Q+");
    f_1->GetParameters(&par[0]);
    f_2->GetParameters(&par[2]);
    TF1* f_3 = new TF1("f_3","expo(0) + gaus(2)",1.0,13);
    f_3->SetParameters(par);
    h_py_mm->Fit(f_3,"R+");
    f_3->SetLineColor(kRed);
    f_3->SetLineWidth(3);
    f_3->Draw("same");
    cout<<"--> chi2/NDF = "<<f_3->GetChisquare()<<"/"<<f_3->GetNDF()<<endl;
}

void function_10()
{
    //---------------------------------------------------------------------------------------//

    TChain* fChain_1;
    TChain* fChain_3;
    Long64_t nEntries;

    Long64_t _COUNTS_1[N_PIXELS*2][N_PIXELS*2];
    Long64_t _COUNTS_3[N_PIXELS*2][N_PIXELS*2];

    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L1_RUN_8.root";
    TString outFileName     = "output_L1_calibration.root";

    fChain_1 = new TChain("Tree_1");
    fChain_3 = new TChain("Tree_3");

    fChain_1->Add(inFileName1.Data());
    fChain_3->Add(inFileName1.Data());

    fChain_1->SetBranchAddress("COUNTS",      _COUNTS_1);
    fChain_3->SetBranchAddress("COUNTS",      _COUNTS_3);

    if(fChain_1->GetEntries() != fChain_3->GetEntries())
    {
        cout<<endl<<endl<<"### ERROR::Different number of entries!!! ###"<<endl<<endl<<endl;
        assert(0);
    }

    nEntries = fChain_1->GetEntries();
    cout<<"--> InputFileName: "<<inFileName1<<endl;
    cout<<"--> Number of Frames: "<<nEntries<<endl;
    //------------------------------------------------------------------------------//
    //------------------------------- HISTOGRAMS -----------------------------------//
    //------------------------------------------------------------------------------//

    TH1D* h_1 = new TH1D("h_1","RP0 Int.",1e3,0,1e6);
    TH1D* h_2 = new TH1D("h_2","RP1 Int.",1e3,0,1e8);
    TH2D* h_3 = new TH2D("h_3","RP1 Int. vs RP0 Int.",1e3,0,1e6,1e3,0,1e8);

    //---------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------- MAIN LOOP ------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------//
    for(Long64_t i = 0; i < nEntries; i++)
    {
        TH2D* h_1_1 = new TH2D("h_1_1","Y vs X vs C (in pixels) [RP1I]",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
        TH2D* h_2_1 = new TH2D("h_2_1","Y vs X vs C (in pixels) [RP1I]",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

        TH2D* h_1_3 = new TH2D("h_1_3","Y vs X vs C (in pixels) [RP0I]",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
        TH2D* h_2_3 = new TH2D("h_2_3","Y vs X vs C (in pixels) [RP0I]",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);

        if(i%1 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntries);
            fflush(stdout);
        }

        fChain_1->GetEntry(i);
        fChain_3->GetEntry(i);

        for(Int_t xi = 1; xi < N_PIXELS - 1; xi++)
        {
            for(Int_t yi = 1; yi < N_PIXELS - 1; yi++)
            {
                if(_COUNTS_1[xi][yi] > 0)
                {
                    h_1_1->Fill(xi,yi,_COUNTS_1[xi][yi]);
                    h_2_1->Fill(xi,yi,_COUNTS_1[xi][yi]);
                }
                if(_COUNTS_3[xi][yi] > 0)
                {
                    h_1_3->Fill(xi,yi,_COUNTS_3[xi][yi]);
                    h_2_3->Fill(xi,yi,_COUNTS_3[xi][yi]);
                }
            }
        }

        median_filter(h_1_1,h_2_1,2.0);
        median_filter(h_1_3,h_2_3,2.0);

        // To check the algorythm

//        TCanvas* c_3_1 = new TCanvas("c_3_1","c_3_1",1800,900);
//        c_3_1->Divide(2,1);

//        c_3_1->cd(1);
//        gPad->SetGrid();
//        h_1_1->GetXaxis()->SetRange(1,256);
//        h_1_1->GetYaxis()->SetRange(1,256);
//        h_1_1->Draw("colz");

//        c_3_1->cd(2);
//        gPad->SetGrid();
//        h_2_1->GetXaxis()->SetRange(1,256);
//        h_2_1->GetYaxis()->SetRange(1,256);
//        h_2_1->Draw("colz");

//        TCanvas* c_3_2 = new TCanvas("c_3_2","c_3_2",1800,900);
//        c_3_2->Divide(2,1);

//        c_3_2->cd(1);
//        gPad->SetGrid();
//        h_1_3->GetXaxis()->SetRange(1,256);
//        h_1_3->GetYaxis()->SetRange(1,256);
//        h_1_3->Draw("colz");

//        c_3_2->cd(2);
//        gPad->SetGrid();
//        h_2_3->GetXaxis()->SetRange(1,256);
//        h_2_3->GetYaxis()->SetRange(1,256);
//        h_2_3->Draw("colz");

        Double_t integral_1, integral_3;
        integral_1= h_2_1->Integral(1,h_2_1->GetNbinsX(),1,h_2_1->GetNbinsY());
        integral_3= h_2_3->Integral(1,h_2_3->GetNbinsX(),1,h_2_3->GetNbinsY());

        h_1->Fill(integral_3);
        h_2->Fill(integral_1);
        h_3->Fill(integral_3,integral_1);

        h_1_1->Delete();
        h_2_1->Delete();

        h_1_3->Delete();
        h_2_3->Delete();
    }
    cout<<endl;

    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File '"<<_file->GetName()<<"' opened successfully!"<<endl;
    h_1->Write();
    h_2->Write();
    h_3->Write();
    _file->Close();
    //---------------------------------------------------------------------------------------//
}

void function_11()
{
//    gStyle->SetOptStat(0);
    TFile *_file_chipid_chip1 = TFile::Open("output_L1_chipid_1.root");
    TGraph *_cry3_angularscan = new TGraph("crystal3_angularscan.dat","%lg %lg");
    _cry3_angularscan->SetName("_cry3_angularscan");
    _cry3_angularscan->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    _cry3_angularscan->GetXaxis()->CenterTitle();
    _cry3_angularscan->GetXaxis()->SetTitle("Time");

    Double_t time_start = 1537238370;   // Tuesday, 18 September 2018, 04:39:30
    Double_t time_stop  = 1537238880;   // Tuesday, 18 September 2018, 04:48:00

//    to fit with pol1
    TF1 *_cry3_angularscan_spline = new TF1("_cry3_angularscan_spline","pol1",time_start,time_stop);
    _cry3_angularscan->Fit(_cry3_angularscan_spline,"R+");

    TH2D* h_chip1_1 = (TH2D*)_file_chipid_chip1->Get("h_rp1_6");
    TH2D* h_chip1_2 = new TH2D("h_chip1_HvsA","h_chip1_HvsA",600,-2000,-1400,N_PIXELS,0,N_PIXELS*0.055);

    vector<Double_t> ratio;
    vector<Double_t> ratio_err;
    vector<Double_t> angle;
    vector<Double_t> angle_err;

    Double_t fit_bkg_lim_min = 256*0.055 - 12.7;
    Double_t fit_bkg_lim_max = 256*0.055 - 13.2;
    Double_t fit_singlech_lim_min = 256*0.055 - 2.0;
    Double_t fit_singlech_lim_max = 256*0.055 - 4.5;
    Double_t par_bkg[2];

    for(Int_t i = 1; i <= h_chip1_1->GetNbinsX(); i++)
    {
        Double_t _time      = h_chip1_1->GetXaxis()->GetBinCenter(i);
        Double_t _angle     = _cry3_angularscan_spline->Eval(_time);

        for(Int_t j = 1; j <= h_chip1_1->GetNbinsY(); j++)
        {
            Double_t _position  = h_chip1_1->GetYaxis()->GetBinCenter(j);
            Double_t _value     = h_chip1_1->GetBinContent(i,j);
            Double_t _value_err = h_chip1_1->GetBinError(i,j);

            if(_time > time_start && _time < time_stop)
            {
                h_chip1_2->SetBinContent(h_chip1_2->GetXaxis()->FindBin(_angle),h_chip1_2->GetYaxis()->FindBin(_position),_value);
                h_chip1_2->SetBinError(h_chip1_2->GetXaxis()->FindBin(_angle),h_chip1_2->GetYaxis()->FindBin(_position),_value_err);
            }
        }

        TH1D* h_1 = (TH1D*)h_chip1_1->ProjectionY("h_1",i,i);
        TH1D* h_1_s = (TH1D*)h_chip1_1->ProjectionY("h_1_s",i,i);
        TH1D* h_1_s_s = (TH1D*)h_chip1_1->ProjectionY("h_1_s_s",i,i);

        if(h_1->Integral(1,h_1->GetNbinsX()) != 1) continue;
        //--------------------------------------------------------------//
        // To count the ratio
        //--------------------------------------------------------------//

        TF1* fit_funct_ch_bkg = new TF1("fit_funct_ch_bkg","pol0",fit_bkg_lim_min,fit_bkg_lim_max);
        fit_funct_ch_bkg->SetParameter(0, 2e-3);
        h_1->Fit(fit_funct_ch_bkg,"R0Q+");

        fit_funct_ch_bkg->GetParameters(&par_bkg[0]);
        TF1* bkg_func = new TF1("bkg_func","pol0",h_1->GetBinCenter(1),h_1->GetBinCenter(h_1->GetNbinsX()));
        bkg_func->SetParameters(par_bkg);

        for(Int_t j = 1; j <= h_1_s_s->GetNbinsX(); j++)
        {
            Double_t val = h_1_s->GetBinContent(j) - bkg_func->Eval(h_1_s->GetBinCenter(j));

            if(val > 0)
                h_1_s_s->SetBinContent(j,val);
            else
                h_1_s_s->SetBinContent(j,0);
        }

        TF1* fit_funct_1 = new TF1("fit_funct_1","gaus",fit_singlech_lim_min,fit_singlech_lim_max);
        h_1_s_s->Fit(fit_funct_1,"R0Q+");

        Double_t min_count_am = fit_funct_1->GetParameter(1)-3.0*fit_funct_1->GetParameter(2);
        Double_t max_count_am = 256*0.055;

        Double_t I1, I1_err;
        Double_t I2, I2_err;
        I1 = h_1_s_s->IntegralAndError(1,h_1_s_s->GetNbinsX(),I1_err);
        I2 = h_1_s_s->IntegralAndError(h_1_s_s->GetXaxis()->FindBin(min_count_am),h_1_s_s->GetXaxis()->FindBin(max_count_am),I2_err);

        ratio.push_back(1.0 - I2/I1);
        ratio_err.push_back(TMath::Sqrt(TMath::Power(I2_err/I1,2) + TMath::Power(I2*I1_err/(I1*I1),2)));
        angle.push_back(_angle);
        angle_err.push_back(1.0);

        //--------------------------------------------------------------//
        // To plot
        //--------------------------------------------------------------//

//        TCanvas* c_2_temp = new TCanvas("c_2_temp","c_2_temp",1000,500);
//        c_2_temp->cd();
//        h_1->SetMinimum(0);
//        h_1->Draw("hist");
//        fit_funct_1->SetLineColor(kBlue);
//        fit_funct_1->SetLineWidth(4);
//        fit_funct_1->Draw("same");

//        TString name = "#chi^{2}/NDF = [";
//        name += (Int_t)fit_funct_ch_bkg->GetChisquare();
//        name += "/";
//        name += fit_funct_ch_bkg->GetNDF();
//        name += "] = ";
//        name += (Int_t)fit_funct_ch_bkg->GetChisquare()/fit_funct_ch_bkg->GetNDF();
//        h_1->SetTitle(name.Data());
//        fit_funct_ch_bkg->SetLineWidth(6);
//        fit_funct_ch_bkg->Draw("same");

//        bkg_func->SetLineColor(kGreen+2);
//        bkg_func->Draw("same");

//        TCanvas* c_3_temp = new TCanvas("c_3_temp","c_3_temp",1000,500);
//        c_3_temp->cd();
//        h_1_s->SetMinimum(0);
//        h_1_s->SetLineColor(kBlue);
//        h_1_s->SetLineWidth(2);
//        h_1_s->Draw("e1");

//        h_1_s_s->SetLineColor(kRed);
//        h_1_s_s->SetLineWidth(2);
//        h_1_s_s->Draw("same & e1");

//        TCanvas* c_4_temp = new TCanvas("c_4_temp","c_4_temp",1000,500);
//        c_4_temp->cd();
//        gPad->SetGrid();
//        h_1_s_s->SetLineWidth(2);
//        h_1_s_s->SetTitle("");
//        h_1_s_s->GetXaxis()-> CenterTitle();
//        h_1_s_s->GetYaxis()-> CenterTitle();
//        h_1_s_s->Draw("e1");

//        TLine *l_min_ch=new TLine(min_count_am,0,min_count_am,h_1_s_s->GetMaximum());
//        l_min_ch->SetLineColor(kMagenta);
//        l_min_ch->SetLineWidth(3);
//        l_min_ch->SetLineStyle(7);
//        l_min_ch->Draw();
//        TLine *l_mean_ch=new TLine(fit_funct_1->GetParameter(1),0,fit_funct_1->GetParameter(1),h_1_s_s->GetMaximum());
//        l_mean_ch->SetLineColor(kMagenta-2);
//        l_mean_ch->SetLineWidth(3);
//        l_mean_ch->SetLineStyle(7);
//        l_mean_ch->Draw();
//        TLine *l_max_ch=new TLine(max_count_am,0,max_count_am,h_1_s_s->GetMaximum());
//        l_max_ch->SetLineColor(kMagenta);
//        l_max_ch->SetLineWidth(3);
//        l_max_ch->SetLineStyle(7);
//        l_max_ch->Draw();

//        break;
        //--------------------------------------------------------------//
    }
    Int_t dim = ratio.size();

    TCanvas* c_0 = new TCanvas("c_0","c_0",1800,900);
    c_0->cd();
    gPad->SetGrid();
    h_chip1_1->Draw("colz");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    TGraphErrors *gr_1 = new TGraphErrors(dim, &angle[0], &ratio[0], &angle_err[0], &ratio_err[0]);
    gr_1->SetMarkerStyle(21);
    gr_1->SetMarkerColor(kBlack);
    gr_1->SetLineWidth(2);
    gr_1->SetLineColor(kBlue);
    gr_1->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    gr_1->GetYaxis()->SetTitle("Ratio (1 - AM/TOTAL)");
    gr_1->GetYaxis()->CenterTitle();
    gr_1->GetXaxis()->CenterTitle();
    gr_1->GetXaxis()->SetLimits(-1900,-1600);
    gr_1->Draw("AP");

    TCanvas* c_2 = new TCanvas("c_2","c_2",1000,1000);
    c_2->cd();
    gPad->SetGrid();
    _cry3_angularscan->SetMarkerStyle(21);
    _cry3_angularscan->Draw("AP");
    _cry3_angularscan_spline->SetLineColor(kRed);
    _cry3_angularscan_spline->SetLineWidth(2);
    _cry3_angularscan_spline->Draw("same & L");

    TCanvas* c_3 = new TCanvas("c_3","c_3",1800,900);
    c_3->cd();
    gPad->SetGrid();
    h_chip1_2->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(-1900),h_chip1_2->GetXaxis()->FindBin(-1600));
    h_chip1_2->SetTitle("");
    h_chip1_2->GetXaxis()->CenterTitle();
    h_chip1_2->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    h_chip1_2->GetYaxis()->CenterTitle();
    h_chip1_2->GetYaxis()->SetTitle("Projection on Horizontal Axis [mm]");
    h_chip1_2->SetMinimum(0.004);
    h_chip1_2->Draw("colz");

    //-----------//
    // METHOD I
    //-----------//

    TCanvas* c_4= new TCanvas("c_4","c_4",1800,900);
    c_4->cd();
    gPad->SetGrid();
    TGraphErrors *gr_2 = new TGraphErrors(dim, &angle[0], &ratio[0], &angle_err[0], &ratio_err[0]);
    gr_2->SetMarkerStyle(21);
    gr_2->SetMarkerColor(kBlack);
    gr_2->SetLineWidth(2);
    gr_2->SetLineColor(kBlue);
    gr_2->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    gr_2->GetYaxis()->SetTitle("Ratio (1 - AM/TOTAL)");
    gr_2->GetYaxis()->CenterTitle();
    gr_2->GetXaxis()->CenterTitle();
    gr_2->GetXaxis()->SetLimits(-1900,-1600);
    gr_2->Draw("AP");

    TF1* fit_landau = new TF1("fit_landau","landau",-1900,-1600);
    TFitResultPtr fit_results_landau = gr_2->Fit(fit_landau,"SR+");

    double angle_fix[1] = {-1820.0};
    double angle_fix_err[1];

    fit_results_landau->GetConfidenceIntervals(1, 1, 1, angle_fix, angle_fix_err, 0.683, false);

    /*
     * void ROOT::Fit::FitResult::GetConfidenceIntervals	(	unsigned int 	n,   unsigned int 	stride1,    unsigned int 	stride2,    const double * 	x,    double * 	ci,    double 	cl = 0.95,    bool 	norm = true     )		const
     *get confidence intervals for an array of n points x.
     * stride1 indicates the stride in the coordinate space while stride2 the stride in dimension space.
     * For 1-dim points : stride1=1, stride2=1
     * for multi-dim points arranged as (x0,x1,...,xN,y0,....yN) stride1=1 stride2=n
     * for multi-dim points arraged as (x0,y0,..,x1,y1,...,xN,yN,..) stride1=ndim, stride2=1
     * the confidence interval are returned in the array ci
     * cl is the desired confidedence interval
     * value norm is a flag to control if the intervals need to be normalized to the chi2/ndf value
     * By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
    */

    Double_t calibr = 0.195/fit_landau->Eval(angle_fix[0]);
    Double_t calibr_err = TMath::Sqrt(TMath::Power(3.1e-5/fit_landau->Eval(angle_fix[0]),2) +
            TMath::Power(0.195*angle_fix_err[0]/(fit_landau->Eval(angle_fix[0])*fit_landau->Eval(angle_fix[0])),2));

    TH1D* h_dist_landau = new TH1D("h_dist_landau","CH Efficiency",13,-1900,-1614);

    for(Int_t i = 1; i <= h_dist_landau->GetNbinsX(); i++)
    {
        Double_t x = h_dist_landau->GetBinCenter(i);

        angle_fix[0]        = x;
        angle_fix_err[0]    = 0;
        fit_results_landau->GetConfidenceIntervals(1, 1, 1, angle_fix, angle_fix_err, 0.683, false);

        h_dist_landau->SetBinContent(i,calibr*fit_landau->Eval(x));
        h_dist_landau->SetBinError(i,TMath::Sqrt(TMath::Power(calibr_err*fit_landau->Eval(x),2) +
                                                 TMath::Power(calibr*angle_fix_err[0],2)));
    }

    TCanvas* c_5= new TCanvas("c_5","c_5",1800,900);
    c_5->cd();
    gPad->SetGrid();
    h_dist_landau->SetTitle("1 bin = 2x #Theta_{c}");
    h_dist_landau->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    h_dist_landau->GetYaxis()->SetTitle("Ratio (1 - AM/TOTAL)");
    h_dist_landau->GetYaxis()->CenterTitle();
    h_dist_landau->GetXaxis()->CenterTitle();
    h_dist_landau->SetMarkerStyle(21);
    h_dist_landau->SetLineWidth(2);
    h_dist_landau->Draw("e1");

    Double_t ch_eff, ch_eff_err;
    ch_eff = h_dist_landau->IntegralAndError(h_dist_landau->GetXaxis()->FindBin(fit_landau->GetParameter(1)-3.0*fit_landau->GetParameter(2)),
                                      h_dist_landau->GetXaxis()->FindBin(fit_landau->GetParameter(1)+3.0*fit_landau->GetParameter(2)),
                                      ch_eff_err);
    cout<<"--> [METHOD I] CH EFF. = "<<ch_eff<<" +/- "<<ch_eff_err<<endl;

    //-----------//
    // METHOD II
    //-----------//

    TCanvas* c_6= new TCanvas("c_6","c_6",1800,900);
    c_6->cd();
    gPad->SetGrid();
    TGraphErrors *gr_3 = new TGraphErrors(dim, &angle[0], &ratio[0], &angle_err[0], &ratio_err[0]);
    gr_3->SetMarkerStyle(21);
    gr_3->SetMarkerColor(kBlack);
    gr_3->SetLineWidth(2);
    gr_3->SetLineColor(kBlue);
    gr_3->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    gr_3->GetYaxis()->SetTitle("Ratio (1 - AM/TOTAL)");
    gr_3->GetYaxis()->CenterTitle();
    gr_3->GetXaxis()->CenterTitle();
    gr_3->GetXaxis()->SetLimits(-1900,-1600);
    gr_3->Draw("AP");

    TF1* fit_gaus = new TF1("fit_gaus","gaus",fit_landau->GetParameter(1)-2.0*fit_landau->GetParameter(2),fit_landau->GetParameter(1)+2.0*fit_landau->GetParameter(2));
    TFitResultPtr fit_results_gaus = gr_3->Fit(fit_gaus,"SR+");

    TF1* gaus_func = new TF1("gaus_func","gaus",-1900,-1600);
    gaus_func->SetParameter(0,fit_gaus->GetParameter(0));
    gaus_func->SetParameter(1,fit_gaus->GetParameter(1));
    gaus_func->SetParameter(2,fit_gaus->GetParameter(2));
    gaus_func->Draw("same");

    angle_fix[0]        = -1820.0;
    angle_fix_err[0]    = 0;

    fit_results_gaus->GetConfidenceIntervals(1, 1, 1, angle_fix, angle_fix_err, 0.683, false);

    calibr = 0.195/gaus_func->Eval(angle_fix[0]);
    calibr_err = TMath::Sqrt(TMath::Power(3.1e-5/gaus_func->Eval(angle_fix[0]),2) +
            TMath::Power(0.195*angle_fix_err[0]/(gaus_func->Eval(angle_fix[0])*gaus_func->Eval(angle_fix[0])),2));

    TH1D* h_dist_gaus = new TH1D("h_dist_gaus","CH Efficiency",13,-1900,-1614);

    for(Int_t i = 1; i <= h_dist_gaus->GetNbinsX(); i++)
    {
        Double_t x = h_dist_gaus->GetBinCenter(i);

        angle_fix[0]        = x;
        angle_fix_err[0]    = 0;
        fit_results_gaus->GetConfidenceIntervals(1, 1, 1, angle_fix, angle_fix_err, 0.683, false);

        h_dist_gaus->SetBinContent(i,calibr*gaus_func->Eval(x));
        h_dist_gaus->SetBinError(i,TMath::Sqrt(TMath::Power(calibr_err*gaus_func->Eval(x),2) +
                                               TMath::Power(calibr*angle_fix_err[0],2)));
    }

    TCanvas* c_7= new TCanvas("c_7","c_7",1800,900);
    c_7->cd();
    gPad->SetGrid();
    h_dist_gaus->SetTitle("1 bin = 2x #Theta_{c}");
    h_dist_gaus->GetXaxis()->SetTitle("CRYSTAL3 LVDT Angle [#murad]");
    h_dist_gaus->GetYaxis()->SetTitle("Ratio (1 - AM/TOTAL)");
    h_dist_gaus->GetYaxis()->CenterTitle();
    h_dist_gaus->GetXaxis()->CenterTitle();
    h_dist_gaus->SetMarkerStyle(21);
    h_dist_gaus->SetLineWidth(2);
    h_dist_gaus->Draw("e1");

    ch_eff = h_dist_gaus->IntegralAndError(h_dist_gaus->GetXaxis()->FindBin(gaus_func->GetParameter(1)-3.0*gaus_func->GetParameter(2)),
                                           h_dist_gaus->GetXaxis()->FindBin(gaus_func->GetParameter(1)+3.0*gaus_func->GetParameter(2)),
                                           ch_eff_err);
    cout<<"--> [METHOD II] CH EFF. = "<<ch_eff<<" +/- "<<ch_eff_err<<endl;
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

void rewriteHisto(TH2D* h_in, TH2D* h_out)// to rewrite a histogram
{
    for(Int_t i = 1; i <= h_in->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h_in->GetNbinsY(); j++)
        {
            h_out->SetBinContent(i,j,h_in->GetBinContent(i,j));
            h_out->SetBinError(i,j,h_in->GetBinError(i,j));
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
