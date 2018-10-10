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

void repairnoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels);
void repairnoizypixelsYT(TH2D* histo, Float_t factor_bad_pixels);
void scale2Dhisto(TH2D* histo, Double_t integral, Double_t integral_err);
void scale2Dhisto_YT(TH2D* histo);
void rotateNscale(TH2D* h_in, TH2D* h_out);
void rewriteHisto(TH2D* h_in, TH2D* h_out);
double getAverage2D(TH2D* histo);
void median_filter(TH2D* h_in, TH2D* h_out, Double_t factor_bad_pixels);
double findMedian(double a[], int n);

int sps_md_analysis_17_09_2018_double_channeling()
{
    cout<<"--> function_1() -- [L1,L2,L3] to plot rp_i 2d beam image normalized by rp_i counts"<<endl;
    cout<<"--> function_2() -- [L1] to plot Timepix data vs CRY3 angle"<<endl;
    cout<<"--> function_3() -- [L3] to plot Timepix data vs CRY2 angle"<<endl;
    cout<<"--> function_4() -- to plot Timepix projection"<<endl;
    cout<<"--> function_5() -- [L2] double CH stability"<<endl;
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
    Long64_t _COUNTS[N_PIXELS][N_PIXELS];

    TString inFileName1     = "/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L3_RUN_8.root";
    TString outFileName     = "output_L3_chipid_1.root";

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
    TH2D* h_4 = new TH2D("h_4","RP1 Internal (normalized)",N_PIXELS,0,N_PIXELS*0.055,N_PIXELS,0,N_PIXELS*0.055);
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
        TH2D* h_1   = new TH2D("h_1","Y vs X vs C (in pixels)",N_PIXELS,0,N_PIXELS,N_PIXELS,0,N_PIXELS);
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

//        repairnoizypixelsXY(h_2,2.0);
        median_filter(h_1,h_2,2.0);

        /*TCanvas* c_3 = new TCanvas("c_3","c_3",1800,900);
        c_3->Divide(2,1);

        c_3->cd(1);
        gPad->SetGrid();
        h_1->GetXaxis()->SetRange(1,256);
        h_1->GetYaxis()->SetRange(1,256);
        h_1->Draw("colz");

        c_3->cd(2);
        gPad->SetGrid();
        h_2->GetXaxis()->SetRange(1,256);
        h_2->GetYaxis()->SetRange(1,256);
        h_2->Draw("colz");*/

        Double_t integral, integral_err;
        integral= h_2->IntegralAndError(1,h_2->GetNbinsX(),1,h_2->GetNbinsY(),integral_err);

        for(Int_t xi = 0; xi < N_PIXELS; xi++)
        {
            for(Int_t yi = 0; yi < N_PIXELS; yi++)
            {
                h_3->Fill(xi,yi,h_2->GetBinContent(xi+1,yi+1));
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

        Int_t dtime = ((event_time - UT_min)/1000 + 0.5);
        TString name = "./frames_l3_chip1/rp1int_";
        name += dtime;
        TString title = "dTime: ";
        title += dtime;
        title += " [sec]";
        h_2->SetTitle(title.Data());
        h_2->SetMaximum(2000);
        TCanvas* cc = new TCanvas(name.Data(),"cc",1000,1000);
        cc->cd();
        gPad->SetGrid();
        h_2->GetXaxis()->SetRange(1,256);
        h_2->GetYaxis()->SetRange(1,256);
        h_2->Draw("colz");
        name += ".png";
        cc->SaveAs(name.Data());
        h_1->Delete();
        h_2->Delete();
    }
    //---------------------------------------------------------------------------------------//

    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File '"<<_file->GetName()<<"' opened successfully!"<<endl;

    TH2D* h_rp1_6 = new TH2D("h_rp1_6","RP1 Internal (normalized)",h_8->GetNbinsX(),h_8->GetXaxis()->GetBinCenter(1),h_8->GetXaxis()->GetBinCenter(h_8->GetNbinsX())+h_8->GetXaxis()->GetBinWidth(1),h_8->GetNbinsY(),0,h_8->GetNbinsY()*0.055);
    h_rp1_6->GetXaxis()->SetTitle("Time");
    h_rp1_6->GetXaxis()->SetTimeOffset(1.2);
    h_rp1_6->GetXaxis()->CenterTitle();
    h_rp1_6->GetYaxis()->SetTitle("Projection on Horizontal Axis [mm]");
    h_rp1_6->GetYaxis()->SetTimeOffset(1.2);
    h_rp1_6->GetYaxis()->CenterTitle();

    Double_t integral, integral_err;
    integral = h_3->IntegralAndError(1,h_3->GetNbinsX(),1,h_3->GetNbinsY(),integral_err);
    scale2Dhisto(h_3,integral,integral_err);

    rewriteHisto(h_8,h_rp1_6);

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

    TCanvas* c_2 = new TCanvas("c_2","c_2",1800,900);
    c_2->cd();
    gPad->SetGrid();
    h_rp1_6->GetYaxis()->SetRange(1,256);
    h_rp1_6->GetXaxis()->SetTimeDisplay(1);
    h_rp1_6->Draw("colz");

    c_0->Write();
    c_00->Write();
    c_1->Write();
    c_2->Write();
    h_3->Write();
    h_4->Write();
    h_8->Write();
    h_rp1_6->Write();

    _file->Close();
    //---------------------------------------------------------------------------------------//
}

void function_2()
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

//    to fit with splines
//    TSpline3 *_cry3_angularscan_spline = new TSpline3("_cry3_angularscan_spline",_cry3_angularscan);

//    to fit with pol1
    TF1 *_cry3_angularscan_spline = new TF1("_cry3_angularscan_spline","pol1",time_start,time_stop);
    _cry3_angularscan->Fit(_cry3_angularscan_spline,"R+");

    TH2D* h_chip1_1 = (TH2D*)_file_chipid_chip1->Get("h_rp1_6");
    TH2D* h_chip1_2 = new TH2D("h_chip1_HvsA","h_chip1_HvsA",600,1400,2000,N_PIXELS/2,0,N_PIXELS*0.055/2);

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

    TString outFileName     = "output_function_2.root";
    TFile *_file = new TFile(outFileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> File opened successfully!"<<endl;

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
//    h_chip1_2->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(1760),h_chip1_2->GetXaxis()->FindBin(1880));
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

    TCanvas* c_3 = new TCanvas("c_3","c_3",900,900);
    c_3->cd();
    gPad->SetGrid();
    h_chip1_3_projection->Draw("hist");
    h_chip1_3_projection->GetXaxis()->SetRange(h_chip1_2->GetXaxis()->FindBin(1760),h_chip1_2->GetXaxis()->FindBin(1880));

    c_0->Write();
    _cry3_angularscan->Write();
    _file->Close();
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
    TH2D* h_chip1_2 = new TH2D("h_chip1_HvsA","h_chip1_HvsA",7000,-6700,-6000,N_PIXELS/2,0,N_PIXELS*0.055/2);

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
    Double_t fit_doublech_lim_min = 10.4;
    Double_t fit_doublech_lim_max = 12.2;
    Double_t fit_singlech_lim_min =  2.0;
    Double_t fit_singlech_lim_max =  4.5;
    Double_t par_bkg[1];

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
    par_bkg[0] = fit_funct_ch_bkg->GetParameter(0);

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
    bkg_func->SetParameter(0,par_bkg[0]);
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
    for(Int_t i = 1; i <= h_in->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h_in->GetNbinsY(); j++)
        {
            h_out->SetBinContent(h_in->GetNbinsX()/2-j+1,i,h_in->GetBinContent(i,j));
            h_out->SetBinError(h_in->GetNbinsX()/2-j+1,i,h_in->GetBinError(i,j));
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
}

double findMedian(double a[], int n)// Function for calculating median
{
    sort(a, a+n);
    return a[n/2];
}
