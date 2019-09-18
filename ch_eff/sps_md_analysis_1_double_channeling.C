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

const Int_t N_PIXELS                = 512;
const Double_t PIXEL_SIZE           = 0.055;

int sps_md_analysis_1_double_channeling()
{

    //----------------------------------------------//
    //---------------- 17.10.2017 ------------------//
    //----------------------------------------------//
    TString file_name_am = "MD_2017_10_17_CRY2CHCRY4AMP1_HISTO_RP1I_RUN_1.root";
    TString file_name_bk = "MD_2017_10_17_CRY2CH_HISTO_RP1I_RUN_1.root";
    TString file_name_ch = "MD_2017_10_17_CRY2CHCRY4CHP1_HISTO_RP1I_RUN_1.root";
    TString fit_ch_func_name = "gaus";
    Int_t min_x_ch_fit = 50;
    Int_t max_x_ch_fit = 150;
    Int_t min_x_dch_ln = 120;
    Int_t max_x_dch_ln = 145;
    Int_t num_sigma_x_fit_ch = 3;
    Float_t factor_bad_pixels = 10000.0;
    TString crystal_info = "CRYSTAL2 & CRYSTAL4 [17.10.2017]";
    //----------------------------------------------//
    TFile *_file_am = TFile::Open(file_name_am.Data());
    TFile *_file_bk = TFile::Open(file_name_bk.Data());
    TFile *_file_ch = TFile::Open(file_name_ch.Data());

    TH2D* h_evn_am = (TH2D*)_file_am->Get("h_10");
    TH2D* h_evn_bk = (TH2D*)_file_bk->Get("h_10");
    TH2D* h_evn_ch = (TH2D*)_file_ch->Get("h_10");

    Double_t time_am = h_evn_am->GetEntries()*0.000246; // interated time
    Double_t time_bk = h_evn_bk->GetEntries()*0.000246; // interated time
    Double_t time_ch = h_evn_ch->GetEntries()*0.000246; // interated time

    cout<<"--> time_am = "<<time_am<<" [sec]"<<endl;
    cout<<"--> time_bk = "<<time_bk<<" [sec]"<<endl;
    cout<<"--> time_ch = "<<time_ch<<" [sec]"<<endl;

    TH2D* h_ch_am = (TH2D*)_file_am->Get("h_1");
    TH2D* h_ch_bk = (TH2D*)_file_bk->Get("h_1");
    TH2D* h_ch_ch = (TH2D*)_file_ch->Get("h_1");

    TH2D* h_ch_clone_am = new TH2D("h_ch_clone_am","h_ch_clone_am",N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);
    TH2D* h_ch_clone_bk = new TH2D("h_ch_clone_bk","h_ch_clone_bk",N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);
    TH2D* h_ch_clone_ch = new TH2D("h_ch_clone_ch","h_ch_clone_ch",N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);

    for(Int_t xi = 1; xi <= N_PIXELS/2; xi++)
    {
        for(Int_t yi = 1; yi <= N_PIXELS/2; yi++)
        {
            if(h_ch_am->GetBinContent(xi,yi) - factor_bad_pixels*h_ch_am->GetBinError(xi,yi)
                    < factor_bad_pixels*h_ch_am->GetBinError(xi,yi+1) + h_ch_am->GetBinContent(xi,yi+1))
            {
                h_ch_clone_am->SetBinContent(xi,yi,h_ch_am->GetBinContent(xi,yi)/time_am);
                h_ch_clone_am->SetBinError(xi,yi,h_ch_am->GetBinError(xi,yi)/time_am);
            }

            if(h_ch_bk->GetBinContent(xi,yi) - factor_bad_pixels*h_ch_bk->GetBinError(xi,yi)
                    < factor_bad_pixels*h_ch_bk->GetBinError(xi,yi+1) + h_ch_bk->GetBinContent(xi,yi+1))
            {
                h_ch_clone_bk->SetBinContent(xi,yi,h_ch_bk->GetBinContent(xi,yi)/time_bk);
                h_ch_clone_bk->SetBinError(xi,yi,h_ch_bk->GetBinError(xi,yi)/time_bk);
            }

            if(h_ch_ch->GetBinContent(xi,yi) - factor_bad_pixels*h_ch_ch->GetBinError(xi,yi)
                    < factor_bad_pixels*h_ch_ch->GetBinError(xi,yi+1) + h_ch_ch->GetBinContent(xi,yi+1))
            {
                h_ch_clone_ch->SetBinContent(xi,yi,h_ch_ch->GetBinContent(xi,yi)/time_ch);
                h_ch_clone_ch->SetBinError(xi,yi,h_ch_ch->GetBinError(xi,yi)/time_ch);
            }
        }
    }

    TH1D* h_ch_clone_projectionY_am = h_ch_clone_am->ProjectionY("h_ch_clone_projectionY_am");
    TH1D* h_ch_clone_projectionY_bk = h_ch_clone_bk->ProjectionY("h_ch_clone_projectionY_bk");
    TH1D* h_ch_clone_projectionY_ch = h_ch_clone_ch->ProjectionY("h_ch_clone_projectionY_ch");

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1","h_ch_clone",1500,500);
    c1->Divide(3,1);

    c1->cd(1);
    h_ch_clone_bk->SetTitle("2D interated image, norm. on 1 sec. CRY2 in CH");
    h_ch_clone_bk->SetMaximum(1000);
    h_ch_clone_bk->GetXaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_bk->GetYaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_bk->Draw("colz");

    c1->cd(2);
    h_ch_clone_am->SetTitle("2D interated image, norm. on 1 sec. CRY2 in CH & CRY4 in AM");
    h_ch_clone_am->SetMaximum(1000);
    h_ch_clone_am->GetXaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_am->GetYaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_am->Draw("colz");

    c1->cd(3);
    h_ch_clone_ch->SetTitle("2D interated image, norm. on 1 sec. CRY2 in CH & CRY4 in CH");
    h_ch_clone_ch->SetMaximum(1000);
    h_ch_clone_ch->GetXaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_ch->GetYaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_ch->Draw("colz");

    Double_t integral_am, integral_err_am;
    Double_t integral_bk, integral_err_bk;
    Double_t integral_ch, integral_err_ch;

    integral_am = h_ch_clone_projectionY_am->IntegralAndError(1,N_PIXELS/2,integral_err_am);
    integral_bk = h_ch_clone_projectionY_bk->IntegralAndError(1,N_PIXELS/2,integral_err_bk);
    integral_ch = h_ch_clone_projectionY_ch->IntegralAndError(1,N_PIXELS/2,integral_err_ch);

    cout<<"--> integral_am = "<<integral_am<<" +/- "<<integral_err_am<<endl;
    cout<<"--> integral_bk = "<<integral_bk<<" +/- "<<integral_err_bk<<endl;
    cout<<"--> integral_ch = "<<integral_ch<<" +/- "<<integral_err_ch<<endl;

    for(Int_t i = 1; i <= N_PIXELS/2; i++)
    {
        h_ch_clone_projectionY_am->SetBinError(i,TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_am->GetBinError(i)/integral_am,2) +
                                                          TMath::Power(h_ch_clone_projectionY_am->GetBinContent(i)*integral_err_am/(integral_am*integral_am),2)));
        h_ch_clone_projectionY_am->SetBinContent(i,h_ch_clone_projectionY_am->GetBinContent(i)/integral_am);

        h_ch_clone_projectionY_bk->SetBinError(i,TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_bk->GetBinError(i)/integral_bk,2) +
                                                          TMath::Power(h_ch_clone_projectionY_bk->GetBinContent(i)*integral_err_bk/(integral_bk*integral_bk),2)));
        h_ch_clone_projectionY_bk->SetBinContent(i,h_ch_clone_projectionY_bk->GetBinContent(i)/integral_bk);

        h_ch_clone_projectionY_ch->SetBinError(i,TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_ch->GetBinError(i)/integral_ch,2) +
                                                          TMath::Power(h_ch_clone_projectionY_ch->GetBinContent(i)*integral_err_ch/(integral_ch*integral_ch),2)));
        h_ch_clone_projectionY_ch->SetBinContent(i,h_ch_clone_projectionY_ch->GetBinContent(i)/integral_ch);
    }

    TCanvas *c2 = new TCanvas("c2","h_ch_clone_projectionY",1500,1000);
    c2->cd();
    gPad->SetGrid();
    h_ch_clone_projectionY_ch->GetXaxis()->SetRange(1,N_PIXELS/2);

    h_ch_clone_projectionY_am->SetLineWidth(2);
    h_ch_clone_projectionY_bk->SetLineWidth(2);
    h_ch_clone_projectionY_ch->SetLineWidth(2);

    h_ch_clone_projectionY_am->SetLineColor(kRed);
    h_ch_clone_projectionY_bk->SetLineColor(kBlack);
    h_ch_clone_projectionY_ch->SetLineColor(kBlue);

    h_ch_clone_projectionY_am->SetMaximum(0.04);
    h_ch_clone_projectionY_bk->SetMaximum(0.04);
    h_ch_clone_projectionY_ch->SetMaximum(0.04);

    h_ch_clone_projectionY_ch->Draw();
    h_ch_clone_projectionY_bk->Draw("same");
    h_ch_clone_projectionY_am->Draw("same");


    TH1D* h_signal_bk = new TH1D("h_signal_bk","h_signal_bk",N_PIXELS/2,1,N_PIXELS/2);
    TH1D* h_signal_am = new TH1D("h_signal_am","h_signal_am",N_PIXELS/2,1,N_PIXELS/2);

    Double_t val, val_err;
    for(Int_t i = 1; i <= N_PIXELS/2; i++)
    {
        val = h_ch_clone_projectionY_ch->GetBinContent(i) - h_ch_clone_projectionY_bk->GetBinContent(i);
        val_err = TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_ch->GetBinError(i),2) + TMath::Power(h_ch_clone_projectionY_bk->GetBinError(i),2));

//                if(val < 0) val = 0;
        h_signal_bk->SetBinContent(i,val);
        h_signal_bk->SetBinError(i,val_err);

        val = h_ch_clone_projectionY_ch->GetBinContent(i) - h_ch_clone_projectionY_am->GetBinContent(i);
        val_err = TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_ch->GetBinError(i),2) + TMath::Power(h_ch_clone_projectionY_am->GetBinError(i),2));

//                if(val < 0) val = 0;
        h_signal_am->SetBinContent(i,val);
        h_signal_am->SetBinError(i,val_err);
    }

    TCanvas *c3 = new TCanvas("c3","h_ch_clone_projectionY_bk_am",1500,1000);
    c3->cd();
    gPad->SetGrid();
    h_signal_bk->SetLineColor(kBlack);
    h_signal_am->SetLineColor(kRed);
    h_signal_bk->SetMaximum(0.006);
    h_signal_am->SetMaximum(0.006);
    h_signal_bk->SetMinimum(-0.017);
    h_signal_am->SetMinimum(-0.017);
    h_signal_bk->Draw();
    h_signal_am->Draw("same");

    TH1D* h_signal_am_clone = (TH1D*)h_signal_am->Clone("h_signal_am_clone");
    TH1D* h_signal_bk_clone = (TH1D*)h_signal_bk->Clone("h_signal_bk_clone");

    TF1* fit_func_am = new TF1("fit_func_am",fit_ch_func_name,min_x_ch_fit,max_x_ch_fit);
    h_signal_am_clone->Fit(fit_func_am,"RQ0+");
    cout<<"--> Chi2/NDF AM = "<<fit_func_am->GetChisquare()/fit_func_am->GetNDF()<<endl;

    TF1* fit_func_bk = new TF1("fit_func_bk",fit_ch_func_name,min_x_ch_fit,max_x_ch_fit);
    h_signal_bk_clone->Fit(fit_func_bk,"RQ0+");
    cout<<"--> Chi2/NDF BK = "<<fit_func_bk->GetChisquare()/fit_func_bk->GetNDF()<<endl;


    TCanvas *c4 = new TCanvas("c4","h_ch_clone_projectionY_final",1000,1000);
    c4->Divide(1,2);

    h_signal_am_clone->GetXaxis()->SetRange(50,150);
    h_signal_am_clone->SetMaximum(0.003);
    h_signal_am_clone->SetMinimum(-0.0005);
    h_signal_bk_clone->GetXaxis()->SetRange(50,150);
    h_signal_bk_clone->SetMaximum(0.003);
    h_signal_bk_clone->SetMinimum(-0.0005);
    h_signal_am_clone->SetLineWidth(2);
    h_signal_bk_clone->SetLineWidth(2);
    h_signal_am_clone->SetLineColor(kRed);
    h_signal_bk_clone->SetLineColor(kBlack);

    c4->cd(1);
    gPad->SetGrid();
    h_signal_bk_clone->Draw();

    Int_t min_x_bk = fit_func_bk->GetParameter(1) - num_sigma_x_fit_ch*fit_func_bk->GetParameter(2);
    Int_t max_x_bk = fit_func_bk->GetParameter(1) + num_sigma_x_fit_ch*fit_func_bk->GetParameter(2);

    TLine* line_1_bk = new TLine(min_x_bk,0,min_x_bk,h_signal_bk_clone->GetMaximum());
    TLine* line_2_bk = new TLine(max_x_bk,0,max_x_bk,h_signal_bk_clone->GetMaximum());

    line_1_bk->SetLineColor(kGreen+2);
    line_2_bk->SetLineColor(kGreen+2);
    line_1_bk->SetLineWidth(4);
    line_2_bk->SetLineWidth(4);
    line_1_bk->Draw("same");
    line_2_bk->Draw("same");

    min_x_dch_ln = max_x_bk + 5;
    TLine* line_3_bk = new TLine(min_x_dch_ln,0,min_x_dch_ln,h_signal_bk_clone->GetMaximum());
    TLine* line_4_bk = new TLine(max_x_dch_ln,0,max_x_dch_ln,h_signal_bk_clone->GetMaximum());

    line_3_bk->SetLineColor(kRed+2);
    line_4_bk->SetLineColor(kRed+2);
    line_3_bk->SetLineWidth(4);
    line_4_bk->SetLineWidth(4);
    line_3_bk->Draw("same");
    line_4_bk->Draw("same");

    Double_t counts_ch_bk = 0.0, counts_ch_err_bk = 0.0, counts_dch_bk = 0.0, counts_dch_err_bk = 0.0;
    for(Int_t i = min_x_bk; i <= max_x_bk; i++)
    {
        counts_ch_bk += h_signal_bk->GetBinContent(i);
        counts_ch_err_bk += TMath::Power(h_signal_bk->GetBinError(i),2);
    }
    Double_t dch_length_bk = 0.0;
    for(Int_t i = min_x_dch_ln; i <= max_x_dch_ln; i++)
    {
        counts_dch_bk += h_signal_bk->GetBinContent(i);
        counts_dch_err_bk += TMath::Power(h_signal_bk->GetBinError(i),2);
        dch_length_bk += PIXEL_SIZE; // [mm]
    }
    counts_ch_err_bk = TMath::Sqrt(counts_ch_err_bk);
    counts_dch_err_bk = TMath::Sqrt(counts_dch_err_bk)/dch_length_bk;
    counts_dch_bk = counts_dch_bk/dch_length_bk;


    cout<<endl<<"--> Crystal info: "<<crystal_info<<endl<<endl;

    cout<<"--> Counts CH(BKG): "<<counts_ch_bk<<" +/- "<<counts_ch_err_bk<<endl;
    cout<<"--> Counts DCH(BKG): "<<counts_dch_bk<<" +/- "<<counts_dch_err_bk<<" [mm-1]"<<endl;
    cout<<"--> DCH length(BKG): "<<dch_length_bk<<" [mm]"<<endl;
    cout<<"--> Ratio CH/DCH(BKG): "<<counts_ch_bk/counts_dch_bk<<" +/- "<<TMath::Sqrt(TMath::Power(counts_ch_err_bk/counts_dch_bk,2) +
                                                                           TMath::Power(counts_ch_bk*counts_dch_err_bk/(counts_dch_bk*counts_dch_bk),2))<<endl;

    c4->cd(2);
    gPad->SetGrid();
    h_signal_am_clone->Draw();

    Int_t min_x_am = fit_func_am->GetParameter(1) - num_sigma_x_fit_ch*fit_func_am->GetParameter(2);
    Int_t max_x_am = fit_func_am->GetParameter(1) + num_sigma_x_fit_ch*fit_func_am->GetParameter(2);

    TLine* line_1_am = new TLine(min_x_am,0,min_x_am,h_signal_am_clone->GetMaximum());
    TLine* line_2_am = new TLine(max_x_am,0,max_x_am,h_signal_am_clone->GetMaximum());

    line_1_am->SetLineColor(kGreen+2);
    line_2_am->SetLineColor(kGreen+2);
    line_1_am->SetLineWidth(4);
    line_2_am->SetLineWidth(4);
    line_1_am->Draw("same");
    line_2_am->Draw("same");

    min_x_dch_ln = max_x_am + 5;
    TLine* line_3_am = new TLine(min_x_dch_ln,0,min_x_dch_ln,h_signal_am_clone->GetMaximum());
    TLine* line_4_am = new TLine(max_x_dch_ln,0,max_x_dch_ln,h_signal_am_clone->GetMaximum());

    line_3_am->SetLineColor(kRed+2);
    line_4_am->SetLineColor(kRed+2);
    line_3_am->SetLineWidth(4);
    line_4_am->SetLineWidth(4);
    line_3_am->Draw("same");
    line_4_am->Draw("same");

    Double_t counts_ch_am = 0.0, counts_ch_err_am = 0.0, counts_dch_am = 0.0, counts_dch_err_am = 0.0;
    for(Int_t i = min_x_am; i <= max_x_am; i++)
    {
        counts_ch_am += h_signal_am->GetBinContent(i);
        counts_ch_err_am += TMath::Power(h_signal_am->GetBinError(i),2);
    }
    Double_t dch_length_am = 0.0;
    for(Int_t i = min_x_dch_ln; i <= max_x_dch_ln; i++)
    {
        counts_dch_am += h_signal_am->GetBinContent(i);
        counts_dch_err_am += TMath::Power(h_signal_am->GetBinError(i),2);
        dch_length_am += PIXEL_SIZE; // [mm]
    }
    counts_ch_err_am = TMath::Sqrt(counts_ch_err_am);
    counts_dch_err_am = TMath::Sqrt(counts_dch_err_am)/dch_length_am;
    counts_dch_am = counts_dch_am/dch_length_am;


    cout<<"--> Counts CH(AM): "<<counts_ch_am<<" +/- "<<counts_ch_err_am<<endl;
    cout<<"--> Counts DCH(AM): "<<counts_dch_am<<" +/- "<<counts_dch_err_am<<" [mm-1]"<<endl;
    cout<<"--> DCH length(AM): "<<dch_length_am<<" [mm]"<<endl;
    cout<<"--> Ratio CH/DCH(AM): "<<counts_ch_am/counts_dch_am<<" +/- "<<TMath::Sqrt(TMath::Power(counts_ch_err_am/counts_dch_am,2) +
                                                                           TMath::Power(counts_ch_am*counts_dch_err_am/(counts_dch_am*counts_dch_am),2))<<endl;

    return 0;
}
