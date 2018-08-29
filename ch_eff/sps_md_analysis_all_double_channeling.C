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

int sps_md_analysis_all_double_channeling()
{

    //----------------------------------------------//
    //---------------- 17.10.2017 ------------------//
    //----------------------------------------------//
    const Int_t nPosition = 8;
    TString file_name_bk = "MD_2017_10_17_CRY2CH_HISTO_RP1I_RUN_1.root";
    TString file_name_ch[] = {
        "MD_2017_10_17_CRY2CHCRY4CHP1_HISTO_RP1I_RUN_1.root",
        "MD_2017_10_17_CRY2CHCRY4CHP2_HISTO_RP1I_RUN_1.root",
        "MD_2017_10_17_CRY2CHCRY4CHP3_HISTO_RP1I_RUN_1.root",
        "MD_2017_10_17_CRY2CHCRY4CHP4_HISTO_RP1I_RUN_1.root",
        "MD_2017_10_17_CRY2CHCRY4CHP5_HISTO_RP1I_RUN_1.root",
        "MD_2017_10_17_CRY2CHCRY4CHP6_HISTO_RP1I_RUN_1.root",
        "MD_2017_10_17_CRY2CHCRY4CHP7_HISTO_RP1I_RUN_1.root",
        "MD_2017_10_17_CRY2CHCRY4CHP8_HISTO_RP1I_RUN_1.root"
    };

    TString fit_ch_func_name = "gaus";
    Int_t min_x_ch_fit = 50;
    Int_t max_x_ch_fit = 150;
    Int_t min_x_dch_ln = 120;
    Int_t max_x_dch_ln = 145;
    Int_t num_sigma_x_fit_ch = 3;
    Float_t factor_bad_pixels = 10000.0;
    TString crystal_info = "CRYSTAL2 & CRYSTAL4 [17.10.2017]";
    //----------------------------------------------//
    TFile *_file_bk = TFile::Open(file_name_bk.Data());
    TH2D* h_evn_bk = (TH2D*)_file_bk->Get("h_10");
    Double_t time_bk = h_evn_bk->GetEntries()*0.000246; // interated time
    cout<<"--> time_bk = "<<time_bk<<" [sec]"<<endl;
    TH2D* h_ch_bk = (TH2D*)_file_bk->Get("h_1");
    TH2D* h_ch_clone_bk = new TH2D("h_ch_clone_bk","h_ch_clone_bk",N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);

    for(Int_t xi = 1; xi <= N_PIXELS/2; xi++)
    {
        for(Int_t yi = 1; yi <= N_PIXELS/2; yi++)
        {
            if(h_ch_bk->GetBinContent(xi,yi) - factor_bad_pixels*h_ch_bk->GetBinError(xi,yi)
                    < factor_bad_pixels*h_ch_bk->GetBinError(xi,yi+1) + h_ch_bk->GetBinContent(xi,yi+1))
            {
                h_ch_clone_bk->SetBinContent(xi,yi,h_ch_bk->GetBinContent(xi,yi)/time_bk);
                h_ch_clone_bk->SetBinError(xi,yi,h_ch_bk->GetBinError(xi,yi)/time_bk);
            }
        }
    }

    TH1D* h_ch_clone_projectionY_bk = h_ch_clone_bk->ProjectionY("h_ch_clone_projectionY_bk");

    TFile *_file_ch[nPosition];
    TH2D* h_evn_ch[nPosition];
    Double_t time_ch[nPosition];
    TH2D* h_ch_ch[nPosition];
    TH2D* h_ch_clone_ch[nPosition];
    TH1D* h_ch_clone_projectionY_ch[nPosition];

    for(Int_t i = 0; i < nPosition; i++)
    {
        _file_ch[i] = TFile::Open(file_name_ch[i].Data());
        h_evn_ch[i] = (TH2D*)_file_ch[i]->Get("h_10");
        time_ch[i] = h_evn_ch[i]->GetEntries()*0.000246; // interated time
        cout<<"--> time_ch["<<i<<"] = "<<time_ch[i]<<" [sec]"<<endl;
        h_ch_ch[i] = (TH2D*)_file_ch[i]->Get("h_1");
        TString name = "h_ch_clone_ch_"; name += i;
        h_ch_clone_ch[i] = new TH2D(name.Data(),name.Data(),N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);

        for(Int_t xi = 1; xi <= N_PIXELS/2; xi++)
        {
            for(Int_t yi = 1; yi <= N_PIXELS/2; yi++)
            {
                if(h_ch_ch[i]->GetBinContent(xi,yi) - factor_bad_pixels*h_ch_ch[i]->GetBinError(xi,yi)
                        < factor_bad_pixels*h_ch_ch[i]->GetBinError(xi,yi+1) + h_ch_ch[i]->GetBinContent(xi,yi+1))
                {
                    h_ch_clone_ch[i]->SetBinContent(xi,yi,h_ch_ch[i]->GetBinContent(xi,yi)/time_ch[i]);
                    h_ch_clone_ch[i]->SetBinError(xi,yi,h_ch_ch[i]->GetBinError(xi,yi)/time_ch[i]);
                }
            }
        }

        name = "h_ch_clone_projectionY_ch_"; name += i;
        h_ch_clone_projectionY_ch[i] = h_ch_clone_ch[i]->ProjectionY(name.Data());
    }

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1","h_ch_clone",1000,1000);
    c1->Divide(3,3);

    c1->cd(1);
    h_ch_clone_bk->SetTitle("2D interated image, norm. on 1 sec. CRY2 in CH");
    h_ch_clone_bk->SetMaximum(1000);
    h_ch_clone_bk->GetXaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_bk->GetYaxis()->SetRange(1,N_PIXELS/2);
    h_ch_clone_bk->Draw("colz");

    for(Int_t i = 0; i < nPosition; i++)
    {
        c1->cd(i+2);
        TString title = "CRY4 in CH position "; title += i+1;
        h_ch_clone_ch[i]->SetTitle(title.Data());
        h_ch_clone_ch[i]->SetMaximum(1000);
        h_ch_clone_ch[i]->GetXaxis()->SetRange(1,N_PIXELS/2);
        h_ch_clone_ch[i]->GetYaxis()->SetRange(1,N_PIXELS/2);
        h_ch_clone_ch[i]->Draw("colz");
    }

    Double_t integral_bk, integral_err_bk;
    Double_t integral_ch[nPosition], integral_err_ch[nPosition];

    integral_bk = h_ch_clone_projectionY_bk->IntegralAndError(1,N_PIXELS/2,integral_err_bk);
    cout<<"--> integral_bk = "<<integral_bk<<" +/- "<<integral_err_bk<<endl;

    for(Int_t i = 1; i <= N_PIXELS/2; i++)
    {
        h_ch_clone_projectionY_bk->SetBinError(i,TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_bk->GetBinError(i)/integral_bk,2) +
                                                          TMath::Power(h_ch_clone_projectionY_bk->GetBinContent(i)*integral_err_bk/(integral_bk*integral_bk),2)));
        h_ch_clone_projectionY_bk->SetBinContent(i,h_ch_clone_projectionY_bk->GetBinContent(i)/integral_bk);
    }

    for(Int_t i = 0; i < nPosition; i++)
    {
        integral_ch[i] = h_ch_clone_projectionY_ch[i]->IntegralAndError(1,N_PIXELS/2,integral_err_ch[i]);
        cout<<"--> integral_ch["<<i<<"] = "<<integral_ch[i]<<" +/- "<<integral_err_ch[i]<<endl;

        for(Int_t j = 1; j <= N_PIXELS/2; j++)
        {
            h_ch_clone_projectionY_ch[i]->SetBinError(j,TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_ch[i]->GetBinError(j)/integral_ch[i],2) +
                                                              TMath::Power(h_ch_clone_projectionY_ch[i]->GetBinContent(j)*integral_err_ch[i]/(integral_ch[i]*integral_ch[i]),2)));
            h_ch_clone_projectionY_ch[i]->SetBinContent(j,h_ch_clone_projectionY_ch[i]->GetBinContent(j)/integral_ch[i]);
        }
    }

    TCanvas *c2 = new TCanvas("c2","h_ch_clone_projectionY",1600,800);
    c2->Divide(4,2);
    h_ch_clone_projectionY_bk->SetLineColor(kBlack);
    h_ch_clone_projectionY_bk->SetLineWidth(2);

    for(Int_t i = 0; i < nPosition; i++)
    {
        c2->cd(i+1);
        gPad->SetGrid();
        h_ch_clone_projectionY_ch[i]->GetXaxis()->SetRange(1,N_PIXELS/2);
        h_ch_clone_projectionY_ch[i]->SetLineWidth(2);
        h_ch_clone_projectionY_ch[i]->SetLineColor(kBlue);
        h_ch_clone_projectionY_ch[i]->SetMaximum(0.04);
        TString title = "CRY4 in CH position "; title += i+1;
        h_ch_clone_projectionY_ch[i]->SetTitle(title.Data());
        h_ch_clone_projectionY_ch[i]->Draw();
        h_ch_clone_projectionY_bk->Draw("same");
    }

    TCanvas *c3 = new TCanvas("c3","h_ch_clone_projectionY_bk_am",1600,800);
    c3->Divide(4,2);

    TH1D* h_signal_bk[nPosition];
    for(Int_t j = 0; j < nPosition; j++)
    {
        TString name = "h_signal_bk_"; name += j;
        h_signal_bk[j] = new TH1D(name.Data(),name.Data(),N_PIXELS/2,1,N_PIXELS/2);

        Double_t val, val_err;
        for(Int_t i = 1; i <= N_PIXELS/2; i++)
        {
            val = h_ch_clone_projectionY_ch[j]->GetBinContent(i) - h_ch_clone_projectionY_bk->GetBinContent(i);
            val_err = TMath::Sqrt(TMath::Power(h_ch_clone_projectionY_ch[j]->GetBinError(i),2) + TMath::Power(h_ch_clone_projectionY_bk->GetBinError(i),2));

//            if(val < 0) val = 0;

            h_signal_bk[j]->SetBinContent(i,val);
            h_signal_bk[j]->SetBinError(i,val_err);
        }

        c3->cd(j+1);
        gPad->SetGrid();
        h_signal_bk[j]->SetLineColor(kBlack);
        h_signal_bk[j]->SetMaximum(0.010);
        h_signal_bk[j]->SetMinimum(-0.017);
        TString title = "CRY4 in CH position "; title += j+1;
        h_signal_bk[j]->SetTitle(title.Data());
        h_signal_bk[j]->Draw();
    }

    TH1D* h_signal_bk_clone[nPosition];
    TF1* fit_func_bk[nPosition];

    for(Int_t j = 0; j < nPosition; j++)
    {
        TString func_name = "fit_func_bk_"; func_name += j;
        fit_func_bk[j] = new TF1(func_name.Data(),fit_ch_func_name.Data(),min_x_ch_fit,max_x_ch_fit);
        TString name = "h_signal_bk_clone_"; name += j;
        h_signal_bk_clone[j] = (TH1D*)h_signal_bk[j]->Clone(name.Data());
        h_signal_bk_clone[j]->Fit(fit_func_bk[j],"RQ0+");
        cout<<"--> Chi2/NDF BK ["<<j<<"] = "<<fit_func_bk[j]->GetChisquare()/fit_func_bk[j]->GetNDF()<<endl;
    }

    TCanvas *c4 = new TCanvas("c4","h_ch_clone_projectionY_final",1600,800);
    c4->Divide(4,2);

    TLine* line_1_bk[nPosition];
    TLine* line_2_bk[nPosition];
    TLine* line_3_bk[nPosition];
    TLine* line_4_bk[nPosition];

    Double_t counts_ch_bk[nPosition] = {}, counts_ch_err_bk[nPosition] = {}, counts_dch_bk[nPosition] = {}, counts_dch_err_bk[nPosition] = {};
    Double_t dch_length_bk[nPosition] = {};

    cout<<endl<<"--> Crystal info: "<<crystal_info<<endl<<endl;

    for(Int_t j = 0; j < nPosition; j++)
    {
        h_signal_bk_clone[j]->GetXaxis()->SetRange(50,150);
        h_signal_bk_clone[j]->SetMaximum(0.005);
        h_signal_bk_clone[j]->SetMinimum(-0.0005);
        h_signal_bk_clone[j]->SetLineWidth(2);
        h_signal_bk_clone[j]->SetLineColor(kBlack);

        c4->cd(j+1);
        gPad->SetGrid();
        h_signal_bk_clone[j]->Draw();

        Int_t min_x_bk = fit_func_bk[j]->GetParameter(1) - num_sigma_x_fit_ch*fit_func_bk[j]->GetParameter(2);
        Int_t max_x_bk = fit_func_bk[j]->GetParameter(1) + num_sigma_x_fit_ch*fit_func_bk[j]->GetParameter(2);

        line_1_bk[j] = new TLine(min_x_bk,0,min_x_bk,h_signal_bk_clone[j]->GetMaximum());
        line_2_bk[j] = new TLine(max_x_bk,0,max_x_bk,h_signal_bk_clone[j]->GetMaximum());

        line_1_bk[j]->SetLineColor(kGreen+2);
        line_2_bk[j]->SetLineColor(kGreen+2);
        line_1_bk[j]->SetLineWidth(4);
        line_2_bk[j]->SetLineWidth(4);
        line_1_bk[j]->Draw("same");
        line_2_bk[j]->Draw("same");

        min_x_dch_ln = max_x_bk + 5;
        line_3_bk[j] = new TLine(min_x_dch_ln,0,min_x_dch_ln,h_signal_bk_clone[j]->GetMaximum());
        line_4_bk[j] = new TLine(max_x_dch_ln,0,max_x_dch_ln,h_signal_bk_clone[j]->GetMaximum());

        line_3_bk[j]->SetLineColor(kRed+2);
        line_4_bk[j]->SetLineColor(kRed+2);
        line_3_bk[j]->SetLineWidth(4);
        line_4_bk[j]->SetLineWidth(4);
        line_3_bk[j]->Draw("same");
        line_4_bk[j]->Draw("same");

        for(Int_t i = min_x_bk; i <= max_x_bk; i++)
        {
            counts_ch_bk[j] += h_signal_bk[j]->GetBinContent(i);
            counts_ch_err_bk[j] += TMath::Power(h_signal_bk[j]->GetBinError(i),2);
        }

        for(Int_t i = min_x_dch_ln; i <= max_x_dch_ln; i++)
        {
            counts_dch_bk[j] += h_signal_bk[j]->GetBinContent(i);
            counts_dch_err_bk[j] += TMath::Power(h_signal_bk[j]->GetBinError(i),2);
            dch_length_bk[j] += PIXEL_SIZE; // [mm]
        }

        counts_ch_err_bk[j] = TMath::Sqrt(counts_ch_err_bk[j]);
        counts_dch_err_bk[j] = TMath::Sqrt(counts_dch_err_bk[j])/dch_length_bk[j];
        counts_dch_bk[j] = counts_dch_bk[j]/dch_length_bk[j];

        cout<<endl<<"Position "<<j+1<<endl;
        cout<<"Counts CH(BKG): "<<counts_ch_bk[j]<<" +/- "<<counts_ch_err_bk[j]<<endl;
        cout<<"Counts DCH(BKG): "<<counts_dch_bk[j]<<" +/- "<<counts_dch_err_bk[j]<<" [mm-1]"<<endl;
        cout<<"DCH length(BKG): "<<dch_length_bk[j]<<" [mm]"<<endl;
        cout<<"Ratio CH/DCH(BKG): "<<counts_ch_bk[j]/counts_dch_bk[j]<<" +/- "<<TMath::Sqrt(TMath::Power(counts_ch_err_bk[j]/counts_dch_bk[j],2) +
                                                                                          TMath::Power(counts_ch_bk[j]*counts_dch_err_bk[j]/(counts_dch_bk[j]*counts_dch_bk[j]),2))<<endl;
    }

    cout<<endl;
    return 0;
}
