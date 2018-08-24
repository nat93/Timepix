//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     Data analysis for SPS MD data 15/08/2018
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

int sps_md_analysis()
{
    TFile *_file = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_RUN_1.root");
    TH2D* h_ch = (TH2D*)_file->Get("h_1");
    TH2D* h_ch_clone = new TH2D("h_ch_clone","h_ch_clone",256,1,256,256,1,256);

    Double_t factor = 10.0;

    for(Int_t xi = 3; xi < 256-3; xi++)
    {
        for(Int_t yi = 3; yi < 256-3; yi++)
        {
            if(h_ch->GetBinContent(xi,yi) - factor*h_ch->GetBinError(xi,yi) < factor*h_ch->GetBinError(xi,yi+1) + h_ch->GetBinContent(xi,yi+1))
            {
                h_ch_clone->SetBinContent(xi,yi,h_ch->GetBinContent(xi,yi));
                h_ch_clone->SetBinError(xi,yi,h_ch->GetBinError(xi,yi));
            }
        }
    }

    TH1D* h_ch_clone_projectionY = h_ch_clone->ProjectionY("h_ch_clone_projectionY");

    TCanvas *c1 = new TCanvas("c1","h_ch",900,450);
    gStyle->SetOptStat(0);
    c1->Divide(2,1);

    c1->cd(1);
    h_ch->GetXaxis()->SetRange(1,256);
    h_ch->GetYaxis()->SetRange(1,256);
    h_ch->Draw("colz");

    c1->cd(2);
    h_ch_clone->GetXaxis()->SetRange(1,256);
    h_ch_clone->GetYaxis()->SetRange(1,256);
    h_ch_clone->Draw("colz");

    TCanvas *c3 = new TCanvas("c3","h_ch_clone_projectionY",1800,900);
    gStyle->SetOptStat(0);
    c3->Divide(1,2);
    c3->cd(1);
    h_ch_clone_projectionY->GetXaxis()->SetRange(1,256);
    h_ch_clone_projectionY->SetMaximum(0.015);
    h_ch_clone_projectionY->SetMinimum(-0.001);
    Double_t integral, integral_err;
    integral = h_ch_clone_projectionY->IntegralAndError(1,256,integral_err);
    for(Int_t i = 1; i <= 256; i++)
    {
        h_ch_clone_projectionY->SetBinError(i,TMath::Sqrt(TMath::Power(h_ch_clone_projectionY->GetBinError(i)/integral,2) +
                                                          TMath::Power(h_ch_clone_projectionY->GetBinContent(i)*integral_err/(integral*integral),2)));
        h_ch_clone_projectionY->SetBinContent(i,h_ch_clone_projectionY->GetBinContent(i)/integral);
    }
    h_ch_clone_projectionY->SetLineWidth(2);
    h_ch_clone_projectionY->Draw();

    TF1* fit_func = new TF1("fit_func","pol3",20,90);
    fit_func->SetLineColor(kMagenta);
    fit_func->SetLineWidth(4);
    h_ch_clone_projectionY->Fit(fit_func,"R0Q+");
    cout<<"--> Chi2/NDF I = "<<fit_func->GetChisquare()/fit_func->GetNDF()<<endl;
    
    TH1D* h_bkg = new TH1D("h_bkg","h_bkg [pol3]",256,1,256);
    TH1D* h_signal = new TH1D("h_signal","h_signal",256,1,256);
    for(Int_t i = 1; i <= 256; i++)
    {
        h_bkg->SetBinContent(i,fit_func->Eval(i));

        Double_t val = h_ch_clone_projectionY->GetBinContent(i) - fit_func->Eval(i);
        Double_t val_err = h_ch_clone_projectionY->GetBinError(i);
        if(val < 0) val = 0;
        h_signal->SetBinContent(i,val);
        h_signal->SetBinError(i,val_err);
    }
    h_bkg->SetLineColor(kBlack);
    h_bkg->Draw("same");
    c3->cd(2);
    h_signal->GetXaxis()->SetRange(1,256);
    h_signal->SetMaximum(0.015);
    h_signal->SetMinimum(-0.001);
    h_signal->SetLineWidth(2);
    h_signal->Draw();

    TF1* fit_func_2 = new TF1("fit_func_2","gaus",150,170);
    fit_func_2->SetLineColor(kRed);
    fit_func_2->SetLineWidth(4);
    h_signal->Fit(fit_func_2,"R0Q+");
    cout<<"--> Chi2/NDF II = "<<fit_func_2->GetChisquare()/fit_func_2->GetNDF()<<endl;

    TH1D* h_signal_2 = new TH1D("h_signal_2","h_signal_2",256,1,256);
    for(Int_t i = 1; i <= 256; i++)
    {
        h_signal_2->SetBinContent(i,fit_func_2->Eval(i));
    }
    h_signal_2->SetLineColor(kBlack);
    h_signal_2->Draw("same");

    Int_t min_x = fit_func_2->GetParameter(1) - 3.0*fit_func_2->GetParameter(2);
    Int_t max_x = fit_func_2->GetParameter(1) + 3.0*fit_func_2->GetParameter(2);

    TLine* line_1 = new TLine(min_x,0,min_x,h_signal->GetMaximum());
    TLine* line_2 = new TLine(max_x,0,max_x,h_signal->GetMaximum());

    line_1->SetLineColor(kGreen+2);
    line_2->SetLineColor(kGreen+2);
    line_1->SetLineWidth(4);
    line_2->SetLineWidth(4);
    line_1->Draw("same");
    line_2->Draw("same");

    TLine* line_3 = new TLine(205,0,205,h_signal->GetMaximum());
    TLine* line_4 = new TLine(220,0,220,h_signal->GetMaximum());

    line_3->SetLineColor(kBlue+2);
    line_4->SetLineColor(kBlue+2);
    line_3->SetLineWidth(4);
    line_4->SetLineWidth(4);
    line_3->Draw("same");
    line_4->Draw("same");

    Double_t counts_ch = 0.0, counts_ch_err = 0.0, counts_dch = 0.0, counts_dch_err = 0.0;
    for(Int_t i = min_x; i <= max_x; i++)
    {
        counts_ch += h_signal->GetBinContent(i);
        counts_ch_err += TMath::Power(h_signal->GetBinError(i),2);
    }
    Double_t dch_length = 0.0;
    for(Int_t i = 205; i <= 220; i++)
    {
        counts_dch += h_signal->GetBinContent(i);
        counts_dch_err += TMath::Power(h_signal->GetBinError(i),2);
        dch_length += PIXEL_SIZE; // [mm]
    }

    counts_ch_err = TMath::Sqrt(counts_ch_err);
    counts_dch_err = TMath::Sqrt(counts_dch_err)/dch_length;
    counts_dch = counts_dch/dch_length;
    cout<<"--> Counts CH: "<<counts_ch<<" +/- "<<counts_ch_err<<endl;
    cout<<"--> Counts DCH: "<<counts_dch<<" +/- "<<counts_dch_err<<" [mm-1]"<<endl;
    cout<<"--> DCH length: "<<dch_length<<" [mm]"<<endl;
    cout<<"--> Ratio CH/DCH: "<<counts_ch/counts_dch<<" +/- "<<TMath::Sqrt(TMath::Power(counts_ch_err/counts_dch,2) +
                                                                           TMath::Power(counts_ch*counts_dch_err/(counts_dch*counts_dch),2))<<endl;

    return 0;
}
