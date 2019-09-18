//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     Data analysis for SPS MD data (RP0I,RP0E,RP3E calibration)
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

int sps_md_analysis_3()
{
    //----------------------------------------------//
    //---------------- 15.08.2018 ------------------//
    //----------------------------------------------//
    const Int_t nFiles = 8;
    gStyle->SetOptStat(0);

    TString file_name_rp1i[] = {
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_295_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_296_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_297_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_298_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_299_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_300_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_301_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_302_RUN_1.root"     // dTHL = 100
    };
    TString file_name_rp0i[] = {
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_295_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_296_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_297_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_298_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_299_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_300_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_301_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_302_RUN_1.root"     // dTHL = 40
    };
    TString file_name_rp0e[] = {
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_295_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_296_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_297_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_298_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_299_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_300_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_301_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0E_302_RUN_1.root"     // dTHL = 40
    };
    TString file_name_rp3e[] = {
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_295_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_296_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_297_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_298_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_299_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_300_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_301_RUN_1.root",    // dTHL = 100
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP3E_302_RUN_1.root"     // dTHL = 40
    };
    Double_t bias[] = {0.0,10.0,20.0,40.0,60.0,80.0,100.0,40};

    Float_t factor_bad_pixels_rp0i = 10;
    Float_t factor_bad_pixels_rp0e = 10;
    Float_t factor_bad_pixels_rp3e = 10;
    Float_t factor_bad_pixels_rp1i = 10;

    Double_t integral_rp0i[nFiles] = {}, integral_err_rp0i[nFiles] = {};
    Double_t integral_rp0e[nFiles] = {}, integral_err_rp0e[nFiles] = {};
    Double_t integral_rp3e[nFiles] = {}, integral_err_rp3e[nFiles] = {};
    Double_t integral_rp1i[nFiles] = {}, integral_err_rp1i[nFiles] = {};

    Double_t integral_norm_rp0i[nFiles] = {}, integral_err_norm_rp0i[nFiles] = {};
    Double_t integral_norm_rp0e[nFiles] = {}, integral_err_norm_rp0e[nFiles] = {};
    Double_t integral_norm_rp3e[nFiles] = {}, integral_err_norm_rp3e[nFiles] = {};

    TH2D* h_rp0i_clone[nFiles];
    TH2D* h_rp0e_clone[nFiles];
    TH2D* h_rp3e_clone[nFiles];
    TH2D* h_rp1i_clone[nFiles];

    TCanvas* c1_rp0i = new TCanvas("c1_rp0i","c1_rp0i",1000,500);
    TCanvas* c1_rp0e = new TCanvas("c1_rp0e","c1_rp0e",1000,500);
    TCanvas* c1_rp3e = new TCanvas("c1_rp3e","c1_rp3e",1000,500);
    TCanvas* c1_rp1i = new TCanvas("c1_rp1i","c1_rp1i",1000,500);

    c1_rp0i->Divide(4,2);
    c1_rp0e->Divide(4,2);
    c1_rp3e->Divide(4,2);
    c1_rp1i->Divide(4,2);

    for(Int_t ifl = 0; ifl < nFiles; ifl++)
    {
        cout<<endl;
        cout<<"--> Input file RP0I: "<<file_name_rp0i[ifl]<<endl;
        cout<<"--> Input file RP0E: "<<file_name_rp0e[ifl]<<endl;
        cout<<"--> Input file RP3E: "<<file_name_rp3e[ifl]<<endl;
        cout<<"--> Input file RP1I: "<<file_name_rp1i[ifl]<<endl;

        TFile *_file_rp0i = TFile::Open(file_name_rp0i[ifl].Data());
        TFile *_file_rp0e = TFile::Open(file_name_rp0e[ifl].Data());
        TFile *_file_rp3e = TFile::Open(file_name_rp3e[ifl].Data());
        TFile *_file_rp1i = TFile::Open(file_name_rp1i[ifl].Data());

        TH2D* h_rp0i = (TH2D*)_file_rp0i->Get("h_1");
        TH2D* h_rp0e = (TH2D*)_file_rp0e->Get("h_1");
        TH2D* h_rp3e = (TH2D*)_file_rp3e->Get("h_1");
        TH2D* h_rp1i = (TH2D*)_file_rp1i->Get("h_1");

        TString hist_name_rp0i = "h_rp0i_clone_"; hist_name_rp0i += ifl;
        TString hist_name_rp0e = "h_rp0e_clone_"; hist_name_rp0e += ifl;
        TString hist_name_rp3e = "h_rp3e_clone_"; hist_name_rp3e += ifl;
        TString hist_name_rp1i = "h_rp1i_clone_"; hist_name_rp1i += ifl;

        h_rp0i_clone[ifl] = new TH2D(hist_name_rp0i.Data(),hist_name_rp0i.Data(),N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);
        h_rp0e_clone[ifl] = new TH2D(hist_name_rp0e.Data(),hist_name_rp0e.Data(),N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);
        h_rp3e_clone[ifl] = new TH2D(hist_name_rp3e.Data(),hist_name_rp3e.Data(),N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);
        h_rp1i_clone[ifl] = new TH2D(hist_name_rp1i.Data(),hist_name_rp1i.Data(),N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);

        for(Int_t xi = 3; xi < N_PIXELS/2-3; xi++)
        {
            for(Int_t yi = 3; yi < N_PIXELS/2-3; yi++)
            {
                if(h_rp0i->GetBinContent(xi,yi) < h_rp0i->GetBinContent(xi,yi+1) + factor_bad_pixels_rp0i*h_rp0i->GetBinError(xi,yi+1) &&
                        h_rp0i->GetBinContent(xi,yi) > h_rp0i->GetBinContent(xi,yi+1) - factor_bad_pixels_rp0i*h_rp0i->GetBinError(xi,yi+1))
                {
                    h_rp0i_clone[ifl]->SetBinContent(xi,yi,h_rp0i->GetBinContent(xi,yi));
                    h_rp0i_clone[ifl]->SetBinError(xi,yi,h_rp0i->GetBinError(xi,yi));
                }

                if(h_rp0e->GetBinContent(xi,yi) < h_rp0e->GetBinContent(xi,yi+1) + factor_bad_pixels_rp0e*h_rp0e->GetBinError(xi,yi+1) &&
                        h_rp0e->GetBinContent(xi,yi) > h_rp0e->GetBinContent(xi,yi+1) - factor_bad_pixels_rp0e*h_rp0e->GetBinError(xi,yi+1))
                {
                    h_rp0e_clone[ifl]->SetBinContent(xi,yi,h_rp0e->GetBinContent(xi,yi));
                    h_rp0e_clone[ifl]->SetBinError(xi,yi,h_rp0e->GetBinError(xi,yi));
                }

                if(h_rp3e->GetBinContent(xi,yi) < h_rp3e->GetBinContent(xi,yi+1) + factor_bad_pixels_rp3e*h_rp3e->GetBinError(xi,yi+1) &&
                        h_rp3e->GetBinContent(xi,yi) > h_rp3e->GetBinContent(xi,yi+1) - factor_bad_pixels_rp3e*h_rp3e->GetBinError(xi,yi+1))
                {
                    h_rp3e_clone[ifl]->SetBinContent(xi,yi,h_rp3e->GetBinContent(xi,yi));
                    h_rp3e_clone[ifl]->SetBinError(xi,yi,h_rp3e->GetBinError(xi,yi));
                }

                if(h_rp1i->GetBinContent(xi,yi) < h_rp1i->GetBinContent(xi,yi+1) + factor_bad_pixels_rp1i*h_rp1i->GetBinError(xi,yi+1) &&
                        h_rp1i->GetBinContent(xi,yi) > h_rp1i->GetBinContent(xi,yi+1) - factor_bad_pixels_rp1i*h_rp1i->GetBinError(xi,yi+1))
                {
                    h_rp1i_clone[ifl]->SetBinContent(xi,yi,h_rp1i->GetBinContent(xi,yi));
                    h_rp1i_clone[ifl]->SetBinError(xi,yi,h_rp1i->GetBinError(xi,yi));
                }
            }
        }

        integral_rp0i[ifl] = h_rp0i_clone[ifl]->IntegralAndError(1,N_PIXELS/2,1,N_PIXELS/2,integral_err_rp0i[ifl]);
        integral_rp0e[ifl] = h_rp0e_clone[ifl]->IntegralAndError(1,N_PIXELS/2,1,N_PIXELS/2,integral_err_rp0e[ifl]);
        integral_rp3e[ifl] = h_rp3e_clone[ifl]->IntegralAndError(1,N_PIXELS/2,1,N_PIXELS/2,integral_err_rp3e[ifl]);
        integral_rp1i[ifl] = h_rp1i_clone[ifl]->IntegralAndError(1,N_PIXELS/2,1,N_PIXELS/2,integral_err_rp1i[ifl]);
        integral_norm_rp0i[ifl] = integral_rp0i[ifl]/integral_rp1i[ifl];
        integral_norm_rp0e[ifl] = integral_rp0e[ifl]/integral_rp1i[ifl];
        integral_norm_rp3e[ifl] = integral_rp3e[ifl]/integral_rp1i[ifl];
        integral_err_norm_rp0i[ifl] = TMath::Sqrt(TMath::Power(integral_err_rp0i[ifl]/integral_rp1i[ifl],2)
                                             + TMath::Power(integral_rp0i[ifl]*integral_err_rp1i[ifl]/(integral_rp1i[ifl]*integral_rp1i[ifl]),2));
        integral_err_norm_rp0e[ifl] = TMath::Sqrt(TMath::Power(integral_err_rp0e[ifl]/integral_rp1i[ifl],2)
                                             + TMath::Power(integral_rp0e[ifl]*integral_err_rp1i[ifl]/(integral_rp1i[ifl]*integral_rp1i[ifl]),2));
        integral_err_norm_rp3e[ifl] = TMath::Sqrt(TMath::Power(integral_err_rp3e[ifl]/integral_rp1i[ifl],2)
                                             + TMath::Power(integral_rp3e[ifl]*integral_err_rp1i[ifl]/(integral_rp1i[ifl]*integral_rp1i[ifl]),2));


        c1_rp0i->cd(ifl+1);
        h_rp0i_clone[ifl]->SetMaximum(5000);
        h_rp0i_clone[ifl]->Draw("colz");

        c1_rp0e->cd(ifl+1);
//        h_rp0e_clone[ifl]->SetMaximum(5000);
        h_rp0e_clone[ifl]->Draw("colz");

        c1_rp3e->cd(ifl+1);
//        h_rp3e_clone[ifl]->SetMaximum(5000);
        h_rp3e_clone[ifl]->Draw("colz");

        c1_rp1i->cd(ifl+1);
        h_rp1i_clone[ifl]->SetMaximum(500);
        h_rp1i_clone[ifl]->Draw("colz");
    }

    TGraphErrors* gr_rp0i = new TGraphErrors(nFiles-1,bias,integral_rp0i,0,integral_err_rp0i);
    TGraphErrors* gr_rp0e = new TGraphErrors(nFiles-1,bias,integral_rp0e,0,integral_err_rp0e);
    TGraphErrors* gr_rp3e = new TGraphErrors(nFiles-1,bias,integral_rp3e,0,integral_err_rp3e);
    TGraphErrors* gr_rp1i = new TGraphErrors(nFiles-1,bias,integral_rp1i,0,integral_err_rp1i);

    TGraphErrors* gr_norm_rp0i = new TGraphErrors(nFiles-1,bias,integral_norm_rp0i,0,integral_err_norm_rp0i);
    TGraphErrors* gr_norm_rp0e = new TGraphErrors(nFiles-1,bias,integral_norm_rp0e,0,integral_err_norm_rp0e);
    TGraphErrors* gr_norm_rp3e = new TGraphErrors(nFiles-1,bias,integral_norm_rp3e,0,integral_err_norm_rp3e);

    gr_rp0i->SetName("RP0I");   gr_rp0i->SetTitle("RP0I");
    gr_rp0e->SetName("RP0E");   gr_rp0e->SetTitle("RP0E");
    gr_rp3e->SetName("RP3E");   gr_rp3e->SetTitle("RP3E");
    gr_rp1i->SetName("RP1I");   gr_rp1i->SetTitle("RP1I");
    gr_norm_rp0i->SetName("RP0I/RP1I");   gr_norm_rp0i->SetTitle("RP0I/RP1I");
    gr_norm_rp0e->SetName("RP0E/RP1I");   gr_norm_rp0e->SetTitle("RP0E/RP1I");
    gr_norm_rp3e->SetName("RP3E/RP1I");   gr_norm_rp3e->SetTitle("RP3E/RP1I");

    gr_rp0i->SetLineColor(kRed);    gr_rp0i->SetLineWidth(2);
    gr_rp0e->SetLineColor(kRed);    gr_rp0e->SetLineWidth(2);
    gr_rp3e->SetLineColor(kRed);    gr_rp3e->SetLineWidth(2);
    gr_rp1i->SetLineColor(kBlue);    gr_rp1i->SetLineWidth(2);
    gr_norm_rp0i->SetLineColor(kBlack);    gr_norm_rp0i->SetLineWidth(2);
    gr_norm_rp0e->SetLineColor(kBlack);    gr_norm_rp0e->SetLineWidth(2);
    gr_norm_rp3e->SetLineColor(kBlack);    gr_norm_rp3e->SetLineWidth(2);

    TCanvas* c_comm = new TCanvas("c_comm","c_comm",1000,500);
    c_comm->Divide(4,1);
    c_comm->cd(1);
    gr_rp0i->Draw("APL");
    c_comm->cd(2);
    gr_rp0e->Draw("APL");
    c_comm->cd(3);
    gr_rp3e->Draw("APL");
    c_comm->cd(4);
    gr_rp1i->Draw("APL");

    TCanvas* c_norm_rp0i = new TCanvas("c_norm_rp0i","c_norm_rp0i",500,500);
    c_norm_rp0i->cd();
    gr_norm_rp0i->Draw("APL");

    TCanvas* c_norm_rp0e = new TCanvas("c_norm_rp0e","c_norm_rp0e",500,500);
    c_norm_rp0e->cd();
    gr_norm_rp0e->Draw("APL");

    TCanvas* c_norm_rp3e = new TCanvas("c_norm_rp3e","c_norm_rp3e",500,500);
    c_norm_rp3e->cd();
    gr_norm_rp3e->Draw("APL");

    return 0;
}
