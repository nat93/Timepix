//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     Data analysis for SPS MD data (RP1I calibration)
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

int sps_md_analysis_4()
{
    //----------------------------------------------//
    //---------------- 18.06.2018 ------------------//
    //----------------------------------------------//
    const Int_t nFiles = 1;
    gStyle->SetOptStat(0);

    TString file_name_rp1i[] = {
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_RUN_1.root"
    };
    TString file_name_rp0i[] = {
        "/home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_RUN_1.root"
    };

    Float_t factor_bad_pixels_rp0i = 10;
    Float_t factor_bad_pixels_rp1i = 10;
    Double_t integral_rp0i[nFiles] = {}, integral_err_rp0i[nFiles] = {};
    Double_t integral_rp1i[nFiles] = {}, integral_err_rp1i[nFiles] = {};
    Double_t integral_norm[nFiles] = {}, integral_err_norm[nFiles] = {};

    TH2D* h_rp0i_clone[nFiles];
    TH2D* h_rp1i_clone[nFiles];

    TCanvas* c1_rp0i = new TCanvas("c1_rp0i","c1_rp0i",500,500);
    TCanvas* c1_rp1i = new TCanvas("c1_rp1i","c1_rp1i",500,500);

    for(Int_t ifl = 0; ifl < nFiles; ifl++)
    {
        cout<<endl;
        cout<<"--> Input file RP0I: "<<file_name_rp0i[ifl]<<endl;
        cout<<"--> Input file RP1I: "<<file_name_rp1i[ifl]<<endl;

        TFile *_file_rp0i = TFile::Open(file_name_rp0i[ifl].Data());
        TFile *_file_rp1i = TFile::Open(file_name_rp1i[ifl].Data());

        TH2D* h_rp0i = (TH2D*)_file_rp0i->Get("h_1");
        TH2D* h_rp1i = (TH2D*)_file_rp1i->Get("h_1");

        TString hist_name_rp0i = "h_rp0i_clone_"; hist_name_rp0i += ifl;
        TString hist_name_rp1i = "h_rp1i_clone_"; hist_name_rp1i += ifl;

        h_rp0i_clone[ifl] = new TH2D(hist_name_rp0i.Data(),hist_name_rp0i.Data(),N_PIXELS/2,1,N_PIXELS/2,N_PIXELS/2,1,N_PIXELS/2);
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

                if(h_rp1i->GetBinContent(xi,yi) < h_rp1i->GetBinContent(xi,yi+1) + factor_bad_pixels_rp1i*h_rp1i->GetBinError(xi,yi+1) &&
                        h_rp1i->GetBinContent(xi,yi) > h_rp1i->GetBinContent(xi,yi+1) - factor_bad_pixels_rp1i*h_rp1i->GetBinError(xi,yi+1))
                {
                    h_rp1i_clone[ifl]->SetBinContent(xi,yi,h_rp1i->GetBinContent(xi,yi));
                    h_rp1i_clone[ifl]->SetBinError(xi,yi,h_rp1i->GetBinError(xi,yi));
                }
            }
        }

        integral_rp0i[ifl] = h_rp0i_clone[ifl]->IntegralAndError(1,N_PIXELS/2,1,N_PIXELS/2,integral_err_rp0i[ifl]);
        integral_rp1i[ifl] = h_rp1i_clone[ifl]->IntegralAndError(1,N_PIXELS/2,1,N_PIXELS/2,integral_err_rp1i[ifl]);
        integral_norm[ifl] = integral_rp1i[ifl]/integral_rp0i[ifl];
        integral_err_norm[ifl] = TMath::Sqrt(TMath::Power(integral_err_rp1i[ifl]/integral_rp0i[ifl],2)
                                             + TMath::Power(integral_rp1i[ifl]*integral_err_rp0i[ifl]/(integral_rp0i[ifl]*integral_rp0i[ifl]),2));


        c1_rp0i->cd();
//        h_rp0i_clone[ifl]->SetMaximum(5000);
        h_rp0i_clone[ifl]->Draw("colz");

        c1_rp1i->cd();
//        h_rp1i_clone[ifl]->SetMaximum(500);
        h_rp1i_clone[ifl]->Draw("colz");
    }

    cout<<endl<<endl;
    cout<<"--> integral_rp0i[0] = "<<integral_rp0i[0]<<" +/- "<<integral_err_rp0i[0]<<endl;
    cout<<"--> integral_rp1i[0] = "<<integral_rp1i[0]<<" +/- "<<integral_err_rp1i[0]<<endl;
    cout<<"--> integral_norm[0] = "<<integral_norm[0]<<" +/- "<<integral_err_norm[0]<<endl;

    return 0;
}
