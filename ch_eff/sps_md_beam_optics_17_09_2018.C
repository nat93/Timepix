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

Double_t get_beam_position(Double_t b1, Double_t b2, Double_t x1, Double_t m1, Double_t m2, Double_t theta);

int sps_md_beam_optics_17_09_2018()
{
    cout<<"--> function_1() -- to convert DAT to ROOT file"<<endl;
    cout<<"--> function_2() -- "<<endl;
    return 0;
}

void function_1()
{
    TString fileName = "twiss_ua9_2018.dat";
    ifstream inputfile(fileName.Data());
    string line, word;

    int nElements = 0;
    string _keyword, _name;
    double _S, _L, _X, _Y, _BETX, _BETY, _ALFX, _ALFY, _MUX, _MUY, _DX, _DY, _DPX, _DPY;

    fileName += ".root";
    TFile *_file = new TFile(fileName.Data(),"RECREATE");
    if ( _file->IsOpen() ) cout<<endl<<"--> Output file '"<<_file->GetName()<<"' opened successfully!"<<endl;
    TTree *tree = new TTree("tree","Tree with vectors");
    tree->Branch("_keyword", &_keyword);
    tree->Branch("_name",    &_name);
    tree->Branch("_S",       &_S);
    tree->Branch("_L",       &_L);
    tree->Branch("_X",       &_X);
    tree->Branch("_Y",       &_Y);
    tree->Branch("_BETX",    &_BETX);
    tree->Branch("_BETY",    &_BETY);
    tree->Branch("_ALFX",    &_ALFX);
    tree->Branch("_ALFY",    &_ALFY);
    tree->Branch("_MUX",     &_MUX);
    tree->Branch("_MUY",     &_MUY);
    tree->Branch("_DX",      &_DX);
    tree->Branch("_DY",      &_DY);
    tree->Branch("_DPX",     &_DPX);
    tree->Branch("_DPY",     &_DPY);

    if(inputfile.is_open())
    {
        while(getline(inputfile, line))
        {
            if(word == "KEYWORD")
            {
                istringstream iss(line);
                iss>>_keyword>>_name>>_S>>_L>>_X>>_Y>>_BETX>>_BETY>>_ALFX>>_ALFY>>_MUX>>_MUY>>_DX>>_DY>>_DPX>>_DPY;
                tree->Fill();
                nElements++;

                printf("\r--> nElements: %7.d",nElements);
                fflush(stdout);
            }
            else
            {
                istringstream iss(line);
                iss>>word;
                if(word == "KEYWORD")
                {
                    cout<<line<<endl;
                }
            }
        }
        inputfile.close();
    }
    else
    {
        cout<<"--> ERROR: Unable to open the input file!"<<endl;
        assert(0);
    }
    cout<<endl;
    cout<<"--> nEntries = "<<tree->GetEntriesFast()<<endl;

    _file->Write();
    delete _file;
}

void function_2()
{
    TString inFileName = "twiss_ua9_2018.dat.root";

    TChain* fChain = new TChain("tree");
    fChain->Add(inFileName.Data());

    string* _keyword = new string();
    string* _name = new string();
    double _S, _L, _X, _Y, _BETX, _BETY, _ALFX, _ALFY, _MUX, _MUY, _DX, _DY, _DPX, _DPY;
    fChain->SetBranchAddress("_keyword", &_keyword);
    fChain->SetBranchAddress("_name",    &_name);
    fChain->SetBranchAddress("_S",       &_S);
    fChain->SetBranchAddress("_L",       &_L);
    fChain->SetBranchAddress("_X",       &_X);
    fChain->SetBranchAddress("_Y",       &_Y);
    fChain->SetBranchAddress("_BETX",    &_BETX);
    fChain->SetBranchAddress("_BETY",    &_BETY);
    fChain->SetBranchAddress("_ALFX",    &_ALFX);
    fChain->SetBranchAddress("_ALFY",    &_ALFY);
    fChain->SetBranchAddress("_MUX",     &_MUX);
    fChain->SetBranchAddress("_MUY",     &_MUY);
    fChain->SetBranchAddress("_DX",      &_DX);
    fChain->SetBranchAddress("_DY",      &_DY);
    fChain->SetBranchAddress("_DPX",     &_DPX);
    fChain->SetBranchAddress("_DPY",     &_DPY);

    Long64_t nEntries = (fChain->GetEntries());
    cout<<"--> nEntries = "<<nEntries<<endl;

    vector<Double_t> S;
    vector<Double_t> BETX;
    vector<Double_t> BETY;
    vector<Double_t> SIGMAX;
    vector<Double_t> SIGMAY;
    vector<Double_t> BEAMX_1SIGMA_TOP;
    vector<Double_t> BEAMX_1SIGMA_BOT;
    vector<Double_t> BEAMY_1SIGMA_TOP;
    vector<Double_t> BEAMY_1SIGMA_BOT;
    vector<Double_t> BEAMX_4SIGMA_TOP;
    vector<Double_t> BEAMX_4SIGMA_BOT;
    vector<Double_t> BEAMY_4SIGMA_TOP;
    vector<Double_t> BEAMY_4SIGMA_BOT;
    vector<Double_t> BEAMX_C4SIGMA_TOP;
    vector<Double_t> BEAMX_C4SIGMA_BOT;

    vector<Double_t> BEAMX_DEFLECTED_CRY2;
    vector<Double_t> BEAMX_DEFLECTED_CRY2_MCRA; // Minus Critical Angle
    vector<Double_t> BEAMX_DEFLECTED_CRY2_PCRA; // Plus Critical Angle
    vector<Double_t> S_DEFLECTED_CRY2;
    vector<Double_t> BEAMX_DEFLECTED_CRY3;
    vector<Double_t> BEAMX_DEFLECTED_CRY3_MCRA; // Minus Critical Angle
    vector<Double_t> BEAMX_DEFLECTED_CRY3_PCRA; // Plus Critical Angle
    vector<Double_t> S_DEFLECTED_CRY3;

    Double_t S_CRYSTAL2, S_CRYSTAL3, S_LHC_COLL, S_RP1_INT, S_RP0_INT;
    Double_t BEAMX_CRYSTAL2, MUX_CRYSTAL2, BETAX_CRYSTAL2, BEAMX_CRYSTAL3, BEAMX_CRYSTAL3_MCRA, BEAMX_CRYSTAL3_PCRA, MUX_CRYSTAL3, BETAX_CRYSTAL3;
    Int_t I_CRYSTAL2, I_CRYSTAL3, I_LHC_COLL, I_RP1_INT, I_RP0_INT;
    Int_t ITERATOR = 0;

    const Double_t EMITTANCE = 3.4e-9;
    const Double_t S_MIN = 5120;
    const Double_t S_MAX = 5320;
    const Double_t BEAMX_MIN = -0.040*1e3;
    const Double_t BEAMX_MAX =  0.010*1e3;
    const Double_t CRYSTAL2_DEFL = -300.0e-6;
    const Double_t CRYSTAL3_DEFL = -200.0e-6;
    const Double_t CRYSTAL_CRITICAL_ANGLE = 10.9e-6;

    Bool_t CRYSTAL2_CH_STATUS = false;
    Bool_t CRYSTAL3_CH_STATUS = false;

    Double_t sigmax_coef = 1.7e-3/(4.0*TMath::Sqrt(EMITTANCE*35.9873));
    cout<<"--> sigmax_coef = "<<sigmax_coef<<endl;

    for(Long64_t i = 0; i < nEntries; i++)
    {
        if(i%10 == 0)
        {
            printf("\r--> Working: %3.1f %%",100*(Double_t)i/nEntries);
            fflush(stdout);
        }

        fChain->GetEntry(i);

        if(_S < S_MIN || _S > S_MAX) continue;

        string name = _name->c_str();
        if(name.find("CRY2.51652.UA9") == 1)
        {
            S_CRYSTAL2 = _S;
            I_CRYSTAL2 = ITERATOR;
            MUX_CRYSTAL2 = _MUX;
            BETAX_CRYSTAL2 = _BETX;
            BEAMX_CRYSTAL2 = _X-sigmax_coef*4.0*TMath::Sqrt(EMITTANCE*_BETX);
            CRYSTAL2_CH_STATUS = true;
        }
        if(name.find("CRY3.51799.UA9") == 1)
        {
            S_CRYSTAL3 = _S;
            I_CRYSTAL3 = ITERATOR;
            MUX_CRYSTAL3 = _MUX;
            BETAX_CRYSTAL3 = _BETX;
            BEAMX_CRYSTAL3 = get_beam_position(BETAX_CRYSTAL2,BETAX_CRYSTAL3,BEAMX_CRYSTAL2,MUX_CRYSTAL2,MUX_CRYSTAL3,CRYSTAL2_DEFL);
            BEAMX_CRYSTAL3_MCRA = get_beam_position(BETAX_CRYSTAL2,BETAX_CRYSTAL3,BEAMX_CRYSTAL2,MUX_CRYSTAL2,MUX_CRYSTAL3,CRYSTAL2_DEFL-CRYSTAL_CRITICAL_ANGLE);
            BEAMX_CRYSTAL3_PCRA = get_beam_position(BETAX_CRYSTAL2,BETAX_CRYSTAL3,BEAMX_CRYSTAL2,MUX_CRYSTAL2,MUX_CRYSTAL3,CRYSTAL2_DEFL+CRYSTAL_CRITICAL_ANGLE);

            cout<<endl<<"--> Beam size on the CRYSTAL3 <--"<<endl;
            cout<<BEAMX_CRYSTAL3_PCRA-BEAMX_CRYSTAL3_MCRA<<" [m]"<<endl;
            cout<<(BEAMX_CRYSTAL3_PCRA-BEAMX_CRYSTAL3_MCRA)*_ALFX/_BETX<<" [rad]"<<endl;

            CRYSTAL3_CH_STATUS = true;
        }
        if(name.find("TCSMS.51934.UA9") == 1)
        {
            S_LHC_COLL = _S;
            I_LHC_COLL = ITERATOR;
        }
        if(name.find("XRPH0.51779.UA9") == 1)
        {
            S_RP0_INT = _S;
            I_RP0_INT = ITERATOR;
        }
        if(name.find("XRPH.51937.UA9") == 1)
        {
            S_RP1_INT = _S;
            I_RP1_INT = ITERATOR;
        }

        S.push_back(_S);
        BETX.push_back(_BETX);
        BETY.push_back(_BETY);
        SIGMAX.push_back(TMath::Sqrt(EMITTANCE*_BETX)*1e3);
        SIGMAY.push_back(TMath::Sqrt(EMITTANCE*_BETY)*1e3);
        BEAMX_1SIGMA_TOP.push_back(_X+1.0*TMath::Sqrt(EMITTANCE*_BETX)*1e3);
        BEAMX_1SIGMA_BOT.push_back(_X-1.0*TMath::Sqrt(EMITTANCE*_BETX)*1e3);
        BEAMY_1SIGMA_TOP.push_back(_Y+1.0*TMath::Sqrt(EMITTANCE*_BETY)*1e3);
        BEAMY_1SIGMA_BOT.push_back(_Y-1.0*TMath::Sqrt(EMITTANCE*_BETY)*1e3);
        BEAMX_4SIGMA_TOP.push_back(_X+4.0*TMath::Sqrt(EMITTANCE*_BETX)*1e3);
        BEAMX_4SIGMA_BOT.push_back(_X-4.0*TMath::Sqrt(EMITTANCE*_BETX)*1e3);
        BEAMY_4SIGMA_TOP.push_back(_Y+4.0*TMath::Sqrt(EMITTANCE*_BETY)*1e3);
        BEAMY_4SIGMA_BOT.push_back(_Y-4.0*TMath::Sqrt(EMITTANCE*_BETY)*1e3);
        BEAMX_C4SIGMA_TOP.push_back(_X+sigmax_coef*4.0*TMath::Sqrt(EMITTANCE*_BETX)*1e3);
        BEAMX_C4SIGMA_BOT.push_back(_X-sigmax_coef*4.0*TMath::Sqrt(EMITTANCE*_BETX)*1e3);

        if(CRYSTAL2_CH_STATUS)
        {
            S_DEFLECTED_CRY2.push_back(_S);
            BEAMX_DEFLECTED_CRY2.push_back(get_beam_position(BETAX_CRYSTAL2,_BETX,BEAMX_CRYSTAL2,MUX_CRYSTAL2,_MUX,CRYSTAL2_DEFL)*1e3);
            BEAMX_DEFLECTED_CRY2_MCRA.push_back(get_beam_position(BETAX_CRYSTAL2,_BETX,BEAMX_CRYSTAL2,MUX_CRYSTAL2,_MUX,CRYSTAL2_DEFL-CRYSTAL_CRITICAL_ANGLE)*1e3);
            BEAMX_DEFLECTED_CRY2_PCRA.push_back(get_beam_position(BETAX_CRYSTAL2,_BETX,BEAMX_CRYSTAL2,MUX_CRYSTAL2,_MUX,CRYSTAL2_DEFL+CRYSTAL_CRITICAL_ANGLE)*1e3);
        }
        if(CRYSTAL3_CH_STATUS)
        {
            S_DEFLECTED_CRY3.push_back(_S);
            BEAMX_DEFLECTED_CRY3.push_back(get_beam_position(BETAX_CRYSTAL3,_BETX,BEAMX_CRYSTAL3,MUX_CRYSTAL3,_MUX,CRYSTAL3_DEFL)*1e3);
            BEAMX_DEFLECTED_CRY3_MCRA.push_back(get_beam_position(BETAX_CRYSTAL3,_BETX,BEAMX_CRYSTAL3_MCRA,MUX_CRYSTAL3,_MUX,CRYSTAL3_DEFL-CRYSTAL_CRITICAL_ANGLE)*1e3);
            BEAMX_DEFLECTED_CRY3_PCRA.push_back(get_beam_position(BETAX_CRYSTAL3,_BETX,BEAMX_CRYSTAL3_PCRA,MUX_CRYSTAL3,_MUX,CRYSTAL3_DEFL+CRYSTAL_CRITICAL_ANGLE)*1e3);
        }

        ITERATOR++;
    }
    cout<<endl;

    Int_t dim = S.size();

    // BETA
    TCanvas* c_0 = new TCanvas("c_0","c_0",1800,900);
    c_0->cd();
    gPad->SetGrid();
    TGraph *gr_betax = new TGraph(dim, &S[0], &BETX[0]);
    TGraph *gr_betay = new TGraph(dim, &S[0], &BETY[0]);
    TMultiGraph* mg_beta = new TMultiGraph();
    mg_beta->SetTitle("#beta_{x}(s) , #beta_{y}(s)");

    gr_betax->SetLineWidth(2);
    gr_betax->SetLineColor(kBlue);
    gr_betax->SetName("#beta_{x}");
    gr_betax->SetTitle("#beta_{x}");
    gr_betay->SetLineWidth(2);
    gr_betay->SetLineColor(kRed);
    gr_betay->SetName("#beta_{y}");
    gr_betay->SetTitle("#beta_{y}");

    mg_beta->Add(gr_betax);
    mg_beta->Add(gr_betay);

    mg_beta->Draw("AL");

    mg_beta->GetXaxis()->SetTitle("s [m]");
    mg_beta->GetYaxis()->SetTitle("#beta [m]");
    mg_beta->GetYaxis()->CenterTitle();
    mg_beta->GetXaxis()->CenterTitle();
    c_0->BuildLegend();

    // SIGMA
    TCanvas* c_1 = new TCanvas("c_1","c_1",1800,900);
    c_1->cd();
    gPad->SetGrid();
    TGraph *gr_sigmax = new TGraph(dim, &S[0], &SIGMAX[0]);
    TGraph *gr_sigmay = new TGraph(dim, &S[0], &SIGMAY[0]);
    TMultiGraph* mg_sigma = new TMultiGraph();
    mg_sigma->SetTitle("#sigma_{x}(s) , #sigma_{y}(s)");

    gr_sigmax->SetLineWidth(2);
    gr_sigmax->SetLineColor(kBlue);
    gr_sigmax->SetName("#sigma_{x}");
    gr_sigmax->SetTitle("#sigma_{x}");
    gr_sigmay->SetLineWidth(2);
    gr_sigmay->SetLineColor(kRed);
    gr_sigmay->SetName("#sigma_{y}");
    gr_sigmay->SetTitle("#sigma_{y}");

    mg_sigma->Add(gr_sigmax);
    mg_sigma->Add(gr_sigmay);

    mg_sigma->Draw("AL");

    mg_sigma->GetXaxis()->SetTitle("s [m]");
    mg_sigma->GetYaxis()->SetTitle("#sigma [mm]");
    mg_sigma->GetYaxis()->CenterTitle();
    mg_sigma->GetXaxis()->CenterTitle();
    c_1->BuildLegend();

    // BEAM X Experiment
    TCanvas* c_2 = new TCanvas("c_2","c_2",1800,900);
    c_2->cd();
    gPad->SetGrid();
    TGraph *gr_beamx_1s_top = new TGraph(dim, &S[0], &BEAMX_1SIGMA_TOP[0]);
    TGraph *gr_beamx_1s_bot = new TGraph(dim, &S[0], &BEAMX_1SIGMA_BOT[0]);
    TGraph *gr_beamx_4s_top = new TGraph(dim, &S[0], &BEAMX_4SIGMA_TOP[0]);
    TGraph *gr_beamx_4s_bot = new TGraph(dim, &S[0], &BEAMX_4SIGMA_BOT[0]);
    TGraph *gr_beamx_c4s_top = new TGraph(dim, &S[0], &BEAMX_C4SIGMA_TOP[0]);
    TGraph *gr_beamx_c4s_bot = new TGraph(dim, &S[0], &BEAMX_C4SIGMA_BOT[0]);

    Int_t dim_defl_cry2 = S_DEFLECTED_CRY2.size();
    TGraph *gr_beamx_defl_cry2 = new TGraph(dim_defl_cry2, &S_DEFLECTED_CRY2[0], &BEAMX_DEFLECTED_CRY2[0]);
    TGraph *gr_beamx_defl_cry2_mcra = new TGraph(dim_defl_cry2, &S_DEFLECTED_CRY2[0], &BEAMX_DEFLECTED_CRY2_MCRA[0]);
    TGraph *gr_beamx_defl_cry2_pcra = new TGraph(dim_defl_cry2, &S_DEFLECTED_CRY2[0], &BEAMX_DEFLECTED_CRY2_PCRA[0]);

    Int_t dim_defl_cry3 = S_DEFLECTED_CRY3.size();
    TGraph *gr_beamx_defl_cry3 = new TGraph(dim_defl_cry3, &S_DEFLECTED_CRY3[0], &BEAMX_DEFLECTED_CRY3[0]);
    TGraph *gr_beamx_defl_cry3_mcra = new TGraph(dim_defl_cry3, &S_DEFLECTED_CRY3[0], &BEAMX_DEFLECTED_CRY3_MCRA[0]);
    TGraph *gr_beamx_defl_cry3_pcra = new TGraph(dim_defl_cry3, &S_DEFLECTED_CRY3[0], &BEAMX_DEFLECTED_CRY3_PCRA[0]);

    TMultiGraph* mg_beamx = new TMultiGraph();
    mg_beamx->SetTitle("beam_{x}");

    gr_beamx_1s_top->SetLineWidth(2);
    gr_beamx_1s_top->SetLineColor(kBlue+2);
    gr_beamx_1s_top->SetName("gr_beamx_1s_top");
    gr_beamx_1s_top->SetTitle("1#sigma");
    gr_beamx_1s_bot->SetLineWidth(2);
    gr_beamx_1s_bot->SetLineColor(kBlue+2);
    gr_beamx_1s_bot->SetName("gr_beamx_1s_bot");
    gr_beamx_4s_top->SetLineWidth(2);
    gr_beamx_4s_top->SetLineColor(kBlue+1);
    gr_beamx_4s_top->SetName("gr_beamx_4s_top");
    gr_beamx_4s_top->SetTitle("4#sigma");
    gr_beamx_4s_bot->SetLineWidth(2);
    gr_beamx_4s_bot->SetLineColor(kBlue+1);
    gr_beamx_4s_bot->SetName("gr_beamx_4s_bot");
    gr_beamx_c4s_top->SetLineWidth(2);
    gr_beamx_c4s_top->SetLineColor(kBlue);
    gr_beamx_c4s_top->SetName("gr_beamx_c4s_top");
    gr_beamx_c4s_top->SetTitle("c4#sigma");
    gr_beamx_c4s_bot->SetLineWidth(2);
    gr_beamx_c4s_bot->SetLineColor(kBlue);
    gr_beamx_c4s_bot->SetName("gr_beamx_c4s_bot");
    gr_beamx_defl_cry2->SetTitle("DEFL. CRYSTAL2");
    gr_beamx_defl_cry2->SetLineWidth(2);
    gr_beamx_defl_cry2->SetLineColor(kGreen);
    gr_beamx_defl_cry2->SetName("gr_beamx_defl_cry2");
    gr_beamx_defl_cry2_mcra->SetTitle("DEFL. CRYSTAL2 -#Theta_{c}");
    gr_beamx_defl_cry2_mcra->SetLineWidth(2);
    gr_beamx_defl_cry2_mcra->SetLineColor(kGreen);
    gr_beamx_defl_cry2_mcra->SetName("gr_beamx_defl_cry2_mcra");
    gr_beamx_defl_cry2_pcra->SetTitle("DEFL. CRYSTAL2 +#Theta_{c}");
    gr_beamx_defl_cry2_pcra->SetLineWidth(2);
    gr_beamx_defl_cry2_pcra->SetLineColor(kGreen);
    gr_beamx_defl_cry2_pcra->SetName("gr_beamx_defl_cry2_pcra");
    gr_beamx_defl_cry3->SetTitle("DEFL. CRYSTAL3");
    gr_beamx_defl_cry3->SetLineWidth(2);
    gr_beamx_defl_cry3->SetLineColor(kGreen+2);
    gr_beamx_defl_cry3->SetName("gr_beamx_defl_cry3");
    gr_beamx_defl_cry3_mcra->SetTitle("DEFL. CRYSTAL3 -#Theta_{c}");
    gr_beamx_defl_cry3_mcra->SetLineWidth(2);
    gr_beamx_defl_cry3_mcra->SetLineColor(kGreen+2);
    gr_beamx_defl_cry3_mcra->SetName("gr_beamx_defl_cry3_mcra");
    gr_beamx_defl_cry3_pcra->SetTitle("DEFL. CRYSTAL3 +#Theta_{c}");
    gr_beamx_defl_cry3_pcra->SetLineWidth(2);
    gr_beamx_defl_cry3_pcra->SetLineColor(kGreen+2);
    gr_beamx_defl_cry3_pcra->SetName("gr_beamx_defl_cry3_pcra");

//    enum EColor { kWhite =0,   kBlack =1,   kGray=920,
//                  kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
//                  kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };

    TLine* CRYSTAL2_linex = new TLine(S_CRYSTAL2,BEAMX_MIN,S_CRYSTAL2,BEAMX_C4SIGMA_BOT[I_CRYSTAL2]);
    CRYSTAL2_linex->SetLineColor(kRed);
    CRYSTAL2_linex->SetLineWidth(2);

    TLine* CRYSTAL3_linex = new TLine(S_CRYSTAL3,BEAMX_MIN,S_CRYSTAL3,BEAMX_C4SIGMA_BOT[I_CRYSTAL3]);
    CRYSTAL3_linex->SetLineColor(kMagenta);
    CRYSTAL3_linex->SetLineWidth(2);

    TLine* LHC_COLL_linex = new TLine(S_LHC_COLL,BEAMX_MIN,S_LHC_COLL,BEAMX_C4SIGMA_BOT[I_LHC_COLL]);
    LHC_COLL_linex->SetLineColor(kOrange);
    LHC_COLL_linex->SetLineWidth(2);

    TLine* RP0_INT_linex = new TLine(S_RP0_INT,BEAMX_MIN,S_RP0_INT,BEAMX_C4SIGMA_BOT[I_RP0_INT]);
    RP0_INT_linex->SetLineColor(kBlack);
    RP0_INT_linex->SetLineWidth(2);

    TLine* RP1_INT_linex = new TLine(S_RP1_INT,BEAMX_MIN,S_RP1_INT,BEAMX_C4SIGMA_BOT[I_RP1_INT]);
    RP1_INT_linex->SetLineColor(kGray+1);
    RP1_INT_linex->SetLineWidth(2);

    mg_beamx->Add(gr_beamx_1s_top);
    mg_beamx->Add(gr_beamx_1s_bot);
    mg_beamx->Add(gr_beamx_4s_top);
    mg_beamx->Add(gr_beamx_4s_bot);
    mg_beamx->Add(gr_beamx_c4s_top);
    mg_beamx->Add(gr_beamx_c4s_bot);
    mg_beamx->Add(gr_beamx_defl_cry2);
    mg_beamx->Add(gr_beamx_defl_cry2_mcra);
    mg_beamx->Add(gr_beamx_defl_cry2_pcra);
    mg_beamx->Add(gr_beamx_defl_cry3);
    mg_beamx->Add(gr_beamx_defl_cry3_mcra);
    mg_beamx->Add(gr_beamx_defl_cry3_pcra);

    mg_beamx->Draw("AL");
    CRYSTAL2_linex->Draw("same");
    CRYSTAL3_linex->Draw("same");
    LHC_COLL_linex->Draw("same");
    RP0_INT_linex->Draw("same");
    RP1_INT_linex->Draw("same");

    mg_beamx->GetXaxis()->SetTitle("s [m]");
    mg_beamx->GetYaxis()->SetTitle("x [mm]");
    mg_beamx->GetYaxis()->CenterTitle();
    mg_beamx->GetXaxis()->CenterTitle();
    mg_beamx->SetMaximum(BEAMX_MAX);
    mg_beamx->SetMinimum(BEAMX_MIN);

    fChain->Delete();
}

Double_t get_beam_position(Double_t b1, Double_t b2, Double_t x1, Double_t m1, Double_t m2, Double_t theta)
{
    return (TMath::Sqrt(b2/b1)*TMath::Cos((m2-m1)*TMath::TwoPi())*x1+
            theta*TMath::Sqrt(b2*b1)*TMath::Sin((m2-m1)*TMath::TwoPi()));
}
