//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     To convert multiple ROOT files to a single ROOT file
//                  for H8 and SPS.
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
#include "TCanvas.h"
#include "TStyle.h"
//C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const Int_t N_PIXELS    = 512;  // The maximum number of pixels per axis per chip
const Int_t N_MAX_CHIP  = 5;    // The maximum number of chips

void epochtime2date(time_t in_time, string &out_time);
bool plot_st = false;

int main(int argc, char *argv[])
{
    cout<<endl<<"--> CONVERT COMMON"<<endl<<endl;
    if(argc != 5)
    {
        cout<<"--> ERROR: Wrong number of input parameters("<<argc<<")!"<<endl;
        cout<<"--> [0] ./script_name"<<endl;
        cout<<"--> [1] Input file name (./Filename_)"<<endl;
        cout<<"--> [2] Initial index"<<endl;
        cout<<"--> [3] Final index"<<endl;
        cout<<"--> [4] Output file name"<<endl;
        assert(0);
    }

    TString fileName0 = argv[1];
    TString fileName1;
    TString fileName2 = argv[4];
    Int_t initial_index = atoi(argv[2]);
    Int_t final_index = atoi(argv[3]);

    TChain* fChain = new TChain("Tree");

    for(Int_t i = initial_index; i <= final_index; i++)
    {
        fileName1 = fileName0;
        fileName1 += i;
        fileName1 += ".root";
        cout<<"--> "<<fileName1<<endl;
        fChain->Add(fileName1.Data());
    }

    Double_t _UnixTime, _DeltaTHR, _Ikrum, _Bias, _Clock, _Gate, _Time_ms;
    UInt_t _AcquisType = -1, _TrigType = -1, _chip;
    Long64_t _event;
    Long64_t _COUNTS[N_PIXELS][N_PIXELS];

    fChain->SetBranchAddress("UnixTime",    &_UnixTime);
    fChain->SetBranchAddress("TrigType",    &_TrigType);
    fChain->SetBranchAddress("AcquisType",  &_AcquisType);
    fChain->SetBranchAddress("DeltaTHR",    &_DeltaTHR);
    fChain->SetBranchAddress("Ikrum",       &_Ikrum);
    fChain->SetBranchAddress("Bias",        &_Bias);
    fChain->SetBranchAddress("Clock",       &_Clock);
    fChain->SetBranchAddress("Gate",        &_Gate);

    fChain->SetBranchAddress("Time_ms",     &_Time_ms);
    fChain->SetBranchAddress("event",       &_event);
    fChain->SetBranchAddress("chip",        &_chip);
    fChain->SetBranchAddress("COUNTS",      _COUNTS);

    Long64_t nEntries = fChain->GetEntries();
    cout<<"--> Number of nEntries: "<<nEntries<<endl;

    TFile* outputfile = new TFile(fileName2.Data(),"RECREATE");
    if (outputfile->IsZombie())
    {
        cout<<"--> ERROR: The output file is zombie!"<<endl;
        assert(0);
    }

    TTree* tree[N_MAX_CHIP];

    Double_t Timems;

    TString Timems_cs   = "Timems[";
    Timems_cs   += N_MAX_CHIP;
    Timems_cs   += "]/D";

    TString COUNTS_cs   = "COUNTS[";
    COUNTS_cs   += N_PIXELS;
    COUNTS_cs   += "][";
    COUNTS_cs   += N_PIXELS;
    COUNTS_cs   += "]/L";

    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        TString tree_name = "Tree_";
        tree_name += i;
        tree[i] = new TTree(tree_name.Data(), "Medipix data");

        tree[i]->Branch("UnixTime",    &_UnixTime,      "UnixTime/D");
        tree[i]->Branch("TrigType",    &_TrigType,      "TrigType/i");
        tree[i]->Branch("AcquisType",  &_AcquisType,    "AcquisType/i");
        tree[i]->Branch("DeltaTHR",    &_DeltaTHR,      "DeltaTHR/D");
        tree[i]->Branch("Ikrum",       &_Ikrum,         "Ikrum/D");
        tree[i]->Branch("Bias",        &_Bias,          "Bias/D");
        tree[i]->Branch("Clock",       &_Clock,         "Clock/D");
        tree[i]->Branch("Gate",        &_Gate,          "Gate/D");
        tree[i]->Branch("Event",       &_event,         "Event/L");
        tree[i]->Branch("Time_ms",     &Timems,         "Time_ms/D");
        tree[i]->Branch("COUNTS",      _COUNTS,          COUNTS_cs.Data());
    }

    fChain->GetEntry(0);

    Double_t zero_time = _UnixTime*1000 + _Time_ms;
    Double_t last_time = _UnixTime*1000 + _Time_ms;

    for(Long64_t i = 0; i < nEntries; i++)
    {
        fChain->GetEntry(i);
        if(i%10 == 0)
        {
            printf("\r--> Working: %3.1f %%", 100*(Double_t)i/nEntries);
            fflush(stdout);
        }

        Timems = _UnixTime*1000 + _Time_ms;
        last_time = Timems;

        tree[_chip]->Fill();
    }

    cout<<endl;

    string zero_time_char;
    string last_time_char;
    epochtime2date(round(zero_time/1000),zero_time_char);
    epochtime2date(round(last_time/1000),last_time_char);

    cout<<"--> Initial time: "<<zero_time_char<<" | Final time: "<<last_time_char<<endl;

    for(Int_t i = 0; i < N_MAX_CHIP; i++)
    {
        tree[i]->Write();
        cout<<"--> The tree["<<i<<"] with "<<tree[i]->GetEntriesFast()<<" entries is written to the file: "<<outputfile->GetName()<<endl;
    }

    outputfile->Close();

    fChain->Delete();

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
    //(void) printf("%u -> '%s'\n", (unsigned) t, res);
}
