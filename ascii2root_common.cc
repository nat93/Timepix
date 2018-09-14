//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     To convert ASCII file to ROOT file for H8 and SPS
//                  The format of the input data LOCAL_PC/PYTHON_SCRIPT type
//--------------------------------------------------------------------------//

//ROOT
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
//C++
#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

const Int_t N_PIXELS    = 512;  // The maximum number of pixels per axis per chip

int main(int argc, char *argv[])
{
    cout<<endl<<"--> ASCII2ROOT COMMON"<<endl<<endl;
    if(argc < 3)
    {
        cout<<"--> ERROR: Wrong number of input parameters("<<argc<<")!"<<endl;
        cout<<"--> [0] ./script_name"<<endl;
        cout<<"--> [1] Input file name"<<endl;
        cout<<"--> [2] Output file name"<<endl;
        assert(0);
    }

    TString fileName1 = argv[1];
    TString fileName2 = argv[2];

    cout<<"--> InputFileName: "<<fileName1<<endl;
    ifstream inputfile(fileName1.Data());
    string word;

    cout<<"--> OutputFileName: "<<fileName2<<endl;
    TFile* outputfile = new TFile(fileName2.Data(), "RECREATE");
    if (outputfile->IsZombie())
    {
        cout<<"--> ERROR: The output file is zombie!"<<endl;
        assert(0);
    }

    Double_t UnixTime = -1, DeltaTHR = -1, Ikrum = -1, Bias = -1, Clock = -1, Gate = -1;
    UInt_t AcquisType = 999, TrigType = 999, chip;

    Long64_t event;
    Double_t Time_ms;
    Long64_t COUNTS[N_PIXELS][N_PIXELS];

    TString COUNTS_cs   = "COUNTS[";
    COUNTS_cs   += N_PIXELS;
    COUNTS_cs   += "][";
    COUNTS_cs   += N_PIXELS;
    COUNTS_cs   += "]/L";

    TTree* tree = new TTree("Tree", "Medipix data");

    // HEADER
    tree->Branch("UnixTime",    &UnixTime,      "UnixTime/D");
    tree->Branch("TrigType",    &TrigType,      "TrigType/i");
    tree->Branch("AcquisType",  &AcquisType,    "AcquisType/i");
    tree->Branch("DeltaTHR",    &DeltaTHR,      "DeltaTHR/D");
    tree->Branch("Ikrum",       &Ikrum,         "Ikrum/D");
    tree->Branch("Bias",        &Bias,          "Bias/D");
    tree->Branch("Clock",       &Clock,         "Clock/D");
    tree->Branch("Gate",        &Gate,          "Gate/D");

    tree->Branch("Time_ms",     &Time_ms,      "Time_ms/D");
    tree->Branch("event",       &event,        "event/L");
    tree->Branch("chip",        &chip,         "chip/i");
    tree->Branch("COUNTS",      COUNTS,        COUNTS_cs.Data());

    if(inputfile.is_open())
    {
        // READ HEADER
        while(1)
        {
            inputfile>>word;
            if(word == "=Time")
            {
                inputfile>>word;
                inputfile>>UnixTime;
            }
            if(word == "=Trigger")
            {
                inputfile>>word;
                inputfile>>word;
                if(word == "HW")
                {
                    inputfile>>word;
                    if(word == "start") {TrigType = 3;}
                    if(word == "stop") {TrigType = 2;}
                }
                if(word == "SW") {TrigType = 1;}
            }
            if(word == "=Acquis.")
            {
                inputfile>>word;
                inputfile>>word;
                if(word == "TOA") {AcquisType = 3;}
                if(word == "TOT") {AcquisType = 1;}
                if(word == "MPX") {AcquisType = 0;}
            }
            if(word == "=DeltaTHR")
            {
                inputfile>>word;
                inputfile>>DeltaTHR;
            }
            if(word == "=Ikrum")
            {
                inputfile>>word;
                inputfile>>Ikrum;
            }
            if(word == "=Bias")
            {
                inputfile>>word;
                inputfile>>Bias;
            }
            if(word == "=Clock")
            {
                inputfile>>word;
                inputfile>>Clock;
            }
            if(word == "=Gate")
            {
                inputfile>>word;
                inputfile>>Gate;
            }
            if(word == "HEADER_END") {break;}
        }

        while(1)
        {
            // READ DATA FROM THE MATRIX
            inputfile>>Time_ms;
            inputfile>>event;
            inputfile>>chip;
            inputfile>>word; // 'START#'

            if(inputfile.eof()) {break;}
            //cout<<Time_ms<<"  "<<event<<"  "<<chip<<endl;

            for(Int_t xi = 0; xi < N_PIXELS; xi++)
            {
                for(Int_t yi = 0; yi < N_PIXELS; yi++)
                {
                    COUNTS[xi][yi] = 0;
                }
            }

            while(inputfile>>word)
            {
                if(word == "END#") {break;}
                Int_t _x, _y;

//                _x = atoi(word.c_str());
//                inputfile>>_y;

                //---- For the rotated Quadpix at H8 ----//
                _y = atoi(word.c_str());
                inputfile>>_x;
                _x = N_PIXELS-_x-1;
                //---------------------------------------//

                inputfile>>COUNTS[_x][_y];
            }

            tree->Fill();
        }
        inputfile.close();
    }
    else
    {
        cout<<"--> ERROR: Unable to open the input file!"<<endl;
        assert(0);
    }
    //tree->Print();
    tree->Write();
    outputfile->Close();

    return 0;
}
