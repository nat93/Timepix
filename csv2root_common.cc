//--------------------------------------------------------------------------//
// Author:          Andrii Natochii
// CERN-LAL-TSNUK:  2018
// Description:     To convert new CSV file to ROOT file for SPS
//                  The format of the new UA9Publisher
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

const Int_t N_PIXELS    = 512;  // The maximum number of pixels per axis per chip

UInt_t chipname2id(TString chipName);
UInt_t getChipIdFromTheNextLine(string nextLine);

int main(int argc, char *argv[])
{
    cout<<endl<<"--> CSV2ROOT COMMON"<<endl<<endl;
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
    UInt_t AcquisType = -1, TrigType = -1, chip, chip_new;

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

            if(word == "=Acq")
            {
                inputfile>>word;
                if(word == "time")
                {
                    inputfile>>word;
                    inputfile>>Gate;
                }
            }
            if(word == "=Timepix")
            {
                inputfile>>word;
                inputfile>>word;
                inputfile>>Clock;
            }
            if(word == "=Bias")
            {
                inputfile>>word;
                inputfile>>Bias;
            }
            if(word == "=Tpx")
            {
                inputfile>>word;
                if(word == "mode")
                {
                    inputfile>>word;
                    inputfile>>AcquisType;
                }
            }
            if(word == "ID,PosX,PosY,Value,")
            {
                inputfile>>word;
                istringstream csvStream(word);
                string csvElement;

                // Unix_time
                getline(csvStream, csvElement, ',');
                UnixTime = atof(csvElement.c_str());

                break;
            }
            if(inputfile.eof()) {break;}
        }

        string csvLine = word;

        while(1)
        {
            if(inputfile.eof()) {break;}
            // Frame time,Event number,Frame ID,Chip ID,PosX,PosY,Value,

            for(Int_t xi = 0; xi < N_PIXELS; xi++)
            {
                for(Int_t yi = 0; yi < N_PIXELS; yi++)
                {
                    COUNTS[xi][yi] = 0;
                }
            }

            while(1)
            {
                istringstream csvStream(csvLine);
                string csvElement;
                // Frame time
                getline(csvStream, csvElement, ',');
                Time_ms = (atof(csvElement.c_str()) - UnixTime)*1000.0;

                // Event number
                getline(csvStream, csvElement, ',');
                event = atof(csvElement.c_str());

                // Frame ID
                getline(csvStream, csvElement, ',');

                // Chip ID
                getline(csvStream, csvElement, ',');
                chip = chipname2id(csvElement);

                // PosX
                getline(csvStream, csvElement, ',');
                Int_t _x = atoi(csvElement.c_str());

                // PosY
                getline(csvStream, csvElement, ',');
                Int_t _y = atoi(csvElement.c_str());

                // Value
                getline(csvStream, csvElement, ',');
                COUNTS[_x][_y] = atoi(csvElement.c_str());

                inputfile>>csvLine;
                if(inputfile.eof()) {break;}

                // Check the next line for the same chip
                chip_new = getChipIdFromTheNextLine(csvLine);
                if(chip != chip_new){break;}
            }
//            cout<<"--> END:: "<<Time_ms<<" "<<event<<" "<<chip<<endl;

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

UInt_t chipname2id(TString chipName)
{
    UInt_t chipID = -999;

    if(chipName == "G02-W0108") chipID = 0;
    else if(chipName == "K09-W0255") chipID = 1;
    else if(chipName == "C08-W0255") chipID = 2;
    else if(chipName == "F04-W0108") chipID = 3;
    else if(chipName == "I02-W0108") chipID = 4;
    else chipID = -1;

    return chipID;
}

UInt_t getChipIdFromTheNextLine(string nextLine)
{
    istringstream csvStream(nextLine);
    string csvElement;
    getline(csvStream, csvElement, ',');
    getline(csvStream, csvElement, ',');
    getline(csvStream, csvElement, ',');
    getline(csvStream, csvElement, ',');
    return chipname2id(csvElement);
}

// From SPS:
//   devRP0E=0 #G02-W0108    FITpix 0384
//   devRP3I=1 #K09-W0255    FITpix 0393
//   devRP3E=2 #C08-W0255    FITpix 0399
//   devRP1I=3 #F04-W0108    FITpix 0409
//   devRP0I=4 #I02-W0108    FITpix 0415
