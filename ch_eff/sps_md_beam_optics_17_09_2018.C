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

int sps_md_beam_optics_17_09_2018()
{
    cout<<"--> function_1() -- to convert DAT to ROOT file"<<endl;
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
