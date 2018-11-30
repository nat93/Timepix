//ROOT
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TH2Poly.h"
#include "TMultiGraph.h"
//C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int plot_tree()
{
    TChain* fChain = new TChain("Tree");
    fChain->Add("/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_L3_RP1I_RUN_1.root");

    cout<<fChain->GetEntries()<<endl;

    Double_t unix_time_clinfo, pos_x_clinfo, pos_x_err_clinfo, pos_y_clinfo, pos_y_err_clinfo, clocks_clinfo;
    Int_t size_clinfo;
    Long64_t event_id_clinfo;
    Double_t _Clock, _Gate;

    fChain->SetBranchAddress("UnixTime",        &unix_time_clinfo);
    fChain->SetBranchAddress("ClusterClocks",   &clocks_clinfo);
    fChain->SetBranchAddress("ClusterSize",     &size_clinfo);
    fChain->SetBranchAddress("ClusterPosX",     &pos_x_clinfo);
    fChain->SetBranchAddress("ClusterPosXerr",  &pos_x_err_clinfo);
    fChain->SetBranchAddress("ClusterPosY",     &pos_y_clinfo);
    fChain->SetBranchAddress("ClusterPosYerr",  &pos_y_err_clinfo);
    fChain->SetBranchAddress("EventID",         &event_id_clinfo);
    fChain->SetBranchAddress("Clock",           &_Clock);
    fChain->SetBranchAddress("Gate",            &_Gate);

    TH2D* h1 = new TH2D("h1","Volume vs Size;Cluster Size [pixels]; Cluster volume [ToT counts]",1000,0,1000,2e4,0,2e5);

    for(Int_t i = 0; i < fChain->GetEntries(); i++)
    {
      fChain->GetEntry(i);

      h1->Fill(size_clinfo,size_clinfo*clocks_clinfo);
    }
    h1->Draw("colz");

    return 0;
}
