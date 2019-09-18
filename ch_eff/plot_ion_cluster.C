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
#include "TEllipse.h"
#include "TArrow.h"
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

int plot_ion_cluster()
{
    TFile *_file0 = TFile::Open("/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_RUN_1_outputfile_function_2_TEST.root");
//    TFile *_file0 = TFile::Open("/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_RUN_1_outputfile_function_2_TEST2.root");
//    TFile *_file0 = TFile::Open("/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_RUN_1_outputfile_function_2_TEST3.root");
//    TFile *_file0 = TFile::Open("/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_RUN_1_outputfile_function_2_TEST4.root");
//    TFile *_file0 = TFile::Open("/media/andrii/F492773C92770302/MedipixData/ROOT_FILES/SPS_2018_11_26_RUN_1_outputfile_function_2_TEST5.root");

    TChain* fChain = new TChain("Tree");
    fChain->Add(_file0->GetName());

    Double_t unix_time_clinfo, pos_x_clinfo, pos_x_err_clinfo, pos_y_clinfo, pos_y_err_clinfo, clocks_clinfo, size_x_clinfo, size_y_clinfo, azimut_clinfo;
    Int_t size_clinfo;
    Long64_t event_id_clinfo;

    fChain->SetBranchAddress("UnixTime",        &unix_time_clinfo);
    fChain->SetBranchAddress("ClusterClocks",   &clocks_clinfo);
    fChain->SetBranchAddress("ClusterSize",     &size_clinfo);
    fChain->SetBranchAddress("ClusterPosX",     &pos_x_clinfo);
    fChain->SetBranchAddress("ClusterPosXerr",  &pos_x_err_clinfo);
    fChain->SetBranchAddress("ClusterPosY",     &pos_y_clinfo);
    fChain->SetBranchAddress("ClusterPosYerr",  &pos_y_err_clinfo);
    fChain->SetBranchAddress("ClusterSizeX",    &size_x_clinfo);
    fChain->SetBranchAddress("ClusterSizeY",    &size_y_clinfo);
    fChain->SetBranchAddress("ClusterAzimut",   &azimut_clinfo);
    fChain->SetBranchAddress("EventID",         &event_id_clinfo);

    Int_t nEntries = fChain->GetEntries();
    cout<<"--> nEntries="<<nEntries<<endl;
    fChain->GetEntry(0);

    cout<<endl<<"event_id_clinfo="<<event_id_clinfo<<endl;
    cout<<"size_clinfo="<<size_clinfo<<" pos_x_clinfo="<<pos_x_clinfo<<" pos_y_clinfo="<<pos_y_clinfo<<endl;
    cout<<"size_x_clinfo="<<size_x_clinfo<<" size_y_clinfo="<<size_y_clinfo<<" azimut_clinfo="<<azimut_clinfo*180.0/TMath::Pi()<<endl;

    TH2D* h_1 = (TH2D*)_file0->Get("h_1");
    TH2D* h_com_4 = (TH2D*)_file0->Get("h_com_4");

//    gStyle->SetOptStat(0);

    TCanvas* cc0 = new TCanvas("cc0","cc0",1000,1000);
    cc0->cd();    
    gPad->SetGrid();
    h_1->GetXaxis()->SetRange(1,256);
    h_1->GetYaxis()->SetRange(1,256);
    h_1->Draw("colz");

    TCanvas* cc1 = new TCanvas("cc1","cc1",1000,1000);
    cc1->cd();
    gPad->SetGrid();
    h_com_4->GetXaxis()->SetRange(pos_x_clinfo-size_x_clinfo,pos_x_clinfo+size_x_clinfo);
    h_com_4->GetYaxis()->SetRange(pos_y_clinfo-size_y_clinfo,pos_y_clinfo+size_y_clinfo);
    h_com_4->Draw("colz");

    TLine* lineX = new TLine(pos_x_clinfo,pos_y_clinfo-size_y_clinfo/2,pos_x_clinfo,pos_y_clinfo+size_y_clinfo/2);
    TLine* lineY = new TLine(pos_x_clinfo-size_x_clinfo/2,pos_y_clinfo,pos_x_clinfo+size_x_clinfo/2,pos_y_clinfo);
    TLine* lineA = new TLine(pos_x_clinfo-size_x_clinfo/2,(-size_x_clinfo/2)*TMath::Tan(azimut_clinfo)+pos_y_clinfo,
                             pos_x_clinfo+size_x_clinfo/2,(+size_x_clinfo/2)*TMath::Tan(azimut_clinfo)+pos_y_clinfo);
    lineX->Draw();
    lineY->Draw();
    lineA->Draw();

    //_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-//

    TH2D* h1 = new TH2D("h1","radius vs theta",160,-8,8,500,0,500);
    TH1D* h2 = new TH1D("h2","theta",160,-8,8);
    TH1D* h3 = new TH1D("h3","radius",500,0,500);

    Double_t x, y, w;
    Double_t x0, y0, r, theta;
    Int_t locmax, locmay, locmaz;
    h_com_4->GetMaximumBin(locmax,locmay,locmaz);
    x0 = h_com_4->GetXaxis()->GetBinCenter(locmax) + h_com_4->GetXaxis()->GetBinWidth(locmax)/2;
    y0 = h_com_4->GetYaxis()->GetBinCenter(locmay) + h_com_4->GetYaxis()->GetBinWidth(locmay)/2;

    for(Int_t binxi = 1; binxi <= h_com_4->GetNbinsX(); binxi++)
    {
        for(Int_t binyi = 1; binyi <= h_com_4->GetNbinsY(); binyi++)
        {
            x = h_com_4->GetXaxis()->GetBinCenter(binxi);
            y = h_com_4->GetYaxis()->GetBinCenter(binyi);
            w = h_com_4->GetBinContent(binxi,binyi);

            getNewRadiusTheta(x0,y0,x,y,r,theta);

            h1->Fill(theta,r,w);
            h2->Fill(theta,w);
            h3->Fill(r,w);
        }
    }

    TCanvas* cc2 = new TCanvas("cc2","cc2",1000,1000);
    cc2->cd();
    h1->Draw("colz");

    TCanvas* cc3 = new TCanvas("cc3","cc3",1000,1000);
    cc3->cd();
    h2->Draw();

    TCanvas* cc4 = new TCanvas("cc4","cc4",1000,1000);
    cc4->cd();
    h3->Draw();

    TCanvas* cc5 = new TCanvas("cc5","cc5",1000,1000);
    cc5->cd();
    gPad->SetGrid();
    h_com_4->GetXaxis()->SetRange(pos_x_clinfo-size_x_clinfo,pos_x_clinfo+size_x_clinfo);
    h_com_4->GetYaxis()->SetRange(pos_y_clinfo-size_y_clinfo,pos_y_clinfo+size_y_clinfo);
    h_com_4->Draw("colz");

    Double_t rc = h3->GetMean();
    Double_t rm = h3->FindLastBinAbove(0);
    Double_t tc = h2->GetBinCenter(h2->GetMaximumBin());
    TLine* lineX2 = new TLine(x0,y0-rc,x0,y0+rc);
    TLine* lineY2 = new TLine(x0-rc,y0,x0+rc,y0);
    TEllipse* circle = new TEllipse(x0,y0,rc);
    circle->SetFillStyle(0);

    TArrow* lineA2 = new TArrow(x0,y0,rm*TMath::Cos(tc)+x0,rm*TMath::Sin(tc)+y0,0.05,"|>");
    lineA2->SetAngle(45);
    lineA2->SetFillColor(0);

    lineX2->Draw();
    lineY2->Draw();
    circle->Draw();
    lineA2->Draw();
    //_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-//
/**/
    return 0;
}

void getNewRadiusTheta(Double_t x0, Double_t y0, Double_t x, Double_t y, Double_t &rnew, Double_t &thetanew)
{
    Double_t r0 = TMath::Sqrt(x0*x0 + y0*y0);
    Double_t sinTheta0 = y0/r0;
    Double_t cosTheta0 = x0/r0;
    Double_t r = TMath::Sqrt(x*x + y*y);
    Double_t sinTheta = y/r;
    Double_t cosTheta = x/r;
    Double_t cosThetaTheta0 = sinTheta*sinTheta0 + cosTheta*cosTheta0;
    rnew = TMath::Sqrt(r*r + r0*r0 - 2.0*r*r0*cosThetaTheta0);

    if(x < x0)
        thetanew = TMath::Pi() - TMath::ASin((r*sinTheta - r0*sinTheta0)/rnew);
    else
        thetanew = TMath::ASin((r*sinTheta - r0*sinTheta0)/rnew);
}
