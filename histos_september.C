// For H8 test beam with septum. Plot histogram difference for diffuser and septum
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

bool checkpixel(TH1D *h, Int_t bin);
void removenoizypixelsX(TH1D* h_1);
void removenoizypixelsY(TH1D* h_1);
void removenoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels);
int findcenter(TH1D* h_profile);
int getMinPositionWithinRangeTH1(Int_t min_x, Int_t max_x, TH1* histo);

int histos_september()
{
    cout<<"--> histos_0() -- to plot beam profiles and image for all runs"<<endl;
    cout<<"--> histos_1() -- to plot all alignment runs on the same plots"<<endl;
    cout<<"--> histos_2() -- to plot all runs with a septum runs and with a align run"<<endl;
    cout<<"--> histos_3() -- to calculate areas above and bellow the zero for all runs"<<endl;
    cout<<"--> histos_4() -- to calculate areas above and bellow the zero for the defined runs with some obstacles on the beam"<<endl;
    cout<<"--> histos_5() -- to calculate areas above left/right and bellow the zero for all runs"<<endl;
    cout<<"--> histos_6() -- to calculate areas above left/right and bellow the zero for the defined runs with some obstacles on the beam"<<endl;
    cout<<"--> histos_7() -- polar coordinate convertor"<<endl;
    cout<<"--> histos_8() -- to plot the difference between runs in 3D"<<endl;
    cout<<"--> histos_9() -- to calculate specified areas for the defined runs with some obstacles on the beam"<<endl;
    cout<<"--> histos_10() -- "<<endl;
    cout<<"--> histos_11() -- to plot the beam movement"<<endl;

    return 0;
}

int histos_0()//to plot beam profiles and image for all runs
{
    for(Int_t j = 1; j <= 26; j++)
    {
        if(j == 18 || j == 19) continue;
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_";
        filename += j;
        filename += ".root";
        TFile *_file0 = TFile::Open(filename.Data());

        filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_";
        filename += j;
        filename += ".root";
        TFile *_file1 = TFile::Open(filename.Data());

        TH2D* h_norm = (TH2D*)_file0->Get("h_1");
        removenoizypixelsXY(h_norm,10);

        TH2D* h_0 = (TH2D*)_file1->Get("h_1");
        removenoizypixelsXY(h_0,10);

        TH1D* h_1 = h_0->ProjectionX("h_2");
        TH1D* h_2 = h_0->ProjectionY("h_3");

        Double_t integral_norm, errintegral_norm;
        Double_t val_1, errval_1, val_2, errval_2;

        integral_norm = h_norm->IntegralAndError(1,h_norm->GetNbinsX(),1,h_norm->GetNbinsY(),errintegral_norm);

        for(Int_t i = 1; i <= h_1->GetNbinsX(); i++)
        {
            val_1 = (h_1->GetBinContent(i))/integral_norm;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_norm,2) +
                                   TMath::Power(h_1->GetBinContent(i)*errintegral_norm/(integral_norm*integral_norm),2));

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);
        }

        for(Int_t i = 1; i <= h_2->GetNbinsX(); i++)
        {
            val_2 = (h_2->GetBinContent(i))/integral_norm;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/integral_norm,2) +
                                   TMath::Power(h_2->GetBinContent(i)*errintegral_norm/(integral_norm*integral_norm),2));

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);
        }

        h_norm->Scale(1.0/integral_norm);
        h_norm->SetTitle("Upstream Timepix");
        h_norm->GetXaxis()->SetTitle("Horizontal Axis [pixels]");
        h_norm->GetYaxis()->SetTitle("Vertical Axis [pixels]");
//        h_norm->SetMaximum(60e-6);

        h_0->Scale(1.0/integral_norm);
        h_0->SetTitle("Downstream Timepix");
        h_0->GetXaxis()->SetTitle("Horizontal Axis [pixels]");
        h_0->GetYaxis()->SetTitle("Vertical Axis [pixels]");
//        h_0->SetMaximum(60e-6);

        h_1->SetLineWidth(4);
        h_1->SetTitle("Downstream Timepix (Projection on X axis)");
        h_1->GetXaxis()->SetTitle("Horizontal Axis [pixels]");
//        h_1->SetMaximum(10e-3);

        h_2->SetLineWidth(4);
        h_2->SetTitle("Downstream Timepix (Projection on Y axis)");
        h_2->GetXaxis()->SetTitle("Horizontal Axis [pixels]");
//        h_2->SetMaximum(10e-3);

        gStyle->SetOptStat(0);

        const Int_t NRGBs = 5;
        const Int_t NCont = 256;
        Double_t stops[NRGBs] = { 0.00, 0.30, 0.61, 0.84, 1.00 };
        Double_t red[NRGBs] = { 0.00, 0.00, 0.57, 0.90, 0.51 };
        Double_t green[NRGBs] = { 0.00, 0.65, 0.95, 0.20, 0.00 };
        Double_t blue[NRGBs] = { 0.51, 0.55, 0.15, 0.00, 0.10 };
        TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        gStyle->SetNumberContours(NCont);

        TString canvasname = "./canvas/run_";
        canvasname += j;
        TCanvas* c1 = new TCanvas(canvasname.Data(),canvasname.Data(),1000,1000);
        c1->Divide(2,2);

        c1->cd(3);
        h_norm->Draw("colz");
        c1->Update();
        TPaletteAxis *palette_norm = (TPaletteAxis*)h_norm->GetListOfFunctions()->FindObject("palette");
        palette_norm->SetX1NDC(0.905);
        palette_norm->SetX2NDC(0.925);
        palette_norm->SetY1NDC(0.1);
        palette_norm->SetY2NDC(0.9);
        c1->Modified();
        c1->Update();

        c1->cd(2);
        h_0->Draw("colz");
        c1->Update();
        TPaletteAxis *palette = (TPaletteAxis*)h_0->GetListOfFunctions()->FindObject("palette");
        palette->SetX1NDC(0.905);
        palette->SetX2NDC(0.925);
        palette->SetY1NDC(0.1);
        palette->SetY2NDC(0.9);
        c1->Modified();
        c1->Update();

        c1->cd(4);
        h_1->Draw();

        c1->cd(1);
        h_2->Draw();

        canvasname += ".png";
        c1->SaveAs(canvasname.Data());
    }

    return 0;
}

int histos_1()//to plot all alignment runs on the same plots
{
    TString canvasname_1 = "./canvas/prof_all_align_runs";
    TCanvas* c_1 = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
    c_1->cd();
    gStyle->SetOptStat(0);
    TLegend *legend = new TLegend(0.15,0.55,0.45,0.85,"","brNDC");

    for(Int_t j = 1; j <= 23; j++)
    {
        if(j == 1 || j == 17)
        {
            TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_";
            filename += j;
            filename += ".root";
            TFile *_file0 = TFile::Open(filename.Data());

            filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_";
            filename += j;
            filename += ".root";
            TFile *_file1 = TFile::Open(filename.Data());

            TH2D* h_norm = (TH2D*)_file0->Get("h_1");
            removenoizypixelsXY(h_norm,10);

            TH2D* h_0 = (TH2D*)_file1->Get("h_1");
            removenoizypixelsXY(h_0,10);

            TH1D* h_1 = h_0->ProjectionX("h_2");

            Double_t integral_norm, errintegral_norm;
            Double_t val_1, errval_1;

            integral_norm = h_norm->IntegralAndError(1,h_norm->GetNbinsX(),1,h_norm->GetNbinsY(),errintegral_norm);

            for(Int_t i = 1; i <= h_1->GetNbinsX(); i++)
            {
                val_1 = (h_1->GetBinContent(i))/integral_norm;
                errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_norm,2) +
                                       TMath::Power(h_1->GetBinContent(i)*errintegral_norm/(integral_norm*integral_norm),2));

                h_1->SetBinContent(i,val_1);
                h_1->SetBinError(i,errval_1);

            }
            h_1->SetLineWidth(4);
            h_1->SetLineColor(1+2*j);
//            h_1->SetMaximum(0.01);
            h_1->SetTitle("Downstream Timepix (Projection on X axis)");
            h_1->GetXaxis()->SetTitle("Horizontal [pixels]");

            h_1->Draw("same");
            TString name = "RUN_";
            name += j;
            legend->AddEntry(h_1,name.Data(),"ple");
        }
    }

    legend->Draw();

    canvasname_1 += ".png";
    c_1->SaveAs(canvasname_1.Data());

    return 0;
}

int histos_2()//to plot all runs with a septum runs and with a align run
{
    for(Int_t j = 1; j <= 23; j++)
    {
        if(j == 18 || j == 19) continue;
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_";
        filename += j;
        filename += ".root";
        TFile *_file_norm1 = TFile::Open(filename.Data());

        filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_";
        filename += j;
        filename += ".root";

        TFile *_file0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_3.root");
        TFile *_file_norm0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_3.root");
        TFile *_file1 = TFile::Open(filename.Data());
        TFile *_file2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_1.root");
        TFile *_file_norm2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_1.root");

        TH2D* h_norm_0 = (TH2D*)_file_norm0->Get("h_1");
        TH2D* h_norm_1 = (TH2D*)_file_norm1->Get("h_1");
        TH2D* h_norm_2 = (TH2D*)_file_norm2->Get("h_1");
        removenoizypixelsXY(h_norm_0,10);
        removenoizypixelsXY(h_norm_1,10);
        removenoizypixelsXY(h_norm_2,10);

        TH2D* h_0_0 = (TH2D*)_file0->Get("h_1");
        TH2D* h_0_1 = (TH2D*)_file1->Get("h_1");
        TH2D* h_0_2 = (TH2D*)_file2->Get("h_1");

        removenoizypixelsXY(h_0_0,10);
        removenoizypixelsXY(h_0_1,10);
        removenoizypixelsXY(h_0_2,10);

        TH1D* h_0 = h_0_0->ProjectionX("h_2_0");
        TH1D* h_1 = h_0_1->ProjectionX("h_2_1");
        TH1D* h_2 = h_0_2->ProjectionX("h_2_2");

        Double_t val_0, errval_0, integral_norm_0, errintegral_norm_0;
        Double_t val_1, errval_1, integral_norm_1, errintegral_norm_1;
        Double_t val_2, errval_2, integral_norm_2, errintegral_norm_2;

        integral_norm_0 = h_norm_0->IntegralAndError(1,h_norm_0->GetNbinsX(),1,h_norm_0->GetNbinsY(),errintegral_norm_0);
        integral_norm_1 = h_norm_1->IntegralAndError(1,h_norm_1->GetNbinsX(),1,h_norm_1->GetNbinsY(),errintegral_norm_1);
        integral_norm_2 = h_norm_2->IntegralAndError(1,h_norm_2->GetNbinsX(),1,h_norm_2->GetNbinsY(),errintegral_norm_2);

        for(Int_t i = 1; i <= h_0->GetNbinsX(); i++)
        {
            val_0 = (h_0->GetBinContent(i))/integral_norm_0;
            errval_0 = TMath::Sqrt(TMath::Power(h_0->GetBinError(i)/integral_norm_0,2) +
                                   TMath::Power(h_0->GetBinContent(i)*errintegral_norm_0/(integral_norm_0*integral_norm_0),2));

            val_1 = (h_1->GetBinContent(i))/integral_norm_1;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_norm_1,2) +
                                   TMath::Power(h_1->GetBinContent(i)*errintegral_norm_1/(integral_norm_1*integral_norm_1),2));

            val_2 = (h_2->GetBinContent(i))/integral_norm_2;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/integral_norm_2,2) +
                                   TMath::Power(h_2->GetBinContent(i)*errintegral_norm_2/(integral_norm_2*integral_norm_2),2));

            h_0->SetBinContent(i,val_0);
            h_0->SetBinError(i,errval_0);

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);
        }
        h_0->SetLineWidth(4);
        h_0->SetLineColor(kBlue);
        h_0->SetMaximum(10e-3);
        h_0->GetXaxis()->SetTitle("Horizontal [pixels]");
        h_0->SetTitle("Beam Profile along Vertical Axis");

        h_1->SetLineWidth(4);
        h_1->SetLineColor(kRed);
        h_1->SetMaximum(10e-3);
        h_1->GetXaxis()->SetTitle("Horizontal [pixels]");
        h_1->SetTitle("Beam Profile along Vertical Axis");

        h_2->SetLineWidth(4);
        h_2->SetLineColor(kBlack);
        h_2->SetMaximum(10e-3);
        h_2->GetXaxis()->SetTitle("Horizontal [pixels]");
        h_2->SetTitle("Beam Profile along Vertical Axis");

        gStyle->SetOptStat(0);

        TString canvasname_1 = "./canvas/prof_run3_run1_run";
        canvasname_1 += j;
        TCanvas* c_1 = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
        c_1->cd();
        h_1->Draw("hist");
        h_0->Draw("same & hist");
        h_2->Draw("same & hist");

        TLegend *legend = new TLegend(0.15,0.55,0.45,0.85,"","brNDC");
        TString name = "RUN_";
        name += j;
        legend->AddEntry(h_0,"RUN_3","ple");
        legend->AddEntry(h_1,name.Data(),"ple");
        legend->AddEntry(h_2,"RUN_1","ple");
        legend->Draw();

        canvasname_1 += ".png";
        c_1->SaveAs(canvasname_1.Data());
    }

    return 0;
}

int histos_3()// to calculate areas above and bellow the zero for all runs
{
    Int_t minCountBinX, maxCountBinX;
    Int_t rangeCountBinX = 50;

    // Diffuser
    const Int_t nn = 19;
    Int_t runid[] = {21, 16, 14, 26/*23*/, 12, 10, 8, 6, 25, 4, 24, 5, 7, 9, 11, 22, 13, 15, 20};
    Double_t position[] = {-2000, -1000, -500, -350, -200, -100, -50, -25, -10, 0, 10, 25, 50, 100, 200, 350, 500, 1000, 2000};
    TString canvasname_2 = "./canvas/diffuser";

    Double_t N1[nn] = {}, err_N1[nn] = {};
    Double_t N2[nn] = {}, err_N2[nn] = {};
    Double_t N12[nn] = {}, err_N12[nn] = {};
    // Difference
    Double_t maxY = 10.0e-4;
    Double_t minY = -10.0e-4;

    TCanvas* c_1[nn];

    for(Int_t j = 0; j < nn; j++)
    {
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_";
        filename += runid[j];
        filename += ".root";
        TFile *_file_norm1 = TFile::Open(filename.Data());

        filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_";
        filename += runid[j];
        filename += ".root";

        TFile *_file0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_3.root");
        TFile *_file_norm0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_3.root");
        TFile *_file1 = TFile::Open(filename.Data());
        TFile *_file2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_1.root");
        TFile *_file_norm2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_1.root");

        TH2D* h_norm_0 = (TH2D*)_file_norm0->Get("h_1");
        TH2D* h_norm_1 = (TH2D*)_file_norm1->Get("h_1");
        TH2D* h_norm_2 = (TH2D*)_file_norm2->Get("h_1");

        removenoizypixelsXY(h_norm_0,10);
        removenoizypixelsXY(h_norm_1,10);
        removenoizypixelsXY(h_norm_2,10);

        TH2D* h_0_0 = (TH2D*)_file0->Get("h_1");
        TH2D* h_0_1 = (TH2D*)_file1->Get("h_1");
        TH2D* h_0_2 = (TH2D*)_file2->Get("h_1");

        removenoizypixelsXY(h_0_0,10);
        removenoizypixelsXY(h_0_1,10);
        removenoizypixelsXY(h_0_2,10);

        TH1D* h_0 = h_0_0->ProjectionX("h_2_0");
        TH1D* h_1 = h_0_1->ProjectionX("h_2_1");
        TH1D* h_2 = h_0_2->ProjectionX("h_2_2");
        TH1D* h_01 = h_0_0->ProjectionX("h_2_01");

        Double_t max_h_2, errmax_h_2;

        Double_t val_0, errval_0, integral_norm_0, errintegral_norm_0;
        Double_t val_1, errval_1, integral_norm_1, errintegral_norm_1;
        Double_t val_2, errval_2;

        integral_norm_0 = h_norm_0->IntegralAndError(1,h_norm_0->GetNbinsX(),1,h_norm_0->GetNbinsY(),errintegral_norm_0);
        integral_norm_1 = h_norm_1->IntegralAndError(1,h_norm_1->GetNbinsX(),1,h_norm_1->GetNbinsY(),errintegral_norm_1);

        max_h_2 = h_2->GetBinContent(h_2->GetMaximumBin());
        errmax_h_2 = h_2->GetBinError(h_2->GetMaximumBin());

        for(Int_t i = 1; i <= h_0->GetNbinsX(); i++)
        {
            val_0 = (h_0->GetBinContent(i))/integral_norm_0;
            errval_0 = TMath::Sqrt(TMath::Power(h_0->GetBinError(i)/integral_norm_0,2) +
                                   TMath::Power(h_0->GetBinContent(i)*errintegral_norm_0/(integral_norm_0*integral_norm_0),2));

            val_1 = (h_1->GetBinContent(i))/integral_norm_1;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_norm_1,2) +
                                   TMath::Power(h_1->GetBinContent(i)*errintegral_norm_1/(integral_norm_1*integral_norm_1),2));

            val_2 = (h_2->GetBinContent(i))/max_h_2;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/max_h_2,2) + TMath::Power(h_2->GetBinContent(i)*errmax_h_2/(max_h_2*max_h_2),2));

            h_0->SetBinContent(i,val_0);
            h_0->SetBinError(i,errval_0);

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);

            // Difference
            if(val_2 != 0)
            {
//                val_2 = 1;
//                errval_2 = 0;
                h_01->SetBinContent(i,(val_1 - val_0)/val_2);
                h_01->SetBinError(i,TMath::Sqrt(TMath::Power(errval_0/val_2,2) +
                                                TMath::Power(errval_1/val_2,2) +
                                                TMath::Power((val_1 - val_0)*errval_2/(val_2*val_2),2)));
            }
            else
            {
	            h_01->SetBinContent(i,0);
	            h_01->SetBinError(i,0);
            }
        }

        h_01->SetLineWidth(2);
        h_01->SetLineColor(kBlue);
        h_01->SetMaximum(maxY);
        h_01->SetMinimum(minY);
        h_01->GetXaxis()->SetTitle("Horizontal [pixels]");

        gStyle->SetOptStat(0);

        if(j == 0)
        {
            TCanvas* c_3 = new TCanvas("Normalized_Alignment_Plot","Normalized_Alignment_Plot",1920,1080);
            c_3->cd();
            h_2->Draw();
            c_3->SaveAs("./canvas/Normalized_Alignment_Plot.png");
        }

        TString canvasname_1 = "./canvas/run3_run1_run";
        canvasname_1 += runid[j];
        c_1[j] = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
        c_1[j]->cd();
        h_01->SetTitle("Beam Profile along Vertical Axis");
        h_01->Draw();

        TLegend *legend = new TLegend(0.15,0.15,0.45,0.35,"","brNDC");

        // Difference
        TString name = "(RUN_";
        name += runid[j];
        name += " -- RUN_3)/RUN_1*";
//        name += " -- RUN_3)";

        gPad->SetGrid();
        legend->AddEntry(h_01,name.Data(),"ple");
        legend->Draw();

        canvasname_1 += ".png";

        //-----------------------------------------------------------------//
        minCountBinX = 350;
        maxCountBinX = 500;
        Int_t min_val_pos = getMinPositionWithinRangeTH1(minCountBinX,maxCountBinX,h_01);
        minCountBinX = min_val_pos - rangeCountBinX;
        maxCountBinX = min_val_pos + rangeCountBinX;
        if(maxCountBinX > h_01->GetNbinsX()) maxCountBinX = h_01->GetNbinsX();
        if(minCountBinX < 1) minCountBinX = 1;
        TLine *l_min=new TLine(minCountBinX,minY,minCountBinX,maxY);
        l_min->SetLineColor(kRed);
        l_min->SetLineWidth(3);
        l_min->SetLineStyle(7);
        l_min->Draw();
        TLine *l_mean=new TLine(min_val_pos,minY,min_val_pos,maxY);
        l_mean->SetLineColor(kBlack);
        l_mean->SetLineWidth(3);
        l_mean->SetLineStyle(7);
        l_mean->Draw();
        TLine *l_max=new TLine(maxCountBinX,minY,maxCountBinX,maxY);
        l_max->SetLineColor(kGreen+2);
        l_max->SetLineWidth(3);
        l_max->SetLineStyle(7);
        l_max->Draw();
        c_1[j]->SaveAs(canvasname_1.Data());

        for(Int_t i = minCountBinX; i <= maxCountBinX; i++)
        {
            if(h_01->GetBinContent(i) > 0)
            {
                N1[j] += h_01->GetBinContent(i);
                err_N1[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(h_01->GetBinContent(i) < 0)
            {
                N2[j] += TMath::Abs(h_01->GetBinContent(i));
                err_N2[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            N12[j] += h_01->GetBinContent(i);
            err_N12[j] += TMath::Power(h_01->GetBinError(i),2);
        }
        err_N1[j] = TMath::Sqrt(err_N1[j]);
        err_N2[j] = TMath::Sqrt(err_N2[j]);
        err_N12[j] = TMath::Sqrt(err_N12[j]);
        //-----------------------------------------------------------------//
    }

    TCanvas* c_2 = new TCanvas(canvasname_2.Data(),canvasname_2.Data(),1920,1080);
    TGraphErrors* gr1 = new TGraphErrors(nn,position,N1,0,err_N1);
    TGraphErrors* gr1_line = new TGraphErrors(nn,position,N1,0,err_N1);
    TGraphErrors* gr2 = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr2_line = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr12 = new TGraphErrors(nn,position,N12,0,err_N12);
    TGraphErrors* gr12_line = new TGraphErrors(nn,position,N12,0,err_N12);
    gr1->SetLineColor(kRed);
    gr1->SetMarkerColor(kRed);
    gr1_line->SetLineColor(kRed);
    gr2->SetLineColor(kBlue);
    gr2->SetMarkerColor(kBlue);
    gr2_line->SetLineColor(kBlue);
    gr12->SetLineColor(kBlack);
    gr12->SetMarkerColor(kBlack);
    gr12_line->SetLineColor(kBlack);
    gr1_line->SetLineStyle(7);
    gr2_line->SetLineStyle(7);
    gr12_line->SetLineStyle(7);
    gr1->SetLineWidth(4);
    gr1_line->SetLineWidth(4);
    gr2->SetLineWidth(4);
    gr2_line->SetLineWidth(4);
    gr12->SetLineWidth(4);
    gr12_line->SetLineWidth(4);
    gr1->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    gr12->SetMarkerStyle(22);
    gr1->SetMarkerSize(2);
    gr2->SetMarkerSize(2);
    gr12->SetMarkerSize(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gr1,"AP");
    mg->Add(gr1_line,"AL");
    mg->Add(gr2,"AP");
    mg->Add(gr2_line,"AL");
    mg->Add(gr12,"AP");
    mg->Add(gr12_line,"AL");
    mg->SetMaximum(0.04);
    mg->SetMinimum(-0.01);
    mg->Draw("A");

    TLegend *legend = new TLegend(0.80,0.80,0.95,0.95,"","brNDC");
    gPad->SetGrid();
    legend->AddEntry(gr1,"N1","ple");
    legend->AddEntry(gr2,"N2","ple");
    legend->AddEntry(gr12,"N1 -- N2","ple");
    legend->Draw();

    canvasname_2 += ".png";
    c_2->SaveAs(canvasname_2.Data());

    return 0;
}

int histos_4()// to calculate areas above and bellow the zero for the defined runs with some obstacles on the beam
{
    Int_t minCountBinX = 300;
    Int_t maxCountBinX = 500;

    const Int_t nn = 6;
    Int_t runid[] = {
        15,//only crystal am
        59,//only septum
        8,//septum + diffuser
        16,//septum + crystal am
        23,//septum + crystal ch
        31//septum + crystal vr
    };
    Double_t position[] = {0,1,2,3,4,5,6};
    TString canvasname_2 = "./canvas/comparison";

    Double_t N1[nn] = {}, err_N1[nn] = {};
    Double_t N2[nn] = {}, err_N2[nn] = {};
    Double_t N12[nn] = {}, err_N12[nn] = {};
    // Difference
    Double_t maxY = 20.0e-4;
    Double_t minY = -20.0e-4;

    for(Int_t j = 0; j < nn; j++)
    {
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_";
        filename += runid[j];
        filename += ".root";

        TFile *_file0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_28.root");
        TFile *_file1 = TFile::Open(filename.Data());
        TFile *_file2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_28.root");

        TH2D* h_0_0 = (TH2D*)_file0->Get("h_91");
        TH2D* h_0_1 = (TH2D*)_file1->Get("h_91");
        TH2D* h_0_2 = (TH2D*)_file2->Get("h_91");

        removenoizypixelsXY(h_0_0,10);
        removenoizypixelsXY(h_0_1,10);
        removenoizypixelsXY(h_0_2,10);

        TH1D* h_0 = h_0_0->ProjectionX("h_92_0");
        TH1D* h_1 = h_0_1->ProjectionX("h_92_1");
        TH1D* h_2 = h_0_2->ProjectionX("h_92_2");
        TH1D* h_01 = h_0_0->ProjectionX("h_92_01");

        Double_t val_0, errval_0, integral_0, errintegral_0;
        Double_t val_1, errval_1, integral_1, errintegral_1;
        Double_t val_2, errval_2, max_h_2, errmax_h_2;

        integral_0 = h_0->IntegralAndError(1,h_0->GetNbinsX(),errintegral_0);
        integral_1 = h_1->IntegralAndError(1,h_1->GetNbinsX(),errintegral_1);

        max_h_2 = h_2->GetBinContent(h_2->GetMaximumBin());
        errmax_h_2 = h_2->GetBinError(h_2->GetMaximumBin());

        for(Int_t i = 1; i <= h_0->GetNbinsX(); i++)
        {
            val_0 = (h_0->GetBinContent(i))/integral_0;
            errval_0 = TMath::Sqrt(TMath::Power(h_0->GetBinError(i)/integral_0,2) + TMath::Power(h_0->GetBinContent(i)*errintegral_0/(integral_0*integral_0),2));

            val_1 = (h_1->GetBinContent(i))/integral_1;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_1,2) + TMath::Power(h_1->GetBinContent(i)*errintegral_1/(integral_1*integral_1),2));

            val_2 = (h_2->GetBinContent(i))/max_h_2;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/max_h_2,2) + TMath::Power(h_2->GetBinContent(i)*errmax_h_2/(max_h_2*max_h_2),2));

            h_0->SetBinContent(i,val_0);
            h_0->SetBinError(i,errval_0);

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);

            // Difference
            if(val_2 != 0)
            {
                h_01->SetBinContent(i,(val_1 - val_0)/val_2);
                h_01->SetBinError(i,TMath::Sqrt(TMath::Power(errval_0/val_2,2) + TMath::Power(errval_1/val_2,2) + TMath::Power((val_1 - val_0)*errval_2/(val_2*val_2),2)));
            }
            else
            {
                h_01->SetBinContent(i,0);
                h_01->SetBinError(i,0);
            }
        }

        h_01->SetLineWidth(2);
        h_01->SetLineColor(kBlue);
        h_01->SetMaximum(maxY);
        h_01->SetMinimum(minY);
        h_01->GetXaxis()->SetTitle("Horizontal [pixels]");

        gStyle->SetOptStat(0);

        TString canvasname_1 = "./canvas/run28_run_28_run";
        canvasname_1 += runid[j];
        TCanvas* c_1 = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
        c_1->cd();
        h_01->SetTitle("Beam Profile along Vertical Axis");
        h_01->Draw();

        TLegend *legend = new TLegend(0.15,0.15,0.45,0.35,"","brNDC");

        // Difference
        TString name = "(RUN_";
        name += runid[j];
        name += " -- RUN_28)/RUN_28*";

        gPad->SetGrid();
        legend->AddEntry(h_01,name.Data(),"ple");
        legend->Draw();

        canvasname_1 += ".png";
        c_1->SaveAs(canvasname_1.Data());

        //-----------------------------------------------------------------//
        for(Int_t i = minCountBinX; i <= maxCountBinX; i++)
        {
            if(h_01->GetBinContent(i) > 0)
            {
                N1[j] += h_01->GetBinContent(i);
                err_N1[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(h_01->GetBinContent(i) < 0)
            {
                N2[j] += TMath::Abs(h_01->GetBinContent(i));
                err_N2[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            N12[j] += h_01->GetBinContent(i);
            err_N12[j] += TMath::Power(h_01->GetBinError(i),2);
        }
        err_N1[j] = TMath::Sqrt(err_N1[j]);
        err_N2[j] = TMath::Sqrt(err_N2[j]);
        err_N12[j] = TMath::Sqrt(err_N12[j]);
        //-----------------------------------------------------------------//
    }

    TCanvas* c_2 = new TCanvas(canvasname_2.Data(),canvasname_2.Data(),1920,1080);
    TGraphErrors* gr1 = new TGraphErrors(nn,position,N1,0,err_N1);
    TGraphErrors* gr1_line = new TGraphErrors(nn,position,N1,0,err_N1);
    TGraphErrors* gr2 = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr2_line = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr12 = new TGraphErrors(nn,position,N12,0,err_N12);
    TGraphErrors* gr12_line = new TGraphErrors(nn,position,N12,0,err_N12);
    gr1->SetLineColor(kRed);
    gr1->SetMarkerColor(kRed);
    gr1_line->SetLineColor(kRed);
    gr2->SetLineColor(kBlue);
    gr2->SetMarkerColor(kBlue);
    gr2_line->SetLineColor(kBlue);
    gr12->SetLineColor(kBlack);
    gr12->SetMarkerColor(kBlack);
    gr12_line->SetLineColor(kBlack);
    gr1_line->SetLineStyle(7);
    gr2_line->SetLineStyle(7);
    gr12_line->SetLineStyle(7);
    gr1->SetLineWidth(4);
    gr1_line->SetLineWidth(4);
    gr2->SetLineWidth(4);
    gr2_line->SetLineWidth(4);
    gr12->SetLineWidth(4);
    gr12_line->SetLineWidth(4);
    gr1->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    gr12->SetMarkerStyle(22);
    gr1->SetMarkerSize(2);
    gr2->SetMarkerSize(2);
    gr12->SetMarkerSize(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gr1,"AP");
    mg->Add(gr1_line,"AL");
    mg->Add(gr2,"AP");
    mg->Add(gr2_line,"AL");
    mg->Add(gr12,"AP");
    mg->Add(gr12_line,"AL");
    mg->SetMaximum(0.04);
    mg->SetMinimum(-0.01);
    mg->Draw("A");

    TLegend *legend = new TLegend(0.20,0.80,0.35,0.95,"","brNDC");
    gPad->SetGrid();
    legend->AddEntry(gr1,"N1","ple");
    legend->AddEntry(gr2,"N2","ple");
    legend->AddEntry(gr12,"N1 -- N2","ple");
    legend->Draw();

    canvasname_2 += ".png";
    c_2->SaveAs(canvasname_2.Data());

    return 0;
}

int histos_5()// to calculate areas above left/right and bellow the zero for all runs
{
    Int_t minCountBinX = 300;
    Int_t maxCountBinX = 500;

    // Diffuser
    const Int_t nn = 9;
    Int_t runid[] = {12, 10, 8, 6, 7, 9, 11, 13, 14};
    Double_t position[] = {66.970, 67.020, 67.045, 67.070, 67.095, 67.120, 67.170, 67.270, 67.570};
    TString canvasname_2 = "./canvas/diffuser";
    Int_t midCountBinX = 420;

    // Crystal Angular scan, fixed position
//    const Int_t nn = 14;
//    Int_t runid[] = {57,56,53,31,55,30,29,33,23,24,32,25,26,27};
//    Double_t position[] = {-170.0,-140.0,-110.0,-80.0,-60.0,-40.0,-20.0,-10.0,0.0,5.0,10.0,20.0,40.0,80.0};
//    TString canvasname_2 = "./canvas/angularscan";
//    Int_t midCountBinX = 400;

    // Crystal AM linear scan
//    const Int_t nn = 11;
//    Int_t runid[] = {22,52,21,49,20,16,17,47,18,50,19};
//    Double_t position[] = {51.750,52.200,52.250,52.300,52.550,52.750,52.950,53.200,53.250,53.300,53.750};
//    TString canvasname_2 = "./canvas/cryam";
//    Int_t midCountBinX = 420; // 420 for  j < 5 and 400 for j > 5

    // Crystal CH linear scan
//    const Int_t nn = 11;
//    Int_t runid[] = {46,43,38,41,36,23,35,40,37,42,45};
//    Double_t position[] = {51.750,52.200,52.250,52.300,52.550,52.750,52.950,53.200,53.250,53.300,53.750};
//    TString canvasname_2 = "./canvas/crych";
//    Int_t midCountBinX = 400;

    Double_t N1l[nn] = {}, err_N1l[nn] = {};
    Double_t N1r[nn] = {}, err_N1r[nn] = {};
    Double_t N1lr[nn] = {}, err_N1lr[nn] = {};
    Double_t N2[nn] = {}, err_N2[nn] = {};
    Double_t N12[nn] = {}, err_N12[nn] = {};
    // Difference
    Double_t maxY = 25.0e-2;
    Double_t minY = -25.0e-2;

    for(Int_t j = 0; j < nn; j++)
    {
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_";
        filename += runid[j];
        filename += ".root";

        TFile *_file0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_59.root");
        TFile *_file1 = TFile::Open(filename.Data());
        TFile *_file2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_28.root");

        TH2D* h_0_0 = (TH2D*)_file0->Get("h_91");
        TH2D* h_0_1 = (TH2D*)_file1->Get("h_91");
        TH2D* h_0_2 = (TH2D*)_file2->Get("h_91");

        removenoizypixelsXY(h_0_0,10);
        removenoizypixelsXY(h_0_1,10);
        removenoizypixelsXY(h_0_2,10);

        TH1D* h_0 = h_0_0->ProjectionX("h_92_0");
        TH1D* h_1 = h_0_1->ProjectionX("h_92_1");
        TH1D* h_2 = h_0_2->ProjectionX("h_92_2");
        TH1D* h_01 = h_0_0->ProjectionX("h_92_01");

        Double_t val_0, errval_0, integral_0, errintegral_0;
        Double_t val_1, errval_1, integral_1, errintegral_1;
        Double_t val_2, errval_2, integral_2, errintegral_2;

        integral_0 = h_0->IntegralAndError(1,h_0->GetNbinsX(),errintegral_0);
        integral_1 = h_1->IntegralAndError(1,h_1->GetNbinsX(),errintegral_1);
        integral_2 = h_2->IntegralAndError(1,h_2->GetNbinsX(),errintegral_2);

        for(Int_t i = 1; i <= h_0->GetNbinsX(); i++)
        {
            val_0 = (h_0->GetBinContent(i))/integral_0;
            errval_0 = TMath::Sqrt(TMath::Power(h_0->GetBinError(i)/integral_0,2) + TMath::Power(h_0->GetBinContent(i)*errintegral_0/(integral_0*integral_0),2));

            val_1 = (h_1->GetBinContent(i))/integral_1;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_1,2) + TMath::Power(h_1->GetBinContent(i)*errintegral_1/(integral_1*integral_1),2));

            val_2 = (h_2->GetBinContent(i))/integral_2;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/integral_2,2) + TMath::Power(h_2->GetBinContent(i)*errintegral_2/(integral_2*integral_2),2));

            h_0->SetBinContent(i,val_0);
            h_0->SetBinError(i,errval_0);

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);

            // Difference
            if(val_2 != 0)
            {
                h_01->SetBinContent(i,(val_1 - val_0)/val_2);
                h_01->SetBinError(i,TMath::Sqrt(TMath::Power(errval_0/val_2,2) + TMath::Power(errval_1/val_2,2) + TMath::Power((val_1 - val_0)*errval_2/(val_2*val_2),2)));
            }
            else
            {
                h_01->SetBinContent(i,0);
                h_01->SetBinError(i,0);
            }
        }

        h_01->SetLineWidth(2);
        h_01->SetLineColor(kBlue);
        h_01->SetMaximum(maxY);
        h_01->SetMinimum(minY);
        h_01->GetXaxis()->SetTitle("Horizontal [pixels]");

        gStyle->SetOptStat(0);

        TString canvasname_1 = "./canvas/run59_run_28_run";
        canvasname_1 += runid[j];
        TCanvas* c_1 = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
        c_1->cd();
        h_01->SetTitle("Beam Profile along Vertical Axis");
        h_01->Draw();

        TLegend *legend = new TLegend(0.15,0.15,0.45,0.35,"","brNDC");

        // Difference
        TString name = "(RUN_";
        name += runid[j];
        name += " -- RUN_59)/RUN_28";

        gPad->SetGrid();
        legend->AddEntry(h_01,name.Data(),"ple");
        legend->Draw();

        canvasname_1 += ".png";
        c_1->SaveAs(canvasname_1.Data());

        //-----------------------------------------------------------------//
        for(Int_t i = minCountBinX; i <= maxCountBinX; i++)
        {
            if(h_01->GetBinContent(i) > 0)
            {
                //if(j < 5){
                    if(i < midCountBinX)
                    {
                        N1l[j] += h_01->GetBinContent(i);
                        err_N1l[j] += TMath::Power(h_01->GetBinError(i),2);
                    }
                    else
                    {
                        N1r[j] += h_01->GetBinContent(i);
                        err_N1r[j] += TMath::Power(h_01->GetBinError(i),2);
                    }
                //}
                /*else
                {
                    if(i < midCountBinX-20)
                    {
                        N1l[j] += h_01->GetBinContent(i);
                        err_N1l[j] += TMath::Power(h_01->GetBinError(i),2);
                    }
                    else
                    {
                        N1r[j] += h_01->GetBinContent(i);
                        err_N1r[j] += TMath::Power(h_01->GetBinError(i),2);
                    }
                }*/
                N1lr[j] += h_01->GetBinContent(i);
                err_N1lr[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(h_01->GetBinContent(i) < 0)
            {
                N2[j] += TMath::Abs(h_01->GetBinContent(i));
                err_N2[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            N12[j] += h_01->GetBinContent(i);
            err_N12[j] += TMath::Power(h_01->GetBinError(i),2);
        }
        err_N1l[j] = TMath::Sqrt(err_N1l[j]);
        err_N1r[j] = TMath::Sqrt(err_N1r[j]);
        err_N1lr[j] = TMath::Sqrt(err_N1lr[j]);
        err_N2[j] = TMath::Sqrt(err_N2[j]);
        err_N12[j] = TMath::Sqrt(err_N12[j]);
        //-----------------------------------------------------------------//
    }

    TCanvas* c_2 = new TCanvas(canvasname_2.Data(),canvasname_2.Data(),1920,1080);
    TGraphErrors* gr1l = new TGraphErrors(nn,position,N1l,0,err_N1l);
    TGraphErrors* gr1l_line = new TGraphErrors(nn,position,N1l,0,err_N1l);
    TGraphErrors* gr1r = new TGraphErrors(nn,position,N1r,0,err_N1r);
    TGraphErrors* gr1r_line = new TGraphErrors(nn,position,N1r,0,err_N1r);
    TGraphErrors* gr1lr = new TGraphErrors(nn,position,N1lr,0,err_N1lr);
    TGraphErrors* gr1lr_line = new TGraphErrors(nn,position,N1lr,0,err_N1lr);
    TGraphErrors* gr2 = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr2_line = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr12 = new TGraphErrors(nn,position,N12,0,err_N12);
    TGraphErrors* gr12_line = new TGraphErrors(nn,position,N12,0,err_N12);
    gr1l->SetLineColor(kRed);
    gr1l->SetMarkerColor(kRed);
    gr1l_line->SetLineColor(kRed);
    gr1r->SetLineColor(kMagenta);
    gr1r->SetMarkerColor(kMagenta);
    gr1r_line->SetLineColor(kMagenta);
    gr1lr->SetLineColor(kGreen+1);
    gr1lr->SetMarkerColor(kGreen+1);
    gr1lr_line->SetLineColor(kGreen+1);
    gr2->SetLineColor(kBlue);
    gr2->SetMarkerColor(kBlue);
    gr2_line->SetLineColor(kBlue);
    gr12->SetLineColor(kBlack);
    gr12->SetMarkerColor(kBlack);
    gr12_line->SetLineColor(kBlack);
    gr1l_line->SetLineStyle(7);
    gr1r_line->SetLineStyle(7);
    gr1lr_line->SetLineStyle(7);
    gr2_line->SetLineStyle(7);
    gr12_line->SetLineStyle(7);
    gr1l->SetLineWidth(4);
    gr1l_line->SetLineWidth(4);
    gr1r->SetLineWidth(4);
    gr1r_line->SetLineWidth(4);
    gr1lr->SetLineWidth(4);
    gr1lr_line->SetLineWidth(4);
    gr2->SetLineWidth(4);
    gr2_line->SetLineWidth(4);
    gr12->SetLineWidth(4);
    gr12_line->SetLineWidth(4);
    gr1l->SetMarkerStyle(20);
    gr1r->SetMarkerStyle(20);
    gr1lr->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    gr12->SetMarkerStyle(22);
    gr1l->SetMarkerSize(2);
    gr1r->SetMarkerSize(2);
    gr1lr->SetMarkerSize(2);
    gr2->SetMarkerSize(2);
    gr12->SetMarkerSize(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gr1l,"AP");
    mg->Add(gr1l_line,"AL");
    mg->Add(gr1r,"AP");
    mg->Add(gr1r_line,"AL");
    mg->Add(gr1lr,"AP");
    mg->Add(gr1lr_line,"AL");
    mg->Add(gr2,"AP");
    mg->Add(gr2_line,"AL");
    mg->Add(gr12,"AP");
    mg->Add(gr12_line,"AL");
    mg->SetMaximum(5.0);
    mg->SetMinimum(-1.0);
    mg->Draw("A");

    TLegend *legend = new TLegend(0.80,0.80,0.95,0.95,"","brNDC");
    gPad->SetGrid();
    legend->AddEntry(gr1l,"N1L","ple");
    legend->AddEntry(gr1r,"N1R","ple");
    legend->AddEntry(gr1lr,"N1L + N1R","ple");
    legend->AddEntry(gr2,"N2","ple");
    legend->AddEntry(gr12,"N1 -- N2","ple");
    legend->Draw();

    canvasname_2 += ".png";
    c_2->SaveAs(canvasname_2.Data());

    return 0;
}

int histos_6()// to calculate areas above left/right and bellow the zero for the defined runs with some obstacles on the beam
{
    Int_t minCountBinX = 300;
    Int_t midCountBinX = 410;
    Int_t maxCountBinX = 500;

    const Int_t nn = 6;
    Int_t runid[] = {
        15,//only crystal am
        59,//only septum
        8,//septum + diffuser
        16,//septum + crystal am
        23,//septum + crystal ch
        31//septum + crystal vr
    };
    Double_t position[] = {0,1,2,3,4,5,6};
    TString canvasname_2 = "./canvas/comparison";

    Double_t N1l[nn] = {}, err_N1l[nn] = {};
    Double_t N1r[nn] = {}, err_N1r[nn] = {};
    Double_t N1lr[nn] = {}, err_N1lr[nn] = {};
    Double_t N2[nn] = {}, err_N2[nn] = {};
    Double_t N12[nn] = {}, err_N12[nn] = {};
    // Difference
    Double_t maxY = 25.0e-2;
    Double_t minY = -25.0e-2;

    for(Int_t j = 0; j < nn; j++)
    {
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_";
        filename += runid[j];
        filename += ".root";

        TFile *_file0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_28.root");
        TFile *_file1 = TFile::Open(filename.Data());
        TFile *_file2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_28.root");

        TH2D* h_0_0 = (TH2D*)_file0->Get("h_91");
        TH2D* h_0_1 = (TH2D*)_file1->Get("h_91");
        TH2D* h_0_2 = (TH2D*)_file2->Get("h_91");

        removenoizypixelsXY(h_0_0,10);
        removenoizypixelsXY(h_0_1,10);
        removenoizypixelsXY(h_0_2,10);

        TH1D* h_0 = h_0_0->ProjectionX("h_92_0");
        TH1D* h_1 = h_0_1->ProjectionX("h_92_1");
        TH1D* h_2 = h_0_2->ProjectionX("h_92_2");
        TH1D* h_01 = h_0_0->ProjectionX("h_92_01");

        Double_t val_0, errval_0, integral_0, errintegral_0;
        Double_t val_1, errval_1, integral_1, errintegral_1;
        Double_t val_2, errval_2, integral_2, errintegral_2;

        integral_0 = h_0->IntegralAndError(1,h_0->GetNbinsX(),errintegral_0);
        integral_1 = h_1->IntegralAndError(1,h_1->GetNbinsX(),errintegral_1);
        integral_2 = h_2->IntegralAndError(1,h_2->GetNbinsX(),errintegral_2);

        for(Int_t i = 1; i <= h_0->GetNbinsX(); i++)
        {
            val_0 = (h_0->GetBinContent(i))/integral_0;
            errval_0 = TMath::Sqrt(TMath::Power(h_0->GetBinError(i)/integral_0,2) + TMath::Power(h_0->GetBinContent(i)*errintegral_0/(integral_0*integral_0),2));

            val_1 = (h_1->GetBinContent(i))/integral_1;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_1,2) + TMath::Power(h_1->GetBinContent(i)*errintegral_1/(integral_1*integral_1),2));

            val_2 = (h_2->GetBinContent(i))/integral_2;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/integral_2,2) + TMath::Power(h_2->GetBinContent(i)*errintegral_2/(integral_2*integral_2),2));

            h_0->SetBinContent(i,val_0);
            h_0->SetBinError(i,errval_0);

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);

            // Difference
            if(val_2 != 0)
            {
                h_01->SetBinContent(i,(val_1 - val_0)/val_2);
                h_01->SetBinError(i,TMath::Sqrt(TMath::Power(errval_0/val_2,2) + TMath::Power(errval_1/val_2,2) + TMath::Power((val_1 - val_0)*errval_2/(val_2*val_2),2)));
            }
            else
            {
                h_01->SetBinContent(i,0);
                h_01->SetBinError(i,0);
            }
        }

        h_01->SetLineWidth(2);
        h_01->SetLineColor(kBlue);
        h_01->SetMaximum(maxY);
        h_01->SetMinimum(minY);
        h_01->GetXaxis()->SetTitle("Horizontal [pixels]");

        gStyle->SetOptStat(0);

        TString canvasname_1 = "./canvas/run28_run_28_run";
        canvasname_1 += runid[j];
        TCanvas* c_1 = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
        c_1->cd();
        h_01->SetTitle("Beam Profile along Vertical Axis");
        h_01->Draw();

        TLegend *legend = new TLegend(0.15,0.15,0.45,0.35,"","brNDC");

        // Difference
        TString name = "(RUN_";
        name += runid[j];
        name += " -- RUN_28)/RUN_28";

        gPad->SetGrid();
        legend->AddEntry(h_01,name.Data(),"ple");
        legend->Draw();

        canvasname_1 += ".png";
        c_1->SaveAs(canvasname_1.Data());

        //-----------------------------------------------------------------//
        for(Int_t i = minCountBinX; i <= maxCountBinX; i++)
        {
            if(h_01->GetBinContent(i) > 0)
            {
                if(i < midCountBinX)
                {
                    N1l[j] += h_01->GetBinContent(i);
                    err_N1l[j] += TMath::Power(h_01->GetBinError(i),2);
                }
                else
                {
                    N1r[j] += h_01->GetBinContent(i);
                    err_N1r[j] += TMath::Power(h_01->GetBinError(i),2);
                }
                N1lr[j] += h_01->GetBinContent(i);
                err_N1lr[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(h_01->GetBinContent(i) < 0)
            {
                N2[j] += TMath::Abs(h_01->GetBinContent(i));
                err_N2[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            N12[j] += h_01->GetBinContent(i);
            err_N12[j] += TMath::Power(h_01->GetBinError(i),2);
        }
        err_N1l[j] = TMath::Sqrt(err_N1l[j]);
        err_N1r[j] = TMath::Sqrt(err_N1r[j]);
        err_N1lr[j] = TMath::Sqrt(err_N1lr[j]);
        err_N2[j] = TMath::Sqrt(err_N2[j]);
        err_N12[j] = TMath::Sqrt(err_N12[j]);
        //-----------------------------------------------------------------//
    }

    TCanvas* c_2 = new TCanvas(canvasname_2.Data(),canvasname_2.Data(),1920,1080);
    TGraphErrors* gr1l = new TGraphErrors(nn,position,N1l,0,err_N1l);
    TGraphErrors* gr1l_line = new TGraphErrors(nn,position,N1l,0,err_N1l);
    TGraphErrors* gr1r = new TGraphErrors(nn,position,N1r,0,err_N1r);
    TGraphErrors* gr1r_line = new TGraphErrors(nn,position,N1r,0,err_N1r);
    TGraphErrors* gr1lr = new TGraphErrors(nn,position,N1lr,0,err_N1lr);
    TGraphErrors* gr1lr_line = new TGraphErrors(nn,position,N1lr,0,err_N1lr);
    TGraphErrors* gr2 = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr2_line = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr12 = new TGraphErrors(nn,position,N12,0,err_N12);
    TGraphErrors* gr12_line = new TGraphErrors(nn,position,N12,0,err_N12);
    gr1l->SetLineColor(kRed);
    gr1l->SetMarkerColor(kRed);
    gr1l_line->SetLineColor(kRed);
    gr1r->SetLineColor(kMagenta);
    gr1r->SetMarkerColor(kMagenta);
    gr1r_line->SetLineColor(kMagenta);
    gr1lr->SetLineColor(kGreen+1);
    gr1lr->SetMarkerColor(kGreen+1);
    gr1lr_line->SetLineColor(kGreen+1);
    gr2->SetLineColor(kBlue);
    gr2->SetMarkerColor(kBlue);
    gr2_line->SetLineColor(kBlue);
    gr12->SetLineColor(kBlack);
    gr12->SetMarkerColor(kBlack);
    gr12_line->SetLineColor(kBlack);
    gr1l_line->SetLineStyle(7);
    gr1r_line->SetLineStyle(7);
    gr1lr_line->SetLineStyle(7);
    gr2_line->SetLineStyle(7);
    gr12_line->SetLineStyle(7);
    gr1l->SetLineWidth(4);
    gr1l_line->SetLineWidth(4);
    gr1r->SetLineWidth(4);
    gr1r_line->SetLineWidth(4);
    gr1lr->SetLineWidth(4);
    gr1lr_line->SetLineWidth(4);
    gr2->SetLineWidth(4);
    gr2_line->SetLineWidth(4);
    gr12->SetLineWidth(4);
    gr12_line->SetLineWidth(4);
    gr1l->SetMarkerStyle(20);
    gr1r->SetMarkerStyle(20);
    gr1lr->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    gr12->SetMarkerStyle(22);
    gr1l->SetMarkerSize(2);
    gr1r->SetMarkerSize(2);
    gr1lr->SetMarkerSize(2);
    gr2->SetMarkerSize(2);
    gr12->SetMarkerSize(2);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gr1l,"AP");
    mg->Add(gr1l_line,"AL");
    mg->Add(gr1r,"AP");
    mg->Add(gr1r_line,"AL");
    mg->Add(gr1lr,"AP");
    mg->Add(gr1lr_line,"AL");
    mg->Add(gr2,"AP");
    mg->Add(gr2_line,"AL");
    mg->Add(gr12,"AP");
    mg->Add(gr12_line,"AL");
    mg->SetMaximum(6.0);
    mg->SetMinimum(-1.0);
    mg->Draw("A");

    TLegend *legend = new TLegend(0.20,0.80,0.35,0.95,"","brNDC");
    gPad->SetGrid();
    legend->AddEntry(gr1l,"N1L","ple");
    legend->AddEntry(gr1r,"N1R","ple");
    legend->AddEntry(gr1lr,"N1L + N1R","ple");
    legend->AddEntry(gr2,"N2","ple");
    legend->AddEntry(gr12,"N1 -- N2","ple");
    legend->Draw();

    canvasname_2 += ".png";
    c_2->SaveAs(canvasname_2.Data());

    return 0;
}

int histos_7()// polar coordinate convertor
{
    for(Int_t j = 1; j <= 66; j++)
    {
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_";
        filename += j;
        filename += ".root";
        TFile *_file_1 = TFile::Open(filename.Data());

        TH2D* h_0 = (TH2D*)_file_1->Get("h_91");
        removenoizypixelsXY(h_0,10);
        TH1D* h_1 = h_0->ProjectionX("h_92");
        TH1D* h_2 = h_0->ProjectionY("h_93");
        TH1D* h_3 = new TH1D("h_3","Polar Angle #Phi [rad]",2.2*TMath::Pi()*1000,-1.1*TMath::Pi(),1.1*TMath::Pi());
        TH2D* h_4 = new TH2D("h_4","Y vs X vs COUNTS",300,-150,150,300,-150,150);

        Double_t integral_0;
        Double_t val_1, errval_1, integral_1, errintegral_1;
        Double_t val_2, errval_2, integral_2, errintegral_2;

        integral_0 = h_0->Integral();
        integral_1 = h_1->IntegralAndError(1,h_1->GetNbinsX(),errintegral_1);
        integral_2 = h_2->IntegralAndError(1,h_2->GetNbinsX(),errintegral_2);

        for(Int_t i = 1; i <= h_1->GetNbinsX(); i++)
        {
            val_1 = (h_1->GetBinContent(i))/integral_1;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_1,2) + TMath::Power(h_1->GetBinContent(i)*errintegral_1/(integral_1*integral_1),2));

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);
        }

        for(Int_t i = 1; i <= h_2->GetNbinsX(); i++)
        {
            val_2 = (h_2->GetBinContent(i))/integral_2;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/integral_2,2) + TMath::Power(h_2->GetBinContent(i)*errintegral_2/(integral_2*integral_2),2));

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);
        }
        h_0->Scale(1.0/integral_0);

        //---------------------------------------//
        TF1* fit_funct = new TF1("fit_funct","gaus",150,250);
        h_2->Fit(fit_funct,"RQ+");
        Int_t meanY = fit_funct->GetParameter(1);
        Int_t meanX = 412;
        Int_t radius = 100;
        for(Int_t xi = 1; xi <= h_0->GetNbinsX(); xi++)
        {
            for(Int_t yi = 1; yi <= h_0->GetNbinsY(); yi++)
            {
                if(xi-meanX == 0 || yi-meanY == 0) continue;
                if(TMath::Abs(xi-meanX) == TMath::Abs(yi-meanY)) continue;
                if(((xi-meanX)*(xi-meanX) + (yi-meanY)*(yi-meanY)) < radius*radius)
                {
                    h_3->Fill(TMath::ATan2((double)(yi-meanY),(double)(xi-meanX)));
                    h_4->Fill(xi-meanX,yi-meanY,h_0->GetBinContent(xi,yi));
                }
            }
        }
        //---------------------------------------//
        h_0->SetMaximum(60e-6);

        h_1->SetLineWidth(4);
        h_1->SetMaximum(10e-3);

        h_2->SetLineWidth(4);
        h_2->SetMaximum(10e-3);

        h_3->SetLineWidth(4);

        h_4->SetMaximum(60e-6);

        gStyle->SetOptStat(0);

        const Int_t NRGBs = 5;
        const Int_t NCont = 256;
        Double_t stops[NRGBs] = { 0.00, 0.30, 0.61, 0.84, 1.00 };
        Double_t red[NRGBs] = { 0.00, 0.00, 0.57, 0.90, 0.51 };
        Double_t green[NRGBs] = { 0.00, 0.65, 0.95, 0.20, 0.00 };
        Double_t blue[NRGBs] = { 0.51, 0.55, 0.15, 0.00, 0.10 };
        TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        gStyle->SetNumberContours(NCont);

        TString canvasname_0 = "./canvas/hxy_run_";
        canvasname_0 += j;
        TCanvas* c_0 = new TCanvas(canvasname_0.Data(),canvasname_0.Data(),1080,1080);
        c_0->cd();
        h_0->Draw("colz");
        canvasname_0 += ".png";
        c_0->Update();
        TPaletteAxis *palette = (TPaletteAxis*)h_0->GetListOfFunctions()->FindObject("palette");
        palette->SetX1NDC(0.905);
        palette->SetX2NDC(0.925);
        palette->SetY1NDC(0.1);
        palette->SetY2NDC(0.9);
        c_0->Modified();
        c_0->Update();
        c_0->SaveAs(canvasname_0.Data());

        TString canvasname_1 = "./canvas/hx_run_";
        canvasname_1 += j;
        TCanvas* c_1 = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
        c_1->cd();
        h_1->Draw();
        canvasname_1 += ".png";
        c_1->SaveAs(canvasname_1.Data());

        TString canvasname_2 = "./canvas/hy_run_";
        canvasname_2 += j;
        TCanvas* c_2 = new TCanvas(canvasname_2.Data(),canvasname_2.Data(),1920,1080);
        c_2->cd();
        h_2->Draw();
        canvasname_2 += ".png";
        c_2->SaveAs(canvasname_2.Data());

        TString canvasname_3 = "./canvas/pol_run_";
        canvasname_3 += j;
        TCanvas* c_3 = new TCanvas(canvasname_3.Data(),canvasname_3.Data(),1920,1080);
        c_3->cd();
        h_3->Draw();
        canvasname_3 += ".png";
        c_3->SaveAs(canvasname_3.Data());

        TString canvasname_4 = "./canvas/hxy_cut_run_";
        canvasname_4 += j;
        TCanvas* c_4 = new TCanvas(canvasname_4.Data(),canvasname_4.Data(),1080,1080);
        c_4->cd();
        h_4->Draw("colz");
        canvasname_4 += ".png";
        c_4->Update();
        TPaletteAxis *palette_4 = (TPaletteAxis*)h_4->GetListOfFunctions()->FindObject("palette");
        palette_4->SetX1NDC(0.905);
        palette_4->SetX2NDC(0.925);
        palette_4->SetY1NDC(0.1);
        palette_4->SetY2NDC(0.9);
        c_4->Modified();
        c_4->Update();
        c_4->SaveAs(canvasname_4.Data());
    }

    return 0;
}

int histos_8()// to plot the difference between runs in 3D
{
    TFile *_file_1 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_17.root");
    TFile *_file_2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_1.root");

    TH2D* h_1 = (TH2D*)_file_1->Get("h_1");
    TH2D* h_2 = (TH2D*)_file_2->Get("h_1");

    removenoizypixelsXY(h_1,10);
    removenoizypixelsXY(h_2,10);

    TH1D* h_1y = h_1->ProjectionY("h_3_1");
    TH1D* h_2y = h_2->ProjectionY("h_3_2");

    Double_t integral_1 = h_1->Integral();
    Double_t integral_2 = h_2->Integral();

    TH2D* h_dif = new TH2D("h_dif","Difference",512,0,512,512,0,512);

    h_1->Scale(1.0/integral_1);
    h_2->Scale(1.0/integral_2);

    TF1* fit_funct = new TF1("fit_funct","gaus",150,250);
    h_1y->Fit(fit_funct,"RQ+");
    Int_t mean1Y = 0;//fit_funct->GetParameter(1);
    h_2y->Fit(fit_funct,"RQ+");
    Int_t mean2Y = 0;//fit_funct->GetParameter(1);

    for(Int_t xi = 1; xi <= h_1->GetNbinsX(); xi++)
    {
        for(Int_t yi = 1; yi <= h_1->GetNbinsY(); yi++)
        {
            Double_t val = h_1->GetBinContent(xi,yi) - h_2->GetBinContent(xi,yi+(mean2Y-mean1Y));
            if(val > 0 && h_1->GetBinContent(xi,yi) > 0 && h_2->GetBinContent(xi,yi+(mean2Y-mean1Y)) > 0)
                h_dif->SetBinContent(xi,yi,val);
        }
    }

    gStyle->SetOptStat(0);

    TCanvas* c_0 = new TCanvas("c_0","c_0",1080,1080);
    c_0->cd();
    gPad->SetGrid();
    h_dif->SetMaximum(1e-4);
    h_dif->Draw("colz");

    TCanvas* c_1 = new TCanvas("c_1","c_1",1080,1080);
    c_1->cd();
    gPad->SetGrid();
    h_1y->Draw();

    TCanvas* c_2 = new TCanvas("c_2","c_2",1080,1080);
    c_2->cd();
    gPad->SetGrid();
    h_2y->Draw();

    TCanvas* c_3 = new TCanvas("c_3","c_3",1080,1080);
    c_3->cd();
    gPad->SetGrid();
    h_1->SetMaximum(1e-4);
    h_1->Draw("colz");

    TCanvas* c_4 = new TCanvas("c_4","c_4",1080,1080);
    c_4->cd();
    gPad->SetGrid();
    h_2->SetMaximum(1e-4);
    h_2->Draw("colz");

    return 0;
}

int histos_9()// to calculate specified areas for the defined runs with some obstacles on the beam
{
    const Int_t nn = 6;
    Int_t runid[] = {
        15,//only crystal am
        59,//only septum
        8,//septum + diffuser
        16,//septum + crystal am
        23,//septum + crystal ch
        31//septum + crystal vr
    };
    Double_t position[] = {0,1,2,3,4,5,6};

    Int_t CountBinX1[] = {370,370,370,370,333,360};
    Int_t CountBinX2[] = {400,401,399,401,391,391};
    Int_t CountBinX3[] = {425,425,426,425,425,422};
    Int_t CountBinX4[] = {450,470,470,460,460,460};

    TLine* line1;
    TLine* line2;
    TLine* line3;
    TLine* line4;

    TString canvasname_2 = "./canvas/comparison";

    Double_t N1[nn] = {}, err_N1[nn] = {};
    Double_t N2[nn] = {}, err_N2[nn] = {};
    Double_t N3[nn] = {}, err_N3[nn] = {};
    Double_t N4[nn] = {}, err_N4[nn] = {};
    Double_t N5[nn] = {}, err_N5[nn] = {};
    Double_t N6[nn] = {}, err_N6[nn] = {};
    // Difference
    Double_t maxY = 20.0e-4;
    Double_t minY = -20.0e-4;

    for(Int_t j = 0; j < nn; j++)
    {
        TString filename = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_";
        filename += runid[j];
        filename += ".root";

        TFile *_file0 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_28.root");
        TFile *_file1 = TFile::Open(filename.Data());
        TFile *_file2 = TFile::Open("/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_28.root");

        TH2D* h_0_0 = (TH2D*)_file0->Get("h_91");
        TH2D* h_0_1 = (TH2D*)_file1->Get("h_91");
        TH2D* h_0_2 = (TH2D*)_file2->Get("h_91");

        removenoizypixelsXY(h_0_0,10);
        removenoizypixelsXY(h_0_1,10);
        removenoizypixelsXY(h_0_2,10);

        TH1D* h_0 = h_0_0->ProjectionX("h_92_0");
        TH1D* h_1 = h_0_1->ProjectionX("h_92_1");
        TH1D* h_2 = h_0_2->ProjectionX("h_92_2");
        TH1D* h_01 = h_0_0->ProjectionX("h_92_01");

        Double_t val_0, errval_0, integral_0, errintegral_0;
        Double_t val_1, errval_1, integral_1, errintegral_1;
        Double_t val_2, errval_2, max_h_2, errmax_h_2;

        integral_0 = h_0->IntegralAndError(1,h_0->GetNbinsX(),errintegral_0);
        integral_1 = h_1->IntegralAndError(1,h_1->GetNbinsX(),errintegral_1);

        max_h_2 = h_2->GetBinContent(h_2->GetMaximumBin());
        errmax_h_2 = h_2->GetBinError(h_2->GetMaximumBin());

        Int_t min_bin = -1;
        Double_t min_val = 9999.9999, min_val_err = 9999.9999;

        for(Int_t i = 1; i <= h_0->GetNbinsX(); i++)
        {
            val_0 = (h_0->GetBinContent(i))/integral_0;
            errval_0 = TMath::Sqrt(TMath::Power(h_0->GetBinError(i)/integral_0,2) + TMath::Power(h_0->GetBinContent(i)*errintegral_0/(integral_0*integral_0),2));

            val_1 = (h_1->GetBinContent(i))/integral_1;
            errval_1 = TMath::Sqrt(TMath::Power(h_1->GetBinError(i)/integral_1,2) + TMath::Power(h_1->GetBinContent(i)*errintegral_1/(integral_1*integral_1),2));

            val_2 = (h_2->GetBinContent(i))/max_h_2;
            errval_2 = TMath::Sqrt(TMath::Power(h_2->GetBinError(i)/max_h_2,2) + TMath::Power(h_2->GetBinContent(i)*errmax_h_2/(max_h_2*max_h_2),2));

            h_0->SetBinContent(i,val_0);
            h_0->SetBinError(i,errval_0);

            h_1->SetBinContent(i,val_1);
            h_1->SetBinError(i,errval_1);

            h_2->SetBinContent(i,val_2);
            h_2->SetBinError(i,errval_2);

            // Difference
            if(val_2 != 0)
            {
                h_01->SetBinContent(i,(val_1 - val_0)/val_2);
                h_01->SetBinError(i,TMath::Sqrt(TMath::Power(errval_0/val_2,2) + TMath::Power(errval_1/val_2,2) + TMath::Power((val_1 - val_0)*errval_2/(val_2*val_2),2)));

                if(min_val > (val_1 - val_0)/val_2 && i > CountBinX1[j] && i < CountBinX4[j])
                {
                    min_val = (val_1 - val_0)/val_2;
                    min_val_err = TMath::Sqrt(TMath::Power(errval_0/val_2,2) + TMath::Power(errval_1/val_2,2) + TMath::Power((val_1 - val_0)*errval_2/(val_2*val_2),2));
                    min_bin = i;
                }
            }
            else
            {
                h_01->SetBinContent(i,0);
                h_01->SetBinError(i,0);
            }
        }

        h_01->SetLineWidth(2);
        h_01->SetLineColor(kBlue);
        h_01->SetMaximum(maxY);
        h_01->SetMinimum(minY);
        h_01->GetXaxis()->SetTitle("Horizontal [pixels]");

        cout<<"--> Global minimum: "<<min_val<<" +/- "<<min_val_err<<" is at "<<min_bin<<endl;

        gStyle->SetOptStat(0);

        TString canvasname_1 = "./canvas/run28_run_28_run";
        canvasname_1 += runid[j];
        TCanvas* c_1 = new TCanvas(canvasname_1.Data(),canvasname_1.Data(),1920,1080);
        c_1->cd();
        h_01->SetTitle("Beam Profile along Vertical Axis");
        h_01->Draw();
        line1 = new TLine(CountBinX1[j],minY,CountBinX1[j],maxY);
        line2 = new TLine(CountBinX2[j],minY,CountBinX2[j],maxY);
        line3 = new TLine(CountBinX3[j],minY,CountBinX3[j],maxY);
        line4 = new TLine(CountBinX4[j],minY,CountBinX4[j],maxY);
        line1->SetLineColor(kRed);
        line2->SetLineColor(kBlack);
        line3->SetLineColor(kGreen+2);
        line4->SetLineColor(kMagenta);
        line1->SetLineWidth(4);
        line2->SetLineWidth(4);
        line3->SetLineWidth(4);
        line4->SetLineWidth(4);
        line1->Draw("same");
        line2->Draw("same");
        line3->Draw("same");
        line4->Draw("same");

        TLegend *legend = new TLegend(0.15,0.15,0.45,0.35,"","brNDC");

        // Difference
        TString name = "(RUN_";
        name += runid[j];
        name += " -- RUN_28)/RUN_28*";

        gPad->SetGrid();
        legend->AddEntry(h_01,name.Data(),"ple");
        legend->Draw();

        canvasname_1 += ".png";
        c_1->SaveAs(canvasname_1.Data());

        //-----------------------------------------------------------------//
        for(Int_t i = 1; i <= h_01->GetNbinsX(); i++)
        {
            if(i < CountBinX1[j])
            {
                N1[j] += h_01->GetBinContent(i);
                err_N1[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(i > CountBinX1[j] && i < CountBinX2[j])
            {
                N2[j] += h_01->GetBinContent(i);
                err_N2[j] += TMath::Power(h_01->GetBinError(i),2);

                N6[j] += h_01->GetBinContent(i);
                err_N6[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(i > CountBinX2[j] && i < CountBinX3[j])
            {
                N3[j] += TMath::Abs(h_01->GetBinContent(i));
                err_N3[j] += TMath::Power(h_01->GetBinError(i),2);

                N6[j] += h_01->GetBinContent(i);
                err_N6[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(i > CountBinX3[j] && i < CountBinX4[j])
            {
                N4[j] += h_01->GetBinContent(i);
                err_N4[j] += TMath::Power(h_01->GetBinError(i),2);

                N6[j] += h_01->GetBinContent(i);
                err_N6[j] += TMath::Power(h_01->GetBinError(i),2);
            }
            if(i > CountBinX4[j])
            {
                N5[j] += h_01->GetBinContent(i);
                err_N5[j] += TMath::Power(h_01->GetBinError(i),2);
            }
        }

        err_N1[j] = TMath::Sqrt(err_N1[j]);
        err_N2[j] = TMath::Sqrt(err_N2[j]);
        err_N3[j] = TMath::Sqrt(err_N3[j]);
        err_N4[j] = TMath::Sqrt(err_N4[j]);
        err_N5[j] = TMath::Sqrt(err_N5[j]);
        err_N6[j] = TMath::Sqrt(err_N6[j]);
        //-----------------------------------------------------------------//
    }

    TCanvas* c_2 = new TCanvas(canvasname_2.Data(),canvasname_2.Data(),1920,1080);

    TGraphErrors* gr1 = new TGraphErrors(nn,position,N1,0,err_N1);
    TGraphErrors* gr1_line = new TGraphErrors(nn,position,N1,0,err_N1);
    TGraphErrors* gr2 = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr2_line = new TGraphErrors(nn,position,N2,0,err_N2);
    TGraphErrors* gr3 = new TGraphErrors(nn,position,N3,0,err_N3);
    TGraphErrors* gr3_line = new TGraphErrors(nn,position,N3,0,err_N3);
    TGraphErrors* gr4 = new TGraphErrors(nn,position,N4,0,err_N4);
    TGraphErrors* gr4_line = new TGraphErrors(nn,position,N4,0,err_N4);
    TGraphErrors* gr5 = new TGraphErrors(nn,position,N5,0,err_N5);
    TGraphErrors* gr5_line = new TGraphErrors(nn,position,N5,0,err_N5);
    TGraphErrors* gr6 = new TGraphErrors(nn,position,N6,0,err_N6);
    TGraphErrors* gr6_line = new TGraphErrors(nn,position,N6,0,err_N6);

    gr1->SetLineColor(kRed);
    gr1->SetMarkerColor(kRed);
    gr1_line->SetLineColor(kRed);
    gr2->SetLineColor(kMagenta);
    gr2->SetMarkerColor(kMagenta);
    gr2_line->SetLineColor(kMagenta);
    gr3->SetLineColor(kGreen+1);
    gr3->SetMarkerColor(kGreen+1);
    gr3_line->SetLineColor(kGreen+1);
    gr4->SetLineColor(kBlue);
    gr4->SetMarkerColor(kBlue);
    gr4_line->SetLineColor(kBlue);
    gr5->SetLineColor(kOrange);
    gr5->SetMarkerColor(kOrange);
    gr5_line->SetLineColor(kOrange);
    gr6->SetLineColor(kBlack);
    gr6->SetMarkerColor(kBlack);
    gr6_line->SetLineColor(kBlack);

    gr1_line->SetLineStyle(7);
    gr2_line->SetLineStyle(7);
    gr3_line->SetLineStyle(7);
    gr4_line->SetLineStyle(7);
    gr5_line->SetLineStyle(7);
    gr6_line->SetLineStyle(7);

    gr1->SetLineWidth(4);
    gr1_line->SetLineWidth(4);
    gr2->SetLineWidth(4);
    gr2_line->SetLineWidth(4);
    gr3->SetLineWidth(4);
    gr3_line->SetLineWidth(4);
    gr4->SetLineWidth(4);
    gr4_line->SetLineWidth(4);
    gr5->SetLineWidth(4);
    gr5_line->SetLineWidth(4);
    gr6->SetLineWidth(4);
    gr6_line->SetLineWidth(4);

    gr1->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    gr3->SetMarkerStyle(22);
    gr4->SetMarkerStyle(23);
    gr5->SetMarkerStyle(24);
    gr6->SetMarkerStyle(25);

    gr1->SetMarkerSize(2);
    gr2->SetMarkerSize(2);
    gr3->SetMarkerSize(2);
    gr4->SetMarkerSize(2);
    gr5->SetMarkerSize(2);
    gr6->SetMarkerSize(2);

    TMultiGraph* mg = new TMultiGraph();
//    mg->Add(gr1,"AP");
//    mg->Add(gr1_line,"AL");
    mg->Add(gr2,"AP");
    mg->Add(gr2_line,"AL");
    mg->Add(gr3,"AP");
    mg->Add(gr3_line,"AL");
    mg->Add(gr4,"AP");
    mg->Add(gr4_line,"AL");
//    mg->Add(gr5,"AP");
//    mg->Add(gr5_line,"AL");
    mg->Add(gr6,"AP");
    mg->Add(gr6_line,"AL");
//    mg->SetMaximum(0.04);
//    mg->SetMinimum(-0.01);
    mg->Draw("A");

    TLegend *legend = new TLegend(0.20,0.70,0.45,0.95,"","brNDC");
    gPad->SetGrid();
//    legend->AddEntry(gr1,"N_{1}","ple");
    legend->AddEntry(gr2,"N_{2} (Left bump)","ple");
    legend->AddEntry(gr3,"N_{3} (Valley)","ple");
    legend->AddEntry(gr4,"N_{4} (Right bump)","ple");
//    legend->AddEntry(gr5,"N_{5}","ple");
    legend->AddEntry(gr6,"N_{2} + N_{4} - N_{3}","ple");
    legend->Draw();

    canvasname_2 += ".png";
    c_2->SaveAs(canvasname_2.Data());

    return 0;
}

int histos_10()
{
    const Int_t nn = 6;

    Double_t pos_left[]     = {388,389,389,389,367,375};
    Double_t pos_right[]    = {433,434,437,436,433,433};
    Double_t pos_mid[]      = {414,413,414,414,413,408};

    Double_t position[] = {0,1,2,3,4,5,6};

    for(Int_t i = 0; i < nn; i++)
    {
        pos_left[i]     = (pos_left[i] - 414)*0.055;
        pos_right[i]    = (pos_right[i] - 414)*0.055;
        pos_mid[i]      = (pos_mid[i] - 414)*0.055;
    }

    TGraph* gr_1 = new TGraph(nn,position,pos_left);
    TGraph* gr_2 = new TGraph(nn,position,pos_right);
    TGraph* gr_3 = new TGraph(nn,position,pos_mid);

    TCanvas* c1 = new TCanvas();
    c1->cd();
    TMultiGraph* mg = new TMultiGraph();

    gr_1->SetLineColor(kRed);
    gr_2->SetLineColor(kBlue);
    gr_3->SetLineColor(kGreen+2);

    gr_1->SetLineWidth(4);
    gr_2->SetLineWidth(4);
    gr_3->SetLineWidth(4);

    gr_1->SetMarkerStyle(20);
    gr_2->SetMarkerStyle(21);
    gr_3->SetMarkerStyle(22);

    gr_1->SetMarkerSize(2);
    gr_2->SetMarkerSize(2);
    gr_3->SetMarkerSize(2);

    gr_1->SetMarkerColor(kRed);
    gr_2->SetMarkerColor(kBlue);
    gr_3->SetMarkerColor(kGreen+2);

    mg->Add(gr_1);
    mg->Add(gr_2);
    mg->Add(gr_3);

    mg->Draw("APL");

    TLegend *legend = new TLegend(0.20,0.70,0.45,0.95,"","brNDC");
    gPad->SetGrid();
    legend->AddEntry(gr_2,"Right Bump","pl");
    legend->AddEntry(gr_3,"Valley","pl");
    legend->AddEntry(gr_1,"Left Bump","pl");
    legend->Draw();

    c1->SaveAs("positions.png");

    return 0;
}

int histos_11()// to plot the beam movement
{
    Int_t i = 0;
    const Int_t nRuns = 24;
    Double_t mean_x[nRuns], mean_y[nRuns], mean_err_x[nRuns], mean_err_y[nRuns];
    TCanvas* c_x[nRuns];
    TCanvas* c_y[nRuns];
    Double_t ttime[nRuns], ttime_err[nRuns];

    gStyle->SetOptStat(0);

    for(Int_t run_id = 1; run_id <= 26; run_id++)
    {
        if(run_id == 18 || run_id == 19) continue; // runs in ToA mode

        TString file_name = "/home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_";
        file_name += run_id;
        file_name += ".root";
        TFile *_file = TFile::Open(file_name.Data());
        cout<<_file->GetName()<<endl;
        TH2D* h_1 = (TH2D*)_file->Get("h_1");
        TH2D* h_10 = (TH2D*)_file->Get("h_10");

        ttime_err[i] = (h_10->GetBinCenter(h_10->GetNbinsX()) - h_10->GetBinCenter(1))/2.0;
        ttime[i] = h_10->GetBinCenter(1) + ttime_err[i];

        removenoizypixelsXY(h_1,4);

        TH1D* h_x = h_1->ProjectionX("h_x");
        TH1D* h_y = h_1->ProjectionY("h_y");

        TF1* fit_x = new TF1("fit_x","gaus",400,450);
        TF1* fit_y = new TF1("fit_y","gaus",110,160);

        h_x->Fit(fit_x,"RQ+");
        h_y->Fit(fit_y,"RQ+");

        mean_x[i] = fit_x->GetParameter(1);
        mean_err_x[i] = fit_x->GetParError(1);
        mean_y[i] = fit_y->GetParameter(1);
        mean_err_y[i] = fit_y->GetParError(1);

        TString c_name_x = "c_x_";
        c_name_x += run_id;
        c_x[i] = new TCanvas(c_name_x.Data(),c_name_x.Data(),1000,1000);
        c_x[i]->cd();
        gPad->SetGrid();
        h_x->Draw();

        TString c_name_y = "c_y_";
        c_name_y += run_id;
        c_y[i] = new TCanvas(c_name_y.Data(),c_name_y.Data(),1000,1000);
        c_y[i]->cd();
        gPad->SetGrid();
        h_y->Draw();

        i++;
    }

    TCanvas* c_x_t = new TCanvas("c_x_t","c_x_t",1000,1000);
    c_x_t->cd();
    TGraphErrors* gr_x = new TGraphErrors(nRuns,ttime,mean_x,ttime_err,mean_err_x);
    gPad->SetGrid();
    gr_x->SetMarkerStyle(21);
    gr_x->SetMarkerColor(kRed);
    gr_x->SetLineStyle(9);
    gr_x->SetLineColor(kRed);
    gr_x->Draw("APL");

    TCanvas* c_y_t = new TCanvas("c_y_t","c_y_t",1000,1000);
    c_y_t->cd();
    TGraphErrors* gr_y = new TGraphErrors(nRuns,ttime,mean_y,ttime_err,mean_err_y);
    gPad->SetGrid();
    gr_y->SetMarkerStyle(21);
    gr_y->SetMarkerColor(kBlue);
    gr_y->SetLineStyle(9);
    gr_y->SetLineColor(kBlue);
    gr_y->Draw("APL");

    return 0;
}

bool checkpixel(TH1D *h, Int_t bin)// to check the neighbor pixels
{
    Double_t val_ref = h->GetBinContent(bin);
    Double_t val_ref_max = h->GetBinContent(bin) + h->GetBinError(bin);
    Double_t val_ref_min = h->GetBinContent(bin) - h->GetBinError(bin);

    Double_t val_left_max = h->GetBinContent(bin-1) + h->GetBinError(bin-1);
    Double_t val_left_min = h->GetBinContent(bin-1) - h->GetBinError(bin-1);

    Double_t val_right_max = h->GetBinContent(bin+1) + h->GetBinError(bin+1);
    Double_t val_right_min = h->GetBinContent(bin+1) - h->GetBinError(bin+1);

    if(bin > 1)
    {
        if(val_ref < val_left_max && val_ref > val_left_min) return true;
        if(val_ref_max < val_left_max && val_ref_max > val_left_min) return true;
        if(val_ref_min < val_left_max && val_ref_max > val_left_min) return true;
    }
    if(bin < h->GetNbinsX())
    {
        if(val_ref < val_right_max && val_ref > val_right_min) return true;
        if(val_ref_max < val_right_max && val_ref_max > val_right_min) return true;
        if(val_ref_min < val_right_max && val_ref_max > val_right_min) return true;
    }
    return false;
}

void removenoizypixelsX(TH1D* h_1)// to remove noizy pixels for X profile
{
    for(Int_t i = 1; i <= h_1->GetNbinsX(); i++)
    {
        if(i == 49 || i == 131 || i == 145 || i == 210 || i == 211 || i == 242
                || i == 245 || i == 362 || i == 372 || i == 386 || i == 436 || i == 473
                || i == 93 || i == 239 || i == 455 || i == 382 || i == 334 || i == 383
                || i == 402 || i == 411)
        {
            h_1->SetBinContent(i,0);
            h_1->SetBinError(i,0);
        }
    }
}

void removenoizypixelsY(TH1D* h_1)// to remove noizy pixels for Y profile
{
    for(Int_t i = 1; i <= h_1->GetNbinsX(); i++)
    {
        if(i == 63 || i == 74 || i == 100 || i == 111 || i == 125 || i == 211
                || i == 210 || i == 223 || i == 250 || i == 255 || i == 270 || i == 339
                || i == 425 || i == 239 || i == 220 || i == 216 || i == 240 || i == 167
                || i == 224 || i == 426)
        {
            h_1->SetBinContent(i,0);
            h_1->SetBinError(i,0);
        }
    }
}

void removenoizypixelsXY(TH2D* histo, Float_t factor_bad_pixels)// to remove noizy pixels for XY image
{
    for(Int_t xi = 1; xi < histo->GetNbinsX(); xi++)
    {
        for(Int_t yi = 1; yi < histo->GetNbinsX(); yi++)
        {
            if(histo->GetBinContent(xi,yi) - factor_bad_pixels*histo->GetBinError(xi,yi) >
                    histo->GetBinContent(xi,yi+1) + factor_bad_pixels*histo->GetBinError(xi,yi+1))
            {
                histo->SetBinContent(xi,yi,0);
                histo->SetBinError(xi,yi,0);
            }
        }
    }
}

int getMinPositionWithinRangeTH1(Int_t min_x, Int_t max_x, TH1* histo)
{
    Double_t min_val = 1e9;
    Int_t min_val_pos = -1;
    if(min_x < 1) min_x = 1;
    if(max_x > histo->GetNbinsX()) max_x = histo->GetNbinsX();
    for(Int_t i = min_x; i <= max_x; i++)
    {
        if(min_val > histo->GetBinContent(i))
        {
            min_val = histo->GetBinContent(i);
            min_val_pos = i;
        }
    }
    return min_val_pos;
}
