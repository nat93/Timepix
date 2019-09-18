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

void function_1();
void function_2();

int plot()
{
	cout<<endl<<"Functions:"<<endl;
        cout<<"--> function_1(): plot EXP data distributions"	<<endl;
        cout<<"--> function_2(): plot EXP vs SIM results"	<<endl;

	return 0;
}

void function_1()
{
    cout<<"### function_1() ###"<<endl;

    TString fileName_rp0 = "output_function_7_RP0.root";
    TString fileName_rp1 = "output_function_7_RP1.root";

    cout<<endl<<"--> input file RP0: "<<fileName_rp0<<endl;
    TFile* _file_rp0 = TFile::Open(fileName_rp0.Data());
    TH2D* h1_rp0 = (TH2D*)_file_rp0->Get("hh_1");
    TH2D* h2_rp0 = (TH2D*)_file_rp0->Get("hh_2");
    h1_rp0->SetName("h1_rp0");
    h2_rp0->SetName("h2_rp0");

    cout<<endl<<"--> input file RP1: "<<fileName_rp1<<endl;
    TFile* _file_rp1 = TFile::Open(fileName_rp1.Data());
    TH2D* h1_rp1 = (TH2D*)_file_rp1->Get("hh_1");
    TH2D* h2_rp1 = (TH2D*)_file_rp1->Get("hh_2");
    h1_rp1->SetName("h1_rp1");
    h2_rp1->SetName("h2_rp1");

    Int_t locmax,locmay,locmaz;
    h1_rp1->GetMaximumBin(locmax,locmay,locmaz);
    scale2Dhisto(h1_rp0, h1_rp1->GetBinContent(locmax,locmay),h1_rp1->GetBinError(locmax,locmay));
    scale2Dhisto(h1_rp1, h1_rp1->GetBinContent(locmax,locmay),h1_rp1->GetBinError(locmax,locmay));

    h2_rp1->GetMaximumBin(locmax,locmay,locmaz);
    scale2Dhisto(h2_rp0, h2_rp1->GetBinContent(locmax,locmay),h2_rp1->GetBinError(locmax,locmay));
    scale2Dhisto(h2_rp1, h2_rp1->GetBinContent(locmax,locmay),h2_rp1->GetBinError(locmax,locmay));

    TCanvas* c1 = new TCanvas("c1","c1",700,700);
    c1->cd();
    TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
    pad1->SetRightMargin(0.2);
    pad1->SetBottomMargin(0.2);
    pad1->Draw();
    pad1->cd();
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    h1_rp0->SetTitle("RP0 Timepix");
    h1_rp0->SetMaximum(1.0);
    h1_rp0->SetMinimum(0.0);
    h1_rp0->GetXaxis()->SetTitle("Horizontal Axis [mm]");
    h1_rp0->GetYaxis()->SetTitle("Vertical Axis [mm]");
    h1_rp0->GetZaxis()->SetTitle("Hits Density");
    h1_rp0->GetXaxis()->CenterTitle();
    h1_rp0->GetYaxis()->CenterTitle();
    h1_rp0->GetZaxis()->CenterTitle();
    h1_rp0->GetYaxis()->SetTitleOffset(0.7);
    h1_rp0->GetZaxis()->SetTitleOffset(1.0);
    h1_rp0->Draw("colz");

    TCanvas* c2 = new TCanvas("c2","c2",700,700);
    c2->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,1);
    pad2->SetRightMargin(0.2);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    h1_rp1->SetTitle("RP1 Timepix");
    h1_rp1->SetMaximum(1.0);
    h1_rp1->SetMinimum(0.0);
    h1_rp1->GetXaxis()->SetTitle("Horizontal Axis [mm]");
    h1_rp1->GetYaxis()->SetTitle("Vertical Axis [mm]");
    h1_rp1->GetZaxis()->SetTitle("Hits Density");
    h1_rp1->GetXaxis()->CenterTitle();
    h1_rp1->GetYaxis()->CenterTitle();
    h1_rp1->GetZaxis()->CenterTitle();
    h1_rp1->GetYaxis()->SetTitleOffset(0.7);
    h1_rp1->GetZaxis()->SetTitleOffset(1.0);
    h1_rp1->Draw("colz");

    TCanvas* c3 = new TCanvas("c3","c3",1000,500);
    c3->cd();
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    h2_rp0->SetTitle("RP0 Timepix");
    h2_rp0->SetMaximum(1.0);
    h2_rp0->SetMinimum(0.0);
    h2_rp0->GetXaxis()->SetTitle("Time [s]");
    h2_rp0->GetYaxis()->SetTitle("Horizontal Axis [mm]");
    h2_rp0->GetZaxis()->SetTitle("Hits Density");
    h2_rp0->GetXaxis()->CenterTitle();
    h2_rp0->GetYaxis()->CenterTitle();
    h2_rp0->GetZaxis()->CenterTitle();
    h2_rp0->Draw("colz");

    TCanvas* c4 = new TCanvas("c4","c4",1000,500);
    c4->cd();
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    h2_rp1->SetTitle("RP1 Timepix");
    h2_rp1->SetMaximum(1.0);
    h2_rp1->SetMinimum(0.0);
    h2_rp1->GetXaxis()->SetTitle("Time [s]");
    h2_rp1->GetYaxis()->SetTitle("Horizontal Axis [mm]");
    h2_rp1->GetZaxis()->SetTitle("Hits Density");
    h2_rp1->GetXaxis()->CenterTitle();
    h2_rp1->GetYaxis()->CenterTitle();
    h2_rp1->GetZaxis()->CenterTitle();
    h2_rp1->Draw("colz");

}

void function_2()
{
        cout<<"### function_2() ###"<<endl;

        TString fileName_exp = "output_function_7_RP1.root";
        TString fileName_sim1 = "output_crystal1_accelerator1_rp0_null_accelerator2_crystal2_null_accelerator3.root";
        TString fileName_sim2 = "output_crystal1_accelerator1_rp0_layers_accelerator2_crystal2_null_accelerator3.root";
        TString fileName_sim3 = "output_crystal1_accelerator1_rp0_simple_accelerator2_crystal2_null_accelerator3.root";

        cout<<endl<<"--> input file exp: "<<fileName_exp<<endl;
        TFile* _file_exp = TFile::Open(fileName_exp.Data());
        TH2D* h1_exp = (TH2D*)_file_exp->Get("hh_1");
        h1_exp->SetName("h1_exp");

        cout<<endl<<"--> input file sim1: "<<fileName_sim1<<endl;
        TFile* _file_sim1 = TFile::Open(fileName_sim1.Data());
        TH2D* h1_sim1 = (TH2D*)_file_sim1->Get("h_1");
        h1_sim1->SetName("h1_sim1");

        cout<<endl<<"--> input file sim2: "<<fileName_sim2<<endl;
        TFile* _file_sim2 = TFile::Open(fileName_sim2.Data());
        TH2D* h1_sim2 = (TH2D*)_file_sim2->Get("h_1");
        h1_sim2->SetName("h1_sim2");

        cout<<endl<<"--> input file sim3: "<<fileName_sim3<<endl;
        TFile* _file_sim3 = TFile::Open(fileName_sim3.Data());
        TH2D* h1_sim3 = (TH2D*)_file_sim3->Get("h_1");
        h1_sim3->SetName("h1_sim3");

        //====================//
        // Experimental
        //====================//

        TH1D* h2_exp = h1_exp->ProjectionX("h2_exp");
        TH1D* h2_tmp = h1_exp->ProjectionX("h2_tmp");

        h2_exp->Rebin(2);
        h2_tmp->Rebin(2);

        TF1* fit1 = new TF1("fit1","pol0",13,14);

        TCanvas* c1 = new TCanvas("c1","c1",1000,700);
        c1->cd();
        gStyle->SetOptStat(0);
        gPad->SetGrid();
        h2_tmp->SetTitle("");
        h2_tmp->SetLineColor(kRed);
        h2_tmp->Draw("HIST");
        h2_tmp->Fit(fit1,"R0Q");

        TF1* draw1 = new TF1("draw1","pol0",0,14);
        draw1->SetParameter(0,fit1->GetParameter(0));
        draw1->SetLineColor(kBlack);
        draw1->SetLineStyle(9);
        draw1->Draw("SAME L");

        for(Int_t i = 1; i <= h2_exp->GetNbinsX(); i++)
        {
            Double_t val = h2_tmp->GetBinContent(i) - fit1->GetParameter(0);
            if(val < 0) val = 0;
            h2_exp->SetBinContent(i,val);
        }
        h2_exp->SetTitle("");
        h2_exp->Draw("SAME HIST");

        //====================//
        // Simulation
        //====================//

        TH1D* h2_tmp21 = h1_sim1->ProjectionX("h2_tmp21");
        TH1D* h2_tmp22 = h1_sim2->ProjectionX("h2_tmp22");
        TH1D* h2_tmp23 = h1_sim3->ProjectionX("h2_tmp23");

        TH1D* h2_sim1 = new TH1D("h2_sim1","h2_sim1",h2_exp->GetNbinsX(),0,256*0.055);
        TH1D* h2_sim2 = new TH1D("h2_sim2","h2_sim2",h2_exp->GetNbinsX(),0,256*0.055);
        TH1D* h2_sim3 = new TH1D("h2_sim3","h2_sim3",h2_exp->GetNbinsX(),0,256*0.055);

        TH1D* h2_tmp3 = new TH1D("h2_tmp3","h2_tmp3",h2_exp->GetNbinsX(),0,256*0.055);

        TCanvas* c2 = new TCanvas("c2","c2",1000,700);
        c2->cd();
        gStyle->SetOptStat(0);
        gPad->SetGrid();
        h2_tmp21->SetTitle("");
        h2_tmp22->SetTitle("");
        h2_tmp23->SetTitle("");
        h2_tmp21->SetLineColor(40);
        h2_tmp22->SetLineColor(41);
        h2_tmp23->SetLineColor(42);
        h2_tmp21->Draw("");
        h2_tmp22->Draw("SAME");
        h2_tmp23->Draw("SAME");

        Double_t shiftX = /*MADX-OPTICS*/1.79604 + /*ALIGNMENT*/(36.76 - 33.76) + /*INTERNAL-OFFSET*/2.97 - 1.2; // [mm]

        for(Int_t i = 1; i <= h2_tmp21->GetNbinsX(); i++)
        {
            Double_t posX = h2_tmp21->GetBinCenter(i);
            posX = (-1e3)*posX - shiftX;

            if(posX >= 0)
            {
                h2_sim1->Fill(posX,h2_tmp21->GetBinContent(i));
                h2_sim2->Fill(posX,h2_tmp22->GetBinContent(i));
                h2_sim3->Fill(posX,h2_tmp23->GetBinContent(i));

                h2_tmp3->Fill(posX,1);
            }
        }
        for(Int_t i = 1; i <= h2_tmp3->GetNbinsX(); i++)
        {
            Double_t nFills = h2_tmp3->GetBinContent(i);

            if(nFills > 0)
            {
                h2_sim1->SetBinContent(i,h2_sim1->GetBinContent(i)/nFills);
                h2_sim2->SetBinContent(i,h2_sim2->GetBinContent(i)/nFills);
                h2_sim3->SetBinContent(i,h2_sim3->GetBinContent(i)/nFills);
            }
        }

        TCanvas* c3 = new TCanvas("c3","c3",1000,700);
        c3->cd();
        gStyle->SetOptStat(0);
        gPad->SetGrid();
        h2_sim1->SetTitle("");
        h2_sim2->SetTitle("");
        h2_sim3->SetTitle("");
        h2_sim1->Draw("");
        h2_sim2->Draw("SAME");
        h2_sim3->Draw("SAME");

        //====================//
        // Comparison
        //====================//

        cout<<"--> bin size : Exp. = "<<h2_exp->GetBinWidth(1)<<" [mm] ; Sim. = "<<h2_sim1->GetBinWidth(1)<<" [mm]"<<endl;

        Double_t integral, integralErr;

        integral = h2_exp->IntegralAndError(1,h2_exp->GetNbinsX(),integralErr);
        scale1Dhisto(h2_exp,integral,integralErr);

        integral = h2_sim1->IntegralAndError(1,h2_sim1->GetNbinsX(),integralErr);
        scale1Dhisto(h2_sim1,integral,integralErr);

        integral = h2_sim2->IntegralAndError(1,h2_sim2->GetNbinsX(),integralErr);
        scale1Dhisto(h2_sim2,integral,integralErr);

        integral = h2_sim3->IntegralAndError(1,h2_sim3->GetNbinsX(),integralErr);
        scale1Dhisto(h2_sim3,integral,integralErr);

        TCanvas* c4 = new TCanvas("c4","c4",1000,700);
        c4->cd();
        gStyle->SetOptStat(0);
        gPad->SetGrid();
        h2_sim1->SetLineColor(1);
        h2_sim2->SetLineColor(2);
        h2_sim3->SetLineColor(3);
        h2_exp->SetLineColor(4);
        h2_sim1->GetYaxis()->SetTitle("Normalized Counts");
        h2_sim1->GetXaxis()->SetTitle("Horizontal Axis Projection [mm]");
        h2_sim1->GetXaxis()->CenterTitle();
        h2_sim1->GetYaxis()->CenterTitle();
        h2_sim1->Draw("HIST");
        h2_sim2->Draw("HIST SAME");
        h2_sim3->Draw("HIST SAME");
        h2_exp->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.3,0.7,0.8,0.85);
        leg->AddEntry(h2_exp,"Exp. (Timepix, SPS-MD)","f");
        leg->AddEntry(h2_sim1,"Sim. (Geant4*): RP0 is OUT","f");
        leg->AddEntry(h2_sim2,"Sim. (Geant4*): RP0 is Timepix+Window","f");
        leg->AddEntry(h2_sim3,"Sim. (Geant4*): RP0 is 1mm of W","f");
        leg->Draw();

/**/
}

void scale2Dhisto(TH2D* histo, Double_t integral, Double_t interal_err)// to scale the 2D histo
{
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= histo->GetNbinsY(); j++)
        {
            Double_t val = histo->GetBinContent(i,j);
            Double_t val_err = histo->GetBinError(i,j);
            val_err = TMath::Sqrt(TMath::Power(val_err/integral,2) + TMath::Power(val*interal_err/(integral*integral),2));
            val = val/integral;

            histo->SetBinContent(i,j,val);
            histo->SetBinError(i,j,val_err);
        }
    }
}

void scale1Dhisto(TH1D* histo, Double_t integral, Double_t integral_err)// to scale the 1D histo
{
    for(Int_t i = 1; i <= histo->GetNbinsX(); i++)
    {
        Double_t val = histo->GetBinContent(i);
        Double_t val_err = histo->GetBinError(i);
        val_err = TMath::Sqrt(TMath::Power(val_err/integral,2) + TMath::Power(val*integral_err/(integral*integral),2));
        val = val/integral;

        histo->SetBinContent(i,val);
        histo->SetBinError(i,val_err);
    }
}
