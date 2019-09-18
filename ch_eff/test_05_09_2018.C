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

void Normalized_Array(Double_t norm, Double_t *array, Int_t nn)
{
	if(norm > 0)
	{
		for(Int_t i = 0; i < nn; i++)
		{
			array[i] = array[i]/norm;
		}
	}
	else
	{
		cout<<"ERROR:: The first input argument must be > 0!!!"<<endl;
		assert(0);
	}
}

int test_05_09_2018()
{
    //----------------------------------------------//
    //---------------- 05.09.2018 ------------------//
    //----------------------------------------------//
    gStyle->SetOptStat(0);
	
	const Int_t nTrials = 5;
	Double_t Bias_Verif_Err_5[] = {1.0, 1.0, 1.0, 1.0, 1.0};
	Double_t Bias_Verif_Err_6[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	Double_t Bias_Verif_Err_7[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	Double_t Bias_Verif_Err_9[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	
	//--- RP0I ---//
	const Int_t nBias_RP0I = 6;
	Double_t Bias_RP0I[] = {0, 10, 20, 40, 80, 100};
	Double_t Bias_Verif_RP0I[] = {4.2, 9.3, 18.8, 37.5, 75.2, 94.1};	
	Double_t Total_Counts_RP0I[nBias_RP0I][nTrials];
	Double_t Average_Counts_RP0I[nBias_RP0I];
	Double_t Average_Counts_Err_RP0I[nBias_RP0I];
	Double_t Max_Average_Counts_RP0I = 0.0;
	Total_Counts_RP0I[0][0] = 110; Total_Counts_RP0I[0][1] = 108; Total_Counts_RP0I[0][2] = 120; Total_Counts_RP0I[0][3] = 106; Total_Counts_RP0I[0][4] = 126;
	Total_Counts_RP0I[1][0] = 107; Total_Counts_RP0I[1][1] = 110; Total_Counts_RP0I[1][2] = 113; Total_Counts_RP0I[1][3] = 112; Total_Counts_RP0I[1][4] = 63;
	Total_Counts_RP0I[2][0] = 102; Total_Counts_RP0I[2][1] = 96; Total_Counts_RP0I[2][2] = 90; Total_Counts_RP0I[2][3] = 130; Total_Counts_RP0I[2][4] = 102;
	Total_Counts_RP0I[3][0] = 124; Total_Counts_RP0I[3][1] = 86; Total_Counts_RP0I[3][2] = 127; Total_Counts_RP0I[3][3] = 78; Total_Counts_RP0I[3][4] = 95;
	Total_Counts_RP0I[4][0] = 1105; Total_Counts_RP0I[4][1] = 1035; Total_Counts_RP0I[4][2] = 923; Total_Counts_RP0I[4][3] = 1061; Total_Counts_RP0I[4][4] = 801;
	Total_Counts_RP0I[5][0] = 1413; Total_Counts_RP0I[5][1] = 1268; Total_Counts_RP0I[5][2] = 1234; Total_Counts_RP0I[5][3] = 1266; Total_Counts_RP0I[5][4] = 1339;
	
	for(Int_t i = 0; i < nBias_RP0I; i++)
	{
		Average_Counts_RP0I[i] = 0.0;
		Average_Counts_Err_RP0I[i] = 0.0;
		for(Int_t j = 0; j < nTrials; j++)
		{
			Average_Counts_RP0I[i] += Total_Counts_RP0I[i][j];			
		}
		Average_Counts_Err_RP0I[i] = TMath::Sqrt(Average_Counts_RP0I[i])/nTrials;
		Average_Counts_RP0I[i] = Average_Counts_RP0I[i]/nTrials;
		
		if(Max_Average_Counts_RP0I < Average_Counts_RP0I[i]) {Max_Average_Counts_RP0I = Average_Counts_RP0I[i];}
	}
	Normalized_Array(Max_Average_Counts_RP0I,Average_Counts_RP0I,nBias_RP0I);
	Normalized_Array(Max_Average_Counts_RP0I,Average_Counts_Err_RP0I,nBias_RP0I);
	
	TGraphErrors* gr_bias_diff_RP0I = new TGraphErrors(nBias_RP0I,Bias_RP0I,Bias_Verif_RP0I,0,Bias_Verif_Err_6);
	gr_bias_diff_RP0I->SetTitle("Verification Bias vs Applied Bias");
	gr_bias_diff_RP0I->SetName("gr_bias_diff_RP0I");
	gr_bias_diff_RP0I->SetLineWidth(2);
	gr_bias_diff_RP0I->SetMarkerStyle(20);
	gr_bias_diff_RP0I->SetLineColor(kRed);
	gr_bias_diff_RP0I->SetLineStyle(1);
	gr_bias_diff_RP0I->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_bias_diff_RP0I->GetYaxis()->SetTitle("Verification Bias [V]");
	
	TGraphErrors* gr_counts_app_bias_RP0I = new TGraphErrors(nBias_RP0I,Bias_RP0I,Average_Counts_RP0I,0,Average_Counts_Err_RP0I);
	gr_counts_app_bias_RP0I->SetTitle("Average Counts vs Applied Bias");
	gr_counts_app_bias_RP0I->SetName("gr_counts_app_bias_RP0I");
	gr_counts_app_bias_RP0I->SetLineWidth(2);
	gr_counts_app_bias_RP0I->SetMarkerStyle(20);
	gr_counts_app_bias_RP0I->SetLineColor(kRed);
	gr_counts_app_bias_RP0I->SetLineStyle(1);
	gr_counts_app_bias_RP0I->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_counts_app_bias_RP0I->GetYaxis()->SetTitle("Average Counts");
	
	TGraphErrors* gr_counts_ver_bias_RP0I = new TGraphErrors(nBias_RP0I,Bias_Verif_RP0I,Average_Counts_RP0I,Bias_Verif_Err_6,Average_Counts_Err_RP0I);
	gr_counts_ver_bias_RP0I->SetTitle("Average Counts vs Verification Bias");
	gr_counts_ver_bias_RP0I->SetName("gr_counts_ver_bias_RP0I");
	gr_counts_ver_bias_RP0I->SetLineWidth(2);
	gr_counts_ver_bias_RP0I->SetMarkerStyle(20);
	gr_counts_ver_bias_RP0I->SetLineColor(kRed);
	gr_counts_ver_bias_RP0I->SetLineStyle(1);
	gr_counts_ver_bias_RP0I->GetXaxis()->SetTitle("Verification Bias [V]");
	gr_counts_ver_bias_RP0I->GetYaxis()->SetTitle("Average Counts");
	
	//--- RP0E ---//
	const Int_t nBias_RP0E = 5;
	Double_t Bias_RP0E[] = {0, 10, 20, 40, 80};
	Double_t Bias_Verif_RP0E[] = {4, 9.1, 18.6, 37.7, 75.8};
	Double_t Total_Counts_RP0E[nBias_RP0E][nTrials];
	Double_t Average_Counts_RP0E[nBias_RP0E];
	Double_t Average_Counts_Err_RP0E[nBias_RP0E];
	Double_t Max_Average_Counts_RP0E = 0.0;
	Total_Counts_RP0E[0][0] = 691; Total_Counts_RP0E[0][1] = 586; Total_Counts_RP0E[0][2] = 610; Total_Counts_RP0E[0][3] = 625; Total_Counts_RP0E[0][4] = 556;
	Total_Counts_RP0E[1][0] = 849; Total_Counts_RP0E[1][1] = 1106; Total_Counts_RP0E[1][2] = 811; Total_Counts_RP0E[1][3] = 865; Total_Counts_RP0E[1][4] = 897;
	Total_Counts_RP0E[2][0] = 1286; Total_Counts_RP0E[2][1] = 1292; Total_Counts_RP0E[2][2] = 1347; Total_Counts_RP0E[2][3] = 1458; Total_Counts_RP0E[2][4] = 1425;
	Total_Counts_RP0E[3][0] = 1527; Total_Counts_RP0E[3][1] = 1770; Total_Counts_RP0E[3][2] = 1831; Total_Counts_RP0E[3][3] = 1488; Total_Counts_RP0E[3][4] = 1438;
	Total_Counts_RP0E[4][0] = 1470; Total_Counts_RP0E[4][1] = 1485; Total_Counts_RP0E[4][2] = 1751; Total_Counts_RP0E[4][3] = 1726; Total_Counts_RP0E[4][4] = 1615;
	
	for(Int_t i = 0; i < nBias_RP0E; i++)
	{
		Average_Counts_RP0E[i] = 0.0;
		Average_Counts_Err_RP0E[i] = 0.0;
		for(Int_t j = 0; j < nTrials; j++)
		{
			Average_Counts_RP0E[i] += Total_Counts_RP0E[i][j];			
		}
		Average_Counts_Err_RP0E[i] = TMath::Sqrt(Average_Counts_RP0E[i])/nTrials;
		Average_Counts_RP0E[i] = Average_Counts_RP0E[i]/nTrials;
		
		if(Max_Average_Counts_RP0E < Average_Counts_RP0E[i]) {Max_Average_Counts_RP0E = Average_Counts_RP0E[i];}
	}
	Normalized_Array(Max_Average_Counts_RP0E,Average_Counts_RP0E,nBias_RP0E);
	Normalized_Array(Max_Average_Counts_RP0E,Average_Counts_Err_RP0E,nBias_RP0E);
	
	TGraphErrors* gr_bias_diff_RP0E = new TGraphErrors(nBias_RP0E,Bias_RP0E,Bias_Verif_RP0E,0,Bias_Verif_Err_5);
	gr_bias_diff_RP0E->SetTitle("Verification Bias vs Applied Bias");
	gr_bias_diff_RP0E->SetName("gr_bias_diff_RP0E");
	gr_bias_diff_RP0E->SetLineWidth(2);
	gr_bias_diff_RP0E->SetMarkerStyle(21);
	gr_bias_diff_RP0E->SetLineColor(kGreen);
	gr_bias_diff_RP0E->SetLineStyle(2);	
	gr_bias_diff_RP0E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_bias_diff_RP0E->GetYaxis()->SetTitle("Verification Bias [V]");
	
	TGraphErrors* gr_counts_app_bias_RP0E = new TGraphErrors(nBias_RP0E,Bias_RP0E,Average_Counts_RP0E,0,Average_Counts_Err_RP0E);
	gr_counts_app_bias_RP0E->SetTitle("Average Counts vs Applied Bias");
	gr_counts_app_bias_RP0E->SetName("gr_counts_app_bias_RP0E");
	gr_counts_app_bias_RP0E->SetLineWidth(2);
	gr_counts_app_bias_RP0E->SetMarkerStyle(21);
	gr_counts_app_bias_RP0E->SetLineColor(kGreen);
	gr_counts_app_bias_RP0E->SetLineStyle(2);
	gr_counts_app_bias_RP0E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_counts_app_bias_RP0E->GetYaxis()->SetTitle("Average Counts");
	
	TGraphErrors* gr_counts_ver_bias_RP0E = new TGraphErrors(nBias_RP0E,Bias_Verif_RP0E,Average_Counts_RP0E,Bias_Verif_Err_5,Average_Counts_Err_RP0E);
	gr_counts_ver_bias_RP0E->SetTitle("Average Counts vs Verification Bias");
	gr_counts_ver_bias_RP0E->SetName("gr_counts_ver_bias_RP0E");
	gr_counts_ver_bias_RP0E->SetLineWidth(2);
	gr_counts_ver_bias_RP0E->SetMarkerStyle(21);
	gr_counts_ver_bias_RP0E->SetLineColor(kGreen);
	gr_counts_ver_bias_RP0E->SetLineStyle(2);
	gr_counts_ver_bias_RP0E->GetXaxis()->SetTitle("Verification Bias [V]");
	gr_counts_ver_bias_RP0E->GetYaxis()->SetTitle("Average Counts");
	
	//--- RP1I ---//
	const Int_t nBias_RP1I = 5;
	Double_t Bias_RP1I[] = {0, 10, 20, 40, 80};
	Double_t Bias_Verif_RP1I[] = {5, 8.9, 18.5, 37.4, 75.5};
	Double_t Total_Counts_RP1I[nBias_RP1I][nTrials];
	Double_t Average_Counts_RP1I[nBias_RP1I];
	Double_t Average_Counts_Err_RP1I[nBias_RP1I];
	Double_t Max_Average_Counts_RP1I = 0.0;
	Total_Counts_RP1I[0][0] = 190; Total_Counts_RP1I[0][1] = 177; Total_Counts_RP1I[0][2] = 236; Total_Counts_RP1I[0][3] = 152; Total_Counts_RP1I[0][4] = 158;
	Total_Counts_RP1I[1][0] = 192; Total_Counts_RP1I[1][1] = 178; Total_Counts_RP1I[1][2] = 147; Total_Counts_RP1I[1][3] = 194; Total_Counts_RP1I[1][4] = 200;
	Total_Counts_RP1I[2][0] = 174; Total_Counts_RP1I[2][1] = 168; Total_Counts_RP1I[2][2] = 186; Total_Counts_RP1I[2][3] = 159; Total_Counts_RP1I[2][4] = 176;
	Total_Counts_RP1I[3][0] = 191; Total_Counts_RP1I[3][1] = 180; Total_Counts_RP1I[3][2] = 161; Total_Counts_RP1I[3][3] = 189; Total_Counts_RP1I[3][4] = 172;
	Total_Counts_RP1I[4][0] = 175; Total_Counts_RP1I[4][1] = 192; Total_Counts_RP1I[4][2] = 198; Total_Counts_RP1I[4][3] = 196; Total_Counts_RP1I[4][4] = 166;
	
	for(Int_t i = 0; i < nBias_RP1I; i++)
	{
		Average_Counts_RP1I[i] = 0.0;
		Average_Counts_Err_RP1I[i] = 0.0;
		for(Int_t j = 0; j < nTrials; j++)
		{
			Average_Counts_RP1I[i] += Total_Counts_RP1I[i][j];			
		}
		Average_Counts_Err_RP1I[i] = TMath::Sqrt(Average_Counts_RP1I[i])/nTrials;
		Average_Counts_RP1I[i] = Average_Counts_RP1I[i]/nTrials;
		
		if(Max_Average_Counts_RP1I < Average_Counts_RP1I[i]) {Max_Average_Counts_RP1I = Average_Counts_RP1I[i];}
	}
	Normalized_Array(Max_Average_Counts_RP1I,Average_Counts_RP1I,nBias_RP1I);
	Normalized_Array(Max_Average_Counts_RP1I,Average_Counts_Err_RP1I,nBias_RP1I);
	
	TGraphErrors* gr_bias_diff_RP1I = new TGraphErrors(nBias_RP1I,Bias_RP1I,Bias_Verif_RP1I,0,Bias_Verif_Err_5);
	gr_bias_diff_RP1I->SetTitle("Verification Bias vs Applied Bias");
	gr_bias_diff_RP1I->SetName("gr_bias_diff_RP1I");
	gr_bias_diff_RP1I->SetLineWidth(2);
	gr_bias_diff_RP1I->SetMarkerStyle(22);
	gr_bias_diff_RP1I->SetLineColor(kBlue);
	gr_bias_diff_RP1I->SetLineStyle(9);
	gr_bias_diff_RP1I->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_bias_diff_RP1I->GetYaxis()->SetTitle("Verification Bias [V]");
	
	TGraphErrors* gr_counts_app_bias_RP1I = new TGraphErrors(nBias_RP1I,Bias_RP1I,Average_Counts_RP1I,0,Average_Counts_Err_RP1I);
	gr_counts_app_bias_RP1I->SetTitle("Average Counts vs Applied Bias");
	gr_counts_app_bias_RP1I->SetName("gr_counts_app_bias_RP1I");
	gr_counts_app_bias_RP1I->SetLineWidth(2);
	gr_counts_app_bias_RP1I->SetMarkerStyle(22);
	gr_counts_app_bias_RP1I->SetLineColor(kBlue);
	gr_counts_app_bias_RP1I->SetLineStyle(9);
	gr_counts_app_bias_RP1I->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_counts_app_bias_RP1I->GetYaxis()->SetTitle("Average Counts");
	
	TGraphErrors* gr_counts_ver_bias_RP1I = new TGraphErrors(nBias_RP1I,Bias_Verif_RP1I,Average_Counts_RP1I,Bias_Verif_Err_5,Average_Counts_Err_RP1I);
	gr_counts_ver_bias_RP1I->SetTitle("Average Counts vs Verification Bias");
	gr_counts_ver_bias_RP1I->SetName("gr_counts_ver_bias_RP1I");
	gr_counts_ver_bias_RP1I->SetLineWidth(2);
	gr_counts_ver_bias_RP1I->SetMarkerStyle(22);
	gr_counts_ver_bias_RP1I->SetLineColor(kBlue);
	gr_counts_ver_bias_RP1I->SetLineStyle(9);
	gr_counts_ver_bias_RP1I->GetXaxis()->SetTitle("Verification Bias [V]");
	gr_counts_ver_bias_RP1I->GetYaxis()->SetTitle("Average Counts");
	
	//--- RP1E ---//
	const Int_t nBias_RP1E = 9;
	Double_t Bias_RP1E[] = {0, 10, 20, 40, 60, 80, 100, 120, 140};
	Double_t Bias_Verif_RP1E[] = {-4.7, 9.8, 19.9, 39.8, 59.8, 79.8, 99.7, 119.9, 139.7};
	Double_t Total_Counts_RP1E[nBias_RP1E][nTrials];
	Double_t Average_Counts_RP1E[nBias_RP1E];
	Double_t Average_Counts_Err_RP1E[nBias_RP1E];
	Double_t Max_Average_Counts_RP1E = 0.0;
	Total_Counts_RP1E[0][0] = 1584; Total_Counts_RP1E[0][1] = 1647; Total_Counts_RP1E[0][2] = 1624; Total_Counts_RP1E[0][3] = 1621; Total_Counts_RP1E[0][4] = 1558;
	Total_Counts_RP1E[1][0] = 1208; Total_Counts_RP1E[1][1] = 1215; Total_Counts_RP1E[1][2] = 1213; Total_Counts_RP1E[1][3] = 1203; Total_Counts_RP1E[1][4] = 1172;
	Total_Counts_RP1E[2][0] = 1234; Total_Counts_RP1E[2][1] = 1311; Total_Counts_RP1E[2][2] = 1413; Total_Counts_RP1E[2][3] = 1183; Total_Counts_RP1E[2][4] = 1291;
	Total_Counts_RP1E[3][0] = 1589; Total_Counts_RP1E[3][1] = 1743; Total_Counts_RP1E[3][2] = 1623; Total_Counts_RP1E[3][3] = 1635; Total_Counts_RP1E[3][4] = 1673;
	Total_Counts_RP1E[4][0] = 1737; Total_Counts_RP1E[4][1] = 1883; Total_Counts_RP1E[4][2] = 1917; Total_Counts_RP1E[4][3] = 1945; Total_Counts_RP1E[4][4] = 1885;
	Total_Counts_RP1E[5][0] = 2094; Total_Counts_RP1E[5][1] = 2052; Total_Counts_RP1E[5][2] = 2094; Total_Counts_RP1E[5][3] = 2234; Total_Counts_RP1E[5][4] = 1826;
	Total_Counts_RP1E[6][0] = 2352; Total_Counts_RP1E[6][1] = 2306; Total_Counts_RP1E[6][2] = 2254; Total_Counts_RP1E[6][3] = 2121; Total_Counts_RP1E[6][4] = 2237;
	Total_Counts_RP1E[7][0] = 2193; Total_Counts_RP1E[7][1] = 2416; Total_Counts_RP1E[7][2] = 2434; Total_Counts_RP1E[7][3] = 2386; Total_Counts_RP1E[7][4] = 2549;
	Total_Counts_RP1E[8][0] = 2449; Total_Counts_RP1E[8][1] = 2617; Total_Counts_RP1E[8][2] = 2691; Total_Counts_RP1E[8][3] = 2232; Total_Counts_RP1E[8][4] = 2395;
	
	for(Int_t i = 0; i < nBias_RP1E; i++)
	{
		Average_Counts_RP1E[i] = 0.0;
		Average_Counts_Err_RP1E[i] = 0.0;
		for(Int_t j = 0; j < nTrials; j++)
		{
			Average_Counts_RP1E[i] += Total_Counts_RP1E[i][j];			
		}
		Average_Counts_Err_RP1E[i] = TMath::Sqrt(Average_Counts_RP1E[i])/nTrials;
		Average_Counts_RP1E[i] = Average_Counts_RP1E[i]/nTrials;
		
		if(Max_Average_Counts_RP1E < Average_Counts_RP1E[i]) {Max_Average_Counts_RP1E = Average_Counts_RP1E[i];}
	}
	Normalized_Array(Max_Average_Counts_RP1E,Average_Counts_RP1E,nBias_RP1E);
	Normalized_Array(Max_Average_Counts_RP1E,Average_Counts_Err_RP1E,nBias_RP1E);
	
	TGraphErrors* gr_bias_diff_RP1E = new TGraphErrors(nBias_RP1E,Bias_RP1E,Bias_Verif_RP1E,0,Bias_Verif_Err_9);
	gr_bias_diff_RP1E->SetTitle("Verification Bias vs Applied Bias");
	gr_bias_diff_RP1E->SetName("gr_bias_diff_RP1E");
	gr_bias_diff_RP1E->SetLineWidth(2);
	gr_bias_diff_RP1E->SetMarkerStyle(24);
	gr_bias_diff_RP1E->SetLineColor(kMagenta);
	gr_bias_diff_RP1E->SetLineStyle(5);
	gr_bias_diff_RP1E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_bias_diff_RP1E->GetYaxis()->SetTitle("Verification Bias [V]");
	
	TGraphErrors* gr_counts_app_bias_RP1E = new TGraphErrors(nBias_RP1E,Bias_RP1E,Average_Counts_RP1E,0,Average_Counts_Err_RP1E);
	gr_counts_app_bias_RP1E->SetTitle("Average Counts vs Applied Bias");
	gr_counts_app_bias_RP1E->SetName("gr_counts_app_bias_RP1E");
	gr_counts_app_bias_RP1E->SetLineWidth(2);
	gr_counts_app_bias_RP1E->SetMarkerStyle(24);
	gr_counts_app_bias_RP1E->SetLineColor(kMagenta);
	gr_counts_app_bias_RP1E->SetLineStyle(5);
	gr_counts_app_bias_RP1E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_counts_app_bias_RP1E->GetYaxis()->SetTitle("Average Counts");
	
	TGraphErrors* gr_counts_ver_bias_RP1E = new TGraphErrors(nBias_RP1E,Bias_Verif_RP1E,Average_Counts_RP1E,Bias_Verif_Err_9,Average_Counts_Err_RP1E);
	gr_counts_ver_bias_RP1E->SetTitle("Average Counts vs Verification Bias");
	gr_counts_ver_bias_RP1E->SetName("gr_counts_ver_bias_RP1E");
	gr_counts_ver_bias_RP1E->SetLineWidth(2);
	gr_counts_ver_bias_RP1E->SetMarkerStyle(24);
	gr_counts_ver_bias_RP1E->SetLineColor(kMagenta);
	gr_counts_ver_bias_RP1E->SetLineStyle(5);
	gr_counts_ver_bias_RP1E->GetXaxis()->SetTitle("Verification Bias [V]");
	gr_counts_ver_bias_RP1E->GetYaxis()->SetTitle("Average Counts");
	
	//--- RP3E ---//
	/*const Int_t nBias_RP3E = 6;
	Double_t Bias_RP3E[] = {0, 10, 20, 40, 80, 100};
	Double_t Bias_Verif_RP3E[] = {2.3, 5.3, 11.3, 26.6, 63.0, 82.4};
	Double_t Total_Counts_RP3E[nBias_RP3E][nTrials];
	Double_t Average_Counts_RP3E[nBias_RP3E];
	Double_t Average_Counts_Err_RP3E[nBias_RP3E];
	Double_t Max_Average_Counts_RP3E = 0.0;
	Total_Counts_RP3E[0][0] = 69; Total_Counts_RP3E[0][1] = 70; Total_Counts_RP3E[0][2] = 57; Total_Counts_RP3E[0][3] = 59; Total_Counts_RP3E[0][4] = 96;
	Total_Counts_RP3E[1][0] = 105; Total_Counts_RP3E[1][1] = 158; Total_Counts_RP3E[1][2] = 84; Total_Counts_RP3E[1][3] = 94; Total_Counts_RP3E[1][4] = 81;
	Total_Counts_RP3E[2][0] = 523; Total_Counts_RP3E[2][1] = 561; Total_Counts_RP3E[2][2] = 468; Total_Counts_RP3E[2][3] = 334; Total_Counts_RP3E[2][4] = 443;
	Total_Counts_RP3E[3][0] = 1021; Total_Counts_RP3E[3][1] = 898; Total_Counts_RP3E[3][2] = 806; Total_Counts_RP3E[3][3] = 883; Total_Counts_RP3E[3][4] = 1061;
	Total_Counts_RP3E[4][0] = 1403; Total_Counts_RP3E[4][1] = 1171; Total_Counts_RP3E[4][2] = 1273; Total_Counts_RP3E[4][3] = 1252; Total_Counts_RP3E[4][4] = 1019;
	Total_Counts_RP3E[5][0] = 1391; Total_Counts_RP3E[5][1] = 1259; Total_Counts_RP3E[5][2] = 1154; Total_Counts_RP3E[5][3] = 1232; Total_Counts_RP3E[5][4] = 1005;
	
	for(Int_t i = 0; i < nBias_RP3E; i++)
	{
		Average_Counts_RP3E[i] = 0.0;
		Average_Counts_Err_RP3E[i] = 0.0;
		for(Int_t j = 0; j < nTrials; j++)
		{
			Average_Counts_RP3E[i] += Total_Counts_RP3E[i][j];			
		}
		Average_Counts_Err_RP3E[i] = TMath::Sqrt(Average_Counts_RP3E[i])/nTrials;
		Average_Counts_RP3E[i] = Average_Counts_RP3E[i]/nTrials;
		
		if(Max_Average_Counts_RP3E < Average_Counts_RP3E[i]) {Max_Average_Counts_RP3E = Average_Counts_RP3E[i];}
	}
	Normalized_Array(Max_Average_Counts_RP3E,Average_Counts_RP3E,nBias_RP3E);
	Normalized_Array(Max_Average_Counts_RP3E,Average_Counts_Err_RP3E,nBias_RP3E);
	
	TGraphErrors* gr_bias_diff_RP3E = new TGraphErrors(nBias_RP3E,Bias_RP3E,Bias_Verif_RP3E,0,Bias_Verif_Err_6);
	gr_bias_diff_RP3E->SetTitle("Verification Bias vs Applied Bias");
	gr_bias_diff_RP3E->SetName("gr_bias_diff_RP3E");
	gr_bias_diff_RP3E->SetLineWidth(2);
	gr_bias_diff_RP3E->SetMarkerStyle(23);
	gr_bias_diff_RP3E->SetLineColor(kBlack);
	gr_bias_diff_RP3E->SetLineStyle(10);
	gr_bias_diff_RP3E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_bias_diff_RP3E->GetYaxis()->SetTitle("Verification Bias [V]");
	
	TGraphErrors* gr_counts_app_bias_RP3E = new TGraphErrors(nBias_RP3E,Bias_RP3E,Average_Counts_RP3E,0,Average_Counts_Err_RP3E);
	gr_counts_app_bias_RP3E->SetTitle("Average Counts vs Applied Bias");
	gr_counts_app_bias_RP3E->SetName("gr_counts_app_bias_RP3E");
	gr_counts_app_bias_RP3E->SetLineWidth(2);
	gr_counts_app_bias_RP3E->SetMarkerStyle(23);
	gr_counts_app_bias_RP3E->SetLineColor(kBlack);
	gr_counts_app_bias_RP3E->SetLineStyle(10);
	gr_counts_app_bias_RP3E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_counts_app_bias_RP3E->GetYaxis()->SetTitle("Average Counts");
	
	TGraphErrors* gr_counts_ver_bias_RP3E = new TGraphErrors(nBias_RP3E,Bias_Verif_RP3E,Average_Counts_RP3E,Bias_Verif_Err_6,Average_Counts_Err_RP3E);
	gr_counts_ver_bias_RP3E->SetTitle("Average Counts vs Verification Bias");
	gr_counts_ver_bias_RP3E->SetName("gr_counts_ver_bias_RP3E");
	gr_counts_ver_bias_RP3E->SetLineWidth(2);
	gr_counts_ver_bias_RP3E->SetMarkerStyle(23);
	gr_counts_ver_bias_RP3E->SetLineColor(kBlack);
	gr_counts_ver_bias_RP3E->SetLineStyle(10);
	gr_counts_ver_bias_RP3E->GetXaxis()->SetTitle("Verification Bias [V]");
	gr_counts_ver_bias_RP3E->GetYaxis()->SetTitle("Average Counts");*/
	
	//--- RP3E ---//
	const Int_t nBias_RP3E = 7;
	Double_t Bias_RP3E[] = {0, 10, 20, 40, 60, 80, 100};
	Double_t Bias_Verif_RP3E[] = {2.3, 5.8, 13, 29.8, 48.1, 66.5, 85};
	Double_t Total_Counts_RP3E[nBias_RP3E][nTrials];
	Double_t Average_Counts_RP3E[nBias_RP3E];
	Double_t Average_Counts_Err_RP3E[nBias_RP3E];
	Double_t Max_Average_Counts_RP3E = 0.0;
	Total_Counts_RP3E[0][0] =  23; Total_Counts_RP3E[0][1] =  40; Total_Counts_RP3E[0][2] =  33; Total_Counts_RP3E[0][3] =  37; Total_Counts_RP3E[0][4] =  42;
	Total_Counts_RP3E[1][0] =  99; Total_Counts_RP3E[1][1] = 109; Total_Counts_RP3E[1][2] =  80; Total_Counts_RP3E[1][3] =  76; Total_Counts_RP3E[1][4] =  78;
	Total_Counts_RP3E[2][0] = 243; Total_Counts_RP3E[2][1] = 217; Total_Counts_RP3E[2][2] = 152; Total_Counts_RP3E[2][3] = 197; Total_Counts_RP3E[2][4] = 251;
	Total_Counts_RP3E[3][0] = 444; Total_Counts_RP3E[3][1] = 417; Total_Counts_RP3E[3][2] = 359; Total_Counts_RP3E[3][3] = 398; Total_Counts_RP3E[3][4] = 414;
	Total_Counts_RP3E[4][0] = 301; Total_Counts_RP3E[4][1] = 541; Total_Counts_RP3E[4][2] = 375; Total_Counts_RP3E[4][3] = 397; Total_Counts_RP3E[4][4] = 376;
	Total_Counts_RP3E[5][0] = 400; Total_Counts_RP3E[5][1] = 342; Total_Counts_RP3E[5][2] = 425; Total_Counts_RP3E[5][3] = 353; Total_Counts_RP3E[5][4] = 416;
	Total_Counts_RP3E[6][0] = 337; Total_Counts_RP3E[6][1] = 443; Total_Counts_RP3E[6][2] = 403; Total_Counts_RP3E[6][3] = 328; Total_Counts_RP3E[6][4] = 449;
	
	for(Int_t i = 0; i < nBias_RP3E; i++)
	{
		Average_Counts_RP3E[i] = 0.0;
		Average_Counts_Err_RP3E[i] = 0.0;
		for(Int_t j = 0; j < nTrials; j++)
		{
			Average_Counts_RP3E[i] += Total_Counts_RP3E[i][j];			
		}
		Average_Counts_Err_RP3E[i] = TMath::Sqrt(Average_Counts_RP3E[i])/nTrials;
		Average_Counts_RP3E[i] = Average_Counts_RP3E[i]/nTrials;
		
		if(Max_Average_Counts_RP3E < Average_Counts_RP3E[i]) {Max_Average_Counts_RP3E = Average_Counts_RP3E[i];}
	}
	Normalized_Array(Max_Average_Counts_RP3E,Average_Counts_RP3E,nBias_RP3E);
	Normalized_Array(Max_Average_Counts_RP3E,Average_Counts_Err_RP3E,nBias_RP3E);
	
	TGraphErrors* gr_bias_diff_RP3E = new TGraphErrors(nBias_RP3E,Bias_RP3E,Bias_Verif_RP3E,0,Bias_Verif_Err_7);
	gr_bias_diff_RP3E->SetTitle("Verification Bias vs Applied Bias");
	gr_bias_diff_RP3E->SetName("gr_bias_diff_RP3E");
	gr_bias_diff_RP3E->SetLineWidth(2);
	gr_bias_diff_RP3E->SetMarkerStyle(23);
	gr_bias_diff_RP3E->SetLineColor(kBlack);
	gr_bias_diff_RP3E->SetLineStyle(10);
	gr_bias_diff_RP3E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_bias_diff_RP3E->GetYaxis()->SetTitle("Verification Bias [V]");
	
	TGraphErrors* gr_counts_app_bias_RP3E = new TGraphErrors(nBias_RP3E,Bias_RP3E,Average_Counts_RP3E,0,Average_Counts_Err_RP3E);
	gr_counts_app_bias_RP3E->SetTitle("Average Counts vs Applied Bias");
	gr_counts_app_bias_RP3E->SetName("gr_counts_app_bias_RP3E");
	gr_counts_app_bias_RP3E->SetLineWidth(2);
	gr_counts_app_bias_RP3E->SetMarkerStyle(23);
	gr_counts_app_bias_RP3E->SetLineColor(kBlack);
	gr_counts_app_bias_RP3E->SetLineStyle(10);
	gr_counts_app_bias_RP3E->GetXaxis()->SetTitle("Applied Bias [V]");
	gr_counts_app_bias_RP3E->GetYaxis()->SetTitle("Average Counts");
	
	TGraphErrors* gr_counts_ver_bias_RP3E = new TGraphErrors(nBias_RP3E,Bias_Verif_RP3E,Average_Counts_RP3E,Bias_Verif_Err_7,Average_Counts_Err_RP3E);
	gr_counts_ver_bias_RP3E->SetTitle("Average Counts vs Verification Bias");
	gr_counts_ver_bias_RP3E->SetName("gr_counts_ver_bias_RP3E");
	gr_counts_ver_bias_RP3E->SetLineWidth(2);
	gr_counts_ver_bias_RP3E->SetMarkerStyle(23);
	gr_counts_ver_bias_RP3E->SetLineColor(kBlack);
	gr_counts_ver_bias_RP3E->SetLineStyle(10);
	gr_counts_ver_bias_RP3E->GetXaxis()->SetTitle("Verification Bias [V]");
	gr_counts_ver_bias_RP3E->GetYaxis()->SetTitle("Average Counts");

	//--- PLOT ---//
	TCanvas* c1 = new TCanvas("c1","Bias Diff",1000,1000);
	c1->cd();
	gPad->SetGrid();
	TMultiGraph* mg1 = new TMultiGraph();
	mg1->Add(gr_bias_diff_RP0I);
	mg1->Add(gr_bias_diff_RP0E);
	mg1->Add(gr_bias_diff_RP1I);
	mg1->Add(gr_bias_diff_RP1E);
	mg1->Add(gr_bias_diff_RP3E);
	mg1->Draw("APL");	
	mg1->GetXaxis()->SetTitle("Applied Bias [V]");
	mg1->GetYaxis()->SetTitle("Verification Bias [V]");	
	mg1->GetXaxis()->SetLimits(-10,150);
	mg1->SetMinimum(-10);
	mg1->SetMaximum(150);
	auto legend1 = new TLegend(0.1,0.7,0.48,0.9);
	legend1->SetHeader("Timepix","C"); // option "C" allows to center the header
	legend1->AddEntry(gr_bias_diff_RP0I,"RP0 Internal (I02-W0108 FITPix-0415)","PLE");
	legend1->AddEntry(gr_bias_diff_RP0E,"RP0 External (G02-W0108 FITPix-0384)","PLE");
	legend1->AddEntry(gr_bias_diff_RP1I,"RP1 Internal (F04-W0108 FITPix-0409)","PLE");
	legend1->AddEntry(gr_bias_diff_RP1E,"RP1 External (H7-W00018 Katherine)","PLE");
	legend1->AddEntry(gr_bias_diff_RP3E,"RP3 External (C08-W0255 FITPix-0399)","PLE");
	legend1->Draw();		
	
	TCanvas* c2 = new TCanvas("c2","Average Counts vs Applied Bias",1500,1000);
	c2->cd();
	gPad->SetGrid();
	TMultiGraph* mg2 = new TMultiGraph();
	mg2->Add(gr_counts_app_bias_RP0I);
	mg2->Add(gr_counts_app_bias_RP0E);
	mg2->Add(gr_counts_app_bias_RP1I);
	mg2->Add(gr_counts_app_bias_RP1E);
	mg2->Add(gr_counts_app_bias_RP3E);
	mg2->Draw("APL");	
	mg2->GetXaxis()->SetTitle("Applied Bias [V]");
	mg2->GetYaxis()->SetTitle("Normalized Counts");	
	mg2->GetXaxis()->SetLimits(-10,150);
	mg2->SetMinimum(-0.1);
	mg2->SetMaximum(1.1);
	auto legend2 = new TLegend(0.52,0.1,0.9,0.3);
	legend2->SetHeader("Timepix","C"); // option "C" allows to center the header
	TString name_str = "RP0 Internal (I02-W0108 FITPix-0415) Max=";
	name_str += (Int_t)Max_Average_Counts_RP0I;
	legend2->AddEntry(gr_counts_app_bias_RP0I,name_str.Data(),"PLE");
	name_str = "RP0 External (G02-W0108 FITPix-0384) Max=";
	name_str += (Int_t)Max_Average_Counts_RP0E;
	legend2->AddEntry(gr_counts_app_bias_RP0E,name_str.Data(),"PLE");
	name_str = "RP1 Internal (F04-W0108 FITPix-0409) Max=";
	name_str += (Int_t)Max_Average_Counts_RP1I;
	legend2->AddEntry(gr_counts_app_bias_RP1I,name_str.Data(),"PLE");
	name_str = "RP1 External (H7-W00018 Katherin) Max=";
	name_str += (Int_t)Max_Average_Counts_RP1E;
	legend2->AddEntry(gr_counts_app_bias_RP1E,name_str.Data(),"PLE");
	name_str = "RP3 External (C08-W0255 FITPix-0399) Max=";
	name_str += (Int_t)Max_Average_Counts_RP3E;
	legend2->AddEntry(gr_counts_app_bias_RP3E,name_str.Data(),"PLE");
	legend2->Draw();
	
	TCanvas* c3 = new TCanvas("c3","Average Counts vs Verification Bias",1500,1000);
	c3->cd();
	gPad->SetGrid();
	TMultiGraph* mg3 = new TMultiGraph();
	mg3->Add(gr_counts_ver_bias_RP0I);
	mg3->Add(gr_counts_ver_bias_RP0E);
	mg3->Add(gr_counts_ver_bias_RP1I);
	mg3->Add(gr_counts_ver_bias_RP1E);
	mg3->Add(gr_counts_ver_bias_RP3E);
	mg3->Draw("APL");	
	mg3->GetXaxis()->SetTitle("Verification Bias [V]");
	mg3->GetYaxis()->SetTitle("Normalized Counts");	
	mg3->GetXaxis()->SetLimits(-10,150);
	mg3->SetMinimum(-0.1);
	mg3->SetMaximum(1.1);
	auto legend3 = new TLegend(0.52,0.1,0.9,0.3);
	legend3->SetHeader("Timepix","C"); // option "C" allows to center the header
	name_str = "RP0 Internal (I02-W0108 FITPix-0415) Max=";
	name_str += (Int_t)Max_Average_Counts_RP0I;
	legend3->AddEntry(gr_counts_ver_bias_RP0I,name_str.Data(),"PLE");
	name_str = "RP0 External (G02-W0108 FITPix-0384) Max=";
	name_str += (Int_t)Max_Average_Counts_RP0E;
	legend3->AddEntry(gr_counts_ver_bias_RP0E,name_str.Data(),"PLE");
	name_str = "RP1 Internal (F04-W0108 FITPix-0409) Max=";
	name_str += (Int_t)Max_Average_Counts_RP1I;
	legend3->AddEntry(gr_counts_ver_bias_RP1I,name_str.Data(),"PLE");
	name_str = "RP1 External (H7-W00018 Katherin) Max=";
	name_str += (Int_t)Max_Average_Counts_RP1E;
	legend3->AddEntry(gr_counts_ver_bias_RP1E,name_str.Data(),"PLE");
	name_str = "RP3 External (C08-W0255 FITPix-0399) Max=";
	name_str += (Int_t)Max_Average_Counts_RP3E;
	legend3->AddEntry(gr_counts_ver_bias_RP3E,name_str.Data(),"PLE");
	legend3->Draw();
	
    return 0;
}
