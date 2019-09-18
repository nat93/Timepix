int plot2()
{
    const Int_t N_PIXELS        = 256;
    const Double_t PIXEL_SIZE   = 0.055;
    const Int_t nFiles          = 4;

    TString input_file_rp0[] = {"output_function_7_RP0_00.0mm.root", "output_function_7_RP0_62.0mm.root", "output_function_7_RP0_64.7mm.root", "output_function_7_RP0_68.5mm.root"};
    TString input_file_rp1[] = {"output_function_7_RP1_00.0mm.root", "output_function_7_RP1_62.0mm.root", "output_function_7_RP1_64.7mm.root", "output_function_7_RP1_68.5mm.root"};
    Double_t position_rp0 [] = {0.00, 62.00, 64.70, 68.50};
    Double_t position_rp1 [] = {27.3, 27.3, 28.6, 33.8};

    TFile* _file_rp0[nFiles];
    TFile* _file_rp1[nFiles];

    TH2D* hh2_rp0[nFiles];
    TH2D* hh2_rp1[nFiles];

    TCanvas* canva1[nFiles];
    TPad *pad0[nFiles];
    TPad *pad1[nFiles];

    TString hName;

    for(Int_t i = 0; i < nFiles; i++)
    {
        // RP0
        cout<<"--> Input file with rp0: "<<input_file_rp0[i]<<endl;
        _file_rp0[i]    = TFile::Open(input_file_rp0[i].Data());
        hh2_rp0[i]       = (TH2D*)_file_rp0[i]->Get("hh_2");

        hName = "hh2_rp0["; hName += i; hName += "]";
        hh2_rp0[i]->SetName(hName.Data());
        hh2_rp0[i]->SetTitle("RP0 Timepix Hor.Axis vs Time");

        hh2_rp0[i]->Rebin2D(2);
        cout<<"--> Binsize XY RP0 = "<<hh2_rp0[i]->GetXaxis()->GetBinWidth(1)<<" [s] & "<<hh2_rp0[i]->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;

        // RP1
        cout<<"--> Input file with rp1: "<<input_file_rp1[i]<<endl;
        _file_rp1[i]    = TFile::Open(input_file_rp1[i].Data());
        hh2_rp1[i]       = (TH2D*)_file_rp1[i]->Get("hh_2");

        hName = "hh2_rp1["; hName += i; hName += "]";
        hh2_rp1[i]->SetName(hName.Data());
        hh2_rp1[i]->SetTitle("RP1 Timepix Hor.Axis vs Time");

        hh2_rp1[i]->Rebin2D(2);
        cout<<"--> Binsize XY RP1 = "<<hh2_rp1[i]->GetXaxis()->GetBinWidth(1)<<" [s] & "<<hh2_rp1[i]->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;


        gStyle->SetOptStat(0);

        hName = "canva1["; hName += i; hName += "]";
        canva1[i] = new TCanvas(hName.Data(),hName.Data(),1200,1000);
        canva1[i]->Divide(1,2);

        canva1[i]->cd(1);
        hName = "pad0["; hName += i; hName += "]";
        pad0[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad0[i]->SetRightMargin(0.2);
        pad0[i]->SetBottomMargin(0.2);
        pad0[i]->Draw();
        pad0[i]->cd();
        gPad->SetGrid();
        hh2_rp0[i]->Draw("colz");

        canva1[i]->cd(2);
        hName = "pad1["; hName += i; hName += "]";
        pad1[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad1[i]->SetRightMargin(0.2);
        pad1[i]->SetBottomMargin(0.2);
        pad1[i]->Draw();
        pad1[i]->cd();
        gPad->SetGrid();
        hh2_rp1[i]->Draw("colz");
    }

    TF1* fit1_rp0[nFiles];
    TF1* fit1_rp1[nFiles];

    Double_t xMin_rp0[] = {0.0,4.0,5.0,8.0};
    Double_t xMax_rp0[] = {0.0,5.0,8.0,12.0};
    Double_t xMin_rp1[] = {1.0,1.0,1.0,7.0};
    Double_t xMax_rp1[] = {2.5,2.5,3.5,9.0};

    TH1D* hMean_rp0[nFiles];
    TH1D* hMean_rp1[nFiles];

    TH1D* hStd_rp0[nFiles];
    TH1D* hStd_rp1[nFiles];

    TH1D* hInt_rp0[nFiles];
    TH1D* hInt_rp1[nFiles];

    TCanvas* canva2[nFiles];
    TPad *pad00[nFiles];
    TPad *pad11[nFiles];

    TCanvas* canva3[nFiles];
    TPad *pad000[nFiles];
    TPad *pad111[nFiles];

    TCanvas* canva4[nFiles];
    TPad *pad0000[nFiles];
    TPad *pad1111[nFiles];

    TCanvas* canva5[nFiles];
    TPad *pad00000[nFiles];
    TPad *pad11111[nFiles];

    Double_t integral, integralErr;

    TGraphErrors * gr_rp0[nFiles];
    TGraphErrors * gr_rp1[nFiles];

    for(Int_t i = 0; i < nFiles; i++)
    {
        hName = "fit1_rp0["; hName += i; hName += "]";
        fit1_rp0[i] = new TF1(hName.Data(),"gaus",xMin_rp0[i],xMax_rp0[i]);

        hName = "fit1_rp1["; hName += i; hName += "]";
        fit1_rp1[i] = new TF1(hName.Data(),"gaus",xMin_rp1[i],xMax_rp1[i]);

        hName = "hMean_rp0["; hName += i; hName += "]";
        hMean_rp0[i] = hh2_rp0[i]->ProjectionX(hName.Data());
        hMean_rp0[i]->SetTitle("RP0 Timepix Mean Peak Position vs Time");

        hName = "hMean_rp1["; hName += i; hName += "]";
        hMean_rp1[i] = hh2_rp1[i]->ProjectionX(hName.Data());
        hMean_rp1[i]->SetTitle("RP1 Timepix Mean Peak Position vs Time");

        hName = "hStd_rp0["; hName += i; hName += "]";
        hStd_rp0[i] = hh2_rp0[i]->ProjectionX(hName.Data());
        hStd_rp0[i]->SetTitle("RP0 Timepix Std Peak Position vs Time");

        hName = "hStd_rp1["; hName += i; hName += "]";
        hStd_rp1[i] = hh2_rp1[i]->ProjectionX(hName.Data());
        hStd_rp1[i]->SetTitle("RP1 Timepix Std Peak Position vs Time");

        hName = "hInt_rp0["; hName += i; hName += "]";
        hInt_rp0[i] = hh2_rp0[i]->ProjectionX(hName.Data());
        hInt_rp0[i]->SetTitle("RP0 Timepix Int Peak Position vs Time");

        hName = "hInt_rp1["; hName += i; hName += "]";
        hInt_rp1[i] = hh2_rp1[i]->ProjectionX(hName.Data());
        hInt_rp1[i]->SetTitle("RP1 Timepix Int Peak Position vs Time");

        gr_rp0[i] = new TGraphErrors();
        hName = "gr_rp0["; hName += i; hName += "]";
        gr_rp0[i]->SetName(hName.Data());

        gr_rp1[i] = new TGraphErrors();
        hName = "gr_rp1["; hName += i; hName += "]";
        gr_rp1[i]->SetName(hName.Data());

        for(Int_t bini = 1; bini <= hh2_rp0[i]->GetNbinsX(); bini++)
        {
            TH1D* hTemp = hh2_rp0[i]->ProjectionY("hTemp_rp0",bini,bini);
            hTemp->Fit(fit1_rp0[i],"R0Q");

            hMean_rp0[i]->SetBinContent(bini,fit1_rp0[i]->GetParameter(1));
            hMean_rp0[i]->SetBinError(bini,fit1_rp0[i]->GetParError(1));

            hStd_rp0[i]->SetBinContent(bini,fit1_rp0[i]->GetParameter(2));
            hStd_rp0[i]->SetBinError(bini,fit1_rp0[i]->GetParError(2));

            integral = hTemp->IntegralAndError(1,hTemp->GetNbinsX(),integralErr);
            hInt_rp0[i]->SetBinContent(bini,integral);
            hInt_rp0[i]->SetBinError(bini,integralErr);

            gr_rp0[i]->SetPoint(gr_rp0[i]->GetN(),integral,fit1_rp0[i]->GetParameter(2));
            gr_rp0[i]->SetPointError(gr_rp0[i]->GetN()-1,integralErr,fit1_rp0[i]->GetParError(2));
        }

        for(Int_t bini = 1; bini <= hh2_rp1[i]->GetNbinsX(); bini++)
        {
            TH1D* hTemp = hh2_rp1[i]->ProjectionY("hTemp_rp1",bini,bini);
            hTemp->Fit(fit1_rp1[i],"R0Q");

            hMean_rp1[i]->SetBinContent(bini,fit1_rp1[i]->GetParameter(1));
            hMean_rp1[i]->SetBinError(bini,fit1_rp1[i]->GetParError(1));

            hStd_rp1[i]->SetBinContent(bini,fit1_rp1[i]->GetParameter(2));
            hStd_rp1[i]->SetBinError(bini,fit1_rp1[i]->GetParError(2));

            integral = hTemp->IntegralAndError(1,hTemp->GetNbinsX(),integralErr);
            hInt_rp1[i]->SetBinContent(bini,integral);
            hInt_rp1[i]->SetBinError(bini,integralErr);

            gr_rp1[i]->SetPoint(gr_rp1[i]->GetN(),integral,fit1_rp1[i]->GetParameter(2));
            gr_rp1[i]->SetPointError(gr_rp1[i]->GetN()-1,integralErr,fit1_rp1[i]->GetParError(2));
        }

        gStyle->SetOptStat(0);

        hName = "canva2["; hName += i; hName += "]";
        canva2[i] = new TCanvas(hName.Data(),hName.Data(),1200,1000);
        canva2[i]->Divide(1,2);

        canva2[i]->cd(1);
        hName = "pad00["; hName += i; hName += "]";
        pad00[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad00[i]->Draw();
        pad00[i]->cd();
        gPad->SetGrid();
        hMean_rp0[i]->SetMinimum(0.0);
        hMean_rp0[i]->SetMaximum(14.0);
        hMean_rp0[i]->Draw();

        canva2[i]->cd(2);
        hName = "pad11["; hName += i; hName += "]";
        pad11[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad11[i]->Draw();
        pad11[i]->cd();
        gPad->SetGrid();
        hMean_rp1[i]->SetMinimum(0.0);
        hMean_rp1[i]->SetMaximum(14.0);
        hMean_rp1[i]->Draw();

        hName = "canva3["; hName += i; hName += "]";
        canva3[i] = new TCanvas(hName.Data(),hName.Data(),1200,1000);
        canva3[i]->Divide(1,2);

        canva3[i]->cd(1);
        hName = "pad000["; hName += i; hName += "]";
        pad000[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad000[i]->Draw();
        pad000[i]->cd();
        gPad->SetGrid();
        hStd_rp0[i]->SetMinimum(0.0);
        hStd_rp0[i]->SetMaximum(2.0);
        hStd_rp0[i]->Draw();

        canva3[i]->cd(2);
        hName = "pad111["; hName += i; hName += "]";
        pad111[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad111[i]->Draw();
        pad111[i]->cd();
        gPad->SetGrid();
        hStd_rp1[i]->SetMinimum(0.0);
        hStd_rp1[i]->SetMaximum(2.0);
        hStd_rp1[i]->Draw();

        hName = "canva4["; hName += i; hName += "]";
        canva4[i] = new TCanvas(hName.Data(),hName.Data(),1200,1000);
        canva4[i]->Divide(1,2);

        canva4[i]->cd(1);
        hName = "pad0000["; hName += i; hName += "]";
        pad0000[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad0000[i]->Draw();
        pad0000[i]->cd();
        gPad->SetGrid();
        hInt_rp0[i]->SetMinimum(0.0);
        hInt_rp0[i]->SetMaximum(1.0e6);
        hInt_rp0[i]->Draw();

        canva4[i]->cd(2);
        hName = "pad1111["; hName += i; hName += "]";
        pad1111[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad1111[i]->Draw();
        pad1111[i]->cd();
        gPad->SetGrid();
        hInt_rp1[i]->SetMinimum(0.0);
        hInt_rp1[i]->SetMaximum(50.0e6);
        hInt_rp1[i]->Draw();

        hName = "canva5["; hName += i; hName += "]";
        canva5[i] = new TCanvas(hName.Data(),hName.Data(),1200,1000);
        canva5[i]->Divide(1,2);

        canva5[i]->cd(1);
        hName = "pad0000["; hName += i; hName += "]";
        pad00000[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad00000[i]->Draw();
        pad00000[i]->cd();
        gPad->SetGrid();
        gr_rp0[i]->Draw("AP");

        canva5[i]->cd(2);
        hName = "pad11111["; hName += i; hName += "]";
        pad11111[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad11111[i]->Draw();
        pad11111[i]->cd();
        gPad->SetGrid();
        gr_rp1[i]->Draw("AP");
    }

    return 0;
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

    void calculate_angular_error(Double_t p /*MeV/c*/)
    {
            Double_t z          = 1;            // charge of particle
            Double_t s          = 0.03;         // [cm] thickness of the material
            Double_t beta_c     = 1;            // beta*c velocity of the particle
//            Double_t s_Xl       = 0.026;      // s/Xl Timepix;
            Double_t s_Xl       = 4.0/17.60;   // s/Xl Stainless Steel;
            Double_t Theta_0    = 13.6*z*TMath::Sqrt(s_Xl)*(1.0 + 0.038*TMath::Log(s_Xl))/(p*beta_c);
            printf("--> Theta = %10.10f [urad]\n",Theta_0*1e6);
    }
