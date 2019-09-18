int plot3()
{
    const Int_t N_PIXELS        = 256;
    const Double_t PIXEL_SIZE   = 0.055;

    TString input_file_rp0 = "output_function_7_RP0_crystal2_angular_scan.root";
    TString input_file_rp1 = "output_function_7_RP1_crystal2_angular_scan.root";

    cout<<"--> Input file with rp0: "<<input_file_rp0<<endl;
    TFile* _file_rp0    = TFile::Open(input_file_rp0.Data());

    cout<<"--> Input file with rp1: "<<input_file_rp1<<endl;
    TFile* _file_rp1    = TFile::Open(input_file_rp1.Data());

    TH2D* h1            = (TH2D*)_file_rp0->Get("hh_2");
    h1->SetName("h1");

    TH2D* h2            = (TH2D*)_file_rp0->Get("hh_2");
    h2->SetName("h2");

    TH2D* hh1            = (TH2D*)_file_rp1->Get("hh_2");
    hh1->SetName("hh1");

    TH2D* hh2            = (TH2D*)_file_rp1->Get("hh_2");
    hh2->SetName("hh2");

    h1->SetTitle("RP0 Timepix Hor.Axis vs Time");
    h2->SetTitle("RP0 Timepix Hor.Axis vs Time (normalized)");

    hh1->SetTitle("RP1 Timepix Hor.Axis vs Time");
    hh2->SetTitle("RP1 Timepix Hor.Axis vs Time (normalized)");

    h1->Rebin2D(2);
    h2->Rebin2D(2);
    hh1->Rebin2D(2);
    hh2->Rebin2D(2);

    cout<<"--> h1 Binsize XY RP0 = "<<h1->GetXaxis()->GetBinWidth(1)<<" [s] & "<<h1->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;
    cout<<"--> h2 Binsize XY RP0 = "<<h2->GetXaxis()->GetBinWidth(1)<<" [s] & "<<h2->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;
    cout<<"--> hh1 Binsize XY RP1 = "<<hh1->GetXaxis()->GetBinWidth(1)<<" [s] & "<<hh1->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;
    cout<<"--> hh2 Binsize XY RP1 = "<<hh2->GetXaxis()->GetBinWidth(1)<<" [s] & "<<hh2->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;

    Double_t integral, integralErr;

    for(Int_t i = 1; i <= h1->GetNbinsX(); i++)
    {
        integral = h1->IntegralAndError(i,i,1,h1->GetNbinsY(),integralErr);
        for(Int_t j = 1; j <= h2->GetNbinsY(); j++)
        {
            if(integral > 0)
            {
                h2->SetBinError(i,j,TMath::Sqrt( TMath::Power(h2->GetBinError(i,j)/integral,2) + TMath::Power(h2->GetBinContent(i,j)*integralErr/(integral*integral),2) ));
                h2->SetBinContent(i,j,h2->GetBinContent(i,j)/integral);
            }
        }
    }

    for(Int_t i = 1; i <= hh1->GetNbinsX(); i++)
    {
        integral = hh1->IntegralAndError(i,i,1,hh1->GetNbinsY(),integralErr);
        for(Int_t j = 1; j <= hh2->GetNbinsY(); j++)
        {
            if(integral > 0)
            {
                hh2->SetBinError(i,j,TMath::Sqrt( TMath::Power(hh2->GetBinError(i,j)/integral,2) + TMath::Power(hh2->GetBinContent(i,j)*integralErr/(integral*integral),2) ));
                hh2->SetBinContent(i,j,hh2->GetBinContent(i,j)/integral);
            }
        }
    }

    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1","c1",1200,1000);
    c1->Divide(1,2);
    c1->cd(1);
    h1->Draw("colz");
    c1->cd(2);
    h2->Draw("colz");

    TCanvas* c2 = new TCanvas("c2","c2",1200,1000);
    c2->Divide(1,2);
    c2->cd(1);
    hh1->Draw("colz");
    c2->cd(2);
    hh2->Draw("colz");

    TCanvas* c3 = new TCanvas("c3","c3",1200,1000);
    c3->Divide(1,2);
    c3->cd(1);
    h2->Draw("colz");
    c3->cd(2);
    hh2->Draw("colz");

    //==================================================================================================================//
    // Kick Angle CRYSTAL2 at RP0
    //==================================================================================================================//
    //Minimizer is Linear
    //Chi2                      =  5.04791e-36
    //NDf                       =          198
    //p0                        =  2.30877e-05   +/-   1.12908e-20
    //p1                        =  1.62293e-05   +/-   1.95557e-22
    //==================================================================================================================//
    // Kick Angle CRYSTAL2 at RP1
    //==================================================================================================================//
    //Minimizer is Linear
    //Chi2                      =   5.0221e-35
    //NDf                       =           98
    //p0                        = -9.16228e-06   +/-   1.42105e-19
    //p1                        =  1.89575e-05   +/-   2.47994e-21

    TF1* cry2Kick_rp0 = new TF1("cry2Kick_rp0","pol1",0,100);
    cry2Kick_rp0->SetParameters(2.30877e-05,1.62293e-05);

    TF1* cry2Kick_rp1 = new TF1("cry2Kick_rp1","pol1",0,100);
    cry2Kick_rp1->SetParameters(-9.16228e-06,1.89575e-05);

    Double_t alignmentRelPos_rp0    = 2.47114;
    Double_t alignmentAbsPos_rp0    = 71.45;
    Double_t actualAbsPos_rp0       = 68.46;
    Double_t timepixOffset_rp0      = 3.0;

    Double_t alignmentRelPos_rp1    = 1.79604;
    Double_t alignmentAbsPos_rp1    = 36.76;
    Double_t actualAbsPos_rp1       = 33.76;
    Double_t timepixOffset_rp1      = 3.0;

    TGaxis *axisCry2_rp0 = new TGaxis(-20,0,-20,0.055*256,(-1)*( cry2Kick_rp0->Eval( (-1)*(alignmentRelPos_rp0 + timepixOffset_rp0 + (alignmentAbsPos_rp0 - actualAbsPos_rp0)) )*1e6 ),
                                                          (-1)*( cry2Kick_rp0->Eval( (-1)*(alignmentRelPos_rp0 + timepixOffset_rp0 + (alignmentAbsPos_rp0 - actualAbsPos_rp0 + 0.055*256)) )*1e6 ),50510,"-");
    axisCry2_rp0->SetName("axisCry2_rp0");
    axisCry2_rp0->SetLabelOffset(0.01);
    axisCry2_rp0->SetLineColor(kBlue);
    axisCry2_rp0->SetLabelColor(kBlue);
    axisCry2_rp0->SetLabelSize(0.03);
    axisCry2_rp0->SetLabelFont(42);
    axisCry2_rp0->SetNdivisions(0*10000 + 5*100 + 20*1);
    axisCry2_rp0->SetTitle("Crystal kick [#murad]");
    axisCry2_rp0->CenterTitle();
    axisCry2_rp0->SetTitleOffset(1.0);
    axisCry2_rp0->SetTitleSize(0.03);
    axisCry2_rp0->SetTitleFont(42);
    axisCry2_rp0->SetTitleColor(kBlue);
    c3->cd(1);
    axisCry2_rp0->Draw();

    TGaxis *axisCry2_rp1 = new TGaxis(-20,0,-20,0.055*256,(-1)*( cry2Kick_rp1->Eval( (-1)*(alignmentRelPos_rp1 + timepixOffset_rp1 + (alignmentAbsPos_rp1 - actualAbsPos_rp1)) )*1e6 ),
                                                          (-1)*( cry2Kick_rp1->Eval( (-1)*(alignmentRelPos_rp1 + timepixOffset_rp1 + (alignmentAbsPos_rp1 - actualAbsPos_rp1 + 0.055*256)) )*1e6 ),50510,"-");
    axisCry2_rp1->SetName("axisCry2_rp1");
    axisCry2_rp1->SetLabelOffset(0.01);
    axisCry2_rp1->SetLineColor(kBlue);
    axisCry2_rp1->SetLabelColor(kBlue);
    axisCry2_rp1->SetLabelSize(0.03);
    axisCry2_rp1->SetLabelFont(42);
    axisCry2_rp1->SetNdivisions(0*10000 + 5*100 + 20*1);
    axisCry2_rp1->SetTitle("Crystal kick [#murad]");
    axisCry2_rp1->CenterTitle();
    axisCry2_rp1->SetTitleOffset(1.0);
    axisCry2_rp1->SetTitleSize(0.03);
    axisCry2_rp1->SetTitleFont(42);
    axisCry2_rp1->SetTitleColor(kBlue);
    c3->cd(2);
    axisCry2_rp1->Draw();


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
