int plot()
{
    const Int_t N_PIXELS        = 256;
    const Double_t PIXEL_SIZE   = 0.055;
    const Int_t nFiles          = 4;

    TString input_file_rp0[] = {"output_function_1_RP0_00.0mm.root", "output_function_1_RP0_62.0mm.root", "output_function_1_RP0_64.7mm.root", "output_function_1_RP0_68.5mm.root"};
    TString input_file_rp1[] = {"output_function_1_RP1_00.0mm.root", "output_function_1_RP1_62.0mm.root", "output_function_1_RP1_64.7mm.root", "output_function_1_RP1_68.5mm.root"};
    Double_t position_rp0 [] = {0.00, 62.00, 64.72, 68.50};
    Double_t position_rp1 [] = {27.3, 27.3, 28.6, 33.8};

    TFile* _file_rp0[nFiles];
    TFile* _file_rp1[nFiles];

    TH2D* h1_rp0[nFiles];
    TH2D* h1_rp1[nFiles];

    TH2D* hh1_rp0[nFiles];
    TH2D* hh1_rp1[nFiles];

    TH1D* projhh1_rp0[nFiles];
    TH1D* projhh1_rp1[nFiles];

    TH1D* projhhh1_rp0[nFiles];
    TH1D* projhhh1_rp1[nFiles];

    TCanvas* canva1[nFiles];
    TPad *pad0[nFiles];
    TPad *pad1[nFiles];

    TString hName;
    Double_t integral, integralErr;
    Int_t locmax, locmay, locmaz;

    for(Int_t i = 0; i < nFiles; i++)
    {
        // RP0
        cout<<"--> Input file with rp0: "<<input_file_rp0[i]<<endl;
        _file_rp0[i]    = TFile::Open(input_file_rp0[i].Data());
        h1_rp0[i]       = (TH2D*)_file_rp0[i]->Get("h_1");

        hName = "h1_rp0["; hName += i; hName += "]";
        h1_rp0[i]->SetName(hName.Data());

        hName = "hh1_rp0["; hName += i; hName += "]";
        hh1_rp0[i] = new TH2D(hName.Data(),"RP0 Timepix YX [mm]",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);

        // RP1
        cout<<"--> Input file with rp1: "<<input_file_rp1[i]<<endl;
        _file_rp1[i]    = TFile::Open(input_file_rp1[i].Data());
        h1_rp1[i]       = (TH2D*)_file_rp1[i]->Get("h_1");

        hName = "h1_rp1["; hName += i; hName += "]";
        h1_rp1[i]->SetName(hName.Data());

        hName = "hh1_rp1["; hName += i; hName += "]";
        hh1_rp1[i] = new TH2D(hName.Data(),"RP1 Timepix YX [mm]",N_PIXELS,0,N_PIXELS*PIXEL_SIZE,N_PIXELS,0,N_PIXELS*PIXEL_SIZE);

        // Loop
        for(Int_t binx = 1; binx <= N_PIXELS; binx++)
        {
            for(Int_t biny = 1; biny <= N_PIXELS; biny++)
            {
                hh1_rp0[i]->SetBinContent(N_PIXELS-biny+1,binx,h1_rp0[i]->GetBinContent(binx,biny));
                hh1_rp0[i]->SetBinError(N_PIXELS-biny+1,binx,h1_rp0[i]->GetBinError(binx,biny));

                hh1_rp1[i]->SetBinContent(N_PIXELS-biny+1,binx,h1_rp1[i]->GetBinContent(binx,biny));
                hh1_rp1[i]->SetBinError(N_PIXELS-biny+1,binx,h1_rp1[i]->GetBinError(binx,biny));
            }
        }

        hh1_rp0[i]->Rebin2D(2);
        cout<<"--> Binsize XY RP0 = "<<hh1_rp0[i]->GetXaxis()->GetBinWidth(1)<<" [mm] & "<<hh1_rp0[i]->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;

        hh1_rp1[i]->Rebin2D(2);
        cout<<"--> Binsize XY RP1 = "<<hh1_rp1[i]->GetXaxis()->GetBinWidth(1)<<" [mm] & "<<hh1_rp1[i]->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;

//        hh1_rp0[i]->GetMaximumBin(locmax,locmay,locmaz);
//        scale2Dhisto(hh1_rp0[i], hh1_rp0[i]->GetBinContent(locmax,locmay),hh1_rp0[i]->GetBinError(locmax,locmay));

//        hh1_rp1[i]->GetMaximumBin(locmax,locmay,locmaz);
//        scale2Dhisto(hh1_rp1[i], hh1_rp1[i]->GetBinContent(locmax,locmay),hh1_rp1[i]->GetBinError(locmax,locmay));

        integral = hh1_rp0[i]->IntegralAndError(1,hh1_rp0[i]->GetNbinsX(),1,hh1_rp0[i]->GetNbinsY(),integralErr);
        scale2Dhisto(hh1_rp0[i], integral,integralErr);

        integral = hh1_rp1[i]->IntegralAndError(1,hh1_rp1[i]->GetNbinsX(),1,hh1_rp1[i]->GetNbinsY(),integralErr);
        scale2Dhisto(hh1_rp1[i], integral,integralErr);

        gStyle->SetOptStat(0);

        hName = "canva["; hName += i; hName += "]";
        canva1[i] = new TCanvas(hName.Data(),hName.Data(),1200,600);
        canva1[i]->Divide(2,1);

        canva1[i]->cd(1);
        hName = "pad0["; hName += i; hName += "]";
        pad0[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad0[i]->SetRightMargin(0.2);
        pad0[i]->SetBottomMargin(0.2);
        pad0[i]->Draw();
        pad0[i]->cd();
        gPad->SetGrid();
        hh1_rp0[i]->Draw("colz");

        canva1[i]->cd(2);
        hName = "pad1["; hName += i; hName += "]";
        pad1[i] = new TPad(hName.Data(),hName.Data(),0,0,1,1);
        pad1[i]->SetRightMargin(0.2);
        pad1[i]->SetBottomMargin(0.2);
        pad1[i]->Draw();
        pad1[i]->cd();
        gPad->SetGrid();
        hh1_rp1[i]->Draw("colz");

        hName = "projhh1_rp0["; hName += i; hName += "]";
        projhh1_rp0[i] = hh1_rp0[i]->ProjectionX(hName.Data());

        hName = "projhhh1_rp0["; hName += i; hName += "]";
        projhhh1_rp0[i] = new TH1D(hName.Data(),hName.Data(),projhh1_rp0[i]->GetNbinsX()*2,-N_PIXELS*PIXEL_SIZE,N_PIXELS*PIXEL_SIZE);

        for(Int_t j = 1; j <= projhh1_rp0[i]->GetNbinsX(); j++)
        {
            Double_t x = projhh1_rp0[i]->GetBinCenter(j);
            Double_t y = projhh1_rp0[i]->GetBinContent(j);
            Double_t ey = projhh1_rp0[i]->GetBinError(j);

            x = x - (position_rp0[i] - position_rp0[1]);

            projhhh1_rp0[i]->SetBinContent(projhhh1_rp0[i]->GetXaxis()->FindBin(x),y);
            projhhh1_rp0[i]->SetBinError(projhhh1_rp0[i]->GetXaxis()->FindBin(x),ey);
        }

        hName = "projhh1_rp1["; hName += i; hName += "]";
        projhh1_rp1[i] = hh1_rp1[i]->ProjectionX(hName.Data());

        hName = "projhhh1_rp1["; hName += i; hName += "]";
        projhhh1_rp1[i] = new TH1D(hName.Data(),hName.Data(),projhh1_rp1[i]->GetNbinsX()*2,-N_PIXELS*PIXEL_SIZE,N_PIXELS*PIXEL_SIZE);

        for(Int_t j = 1; j <= projhh1_rp1[i]->GetNbinsX(); j++)
        {
            Double_t x = projhh1_rp1[i]->GetBinCenter(j);
            Double_t y = projhh1_rp1[i]->GetBinContent(j);
            Double_t ey = projhh1_rp1[i]->GetBinError(j);

            x = x - (position_rp1[i] - position_rp1[0]);

            projhhh1_rp1[i]->SetBinContent(projhhh1_rp1[i]->GetXaxis()->FindBin(x),y);
            projhhh1_rp1[i]->SetBinError(projhhh1_rp1[i]->GetXaxis()->FindBin(x),ey);
        }
    }

    TCanvas* canva2 = new TCanvas("canva2","canva2",1200,600);
    canva2->cd();
    TLegend* legend2 = new TLegend(0.5,0.5,0.8,0.8);

    for(Int_t i = 0; i < nFiles; i++)
    {
        projhh1_rp0[i]->SetMaximum(0.1);
        projhh1_rp0[i]->SetMinimum(0.0);
        projhh1_rp0[i]->SetLineColor(1+i);
        projhh1_rp0[i]->SetLineWidth(2);

        projhh1_rp0[i]->Draw("HIST SAME");

        hName.Form("RP0 is at %.2f mm",position_rp0[i]);
        legend2->AddEntry(projhh1_rp0[i],hName.Data(),"l");
    }
    legend2->Draw();

    TCanvas* canva3 = new TCanvas("canva3","canva3",1200,600);
    canva3->cd();
    TLegend* legend3 = new TLegend(0.5,0.5,0.8,0.8);
    TF1* fit_rp1[nFiles];
    Double_t xMin[] = {1.0,1.0,1.0,7.0};
    Double_t xMax[] = {2.5,2.5,3.5,9.0};

    TF1* rp0Kick = new TF1("rp0Kick","pol1",0,100);
    rp0Kick->SetParameters(2.65441e-05,2.33628e-05);

    Double_t alignmentRelPos    = 1.79604;
    Double_t alignmentAbsPos    = 36.8;
    Double_t timepixOffset      = 2.97;
    Double_t meanKickcr2[nFiles], meanKickcr2err[nFiles], sigmaKickcr2[nFiles], sigmaKickcr2err[nFiles], std3[nFiles], std3err[nFiles];

    for(Int_t i = 0; i < nFiles; i++)
    {
        projhh1_rp1[i]->SetMaximum(0.1);
        projhh1_rp1[i]->SetMinimum(0.0);
        projhh1_rp1[i]->SetLineColor(1+i);
        projhh1_rp1[i]->SetLineWidth(2);

        projhh1_rp1[i]->Draw("HIST SAME");
        hName = "fit_rp1["; hName += i; hName += "]";
        fit_rp1[i] = new TF1(hName.Data(),"gaus",xMin[i],xMax[i]);
        projhh1_rp1[i]->Fit(fit_rp1[i],"R+");
        fit_rp1[i]->SetLineStyle(9);
        fit_rp1[i]->SetLineWidth(3);
        fit_rp1[i]->SetLineColor(kMagenta);
        fit_rp1[i]->Draw("SAME L");

        meanKickcr2[i]         = (-1)*( rp0Kick->Eval( (-1)*(alignmentRelPos + timepixOffset + (alignmentAbsPos-position_rp1[i]+fit_rp1[i]->GetParameter(1))) )*1e6 );
        meanKickcr2err[i]      = TMath::Abs((-1)*( rp0Kick->Eval( (-1)*(alignmentRelPos + timepixOffset + (alignmentAbsPos-position_rp1[i]+fit_rp1[i]->GetParameter(1)+fit_rp1[i]->GetParError(1))) )*1e6 ) - meanKickcr2[i]);
        sigmaKickcr2[i]        = TMath::Abs((-1)*( rp0Kick->Eval( (-1)*(alignmentRelPos + timepixOffset + (alignmentAbsPos-position_rp1[i]+fit_rp1[i]->GetParameter(1)+fit_rp1[i]->GetParameter(2))) )*1e6 ) - meanKickcr2[i]);
        sigmaKickcr2err[i]     = TMath::Sqrt( TMath::Power(1.89575e-05*TMath::Sqrt(fit_rp1[i]->GetParError(1)*fit_rp1[i]->GetParError(1) + fit_rp1[i]->GetParError(2)*fit_rp1[i]->GetParError(2))*1e6,2) + TMath::Power(meanKickcr2err[i]*meanKickcr2err[i],2));

        cout<<"--> mean = "<<meanKickcr2[i]<<" +/- "<<meanKickcr2err[i]<<endl;
        cout<<"--> sigma = "<<sigmaKickcr2[i]<<" +/- "<<sigmaKickcr2err[i]<<endl;

        hName.Form("RP0 is at %.2f mm",position_rp0[i]);
        legend3->AddEntry(projhh1_rp1[i],hName.Data(),"l");

        std3[i]     = TMath::Sqrt(sigmaKickcr2[i]*sigmaKickcr2[i] - sigmaKickcr2[0]*sigmaKickcr2[0]);
        std3err[i]  = TMath::Sqrt(  TMath::Power(sigmaKickcr2err[i]*sigmaKickcr2[i]/TMath::Sqrt(sigmaKickcr2[i]*sigmaKickcr2[i] - sigmaKickcr2[0]*sigmaKickcr2[0]),2) +
                                    TMath::Power(sigmaKickcr2err[0]*sigmaKickcr2[0]/TMath::Sqrt(sigmaKickcr2[i]*sigmaKickcr2[i] - sigmaKickcr2[0]*sigmaKickcr2[0]),2));

        cout<<"--> Std MCS RP0: "<<std3[i]<<" +/- "<<std3err[i]<<" urad"<<endl;
    }
    legend3->Draw();

    Double_t z = 1;             // charge of particle
    Double_t p = 270.0e3;	// [MeV] particle momentum
    Double_t beta_c = 1; 	// beta*c velocity of the particle

    const Int_t nnn = 50000;
    Double_t sXL[nnn], theta0[nnn];
    for(Int_t i = 0; i < nnn; i++)
    {
            sXL[i] 		= 0.000001 + i*0.000005;
            theta0[i] 	= (1e6)*13.6*z*TMath::Sqrt(sXL[i])*(1.0 + 0.038*TMath::Log(sXL[i]))/(p*beta_c);
    }

    TGraph* grRadLen = new TGraph(nnn,sXL,theta0);

    Double_t x1,y1, x2,y2, sXLmean[nFiles];
    for(Int_t j = 0; j < nFiles; j++)
    {
        for(Int_t i = 1; i < grRadLen->GetN(); i++)
        {
                grRadLen->GetPoint(i,x2,y2);
                grRadLen->GetPoint(i-1,x1,y1);

                if(std3[j] < y2 && std3[j] >= y1)
                {
                        sXLmean[j] = (x2 + x1)/2.0;
                }
        }
    }


    TCanvas* c_7 = new TCanvas("c_7","c_7",1000,1000);
    c_7->cd();
    TPad *pad_7 = new TPad("pad_7","pad_7",0,0,1,1);
    pad_7->SetLeftMargin(0.15);
    pad_7->Draw();
    pad_7->cd();
    gStyle->SetOptStat(0);
    gPad->SetGrid();

    grRadLen->SetTitle("#theta_{MCS} = #frac{13.6 MeV}{#betacp}z#sqrt{s/X_{0}}[1+0.038ln(s/X_{0})]");
    grRadLen->GetYaxis()->SetTitle("Scattering (#theta_{MCS}) [#murad]");
    grRadLen->GetXaxis()->SetTitle("Thickness (s) [X_{0}^{-1}]");
    grRadLen->GetXaxis()->CenterTitle();
    grRadLen->GetYaxis()->CenterTitle();
    grRadLen->GetYaxis()->SetTitleOffset(1.8);
    grRadLen->GetXaxis()->SetTitleOffset(1);
    grRadLen->Draw("APC");

    TLine* line1[nFiles];
    TLine* line2[nFiles];

    cout<<endl;
    for(Int_t i = 0; i < nFiles; i++)
    {
        line1[i] = new TLine(0,std3[i],sXLmean[i],std3[i]);
        line2[i] = new TLine(sXLmean[i],0,sXLmean[i],std3[i]);

        line1[i]->SetLineWidth(2);
        line2[i]->SetLineWidth(2);

        line1[i]->SetLineColor(1+i);
        line2[i]->SetLineColor(1+i);

        line1[i]->Draw();
        line2[i]->Draw();

        cout<<"--> s/X_{0} = "<<sXLmean[i]<<endl;
    }
    cout<<endl;

    Double_t X0_W 	= 3.504; // Radiation length [mm]
    Double_t X0_SS      = 17.60; // Radiation length [mm]

    for(Int_t i = 0; i < nFiles; i++)
    {
        cout<<"--> s_W = "<<sXLmean[i]*X0_W<<endl;
        cout<<"--> s_SS = "<<sXLmean[i]*X0_SS<<endl;
    }
    cout<<endl;

    calculate_angular_error(270e3);

    TCanvas* canva22 = new TCanvas("canva22","canva22",1200,600);
    canva22->cd();
    TLegend* legend22 = new TLegend(0.5,0.5,0.8,0.8);

    for(Int_t i = 0; i < nFiles; i++)
    {
        projhhh1_rp0[i]->SetMaximum(0.1);
        projhhh1_rp0[i]->SetMinimum(0.0);
        projhhh1_rp0[i]->SetLineColor(1+i);
        projhhh1_rp0[i]->SetLineWidth(2);

        projhhh1_rp0[i]->Draw("HIST SAME");

        hName.Form("RP0 is at %.2f mm",position_rp0[i]);
        legend22->AddEntry(projhhh1_rp0[i],hName.Data(),"l");
    }
    legend22->Draw();

    TCanvas* canva33 = new TCanvas("canva33","canva33",1200,600);
    canva33->cd();
    TLegend* legend33 = new TLegend(0.5,0.5,0.8,0.8);

    for(Int_t i = 0; i < nFiles; i++)
    {
        projhhh1_rp1[i]->SetMaximum(0.1);
        projhhh1_rp1[i]->SetMinimum(0.0);
        projhhh1_rp1[i]->SetLineColor(1+i);
        projhhh1_rp1[i]->SetLineWidth(2);

        projhhh1_rp1[i]->Draw("HIST SAME");

        hName.Form("RP1 is at %.2f mm",position_rp1[i]);
        legend33->AddEntry(projhhh1_rp1[i],hName.Data(),"l");
    }
    legend33->Draw();

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
