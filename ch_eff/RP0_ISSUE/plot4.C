int plot4()
{
    const Int_t N_PIXELS        = 256;
    const Double_t PIXEL_SIZE   = 0.055;

    TString input_file_rp0 = "output_function_7_RP0_move.root";
    TString input_file_rp1 = "output_function_7_RP1_move.root";

    cout<<"--> Input file with rp0: "<<input_file_rp0<<endl;
    TFile* _file_rp0    = TFile::Open(input_file_rp0.Data());

    cout<<"--> Input file with rp1: "<<input_file_rp1<<endl;
    TFile* _file_rp1    = TFile::Open(input_file_rp1.Data());

    TH2D* h1            = (TH2D*)_file_rp0->Get("hh_2");
    h1->SetName("h1");

    TH2D* h2            = (TH2D*)_file_rp1->Get("hh_2");
    h2->SetName("h2");

    h1->SetTitle("RP0 Timepix Hor.Axis vs Time (normalized)");
    h2->SetTitle("RP1 Timepix Hor.Axis vs Time (normalized)");

    h1->Rebin2D(2);
    h2->Rebin2D(2);

    cout<<"--> h1 Binsize XY RP1 = "<<h1->GetXaxis()->GetBinWidth(1)<<" [s] & "<<h1->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;
    cout<<"--> h2 Binsize XY RP1 = "<<h2->GetXaxis()->GetBinWidth(1)<<" [s] & "<<h2->GetYaxis()->GetBinWidth(1)<<" [mm]"<<endl;

    Double_t integral, integralErr;

    for(Int_t i = 1; i <= h1->GetNbinsX(); i++)
    {
        integral = h1->IntegralAndError(i,i,1,h1->GetNbinsY(),integralErr);
        for(Int_t j = 1; j <= h1->GetNbinsY(); j++)
        {
            if(integral > 0)
            {
                h1->SetBinError(i,j,TMath::Sqrt( TMath::Power(h1->GetBinError(i,j)/integral,2) + TMath::Power(h1->GetBinContent(i,j)*integralErr/(integral*integral),2) ));
                h1->SetBinContent(i,j,h1->GetBinContent(i,j)/integral);
            }
        }
    }

    for(Int_t i = 1; i <= h2->GetNbinsX(); i++)
    {
        integral = h2->IntegralAndError(i,i,1,h2->GetNbinsY(),integralErr);
        for(Int_t j = 1; j <= h2->GetNbinsY(); j++)
        {
            if(integral > 0)
            {
                h2->SetBinError(i,j,TMath::Sqrt( TMath::Power(h2->GetBinError(i,j)/integral,2) + TMath::Power(h2->GetBinContent(i,j)*integralErr/(integral*integral),2) ));
                h2->SetBinContent(i,j,h2->GetBinContent(i,j)/integral);
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

    return 0;
}
