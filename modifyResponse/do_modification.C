#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

Double_t SingleGaussian(Double_t *x, Double_t *par){
  return par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));
}

double ComputeGaussianValue(double scale, double mean, double sigma, float x){
  return scale*TMath::Exp(-0.5*((x - mean)/sigma)*((x - mean)/sigma));
}

void do_modification(){

  //core + 25%
  float core_mod = 1.25;
  float tail_mod = 1.0;
  std::string plotdir = "core_plus_25percent";

  //core - 25%
  //float core_mod = 0.75;
  //float tail_mod = 1.0;
  //std::string plotdir = "core_minus_25percent";

  //tail + 25%
  //float core_mod = 1.0;
  //float tail_mod = 1.25;
  //std::string plotdir = "tail_plus_25percent";

  //tail - 25%
  //float core_mod = 1.0;
  //float tail_mod = 0.75;
  //std::string plotdir = "tail_minus_25percent";

  TFile* outfile_fits = new TFile(Form("template_fits_%s.root", plotdir.c_str()), "RECREATE");
  TFile* outfile_templates = new TFile(Form("templates_%s.root", plotdir.c_str()), "RECREATE");

  TFile* responseFile = new TFile("input_templates.root", "READ");
  TCanvas* c1 = new TCanvas();
  TIter it(responseFile->GetListOfKeys());
  TKey* key;
  std::string keep = "h_tot_JetAll_ResponsePt";
  //std::string keep = "h_tot_JetAll_ResponsePt_Pt0_Eta0";
  while ( (key = (TKey *)it()) ) {
    if (strncmp (key->GetTitle(), keep.c_str(), keep.length()) != 0) continue;
    std::string hist_name = (key->GetTitle());
    TH1F* h = (TH1F*) responseFile->Get(TString(hist_name));
    if(h->GetEntries() < 100) continue;
    h->Scale(1.0/h->Integral());

    TH1F* h_tail = (TH1F*) h->Clone();
    TH1F* h_core = (TH1F*) h->Clone();
    TH1F* h_new_core = (TH1F*) h->Clone();
    TH1F* h_new_tail = (TH1F*) h->Clone();

    TF1* fit = new TF1(Form("fit_%s", hist_name.c_str()), SingleGaussian, 0.03, 1.0, 3);
    fit->SetTitle(Form("fit_%s", hist_name.c_str()));
    float mean = h->GetBinCenter(h->GetMaximumBin());
    float max = h->GetBinContent(h->GetMaximumBin());
    
    fit->SetParameters(max, mean, 0.1);
    fit->SetParLimits(0, max - 0.0005, max + 0.0005);
    fit->SetParLimits(1, mean - 0.02, mean + 0.02);
    float side = h->GetRMS();
    if(side > 0.1) side *= 0.5;
    fit->SetRange(mean - side, mean + side);
    h->Fit(Form("fit_%s", hist_name.c_str()), "R");
    fit->SetRange(0, 3);
    fit->Print();
    fit->Update();
    c1->cd();
    c1->Print(Form("fits/fit_%s.pdf", hist_name.c_str()));

    double g_scale = fit->GetParameter(0);
    double g_mean  = fit->GetParameter(1);
    double g_sigma = fit->GetParameter(2);

    for(int i=1; i<=h->GetNbinsX(); i++){
        float gaussian_contribution = ComputeGaussianValue(g_scale, g_mean, g_sigma, h->GetBinCenter(i));
        float tail_contribution = h->GetBinContent(i) - gaussian_contribution;
        h_tail->SetBinContent(i, std::max(tail_contribution, 0)); 
        h_core->SetBinContent(i, gaussian_contribution); 
    }
    h->Draw();
    h_tail->SetLineColor(kGreen);
    h_tail->Draw("same");
    h_core->SetLineColor(kRed);
    h_core->Draw("same");
    c1->Print(Form("plots/core_mod_%s.pdf", hist_name.c_str()));
    c1->Clear();


    for(int i=1; i<=h_new_core->GetNbinsX(); i++){
      float gaussian_contribution = ComputeGaussianValue(g_scale, g_mean, core_mod*g_sigma, h_new_core->GetBinCenter(i));
      h_new_core->SetBinContent(i, gaussian_contribution); 
    }
    for(int i=1; i<=h_new_tail->GetNbinsX(); i++){
      h_new_tail->SetBinContent(i, tail_mod*h_tail->GetBinContent(i)); 
    }

    //rescale new core histogram to have same area as original core
    float core_integral_before = h_core->Integral();
    float core_integral_after = h_new_core->Integral();
    float integral_scale_factor = core_integral_before/core_integral_after;
    h_new_core->Scale(integral_scale_factor);

    outfile_fits->cd();
    fit->SetParameter(2, core_mod*g_sigma);
    fit->SetParameter(0, integral_scale_factor*g_sigma);
    fit->Write();

    delete fit;
    delete c1;
    TCanvas* c1 = new TCanvas();
    c1->SetLogy();

    h_new_core->Add(h_new_tail, 1.0);
    //h_new_core->Scale(1.0/h_new_core->Integral());
    std::cout << "h_new_core->Integral() = " << h_new_core->Integral() << std::endl;
    h_core->SetMarkerColor(kRed);
    h_core->SetMarkerStyle(20);
    h_core->SetMarkerSize(0.8);
    h_tail->SetMarkerColor(kGreen);
    h_tail->SetMarkerStyle(20);
    h_tail->SetMarkerSize(0.8);
    h_new_core->SetMarkerStyle(20);
    h_new_core->SetMarkerSize(0.8);
    h->SetFillColor(kGray+1);
    h->SetLineColor(kGray+1);
    //h->SetMaximum(1.2 * std::max(h_new_core->GetMaximum(), h->GetMaximum()));
    h->SetMaximum(2.0 * std::max(h_new_core->GetMaximum(), h->GetMaximum()));
    h->SetMinimum(0.0001);
    h->Draw("H");
    h_core->Draw("SAMEL");
    h_tail->Draw("SAMEL");
    h_new_core->Draw("SAME");
    TLegend* legend = new TLegend(0.6, 0.6, 0.85, 0.85);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->AddEntry( h, "Original Template", "f" );
    if(plotdir == "core_plus_25percent") legend->AddEntry( h_new_core, "New Template (width + 25%)", "p" );
    else if(plotdir == "core_minus_25percent") legend->AddEntry( h_new_core, "New Template (width - 25%)", "p" );
    else legend->AddEntry( h_new_core, "New Template", "p" );
    legend->AddEntry( h_core, "Original Core", "p" );
    legend->AddEntry( h_tail, "Original Tail", "p" );
    legend->Draw("SAME");
    c1->Print(Form("%s/new_core_%s.pdf", plotdir.c_str(), hist_name.c_str()));
    outfile_templates->cd();
    h_new_core->Write(hist_name.c_str());
    responseFile->cd();

  }
  responseFile->Close();
  delete responseFile;
  delete c1;

  outfile_fits->Close();
  delete outfile_fits;
  outfile_templates->Close();
  delete outfile_templates;
}
