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

void findTail(){

  TFile* outfile = new TFile("template_fits_MC.root", "RECREATE");
  //TFile* outfile = new TFile("testing.root", "RECREATE");
  //TFile* responseFile = new TFile("/home/users/jgran/response_templates_fitting.root", "READ");
  //TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_madgraph_PHYS14_fineBins_wideRange_eventVeto_newMatching035.root", "READ");
  //TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_MGMLM_Spring15_bestMatching_angles_withNeutrinos.root", "READ");
  TFile* responseFile = new TFile("input_templates.root", "READ");
  TCanvas* c1 = new TCanvas();
  TIter it(responseFile->GetListOfKeys());
  TKey* key;
  std::string keep = "h_tot_JetAll_ResponsePt";
  //std::string keep = "h_tot_JetAll_ResponsePt_Pt6_Eta9";
  //std::string keep = "h_tot_JetAll_ResponsePt_Pt0_Eta0";
  while ( (key = (TKey *)it()) ) {
    if (strncmp (key->GetTitle(), keep.c_str(), keep.length()) != 0) continue;
    std::string hist_name = (key->GetTitle());
    TH1F* h = (TH1F*) responseFile->Get(TString(hist_name));
    if(h->GetEntries() < 100) continue;
    h->Scale(1.0/h->Integral());
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
    fit->Print();
    fit->Update();
    c1->cd();
    c1->Print(Form("fits/fit_%s.pdf", hist_name.c_str()));
    outfile->cd();
    fit->Write();
    responseFile->cd();

    double g_scale = fit->GetParameter(0);
    double g_mean  = fit->GetParameter(1);
    double g_sigma = fit->GetParameter(2);

    TH1F* h_tail = (TH1F*) h->Clone();
    TH1F* h_core = (TH1F*) h->Clone();
    for(int i=1; i<=h->GetNbinsX(); i++){
        float gaussian_contribution = ComputeGaussianValue(g_scale, g_mean, g_sigma, h->GetBinLowEdge(i));
        float tail_contribution = h->GetBinContent(i) - gaussian_contribution;
        h_tail->SetBinContent(i, std::max(tail_contribution, 0)); 
        h_core->SetBinContent(i, gaussian_contribution); 
    }
    h_tail->SetLineColor(kGreen);
    h_tail->Draw();
    h_core->SetLineColor(kRed);
    h_core->Draw("same");
    h->Draw("same");
    c1->Print(Form("plots/core_tail_%s.pdf", hist_name.c_str()));

    delete fit;
  }
  responseFile->Close();
  delete responseFile;
  delete c1;

  outfile->Close();
  delete outfile;
}
