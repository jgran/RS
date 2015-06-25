#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

//fits with three fucking gaussians
Double_t TripleGaussian(Double_t *x, Double_t *par){
  return par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])) + par[3]*TMath::Exp(-0.5*((x[0]-par[4])/par[5])*((x[0]-par[4])/par[5])) + par[6]*TMath::Exp(-0.5*((x[0]-par[7])/par[8])*((x[0]-par[7])/par[8]));
}

Double_t SingleGaussian(Double_t *x, Double_t *par){
  return par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));
}

void ResponseFitter(){

  TFile* outfile = new TFile("template_fits.root", "RECREATE");
  //TFile* responseFile = new TFile("/home/users/jgran/response_templates_fitting.root", "READ");
  TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_madgraph_PHYS14_fineBins_wideRange_eventVeto_newMatching035.root", "READ");
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
    if(h->Integral() < 0.001) continue;
    h->Scale(1.0/h->Integral());
    //TF1* fit = new TF1("fit", TripleGaussian, (h->GetMinimumBin())*0.01, (h->GetMaximumBin())*0.01, 9);
    //TF1* fit = new TF1(Form("fit_%s", hist_name.c_str()), TripleGaussian, (h->GetMinimumBin())*0.01, (h->GetMaximumBin())*0.01, 9);
    //TF1* fit = new TF1(Form("fit_%s", hist_name.c_str()), TripleGaussian, 0, 10, 9);
    TF1* fit = new TF1(Form("fit_%s", hist_name.c_str()), SingleGaussian, 0.2, 2.2, 3);
    fit->SetTitle(Form("fit_%s", hist_name.c_str()));
    float mean = h->GetMaximumBin();
    float max = h->GetBinContent(mean);
    mean *= 0.01;//multiply bin by bin width
/*
    fit->SetParameters(0.00001, mean/1.2, 0.05, max, mean, 0.01, 0.00001, mean*1.2, 0.05);
    fit->FixParameter(3, max);//fix central gaussian peak height
    fit->FixParameter(4, mean);//fix central gaussian mean
    fit->SetParLimits(0, 0.0, 1000);
    fit->SetParLimits(1, 0.1, 1.0);
    fit->SetParLimits(2, 0.0, 1000);
    fit->SetParLimits(6, 0.0, 1000);
    fit->SetParLimits(7, 1.0, 3.0);
    fit->SetParLimits(8, 0.0, 1000);
*/
    fit->SetParameters(max, mean, 0.01);
    fit->FixParameter(0, max);//fix central gaussian peak height
    //fit->FixParameter(1, mean);//fix central gaussian mean
    h->Fit(Form("fit_%s", hist_name.c_str()));
    //fit->SetParameter(0, 0);//set first gaussian to zero
    //fit->SetParameter(6, 0);//set third gaussian to zero
    fit->SetRange(0, 10);
    //fit->SetLineColor(kBlue);
    //fit->Draw("SAME");
    fit->Print();
    fit->Update();
    c1->cd();
    c1->Print(Form("fits/fit_%s.pdf", hist_name.c_str()));
    outfile->cd();
    fit->Write();
    responseFile->cd();
    delete fit;
  }
  responseFile->Close();
  delete responseFile;
  delete c1;

  outfile->Close();
  delete outfile;
}
