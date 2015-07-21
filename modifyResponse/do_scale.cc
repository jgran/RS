#include "do_scale.h"

void scale::do_scale(std::string modType, float modValue)
{

  TFile* outFile = new TFile(Form("response_templates_scale_%s_plus_%f.root", modType.c_str(), modValue), "RECREATE");

  TCanvas *c1 = new TCanvas("c1", "", 300, 400);
  c1->Divide(1, 2);

  int counter = 0;
  TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_madgraph_PHYS14_fineBins_wideRange_eventVeto_newMatching035.root", "READ");
  std::string keep = "h_tot_JetAll_ResponsePt";
  TIter it(responseFile->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey *)it()) ) {
    if (strncmp (key->GetTitle(), keep.c_str(), keep.length()) != 0) continue;
    std::string hist_name = (key->GetTitle());
    TH1F* h = (TH1F*) responseFile->Get(TString(hist_name));
    if(h->GetEntries() < 20){
      std::cout << "Number of entries = " << h->GetEntries() << std::endl;
      std::cout << "skipping " << hist_name << std::endl;
      outFile->cd();
      h->Write(hist_name.c_str());
      continue;
    }
    h->Scale(1.0/h->Integral());

    float tail_integral_old = -999.0;
    float tail_integral_new = -999.0;

    float tail_def = 0.05;
    float sf_tail = 1.5;

    counter++;
    if(counter > 5) break;

    float integral_low = 0;
    int integral_bin_low = -1;
    for(int i=1; i<=h->GetNbinsX(); i++){
      integral_low += h->GetBinContent(i);
      if(integral_low/h->Integral() > tail_def){
        integral_bin_low = i;   
        break;
      }
    }

    float integral_high = 0;
    int integral_bin_high = -1;
    for(int i=h->GetNbinsX(); i>=1; i--){
      integral_high += h->GetBinContent(i);
      if(integral_high/h->Integral() > tail_def){
        integral_bin_high= i;   
        break;
      }
    }

    TH1F *h_result = (TH1F*)h->Clone("h_result");

    if(modType == "low_tail"){
      float sf_bulk = (1 - integral_low*sf_tail)/(1 - integral_low);
      for(int i=1; i<h_result->GetNbinsX(); i++){
        if(i<=integral_bin_low) h_result->SetBinContent(i, sf_tail*h_result->GetBinContent(i));
        else h_result->SetBinContent(i, sf_bulk*h_result->GetBinContent(i));
      }
    } else if (modType == "high_tail"){
      float sf_bulk = (1 - integral_high*sf_tail)/(1 - integral_high);
      for(int i=1; i<h_result->GetNbinsX(); i++){
        if(i>=integral_bin_high) h_result->SetBinContent(i, sf_tail*h_result->GetBinContent(i));
        else h_result->SetBinContent(i, sf_bulk*h_result->GetBinContent(i));
      }
    } else if(modType == "low_and_high_tail"){
      float sf_bulk = (1 - integral_low*sf_tail - integral_high*sf_tail)/(1 - integral_low - integral_high);
      for(int i=1; i<h_result->GetNbinsX(); i++){
        if(i<=integral_bin_low) h_result->SetBinContent(i, sf_tail*h_result->GetBinContent(i));
        else if(i>=integral_bin_high) h_result->SetBinContent(i, sf_tail*h_result->GetBinContent(i));
        else h_result->SetBinContent(i, sf_bulk*h_result->GetBinContent(i));
      }
      h_result->Scale(1.0/h_result->Integral());
    } else std::cout << "unknown modType: " << modType << std::endl;

    //std::cout << "Lower tail old = " << h->Integral(1, integral_bin_low) << std::endl;
    //std::cout << "Lower tail new = " << h_result->Integral(1, integral_bin_low) << std::endl;
    //std::cout << "Upper tail old = " << h->Integral(integral_bin_high, 10000) << std::endl;
    //std::cout << "Upper tail new = " << h_result->Integral(integral_bin_high, 10000) << std::endl;

    c1->cd(1);
    gPad->SetLogy();
    h->SetMinimum(0.00001);
    h->Draw();

    c1->cd(2);
    gPad->SetLogy();
    h_result->SetMinimum(0.00001);
    h_result->Draw();

    c1->Print(Form("convolutions_simpleScale_low_and_high_tail_plus_50percent/simpleScale_high_tail_plus_50percent_%s.pdf", hist_name.c_str()));
    outFile->cd();
    h_result->Write(hist_name.c_str());

    delete h_result;
    
  }

  outFile->Close();
  delete outFile;
  
}
