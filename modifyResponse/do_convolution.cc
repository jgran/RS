#include "do_convolution.h"

void convol::do_convol(std::string modType, float modValue)
{

  TFile* outFile = new TFile(Form("response_templates_convolution_%s_plus_%f.root", modType.c_str(), modValue), "RECREATE");

  TCanvas *c1 = new TCanvas("c1", "", 750, 750);
  c1->Divide(2, 2);

  TH1F* h_CrystalBall_neg = new TH1F("h_CrystalBall_neg", "CrystalBall_neg", 2000, -10, 10);
  TH1F* h_CrystalBall_pos = new TH1F("h_CrystalBall_pos", "CrystalBall_pos", 2000, -10, 10);
  TH1F* h_CrystalBall_neg_and_pos = new TH1F("h_CrystalBall_neg_and_pos", "CrystalBall_neg_and_pos", 2000, -10, 10);
  TH1F* h_Convolution = new TH1F("h_Convolution", "Convolution", 1000, 0, 10);
  
  if(modType != "low_tail") delete h_CrystalBall_neg;
  if(modType != "high_tail") delete h_CrystalBall_pos;
  if(modType != "low_and_high_tail") delete h_CrystalBall_neg_and_pos;

  int counter = 0;
  TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_madgraph_PHYS14_fineBins_wideRange_eventVeto_newMatching035.root", "READ");
  std::string keep = "h_tot_JetAll_ResponsePt";
  TIter it(responseFile->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey *)it()) ) {
    if (strncmp (key->GetTitle(), keep.c_str(), keep.length()) != 0) continue;
    std::string hist_name = (key->GetTitle());
    TH1F* h = (TH1F*) responseFile->Get(TString(hist_name));
    if(h->GetEntries() < 5){
      std::cout << "Number of entries = " << h->GetEntries() << std::endl;
      std::cout << "skipping " << hist_name << std::endl;
      outFile->cd();
      h->Write(hist_name.c_str());
      continue;
    }

    counter++;
    if(counter > 5) break;

    h->Scale(1.0/h->Integral());

    float tail_integral_old = -999.0;
    float tail_integral_new = -999.0;
    bool done = false;
    float alpha_neg = 1.0;
    float alpha_pos = 1.0;
    float cb_sigma = h->GetRMS()/5;
    std::cout << "cb_sigma = " << cb_sigma << std::endl;

    float integral = 0;
    int integral_bin_neg_old = -1;
    int integral_bin_neg_new = -1;
    for(int i=1; i<=h->GetNbinsX(); i++){
      integral += h->GetBinContent(i);
      if(integral/h->Integral() > 0.05){
        integral_bin_neg_old = i;   
        break;
      }
    }
    integral_bin_neg_new = integral_bin_neg_old;

    integral = 0;
    int integral_bin_pos_old = -1;
    int integral_bin_pos_new = -1;
    for(int i=h->GetNbinsX(); i>=1; i--){
      integral += h->GetBinContent(i);
      if(integral/h->Integral() > 0.05){
        integral_bin_pos_old = i;
        break;
      }
    }
    integral_bin_pos_new = integral_bin_pos_old;

    if(modType == "low_tail"){

      bool bad_alpha_range = true;
      done = false;
      while(!done){

        TF1* fneg = new TF1("fneg", CrystalBall_neg, -10, 10, 5);
        fneg->SetParameters(alpha_neg, 3, 0, cb_sigma/2, 1);

        for(int i=1; i<=h_CrystalBall_neg->GetNbinsX()/2; i++){
          h_CrystalBall_neg->SetBinContent(i, fneg->Eval(-10.0 + i*0.01));
        }

        for(int i=1; i<=h_Convolution->GetNbinsX(); i++){
          float sum = 0.0;
          if(i <= 10){
            sum = h->GetBinContent(i);
          }
          else{
            for(int j=1; j<=h->GetNbinsX(); j++){
              sum += h->GetBinContent(j)*h_CrystalBall_neg->GetBinContent(1000 + i-j);
            }
          }
          h_Convolution->SetBinContent(i, sum);
        }
        h_Convolution->Scale(1.0/h_Convolution->Integral());//FIXME

        int shift = h->GetMaximumBin() - h_Convolution->GetMaximumBin();
        tail_integral_old = h->Integral(1, integral_bin_neg_old);
        tail_integral_new = h_Convolution->Integral(1, integral_bin_neg_new - shift);

        if(tail_integral_new > modValue*tail_integral_old){
          if(alpha_neg > 10){
            std::cout << "ALPHA > 10. Stopping." << std::endl;
            done = true;
          }
          else if(bad_alpha_range){
            alpha_neg += 1.0;
          }
          else done = true;
          continue;
        }

        bad_alpha_range = false;
        if(alpha_neg > 1.5) alpha_neg -= 0.1;
        else alpha_neg -= 0.05;

        if(alpha_neg < 0){
          std::cout << "ALPHA WENT NEGATIVE" << std::endl;
          break;
         }

        delete fneg;
      
      }
    }
    else if(modType == "high_tail"){

      bool bad_alpha_range = true;
      done = false;
      while(!done){

        TF1* fpos = new TF1("fpos", CrystalBall_pos, -10, 10, 5);
        fpos->SetParameters(alpha_pos, 3, 0, cb_sigma/2, 1);

        for(int i=h_CrystalBall_pos->GetNbinsX()/2; i<=h_CrystalBall_pos->GetNbinsX(); i++){
          h_CrystalBall_pos->SetBinContent(i, fpos->Eval(-10.0 + i*0.01));
        }

        for(int i=1; i<=h_Convolution->GetNbinsX(); i++){
          float sum = 0.0;
          if(i <= 10){
            sum = h->GetBinContent(i);
          }
          else{
            for(int j=1; j<=h->GetNbinsX(); j++){
              sum += h->GetBinContent(j)*h_CrystalBall_pos->GetBinContent(1000 + i-j);
            }
          }
          h_Convolution->SetBinContent(i, sum);
        }
        h_Convolution->Scale(1.0/h_Convolution->Integral());//FIXME

        int shift = h->GetMaximumBin() - h_Convolution->GetMaximumBin();
        tail_integral_old = h->Integral(integral_bin_pos_old, 10000);
        tail_integral_new = h_Convolution->Integral(integral_bin_pos_new - shift, 10000);

        if(tail_integral_new > modValue*tail_integral_old){
          if(alpha_pos > 10){
            std::cout << "ALPHA > 10. Stopping." << std::endl;
            done = true;
          }
          else if(bad_alpha_range){
            alpha_pos += 1.0;
          }
          else done = true;
          continue;
        }

        bad_alpha_range = false;
        if(alpha_pos > 1.5) alpha_pos -= 0.1;
        else alpha_pos -= 0.05;

        if(alpha_pos < 0){
          std::cout << "ALPHA WENT NEGATIVE" << std::endl;
          break;
         }

        delete fpos;
      
      }
    }
    else if(modType == "low_and_high_tail"){

      bool bad_alpha_range_neg = true;
      bool bad_alpha_range_pos = true;
      bool done = false;
      bool done_neg = false;
      bool done_pos = false;
      while(!done){

        TF1* fneg = new TF1("fneg", CrystalBall_neg, -10, 10, 5);
        fneg->SetParameters(alpha_neg, 3, 0, cb_sigma/2, 1);

        TF1* fpos = new TF1("fpos", CrystalBall_pos, -10, 10, 5);
        fpos->SetParameters(alpha_pos, 3, 0, cb_sigma/2, 1);

        for(int i=1; i<=h_CrystalBall_neg_and_pos->GetNbinsX()/2; i++){
          h_CrystalBall_neg_and_pos->SetBinContent(i, fneg->Eval(-10.0 + i*0.01));
        }
        for(int i=h_CrystalBall_neg_and_pos->GetNbinsX()/2; i<=h_CrystalBall_neg_and_pos->GetNbinsX(); i++){
          h_CrystalBall_neg_and_pos->SetBinContent(i, fpos->Eval(-10.0 + i*0.01));
        }

        for(int i=1; i<=h_Convolution->GetNbinsX(); i++){
          float sum = 0.0;
          if(i <= 10){
            sum = h->GetBinContent(i);
          }
          else{
            for(int j=1; j<=h->GetNbinsX(); j++){
              sum += h->GetBinContent(j)*h_CrystalBall_neg_and_pos->GetBinContent(1000 + i-j);
            }
          }
          h_Convolution->SetBinContent(i, sum);
        }
        h_Convolution->Scale(1.0/h_Convolution->Integral());//FIXME

        int shift = h->GetMaximumBin() - h_Convolution->GetMaximumBin();
        float tail_integral_neg_old = h->Integral(1, integral_bin_neg_old);
        float tail_integral_neg_new = h_Convolution->Integral(1, integral_bin_neg_new - shift);
        float tail_integral_pos_old = h->Integral(integral_bin_pos_old, 10000);
        float tail_integral_pos_new = h_Convolution->Integral(integral_bin_pos_new - shift, 10000);

        if(tail_integral_neg_new > modValue*tail_integral_neg_old){
          if(alpha_neg > 10){
            std::cout << "ALPHA > 10. Stopping." << std::endl;
            done_neg = true;
          }
          else if(bad_alpha_range_neg){
            alpha_neg += 1.0;
          }
          else done_neg = true;
        }
        else bad_alpha_range_neg = false;

        if(tail_integral_pos_new > modValue*tail_integral_pos_old){
          if(alpha_pos > 10){
            std::cout << "ALPHA > 10. Stopping." << std::endl;
            done_pos = true;
          }
          else if(bad_alpha_range_pos){
            alpha_pos += 1.0;
          }
          else done_pos = true;
        }
        else bad_alpha_range_pos = false;


        if(!done_neg){
          if(alpha_neg > 1.5) alpha_neg -= 0.1;
          else alpha_neg -= 0.05;
        }
        if(!done_pos){
          if(alpha_pos > 1.5) alpha_pos -= 0.1;
          else alpha_pos -= 0.05;
        }

        if(alpha_neg < 0){
          std::cout << "ALPHA WENT NEGATIVE" << std::endl;
          break;
         }
        if(alpha_pos < 0){
          std::cout << "ALPHA WENT NEGATIVE" << std::endl;
          break;
         }

        if(done_neg && done_pos) done = true;

        delete fneg;
        delete fpos;
      }
    }
    else if(modType == "core"){
      std::cout << " modType: " << modType << " not implemented yet!" << std::endl;
      return;
    }
    else{
      std::cout << " modType: " << modType << " not recognized!" << std::endl;
      return;
    }
    
    int shift = h->GetMaximumBin() - h_Convolution->GetMaximumBin();
    TH1F* h_Convolution_shifted = new TH1F("h_Convolution_shifted", hist_name.c_str(), 1000, 0, 10);
    for(int i=1; i<=h_Convolution_shifted->GetNbinsX(); i++){
      h_Convolution_shifted->SetBinContent(i, h_Convolution->GetBinContent(i - shift));
    }

    std::cout << "h->GetMaximumBin() = " << h->GetMaximumBin() << std::endl;
    std::cout << "h_Convolution_shifted->GetMaximumBin() = " << h_Convolution_shifted->GetMaximumBin() << std::endl;

    std::cout << "Lower tail old = " << h->Integral(1, integral_bin_neg_old) << std::endl;
    std::cout << "Upper tail old = " << h->Integral(integral_bin_pos_old, 10000) << std::endl;
    std::cout << "Lower tail new = " << h_Convolution_shifted->Integral(1, integral_bin_neg_new) << std::endl;
    std::cout << "Upper tail new = " << h_Convolution_shifted->Integral(integral_bin_pos_new, 10000) << std::endl;

    c1->cd(1);
    gPad->SetLogy();
    if(modType == "low_tail"){
      h_CrystalBall_neg->SetMinimum(0.00001);
      h_CrystalBall_neg->Draw();
    }
    else if(modType == "low_and_high_tail"){
      h_CrystalBall_neg_and_pos->SetMinimum(0.00001);
      h_CrystalBall_neg_and_pos->Draw();
    }

    c1->cd(2);
    gPad->SetLogy();
    h->SetMinimum(0.00001);
    h->Draw();

    c1->cd(3);
    gPad->SetLogy();
    if(modType == "high_tail"){
      h_CrystalBall_pos->SetMinimum(0.00001);
      h_CrystalBall_pos->Draw();
    }

    c1->cd(4);
    gPad->SetLogy();
    h_Convolution_shifted->SetMinimum(0.00001);
    h_Convolution_shifted->Draw();

    c1->Print(Form("convolutions_low_tail_plus_100percent/low_tail_plus_100percent_%s.pdf", hist_name.c_str()));
    outFile->cd();
    h_Convolution_shifted->Write(hist_name.c_str());

    delete h_Convolution_shifted;
    
  }

  delete h_Convolution;
  if(modType == "low_tail") delete h_CrystalBall_neg;
  if(modType == "high_tail") delete h_CrystalBall_pos;
  if(modType == "low_and_high_tail") delete h_CrystalBall_neg_and_pos;

  outFile->Close();
  delete outFile;
  
}
