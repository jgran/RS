#include "TH2.h"

void MakeScaleFactor(){

  TFile* outfile = new TFile("rebalance_scale_factor.root", "RECREATE");
  //TFile* infile = new TFile("/home/users/jgran/temp/V00-00-12/MT2Analysis/SmearLooper/output/test_fit/RS_qcd_pt.root", "READ");
  TFile* infile = new TFile("/home/users/jgran/temp/V00-00-12/MT2Analysis/SmearLooper/output/mht_minimization/RS_qcd_pt.root", "READ");
  TH2F* h = (TH2F*) infile->Get("h2d_recopt_RebalancedOverGen");
  
  TH1F* sf = new TH1F("sf", "sf", 100, 0, 1000);
  //TH2F* sf = new TH2F("sf", "sf", 99, 10, 1000, 40, 0.9, 1.3);

  for(int i=1; i<=100; i++){
    float entry = h->ProjectionY("bin",i,i)->GetMean();
    sf->Fill((i-1)*10, entry);
    //sf->Fill(i*10, entry, 1.0);
  }

  outfile->cd();
  sf->Write();
  outfile->Close();
  delete outfile;
  delete h;
  infile->Close();
  delete infile;

}
