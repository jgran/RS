#include <string>
#include "TChain.h"
#include "TString.h"

void doAll(std::string input_dir = "/nfs-6/userdata/mt2/V00-00-02", std::string sample = "ttall_msdecays", std::string output_dir = ".") {

  //  gSystem->Load("libMiniFWLite.so");
  gSystem->Load("libRebalanceLooperCORE.so");
  gSystem->Load("libRebalanceLooper.so");
  gSystem->Load("libBabymakerTools.so");

  TChain *ch = new TChain("mt2"); 
 
  ch->Add(Form("%s/%s*.root",input_dir.c_str(),sample.c_str()));

  RebalanceLooper *looper = new RebalanceLooper();
  //looper->loop(ch, output_dir + "/" + "rebalance.root"); 
  looper->loop(ch, Form("%s/%s.root", output_dir.c_str(), sample.c_str())); 
}
