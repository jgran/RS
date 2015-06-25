#include <string>
#include "TChain.h"
#include "TString.h"

void doAll(std::string input_dir = "/nfs-6/userdata/mt2/V00-00-02", std::string sample = "ttall_msdecays", std::string output_dir = ".") {

  //  gSystem->Load("libMiniFWLite.so");
  gSystem->Load("libSmearLooperCORE.so");
  gSystem->Load("libSmearLooper.so");
  gSystem->Load("libBabymakerTools.so");

  TChain *ch = new TChain("mt2"); 
 
  ch->Add(Form("%s/%s*.root",input_dir.c_str(),sample.c_str()));

  SmearLooper *looper = new SmearLooper();
  looper->loop(ch, output_dir + "/" + sample + ".root", sample); 
}
