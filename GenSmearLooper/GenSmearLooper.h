#ifndef GenSmearLooper_h
#define GenSmearLooper_h

// C++ includes
//#include <string>
//#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "Math/LorentzVector.h"

//MT2
#include "../MT2CORE/mt2tree.h"
#include "../MT2CORE/sigSelections.h"
#include "../MT2CORE/SR.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
using namespace mt2;
class GenSmearLooper {

 public:

  struct lepcand {
    float pt;
    float eta;
    float phi;
    int pdgId;
    float mt;
    bool isPFCand;
  };

  GenSmearLooper();
  ~GenSmearLooper();

  void setup();
  void SetSignalRegions();
  void loop(TChain* chain, std::string output_name = "test.root", std::string sample = "test.root");
  void fillHistosSRBase();
  void fillHistosInclusive();
  void fillHistosSignalRegion(const std::string& prefix = "", const std::string& suffix = "");
  void fillHistos(std::map<std::string, TH1*>& h_1d, 
		  const std::string& dir = "", const std::string& suffix = ""); 
 private:

  TFile * outfile_;
  mt2tree t;
  float evtweight_;
  int nlepveto_;
  float leppt_;
  float mt_;
  std::map<std::string, TH1*> h_1d_global;
  std::vector<SR> SRVec;
  std::vector<SR> InclusiveRegions;
  SR SRBase;
  SR SRNoCut;
  float met_pt;
  float mht_pt;
  float met_phi;
  float mt2;
  float ht;
  int nJet40;
  int nBJet20;
  float diffMetMht;
  float deltaPhiMin;
  float jet1_pt;
  float jet2_pt;
  
};

#endif

