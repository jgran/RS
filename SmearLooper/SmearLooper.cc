// C++
#include <iostream>
#include <vector>
#include <set>
#include <cmath>

// ROOT
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TList.h"
#include "TCollection.h"
#include "TKey.h"
#include "TLorentzVector.h"

#include "../Tools/utils.h"
#include "../Tools/hemJet.h"
#include "../Tools/MT2/MT2.h"

// header
#include "SmearLooper.h"

//MT2
#include "../MT2CORE/Plotting/PlotUtilities.h"


using namespace std;
using namespace mt2;

class mt2tree;
class SR;

// mt2 binning for results
const int n_mt2bins = 5;
const float mt2bins[n_mt2bins+1] = {200., 300., 400., 600., 1000., 1500.};
const std::string mt2binsname[n_mt2bins+1] = {"200", "300", "400", "600", "1000", "1500"};
bool useDRforGammaQCDMixing = true; // requires GenParticles

//for rebalance and smear
const int numberOfSmears = 100;
std::vector<float> jet_pt;
std::vector<float> jet_eta;
std::vector<float> jet_phi;
std::vector<float> jet_btagCSV;
std::vector<float> PU_passes_id_jet_pt;
std::vector<float> PU_passes_id_jet_eta;
std::vector<float> PU_passes_id_jet_phi;
std::vector<float> PU_passes_id_jet_btagCSV;
std::vector<float> PU_fails_id_jet_pt;
std::vector<float> PU_fails_id_jet_eta;
std::vector<float> PU_fails_id_jet_phi;
float pt_soft_x;
float pt_soft_y;
std::map<std::string, std::vector<float> > responseMap;

inline bool sortByPt(const LorentzVector &vec1, const LorentzVector &vec2 ) {
    return vec1.pt() > vec2.pt();
}


void SmearLooper::setup(){
  //TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_madgraph_PHYS14_fineBins_wideRange_eventVeto_newJetID.root", "READ");
  TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_madgraph_PHYS14_fineBins_wideRange_eventVeto_newMatching035.root", "READ");
  TIter it(responseFile->GetListOfKeys());
  TKey* key;
  std::string keep = "h_tot_JetAll_ResponsePt";
  while ((key = (TKey *)it())) {
    if (strncmp (key->GetTitle(), keep.c_str(), keep.length()) != 0) continue;
    std::string hist_name = (key->GetTitle());
    //std::cout << hist_name << std::endl;
    TH1F* h = (TH1F*) responseFile->Get(TString(hist_name));
    //if(h->Integral() > 0) h->Scale(1.0/h->Integral());
    h->Scale(1.0/h->Integral());
    std::vector<float> temp_vec;
    for(int i=1; i<=h->GetNbinsX(); i++){
      temp_vec.push_back(h->GetBinContent(i));
    } 
    responseMap[hist_name] = temp_vec;
  }
  responseFile->Close();
  delete responseFile;
}

/*
float GetValue(float xvalue, TH1F &hist){
  float xmax = hist.GetXaxis()->GetXmax()-0.01;
  int xbin   = hist.GetXaxis()->FindBin(min(xvalue, xmax));
  return hist.GetBinContent(xbin);
}
*/

std::vector<float> GetResponseVector(float pt, float eta){

  std::string pt_string;
  if(pt < 20.0) pt_string = "Pt0";
  else if(pt < 30.0)   pt_string = "Pt1";
  else if(pt < 50.0)   pt_string = "Pt2";
  else if(pt < 80.0)   pt_string = "Pt3";
  else if(pt < 120.0)  pt_string = "Pt4";
  else if(pt < 170.0)  pt_string = "Pt5";
  else if(pt < 230.0)  pt_string = "Pt6";
  else if(pt < 300.0)  pt_string = "Pt7";
  else if(pt < 380.0)  pt_string = "Pt8";
  else if(pt < 470.0)  pt_string = "Pt9";
  else if(pt < 570.0)  pt_string = "Pt10";
  else if(pt < 680.0)  pt_string = "Pt11";
  else if(pt < 800.0)  pt_string = "Pt12";
  else pt_string = "Pt13";

/*
  else if(pt < 1000.0) pt_string = "Pt13";
  else if(pt < 1300.0) pt_string = "Pt14";
  else if(pt < 1700.0) pt_string = "Pt15";
  else if(pt < 2200.0) pt_string = "Pt16";
  else if(pt < 2800.0) pt_string = "Pt17";
  else if(pt < 3500.0) pt_string = "Pt18";
  else if(pt < 4300.0) pt_string = "Pt19";
  else if(pt < 5200.0) pt_string = "Pt20";
  else if(pt < 6500.0) pt_string = "Pt21";
*/

  std::string eta_string;
  if(eta < 0.3) eta_string = "Eta0";
  else if(eta < 0.5)   eta_string = "Eta1";
  else if(eta < 0.8)   eta_string = "Eta2";
  else if(eta < 1.1)   eta_string = "Eta3";
  else if(eta < 1.4)  eta_string = "Eta4";
  else if(eta < 1.7)  eta_string = "Eta5";
  else if(eta < 2.3)  eta_string = "Eta6";
  else if(eta < 2.8)  eta_string = "Eta7";
  else eta_string = "Eta8";

  //else if(eta < 3.2)  eta_string = "Eta8";
  //else if(eta < 4.1)  eta_string = "Eta9";
  //else if(eta < 5.0)  eta_string = "Eta10";
  //else if(eta < 0.5)  eta_string = "Eta11"; //looks like the binning in the email might not be right

  //std::cout << Form("h_tot_JetAll_ResponsePt_%s_%s", pt_string.c_str(), eta_string.c_str()) << std::endl;

  std::string hist_name = "h_tot_JetAll_ResponsePt_" + pt_string + "_" + eta_string;
  return responseMap.at(hist_name);
}


SmearLooper::SmearLooper(){
}
SmearLooper::~SmearLooper(){

};

void SmearLooper::SetSignalRegions(){

  //SRVec = getSignalRegions2015LowLumi(); //Phys14 AN selection
  //SRVec =  getSignalRegions2015SevenJets_UltraHighHT(); //new selection with additional njet boundary at 7 jets, 4th ht bin
  //SRVec =  getSignalRegionsZurich(); //same as getSignalRegions2015SevenJets_UltraHighHT(), but with minMT binning removed
  SRVec =  getSignalRegionsZurich_jetpt40(); //same as getSignalRegionsZurich(), but with j1pt and j2pt cuts changed to 40 GeV

  //store histograms with cut values for all variables
  for(unsigned int i = 0; i < SRVec.size(); i++){
    std::vector<std::string> vars = SRVec.at(i).GetListOfVariables();
    TDirectory * dir = (TDirectory*)outfile_->Get(("sr"+SRVec.at(i).GetName()).c_str());
    if (dir == 0) {
      dir = outfile_->mkdir(("sr"+SRVec.at(i).GetName()).c_str());
    } 
    dir->cd();
    for(unsigned int j = 0; j < vars.size(); j++){
      plot1D("h_"+vars.at(j)+"_"+"LOW",  1, SRVec.at(i).GetLowerBound(vars.at(j)), SRVec.at(i).srHistMap, "", 1, 0, 2);
      plot1D("h_"+vars.at(j)+"_"+"UP",   1, SRVec.at(i).GetUpperBound(vars.at(j)), SRVec.at(i).srHistMap, "", 1, 0, 2);
    }

    dir = (TDirectory*)outfile_->Get(("crsl"+SRVec.at(i).GetName()).c_str());
    if (dir == 0) {
      dir = outfile_->mkdir(("crsl"+SRVec.at(i).GetName()).c_str());
    } 
    dir->cd();
    for(unsigned int j = 0; j < vars.size(); j++){
      plot1D("h_"+vars.at(j)+"_"+"LOW",  1, SRVec.at(i).GetLowerBound(vars.at(j)), SRVec.at(i).crslHistMap, "", 1, 0, 2);
      plot1D("h_"+vars.at(j)+"_"+"UP",   1, SRVec.at(i).GetUpperBound(vars.at(j)), SRVec.at(i).crslHistMap, "", 1, 0, 2);
    }

    dir = (TDirectory*)outfile_->Get(("crgj"+SRVec.at(i).GetName()).c_str());
    if (dir == 0) {
      dir = outfile_->mkdir(("crgj"+SRVec.at(i).GetName()).c_str());
    } 
    dir->cd();
    for(unsigned int j = 0; j < vars.size(); j++){
      plot1D("h_"+vars.at(j)+"_"+"LOW",  1, SRVec.at(i).GetLowerBound(vars.at(j)), SRVec.at(i).crgjHistMap, "", 1, 0, 2);
      plot1D("h_"+vars.at(j)+"_"+"UP",   1, SRVec.at(i).GetUpperBound(vars.at(j)), SRVec.at(i).crgjHistMap, "", 1, 0, 2);
    }
    outfile_->cd();
  }

  SRBase.SetName("srbase");
  SRBase.SetVar("mt2", 200, -1);
  SRBase.SetVar("j1pt", 40, -1);
  SRBase.SetVar("j2pt", 40, -1);
  SRBase.SetVar("deltaPhiMin", 0.3, -1);
  SRBase.SetVar("diffMetMhtOverMet", 0, 0.5);
  SRBase.SetVar("nlep", 0, 1);
  SRBase.SetVar("passesHtMet", 1, 2);

  std::vector<std::string> vars = SRBase.GetListOfVariables();
  TDirectory * dir = (TDirectory*)outfile_->Get((SRBase.GetName()).c_str());
  if (dir == 0) {
    dir = outfile_->mkdir((SRBase.GetName()).c_str());
  } 
  dir->cd();
  for(unsigned int j = 0; j < vars.size(); j++){
    plot1D("h_"+vars.at(j)+"_"+"LOW",  1, SRBase.GetLowerBound(vars.at(j)), SRBase.srHistMap, "", 1, 0, 2);
    plot1D("h_"+vars.at(j)+"_"+"UP",   1, SRBase.GetUpperBound(vars.at(j)), SRBase.srHistMap, "", 1, 0, 2);
  }
  outfile_->cd();

  //setup inclusive regions
  SR InclusiveHT450to575 = SRBase;
  InclusiveHT450to575.SetName("InclusiveHT450to575");
  InclusiveHT450to575.SetVar("ht", 450, 575);
  InclusiveRegions.push_back(InclusiveHT450to575);

  SR InclusiveHT575to1000 = SRBase;
  InclusiveHT575to1000.SetName("InclusiveHT575to1000");
  InclusiveHT575to1000.SetVar("ht", 575, 1000);
  InclusiveRegions.push_back(InclusiveHT575to1000);

  SR InclusiveHT1000to1500 = SRBase;
  InclusiveHT1000to1500.SetName("InclusiveHT1000to1500");
  InclusiveHT1000to1500.SetVar("ht", 1000, 1500);
  InclusiveRegions.push_back(InclusiveHT1000to1500);

  SR InclusiveHT1500toInf = SRBase;
  InclusiveHT1500toInf.SetName("InclusiveHT1500toInf");
  InclusiveHT1500toInf.SetVar("ht", 1500, -1);
  InclusiveRegions.push_back(InclusiveHT1500toInf);

  SR InclusiveNJets2to3 = SRBase;
  InclusiveNJets2to3.SetName("InclusiveNJets2to3");
  InclusiveNJets2to3.SetVar("njets", 2, 4);
  InclusiveRegions.push_back(InclusiveNJets2to3);

  SR InclusiveNJets4to6 = SRBase;
  InclusiveNJets4to6.SetName("InclusiveNJets4to6");
  InclusiveNJets4to6.SetVar("njets", 4, 7);
  InclusiveRegions.push_back(InclusiveNJets4to6);

  SR InclusiveNJets7toInf = SRBase;
  InclusiveNJets7toInf.SetName("InclusiveNJets7toInf");
  InclusiveNJets7toInf.SetVar("njets", 7, -1);
  InclusiveRegions.push_back(InclusiveNJets7toInf);

  SR InclusiveNBJets0 = SRBase;
  InclusiveNBJets0.SetName("InclusiveNBJets0");
  InclusiveNBJets0.SetVar("nbjets", 0, 1);
  InclusiveRegions.push_back(InclusiveNBJets0);

  SR InclusiveNBJets1 = SRBase;
  InclusiveNBJets1.SetName("InclusiveNBJets1");
  InclusiveNBJets1.SetVar("nbjets", 1, 2);
  InclusiveRegions.push_back(InclusiveNBJets1);

  SR InclusiveNBJets2 = SRBase;
  InclusiveNBJets2.SetName("InclusiveNBJets2");
  InclusiveNBJets2.SetVar("nbjets", 2, 3);
  InclusiveRegions.push_back(InclusiveNBJets2);

  SR InclusiveNBJets3toInf = SRBase;
  InclusiveNBJets3toInf.SetName("InclusiveNBJets3toInf");
  InclusiveNBJets3toInf.SetVar("nbjets", 3, -1);
  InclusiveRegions.push_back(InclusiveNBJets3toInf);

  for(unsigned int i=0; i<InclusiveRegions.size(); i++){
    TDirectory * dir = (TDirectory*)outfile_->Get((InclusiveRegions.at(i).GetName()).c_str());
    if (dir == 0) {
      dir = outfile_->mkdir((InclusiveRegions.at(i).GetName()).c_str());
    } 
  }

}


void SmearLooper::loop(TChain* chain, std::string output_name, std::string sample){

  setup();

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");


  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/all/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/jetpt10/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/jetpt10_PUid/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/jetpt10_PUid_200k/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/jetpt10_PUid_all/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test_simplex/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test_gaussian/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test_fit/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/10k_test_fit/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/mht_minimization/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test_sigmasoft5/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test_sigmasoft5_10percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test_sigmasoft5_50percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/test_sigmasoft5_50percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft20_2percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft5_2percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft30_2percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft10_2percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft15_2percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft_nvert_2percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft_nvert_20percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft_sqrtnvert_10percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft20_100percent/%s.root", sample.c_str()), "READ");
  //TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft15_100percent/%s.root", sample.c_str()), "READ");
  TFile* rebalanceFile = new TFile(Form("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/output/met_minimization_sigmasoft25_100percent/%s.root", sample.c_str()), "READ");

  TTree *rebalanceTree = (TTree*)rebalanceFile->Get("rebalance");
  if(rebalanceTree) std::cout << "Good rebalance tree" << std::endl;
  if(!rebalanceTree) std::cout << "Bad rebalance tree" << std::endl;
  std::vector<float>* rebalanceFactors = 0;
  rebalanceTree->SetBranchAddress("rebalanceFactors", &rebalanceFactors);
  std::vector<int>* useJet = 0;
  rebalanceTree->SetBranchAddress("useJet", &useJet);


  gROOT->cd();
  cout << "[SmearLooper::loop] creating output file: " << output_name << endl;

  outfile_ = new TFile(output_name.c_str(),"RECREATE") ; 

  cout << "[SmearLooper::loop] setting up histos" << endl;

  SetSignalRegions();

  SRNoCut.SetName("nocut");

  // File Loop
  int nDuplicates = 0;
  int nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  cout << "[SmearLooper::loop] running on " << nEventsChain << " events" << endl;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  int rebalanceCounter = 0;
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    cout << "[SmearLooper::loop] running on file: " << currentFile->GetTitle() << endl;

    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("mt2");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);

    t.Init(tree);

    TRandom3 *random = new TRandom3();

    float pass_met = 0;
    float fail_met = 0;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < (nEventsTree/10); ++event) {
      t.GetEntry(event);
      // remove low pt QCD samples 
      if (t.evt_id >= 100 && t.evt_id < 109) continue;
      // remove low HT QCD samples 
      if (t.evt_id >= 120 && t.evt_id < 123) continue;
      if (t.nVert == 0) continue;
      if (t.njet < 2) continue;
      if (t.met_pt <= 20.0) pass_met++;
      else fail_met++;
    }
    float met_eff = pass_met/(pass_met + fail_met);
    std::cout << "pass_met = " << pass_met << std::endl;
    std::cout << "fail_met = " << fail_met << std::endl;
    std::cout << "met_eff = " << met_eff << std::endl;

    // Event Loop
    //unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    //for( unsigned int event = 0; event < (nEventsTree/50); ++event) {
    //for( unsigned int event = 0; event < 10000; ++event) {
    //for( unsigned int event = 0; event < 200000; ++event) {

      t.GetEntry(event);
      rebalanceTree->GetEntry(event);

      //---------------------
      // bookkeeping and progress report
      //---------------------
      ++nEventsTotal;
      if (nEventsTotal%10000==0) {
        ULong64_t i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
        if (isatty(1)) {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
              "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
          fflush(stdout);
        }
        else {
          cout<<i_permille/10.<<" ";
          if (nEventsTotal%100000==0) cout<<endl;
        }
      }

      //---------------------
      // skip duplicates -- will need this eventually
      //---------------------
      //if( isData ) {
      //  DorkyEventIdentifier id = {stopt.run(), stopt.event(), stopt.lumi() };
      //  if (is_duplicate(id, already_seen) ){
      //    continue;
      //  }
      //}
      //

      //---------------------
      // basic event selection and cleaning
      //---------------------
      // remove low pt QCD samples 
      if (t.evt_id >= 100 && t.evt_id < 109) continue;
      // remove low HT QCD samples 
      if (t.evt_id >= 120 && t.evt_id < 123) continue;

      if (t.nVert == 0) continue;
      if (t.njet < 2) continue;

/*
      int rebalanced_jet_counter = 0;
      for(int i=0; i<t.njet; i++){
        if( (t.jet_puId[i] < 1) && (t.jet_pt[i] < 100.0) ) continue;
        else rebalanced_jet_counter++;
      }
      if(rebalanced_jet_counter < 2) continue;
*/

/*
      if(rebalanceFactors->size() != t.njet){
        std::cout << "rebalance vector mismatch!" << std::endl;
        std::cout << "size of rebalance vector = " << rebalanceFactors->size() << std::endl;
        std::cout << "t.njet = " << t.njet << std::endl;
        return;
      }
*/


      //---------------------
      // set weights and start making plots
      //---------------------
      outfile_->cd();
      const float lumi = 4.;
      const float smearNormalization = 1.0/numberOfSmears;
      evtweight_ = t.evt_scale1fb * lumi;


      jet_pt.clear();
      jet_eta.clear();
      jet_phi.clear();
      jet_btagCSV.clear();
      PU_passes_id_jet_pt.clear();
      PU_passes_id_jet_eta.clear();
      PU_passes_id_jet_phi.clear();
      PU_passes_id_jet_btagCSV.clear();
      PU_fails_id_jet_pt.clear();
      PU_fails_id_jet_eta.clear();
      PU_fails_id_jet_phi.clear();

      float new_met_x = t.met_pt*cos(t.met_phi);
      float new_met_y = t.met_pt*sin(t.met_phi);

      float pt_soft_true_x = 0.0;
      float pt_soft_true_y = 0.0;
      float pt_soft_reco_x = - t.met_pt*cos(t.met_phi);
      float pt_soft_reco_y = - t.met_pt*sin(t.met_phi);
      bool dojet = true;//FIXME
      int jetCounter = -1;//FIXME
      float rf = -999;
      for(int i=0; i<t.njet; i++){
        //if(t.jet_puId[i] > 0) dojet = true;

        new_met_x += t.jet_pt[i]*cos(t.jet_phi[i]);//FIXME
        new_met_y += t.jet_pt[i]*sin(t.jet_phi[i]);//FIXME

        pt_soft_reco_x -= t.jet_pt[i]*cos(t.jet_phi[i]);
        pt_soft_reco_y -= t.jet_pt[i]*sin(t.jet_phi[i]);
        
/*
        if( (t.jet_puId[i] < 1) && (t.jet_pt[i] < 100.0) ){
          dojet = false;//FIXME
          PU_jet_pt.push_back(t.jet_pt[i]);
          PU_jet_eta.push_back(t.jet_eta[i]);
          PU_jet_phi.push_back(t.jet_phi[i]);
          PU_jet_btagCSV.push_back(t.jet_btagCSV[i]);
        }
*/
        if((t.jet_pt[i] < 100) && (t.jet_puId[i] < 1)){
          dojet = false;//FIXME
          if(t.jet_id[i] > 0){
            PU_passes_id_jet_pt.push_back(t.jet_pt[i]);
            PU_passes_id_jet_eta.push_back(t.jet_eta[i]);
            PU_passes_id_jet_phi.push_back(t.jet_phi[i]);
            PU_passes_id_jet_btagCSV.push_back(t.jet_btagCSV[i]);
          } else {
            PU_fails_id_jet_pt.push_back(t.jet_pt[i]);
            PU_fails_id_jet_eta.push_back(t.jet_eta[i]);
            PU_fails_id_jet_phi.push_back(t.jet_phi[i]);
          }
        }
        else dojet = true;//FIXME
        
        if(!dojet) continue;

        jetCounter++;
        rf = rebalanceFactors->at(jetCounter);//FIXME
        jet_pt.push_back(t.jet_pt[i]*rf);//FIXME
        //jet_pt.push_back(t.jet_pt[i]);//FIXME
        jet_eta.push_back(t.jet_eta[i]);
        jet_phi.push_back(t.jet_phi[i]);
        jet_btagCSV.push_back(t.jet_btagCSV[i]);

        new_met_x -= jet_pt.at(jetCounter)*cos(jet_phi.at(jetCounter));//FIXME
        new_met_y -= jet_pt.at(jetCounter)*sin(jet_phi.at(jetCounter));//FIXME
        pt_soft_true_x -= jet_pt.at(jetCounter)*cos(jet_phi.at(jetCounter));
        pt_soft_true_y -= jet_pt.at(jetCounter)*sin(jet_phi.at(jetCounter));

        if(dojet) {
          plot1DUnderOverFlow("h_rebalance_factor",       1.0/rf,       evtweight_, h_1d_global, ";rebalance factor", 1000, 0, 10);
        }
        if(t.jet_mcPt[i] > 1){
          if(dojet){
            plot1DUnderOverFlow("h_RecoptOverGenpt", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reco jet p_{T}/gen jet p_{T}", 100, 0, 10);
            if(t.jet_pt[i] < 20.0) plot1DUnderOverFlow("h_RecoptOverGenpt_10to20", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt10to20)", 100, 0, 10);
            if(t.jet_pt[i] < 50.0) plot1DUnderOverFlow("h_RecoptOverGenpt_20to50", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt20to50)", 100, 0, 10);
            else if(t.jet_pt[i] < 100.0) plot1DUnderOverFlow("h_RecoptOverGenpt_50to100", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt50to100)", 100, 0, 10);
            else if(t.jet_pt[i] < 200.0) plot1DUnderOverFlow("h_RecoptOverGenpt_100to200", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt100to200)", 100, 0, 10);
            else if(t.jet_pt[i] < 500.0) plot1DUnderOverFlow("h_RecoptOverGenpt_200to500", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt200to500)", 100, 0, 10);
            else plot1DUnderOverFlow("h_RecoptOverGenpt_pt500", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt500)", 100, 0, 10);

            plot1DUnderOverFlow("h_RebalancedptOverGenpt", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T}", 100, 0, 10);
            if(t.jet_pt[i] < 20.0) plot1DUnderOverFlow("h_RebalancedptOverGenpt_10to20", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt10to20)", 100, 0, 10);
            if(t.jet_pt[i] < 50.0) plot1DUnderOverFlow("h_RebalancedptOverGenpt_20to50", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt20to50)", 100, 0, 10);
            else if(t.jet_pt[i] < 100.0) plot1DUnderOverFlow("h_RebalancedptOverGenpt_50to100", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt50to100)", 100, 0, 10);
            else if(t.jet_pt[i] < 200.0) plot1DUnderOverFlow("h_RebalancedptOverGenpt_100to200", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt100to200)", 100, 0, 10);
            else if(t.jet_pt[i] < 500.0) plot1DUnderOverFlow("h_RebalancedptOverGenpt_200to500", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt200to500)", 100, 0, 10);
            else plot1DUnderOverFlow("h_RebalancedptOverGenpt_pt500", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (pt500)", 100, 0, 10);

            plot2DUnderOverFlow("h2d_RecoptOverGenpt_RebalancedptOverGenpt", t.jet_pt[i]/t.jet_mcPt[i], jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, "; reco jet p_{T}/gen jet p_{T}; reb jet p_{T}/gen jet p_{T}",100,0,10,100,0,10);
            plot2DUnderOverFlow("h2d_mcpt_rebalancedpt", t.jet_mcPt[i], jet_pt.at(jetCounter), evtweight_, h_1d_global, "; gen jet pt [GeV]; rebalanced jet pt [GeV]",100,0,1000,100,0,1000);
            plot2DUnderOverFlow("h2d_mcpt_recopt", t.jet_mcPt[i], t.jet_pt[i], evtweight_, h_1d_global, "; gen jet pt [GeV]; reco jet pt [GeV]",100,0,1000,100,0,1000);
            plot2DUnderOverFlow("h2d_rebalancedpt_recopt", jet_pt.at(jetCounter), t.jet_pt[i], evtweight_, h_1d_global, "; rebalanced jet pt [GeV]; reco jet pt [GeV]",100,0,1000,100,0,1000);
            plot2DUnderOverFlow("h2d_mcpt_RebalancedOverGen", t.jet_mcPt[i], jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, "; gen jet p_{T} [GeV]; reb jet p_{T}/gen jet p_{T}",100,0,1000,100,0,10);
            plot2DUnderOverFlow("h2d_recopt_RebalancedOverGen", t.jet_pt[i], jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, "; gen jet p_{T} [GeV]; reb jet p_{T}/gen jet p_{T}",100,0,1000,100,0,10);
            if(t.nJet40 >= 3){
              plot1DUnderOverFlow("h_RecoptOverGenpt_njet3", t.jet_pt[i]/t.jet_mcPt[i], evtweight_, h_1d_global, ";reco jet p_{T}/gen jet p_{T} (njet3)", 100, 0, 10);
              plot1DUnderOverFlow("h_RebalancedptOverGenpt_njet3", jet_pt.at(jetCounter)/t.jet_mcPt[i], evtweight_, h_1d_global, ";reb jet p_{T}/gen jet p_{T} (njet3)", 100, 0, 10);
            }
          }
        }
      }

      float new_met = sqrt(new_met_x*new_met_x + new_met_y*new_met_y);
      float pt_soft_true = sqrt(pt_soft_true_x*pt_soft_true_x + pt_soft_true_y*pt_soft_true_y);
      float pt_soft_reco = sqrt(pt_soft_reco_x*pt_soft_reco_x + pt_soft_reco_y*pt_soft_reco_y);

      plot1DUnderOverFlow("h_before_met_x", t.met_pt*cos(t.met_phi), evtweight_, h_1d_global, ";met_x", 50, 0, 500);
      plot1DUnderOverFlow("h_before_met_y", t.met_pt*sin(t.met_phi), evtweight_, h_1d_global, ";met_y", 50, 0, 500);
      plot1DUnderOverFlow("h_before_met", t.met_pt, evtweight_, h_1d_global, ";met", 250, 0, 250);
      plot1DUnderOverFlow("h_rebalanced_met_x", new_met_x, evtweight_, h_1d_global, ";rebalanced met_x", 250, 0, 250);
      plot1DUnderOverFlow("h_rebalanced_met_y", new_met_y, evtweight_, h_1d_global, ";rebalanced met_y", 250, 0, 250);
      plot1DUnderOverFlow("h_rebalanced_met", new_met, evtweight_, h_1d_global, ";rebalanced met", 250, 0, 250);
      plot2DUnderOverFlow("h2d_met_rebalancedmet", t.met_pt, new_met, evtweight_, h_1d_global, "; reco met [GeV]; rebalanced met [GeV]",50,0,500,20,0,200);
      plot1DUnderOverFlow("h_pt_soft_true", pt_soft_true, evtweight_, h_1d_global, ";p_{T, soft}^{true}", 20, 0, 200);
      plot1DUnderOverFlow("h_pt_soft_reco", pt_soft_reco, evtweight_, h_1d_global, ";p_{T, soft}^{reco}", 20, 0, 200);
      plot2DUnderOverFlow("h2d_ptsoft_reco_true", pt_soft_reco, pt_soft_true, evtweight_, h_1d_global, "; p_{T, soft}^{reco}; p_{T, soft}^{true}",200,0,200,200,0,200);
      plot1DUnderOverFlow("h_pt_soft_reco_minus_true", pt_soft_reco - pt_soft_true, evtweight_, h_1d_global, ";p_{T, soft}^{reco} - p_{T, soft}^{true}", 100, -50, 50);
  
      //if(new_met > 20.0) continue; //FIXME

      bool do_before_plots = true;
      if(t.nJet40 < 2) do_before_plots = false;
      //if(t.met_pt < 30.0) do_before_plots = false;
      if(t.ht >= 1000.0) do_before_plots = false;
      if(t.ht < 450.0) do_before_plots = false;
      //if(t.ht < 1000.0 && t.met_pt < 200.0) do_before_plots = false;
      //if(t.ht >= 1000.0 && t.met_pt < 30.0) do_before_plots = false;
      //if(t.diffMetMht/t.met_pt > 0.5) do_before_plots = false;
      //if(t.mt2 < 200.0) do_before_plots = false;

      // note: this will double count some leptons, since reco leptons can appear as PFcands
      nlepveto_ = t.nMuons10 + t.nElectrons10 + t.nPFLep5LowMT + t.nPFHad10LowMT;

      //if(t.deltaPhiMin < 0.3) do_before_plots = false;
      //if(nlepveto_ > 0) do_before_plots = false;


      if(do_before_plots){
        plot1DUnderOverFlow("h_mt2",       t.mt2,       evtweight_, h_1d_global, ";M_{T2} [GeV]", 40, 0, 1000);
        plot1DUnderOverFlow("h_met",       t.met_pt,    evtweight_, h_1d_global, ";MET [GeV]", 50, 0, 500);
        plot1DUnderOverFlow("h_mht",       t.mht_pt,    evtweight_, h_1d_global, ";MHT [GeV]", 50, 0, 500);
        plot1DUnderOverFlow("h_met_phi",   t.met_phi,   evtweight_, h_1d_global, ";MET Phi", 64, -3.2, 3.2);
        plot1DUnderOverFlow("h_jet1_pt",   t.jet1_pt,   evtweight_, h_1d_global, ";Jet1 p_{T} [GeV]", 100, 0, 1000);
        plot1DUnderOverFlow("h_jet2_pt",   t.jet2_pt,   evtweight_, h_1d_global, ";Jet2 p_{T} [GeV]", 100, 0, 1000);
        plot1DUnderOverFlow("h_ht",        t.ht,   evtweight_, h_1d_global, ";H_{T} [GeV]", 55, 450, 1000);
        //plot1DUnderOverFlow("h_jet1_pt",   t.jet1_pt,   evtweight_, h_1d_global, ";Jet1 p_{T} [GeV]", 60, 0, 3000);
        //plot1DUnderOverFlow("h_jet2_pt",   t.jet2_pt,   evtweight_, h_1d_global, ";Jet2 p_{T} [GeV]", 60, 0, 3000);
        //plot1DUnderOverFlow("h_ht",        t.ht,   evtweight_, h_1d_global, ";H_{T} [GeV]", 60, 1000, 4000);
        //plot1DUnderOverFlow("h_ht",        t.ht,   evtweight_, h_1d_global, ";H_{T} [GeV]", 80, 0, 4000);
        plot1DUnderOverFlow("h_nJet40",       t.nJet40,   evtweight_, h_1d_global, ";N(jets)", 15, 0, 15);
        plot1DUnderOverFlow("h_nBJet20",      t.nBJet20,   evtweight_, h_1d_global, ";N(bjets)", 6, 0, 6);
        plot1DUnderOverFlow("h_deltaPhiMin",  t.deltaPhiMin,   evtweight_, h_1d_global, ";#Delta#phi_{min}", 32, 0, 3.2);
        plot1DUnderOverFlow("h_diffMetMht",   t.diffMetMht,   evtweight_, h_1d_global, ";|E_{T}^{miss} - MHT| [GeV]", 120, 0, 300);
        plot1DUnderOverFlow("h_diffMetMhtOverMet",   t.diffMetMht/t.met_pt,   evtweight_, h_1d_global, ";|E_{T}^{miss} - MHT| / E_{T}^{miss}", 100, 0, 2.);
      }

      evtweight_ = t.evt_scale1fb * lumi * smearNormalization;//FIXME
      //evtweight_ = t.evt_scale1fb * lumi * smearNormalization * (1.0/met_eff);//FIXME
      random->SetSeed();

      float avg_met_pt = 0.0;
      float avg_mht_pt = 0.0;
      float avg_met_phi = 0.0;
      float avg_ht = 0.0;
      float avg_nJet40 = 0.0;
      float avg_nBJet20 = 0.0;

      //if (t.met_pt > 20) continue;//FIXME

      for(int iSmear=0; iSmear<numberOfSmears; iSmear++){

        std::vector<float> jet_pt_smeared = jet_pt;

        for(int i=0; i<jet_pt_smeared.size(); i++){
          double r = random->Rndm();
          std::vector<float> response = GetResponseVector(jet_pt_smeared.at(i), fabs(jet_eta.at(i)));
          float prob_integral = 0;
          int iBin = 0;
          while( prob_integral < r ){
            prob_integral += response.at(iBin);
            //plot2D("h2d_bin_response", iBin, response.at(iBin), evtweight_, h_1d_global, "; bin; response",1000,0,1000,100,0,0.2);
            iBin++;
            if(iBin==1000) break;
          }
          float smear = float(iBin) * 0.01;
          plot1D("h_smear",      smear,   evtweight_, h_1d_global, ";smear", 10000, 0, 10);
          plot1D("h_bin",        iBin,   evtweight_, h_1d_global, ";bin", 1000, 0, 1000);
          plot1D("h_random",     r,   evtweight_, h_1d_global, ";random number", 100, 0, 1);
          jet_pt_smeared.at(i) *= smear;
        }

        new_met_x = t.met_pt*cos(t.met_phi);
        new_met_y = t.met_pt*sin(t.met_phi);

        float jet_x = 0;
        float jet_y = 0;
        nJet40 = 0;
        nBJet20 = 0;
        ht = 0;
        for(int i=0; i<t.njet; i++){
          new_met_x += t.jet_pt[i]*cos(t.jet_phi[i]);
          new_met_y += t.jet_pt[i]*sin(t.jet_phi[i]);
        }
        for(int i=0; i<jet_pt_smeared.size(); i++){
          new_met_x -= jet_pt_smeared.at(i)*cos(jet_phi.at(i));
          new_met_y -= jet_pt_smeared.at(i)*sin(jet_phi.at(i));
          jet_x += (jet_pt_smeared.at(i))*cos(jet_phi.at(i));
          jet_y += (jet_pt_smeared.at(i))*sin(jet_phi.at(i));
          if( (jet_pt_smeared.at(i) > 40.0) && (fabs(jet_eta.at(i)) < 2.5) ) ht += jet_pt_smeared.at(i);
          if( (jet_pt_smeared.at(i) > 40.0) && (fabs(jet_eta.at(i)) < 2.5) ) nJet40++;
          if( (jet_pt_smeared.at(i) > 20.0) && (fabs(jet_eta.at(i)) < 2.5) && (jet_btagCSV.at(i) > 0.814) ) nBJet20++;
        }
        for(int i=0; i<PU_passes_id_jet_pt.size(); i++){
          new_met_x -= PU_passes_id_jet_pt.at(i)*cos(PU_passes_id_jet_phi.at(i));
          new_met_y -= PU_passes_id_jet_pt.at(i)*sin(PU_passes_id_jet_phi.at(i));
          jet_x += (PU_passes_id_jet_pt.at(i))*cos(PU_passes_id_jet_phi.at(i));
          jet_y += (PU_passes_id_jet_pt.at(i))*sin(PU_passes_id_jet_phi.at(i));
          if( (PU_passes_id_jet_pt.at(i) > 40.0) && (fabs(PU_passes_id_jet_eta.at(i)) < 2.5) ) ht += PU_passes_id_jet_pt.at(i);
          if( (PU_passes_id_jet_pt.at(i) > 40.0) && (fabs(PU_passes_id_jet_eta.at(i)) < 2.5) ) nJet40++;
          if( (PU_passes_id_jet_pt.at(i) > 20.0) && (fabs(PU_passes_id_jet_eta.at(i)) < 2.5) && (PU_passes_id_jet_btagCSV.at(i) > 0.814) ) nBJet20++;
        }
        for(int i=0; i<PU_fails_id_jet_pt.size(); i++){
          new_met_x -= PU_fails_id_jet_pt.at(i)*cos(PU_fails_id_jet_phi.at(i));
          new_met_y -= PU_fails_id_jet_pt.at(i)*sin(PU_fails_id_jet_phi.at(i));
          jet_x += (PU_fails_id_jet_pt.at(i))*cos(PU_fails_id_jet_phi.at(i));
          jet_y += (PU_fails_id_jet_pt.at(i))*sin(PU_fails_id_jet_phi.at(i));
        }


        met_pt = sqrt(new_met_x*new_met_x + new_met_y*new_met_y);
        TLorentzVector met_vec;
        met_vec.SetPxPyPzE(new_met_x, new_met_y, 0, met_pt);
        met_phi = met_vec.Phi();

        if(nJet40 < 2) continue;
        if(ht >= 1000.0) continue;
        if(ht < 450.0) continue;
        //if(met_pt < 30.0) continue;
        //if(ht < 1000.0 && met_pt < 200.0) continue;
        //if(ht >= 1000.0 && met_pt < 30.0) continue;

        std::vector<LorentzVector> p4sForDphi;
        for(int i=0; i<jet_pt_smeared.size(); i++){
          if(jet_pt_smeared.at(i) < 40.0 || fabs(jet_eta.at(i)) > 4.7) continue;
          LorentzVector tempVec;
          tempVec.SetPx((jet_pt_smeared.at(i))*cos(jet_phi.at(i)));
          tempVec.SetPy((jet_pt_smeared.at(i))*sin(jet_phi.at(i)));
          tempVec.SetPz((jet_pt_smeared.at(i))*sinh(jet_eta.at(i)));
          tempVec.SetE(tempVec.P());
          p4sForDphi.push_back(tempVec);
        }
        for(int i=0; i<PU_passes_id_jet_pt.size(); i++){
          if(PU_passes_id_jet_pt.at(i) < 40.0 || fabs(PU_passes_id_jet_eta.at(i)) > 4.7) continue;
          LorentzVector tempVec;
          tempVec.SetPx((PU_passes_id_jet_pt.at(i))*cos(PU_passes_id_jet_phi.at(i)));
          tempVec.SetPy((PU_passes_id_jet_pt.at(i))*sin(PU_passes_id_jet_phi.at(i)));
          tempVec.SetPz((PU_passes_id_jet_pt.at(i))*sinh(PU_passes_id_jet_eta.at(i)));
          tempVec.SetE(tempVec.P());
          p4sForDphi.push_back(tempVec);
        }
        sort(p4sForDphi.begin(), p4sForDphi.end(), sortByPt);

        deltaPhiMin = 999;
        for (unsigned int ip4 = 0; ip4 < p4sForDphi.size(); ++ip4) {
          if(ip4 < 4) deltaPhiMin = min(deltaPhiMin, DeltaPhi( met_phi, p4sForDphi.at(ip4).phi() ));
        }

        //if(deltaPhiMin < 0.3) continue;//FIXME
        //if(nlepveto_ > 0) continue;//FIXME

        std::vector<LorentzVector> p4sForHems;
        for(int i=0; i<jet_pt_smeared.size(); i++){
          if(jet_pt_smeared.at(i) < 40.0) continue;
          if(fabs(jet_eta.at(i)) > 2.5) continue;
          LorentzVector tempVec;
          tempVec.SetPx((jet_pt_smeared.at(i))*cos(jet_phi.at(i)));
          tempVec.SetPy((jet_pt_smeared.at(i))*sin(jet_phi.at(i)));
          tempVec.SetPz((jet_pt_smeared.at(i))*sinh(jet_eta.at(i)));
          tempVec.SetE(tempVec.P());
          p4sForHems.push_back(tempVec);
        }
        for(int i=0; i<PU_passes_id_jet_pt.size(); i++){
          if(PU_passes_id_jet_pt.at(i) < 40.0) continue;
          if(fabs(PU_passes_id_jet_eta.at(i)) > 2.5) continue;
          LorentzVector tempVec;
          tempVec.SetPx((PU_passes_id_jet_pt.at(i))*cos(PU_passes_id_jet_phi.at(i)));
          tempVec.SetPy((PU_passes_id_jet_pt.at(i))*sin(PU_passes_id_jet_phi.at(i)));
          tempVec.SetPz((PU_passes_id_jet_pt.at(i))*sinh(PU_passes_id_jet_eta.at(i)));
          tempVec.SetE(tempVec.P());
          p4sForHems.push_back(tempVec);
        }
        sort(p4sForHems.begin(), p4sForHems.end(), sortByPt);

        LorentzVector sumMhtp4 = LorentzVector(0,0,0,0);
        for (unsigned int ip4 = 0; ip4 < p4sForHems.size(); ++ip4) {
          sumMhtp4 -= p4sForHems.at(ip4);
        }
        float mht_pt  = sumMhtp4.pt();
        float mht_phi = sumMhtp4.phi();
        TVector2 mhtVector = TVector2(mht_pt*cos(mht_phi), mht_pt*sin(mht_phi));
        TVector2 metVector = TVector2(met_pt*cos(met_phi), met_pt*sin(met_phi));
        diffMetMht = (mhtVector - metVector).Mod();

        //if(diffMetMht/met_pt > 0.5) continue;

        vector<LorentzVector> hemJets;
        if(p4sForHems.size() > 1){
          hemJets = getHemJets(p4sForHems); 
        }else{
          std::cout << "Problem with jets for hemispheres!" << std::endl;
          return;
        }

        mt2 = HemMT2(met_pt, met_phi, hemJets.at(0), hemJets.at(1));

        //if(mt2 < 200.0) continue;

        jet1_pt = p4sForHems.at(0).pt();
        jet2_pt = p4sForHems.at(1).pt();

        plot1DUnderOverFlow("h_nvtx",      t.nVert,   evtweight_, h_1d_global, ";N(vtx)", 80, 0, 80);
        plot1DUnderOverFlow("h_smear_mt2",       mt2,       evtweight_, h_1d_global, ";M_{T2} [GeV]", 40, 0, 1000);
        plot1DUnderOverFlow("h_smear_met",       met_pt,    evtweight_, h_1d_global, ";MET [GeV]", 50, 0, 500);
        plot1DUnderOverFlow("h_smear_mht",       mht_pt,    evtweight_, h_1d_global, ";MHT [GeV]", 50, 0, 500);
        plot1DUnderOverFlow("h_smear_met_phi",   met_phi,   evtweight_, h_1d_global, ";MET Phi", 64, -3.2, 3.2);
        plot1DUnderOverFlow("h_smear_jet1_pt",   jet1_pt,   evtweight_, h_1d_global, ";Jet1 p_{T} [GeV]", 100, 0, 1000);
        plot1DUnderOverFlow("h_smear_jet2_pt",   jet2_pt,   evtweight_, h_1d_global, ";Jet2 p_{T} [GeV]", 100, 0, 1000);
        plot1DUnderOverFlow("h_smear_ht",        ht,   evtweight_, h_1d_global, ";H_{T} [GeV]", 55, 450, 1000);
        //plot1DUnderOverFlow("h_smear_jet1_pt",   jet1_pt,   evtweight_, h_1d_global, ";Jet1 p_{T} [GeV]", 60, 0, 3000);
        //plot1DUnderOverFlow("h_smear_jet2_pt",   jet2_pt,   evtweight_, h_1d_global, ";Jet2 p_{T} [GeV]", 60, 0, 3000);
        //plot1DUnderOverFlow("h_smear_ht",        ht,   evtweight_, h_1d_global, ";H_{T} [GeV]", 80, 0, 4000);
        //plot1DUnderOverFlow("h_smear_ht",        ht,   evtweight_, h_1d_global, ";H_{T} [GeV]", 60, 1000, 4000);
        plot1DUnderOverFlow("h_smear_nJet40",       nJet40,   evtweight_, h_1d_global, ";N(jets)", 15, 0, 15);
        plot1DUnderOverFlow("h_smear_nBJet20",      nBJet20,   evtweight_, h_1d_global, ";N(bjets)", 6, 0, 6);
        plot1DUnderOverFlow("h_smear_deltaPhiMin",  deltaPhiMin,   evtweight_, h_1d_global, ";#Delta#phi_{min}", 32, 0, 3.2);
        plot1DUnderOverFlow("h_smear_diffMetMht",   diffMetMht,   evtweight_, h_1d_global, ";|E_{T}^{miss} - MHT| [GeV]", 120, 0, 300);
        plot1DUnderOverFlow("h_smear_diffMetMhtOverMet",   diffMetMht/met_pt,   evtweight_, h_1d_global, ";|E_{T}^{miss} - MHT| / E_{T}^{miss}", 100, 0, 2.);


        ////////////////////////////////////
        /// done with overall selection  /// 
        ////////////////////////////////////
        ///   time to fill histograms    /// 
        ////////////////////////////////////

        fillHistos(SRNoCut.srHistMap, SRNoCut.GetName(), "");
        fillHistosSignalRegion("sr");
        fillHistosSRBase();
        fillHistosInclusive();

        //avg_met_pt += met_pt;
        //avg_mht_pt += mht_pt;
        //avg_met_phi += met_phi;
        //avg_ht += ht;
        //avg_nJet40 += nJet40;
        //avg_nBJet20 += nBJet20;

      }//end loop on smears

/*
      avg_met_pt *= (1.0/numberOfSmears);
      //avg_mht_pt *= (1.0/numberOfSmears);
      //avg_met_phi *= (1.0/numberOfSmears);
      avg_ht *= (1.0/numberOfSmears);
      avg_nJet40 *= (1.0/numberOfSmears);
      avg_nBJet20 *= (1.0/numberOfSmears);

      evtweight_ = t.evt_scale1fb * lumi;

      //plot1D("h_smear_mt2",       avg_mt2,       evtweight_, h_1d_global, ";M_{T2} [GeV]", 200, 0, 2000);
      plot1D("h_smear_met",       avg_met_pt,    evtweight_, h_1d_global, ";MET [GeV]", 50, 0, 500);
      plot1D("h_smear_mht",       avg_mht_pt,    evtweight_, h_1d_global, ";MHT [GeV]", 50, 0, 500);
      plot1D("h_smear_met_phi",   avg_met_phi,   evtweight_, h_1d_global, ";MET Phi", 64, -3.2, 3.2);
      //plot1D("h_smear_jet1_pt",   avg_jet1_pt,   evtweight_, h_1d_global, ";Jet1 p_{T} [GeV]", 100, 0, 1000);
      //plot1D("h_smear_jet2_pt",   avg_jet2_pt,   evtweight_, h_1d_global, ";Jet2 p_{T} [GeV]", 100, 0, 1000);
      plot1D("h_smear_ht",        avg_ht,   evtweight_, h_1d_global, ";H_{T} [GeV]", 120, 0, 3000);
      plot1D("h_smear_nJet40",       avg_nJet40,   evtweight_, h_1d_global, ";N(jets)", 15, 0, 15);
      plot1D("h_smear_nBJet20",      avg_nBJet20,   evtweight_, h_1d_global, ";N(bjets)", 6, 0, 6);
      //plot1D("h_smear_deltaPhiMin",  avg_deltaPhiMin,   evtweight_, h_1d_global, ";#Delta#phi_{min}", 32, 0, 3.2);
      //plot1D("h_smear_diffMetMht",   avg_diffMetMht,   evtweight_, h_1d_global, ";|E_{T}^{miss} - MHT| [GeV]", 120, 0, 300);
      //plot1D("h_smear_diffMetMhtOverMet",   avg_diffMetMht/avg_met_pt,   evtweight_, h_1d_global, ";|E_{T}^{miss} - MHT| / E_{T}^{miss}", 100, 0, 2.);


      ////////////////////////////////////
      /// done with overall selection  /// 
      ////////////////////////////////////
      ///   time to fill histograms    /// 
      ////////////////////////////////////

      //fillHistos(SRNoCut.srHistMap, SRNoCut.GetName(), "");
      //fillHistosSignalRegion("sr");
      //fillHistosSRBase();
      //fillHistosInclusive();
*/

    }//end loop on events in a file

    delete random;

    delete tree;
    f.Close();
  }//end loop on files

  cout << "[SmearLooper::loop] processed " << nEventsTotal << " events" << endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  cout << nDuplicates << " duplicate events were skipped." << endl;

  //---------------------
  // Save Plots
  //---------------------

  outfile_->cd();
  savePlotsDir(h_1d_global,outfile_,"");
  savePlotsDir(SRNoCut.srHistMap,outfile_,SRNoCut.GetName().c_str());
  savePlotsDir(SRBase.srHistMap,outfile_,SRBase.GetName().c_str());

  for(unsigned int srN = 0; srN < SRVec.size(); srN++){
    if(!SRVec.at(srN).srHistMap.empty()){
      savePlotsDir(SRVec.at(srN).srHistMap, outfile_, ("sr"+SRVec.at(srN).GetName()).c_str());
    }
  }

  //---------------------
  // Write and Close file
  //---------------------
  outfile_->Write();
  outfile_->Close();
  delete outfile_;

  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f s", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f s", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;

  return;
  }

  void SmearLooper::fillHistosSRBase() {

    std::map<std::string, float> values;
    values["deltaPhiMin"] = deltaPhiMin;
    values["diffMetMhtOverMet"]  = diffMetMht/met_pt;
    values["nlep"]        = nlepveto_;
    values["j1pt"]        = jet1_pt;
    values["j2pt"]        = jet2_pt;
    values["mt2"]         = mt2;
    values["passesHtMet"] = ( (ht > 450. && met_pt > 200.) || (ht > 1000. && met_pt > 30.) );

    if(SRBase.PassesSelection(values)) fillHistos(SRBase.srHistMap, SRBase.GetName(), "");

    return;
  }

  void SmearLooper::fillHistosInclusive() {

    std::map<std::string, float> values;
    values["deltaPhiMin"] = deltaPhiMin;
    values["diffMetMhtOverMet"]  = diffMetMht/met_pt;
    values["nlep"]        = nlepveto_;
    values["j1pt"]        = jet1_pt;
    values["j2pt"]        = jet2_pt;
    values["mt2"]         = mt2;
    values["passesHtMet"] = ( (ht > 450. && met_pt > 200.) || (ht > 1000. && met_pt > 30.) );

    for(unsigned int srN = 0; srN < InclusiveRegions.size(); srN++){
      std::map<std::string, float> values_temp = values;
      std::vector<std::string> vars = InclusiveRegions.at(srN).GetListOfVariables();
      for(unsigned int iVar=0; iVar<vars.size(); iVar++){
        if(vars.at(iVar) == "ht") values_temp["ht"] = t.ht;
        else if(vars.at(iVar) == "njets") values_temp["njets"] = nJet40;
        else if(vars.at(iVar) == "nbjets") values_temp["nbjets"] = nBJet20;
      }
      if(InclusiveRegions.at(srN).PassesSelection(values_temp)){
        fillHistos(InclusiveRegions.at(srN).srHistMap, InclusiveRegions.at(srN).GetName(), "");
      }
    }

    return;
  }

  void SmearLooper::fillHistosSignalRegion(const std::string& prefix, const std::string& suffix) {

    std::map<std::string, float> values;
    values["deltaPhiMin"] = deltaPhiMin;
    values["diffMetMhtOverMet"]  = diffMetMht/met_pt;
    values["nlep"]        = nlepveto_;
    values["j1pt"]        = jet1_pt;
    values["j2pt"]        = jet2_pt;
    values["njets"]       = nJet40;
    values["nbjets"]      = nBJet20;
    values["mt2"]         = mt2;
    values["ht"]          = ht;
    values["met"]         = met_pt;
    //values["passesHtMet"] = ( (t.ht > 450. && met_pt > 200.) || (t.ht > 1000. && met_pt > 30.) );


    for(unsigned int srN = 0; srN < SRVec.size(); srN++){
      if(SRVec.at(srN).PassesSelection(values)){
        fillHistos(SRVec.at(srN).srHistMap, prefix+SRVec.at(srN).GetName(), suffix);
        break;//signal regions are orthogonal, event cannot be in more than one
      }
    }

    //plot1D("h_SignalRegion",  sr_jets+sr_htmet,   evtweight_, h_1d_global, ";Signal Region", 100, 0, 100);

    return;
  }

  void SmearLooper::fillHistos(std::map<std::string, TH1*>& h_1d, const std::string& dirname, const std::string& s) {
    TDirectory * dir = (TDirectory*)outfile_->Get(dirname.c_str());
    if (dir == 0) {
      dir = outfile_->mkdir(dirname.c_str());
    } 
    dir->cd();

    plot1D("h_Events"+s,  1, 1, h_1d, ";Events, Unweighted", 1, 0, 2);
    plot1D("h_Events_w"+s,  1,   evtweight_, h_1d, ";Events, Weighted", 1, 0, 2);
    plot1D("h_mt2"+s,       mt2,   evtweight_, h_1d, "; M_{T2} [GeV]", 150, 0, 1500);
    plot1D("h_met"+s,       met_pt,   evtweight_, h_1d, ";E_{T}^{miss} [GeV]", 150, 0, 1500);
    plot1D("h_ht"+s,       ht,   evtweight_, h_1d, ";H_{T} [GeV]", 120, 0, 3000);
    plot1D("h_nJet40"+s,       nJet40,   evtweight_, h_1d, ";N(jets)", 15, 0, 15);
    plot1D("h_nBJet20"+s,      nBJet20,   evtweight_, h_1d, ";N(bjets)", 6, 0, 6);
    plot1D("h_deltaPhiMin"+s,  deltaPhiMin,   evtweight_, h_1d, ";#Delta#phi_{min}", 32, 0, 3.2);
    plot1D("h_diffMetMht"+s,   diffMetMht,   evtweight_, h_1d, ";|E_{T}^{miss} - MHT| [GeV]", 120, 0, 300);
    plot1D("h_diffMetMhtOverMet"+s,   diffMetMht/met_pt,   evtweight_, h_1d, ";|E_{T}^{miss} - MHT| / E_{T}^{miss}", 100, 0, 2.);
    plot1D("h_minMTBMet"+s,   t.minMTBMet,   evtweight_, h_1d, ";min M_{T}(b, E_{T}^{miss}) [GeV]", 150, 0, 1500);
    plot1D("h_nlepveto"+s,     nlepveto_,   evtweight_, h_1d, ";N(leps)", 10, 0, 10);
    plot1D("h_J0pt"+s,       jet1_pt,   evtweight_, h_1d, ";p_{T}(jet1) [GeV]", 150, 0, 1500);
    plot1D("h_J1pt"+s,       jet2_pt,   evtweight_, h_1d, ";p_{T}(jet2) [GeV]", 150, 0, 1500);
    plot1D("h_mt2bins"+s,       mt2,   evtweight_, h_1d, "; M_{T2} [GeV]", n_mt2bins, mt2bins);

    outfile_->cd();
    return;
  }
