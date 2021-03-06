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
#include "TMinuit.h"
#include "TList.h"
#include "TCollection.h"
#include "TKey.h"
#include "TF1.h"
#include "TStyle.h"
#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"

#include "../Tools/utils.h"

// header
#include "RebalanceLooper.h"

//MT2
#include "../MT2CORE/Plotting/PlotUtilities.h"


using namespace std;

class mt2tree;

//for rebalance and smear
std::vector<float> jetpt;
std::vector<float> jeteta;
std::vector<float> jetphi;
std::map<std::string, std::vector<double> > responseMap;
std::map<std::string, std::pair<TF1, float> > fitMap;
//float sigma = 10.0;
float sigma_soft;
float pt_soft_x;
float pt_soft_y;
//TH2F* NegLogL_par0_par1 = new TH2F("NegLogL_par0_par1", "NegLogL_par0_par1", 150, 0.495, 1.995, 150, 0.495, 1.995);
//TH2F* NegLogL_par0_par1 = new TH2F("NegLogL_par0_par1", "NegLogL_par0_par1", 20, 0.495, 0.695, 20, 0.495, 0.695);
//TH2F* NegLogL_par0_par1 = new TH2F("NegLogL_par0_par1", "NegLogL_par0_par1", 100, 0.695, 1.695, 100, 0.695, 1.695);
//TH2F* NegLogL_par0_par1 = new TH2F("NegLogL_par0_par1", "NegLogL_par0_par1", 200, 0.495, 2.495, 200, 0.495, 2.495);
TH2F* NegLogL_par0_par1 = new TH2F("NegLogL_par0_par1", "NegLogL_par0_par1", 170, 0.795, 2.495, 170, 0.795, 2.495);
TH2F* par0_par1 = new TH2F("par0_par1", "par0_par1", 100, 0.5, 0.6, 100, 0.5, 0.6);

void setup(){
  TFile* responseFile = new TFile("/home/users/jgran/QCD_13TeV_madgraph_PHYS14_fineBins_wideRange_eventVeto_newJetID.root", "READ");
  TIter it(responseFile->GetListOfKeys());
  TKey* key;
  std::string keep = "h_tot_JetAll_ResponsePt";
  while ( (key = (TKey *)it()) ) {
    if (strncmp (key->GetTitle(), keep.c_str(), keep.length()) != 0) continue;
    std::string hist_name = (key->GetTitle());
    TH1F* h = (TH1F*) responseFile->Get(TString(hist_name));
    h->Scale(1.0/h->Integral());
    std::vector<double> temp_vec;
    for(int i=1; i<=h->GetNbinsX(); i++){
      temp_vec.push_back(h->GetBinContent(i));
    } 
    responseMap[hist_name] = temp_vec;
  }
  responseFile->Close();
  delete responseFile;
}

void setup_fits(){
  TFile* fitFile = new TFile("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/template_fits.root", "READ");
  TIter it(fitFile->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey *)it()) ) {
    std::string fit_name = (key->GetTitle());
    //std::cout << fit_name << std::endl;
    TF1* f = (TF1*) fitFile->Get(TString(fit_name));
    float f_scale = 0.01/f->Integral(0, 10);//hist bin size = 0.01;
    std::cout << "f_scale = " << f_scale << std::endl;
    std::pair<TF1, float> my_pair;
    my_pair.second = f_scale;
    if(!f){
      std::cout << "Could not find fit " << fit_name << "!" << std::endl;
      continue;
    }
    my_pair.first = *f;
    //fitMap[fit_name] = *f;
    fitMap[fit_name] = my_pair;
  }
  fitFile->Close();
  delete fitFile;
}

float InitializeByMHT(std::vector<float>& pt_vec, std::vector<float>& phi_vec, int idx){
  float pt_x = 0.0;
  float pt_y = 0.0;
  pt_x -= pt_soft_x;
  pt_y -= pt_soft_y;
  for(int i=0; i<pt_vec.size(); i++){
    if(i==idx) continue;
    pt_x -= (pt_vec.at(i))*cos(phi_vec.at(i));
    pt_y -= (pt_vec.at(i))*sin(phi_vec.at(i));
  }
  float pt_value = sqrt(pt_x*pt_x + pt_y*pt_y);
  return pt_vec.at(idx)/pt_value;
}

std::vector<float> InitializeByFit(float pt, float eta){

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

  std::string fit_name = "fit_h_tot_JetAll_ResponsePt_" + pt_string + "_" + eta_string;

  float mean = fitMap.at(fit_name).first.GetParameter(4);
  float sigma = fitMap.at(fit_name).first.GetParameter(5);
 
  std::vector<float> result_vec;
  result_vec.push_back(mean);
  result_vec.push_back(sigma);
  return result_vec;
}

float GetValue(float xvalue, TH1F &hist){
  float xmax = hist.GetXaxis()->GetXmax()-0.01;
  int xbin   = hist.GetXaxis()->FindBin(min(xvalue, xmax));
  return hist.GetBinContent(xbin);
}

float GetValueUniformBinning(float xvalue, TH1F &hist){
  float bin_width = hist.GetBinWidth(1);
  //float bin_width = 0.01;
  int xbin = xvalue/bin_width;
  return hist.GetBinContent(xbin);
}

double ComputeGaussianValue(double scale, double mean, double sigma, float x){
  return scale*TMath::Exp(-0.5*((x - mean)/sigma)*((x - mean)/sigma));
}

float GetProb(float pt, float eta, double par){
  //std::cout << "par = " << par << std::endl;
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
  else pt_string = "Pt21";
*/
  //else if(pt < 6500.0) pt_string = "Pt21";

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
/*
  else if(eta < 3.2)  eta_string = "Eta8";
  else if(eta < 4.1)  eta_string = "Eta9";
  else if(eta < 5.0)  eta_string = "Eta10";
  else eta_string = "Eta11";
*/

  //std::cout << Form("h_tot_JetAll_ResponsePt_%s_%s", pt_string.c_str(), eta_string.c_str()) << std::endl;

/*
  if(par < 0) par = 0;
  if(par > 10) par = 10;
  double bin_width = 0.01;
  std::string hist_name = "h_tot_JetAll_ResponsePt_" + pt_string + "_" + eta_string;
  int bin = par/bin_width;
  bin = min(bin, 999);//only 1000 bins in the histogram
  return responseMap.at(hist_name).at(bin);
*/

  std::string fit_name = "fit_h_tot_JetAll_ResponsePt_" + pt_string + "_" + eta_string;
  //std::cout << "Looking up prob of " << par << " using fit " << fit_name << std::endl;
  double g_scale = fitMap.at(fit_name).first.GetParameter(0);
  double g_mean  = fitMap.at(fit_name).first.GetParameter(1);
  double g_sigma = fitMap.at(fit_name).first.GetParameter(2);
  
  double my_prob = ComputeGaussianValue(g_scale, g_mean, g_sigma, par);
  //double my_prob = fitMap.at(fit_name).first.Eval(par);
  double my_scale = fitMap.at(fit_name).second;
  return my_prob*my_scale;
}

void PlotLikelihood(std::vector<float> &mypar) {
  //std::cout << std::endl;
  par0_par1->Fill(mypar.at(0), mypar.at(1), 1.0);
  float likelihood = 0;
  float pt_constrained_x = 0.0;
  float pt_constrained_y = 0.0;
  for(int i=0; i<jetpt.size(); i++){
    std::cout << "par " << i << " = " << mypar.at(i) << std::endl;
    std::cout << "jet pt = " << jetpt.at(i) << std::endl;
    std::cout << "jet eta = " << jeteta.at(i) << std::endl;
    float prob = GetProb(jetpt.at(i)/mypar.at(i), fabs(jeteta.at(i)), mypar.at(i));
    float min_prob = 0.000001;
    prob = max(prob, min_prob);
    std::cout << "prob = " << prob << std::endl;
    likelihood += log(prob);
    pt_constrained_x -= (jetpt.at(i))*cos(jetphi.at(i))/mypar.at(i);
    pt_constrained_y -= (jetpt.at(i))*sin(jetphi.at(i))/mypar.at(i);
  }
  float x1 = (pt_soft_x - pt_constrained_x)/sigma_soft;
  float x2 = (pt_soft_y - pt_constrained_y)/sigma_soft;
  likelihood += -x1*x1/2;
  likelihood += -x2*x2/2;
  std::cout << "pt_constrained_x = " << pt_constrained_x << std::endl;
  std::cout << "pt_constrained_y = " << pt_constrained_y << std::endl;
  std::cout << "pt_soft_x = " << pt_soft_x << std::endl;
  std::cout << "pt_soft_y = " << pt_soft_y << std::endl;
  std::cout << "sigma_soft = " << sigma_soft << std::endl;
  std::cout << "x1 = " << x1 << std::endl;
  std::cout << "x2 = " << x2 << std::endl;
  std::cout << "-log likelihood = " << -likelihood << std::endl;
  std::cout << std::endl;

  float result = -likelihood;
  NegLogL_par0_par1->Fill(mypar.at(0), mypar.at(1), result);
  std::cout << "par.at(0), par.at(1), result = " << mypar.at(0) << ", " << mypar.at(1) << ", " << result << std::endl;
}

/*
//try requiring met = 0;
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  float likelihood = 0;
  float pt_constrained_x = 0.0;
  float pt_constrained_y = 0.0;
  float min_prob = 1E-20;;
  for(int i=0; i<(jetpt.size() -1); i++){
    float prob = GetProb(jetpt.at(i)/par[i], fabs(jeteta.at(i)), par[i]);
    prob = max(prob, min_prob);
    likelihood += log(prob);
    pt_constrained_x -= (jetpt.at(i))*cos(jetphi.at(i))/par[i];
    pt_constrained_y -= (jetpt.at(i))*sin(jetphi.at(i))/par[i];
  }
  pt_constrained_x -= pt_soft_x;
  pt_constrained_y -= pt_soft_y;
  float pt_constrained = sqrt(pt_constrained_x*pt_constrained_x + pt_constrained_y*pt_constrained_y);
  float prob = GetProb(pt_constrained, fabs(jeteta.at(jetpt.size() -1)), jetpt.at(jetpt.size() -1)/pt_constrained);
  prob = max(prob, min_prob);
  likelihood += log(prob);

  result = -likelihood;
}
*/

/*
//requiring mht = 0;
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  float likelihood = 0;
  float fit_mht_x = 0.0;
  float fit_mht_y = 0.0;
  float min_prob = 1E-20;;
  for(int i=0; i<jetpt.size(); i++){
    float prob = GetProb(jetpt.at(i)/par[i], fabs(jeteta.at(i)), par[i]);
    prob = max(prob, min_prob);
    likelihood += log(prob);
    fit_mht_x -= (jetpt.at(i))*cos(jetphi.at(i))/par[i];
    fit_mht_y -= (jetpt.at(i))*sin(jetphi.at(i))/par[i];
  }
  float fit_mht = sqrt(fit_mht_x*fit_mht_x + fit_mht_y*fit_mht_y);
  float fit_res = 5.0;
  float penalty = - (fit_mht*fit_mht)/(2*fit_res*fit_res);
  likelihood += penalty;

  result = -likelihood;
}
*/

//the standard one
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  float likelihood = 0;
  float pt_constrained_x = 0.0;
  float pt_constrained_y = 0.0;
  float min_prob = 1E-20;;
  for(unsigned int i=0; i<jetpt.size(); i++){
    //std::cout << "par " << i << " = " << par[i] << std::endl;
    //std::cout << "jet pt = " << jetpt.at(i) << std::endl;
    //std::cout << "jet eta = " << jeteta.at(i) << std::endl;
    float prob = GetProb((jetpt.at(i))/par[i], fabs(jeteta.at(i)), par[i]);
    //float prob = (1.92151e-02)*exp(-0.5*((par[i]-1.08963)/(1.95611e-0))*((par[i]-1.08963)/(1.95611e-0)));
    prob = max(prob, min_prob);
    //std::cout << "prob = " << prob << std::endl;
    likelihood += log(prob);
    pt_constrained_x -= (jetpt.at(i))*cos(jetphi.at(i))/par[i];
    pt_constrained_y -= (jetpt.at(i))*sin(jetphi.at(i))/par[i];
  }
  float x1 = (pt_soft_x - pt_constrained_x)/sigma_soft;
  float x2 = (pt_soft_y - pt_constrained_y)/sigma_soft;
  likelihood += -x1*x1/2;
  likelihood += -x2*x2/2;
  //std::cout << "pt_constrained_x = " << pt_constrained_x << std::endl;
  //std::cout << "pt_constrained_y = " << pt_constrained_y << std::endl;
  //std::cout << "pt_soft_x = " << pt_soft_x << std::endl;
  //std::cout << "pt_soft_y = " << pt_soft_y << std::endl;
  //std::cout << "sigma_soft = " << sigma_soft << std::endl;
  //std::cout << "x1 = " << x1 << std::endl;
  //std::cout << "x2 = " << x2 << std::endl;
  //std::cout << "-log likelihood = " << -likelihood << std::endl;
  //std::cout << std::endl;

  result = -likelihood;
}

/*
//this one uses chi2 with fixed sigmas
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  float chi2 = 0;
  float pt_constrained_x = 0.0;
  float pt_constrained_y = 0.0;
  for(int i=0; i<jetpt.size(); i++){
    float sigma = 0.05*jetpt.at(i);
    chi2 += ((jetpt.at(i) - ((jetpt.at(i))/par[i]))/sigma)*((jetpt.at(i) - ((jetpt.at(i))/par[i]))/sigma);
    pt_constrained_x -= (jetpt.at(i))*cos(jetphi.at(i))/par[i];
    pt_constrained_y -= (jetpt.at(i))*sin(jetphi.at(i))/par[i];
  }
  chi2 += ((pt_soft_x - pt_constrained_x)/sigma_soft)*((pt_soft_x - pt_constrained_x)/sigma_soft);
  chi2 += ((pt_soft_y - pt_constrained_y)/sigma_soft)*((pt_soft_y - pt_constrained_y)/sigma_soft);

  result = chi2;
}
*/

RebalanceLooper::RebalanceLooper(){
}
RebalanceLooper::~RebalanceLooper(){

};


void RebalanceLooper::loop(TChain* chain, std::string output_name){

  
  //setup();
  setup_fits();

/*
  TFile* sf_file = new TFile("/home/users/jgran/temp/V00-00-12/MT2Analysis/RebalanceLooper/rebalance_scale_factor.root", "READ");
  TH1F* h_sf = (TH1F*) sf_file->Get("sf");
  if(!h_sf){
    std::cout << "Couldn't open rebalance scale factor hist." << std::endl; 
    return;
  }
*/

/*
  float sum = 0.0;
  float min_toAdd = 1E-20;
  //float min_toAdd = 0.000001;
  for(int i=0; i<1000; i++){
    float toAdd = max(GetProb(100, 0, 0.01*i), min_toAdd);
    //std::cout << toAdd << std::endl;
    sum += toAdd;
  }
  std::cout << "sum = " << sum << std::endl;
  sum = 0.0;
  for(int i=0; i<1000; i++){
    float toAdd = max(GetProb(200, 1.0, 0.01*i), min_toAdd);
    sum += toAdd;
  }
  std::cout << "sum = " << sum << std::endl;
  sum = 0.0;
  for(int i=0; i<1000; i++){
    float toAdd = max(GetProb(500, 1.0, 0.01*i), min_toAdd);
    sum += toAdd;
  }
  std::cout << "sum = " << sum << std::endl;
  return;
*/

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  MakeNtuple( Form("%s", output_name.c_str()) );

  // File Loop
  int nDuplicates = 0;
  int nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  cout << "[RebalanceLooper::loop] running on " << nEventsChain << " events" << endl;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    cout << "[RebalanceLooper::loop] running on file: " << currentFile->GetTitle() << endl;

    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("mt2");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    //mt2tree t(tree);

    // Use this to speed things up when not looking at genParticles
    //tree->SetBranchStatus("genPart_*", 0); 
    //tree->SetBranchStatus("*", 0); 
    //tree->SetBranchStatus("njet", 1); 
    //tree->SetBranchStatus("jet_*", 1); 
    //tree->SetBranchStatus("met_*", 1); 
    //tree->SetBranchStatus("mht_*", 1); 

    t.Init(tree);

    // Event Loop
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    //for( unsigned int event = 0; event < (nEventsTree/10); ++event) {
    //for( unsigned int event = 0; event < 100; ++event) {
    //for( unsigned int event = 150000; event < 200000; ++event) {

      t.GetEntry(event);
      //if(t.evt != 1637495) continue;
      //if(t.evt != 477561) continue;

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

      InitNtuple();

/*
      for(int i=0; i<t.njet; i++){
        //if(t.jet_id[i] == 0) std::cout << "t.jet_id[" << i << "] = " << t.jet_id[i] << std::endl;
        std::cout << "t.jet_pt[" << i << "] = " << t.jet_pt[i] << std::endl;
        std::cout << "t.jet_mcPt[" << i << "] = " << t.jet_mcPt[i] << std::endl;
        std::cout << "t.jet_id[" << i << "] = " << t.jet_id[i] << std::endl;
        std::cout << "t.jet_puId[" << i << "] = " << t.jet_puId[i] << std::endl;
        std::cout << std::endl;
      }
*/

      //if(t.njet != 2) continue;//FIXME

      if(t.njet < 2){
        status = -1;
        FillNtuple(); 
        continue;
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


      ////////////////////////// rebalance //////////////////////////

      jetpt.clear();
      jeteta.clear();
      jetphi.clear();
      //responseHists.clear();

      float met_x = (t.met_pt)*cos(t.met_phi);
      float met_y = (t.met_pt)*sin(t.met_phi);

      float jet_x = 0;
      float jet_y = 0;
      for(int i=0; i<t.njet; i++){
/*
        if(t.jet_pt[i] < 20.0){
          useJet.push_back(0);
          continue;
        }
*/
        if(t.jet_pt[i] < 100.0 && t.jet_puId[i] < 1){
          useJet.push_back(0);
          met_x += (t.jet_pt[i])*cos(t.jet_phi[i]);//FIXME
          met_y += (t.jet_pt[i])*sin(t.jet_phi[i]);//FIXME
          continue;
        }
        //std::cout << "jet mcPt = " << t.jet_mcPt[i] << std::endl;
        jet_x += (t.jet_pt[i])*cos(t.jet_phi[i]);
        jet_y += (t.jet_pt[i])*sin(t.jet_phi[i]);
        useJet.push_back(1);
        //float sf = GetValue(t.jet_pt[i], *h_sf);
        //float corr_pt = (t.jet_pt[i])/sf;
        //jetpt.push_back(corr_pt);
        jetpt.push_back(t.jet_pt[i]);
        jeteta.push_back(t.jet_eta[i]);
        jetphi.push_back(t.jet_phi[i]);
      }

      if(jetpt.size() < 2){
        for(int i=0; i<jetpt.size(); i++){
          rebalanceFactors.push_back(1.0);
        }
        status = -1;
        FillNtuple(); 
        continue;
      }

      pt_soft_x = -met_x - jet_x;
      pt_soft_y = -met_y - jet_y;
      float pt_soft = sqrt(pt_soft_x*pt_soft_x + pt_soft_y*pt_soft_y);
      //sigma_soft = 0.1*pt_soft;
      //sigma_soft = 20.0;
      //sigma_soft = 30.0;
      //sigma_soft = 10.0;
      //sigma_soft = 15.0;
      sigma_soft = 25.0;
      //sigma_soft = 1.5*t.nVert;
      //sigma_soft = 5.0*sqrt(t.nVert);


/*
      //FIXME
      std::vector<float> parameters;
      for(int i=0; i<jetpt.size(); i++){
        parameters.push_back(-1);
      }
      for(float i=0; i<170; i++){
        for(float j=0; j<170; j++){
          parameters.clear();
          float temp1 = 0.8 + 0.01*i;
          float temp2 = 0.8 + 0.01*j;
          parameters.push_back(temp1);
          parameters.push_back(temp2);
          PlotLikelihood(parameters);
        }
      }

      gStyle->SetPaintTextFormat("1.1f");
      NegLogL_par0_par1->Draw("COLZ");
      //NegLogL_par0_par1->Draw("TEXTSAME");
      //NegLogL_par0_par1->Draw();
      //par0_par1->Draw("TEXT");
      //break;
      //FIXME
*/
      
      //if(t.met_pt > 50) continue;
      //if(pt_soft > 50) continue;
      
      TMinuit* minimizer = new TMinuit(jetpt.size());
      //TMinuit* minimizer = new TMinuit((jetpt.size()-1));
      minimizer->SetFCN(minuitFunction);
      int iflag = 0;
      Double_t arglist[10];

      arglist[0] = -1;
      minimizer->mnexcm("SET PRI", arglist, 1, iflag);

      arglist[0] = 1;
      minimizer->mnexcm("SET STRATEGY", arglist, 1, iflag);

      minimizer->SetErrorDef(0.5);
      //minimizer->SetErrorDef(10.0);

      for(int i=0; i<jetpt.size(); i++){
      //for(int i=0; i<(jetpt.size() - 1); i++){
        std::string name = Form("c%d", i);
        //minimizer->mnparm(i,name,1,0.01,0.05,10,iflag);
        //minimizer->mnparm(i,name,1.0,0.1,0.1,3.0,iflag);
        //minimizer->mnparm(i,name,1,0.1,0.1,5,iflag);
        //float initial_value = InitializeByMHT(jetpt, jetphi, i);
        //minimizer->mnparm(i,name,initial_value,0.05,0,0,iflag);
        //minimizer->mnparm(i,name,1.0,0.1,0.1,10,iflag);
        minimizer->mnparm(i,name,1.0,0.05,0,0,iflag);
        //std::vector<float> initial_values = InitializeByFit(jetpt.at(i), jeteta.at(i));
        //minimizer->mnparm(i,name,initial_values.at(0),5*initial_values.at(1),0,0,iflag);
      }
/*
        minimizer->mnparm(0,"c0",1.0489,0.01,0,0,iflag);
        minimizer->mnparm(1,"c1",1.0245,0.01,0,0,iflag);
        minimizer->mnparm(2,"c2",1.0428,0.01,0,0,iflag);
        minimizer->mnparm(3,"c3",1.0834,0.01,0,0,iflag);
        minimizer->mnparm(4,"c4",1.0903,0.01,0,0,iflag);
*/


      arglist[0] = 10000;
      arglist[1] = 1.0;
      //arglist[1] = 0;
      //arglist[2] = 1;

/*
      for(int i=0; i<jetpt.size(); i++){
        minimizer->FixParameter(i);
      }

      //for(int i=0; i<jetpt.size(); i++){
      for(int i=0; i<2; i++){
        minimizer->Release(i);
        minimizer->mnexcm("MIGRAD", arglist, 0, iflag);
        minimizer->FixParameter(i);
      }

      for(int i=0; i<jetpt.size(); i++){
        minimizer->Release(i);
      }
*/


      //arglist[0] = 100000;
      //arglist[1] = 20;

      //minimizer->mnexcm("SIMPLEX", arglist, 1, iflag);
      //std::cout << "SIMPLEX iflag = " << iflag << std::endl;
      minimizer->mnexcm("MIGRAD", arglist, 2, iflag);
      std::cout << "MIGRAD iflag = " << iflag << std::endl;
      //minimizer->mnexcm("MINOS", arglist, 3, iflag);
      //std::cout << "MINOS iflag = " << iflag << std::endl;
      //minimizer->mnexcm("MINIMIZE", arglist, 1, iflag);
      //std::cout << "MINIMIZE iflag = " << iflag << std::endl;
      //minimizer->mnexcm("SEEK", arglist, 2, iflag);
      //std::cout << "SEEK iflag = " << iflag << std::endl;
      status = iflag;
      
      if(iflag !=0){
        arglist[1] = 10.0;//easier threshold for convergence
        minimizer->mnexcm("MIGRAD", arglist, 2, iflag);
        std::cout << "second MIGRAD iflag = " << iflag << std::endl;
        status = iflag;
      }

/*
      arglist[0] = 5000;
      arglist[1] = 0;
      arglist[2] = 1;

      if(iflag !=0){
        minimizer->mnexcm("MINOS", arglist, 3, iflag);
        status = iflag;
      }
*/

/*
         std::cout << std::endl;
         std::cout << "iflag = " << iflag << std::endl;
         std::cout << "\t jet_x = " << jet_x << std::endl;
         std::cout << "\t jet_y = " << jet_y << std::endl;
         std::cout << "\t met_x = " << t.met_pt*cos(t.met_phi) << std::endl;
         std::cout << "\t met_y = " << t.met_pt*sin(t.met_phi) << std::endl;
         std::cout << "\t mht_x = " << t.mht_pt*cos(t.mht_phi) << std::endl;
         std::cout << "\t mht_y = " << t.mht_pt*sin(t.mht_phi) << std::endl;
         std::cout << "\t pt_soft_x = " << pt_soft_x << std::endl;
         std::cout << "\t pt_soft_y = " << pt_soft_y << std::endl;

         std::cout << "\t met_pt = " << t.met_pt << std::endl;
         std::cout << "\t mht_pt = " << t.mht_pt << std::endl;
         std::cout << "\t pt_soft = " << pt_soft << std::endl;
         std::cout << std::endl;
*/
         
      for(int i=0; i<jetpt.size(); i++){
        double par_value;
        double par_error;
        minimizer->GetParameter(i, par_value, par_error);
        par_value = 1.0/par_value;
        //std::cout << "\t Before correction, jet pt =  " << jetpt.at(i) << ", phi = " << jetphi.at(i) << ", eta = " << jeteta.at(i) << std::endl;
        //std::cout << "\t The value of the parameter is " << par_value << "+/- " << par_error << std::endl;
        jetpt.at(i) *= par_value;
        rebalanceFactors.push_back(par_value); 
      }

      delete minimizer;

      jet_x = 0;
      jet_y = 0;
      //std::cout << "\t Corrected:" << std::endl;
      for(int i=0; i<jetpt.size(); i++){
        jet_x += (jetpt.at(i))*cos(jetphi.at(i));
        jet_y += (jetpt.at(i))*sin(jetphi.at(i));
        //std::cout << "\t jet" << i << ":  reb pt = " << jetpt.at(i) << ", gen pt = " << t.jet_mcPt[i] << " | reb/gen = " << jetpt.at(i)/t.jet_mcPt[i] << std::endl;
      }
      //std::cout << "\t new jet_x = " << jet_x << std::endl;
      //std::cout << "\t new jet_y = " << jet_y << std::endl;
      float new_met_x = -pt_soft_x - jet_x;
      float new_met_y = -pt_soft_y - jet_y;
      float new_met = sqrt(new_met_x*new_met_x + new_met_y*new_met_y);
      float new_mht = sqrt(jet_x*jet_x + jet_y*jet_y);

      //std::cout << "old met = " << t.met_pt << ", new met = " << new_met << std::endl;
      //std::cout << "old mht = " << t.mht_pt << ", new mht = " << new_mht << std::endl;
      //std::cout << "pt_soft = " << pt_soft << std::endl;
      //std::cout << std::endl;


/*
      plot1D("h_metx",       met_x,       evtweight_, h_1d_global, "; metx)", 40, -100, 100);
      plot1D("h_mety",       met_y,       evtweight_, h_1d_global, "; mety)", 40, -100, 100);
      plot1D("h_met" ,       t.met_pt,   evtweight_, h_1d_global, "; mety)", 40, 0, 100);
      plot1D("h_ptsoftx_plus_jetx",       pt_soft_x + jet_x,       evtweight_, h_1d_global, "; balance in x)", 40, -100, 100);
      plot1D("h_ptsofty_plus_jety",       pt_soft_y + jet_y,       evtweight_, h_1d_global, "; balance in y)", 40, -100, 100);
      plot1D("h_rebalanced_metx",         new_met_x,      evtweight_, h_1d_global, "; rebalanced metx)", 40, -100, 100);
      plot1D("h_rebalanced_mety",         new_met_y,      evtweight_, h_1d_global, "; rebalanced mety)", 40, -100, 100);
      plot1D("h_rebalanced_met" ,         new_met,        evtweight_, h_1d_global, "; rebalanced met)", 40, 0, 100);
*/
      /////////////////////        end          /////////////////////
     
      FillNtuple();

    }//end loop on events in a file

    delete tree;
    f.Close();
    }//end loop on files

    cout << "[RebalanceLooper::loop] processed " << nEventsTotal << " events" << endl;
    if ( nEventsChain != nEventsTotal ) {
      std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    cout << nDuplicates << " duplicate events were skipped." << endl;

    CloseNtuple();

    bmark->Stop("benchmark");
    cout << endl;
    cout << nEventsTotal << " Events Processed" << endl;
    cout << "------------------------------" << endl;
    cout << "CPU  Time:	" << Form( "%.01f s", bmark->GetCpuTime("benchmark")  ) << endl;
    cout << "Real Time:	" << Form( "%.01f s", bmark->GetRealTime("benchmark") ) << endl;
    cout << endl;
    delete bmark;

    //delete h_sf;
    //sf_file->Close();
    //delete sf_file;

    return;
}


void RebalanceLooper::MakeNtuple(const char *Filename){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  outfile_ = new TFile(Form("%s", Filename), "RECREATE");
  outfile_->cd();
  outTree_ = new TTree("rebalance", "rebalance");

  outTree_->Branch("rebalanceFactors", &rebalanceFactors );
  outTree_->Branch("useJet", &useJet );
  outTree_->Branch("status", &status );

  return;
}


void RebalanceLooper::InitNtuple () {
  rebalanceFactors.clear();
  useJet.clear();
  status = -999;
  return;
}


void RebalanceLooper::FillNtuple(){
  outTree_->Fill();
  return;
}

void RebalanceLooper::CloseNtuple(){
  outfile_->cd();
  outTree_->Write();
  outfile_->Close();
  return;
}
