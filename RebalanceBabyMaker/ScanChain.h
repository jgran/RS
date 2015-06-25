#ifndef ScanChain_h
#define ScanChain_h

// C++ includes
//#include <string>
//#include <vector>
#include <stdint.h>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class babyMaker {

 public:

  babyMaker() {};
  ~babyMaker() {
    delete BabyFile_;
    delete BabyTree_;
  };

  void ScanChain(TChain*, std::string = "testSample");

  void MakeBabyNtuple(const char *);
  void InitBabyNtuple();
  void FillBabyNtuple();
  void CloseBabyNtuple();

 private:

  TFile *BabyFile_;
  TTree *BabyTree_;

  TH1D* count_hist_;

  //baby ntuple variables

  Int_t           nPUvertices;

  Int_t           run;
  Int_t           lumi;
  ULong64_t       evt;
  Int_t           isData;

  Float_t         evt_scale1fb;
  Float_t         evt_xsec;
  Float_t         evt_kfactor;
  Float_t         evt_filter;
  ULong64_t       evt_nEvts;
  Int_t           evt_id;
  Float_t         puWeight;
  Int_t           nVert;
  Int_t           nTrueInt;
  Float_t         rho;

  Int_t           nJet40;
  Int_t           nBJet20;
  Int_t           nBJet25;
  Int_t           nBJet40;
  Int_t           nMuons10;
  Int_t           nElectrons10;
  Int_t           nLepLowMT;
  Int_t           nTaus20;
  Int_t           nGammas20;

  Float_t         deltaPhiMin;
  Float_t         diffMetMht;
  Float_t         minMTBMet;
  Float_t         gamma_minMTBMet;
  Float_t         zll_minMTBMet;

  Float_t         ht;
  Float_t         mt2;
  Float_t         mt2_gen;

  Float_t         jet1_pt;
  Float_t         jet2_pt;
  Float_t         gamma_jet1_pt;
  Float_t         gamma_jet2_pt;

  Float_t         pseudoJet1_pt;
  Float_t         pseudoJet1_eta;
  Float_t         pseudoJet1_phi;
  Float_t         pseudoJet1_mass;
  Float_t         pseudoJet2_pt;
  Float_t         pseudoJet2_eta;
  Float_t         pseudoJet2_phi;
  Float_t         pseudoJet2_mass;

  Float_t         mht_pt;
  Float_t         mht_phi;
  Float_t         met_pt;
  Float_t         met_phi;
  Float_t         met_rawPt;
  Float_t         met_rawPhi;
  Float_t         met_caloPt;
  Float_t         met_caloPhi;
  Float_t         met_genPt;
  Float_t         met_genPhi;

//----- MET FILTERS
  Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
  Int_t           Flag_trkPOG_manystripclus53X;
  Int_t           Flag_ecalLaserCorrFilter;
  Int_t           Flag_trkPOG_toomanystripclus53X;
  Int_t           Flag_hcalLaserEventFilter;
  Int_t           Flag_trkPOG_logErrorTooManyClusters;
  Int_t           Flag_trkPOGFilters;
  Int_t           Flag_trackingFailureFilter;
  Int_t           Flag_CSCTightHaloFilter;
  Int_t           Flag_HBHENoiseFilter;
  Int_t           Flag_goodVertices;
  Int_t           Flag_eeBadScFilter;
  Int_t           Flag_METFilters;

//----- TRIGGER (to be better defined)
  Int_t           HLT_HT900;   
  Int_t           HLT_MET170;
  Int_t           HLT_ht350met120;   
  Int_t           HLT_SingleMu;   
  Int_t           HLT_DoubleEl;   
  Int_t           HLT_MuEG;   
  Int_t           HLT_DoubleMu;   
  Int_t           HLT_Photons;   

//----- LEPTONS
  static const int max_nlep = 50;
  Int_t           nlep;
  Float_t         lep_pt[max_nlep];   //[nlep]
  Float_t         lep_eta[max_nlep];   //[nlep]
  Float_t         lep_phi[max_nlep];   //[nlep]
  Float_t         lep_mass[max_nlep];   //[nlep]
  Int_t           lep_charge[max_nlep];   //[nlep]
  Int_t           lep_pdgId[max_nlep];   //[nlep]
  Float_t         lep_dxy[max_nlep];   //[nlep]
  Float_t         lep_dz[max_nlep];   //[nlep]
  Int_t           lep_tightId[max_nlep];   //[nlep]
  Float_t         lep_relIso03[max_nlep];   //[nlep]
  Float_t         lep_relIso04[max_nlep];   //[nlep]
  Float_t         lep_miniRelIso[max_nlep];   //[nlep]
  Int_t           lep_mcMatchId[max_nlep];   //[nlep]
  Int_t           lep_lostHits[max_nlep];   //[nlep]
  Int_t           lep_convVeto[max_nlep];   //[nlep]
  Int_t           lep_tightCharge[max_nlep];   //[nlep]

//----- ISOLATED TRACK
  static const int max_nisoTrack = 50;
  Int_t             nisoTrack;
  Float_t           isoTrack_pt[max_nisoTrack];
  Float_t           isoTrack_eta[max_nisoTrack];
  Float_t           isoTrack_phi[max_nisoTrack];
  Float_t           isoTrack_mass[max_nisoTrack];
  Float_t           isoTrack_absIso[max_nisoTrack];
  Float_t           isoTrack_dz[max_nisoTrack];
  Int_t             isoTrack_pdgId[max_nisoTrack];
  Int_t             isoTrack_mcMatchId[max_nisoTrack];

  Int_t             nPFLep5LowMT;
  Int_t             nPFHad10LowMT;

//----- TAUS
  static const int max_ntau = 50;
  Int_t           ntau;
  Float_t         tau_pt[max_ntau];   //[ntau]
  Float_t         tau_eta[max_ntau];   //[ntau]
  Float_t         tau_phi[max_ntau];   //[ntau]
  Float_t         tau_mass[max_ntau];   //[ntau]
  Int_t           tau_charge[max_ntau];   //[ntau]
  Float_t         tau_dxy[max_ntau];   //[ntau]
  Float_t         tau_dz[max_ntau];   //[ntau]
  Int_t           tau_idCI3hit[max_ntau];   //[ntau]
  Float_t         tau_isoCI3hit[max_ntau];   //[ntau]
  Int_t           tau_mcMatchId[max_ntau];   //[ntau]

//----- PHOTONS
  static const int max_ngamma = 50;
  Int_t           ngamma;
  Float_t         gamma_pt[max_ngamma];   //[ngamma]
  Float_t         gamma_eta[max_ngamma];   //[ngamma]
  Float_t         gamma_phi[max_ngamma];   //[ngamma]
  Float_t         gamma_mass[max_ngamma];   //[ngamma]
  Int_t           gamma_mcMatchId[max_ngamma];   //[ngamma]
  Float_t         gamma_genIso04[max_ngamma];   //[ngamma]
  Float_t         gamma_drMinParton[max_ngamma];   //[ngamma]
  Float_t         gamma_chHadIso[max_ngamma];   //[ngamma]
  Float_t         gamma_neuHadIso[max_ngamma];   //[ngamma]
  Float_t         gamma_phIso[max_ngamma];   //[ngamma]
  Float_t         gamma_sigmaIetaIeta[max_ngamma];   //[ngamma]
  Float_t         gamma_r9[max_ngamma];   //[ngamma]
  Float_t         gamma_hOverE[max_ngamma];   //[ngamma]
  Int_t           gamma_idCutBased[max_ngamma];   //[ngamma]

  // event level vars recalculated for photon+jets control region
  Float_t         gamma_mt2;
  Int_t           gamma_nJet40;
  Int_t           gamma_nBJet20;
  Int_t           gamma_nBJet25;
  Int_t           gamma_nBJet40;
  Float_t         gamma_ht;
  Float_t         gamma_deltaPhiMin;
  Float_t         gamma_diffMetMht;
  Float_t         gamma_mht_pt;
  Float_t         gamma_mht_phi;
  Float_t         gamma_met_pt;
  Float_t         gamma_met_phi;

  // event level vars recalculated for Z-->ll control region
  Float_t         zll_mt2;
  Float_t         zll_deltaPhiMin;
  Float_t         zll_diffMetMht;
  Float_t         zll_met_pt;
  Float_t         zll_met_phi;
  Float_t         zll_mht_pt;
  Float_t         zll_mht_phi;
  Float_t         zll_mass;
  Float_t         zll_pt;
  Float_t         zll_eta;
  Float_t         zll_phi;
  Float_t         zll_ht;

//----- GEN PARTICLES
  static const int max_ngenPart = 300;
  Int_t           ngenPart;
  Float_t         genPart_pt[max_ngenPart];   //[ngenPart]
  Float_t         genPart_eta[max_ngenPart];   //[ngenPart]
  Float_t         genPart_phi[max_ngenPart];   //[ngenPart]
  Float_t         genPart_mass[max_ngenPart];   //[ngenPart]
  Int_t           genPart_pdgId[max_ngenPart];   //[ngenPart]
  Int_t           genPart_status[max_ngenPart];   //[ngenPart]
  Float_t         genPart_charge[max_ngenPart];   //[ngenPart]
  Int_t           genPart_motherId[max_ngenPart];   //[ngenPart]
  Int_t           genPart_grandmotherId[max_ngenPart];   //[ngenPart]

//----- GEN LEPTONS (ELECTRONS/MUONS)
  static const int max_ngenLep = 10;
  Int_t           ngenLep;
  Float_t         genLep_pt[max_ngenLep];   //[ngenLep]
  Float_t         genLep_eta[max_ngenLep];   //[ngenLep]
  Float_t         genLep_phi[max_ngenLep];   //[ngenLep]
  Float_t         genLep_mass[max_ngenLep];   //[ngenLep]
  Int_t           genLep_pdgId[max_ngenLep];   //[ngenLep]
  Int_t           genLep_status[max_ngenLep];   //[ngenLep]
  Float_t         genLep_charge[max_ngenLep];   //[ngenLep]
  Int_t           genLep_sourceId[max_ngenLep];   //[ngenLep]

//----- GEN TAUS
  static const int max_ngenTau = 10;
  Int_t           ngenTau;
  Float_t         genTau_pt[max_ngenTau];   //[ngenTau]
  Float_t         genTau_eta[max_ngenTau];   //[ngenTau]
  Float_t         genTau_phi[max_ngenTau];   //[ngenTau]
  Float_t         genTau_mass[max_ngenTau];   //[ngenTau]
  Int_t           genTau_pdgId[max_ngenTau];   //[ngenTau]
  Int_t           genTau_status[max_ngenTau];   //[ngenTau]
  Float_t         genTau_charge[max_ngenTau];   //[ngenTau]
  Int_t           genTau_sourceId[max_ngenTau];   //[ngenTau]

//----- GEN LEPTONS FROM TAUS
  static const int max_ngenLepFromTau = 10;
  Int_t           ngenLepFromTau;
  Float_t         genLepFromTau_pt[max_ngenLepFromTau];   //[ngenLepFromTau]
  Float_t         genLepFromTau_eta[max_ngenLepFromTau];   //[ngenLepFromTau]
  Float_t         genLepFromTau_phi[max_ngenLepFromTau];   //[ngenLepFromTau]
  Float_t         genLepFromTau_mass[max_ngenLepFromTau];   //[ngenLepFromTau]
  Int_t           genLepFromTau_pdgId[max_ngenLepFromTau];   //[ngenLepFromTau]
  Int_t           genLepFromTau_status[max_ngenLepFromTau];   //[ngenLepFromTau]
  Float_t         genLepFromTau_charge[max_ngenLepFromTau];   //[ngenLepFromTau]
  Int_t           genLepFromTau_sourceId[max_ngenLepFromTau];   //[ngenLepFromTau]

//----- JETS
  static const int max_njet = 100;
  Int_t           njet;
  Float_t         jet_pt[max_njet];   //[njet]
  Float_t         jet_eta[max_njet];   //[njet]
  Float_t         jet_phi[max_njet];   //[njet]
  Float_t         jet_mass[max_njet];   //[njet]
  Float_t         jet_btagCSV[max_njet];   //[njet]
  Float_t         jet_rawPt[max_njet];   //[njet]
  Float_t         jet_mcPt[max_njet];   //[njet]
  Int_t           jet_mcFlavour[max_njet];   //[njet]
  Float_t         jet_qgl[max_njet];   //[njet]
  Float_t         jet_area[max_njet];   //[njet]
  Int_t           jet_id[max_njet];   //[njet]
  Int_t           jet_puId[max_njet];   //[njet]
  std::vector<std::vector<Int_t> > jet_pfcandIndices;

//----- PFCANDS
  static const int max_npfcand = 10000;
  Int_t           npfcand;
  Float_t         pfcand_pt[max_npfcand];
  Float_t         pfcand_eta[max_npfcand];
  Float_t         pfcand_phi[max_npfcand];
  Float_t         pfcand_dz[max_npfcand];
  Int_t           pfcand_charge[max_npfcand];
  Int_t           pfcand_fromPV[max_npfcand];

//FIXME
  static const int max_ngenjet = 100;
  Int_t           ngenjet;
  Float_t         genjet_pt[max_ngenjet];   //[ngenjet]
  Float_t         genjet_eta[max_ngenjet];   //[ngenjet]
  Float_t         genjet_phi[max_ngenjet];   //[ngenjet]

//----- SUSY SIGNALS
  Int_t           GenSusyMScan1;
  Int_t           GenSusyMScan2;
  Int_t           GenSusyMScan3;
  Int_t           GenSusyMScan4;

//----- WEIGHTS AND VARIATIONS
  Float_t         weight_lepsf;
  Float_t         weight_lepsf_UP;
  Float_t         weight_lepsf_DN;
  Float_t         weight_btagsf;
  Float_t         weight_btagsf_UP;
  Float_t         weight_btagsf_DN;
  Float_t         weight_sigtrigsf;
  Float_t         weight_dileptrigsf;
  Float_t         weight_phottrigsf;
  Float_t         weight_pu;
  Float_t         weight_isr;
  Float_t         weight_scales_UP;
  Float_t         weight_scales_DN;
  Float_t         weight_pdfs_UP;
  Float_t         weight_pdfs_DN;

};

#endif

