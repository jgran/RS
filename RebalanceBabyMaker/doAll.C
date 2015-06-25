{

  gSystem->Load("libMiniFWLite.so");
  gSystem->Load("libBabymakerCORE.so");
  gSystem->Load("libBabymakerTools.so");
  gSystem->Load("libBabymakerMT2CORE.so");
  gSystem->Load("libScanChain.so");

  TChain *ch = new TChain("Events"); 
 
  //CMS3 ntuples made from  
  // /hadoop/cms/store/user/dalfonso/MINIAOD/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25.1.root
  // /hadoop/cms/store/user/dalfonso/MINIAOD/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25.2.root 
  //ch->Add("/nfs-3/userdata/jgran/MT2Sync/ntuple1.root");
  //ch->Add("/nfs-3/userdata/jgran/MT2Sync/ntuple2.root");
  //  ch->Add("/home/users/gzevi/miniAOD/CMSSW_7_0_6/src/CMS2/NtupleMaker/test/MT2ntuple1000.root");
  //  ch->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_190.root");
  //ch->Add("/home/users/jgran/CMSSW_7_0_6_patch1/src/CMS2/NtupleMaker/ntuple_first1000.root");
  //ch->Add("/nfs-7/userdata/olivito/mt2/phys14/ttjets_cms3_1000.root");
  //ch->Add("/nfs-7/userdata/olivito/mt2/phys14/gjets_ht200to400_cms3_1000.root");
  //ch->Add("/hadoop/cms/store/group/snt/phys14/QCD_Pt-800to1000_Tune4C_13TeV_pythia8_Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/V07-02-08/merged_ntuple_1.root");
  ch->Add("/hadoop/cms/store/group/snt/phys14/QCD_Pt-170to300_Tune4C_13TeV_pythia8_Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/V07-02-08/merged_ntuple_1.root");
  //ch->Add("/hadoop/cms/store/group/snt/phys14/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8_Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/V07-02-08/merged_ntuple_1.root");

  babyMaker *looper = new babyMaker();
  looper->ScanChain(ch, "sntMT2Baby"); 
}
