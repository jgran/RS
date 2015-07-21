#ifndef DO_SCALE_H
#define DO_SCALE_H

#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TKey.h>
#include <TH1F.h>

class scale{

  public:

  scale() {};
  ~scale() {};
  
  void do_scale(std::string, float);

};
#endif
