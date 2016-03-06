#ifndef __GEANTVFITNESS__
#define __GEANTVFITNESS__

#include <typeinfo>
#include <vector>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

class GeantVFitness {
public:
  GeantVFitness() : hfile(0), hMemVirt(0), hMemRes(0), fMemoryVector(0) {}

  virtual ~GeantVFitness() {
    hfile->cd();
    hfile->Write();
    hfile->Close();
    delete hfile;
  }
  Double_t GetmaxMemResident() const { return maxMemResident; }
  void LogMemoryFitness();
  void LogTimeFitness();
  void HistOutputFitness(std::string file);

private:
  TFile *hfile;
  TH1F *hMemVirt;
  TH1F *hMemRes;
  std::vector<ProcInfo_t> fMemoryVector;
  Double_t maxMemResident;

  ClassDef(GeantVFitness, 1)
};

#endif
