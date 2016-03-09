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
  GeantVFitness() : hfile(0), hMemVirt(0), hMemRes(0), fMemoryVector(0), maxMemResident(0), fMemorySwitch(true) {}

  virtual ~GeantVFitness() {
    hfile->cd();
    hfile->Write();
    hfile->Close();
    delete hfile;
    maxMemResident = 0;
  }
  
  Double_t GetmaxMemResident() const { return maxMemResident; }
  Bool_t SetMemorySwitch() const { return fMemorySwitch; }
  void SetMemorySwitch(Bool_t i){fMemorySwitch = i;}
  void LogMemoryFitness(std::string file);
  void LogTimeFitness();
  void TemporarySolution();

private:
  TFile *hfile;
  TH1F *hMemVirt;
  TH1F *hMemRes;
  std::vector<ProcInfo_t> fMemoryVector;
  Double_t maxMemResident;
  Bool_t fMemorySwitch;

  ClassDef(GeantVFitness, 1)
};

#endif
