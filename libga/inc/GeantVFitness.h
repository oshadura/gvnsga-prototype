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
  GeantVFitness(): hfile(0), hMemVirt(0), hMemRes(0), fMemoryVector(0) {
  	hMemRes = new TH1F("memory_resident", "Resident memory usage", 0, 0, 0);
  	hMemVirt = new TH1F("memory_virtual", "Virtual memory usage", 0, 0, 0);
  }

  virtual ~GeantVFitness() {
  	// think think!!!
  }

  void LogMemoryFitness();
  void LogTimeFitness();
  void HistOutputFitness(std::string file);

private:
  TFile *hfile;
  TH1F *hMemVirt;
  TH1F *hMemRes;
  std::vector<ProcInfo_t> fMemoryVector;

  ClassDef(GeantVFitness, 1)
};

#endif
