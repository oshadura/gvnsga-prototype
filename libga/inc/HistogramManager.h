#ifndef __HISTOGRAMMANAGER__
#define __HISTOGRAMMANAGER__

#include "Population.h"
#include "TObject.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <vector>

class HistogramManager : public TObject {
private:
  TH1F *fPopulationSize;
  TH1F *fAllev;
  TH1F *fBuffev;
  TH1F *fThread;
  TH1F *fPriority;
  TH1F *fSteps;
  TH1F *fVector;
  TH1F *fMemory;
  TH1F *fTime;

public:
  HistogramManager(TDirectory *dir);
  virtual ~HistogramManager() {}
  void HFill(Population<Double_t> *pop, char *file);
  bool CheckValue(ROOT::TTreeReaderValueBase *value);

  ClassDef(HistogramManager, 1)
};

#endif
