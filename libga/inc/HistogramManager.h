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
  TH1F *hAllev;
  TH1F *hBuffev;
  TH1F *hThread;
  TH1F *hPriority;
  TH1F *hSteps;
  TH1F *hVector;
  TH1F *hMemory;
  TH1F *hTime;

public:
  HistogramManager(TDirectory *dir);
  virtual ~HistogramManager() {}
  void HFill(Population<Double_t> &pop, char *file);
  bool CheckValue(ROOT::Internal::TTreeReaderValueBase *value);

  ClassDef(HistogramManager, 1)
};

#endif
