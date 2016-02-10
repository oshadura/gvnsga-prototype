#ifndef __HISTOGRAMMANAGER__
#define __HISTOGRAMMANAGER__

#include "Population.h"
#include "TObject.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <vector>

class HistogramManager : public TObject {
private:
  //========= Individual ========//
  TH1F *fPopulationSize; //#0
  TH1F *hAllev;          //#1
  TH1F *hBuffev;         //#3
  TH1F *hThread;         //#4
  TH1F *hPriority;       //#4
  TH1F *hSteps;          //#5
  TH1F *hVector;         //#6
  TH1F *hMaxVector;      //#7
  //========= Fitness ========//
  TH1F *hMemory;         //#0
  TH1F *hTime;           //#1

public:
  HistogramManager(TDirectory *dir);
  virtual ~HistogramManager() {}
  Bool_t HistoFill(Population<Double_t> &pop, char *file);
  Bool_t CheckValue(ROOT::Internal::TTreeReaderValueBase *value);

  ClassDef(HistogramManager, 1)
};

#endif
