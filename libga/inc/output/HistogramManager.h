#ifndef __HISTOGRAMMANAGER__
#define __HISTOGRAMMANAGER__

#include "generic/Population.h"
#include "TObject.h"
#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <vector>
#include <iostream>

class HistogramManager : public TObject {
private:
//========= Individual ======//
#ifdef ENABLE_GEANTV
  TH1F *hAllev;     //#1
  TH1F *hBuffev;    //#3
  TH1F *hThread;    //#4
  TH1F *hPriority;  //#4
  TH1F *hSteps;     //#5
  TH1F *hVector;    //#6
  TH1F *hMaxVector; //#7
#else
  TH1F *hx;
#endif

public:
  HistogramManager(){}
  HistogramManager *Instance();
  virtual ~HistogramManager(){}
  Bool_t HistoFill(Population<Double_t> &pop, char *hfile);
  Bool_t CheckValue(ROOT::Internal::TTreeReaderValueBase *value);
  void Reset();

private:
  static HistogramManager *HistoInstance;

  HistogramManager(const HistogramManager &) = delete;
  HistogramManager &operator=(const HistogramManager &) =delete;

  ClassDef(HistogramManager, 1)
};

#endif
