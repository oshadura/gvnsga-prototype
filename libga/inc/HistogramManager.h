#ifndef __HISTOGRAMMANAGER__
#define __HISTOGRAMMANAGER__

#include "Population.h"
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
  //========== TFile ===========//
  // TFile *fHisto;
  //========= Individual ======//
  TH1F *hAllev;     //#1
  TH1F *hBuffev;    //#3
  TH1F *hThread;    //#4
  TH1F *hPriority;  //#4
  TH1F *hSteps;     //#5
  TH1F *hVector;    //#6
  TH1F *hMaxVector; //#7
  //========= Fitness ========//
  TH1F *hMemory; //#1
  TH1F *hTime;   //#2

public:
  HistogramManager();
  HistogramManager *Instance();
  virtual ~HistogramManager();
  Bool_t HistoFill(Population<Double_t> &pop, char *file);
  Bool_t CheckValue(ROOT::Internal::TTreeReaderValueBase *value);
  void Reset();

private:
  static HistogramManager *HistoInstance;

  ClassDef(HistogramManager, 1)
};

#endif
