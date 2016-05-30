#ifndef __HISTOGRAMMANAGER__
#define __HISTOGRAMMANAGER__

//===--- HistogramManager.h - LibGA ---------------------------------*- C++
//-*-===//
//
//                     LibGA Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file HistogramManager.h
 * @brief Implementation of CSV class for LibGA
 * prototype
 */
//

#include "generic/Population.h"
#include "TObject.h"
#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <vector>
#include <iostream>

namespace geantvmoop{

template <typename F> class HistogramManager : public TObject {

private:
  F problem; 
public:
    static HistogramManager &GetInstance() {
    static HistogramManager Instance;
    return Instance;
  }

  HistogramManager(HistogramManager const &) = delete;
  void operator=(HistogramManager const &) = delete;
  bool HistoFill(Population<F> &pop, char *hfile);
  void Reset();

private:
  HistogramManager() {}
  ~HistogramManager() {}

//ClassDef(HistogramManager, 1)
ClassDefT(HistogramManager<F>,0)
};

ClassDefT2(HistogramManager,F)

}

#endif
