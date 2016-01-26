#include "TGenes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"
#include "Functions.h"
#include "invoke_cpp11.hpp"

#include <vector>
#include <ostream>
#include <string>
#include <utility>

using namespace std::placeholders;

ClassImp(Functions)

    Functions::Functions(const Functions &func)
    : fNParam(func.fNParam), fNCons(func.fNCons), fInterval(func.fInterval),
      fNObjectives(func.fNObjectives), fPMut(func.fPMut), fEtaMut(func.fPMut) {}

void Functions::SetInterval() {
  for (Int_t i = 0; i < fNParam; ++i) {
    SetIntervalLimit(i, 1, 100);
  }
}

// Implementation that doesnt allow to change number of parameters
void Functions::SetIntervalLimit(Int_t i, Double_t fMin, Double_t fMax) {
  auto value = std::make_pair(fMin, fMax);
  fInterval.emplace(fInterval.begin() + i, value);
}
