#include "generic/TGenes.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "algorithms/AlgorithmNSGA.h"
#include "output/HistogramManager.h"

#include <vector>
#include <ostream>
#include <string>
#include <utility>

//using namespace std::placeholders;

ClassImp(Functions)

    Functions::Functions(const Functions &func)
    : fNParam(func.fNParam), fNCons(func.fNCons), fInterval(func.fInterval),
      fNObjectives(func.fNObjectives), fPMut(func.fPMut), fEtaMut(func.fPMut) {}

void Functions::SetInterval() {
  for (Int_t i = 0; i < fNParam; ++i) {
    SetIntervalLimit(i, 0, 1);
  }
}

#ifdef ENABLE_GEANTV
void Functions::SetIntervalGeantV() {
  // for GetAllev(Genes<T> &ind) const { return ind.GetGene(0); }
  SetIntervalLimit(0, 1, 100);
  // for GetBuffev(Genes<T> &ind) const { return ind.GetGene(1); }
  // FIX -> get max value from (GetValue[0] - 1)
  SetIntervalLimit(1, 1, 99);
  // for GetThread(Genes<T> &ind) const { return ind.GetGene(2); }
  SetIntervalLimit(2, 1, 1);
  // for GetPriority(Genes<T> &ind) const { return ind.GetGene(3); }
  SetIntervalLimit(3, 0, 0.1);
  // for GetSteps(Genes<T> &ind) const { return ind.GetGene(4); }
  SetIntervalLimit(4, 1, 10000);
  // for T GetVector(Genes<T> &ind) const { return ind.GetGene(5); }
  SetIntervalLimit(5, 1, 64);
  // for GetMaxVector(Genes<T> &ind) const { return ind.GetGene(6); }
  SetIntervalLimit(6, 1, 512);
}
#endif

// Implementation that doesnt allow to change number of parameters
void Functions::SetIntervalLimit(Int_t i, Double_t fMin, Double_t fMax) {
  auto value = std::make_pair(fMin, fMax);
  fInterval.emplace(fInterval.begin() + i, value);
}
