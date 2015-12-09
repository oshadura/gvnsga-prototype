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
Functions *Functions::fgFunction = 0;

ClassImp(Functions)

    Functions::Functions(const Functions &func)
    : fNParam(func.fNParam), fNCons(func.fNCons), fInterval(func.fInterval), fNObjectives(func.fNObjectives) {
}

Functions *Functions::Instance() {
  if (!fgFunction)
    Functions *fgFunction = new Functions();
  return fgFunction;
}

void Functions::SetInterval() {
  for (Int_t i = 0; i < GetNParam(); ++i) {
    SetIntervalLimit(i, 1, 100);
  }
}
// Implementation that doesnt allow to change number of parameters
void Functions::SetIntervalLimit(Int_t i, Double_t fMin, Double_t fMax) {
  auto value = std::make_pair(fMin, fMax);
  fInterval.emplace(fInterval.begin() + i, value);
}

/*
void SetFunction(void (*fFunction)()) { std::function<void()> f = Function; }

void SetFunctionGenes(void (*fFunction)(Genes<Double_t> &),
                      Genes<Double_t> &ind) {
  auto f = std::bind(fFunction, ind);
  functional::cpp11::invoke(fFunction, ind);
}
*/
