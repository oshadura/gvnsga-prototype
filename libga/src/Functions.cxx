#include "Genes.h"
#include "Population.h"
#include "AlgorithmNSGA.h"
#include "HistogramManager.h"
#include "Functions.h"

#include <vector>
#include <ostream>
#include <string>
#include <utility>

using namespace std::placeholders;

ClassImp(Functions)

    Functions::Functions(const Functions &func)
    : fNParam(func.fNParam), fNCons(func.fNCons), fInterval(func.fInterval),
      fConstraines(func.fConstraines) {
  /*fConstraines.clear();
  for (std::vector<Double_t*>::const_iterator it = fConstraines.begin(),
  itEnd = (fConstraines.end() - GetNParam()); it != itEnd; ++it){
          Double_t value=*(*it);
          fConstraines.push_back(value);
          ++it;}
  fInterval.clear();
  for (std::vector<Double_t*>::const_iterator it = fInterval.begin(),
  itEnd = fInterval.end() - GetNParam(); it != itEnd; ++it){
          Double_t value=*(*it);
          fInterval.push_back(value);
          ++it;}
          */
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

void Functions::SetConstrain(Int_t i, Double_t value) {
  fConstraines.emplace(fConstraines.begin() + i, value);
}

void SetFunction(void (*fFunction)()) { std::function<void()> f = fFunction; }

void SetFunctionGenes(void (*fFunction)(Genes &), Genes &Individual) {
  auto f = std::bind(fFunction, Individual);
}
