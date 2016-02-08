#ifndef __POPULATION__
#define __POPULATION__

#include "TGenes.h"
#include "ExceptionMessenger.h"

#ifdef ENABLE_GEANTV
#include "GeantPropagator.h"
#endif

#include "TH1.h"
#include "TFile.h"
#include "TObject.h"

#include <vector>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <memory>
#include <algorithm>

#define EPS 1e-14
#define INF 1e+14

template <class T> class Genes;
template <class T> class Population : public Genes<T>, public Functions {
public:
  Population()
      : fFront(), fPopulation(), fCrowdingObj(true), fSizePop(0), fH(0),
        fPopFunction(NULL), setupPop() {}
  Population(Int_t size)
      : fFront(), fPopulation(), fCrowdingObj(true), fSizePop(size), fH(0),
        fPopFunction(NULL), setupPop() {
    fFront.reserve(size);
    fPopulation.reserve(size);
  }

  Population(const Int_t fSizePop, const Int_t fNParam, const Int_t fNCons,
             const Int_t fNObjectives, const Double_t fEpsilonC,
             const Double_t fPMut, const Double_t fEtaMut,
             const std::vector<std::pair<Double_t, Double_t>> fInterval,
             const Functions::functype func) throw(ExceptionMessenger);

  Population(const Population &pop) {}

  Population<T> &operator=(const Population<T> &pop) {
    if (this != &pop) {
      setupPop = pop.setupPop;
      fPopFunction = pop.fPopFunction;
      fFront = pop.fFront;
      fPopulation = pop.fPopulation;
      fSizePop = pop.fSizePop;
      fH = pop.fH;
    }
    return *this;
  }
  std::vector<Genes<T>> operator=(Population<T> pop) { return fPopulation; }

  virtual ~Population() {}
  /////
  Genes<T> &GetGenes(Int_t i) { return fPopulation.at(i); }
  void SetGenes(Int_t i, const Genes<T> &value) {
    fPopulation.emplace(fPopulation.begin() + i, value);
  }
  void PushGenes(const Genes<T> &value) { fPopulation.push_back(value); }
  /////
  void SetPopulationSize(Int_t s) { fPopulation.resize(s); }
  Int_t GetPopulationSize() { return fPopulation.size(); }
  Int_t GetPopulationSetupSize() const { return fSizePop; }
  /////
  std::vector<std::vector<Int_t>> GetFront() { return fFront; }
  std::vector<Int_t> &GetFront(Int_t i) { return fFront.at(i); }
  /////
  Bool_t IsCrowdingObj() { return fCrowdingObj; }
  void SetCrowdingObj(Bool_t co) { fCrowdingObj = co; }
  /////
  void Build() throw(ExceptionMessenger);
  void CrowdingDistanceAll();
  void CrowdingDistanceFront(Int_t i);
  void FastNonDominantSorting();
  void Merge(const Population &population1,
             const Population &population2); // Merging two populations
  Int_t Mutate();
#ifdef ENABLE_GEANTV
  void Evaluate(GeantPropagator* prop);
#else
  void Evaluate();
#endif
  /////
  void SetPopFunction(Functions::popfunctype f) { fPopFunction = f; }
  /////
  void ResetHistogramPointer() {
    fH = 0;
  } // Function that reset histogram pointer
  TH1F *GetHistogram() const { return fH; } // Return histosgrames
  /////
  void WritePopulationTree(Population &pop, const char *file);
  void UpdatePopulationTree(Population &pop, const char *file);
  void ReadPopulationTree(Population &pop, const char *file);
  Int_t PrintTree(const char *file, const char *name);
  /////
  Genes<T> operator[](Int_t i) { return fPopulation.at(i); }
  /////
  friend std::ostream &operator<<(std::ostream &os, Population<T> &pop) {
    os << "Population: [\n";
    std::ostream_iterator<Genes<T>> fGenesOutIt(os, "\n");
    std::copy(pop.fPopulation.begin(), pop.fPopulation.end(), fGenesOutIt);
    for (auto it = pop.begin(); it != pop.end(); ++it) {
      os << *it;
    }
    os << "]";
    return os;
  }
  // Const <-> non const iterator problem
  void printPopulation(const Population<T> &pop) {
    std::cout << "Population: [\n" << std::endl;
    for (int i = 0; i < const_cast<Population<T> &>(pop).GetPopulationSize();
         ++i) {
      for (auto it = const_cast<Population<T> &>(pop).GetGenes(i).begin();
           it != const_cast<Population<T> &>(pop).GetGenes(i).end(); ++it) {
        std::cout << *it << std::endl;
      }
    }
    std::cout << " ]\n" << std::endl;
  }

public:
  Bool_t fCrowdingObj; // true: crowding over objective (default)
                       // false: crowding over real variable
  Int_t GenCounter;
  std::vector<std::vector<Int_t>> fFront;
  std::vector<Genes<T>> fPopulation;

private:
#ifdef ENABLE_GEANTV
  void Evaluation(GeantPropagator* prop);
  void EvaluationOpenMP(GeantPropagator* prop);
#else
  void Evaluation();
  void EvaluationOpenMP();
#endif

private:
#ifdef ENABLE_GEANTV
  GeantPropagator *prop;
#endif
  Functions setupPop;
  Functions::popfunctype fPopFunction;
  Int_t fSizePop;
  TH1F *fH;

  ClassDef(Population, 1)
};

#endif
