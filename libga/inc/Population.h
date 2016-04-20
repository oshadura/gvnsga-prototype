#ifndef __POPULATION__
#define __POPULATION__

#define RESET "\033[0m"
#define BLACK "\033[30m"   /* Black */
#define RED "\033[31m"     /* Red */
#define GREEN "\033[32m"   /* Green */
#define YELLOW "\033[33m"  /* Yellow */
#define BLUE "\033[34m"    /* Blue */
#define MAGENTA "\033[35m" /* Magenta */
#define CYAN "\033[36m"    /* Cyan */

#include "TGenes.h"
#include "ExceptionMessenger.h"
#include "HistogramManager.h"

//#ifdef ENABLE_GEANTV
//#include "GeantPropagator.h"
//#endif

#include "TH1.h"
#include "TFile.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"

#include <vector>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <memory>
#include <algorithm>
#include <ios>

#define EPS 1e-14
#define INF 1e+14

template <class T> class Genes;
template <class T>
class Population : public Genes<T>, public Functions, public HistogramManager {
public:
  Population()
      : fFront(), fPopulation(), fCrowdingObj(true), fSizePop(0), fHisto(0),
        fPopFunction(NULL), setupPop() {}
  Population(Int_t size)
      : fFront(), fPopulation(), fCrowdingObj(true), fSizePop(size), fHisto(0),
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
      fHisto = pop.fHisto;
    }
    return *this;
  }
  std::vector<Genes<T>> operator=(Population<T> pop) { return fPopulation; }

  virtual ~Population() { delete fHisto; }

  Genes<T> &GetGenes(Int_t i) { return fPopulation.at(i); }
  void SetGenes(Int_t i, const Genes<T> &value) {
    fPopulation.emplace(fPopulation.begin() + i, value);
  }
  void PushGenes(const Genes<T> &value) { fPopulation.push_back(value); }
  void SetPopulationSize(Int_t s) { fPopulation.resize(s); }
  Int_t GetPopulationSize() { return fPopulation.size(); }
  Int_t GetPopulationSetupSize() const { return fSizePop; }
  std::vector<std::vector<Int_t>> GetFront() { return fFront; }
  std::vector<Int_t> &GetFront(Int_t i) { return fFront.at(i); }
  Bool_t IsCrowdingObj() { return fCrowdingObj; }
  void SetCrowdingObj(Bool_t co) { fCrowdingObj = co; }
  void Build() throw(ExceptionMessenger);
  void CrowdingDistanceAll();
  void CrowdingDistanceFront(Int_t i);
  void FastNonDominantSorting();
  void Merge(const Population &population1,
             const Population &population2); // Merging two populations
  Int_t Mutate();
  void Print();
  void LoadCSV(const Population<T> &pop);
  std::ofstream &CreateCVS(std::string file);
  void CSVOutput(std::ofstream &populationcvs, const Population<T> &pop);
  //#ifdef ENABLE_GEANTV
  //  void Evaluate(GeantPropagator* prop);
  //#else
  void Evaluate();
  //#endif
  void SetPopFunction(Functions::popfunctype f) { fPopFunction = f; }

  void ResetHistogramPointer() {
    fHisto->Reset();
  } // Function that reset histogram pointer
  HistogramManager *GetHistograms() const {
    return fHisto;
  } // Return histosgrames
  void WritePopulationTree(Population &pop, const char *file);
  void UpdatePopulationTree(Population &pop, const char *file);
  void ReadPopulationTree(Population &pop, const char *file);
  Int_t PrintTree(const char *file, const char *name);
  void Store(const char *file, const Population<T> &pop);
  Genes<T> operator[](Int_t i) { return fPopulation.at(i); }

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
  //
  void push_back(Genes<T> &i) { return fPopulation.push_back(i); }

public:
  Bool_t fCrowdingObj; // true: crowding over objective (default)
                       // false: crowding over real variable
  Int_t GenCounter;
  std::vector<std::vector<Int_t>> fFront;
  std::vector<Genes<T>> fPopulation;

private:
  //#ifdef ENABLE_GEANTV
  //  void Evaluation(GeantPropagator* prop);
  //  void EvaluationOpenMP(GeantPropagator* prop);
  //#else
  void Evaluation();
  void EvaluationOpenMP();
  void EvaluationMPI();
  void EvaluationSLURM();
  //#endif

private:
  //#ifdef ENABLE_GEANTV
  //  GeantPropagator *prop;
  //#endif
  Functions setupPop;
  Functions::popfunctype fPopFunction;
  Int_t fSizePop;
  HistogramManager *fHisto;

  ClassDef(Population, 1)
};

#endif
