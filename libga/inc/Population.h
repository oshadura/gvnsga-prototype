#ifndef __POPULATION__
#define __POPULATION__

#include "TGenes.h"
#include "ExceptionMessenger.h"

#include "TH1.h"
#include "TFile.h"

#include <vector>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <memory>
#include <algorithm>

template <class T> class Genes;
template <class T> class Population : public Genes<T>, public Functions {

  typedef Genes<T> Individual;
  // To be used after
  //template <typename ...Args> using myType = std::function<void(Args...)>;

protected:
public:
  Population()
      : fFront(), fPopulation(), fCrowdingObj(true), fSizePop(0), fH(0), fPopFunction(NULL), setupPop() {}
  Population(Int_t size)
      : fFront(), fPopulation(), fCrowdingObj(true), fSizePop(size), fH(0), fPopFunction(NULL), setupPop() {
    fFront.reserve(fSizePop);
    fPopulation.reserve(fSizePop);
  }

  Population(const Int_t fSizePop,
    const Int_t fNParam,
    const Int_t fNCons,
    const Int_t fNObjectives,
    const Double_t fEpsilonC,
    const Double_t fPMut,
    const Double_t fEtaMut,
    const std::vector<std::pair<Double_t, Double_t>> fInterval,
    const Functions::functype func) 
  throw (ExceptionMessenger);

  Population(const Population& pop){}

  Population<T> &operator=(const Population<T> &pop){
    if (this != &pop) {
      genes = pop.genes;
      setupPop = pop.setupPop;
      fPopFunction = pop.fPopFunction;
      fFront = pop.fFront;
      fPopulation = pop.fPopulation;
      fSizePop = pop.fSizePop;
      fH = pop.fH;
    }
    return *this;
  }

  virtual ~Population(){}
  Individual &GetGenes(Int_t i) { return fPopulation.at(i); }
  void SetGenes(Int_t i, const Genes<T> &value) {
    fPopulation.emplace(fPopulation.begin() + i, value);
  }
  void PushGenes(const Individual &value) { fPopulation.push_back(value); }
  void SetPopulationSize(Int_t s) { fPopulation.resize(s); }
  Int_t GetPopulationSize() const { return fPopulation.size(); }
  Int_t GetPopulationSetupSize() const { return fSizePop; }
  std::vector<Individual> GetIndividuals() const { return fPopulation; }
  std::vector<Individual> GetFront() const { return fFront; }
  Individual GetFront(Int_t i) { return fFront.at(i); }
  Bool_t IsCrowdingObj() { return fCrowdingObj; }
  std::vector<Individual> operator=(Population<T> pop) { return fPopulation; }
  void SetCrowdingObj(Bool_t co){fCrowdingObj = co;}
  void Build() throw (ExceptionMessenger);
  void CrowdingDistanceAll();
  void CrowdingDistanceFront(Int_t i);
  //void FastNonDominantSorting(const Population<T> &population);
  void FastNonDominantSorting();
  void Merge(const Population &population1,
             const Population &population2); // Merging two populations
  Int_t Mutate();
  void Clear(Option_t *option = "");         // Clear function
  static void Reset(Option_t *option = "");  // Reset function
  //void Evaluate(const Population<T>& pop);
  void Evaluate();
  void SetPopFunction(Functions::popfunctype f){
    fPopFunction = f;
  }

  void ResetHistogramPointer() {
    fH = 0;
  } // Function that reset histogram pointer
  TH1F *GetHistogram() const { return fH; } // Return histosgrames
  //////////////////// Playing with ROOT files///////////////
  void WritePopulationTree(Population &pop, const char *file);
  void UpdatePopulationTree(Population &pop, const char *file);
  void ReadPopulationTree(Population &pop, const char *file);
  Int_t PrintTree(const char *file, const char *name);
  /*
  friend std::ostream &operator<<(std::ostream &os, Population<T> &pop){
    os << "Population: [\n";
    //std::ostream_iterator<Genes<T>> fGenesOutIt (os,"\n");
    //std::copy(pop.GetIndividuals().begin(), pop.GetIndividuals().end(),fGenesOutIt);
    for(auto it = pop.begin(); it != pop.end(); ++it){
      os << *it;
    }
    os << "]";
    return os;
  }
  */
  void printPopulation(const Population<T>& pop){
    std::cout << "Population: [\n" << std::endl;
    for(auto it = pop.begin(); it != pop.end(); ++it){
      std::cout << *it << std::endl;
    }
    std::cout << " ]\n" << std::endl;
  }
  
public:
  Bool_t fCrowdingObj; // true: crowding over objective (default)
                       // false: crowding over real variable
  Int_t GenCounter;

private:
  Genes<T> genes;
  Functions setupPop;
  Functions::popfunctype fPopFunction;
  std::vector<Individual> fFront;
  std::vector<Individual> fPopulation;
  Int_t fSizePop;
  TH1F *fH;

  //Counter
  ClassDef(Population, 1)
};

#endif
