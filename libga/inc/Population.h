#ifndef __POPULATION__
#define __POPULATION__

#include "AlgorithmNSGA.h"
#include "TGenes.h"

#include "TH1.h"
#include "TFile.h"

#include <vector>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <memory>
#include <algorithm>

#define EPS 1e-14
#define INF 1e+14

//template<class T> class Genes;

template<class T> class Population : public Genes<T> {

  typedef Genes<T> Individual;

public:
  Population()
      : fFront(), fPopulation(), fCrowdingObj(1), fSizePop(0), fGen(0), fH() {}
  Population(Int_t size)
      : fFront(), fPopulation(), fCrowdingObj(1), fSizePop(size), fGen(0),
        fH() {
    fFront.reserve(fSizePop);
    fPopulation.reserve(fSizePop);
  }
  virtual ~Population() {}

  Individual &GetGenes(Int_t i) { return fPopulation.at(i); }
  void SetGenes(Int_t i, const Genes<T> &value) {
    fPopulation.emplace(fPopulation.begin() + i, value);
  }
  void PushGenes(const Genes<T> &value) { fPopulation.push_back(value); }
  void SetPopulationSize(Int_t s) {
    // std::unique_ptr<Population> pop(new Population());
    fPopulation.resize(s);
  }
  Int_t GetPopulationSize() const { return fPopulation.size(); }
  Int_t GetPopulationSetupSize() const { return fSizePop; }
  std::vector<Individual> GetIndividuals() const { return fPopulation; }
  std::vector<Individual> GetFront() const { return fFront; }
  Individual GetFront(Int_t i) { return fFront.at(i); }
  Bool_t IsCrowdingObj() { return fCrowdingObj; }
  std::vector<Individual> operator=(Population<T> pop) {
    return fPopulation;
  }
  // void SetCrowdingObj();
  void Build();
  void CrowdingDistanceAll();
  void CrowdingDistanceFront(Int_t i);
  void FastNonDominantSorting();
  void Merge(const Population &population1,
             const Population &population2); // Merging two populations
  void Clear(Option_t *option = "");         // Clear function
  static void Reset(Option_t *option = "");  // Reset function

  void ResetHistogramPointer() {
    fH = 0;
  } // Function that reset histogram pointer
  TH1F *GetHistogram() const { return fH; } // Returns historgrames
  //////////////////////////////// Playing with ROOT files
  ///////////////////////////////
  void WritePopulationTree(Population &pop, const char *file);
  void UpdatePopulationTree(Population &pop, const char *file);
  void ReadPopulationTree(Population &pop, const char *file);
  Int_t PrintTree(const char *file, const char *name);
  ////////////////////////////////////////////////////////////
  void Evaluate();
  void SetGenNumber(Int_t i) { fGen = i; }
  Int_t GetGenNumber() const { return fGen; }
  ///////////////////////////////////////////////////////////
  /*
  friend std::ostream &operator<<(std::ostream &os, Population<T> &pop){
    std::ostream_iterator<Genes<T>> fGenesOutIt (os,"\n");
    std::copy(pop.GetIndividuals().begin(), pop.GetIndividuals().end(), fGenesOutIt);
    return os;
  }
  */
  //void printPopulation(const Population<T>& p);

public:
  Bool_t fCrowdingObj; // true: crowding over objective (default)
                       // false: crowding over real variable
private:
  std::vector<Individual> fFront;
  std::vector<Individual> fPopulation;
  Int_t fSizePop;
  Int_t fGen; // Generation counter
  TH1F *fH;

  ClassDef(Population, 1)
};

#endif
