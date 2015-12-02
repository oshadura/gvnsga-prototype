#ifndef __POPULATION__
#define __POPULATION__

#include "AlgorithmNSGA.h"

#include "TH1.h"
#include "TFile.h"

#include <vector>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>

#define EPS 1e-14
#define INF 1e+14

class Genes;

class Population : public Genes {
public:
  Population() : fFront(), fPopulation(), fCrowdingObj(1), fNGen(1) {}
  virtual ~Population() {}

  Genes &GetGenes(Int_t i) { return fPopulation.at(i); }
  void SetGenes(Int_t i, const Genes &value) {
    fPopulation.emplace(fPopulation.begin() + i, value);
  }
  void PushGenes(const Genes &value) { fPopulation.push_back(value); }
  Int_t GetGenerationNumber() { return fNGen; }
  void SetGenerationNumber(Int_t gen) { fNGen = gen; }
  Int_t GetPopulationSize() const { return fPopulation.size(); }
  std::vector<Genes> &GetIndividuals() { return fPopulation; }
  // std::vector<std::vector<Double_t> > GetFront() { return fFront;} //
  // std::vector<Genes>
  // std::vector<Double_t> GetFront(Int_t i) { return fFront.at(i); } // Genes
  std::vector<Genes> GetFront() { return fFront; }
  Genes GetFront(Int_t i) { return fFront.at(i); }
  Bool_t IsCrowdingObj() { return fCrowdingObj; }

  std::vector<Genes> operator=(std::vector<Genes> fPopulation) {
    return fPopulation;
  }
  void SetCrowdingObj();
  void Build();
  void CrowdingDistanceAll();
  void CrowdingDistanceFront(Int_t i);
  void FastNonDominantSorting();
  void Evaluation() { ; }
  void Merge(const Population &population1,
             const Population &population2); // Merging two populations
  void Clear(Option_t *option = "");         // Clear function
  static void Reset(Option_t *option = "");  // Reset function

  void ResetHistogramPointer() {
    fH = 0;
  } // Function that reset histogram pointer
  TH1F *GetHistogram() const { return fH; } // Returns historgrames
  void WritePopulationTree(Population &pop, const char *file);
  void UpdatePopulationTree(Population &pop, const char *file);
  void ReadPopulationTree(Population &pop, const char *file);
  Int_t PrintTree(const char *file, const char *name);
  void Evaluate();

public:
  Bool_t fCrowdingObj; // true: crowding over objective (default) false:
                       // crowding over real variable
private:
  // std::vector<std::vector<Double_t> > fFront; // std::vector<Genes>
  std::vector<Genes> fFront;
  std::vector<Genes> fPopulation;
  Int_t fNGen;
  TH1F *fH;

  ClassDef(Population, 1)
};

#endif
