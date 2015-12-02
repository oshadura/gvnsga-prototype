#ifndef __GENES__
#define __GENES__

#include "TObject.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <vector>
#include "Functions.h"

class Functions;
class Population;
class AlgorithmNSGA;

class Genes : public TObject {
public:
  Genes();
  Genes(const Genes &copy) { ; }   // Copy operator for vector of Individual
  Genes(std::vector<Double_t> &f); // Vector of parameters = individual
  virtual ~Genes() { ; }
  void Set(/*Double_t fAllev, Double_t fBuffev, Double_t fThread,
           Double_t fPriority, Double_t fSteps, Double_t fVector*/);
  void SetIt(Int_t i);

  void Clear(Option_t *option = "");
  Double_t CheckDominance(const Genes *ind2);
  Double_t Mutate();
  void StoreGenesTree(Genes *ind);
  /*operators of comparision = <  > */

  // Double_t *operator[](Int_t i) { return fGenes[i]; }
  // Bool_t *operator<(Genes& ind1, Genes& ind2) {return ind1.fFitness <
  // ind2.fFitness;}

  Genes &operator=(const Genes &gen);
  // virtual const std::vector<Genes>&  operator[] ( Int_t i ) const { return
  // fGenes.at(i); }

  Int_t GetDominatedCounter() { return fDominationCounter; }
  std::vector<Double_t> GetDominated() {
    return fDominated;
  } // Vector of dominanted individuals
  Double_t GetDominated(Int_t i) {
    return fDominated.at(i);
  } // 1 Dominated individual
  Int_t SetDominatedCounter(Int_t dc) { return fDominationCounter = dc; }

  // Double_t GetGene(Int_t i) const { return fGenes.at(i); }
  const Genes &SetGene(Int_t i, Double_t value) {
    fGenes.emplace(fGenes.begin() + i, value);
    return fGenes;
  }

  std::vector<Double_t> GetFitness() const {
    return fFitness;
  } // Get a vector of fTime and fMemory
  Double_t GetFitness(Int_t i) const {
    return fFitness.at(i);
  } // Get Fitness at i position
  void SetFitness(std::vector<Double_t> fitness) {
    fFitness = fitness;
  } // How to add two fTime and fMemory in Vector
  Int_t GetNObjectives() const { return fNObjectives; }
  Double_t GetCrowdingDistance() const { return fCrowdingDistance; }
  void SetCrowdingDistance(Double_t dist) { fCrowdingDistance = dist; }
  Int_t GetRank() const { return fRank; }
  void SetRank(Int_t rank) { fRank = rank; }
  void WriteGenesTree(Genes &ind, Population &pop, const char *file);
  void UpdateGenesTree(Genes &ind1, Genes &ind2, Population &pop,
                       const char *file);
  void ReadGenesTree(Genes &ind, Population &pop, const char *file);

  void SetEpsilonC(Double_t epsc) { fEpsilonC = epsc; }
  Double_t GetEpsilonC() const { return fEpsilonC; }
  void EvaluateGene();

  Int_t size() const { return fGenes.size(); }
  ////// Horrible /////
  std::vector<Double_t>::iterator begin() {
    std::vector<Double_t>::iterator it = fGenes.begin();
    return it;
  }

  std::vector<Double_t>::iterator end() {
    std::vector<Double_t>::iterator it = fGenes.end();
    return it;
  }

  virtual const Double_t operator[](Int_t i) const { return fGenes.at(i); }

  void clear() { fGenes.clear(); }
  void push_back(Int_t i) { return fGenes.push_back(i); }

private:
  std::vector<Double_t> fFitness; // Vector of values of different fitness
                                  // function (objectives)
  Int_t fNObjectives;             // Number of fitness values (objectives)
  Int_t fRank;                    // Rank of Individual
  Int_t fDominationCounter;       // Domination counter for individual (used in
                                  // Non-Dominant sorting)
  Double_t fCrowdingDistance;     // Crowding distance per individual
  Bool_t fEvaluated;              // Evaluated or not
  std::vector<Double_t> fDominated; // Vector of dominanted individuals
  Double_t ConstViol;               // Violation of constrains
  std::vector<Double_t> fGenes;
  Double_t fEpsilonC;

  ClassDef(Genes, 1)
};

#endif
