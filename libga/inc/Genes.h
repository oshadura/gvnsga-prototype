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
  // Constructor
  Genes();
  // Copy Constructor
  Genes(const Genes &copy) {} // Copy operator for vector of Individual
  // Parametrized constructor
  Genes(std::vector<Double_t> &f); // Vector of parameters = individual
  // Destructor
  virtual ~Genes() {}
  // Function building Genes (meved in Population and Functions)
  /*
  void Set();
  void SetIt(Int_t i);
  */

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
  const Genes SetGene(Int_t i, Double_t value) {
    fGenes.emplace(fGenes.begin() + i, value);
    return fGenes;
  }
  std::vector<Double_t> GetFitness() const {
    return fFitness;
  } // Get a vector of fTime and fMemory
  Double_t GetFitness(Int_t i) const {
    return fFitness.at(i);
  } // Get Fitness at i position
  void SetFitness(std::vector<Double_t> fitness) { fFitness = fitness; }
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

  ////////// Parameters definition /////////////
  Double_t GetAllev() const { return fAllev; }
  Double_t GetBuffer() const { return fBuffev; }
  Double_t GetThread() const { return fThread; }
  Double_t GetPriority() const { return fPriority; }
  Double_t GetSteps() const { return fSteps; }
  Double_t GetVector() const { return fVector; }
  Double_t GetTime() const { return fTime; }
  Double_t GetMemory() const { return fMemory; }

private:
  ///////////////////////////////////////////////////
  // Individual parts (Genes)
  Int_t fAllev;  // All events (after will be translated in GeantV namespace) #0
  Int_t fBuffev; // Buffered events (after will be translated in GeantV
                 // namespace) #1
  Int_t fThread; // Number of threads (after will be translated in GeantV
                 // namespace) #3
  Double_t fPriority; // Priority value (after will be translated in GeantV
                      // namespace) #4
  Int_t fSteps;       // Number of steps (after will be translated in GeantV
                      // namespace) #5
  Int_t fVector;      // Vector size (after will be translated in GeantV
                      // namespace) #6
  //////////////////////////////////////////////////
  // Parts of fitness vector
  Double_t fTime;   // RT from GeantV (after will be translated in GeantV
                    // namespace)
  Double_t fMemory; // RT from GeantV (after will be translated in GeantV
                    // namespace)
  ///////////////////////////////////////////////////
  std::vector<Double_t> fFitness; // Vector of values of different fitness
                                  // function (objectives)
  Int_t fNObjectives;             // Number of fitness values (objectives)
  Int_t fRank;                    // Rank of Individual
  Int_t fDominationCounter;       // Domination counter for individual (used in
                                  // Non-Dominant sorting)
  Double_t fCrowdingDistance;     // Crowding distance per individual
  Bool_t fEvaluated;              // Evaluated or not
  std::vector<Double_t> fDominated; // Vector of dominanted individuals
  Double_t ConstViol;               // Violation of constraints
  std::vector<Double_t> fGenes;
  Double_t fEpsilonC;

  ClassDef(Genes, 1)
};

#endif
