#ifndef __GENES__
#define __GENES__

#include "TObject.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include <vector>

#include "Functions.h"
#include "ExceptionMessenger.h"

#define EPS 1e-14
#define INF 1e+14

class Functions;
class AlgorithmNSGA;

template <class T> class Population;

template <class T> class Genes : public TObject {
public:
  // Constructor
  Genes();
  Genes(const Functions& config) throw (ExceptionMessenger);
  // Copy constructor
  Genes(const Genes<T> &copy) {} // Copy operator for vector of Individual
  // Parametrized constructor
  Genes(std::vector<T> &f); // Vector of parameters = individual
  // Destructor
  virtual ~Genes() {}
  // Function building Genes (moved in Population and Functions)
  void Clear(Option_t *option = "");
  T CheckDominance(const Genes<T> *ind2);
  Int_t Mutate();
  void StoreGenesTree(Genes<T> *ind);
  /*operators of comparision = <  > */

  // Double_t *operator[](Int_t i) { return fGenes[i]; }
  // Bool_t *operator<(Genes& ind1, Genes& ind2) {return ind1.fFitness <
  // ind2.fFitness;}

  Genes<T> &operator=(const Genes<T> &gen);
  // virtual const std::vector<Genes>&  operator[] ( Int_t i ) const { return
  // fGenes.at(i); }
  void Set();
  void Evaluate(Genes<T> &ind);
  Int_t GetDominatedCounter() { return fDominationCounter; }
  std::vector<T> GetDominated() {
    return fDominated;
  } // Vector of dominanted individuals
  T GetDominated(Int_t i) { return fDominated.at(i); } // 1 Dominated individual
  T SetDominatedCounter(Int_t dc) { return fDominationCounter = dc; }
  T GetGene(Int_t i) const { return fGenes.at(i); }
  const Genes<T> SetGene(Int_t i, T value) {
    fGenes.emplace(fGenes.begin() + i, value);
    return fGenes;
  }
  std::vector<T> GetFitnessVector() const {
    return fFitness;
  } // Get a vector of fTime and fMemory
  T GetFitness(Int_t i) const {
    return fFitness.at(i);
  } // Get Fitness at i position
  void SetFitness(std::vector<T> fitness) { fFitness = fitness; }
  T GetCrowdingDistance() const { return fCrowdingDistance; }
  void SetCrowdingDistance(T dist) { fCrowdingDistance = dist; }
  Int_t GetRank() const { return fRank; }
  void SetRank(Int_t rank) { fRank = rank; }
  void WriteGenesTree(Genes<T> &ind, Population<T> &pop, const char *file);
  void UpdateGenesTree(Genes<T> &ind1, Genes<T> &ind2, Population<T> &pop,
                       const char *file);
  void ReadGenesTree(Genes<T> &ind, Population<T> &pop, const char *file);

  void SetEpsilonC(T epsc) { fEpsilonC = epsc; }
  T GetEpsilonC() const { return fEpsilonC; }

  void SetConstrain(Int_t i, T value);
  std::vector<Double_t> GetConstraines() const { return fConstraines; }

  Int_t size() { return fGenes.size(); }

  typename std::vector<T>::iterator begin() {
    typename std::vector<T>::iterator it = fGenes.begin();
    return it;
  }

  typename std::vector<T>::iterator end() {
    typename std::vector<T>::iterator it = fGenes.end();
    return it;
  }

  T operator[](Int_t i) const { return fGenes.at(i); }
  void clear() { fGenes.clear(); }
  void push_back(Int_t i) { return fGenes.push_back(i); }
  ////////// Parameters definition /////////////
  T GetAllev() const { return fAllev; }
  T GetBuffev() const { return fBuffev; }
  T GetThread() const { return fThread; }
  T GetPriority() const { return fPriority; }
  T GetSteps() const { return fSteps; }
  T GetVector() const { return fVector; }
  T GetTime() const { return fTime; }
  T GetMemory() const { return fMemory; }
  /////////////////////////////////////////////

  T GetAllev(Genes<T> &ind) const { return ind.GetGene(0); }

  T GetBuffev(Genes<T> &ind) const { return ind.GetGene(1); }

  T GetThread(Genes<T> &ind) const { return ind.GetGene(2); }

  T GetPriority(Genes<T> &ind) const { return ind.GetGene(3); }

  T GetSteps(Genes<T> &ind) const { return ind.GetGene(4); }

  T GetVector(Genes<T> &ind) const { return ind.GetGene(5); }

  T GetTime(Genes<T> &ind) const { return ind.GetFitness(0); }

  T GetMemory(Genes<T> &ind) const { return ind.GetFitness(1); }

  friend std::ostream &operator<<(std::ostream &os, Genes<T> &g){
       std::copy(g.begin(), g.end(), std::ostream_iterator<T>(os, "\n"));
       return os;
  }
  
  /*
  void printGene(const Genes<T>& g) {
    std::copy(g.begin(), g.end(), std::ostream_iterator<typename
T::value_type>(std::cout, ", "));
}
*/
private:
  ///////////////////////////////////////////////////
  // Individual parts (Genes)

  T fAllev;    // All events (after will be translated in GeantV namespace) #0
  T fBuffev;   // Buffered events (after will be translated in GeantV
               // namespace) #1
  T fThread;   // Number of threads (after will be translated in GeantV
               // namespace) #3
  T fPriority; // Priority value (after will be translated in GeantV
               // namespace) #4
  T fSteps;    // Number of steps (after will be translated in GeantV
               // namespace) #5
  T fVector;   // Vector size (after will be translated in GeantV
               // namespace) #6

  //////////////////////////////////////////////////
  // Parts of fitness vector
  T fTime;   // RT from GeantV (after will be translated in GeantV
             // namespace)
  T fMemory; // RT from GeantV (after will be translated in GeantV
             // namespace)

  ///////////////////////////////////////////////////

  std::vector<T> fFitness;    // Vector of values of different fitness
                              // function (objectives)
  Int_t fDominationCounter;   // Domination counter for individual (used in
                              // Non-Dominant sorting)
  Int_t fRank;                // Rank of Individual
  Double_t fCrowdingDistance; // Crowding distance per individual
  Bool_t fEvaluated;          // Evaluated or not
  std::vector<T> fDominated;  // Vector of dominanted individuals
  Double_t ConstViol;         // Violation of constraints
  std::vector<T> fGenes;
  Double_t fEpsilonC;
  std::vector<T> fConstraines; // Vector of constraines for NSGA2


  /////// Trial for Genes() constructor ///////////////
  const Functions* setup;
  ///////////////////////////////////////////////////

  // Should we write a map to be sure about connection between Limits[] and
  // Genes[] || Fitmess[] and Constraint[]?

  // static std::multimap<Genes,Limits> fInputMap;
  // static std::multimap<Genes,Constraint> fOutputMap;

  ClassDef(Genes, 1)
};

#endif
