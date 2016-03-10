#ifndef __GENES__
#define __GENES__

#define RESET "\033[0m"
#define BLACK "\033[30m"   /* Black */
#define RED "\033[31m"     /* Red */
#define GREEN "\033[32m"   /* Green */
#define YELLOW "\033[33m"  /* Yellow */
#define BLUE "\033[34m"    /* Blue */
#define MAGENTA "\033[35m" /* Magenta */
#define CYAN "\033[36m"    /* Cyan */

#include "TObject.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

//#ifdef ENABLE_GEANTV
//#include "GeantPropagator.h"
//#endif

#include <vector>

#include "Functions.h"
#include "ExceptionMessenger.h"

//#ifdef ENABLE_GEANTV
//  class GeantPropagator;
//#endif

class Functions;
template <class T> class Population;
template <class T> class Genes : public TObject {

protected:
  T fGeneType;

public:
  // Constructor
  Genes() throw();
  Genes(const Functions &config) throw(ExceptionMessenger);
  // Copy constructor
  Genes(const Genes<T> &copy) =
      default; // Copy operator for vector of Individual
  // Destructor
  virtual ~Genes() {}
  // Function building Genes (moved in Population and Functions)
  void Clear(Option_t *option = "");
  T CheckDominance(Functions *setup,
                   const Genes<T> *ind2) throw(ExceptionMessenger);
  Int_t Mutate(const Functions *setup);
  void StoreGenesTree(Genes<T> *ind);
  Genes<T> &operator=(const Genes<T> &gen);
  void Set(Functions &setup, Genes<T> &ind) throw(ExceptionMessenger);

#ifdef ENABLE_GEANTV
  void SetGeantV(Functions &setup, Genes<T> &ind) throw(ExceptionMessenger);
#endif

  //#ifdef ENABLE_GEANTV
  //  void Evaluate(GeantPropagator *prop, Functions &setup,
  //                      Genes<T> &ind) throw(ExceptionMessenger);
  //#else
  void Evaluate(Functions &setup, Genes<T> &ind) throw(ExceptionMessenger);
  //#endif

  Int_t GetDominatedCounter() { return fDominationCounter; }
  std::vector<Int_t> GetDominated() {
    return fDominated;
  } // Vector of dominanted counters
  T GetDominated(Int_t i) { return fDominated.at(i); }
  T SetDominatedCounter(Int_t dc) { return fDominationCounter = dc; }
  T GetGene(Int_t i) const { return fGenes.at(i); }
  //void SetGene(Int_t i, T value) { fGenes[i] = value; }
  void SetGene(Int_t i, T value) {fGenes.at(i) = value;}
  std::vector<T> GetFitnessVector() const {
    return fFitness;
  } // Get a vector of fTime and fMemory
  T GetFitness(Int_t i) const {
    return fFitness.at(i);
  } // Get Fitness at i position
  void SetFitness(std::vector<T> fitness) { fFitness = fitness; }
  void SetFitness(Int_t i, T value) { fFitness[i] = value; }
  T GetCrowdingDistance() const { return fCrowdingDistance; }
  void SetCrowdingDistance(T dist) { fCrowdingDistance = dist; }
  Int_t GetRank() const { return fRank; }
  void SetRank(Int_t rank) { fRank = rank; }
  Double_t GetConsViol() const { return ConstViol; }
  void SetConsViol(Double_t consv) { ConstViol = consv; }
  void SetEvaluated(Bool_t value) { fEvaluated = value; }

  void WriteGenesTree(Genes<T> &ind, Population<T> &pop, const char *file);
  void UpdateGenesTree(Genes<T> &ind1, Genes<T> &ind2, Population<T> &pop,
                       const char *file);
  void ReadGenesTree(Genes<T> &ind, Population<T> &pop, const char *file);
  void SetConstrain(Int_t i, T value);
  std::vector<T> GetConstraines() const { return fConstraines; }
  T GetConstrain(Int_t i) const { return fConstraines.at(i); }
  std::vector<T> GetfGenes() const { return fGenes; }
  const Functions *GetSetup() const { return setup; }
  void SetDominated(std::vector<Int_t> &d) { fDominated = d; }
  // faster access
  Int_t capacity() { return fGenes.capacity(); }
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
  void push_back(T i) { return fGenes.push_back(i); }
  void resize(Int_t i) { return fGenes.resize(i); }
  void reserve(Int_t i) { return fGenes.reserve(i); }
  Bool_t empty() { return fGenes.empty(); }
  T at(Int_t i){ return fGenes.at(i);}
  void assign(Int_t i, T value){return fGenes.assign(i, value);}
  // Pointless thing
  T GetAllev() const { return fAllev; }
  T GetBuffev() const { return fBuffev; }
  T GetThread() const { return fThread; }
  T GetPriority() const { return fPriority; }
  T GetSteps() const { return fSteps; }
  T GetVector() const { return fVector; }
  T GetTime() const { return fTime; }
  T GetMemory() const { return fMemory; }
  //
  T GetAllev(Genes<T> &ind) const { return ind.GetGene(0); }
  T GetBuffev(Genes<T> &ind) const { return ind.GetGene(1); }
  T GetThread(Genes<T> &ind) const { return ind.GetGene(2); }
  T GetPriority(Genes<T> &ind) const { return ind.GetGene(3); }
  T GetSteps(Genes<T> &ind) const { return ind.GetGene(4); }
  T GetVector(Genes<T> &ind) const { return ind.GetGene(5); }
  T GetMaxVector(Genes<T> &ind) const { return ind.GetGene(6); }
  T GetTime(Genes<T> &ind) const { return ind.GetFitness(0); }
  T GetMemory(Genes<T> &ind) const { return ind.GetFitness(1); }
  //
  void printGenes(Genes<T> &g) {

    std::cout << "Individual rank = " << g.GetRank() << "\n" << std::endl;

    std::cout << "Available constraint violations = " << g.GetConsViol() << "\n"
              << std::endl;
    std::cout << "Size of Gene<T> is " << g.capacity() << "\n" << std::endl;
    std::cout << "Gene<T> = [";
    for (auto it = g.begin(); it != g.end(); ++it) {
      std::cout << *it << ' ';
    }
    std::cout << "]"
              << "\n" << std::endl;
    std::cout << "Fitness<T> = [";
    for (auto &it : g.GetFitnessVector()) {
      std::cout << it << ' ';
    }
    std::cout << "]"
              << "\n" << std::endl;
    std::cout << "Constraint<T> = [";
    for (auto &it : g.GetConstraines()) {
      std::cout << it << ' ';
    }
    std::cout << "]"
              << "\n" << std::endl;
    std::cout << "fCrowdingDistance = " << g.GetCrowdingDistance() << "\n"
              << std::endl;
  }

private:
  // Individual parts (Genes)
  T fAllev;     // All events (after will be translated in GeantV namespace) #0
  T fBuffev;    // Buffered events (after will be translated in GeantV
                // namespace) #1
  T fThread;    // Number of threads (after will be translated in GeantV
                // namespace) #3
  T fPriority;  // Priority value (after will be translated in GeantV
                // namespace) #4
  T fSteps;     // Number of steps (after will be translated in GeantV
                // namespace) #5
  T fVector;    // Vector size (after will be translated in GeantV
                // namespace) #6
  T fMaxVector; // Max VectorSize
  // Parts of fitness vector
  T fTime;                 // RT from GeantV (after will be translated in GeantV
                           // namespace)
  T fMemory;               // RT from GeantV (after will be translated in GeantV
                           // namespace)
  std::vector<T> fFitness; // Vector of values of different fitness
                           // function (objectives)
  Int_t fDominationCounter;      // Domination counter for individual (used in
                                 // Non-Dominant sorting)
  Int_t fRank;                   // Rank of Individual
  Double_t fCrowdingDistance;    // Crowding distance per individual
  Bool_t fEvaluated;             // Evaluated or not
  std::vector<Int_t> fDominated; // Vector of dominanted counters
  Double_t ConstViol;            // Violation of constraints
  std::vector<T> fGenes;
  std::vector<T> fConstraines; // Vector of constraines for NSGA2
  const Functions *setup;

  //#ifdef ENABLE_GEANTV
  //  GeantPropagator *prop;
  //#endif
  // Should we write a map to be sure about connection between Limits[] and
  // Genes[] || Fitmess[] and Constraint[]?
  // static std::multimap<Genes,Limits> fInputMap;
  // static std::multimap<Genes,Constraint> fOutputMap;

  ClassDef(Genes, 1)
};

#endif
