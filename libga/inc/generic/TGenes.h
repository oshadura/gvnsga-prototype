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

#include <vector>
#include <memory>

#include "generic/ExceptionMessenger.h"

// template <typename F> using GenesPtr = std::shared_ptr<Genes<F> >;

template <typename F> class Genes {

private:
  typename F::Input fGenes;
  typename F::Output fFitness;

public:
private:



/*
T GetAllev(Genes<T> &ind) const { return ind.GetGene(0); }
T GetBuffev(Genes<T> &ind) const { return ind.GetGene(1); }
T GetThread(Genes<T> &ind) const { return ind.GetGene(2); }
T GetPriority(Genes<T> &ind) const { return ind.GetGene(3); }
T GetSteps(Genes<T> &ind) const { return ind.GetGene(4); }
T GetVector(Genes<T> &ind) const { return ind.GetGene(5); }
T GetMaxVector(Genes<T> &ind) const { return ind.GetGene(6); }
T GetTime(Genes<T> &ind) const { return ind.GetFitness(0); }
T GetMemory(Genes<T> &ind) const { return ind.GetFitness(1); }

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
T fTime;                 // RT from GeantV (after will be translated in GeantV
                     // namespace)
T fMemory;               // RT from GeantV (after will be translated in GeantV
                     // namespace)
*/
	ClassDef(Genes, 1)
};

#endif