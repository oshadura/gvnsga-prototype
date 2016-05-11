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

template <typename F> class Genes;
template <typename F> using individual_t = std::shared_ptr<Genes<F>>;

template <typename F> class Genes {

private:
  typename F::Input fGenes;
  typename F::Output fFitness;

public:
  Genes() {}

  ~Genes() {}

  Genes(Genes const &);

  void operator=(Genes const &);

  Genes(const typename F::Input &i, bool fEvaluated = true) : fGenes(i) {
    if (fEvaluated)
      Evaluate();
  };

  bool IsDominating(const Genes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (fFitness[i] > other.fFitness[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsDominated(const Genes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (fFitness[i] < other.fFitness[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsEqual(const Genes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (other.fFitness[i] != fFitness[i])
        return false;
    }
    return true;
  }

  void Evaluate() { fFitness = F::Evaluate(fGenes); }

  const typename F::Input &GetInput() const { return fGenes; }

  const typename F::Output &GetOutput() const { return fFitness; }
};

#endif
