#ifndef __GENES__
#define __GENES__

#include <vector>
#include <memory>

#include "generic/ExceptionMessenger.h"

namespace geantvmoop{

template <typename F> class Genes;
template <typename F> using individual_t = std::shared_ptr<Genes<F>>;

template <typename F> class Genes {

private:
  typename F::Input fGenes;
  typename F::Output fFitness;

public:

  Genes(const typename F::Input &i, bool fEvaluated = true) : fGenes(i) {
    if (fEvaluated)
      Evaluate();
  };
  //~Genes() {}
  //Genes(Genes const &);
  //void operator=(Genes const &);
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

}

#endif
