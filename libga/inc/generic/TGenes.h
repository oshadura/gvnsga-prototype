#ifndef __GENES__
#define __GENES__

#include <vector>
#include <memory>

#include "generic/ExceptionMessenger.h"

namespace geantvmoop {

template <typename F> class Genes;
template <typename F> using individual_t = std::shared_ptr<Genes<F> >;
template <typename F> class Genes {

private:
  typename F::Input fInput;
  typename F::Output fOutput;

public:
  Genes() {}

  Genes(const typename F::Input &i, bool fEvaluated = true) : fInput(i) {
    if (fEvaluated)
      Evaluate();
  };
  //~Genes() {}
  // Genes(Genes const &);
  // void operator=(Genes const &);
  bool IsDominating(const Genes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (fOutput[i] > other.fOutput[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsDominated(const Genes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (fOutput[i] < other.fOutput[i])
        return false;
    }
    return !IsEqual(other);
  }

  bool IsEqual(const Genes &other) const {
    for (unsigned int i = 0; i < GetOutput().size(); ++i) {
      if (other.fOutput[i] != fOutput[i])
        return false;
    }
    return true;
  }

  void Evaluate() { fOutput = F::Evaluate(fInput); }

  const typename F::Input &GetInput() const { return fInput; }

  const typename F::Output &GetOutput() const { return fOutput; }
};
}

#endif
